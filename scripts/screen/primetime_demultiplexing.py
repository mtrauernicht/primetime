# ==============================================================================
# Prime Time: TF reporter pipeline
# Vinícius H. Franceschini-Santos, Max Trauernicht 2026-01-13
# Version 0.2
# ==============================================================================
# Description:
#
# This script extracts barcodes from fastq files and counts the number of valid
# barcodes. It uses regex to match the upstream sequence of the barcode and
# allows for a specified number of mismatches.
#
# The barcodes extracted in this step were added during reverse transcription.
# These barcodes indicate the identity of the sample. So, in this pipeline step,
# we will demultiplex the reads based on these barcodes. In the next step (not in
# this script), the internal barcodes (added during PCR) will be used to identify 
# TF reporters.
#
# To be more precise, we will open all read 2 fastq.gz files provided (each 
# corresponding to a 384-well plate), and for each read, we will extract the barcode
# located right after a constant 28bp upstream sequence at the beginning of the read
# (with some allowed mismatches). We will just extract the 10bp following this 
# upstream sequence, which should correspond to the RT barcode. For each of the 
# 384-well plates, we will generate 384 new fastq.gz files, each corresponding to a well. 
# The rt_barcode.csv file will be used to map the extracted barcodes to well IDs.
# Then, the corresponding read 1 (which will be used in the next step) will be
# written to the appropriate well-specific fastq.gz file. At the end of the script,
# we will provide statistics on the number of reads processed, the number of
# valid barcodes, invalid barcodes, and mismatches. Also, we will output a file
# (SAMPLES_FILE) that indicates which samples (wells) correspond to which file 
# paths (for use in the next step of the pipeline). Important: samples of different
# plates can be of the same condition (replicates across plates).

# We also include basic plots to visualize the number of reads per sample (well)
# per 384-well plate.

# ==============================================================================
# Versions:
# 0.1 - Initial version
# 0.2 - Performance fixes:
#         * Dropped regex.BESTMATCH: it forces an exhaustive search for the
#           optimal fuzzy match rather than stopping at the first match within
#           the allowed error count, and was the dominant per-read cost.
#         * Switched batch processing from ThreadPoolExecutor to
#           ProcessPoolExecutor: the per-batch work is CPU-bound (regex +
#           string ops), so threads were serialized by the GIL and --threads
#           > 1 wasn't buying real parallelism before.
#         * Replaced Biopython SeqIO parsing/formatting with manual 4-line
#           FASTQ parsing, avoiding SeqRecord construction/formatting
#           overhead per read.
# ==============================================================================

import regex
from Bio.Seq import Seq
import argparse
import sys
import gzip
import os
import re
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

# ==============================================================================
# Parse arguments
# ==============================================================================
def parse_arguments():
    """
    Parse command-line arguments
    """
    parser = argparse.ArgumentParser(description="Demultiplex fastq.gz files based on RT barcodes")
    # Read containing the upstream sequence and barcode (typically read 2)
    parser.add_argument("--fastq", "--fastq_r2", dest="fastq_r2", type=str, required=True, help="Input R2 fastq file (contains RT barcode)")
    # Corresponding read 1 file whose records will be written to well-specific fastqs
    # For read1-only demultiplexing, this can be omitted and the --fastq_r2 input will be used.
    parser.add_argument("--read1", type=str, required=False, default=None, help="Input R1 fastq file to split into well-specific fastqs")
    parser.add_argument("--bc_length", type=int, required=True, help="Length of barcode to extract after the upstream sequence")
    parser.add_argument("--rt_bc_upstream_seq", type=str, required=True, help="Sequence upstream of barcode")
    parser.add_argument("--max_mismatch", type=int, default=0, help="Maximum number of mismatches allowed when matching upstream sequence")
    parser.add_argument("--threads", type=int, default=1, help="Number of worker processes to use for parallel batch processing")

    parser.add_argument("--invalid_bc_file", type=str, default=None, help="Optional output file for invalid barcodes and unmapped reads")
    parser.add_argument("--rt_barcode_map", type=str, required=True, help="CSV file mapping RT barcodes to well IDs (columns: well,barcode)")
    parser.add_argument("--samples_file", type=str, required=True, help="CSV file mapping well IDs to sample names (columns: well,sample)")
    parser.add_argument("--out_base", type=str, required=True, help="Base directory where demultiplexed files will be written (a 'demultiplex' dir will be created inside)")
    parser.add_argument("--plate_id", type=str, default=None, help="Plate ID to place demultiplexed well files under (out_base/demultiplex/<plate_id>)")
    parser.add_argument("--sample_id", type=str, default=None, help="Sample ID to filter demultiplexing output (only write fastq for this sample)")
    parser.add_argument("--barcode_ids_file", type=str, default=None, help="Optional output file to write barcode\tread_id\twell records for matched reads")
    parser.add_argument("--samples_output", type=str, default=None, help="Optional CSV output listing wells and sample names (well,sample) for upstream steps")
    parser.add_argument("--plot_dir", type=str, default=None, help="Directory where per-plate heatmap PDF should be written (optional)")
    parser.add_argument("--stats_file", type=str, default=None, help="Optional output file to write demultiplexing stats (in addition to stderr)")
    parser.add_argument("--rt_cluster_file", type=str, default=None, help="Optional starcode cluster file mapping raw RT barcodes to centroid/correct barcodes")

    return parser.parse_args()

# ==============================================================================
# Build regular expression pattern
# ==============================================================================
def build_regexp_pattern(rt_bc_upstream_seq, max_mismatch):
    """
    Build regular expression pattern for the upstream sequence of the barcode
    """
    pattern = "(" + rt_bc_upstream_seq + "){e<" + str(max_mismatch) + "}"
    return pattern


def count_reads_in_fastq(fastq_gz_path):
    """Quick count of reads in gzipped FASTQ by counting @ lines."""
    import subprocess
    try:
        result = subprocess.run(
            f"zcat {fastq_gz_path} | grep -c '^@'",
            shell=True,
            capture_output=True,
            text=True,
            timeout=300
        )
        if result.returncode == 0:
            return int(result.stdout.strip())
    except Exception as e:
        sys.stderr.write(f"Warning: could not count reads in {fastq_gz_path}: {e}\n")
    return None


def fastq_records(fh):
    """Yield (header, seq, plus, qual) 4-line tuples from an open FASTQ text handle.

    Each element still includes its trailing newline, so the four lines can be
    concatenated directly for output without any reformatting (unlike
    Biopython's SeqRecord.format(), which re-serializes the record).
    """
    while True:
        header = fh.readline()
        if not header:
            return
        seq = fh.readline()
        plus = fh.readline()
        qual = fh.readline()
        if not qual:
            return  # truncated/incomplete final record
        yield header, seq, plus, qual

# ==============================================================================
# Demultiplex and write per-well read1 files
# ==============================================================================

def build_barcode_map(mapping_file):
    """Read the RT barcode mapping CSV (columns: well,barcode) and return
    a dict barcode -> well."""
    import csv

    mapping = {}
    with open(mapping_file, "r", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            well = row.get("well") or row.get("Well") or row.get("WELL")
            barcode = row.get("barcode")
            if well and barcode:
                bc = barcode.strip()
                bc_rc = str(Seq(bc).reverse_complement())
                mapping[bc_rc] = well.strip()
    return mapping


def build_well_to_sample_map(samples_file):
    """Read the samples file (columns: well,sample) and return
    a dict well -> sample."""
    import csv

    mapping = {}
    with open(samples_file, "r", encoding="utf-8-sig") as fh:
        reader = csv.reader(fh)
        for row in reader:
            if len(row) >= 2:
                well = row[0].strip()
                sample = row[1].strip()
                if well and sample:
                    mapping[well] = sample
    return mapping


# ------------------------------------------------------------------------------
# Module-level worker state for ProcessPoolExecutor
# ------------------------------------------------------------------------------
# process_batch is defined at module level (rather than nested inside
# get_barcode_counts, as it was before) so it can be pickled and sent to
# worker processes. _init_worker populates _worker_state once per worker
# process (and once in the main process, for the single-process fallback
# path where no executor is used) so the regex pattern and lookup dicts
# don't need to be re-sent with every batch.
_worker_state = {}


def _init_worker(regex_pattern, bc_length, barcode_to_well, well_to_sample, cluster_map):
    _worker_state["compiled_pattern"] = regex.compile(regex_pattern)
    _worker_state["bc_length"] = bc_length
    _worker_state["barcode_to_well"] = barcode_to_well
    _worker_state["well_to_sample"] = well_to_sample
    _worker_state["cluster_map"] = cluster_map


def process_batch(batch):
    """Process a batch of records. batch is list of (seq2, read_id, r1_fastq_str)."""
    compiled_pattern = _worker_state["compiled_pattern"]
    bc_length = _worker_state["bc_length"]
    barcode_to_well = _worker_state["barcode_to_well"]
    well_to_sample = _worker_state["well_to_sample"]
    cluster_map = _worker_state["cluster_map"]

    local_matched_valid = 0
    local_matched_invalid = 0
    local_mismatched = 0
    local_invalids = []
    local_well_reads = {}
    local_barcode_ids = []
    local_well_counts = {}
    local_sample_counts = {}

    for seq, read_id, r1_fastq in batch:
        match = compiled_pattern.search(seq)
        if match is None:
            local_mismatched += 1
            continue

        start_bc = match.span()[1]
        barcode = seq[start_bc:start_bc + bc_length]
        # correct barcode via cluster centroid if available
        if cluster_map:
            barcode = cluster_map.get(barcode, barcode)

        well = barcode_to_well.get(barcode)
        if not well:
            local_matched_invalid += 1
            local_invalids.append((barcode, read_id))
            continue

        sample = well_to_sample.get(well)
        if not sample:
            local_matched_invalid += 1
            local_invalids.append((barcode, read_id))
            continue

        local_matched_valid += 1
        local_well_reads.setdefault(well, []).append(r1_fastq)
        local_well_counts[well] = local_well_counts.get(well, 0) + 1
        local_sample_counts[sample] = local_sample_counts.get(sample, 0) + 1
        local_barcode_ids.append((barcode, read_id))

    return (local_matched_valid, local_matched_invalid, local_mismatched,
            local_invalids, local_well_reads, local_barcode_ids,
            local_well_counts, local_sample_counts)


def get_barcode_counts(fastq_r2,
                       fastq_r1,
                       bc_length,
                       rt_bc_upstream_seq,
                       max_mismatch,
                       num_cores,
                       invalids_file,
                       rt_barcode_map_file,
                       samples_file,
                       out_base,
                       plate_id=None,
                       barcode_ids_file=None,
                       samples_output=None,
                       plot_dir=None,
                       rt_cluster_file=None):
    """Stream through paired R2 (contains RT barcode) and R1 (reads to write).

    Process entire plate in one pass: read R2 once, extract all barcodes,
    and write all sample-specific R1 fastq files simultaneously.

    Special-case: if plate_id == 'pDNA', this file is a single pDNA fastq that
    should NOT be demultiplexed but copied directly to out_base/pDNA.fastq.gz
    and recorded in the samples_output if requested.
    """
    # If read1 is not provided, use fastq_r2 for both barcode extraction and output
    if fastq_r1 is None:
        fastq_r1 = fastq_r2

    # Handle pDNA short-circuit: copy input raw file into pDNA.fastq.gz
    if plate_id == "pDNA":
        os.makedirs(out_base, exist_ok=True)
        # pDNA only has a read1 fastq (specified in plates file) — prefer read1
        src = fastq_r1 if fastq_r1 and os.path.exists(fastq_r1) else fastq_r2
        dst = os.path.join(out_base, "pDNA.fastq.gz")
        if src and os.path.exists(src):
            shutil.copyfile(src, dst)
            if samples_output:
                with open(samples_output, "w") as fh:
                    fh.write(f"pDNA,{dst}\n")
            # pDNA isn't demultiplexed; return zeros and empty counts.
            return 0, 0, 0, 0, {}
        else:
            raise FileNotFoundError("pDNA plate specified but source pDNA fastq not found")

    regex_pattern = build_regexp_pattern(rt_bc_upstream_seq, max_mismatch)

    barcode_to_well = build_barcode_map(rt_barcode_map_file)
    well_to_sample = build_well_to_sample_map(samples_file)

    # Optional: load starcode cluster mapping raw_bc -> centroid_bc
    cluster_map = {}
    if rt_cluster_file and os.path.exists(rt_cluster_file):
        try:
            with open(rt_cluster_file, "r") as fh:
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split("\t")
                    centroid = parts[0]
                    # map centroid to itself
                    cluster_map[centroid] = centroid
                    if len(parts) >= 3:
                        members = parts[2]
                        # members may be comma-separated or space-separated
                        members = members.replace(",", " ")
                        for m in members.split():
                            cluster_map[m] = centroid
        except Exception as e:
            sys.stderr.write(f"Warning: failed to parse cluster file {rt_cluster_file}: {e}\n")

    # Populate worker state in the main process too, so the single-process
    # fallback path (num_cores <= 1) can call process_batch directly.
    _init_worker(regex_pattern, bc_length, barcode_to_well, well_to_sample, cluster_map)

    matched_but_invalid = 0
    matched_and_valid = 0
    mismatched = 0
    total_reads = 0
    invalids = []
    well_counts = {}
    sample_counts = {}

    # Dictionary to hold open file handles for each well
    out_writers = {}
    
    # Create output directory (use out_base directly)
    os.makedirs(out_base, exist_ok=True)

    # Materialize all expected per-well outputs up front so Snakemake can
    # observe the declared files even when a well receives zero reads.
    if plate_id:
        for well in well_to_sample:
            expected_path = os.path.join(out_base, f"{plate_id}_{well}.fastq.gz")
            with gzip.open(expected_path, "wt", compresslevel=1):
                pass

    # Optional barcode -> read id mapping file
    barcode_ids_fh = None
    if barcode_ids_file:
        barcode_ids_fh = open(barcode_ids_file, "w")
        barcode_ids_fh.write("barcode\tread_id\n")

    def handle_result(result):
        nonlocal matched_and_valid, matched_but_invalid, mismatched, total_reads
        (loc_valid, loc_invalid, loc_mismatched, loc_invalids,
         loc_well_reads, loc_barcode_ids, loc_well_counts, loc_sample_counts) = result

        matched_and_valid += loc_valid
        matched_but_invalid += loc_invalid
        mismatched += loc_mismatched
        invalids.extend(loc_invalids)

        for well, reads in loc_well_reads.items():
            if plate_id:
                out_path = os.path.join(out_base, f"{plate_id}_{well}.fastq.gz")
            else:
                out_path = os.path.join(out_base, f"{well}.fastq.gz")
            if out_path not in out_writers:
                out_writers[out_path] = gzip.open(out_path, "wt", compresslevel=1)
            out_writers[out_path].writelines(reads)

        if barcode_ids_fh and loc_barcode_ids:
            for bc, rid in loc_barcode_ids:
                barcode_ids_fh.write(f"{bc}\t{rid}\n")

        for w, c in loc_well_counts.items():
            well_counts[w] = well_counts.get(w, 0) + c
        for s, c in loc_sample_counts.items():
            sample_counts[s] = sample_counts.get(s, 0) + c

    batch_size = 50000
    futures = []
    executor = ProcessPoolExecutor(
        max_workers=num_cores,
        initializer=_init_worker,
        initargs=(regex_pattern, bc_length, barcode_to_well, well_to_sample, cluster_map),
    ) if num_cores and num_cores > 1 else None

    # Skip a second pass over the FASTQ just for ETA estimation.
    total_reads_expected = None

    try:
        with gzip.open(fastq_r2, "rt") as h2, gzip.open(fastq_r1, "rt") as h1:
            iter2 = fastq_records(h2)
            iter1 = fastq_records(h1)

            batch = []
            pbar = tqdm(total=total_reads_expected, desc=f"Demultiplexing {plate_id if plate_id else 'reads'}", unit=" reads", unit_scale=True)
            for (header2, seq2, plus2, qual2), (header1, seq1, plus1, qual1) in zip(iter2, iter1):
                total_reads += 1
                pbar.update(1)
                read_id = header2[1:].split()[0]
                r1_fastq = header1 + seq1 + plus1 + qual1
                batch.append((seq2.rstrip("\n"), read_id, r1_fastq))
                if len(batch) >= batch_size:
                    if executor:
                        futures.append(executor.submit(process_batch, batch))
                    else:
                        handle_result(process_batch(batch))
                    batch = []

            if batch:
                if executor:
                    futures.append(executor.submit(process_batch, batch))
                else:
                    handle_result(process_batch(batch))
            
            pbar.close()
            
            if executor:
                pbar_futures = tqdm(total=len(futures), desc="Processing batches", unit=" batch")
                for fut in as_completed(futures):
                    handle_result(fut.result())
                    pbar_futures.update(1)
                pbar_futures.close()
    finally:
        if executor:
            executor.shutdown(wait=True)
        for fh in out_writers.values():
            fh.close()

        if barcode_ids_fh:
            barcode_ids_fh.close()

    # Write invalids to file (only if specified)
    if invalids_file:
        invalids_dir = os.path.dirname(invalids_file)
        if invalids_dir:
            os.makedirs(invalids_dir, exist_ok=True)
        with open(invalids_file, "w") as fh:
            for barcode, read_id in invalids:
                fh.write(f"{barcode}\t{read_id}\n")

    # Plot per-well read counts as a heatmap (rows A-P, columns numeric)
    try:
        import numpy as np
        import matplotlib.pyplot as plt
    except Exception:
        sys.stderr.write("matplotlib/numpy not available, skipping per-well plot\n")
    else:
        # Build a fixed 384-well plate layout (A-P, 1-24)
        rows = [chr(i) for i in range(ord('A'), ord('P') + 1)]
        maxcol = 24
        mat = np.zeros((len(rows), maxcol), dtype=int)
        for well, count in well_counts.items():
            m = re.match(r"([A-Za-z]+)(\d+)$", well)
            if m:
                row_idx = ord(m.group(1)[0].upper()) - ord('A')
                col_idx = int(m.group(2)) - 1
                if 0 <= row_idx < len(rows) and 0 <= col_idx < maxcol:
                    mat[row_idx, col_idx] = count
        fig, ax = plt.subplots(figsize=(max(8, maxcol / 2), 6))
        im = ax.imshow(mat, aspect='auto', cmap='viridis', origin='upper')
        ax.set_xticks(range(maxcol))
        ax.set_xticklabels([str(i + 1) for i in range(maxcol)], rotation=90)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels(rows)
        ax.set_xlabel("Column")
        ax.set_ylabel("Row")
        title_plate = plate_id if plate_id else ""
        ax.set_title(f"Read counts per well - plate {title_plate}")
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("Valid reads")
        if plot_dir:
            os.makedirs(plot_dir, exist_ok=True)
            out_plot = os.path.join(plot_dir, f"{plate_id}_well_read_counts.pdf") if plate_id else os.path.join(plot_dir, "well_read_counts.pdf")
        else:
            out_plot = os.path.join(out_base, f"{plate_id}_well_read_counts.pdf") if plate_id else os.path.join(out_base, "well_read_counts.pdf")
        plt.tight_layout()
        plt.savefig(out_plot)
        plt.close()
        sys.stderr.write(f"Wrote per-well read-count heatmap: {out_plot}\n")

        # Also create a density plot showing distribution of read counts
        read_count_values = [count for count in well_counts.values() if count > 0]
        if read_count_values:
            fig, ax = plt.subplots(figsize=(10, 6))
            
            # Create density plot (histogram with KDE-like appearance)
            ax.hist(read_count_values, bins=50, density=True, alpha=0.6, color='#377EB8', edgecolor='black')
            
            # Add mean and median lines
            mean_count = np.mean(read_count_values)
            median_count = np.median(read_count_values)
            
            ax.axvline(median_count, color='red', linestyle='dashed', linewidth=2, label=f'Median: {int(median_count):,}')
            ax.axvline(mean_count, color='darkred', linestyle='dotted', linewidth=2, label=f'Mean: {int(mean_count):,}')
            
            ax.set_xlabel("Read Count per Well", fontsize=12)
            ax.set_ylabel("Density", fontsize=12)
            title_plate = plate_id if plate_id else ""
            ax.set_title(f"Distribution of Read Counts per Well - plate {title_plate}", fontsize=14, fontweight='bold')
            ax.legend(fontsize=10)
            ax.grid(axis='y', alpha=0.3)
            
            # Format x-axis with thousands separator
            ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x):,}'))
            
            if plot_dir:
                out_density_plot = os.path.join(plot_dir, f"{plate_id}_read_count_distribution.pdf") if plate_id else os.path.join(plot_dir, "read_count_distribution.pdf")
            else:
                out_density_plot = os.path.join(out_base, f"{plate_id}_read_count_distribution.pdf") if plate_id else os.path.join(out_base, "read_count_distribution.pdf")
            
            plt.tight_layout()
            plt.savefig(out_density_plot)
            plt.close()
            sys.stderr.write(f"Wrote read-count density plot: {out_density_plot}\n")

    return matched_and_valid, matched_but_invalid, mismatched, total_reads, {}

# ==============================================================================
# Write statistics
# ==============================================================================
def write_stats(matched_and_valid, matched_but_invalid, mismatched, total_reads, stats_file=None):
    def _write_message(message):
        sys.stderr.write(message)
        if stats_file:
            stats_dir = os.path.dirname(stats_file)
            if stats_dir:
                os.makedirs(stats_dir, exist_ok=True)
            with open(stats_file, "w") as fh:
                fh.write(message)

    if total_reads == 0:
        _write_message("No reads processed.\n")
        return

    mismatched_pct = mismatched / total_reads * 100
    matched_and_valid_pct = matched_and_valid / total_reads * 100
    matched_but_invalid_pct = matched_but_invalid / total_reads * 100

    message = (
        f"Total reads: {total_reads}, of those: \n"
        f"    {matched_and_valid} ({matched_and_valid_pct:.2f}%) matched the upstream sequence and showed a valid barcode -- KEPT\n"
        f"    {matched_but_invalid} ({matched_but_invalid_pct:.2f}%) matched the upstream sequence but showed an invalid barcode* -- DISCARDED \n"
        f"    {mismatched} ({mismatched_pct:.2f}%) did not match the upstream sequence -- DISCARDED \n\n"
        f"* Invalid barcodes are those that are unmapped\n"
    )
    _write_message(message)

# ==============================================================================
# Main function
# ==============================================================================
def main():
    args = parse_arguments()
    matched_and_valid, matched_but_invalid, mismatched, total_reads, well_counts = (
        get_barcode_counts(
            args.fastq_r2,
            args.read1,
            args.bc_length,
            args.rt_bc_upstream_seq,
            args.max_mismatch,
            args.threads,
            args.invalid_bc_file,
            args.rt_barcode_map,
            args.samples_file,
            args.out_base,
            plate_id=args.plate_id,
            barcode_ids_file=args.barcode_ids_file,
            samples_output=args.samples_output,
            plot_dir=args.plot_dir,
            rt_cluster_file=args.rt_cluster_file,
        )
    )
    write_stats(matched_and_valid, matched_but_invalid, mismatched, total_reads, stats_file=args.stats_file)

    # Report per-well counts
    if well_counts:
        sys.stderr.write("Per-well read counts:\n")
        for well, count in sorted(well_counts.items()):
            sys.stderr.write(f"  {well}: {count}\n")

if __name__ == "__main__":
    main()