# ==============================================================================
# Prime Time: TF reporter pipeline
# Vinícius H. Franceschini-Santos, Max Trauernicht 2024-10-22
# Version 0.1
# ==============================================================================
# Description:
#
# This script extracts barcodes from fastq files and counts the number of valid
# barcodes. It uses regex to match the downstream sequence of the barcode and
# allows for a specified number of mismatches.
#
# ==============================================================================
# Versions:
# 0.1 - Initial version
# ==============================================================================

import regex
from Bio import SeqIO
import argparse
import sys
import gzip
from concurrent.futures import ProcessPoolExecutor, as_completed

# ==============================================================================
# Parse arguments
# ==============================================================================
def parse_arguments():
    """
    Parse command-line arguments
    """
    parser = argparse.ArgumentParser(description="Extract barcode counts from fastq.gz")
    parser.add_argument("--fastq", type=str, help="Input fastq file")
    parser.add_argument("--bc_length", type=int, help="Length of barcode")
    parser.add_argument(
        "--bc_downstream_seq", type=str, help="Sequence downstream of barcode"
    )
    parser.add_argument(
        "--max_mismatch", type=int, help="Maximum number of mismatches of the downstream allowed"
    )
    parser.add_argument(
        "--threads", type=int, default=20, help="Number of threads to use"
    )
    parser.add_argument(
        "--invalid_bc_file", type=str, help="Output file for invalid barcodes"
    )

    return parser.parse_args()

# ==============================================================================
# Build regular expression pattern
# ==============================================================================
def build_regexp_pattern(bc_downstream_seq, max_mismatch):
    """
    Build regular expression pattern for the downstream sequence of the barcode
    """
    pattern = "(" + bc_downstream_seq + "){e<" + str(max_mismatch) + "}"
    return pattern

# ==============================================================================
# Process chunk
# ==============================================================================
def process_chunk(chunk, regex_pattern, bc_length):
    matched_but_invalid = 0
    matched_and_valid = 0
    mismatched = 0
    total_reads = 0
    results = []
    invalids = []

    for record in chunk:
        total_reads += 1
        seq = str(record.seq)
        match = regex.search(regex_pattern, seq, regex.BESTMATCH)
        if match is None:
            mismatched += 1
            continue
        end_bc = match.span()[0]
        barcode = seq[0:end_bc]
        results.append(barcode)
        matched_and_valid += 1
        # if (len(barcode) >= bc_length) and ("N" not in barcode):
        #     results.append(barcode)
        #     matched_and_valid += 1
        # else:
        #     invalids.append(barcode)
        #     matched_but_invalid += 1

    return matched_and_valid, matched_but_invalid, mismatched, total_reads, results, invalids

# ==============================================================================
# Get barcode counts
# ==============================================================================
def get_barcode_counts(fastq, bc_length, bc_downstream_seq, max_mismatch, num_cores, invalids_file):
    regex_pattern = build_regexp_pattern(bc_downstream_seq, max_mismatch)
    matched_but_invalid = 0
    matched_and_valid = 0
    mismatched = 0
    total_reads = 0
    results = []

    with gzip.open(fastq, "rt") as handle:
        records = list(SeqIO.parse(handle, "fastq"))
        chunk_size = 1000  # Set chunk size to 1000
        chunks = [records[i:i + chunk_size] for i in range(0, len(records), chunk_size)]

    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        futures = [executor.submit(process_chunk, chunk, regex_pattern, bc_length) for chunk in chunks]
        for future in as_completed(futures):
            try:
                m_valid, m_invalid, mm, t_reads, res, invalids = future.result()
                matched_and_valid += m_valid
                matched_but_invalid += m_invalid
                mismatched += mm
                total_reads += t_reads
                results.extend(res)
                invalids.extend(invalids)
                # sys.stdout.write(f"Processed chunk: {t_reads} reads\n")
            except Exception as e:
                sys.stderr.write(f"Error processing chunk: {e}\n")

    for barcode in results:
        sys.stdout.write(barcode + "\n")

    with open(invalids_file, "w") as f:
        for barcode in invalids:
            f.write(barcode + "\n")

    return matched_and_valid, matched_but_invalid, mismatched, total_reads

# ==============================================================================
# Write statistics
# ==============================================================================
def write_stats(matched_and_valid, matched_but_invalid, mismatched, total_reads):
    mismatched_pct = mismatched / total_reads * 100
    matched_and_valid_pct = matched_and_valid / total_reads * 100
    matched_but_invalid_pct = matched_but_invalid / total_reads * 100
    sys.stderr.write(f"Total reads: {total_reads}, of those: \n")
    sys.stderr.write(
        f"    {matched_and_valid} ({matched_and_valid_pct:.2f}%) matched the downstream sequence and showed a valid barcode -- KEPT\n"
    )
    sys.stderr.write(
        f"    {matched_but_invalid} ({matched_but_invalid_pct:.2f}%) matched the downstream sequence but showed an invalid barcode* -- DISCARDED \n"
    )
    sys.stderr.write(
        f"    {mismatched} ({mismatched_pct:.2f}%) did not match the downstream sequence -- DISCARDED \n\n"
    )
    sys.stderr.write(
        f"* Invalid barcodes are those that are shorter than the expected length or contain Ns\n"
    )

# ==============================================================================
# Main function
# ==============================================================================
def main():
    args = parse_arguments()
    matched_and_valid, matched_but_invalid, mismatched, total_reads = (
        get_barcode_counts(
            args.fastq,
            args.bc_length, 
            args.bc_downstream_seq, 
            args.max_mismatch, 
            args.threads, 
            args.invalid_bc_file,
        )
    )
    write_stats(matched_and_valid, matched_but_invalid, mismatched, total_reads)

if __name__ == "__main__":
    main()