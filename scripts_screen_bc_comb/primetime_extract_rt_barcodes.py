# ==============================================================================
# Extract RT barcodes from R2 and write counts
# ==============================================================================
# Usage: python primetime_extract_rt_barcodes.py --fastq_r2 <file> --bc_length 10 \
#               --rt_bc_upstream_seq <seq> --max_mismatch 0 --output <counts.txt>
# ==============================================================================

import argparse
import gzip
import os
import sys
from Bio import SeqIO
import regex


def parse_args():
    p = argparse.ArgumentParser(description="Extract RT barcodes from R2 and write counts (barcode\tcount)")
    p.add_argument("--fastq_r2", required=True, help="Input R2 fastq.gz (contains RT barcode)")
    p.add_argument("--bc_length", type=int, required=True, help="Barcode length to extract")
    p.add_argument("--rt_bc_upstream_seq", required=True, help="Sequence upstream of barcode")
    p.add_argument("--max_mismatch", type=int, default=0, help="Max mismatches for upstream matching")
    p.add_argument("--output", required=True, help="Output path for counts (barcode\tcount)")
    return p.parse_args()


def build_pattern(rt_bc_upstream_seq, max_mismatch):
    return "(" + rt_bc_upstream_seq + "){e<" + str(max_mismatch) + "}"


def main():
    args = parse_args()
    pattern = build_pattern(args.rt_bc_upstream_seq, args.max_mismatch)
    compiled = regex.compile(pattern, regex.BESTMATCH)

    counts = {}
    total = 0

    with gzip.open(args.fastq_r2, "rt") as h2:
        for rec in SeqIO.parse(h2, "fastq"):
            total += 1
            seq = str(rec.seq)
            m = compiled.search(seq)
            if m is None:
                continue
            start_bc = m.span()[1]
            bc = seq[start_bc:start_bc + args.bc_length]
            if len(bc) == args.bc_length:
                counts[bc] = counts.get(bc, 0) + 1

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.output, "w") as fh:
        for bc, c in counts.items():
            fh.write(f"{bc}\t{c}\n")

    sys.stderr.write(f"Extracted {len(counts)} unique RT barcodes from {total} reads.\n")


if __name__ == "__main__":
    main()
