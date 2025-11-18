#!/usr/bin/env python3
"""
extract_clipped_seq.py

Read a regional WGS BAM file and extract clipped / inserted sequences with read IDs.
minimum inserted/clipped sequence length 600

Usage:
  python extract_clipped_seq.py insertion_bearing_regional.bam --min-length 600 --fasta -o inserted_reads.fa

"""

import argparse
import sys
import pysam

CIGAR_OPS = {
    0: "M",
    1: "I",
    2: "D",
    3: "S",
    4: "H",
    5: "=",
}

def extract_from_bam(bam_path, out_handle, min_length=1, fasta=False):
    # iterate regional bam and fetch ins seq
    bam = pysam.AlignmentFile(bam_path, "rb")

    for read in bam.fetch(until_eof=True):
        if read.is_unmapped:
            continue

        seq = read.query_sequence
        if seq is None:
            continue

        cigartuples = read.cigartuples
        if not cigartuples:
            continue

        read_pos = 0  # pointer on read.sequence (0-based)
        read_len = len(seq)

        for (op, length) in cigartuples:
            op_char = CIGAR_OPS.get(op, str(op))

            # Deletions
            if op == 2:
                # read through
                continue
            # Hard clip
            if op == 4:
                # read through as well
                continue

            # potential insertion Cigar values I & S
            if op in (1, 3):  # I or S
                if length >= min_length:
                    start = read_pos
                    end = read_pos + length
                    # sanity check bounds
                    if start < 0 or end > read_len:
                        pass
                    else:
                        subseq = seq[start:end]
                        # checking direction
                        # soft clipped at left or right
                        side = "."
                        if op == 3:  # S
                            if start == 0:
                                side = "left"
                            elif end == read_len:
                                side = "right"
                            else:
                                # I in the middle
                                side = "internal"
                        elif op == 1:  # insertion
                            side = "internal"

                        if fasta:
                            header = f">{read.query_name}|{op_char}|{side}|{start}:{end}|len={length}"
                            out_handle.write(header + "\n")
                            out_handle.write(subseq + "\n")
                        else:
                            out_handle.write("\t".join([
                                read.query_name,
                                op_char,
                                side,
                                str(start),
                                str(length),
                                subseq
                            ]) + "\n")

                read_pos += length
            else:
                # exact match count in
                if op in (0, 5):  # M, =
                    read_pos += length
                else:
                    read_pos += length

    bam.close()


def main():
    p = argparse.ArgumentParser(
        description="Extract soft-clipped and inserted sequences from a regional WGS BAM file."
    )
    p.add_argument("bam", help="Input BAM")
    p.add_argument("-o", "--output", help="Output file", default=None)
    p.add_argument("--min-length", type=int, default=600,
                   help="Minimum length of clipped/inserted sequence (default 600).")
    p.add_argument("--fasta", action="store_true", help="Write output in FASTA instead of TSV.")
    args = p.parse_args()

    if args.output:
        out_handle = open(args.output, "w")
    else:
        out_handle = sys.stdout

    try:
        extract_from_bam(args.bam, out_handle, min_length=args.min_length, fasta=args.fasta)
    finally:
        if args.output:
            out_handle.close()


if __name__ == "__main__":
    main()
