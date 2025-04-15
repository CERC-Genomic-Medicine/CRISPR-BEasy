#!/usr/bin/env python3

import argparse
from pathlib import Path

def write_chunk(lines, output_prefix, chunk_id):
    with open(f"chunk{chunk_id:03d}_{output_prefix}.bed", "w") as out:
        out.writelines(lines)

def split_long_entry(chrom, start, end, rest, soft_limit, overlap):
    segments = []
    pos = start
    while pos < end:
        split_end = min(pos + soft_limit, end)
        fields = [chrom, str(pos), str(split_end)] + rest
        segments.append("\t".join(fields) + "\n")
        pos += soft_limit - overlap
    return segments

def split_bed_by_bp(bed_path, soft_limit, overlap):
    chunk, chunk_bp, chunk_id = [], 0, 1
    stem = Path(bed_path).stem
    output_prefix = stem.replace(".bed", "").replace(".BED", "")

    with open(bed_path) as bed:
        for line in bed:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            try:
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                rest = parts[3:]
            except ValueError:
                continue

            length = end - start

            if length > soft_limit:
                sub_entries = split_long_entry(chrom, start, end, rest, soft_limit, overlap)
                for sub in sub_entries:
                    sub_len = int(sub.split()[2]) - int(sub.split()[1])
                    if chunk_bp + sub_len > soft_limit and chunk:
                        write_chunk(chunk, output_prefix, chunk_id)
                        chunk_id += 1
                        chunk, chunk_bp = [], 0
                    chunk.append(sub)
                    chunk_bp += sub_len
                continue

            if chunk_bp + length > soft_limit and chunk:
                write_chunk(chunk, output_prefix, chunk_id)
                chunk_id += 1
                chunk, chunk_bp = [], 0

            chunk.append(line)
            chunk_bp += length

    if chunk:
        write_chunk(chunk, output_prefix, chunk_id)

def main():
    parser = argparse.ArgumentParser(
        description="Split a BED file by soft limit on base pairs. If a record exceeds the limit, it will be split with overlap."
    )
    parser.add_argument(
        "bed_file",
        type=str,
        metavar="BED_FILE",
        help="Path to input BED file"
    )
    parser.add_argument(
        "soft_limit",
        type=int,
        metavar="SOFT_LIMIT",
        help="Soft limit (in base pairs) per output chunk"
    )
    parser.add_argument(
        "--overlap",
        type=int,
        metavar="BP",
        default=30,
        help="Overlap in bp when splitting large intervals (default: 0)"
    )

    args = parser.parse_args()
    split_bed_by_bp(
        bed_path=args.bed_file,
        soft_limit=args.soft_limit,
        overlap=args.overlap
    )

if __name__ == "__main__":
    main()

