#!/usr/bin/env python3
"""
change_header_name.py — Rename FASTA sequence headers from a mapping file

Replaces sequence headers (>old_name) with new names from a two-column TSV
mapping file. Useful for standardizing headers to Genus_species_accession format
before alignment or tree inference.

Usage:
    python change_header_name.py \
        --input <input.fasta> \
        --map <mapping.tsv> \
        --output <output.fasta> \
        [--column1 <old_name_col>] \
        [--column2 <new_name_col>] \
        [--strict]

Arguments:
    --input    Input FASTA file
    --map      Tab-separated mapping file (two columns: old_name, new_name)
               Lines starting with # are ignored.
               Extra columns beyond the first two are ignored.
    --output   Output FASTA file with renamed headers
    --column1  Column index (0-based) for old name (default: 0)
    --column2  Column index (0-based) for new name (default: 1)
    --strict   Exit with error if any sequence in the FASTA has no mapping
               (default: warn and keep original name)

Notes:
    - Only the portion up to the first space in the FASTA header is matched
      against the mapping file (the accession/ID part).
    - The full header line (including description after the space) is replaced
      by the new name only; description is dropped.

Tool version requirements:
    Python >= 3.8  (no external dependencies)
"""

import argparse
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Rename FASTA headers from a mapping file"
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--map", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--column1", type=int, default=0,
                        help="Column index for old name (default: 0)")
    parser.add_argument("--column2", type=int, default=1,
                        help="Column index for new name (default: 1)")
    parser.add_argument("--strict", action="store_true",
                        help="Exit with error on missing mapping")
    return parser.parse_args()


def load_mapping(map_file: str, col1: int, col2: int) -> dict[str, str]:
    mapping = {}
    with open(map_file) as fh:
        for line_num, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) <= max(col1, col2):
                print(f"WARNING: Line {line_num} in mapping file has fewer than "
                      f"{max(col1, col2) + 1} columns — skipped.", file=sys.stderr)
                continue
            old_name = parts[col1].strip()
            new_name = parts[col2].strip()
            if old_name in mapping and mapping[old_name] != new_name:
                print(f"WARNING: Duplicate mapping for '{old_name}' "
                      f"('{mapping[old_name]}' vs '{new_name}') — using last.",
                      file=sys.stderr)
            mapping[old_name] = new_name
    return mapping


def main():
    args = parse_args()
    print(f"# change_header_name.py")
    print(f"# Input:  {args.input}")
    print(f"# Map:    {args.map}")
    print(f"# Output: {args.output}")

    mapping = load_mapping(args.map, args.column1, args.column2)
    print(f"# Loaded {len(mapping)} name mappings")

    n_renamed = 0
    n_kept = 0
    n_missing = 0

    with open(args.input) as fin, open(args.output, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                # Extract ID (portion before first space)
                header_rest = line[1:].rstrip()
                seq_id = header_rest.split()[0]
                if seq_id in mapping:
                    fout.write(f">{mapping[seq_id]}\n")
                    n_renamed += 1
                else:
                    if args.strict:
                        print(f"ERROR: No mapping found for '{seq_id}' (--strict mode)",
                              file=sys.stderr)
                        sys.exit(1)
                    print(f"WARNING: No mapping for '{seq_id}' — keeping original",
                          file=sys.stderr)
                    fout.write(line)
                    n_missing += 1
                    n_kept += 1
            else:
                fout.write(line)

    print(f"# Renamed: {n_renamed}")
    print(f"# Kept (no mapping): {n_kept}")
    print(f"# Output: {args.output}")


if __name__ == "__main__":
    main()
