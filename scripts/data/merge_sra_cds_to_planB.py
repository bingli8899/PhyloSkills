#!/usr/bin/env python3
"""
merge_sra_cds_to_planB.py — Merge SRA-extracted marker sequences into the
Plan B genbank_dedup directory.

Scans data/planA/cds/<marker>/ for FASTA files matching the Plan B naming
convention (Genus_species_<run>_<marker>.fasta) and copies them into
data/planB/genbank_dedup/<marker>/. Skips files that already exist in the
destination (idempotent — safe to re-run).

Usage:
    python merge_sra_cds_to_planB.py \\
        --cds_dir   data/planA/cds \\
        --planB_dir data/planB/genbank_dedup \\
        --markers   matK,rbcL,trnL,psbA-trnH,rpoB,rpoC1,atpB,ndhF,ycf1,ycf2 \\
        [--dry_run] \\
        [--provenance results/planA/provenance]

Arguments:
    --cds_dir     Root directory of extracted marker sequences (planA/cds)
    --planB_dir   Root of Plan B deduplicated dataset (genbank_dedup)
    --markers     Comma-separated marker names to merge
    --dry_run     Print actions without copying files
    --provenance  Directory for JSON provenance log

Tool version requirements:
    Python >= 3.8  (standard library only)
"""

import argparse
import os
import shutil
import sys
import time
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils"))
from provenance_py import write_provenance


def parse_args():
    p = argparse.ArgumentParser(
        description="Copy SRA-extracted marker sequences into Plan B genbank_dedup."
    )
    p.add_argument("--cds_dir",   required=True,
                   help="Root of extracted sequences (data/planA/cds)")
    p.add_argument("--planB_dir", required=True,
                   help="Root of Plan B genbank_dedup")
    p.add_argument("--markers",   required=True,
                   help="Comma-separated marker names")
    p.add_argument("--dry_run",   action="store_true",
                   help="Print actions without writing")
    p.add_argument("--provenance", default=None,
                   help="Provenance directory")
    return p.parse_args()


def main():
    start   = time.time()
    args    = parse_args()
    markers = [m.strip() for m in args.markers.split(",")]

    print("merge_sra_cds_to_planB.py")
    print(f"  Source:  {args.cds_dir}")
    print(f"  Dest:    {args.planB_dir}")
    print(f"  Markers: {', '.join(markers)}")
    if args.dry_run:
        print("  [dry_run] No files will be copied.")
    print()

    n_copied  = 0
    n_skipped = 0   # already exists
    n_missing = 0   # cds_dir marker subdir absent
    by_marker = defaultdict(int)

    for marker in markers:
        src_dir = os.path.join(args.cds_dir, marker)
        dst_dir = os.path.join(args.planB_dir, marker)

        if not os.path.isdir(src_dir):
            print(f"  {marker:<14} — no extracted sequences (source dir absent)")
            n_missing += 1
            continue

        fasta_files = [f for f in os.listdir(src_dir) if f.endswith(".fasta")]
        if not fasta_files:
            print(f"  {marker:<14} — no FASTA files in source dir")
            n_missing += 1
            continue

        if not args.dry_run:
            os.makedirs(dst_dir, exist_ok=True)

        for fname in sorted(fasta_files):
            src_path = os.path.join(src_dir, fname)
            dst_path = os.path.join(dst_dir, fname)

            if os.path.exists(dst_path):
                n_skipped += 1
                continue

            if args.dry_run:
                print(f"  [dry_run] copy {fname} → {dst_dir}/")
            else:
                shutil.copy2(src_path, dst_path)

            n_copied += 1
            by_marker[marker] += 1

        print(f"  {marker:<14} {by_marker.get(marker, 0)} files merged")

    print()
    print(f"Total copied:  {n_copied}")
    print(f"Total skipped (already exist): {n_skipped}")
    print(f"Markers with no source: {n_missing}")

    if args.provenance and not args.dry_run:
        write_provenance(
            script="merge_sra_cds_to_planB",
            tool="python",
            version=sys.version.split()[0],
            provenance_dir=args.provenance,
            parameters={
                "cds_dir":    args.cds_dir,
                "planB_dir":  args.planB_dir,
                "markers":    markers,
                "n_copied":   n_copied,
                "n_skipped":  n_skipped,
                "by_marker":  dict(by_marker),
            },
            input_files={"cds_dir": args.cds_dir},
            output_files={"planB_dir": args.planB_dir},
            runtime_seconds=int(time.time() - start),
            exit_code=0,
        )

    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
