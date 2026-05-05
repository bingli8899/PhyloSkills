#!/usr/bin/env python3
"""
select_best_accession.py — Deduplicate multi-marker GenBank downloads to one accession per species

For each species, selects the single accession that appears in the most markers
(i.e., the accession sourced from the richest dataset, such as a complete plastome).
Ties are broken randomly. The selected accession's sequences are copied into a
cleaned output directory; sequences from other accessions for the same species
are discarded.

Usage:
    python select_best_accession.py \\
        --genbank_dir  <data/planB/genbank> \\
        --output_dir   <data/planB/genbank_dedup> \\
        --markers      matK,rbcL,trnL,psbA-trnH,rpoB,rpoC1,atpB,ndhF,ycf1,ycf2 \\
        [--report      <data/planB/accession_selection.tsv>] \\
        [--seed        42] \\
        [--provenance  <results/planB/provenance>]

Arguments:
    --genbank_dir   Root directory containing one subdir per marker
                    (e.g. genbank/matK/, genbank/rbcL/, ...)
    --output_dir    Output root; same subdir layout as genbank_dir but deduplicated
    --markers       Comma-separated list of markers to process (must match subdir names)
    --report        TSV file recording selection decisions (optional but recommended)
    --seed          Random seed for tie-breaking (default: 42, for reproducibility)
    --provenance    Directory for JSON provenance log (optional)

Selection algorithm:
    1. Parse every FASTA filename: Genus_species_accession_marker.fasta
    2. Build: species → accession → set(markers where that accession has a file)
    3. For each species, choose the accession covering the most markers.
       On tie: randomly select one (controlled by --seed).
    4. Copy the chosen accession's files to output_dir/<marker>/.
    5. Markers where the chosen accession has NO file → left absent (missing data).

Filename convention (input and output):
    <Genus>_<species>_<accession>_<marker>.fasta
    e.g. Zingiber_officinale_NC_037455.1_matK.fasta

Provenance output (JSON):
    {
      "script": "select_best_accession",
      "date": "...",
      "parameters": {...},
      "summary": {
        "total_species": N,
        "species_with_multiple_accessions": N,
        "total_files_selected": N,
        "total_files_discarded": N
      }
    }

Tool version requirements:
    Python >= 3.8  (standard library only)
"""

import argparse
import os
import random
import shutil
import sys
import time
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils"))
from provenance_py import write_provenance


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Deduplicate GenBank downloads: one representative accession per species."
    )
    p.add_argument("--genbank_dir", required=True,
                   help="Root dir with one subdir per marker (e.g. data/planB/genbank)")
    p.add_argument("--output_dir", required=True,
                   help="Output root dir; will be created if absent")
    p.add_argument("--markers", required=True,
                   help="Comma-separated marker names matching subdir names")
    p.add_argument("--report", default=None,
                   help="Path for TSV selection report (recommended)")
    p.add_argument("--seed", type=int, default=42,
                   help="Random seed for reproducible tie-breaking (default: 42)")
    p.add_argument("--provenance", default=None,
                   help="Directory to write JSON provenance log")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def parse_filename(fname, marker):
    """
    Parse Genus_species_accession_marker.fasta → (binomial, accession).
    The accession is everything between species and marker fields.
    Handles accessions with dots and underscores (e.g. NC_050165.1).
    Returns (None, None) if filename is unparseable.
    """
    if not fname.endswith(".fasta"):
        return None, None
    stem = fname[:-6]  # strip .fasta
    # Expected tail: _<marker>
    marker_suffix = f"_{marker}"
    if not stem.endswith(marker_suffix):
        return None, None
    core = stem[: -len(marker_suffix)]          # Genus_species_accession
    parts = core.split("_")
    if len(parts) < 3:
        return None, None
    genus   = parts[0]
    species = parts[1]
    # Accession = everything from parts[2] onward, re-joined with "_"
    accession = "_".join(parts[2:])
    binomial = f"{genus}_{species}"
    return binomial, accession


def scan_genbank_dir(genbank_dir, markers):
    """
    Returns:
        species_data  dict: species → accession → {marker: filepath}
        all_files     list of all FASTA filepaths found
    """
    # species → accession → {marker → filepath}
    species_data = defaultdict(lambda: defaultdict(dict))
    all_files = []

    for marker in markers:
        marker_dir = os.path.join(genbank_dir, marker)
        if not os.path.isdir(marker_dir):
            print(f"  WARNING: marker directory not found: {marker_dir}", file=sys.stderr)
            continue
        for fname in os.listdir(marker_dir):
            if not fname.endswith(".fasta"):
                continue
            fpath = os.path.join(marker_dir, fname)
            all_files.append(fpath)
            binomial, accession = parse_filename(fname, marker)
            if binomial is None:
                print(f"  WARNING: could not parse filename: {fname}", file=sys.stderr)
                continue
            species_data[binomial][accession][marker] = fpath

    return species_data, all_files


def select_best_accessions(species_data, rng):
    """
    For each species, pick the accession covering the most markers.
    Ties broken randomly using rng.

    Returns:
        selections   dict: species → {
            "chosen_accession": str,
            "marker_count": int,
            "all_accessions": {accession: marker_count},
            "files": {marker: filepath},
            "tie": bool
        }
        n_discarded  int: total files belonging to non-representative accessions
    """
    selections  = {}
    n_discarded = 0

    for species, accessions in species_data.items():
        # accession → number of markers it covers
        acc_counts = {acc: len(markers) for acc, markers in accessions.items()}
        max_count  = max(acc_counts.values())
        best_accs  = [acc for acc, cnt in acc_counts.items() if cnt == max_count]
        is_tie     = len(best_accs) > 1

        chosen = rng.choice(best_accs)

        for acc, markers in accessions.items():
            if acc != chosen:
                n_discarded += len(markers)

        selections[species] = {
            "chosen_accession": chosen,
            "marker_count":     max_count,
            "all_accessions":   acc_counts,
            "files":            accessions[chosen],   # {marker → filepath}
            "tie":              is_tie,
        }

    return selections, n_discarded


def copy_selected_files(selections, output_dir):
    """
    Copy selected FASTA files to output_dir/<marker>/<filename>.
    Returns n_copied.
    """
    # Pre-create all marker output directories before copying
    markers_seen = {marker for info in selections.values() for marker in info["files"]}
    for marker in markers_seen:
        os.makedirs(os.path.join(output_dir, marker), exist_ok=True)

    n_copied = 0
    for species, info in selections.items():
        for marker, src_path in info["files"].items():
            dst_path = os.path.join(output_dir, marker, os.path.basename(src_path))
            shutil.copy2(src_path, dst_path)
            n_copied += 1

    return n_copied


def write_report(selections, report_path):
    """Write a TSV report of selection decisions."""
    os.makedirs(os.path.dirname(report_path) or ".", exist_ok=True)
    with open(report_path, "w") as fh:
        fh.write("species\tchosen_accession\tmarkers_covered\ttie\tall_accessions_and_counts\n")
        for species in sorted(selections):
            info = selections[species]
            # Format all accessions: acc1(N),acc2(M),...
            all_acc_str = ",".join(
                f"{acc}({cnt})"
                for acc, cnt in sorted(info["all_accessions"].items(),
                                       key=lambda x: -x[1])
            )
            fh.write(
                f"{species}\t{info['chosen_accession']}\t"
                f"{info['marker_count']}\t{info['tie']}\t{all_acc_str}\n"
            )
    print(f"  Selection report written: {report_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    start  = time.time()
    args   = parse_args()
    markers = [m.strip() for m in args.markers.split(",")]
    rng    = random.Random(args.seed)

    print(f"select_best_accession.py")
    print(f"  Input:   {args.genbank_dir}")
    print(f"  Output:  {args.output_dir}")
    print(f"  Markers: {', '.join(markers)}")
    print(f"  Seed:    {args.seed}")
    print()

    # 1. Scan input
    print("Scanning downloaded files...")
    species_data, all_files = scan_genbank_dir(args.genbank_dir, markers)
    n_species = len(species_data)
    print(f"  Found {len(all_files)} FASTA files across {n_species} species")

    multi_acc = sum(1 for acc in species_data.values() if len(acc) > 1)
    print(f"  Species with multiple accessions (require dedup): {multi_acc}")
    print()

    # 2. Select best accession per species (n_discarded accumulated during selection)
    print("Selecting best accession per species (most markers covered)...")
    selections, n_discarded = select_best_accessions(species_data, rng)

    ties = sum(1 for info in selections.values() if info["tie"])
    print(f"  Tie-breaking applied for {ties} species (random, seed={args.seed})")
    print()

    # 3. Copy selected files
    print(f"Copying selected files to: {args.output_dir}")
    os.makedirs(args.output_dir, exist_ok=True)
    n_copied = copy_selected_files(selections, args.output_dir)
    print(f"  Files copied (selected):    {n_copied}")
    print(f"  Files discarded (non-repr): {n_discarded}")
    print()

    # 4. Per-marker summary (derived from selections dict, no re-scan)
    print("Per-marker sequence counts (deduplicated):")
    marker_counts = defaultdict(int)
    for info in selections.values():
        for marker in info["files"]:
            marker_counts[marker] += 1
    for marker in markers:
        print(f"  {marker:<14} {marker_counts.get(marker, 0)}")
    print()

    # 5. Report
    if args.report:
        write_report(selections, args.report)

    # 6. Provenance
    summary = {
        "total_species":                   n_species,
        "species_with_multiple_accessions": multi_acc,
        "tie_broken_species":              ties,
        "total_files_selected":            n_copied,
        "total_files_discarded":           n_discarded,
    }
    if args.provenance:
        write_provenance(
            script="select_best_accession",
            tool="python",
            version=sys.version.split()[0],
            provenance_dir=args.provenance,
            parameters={
                "genbank_dir": args.genbank_dir,
                "output_dir":  args.output_dir,
                "markers":     markers,
                "seed":        args.seed,
                "report":      args.report,
                **summary,
            },
            input_files={"genbank_dir": args.genbank_dir},
            output_files={"output_dir": args.output_dir, "report": args.report or ""},
            runtime_seconds=int(time.time() - start),
            exit_code=0,
        )

    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
