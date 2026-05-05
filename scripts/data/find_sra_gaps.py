#!/usr/bin/env python3
"""
find_sra_gaps.py — Identify Zingiberaceae species absent from Plan B but present in SRA

Queries NCBI SRA for WGS and Targeted-Capture runs matching the target taxonomic group,
compares the SRA species list against the Plan B GenBank dataset (post-dedup), and
outputs a prioritized list of SRA runs to download for gap-filling.

Gap-filling strategy:
    - A "gap species" is one with SRA genomic data (WGS or Targeted-Capture) but NO
      sequences in the current Plan B marker dataset (data/planB/genbank/).
    - For each gap species, one SRA run is selected:
        Priority 1: WGS run (best for GetOrganelle plastome assembly)
        Priority 2: Targeted-Capture run (acceptable; yields partial plastome)
    - Within a priority tier, prefer the run with the most bases (richest data).
    - The selected runs are saved to data/planA/ for plastome assembly.
    - Assembled plastomes are then processed with extract_cds.py → output merged
      into data/planB/genbank/<marker>/ to fill the Plan B dataset.

Usage:
    python find_sra_gaps.py \\
        --group         "Zingiberaceae" \\
        --genbank_dir   data/planB/genbank \\
        --markers       matK,rbcL,trnL,psbA-trnH,rpoB,rpoC1,atpB,ndhF,ycf1,ycf2 \\
        --sra_out       data/planA/sra_gaps.tsv \\
        [--min_bases    100000000] \\
        [--strategy     WGS,Targeted-Capture] \\
        [--provenance   results/planA/provenance] \\
        [--dry_run]

Arguments:
    --group         Taxonomic group to query in SRA (must match NCBI organism name)
    --genbank_dir   Root dir of Plan B genbank downloads (one subdir per marker)
    --markers       Comma-separated marker names to check coverage for
    --sra_out       Output TSV of gap-filling SRA runs (one run per gap species)
    --min_bases     Minimum base count for a run to be considered (default: 100,000,000)
                    Set lower for targeted-capture runs (~50M bases typical)
    --strategy      Comma-separated SRA library strategies to query (default: WGS,Targeted-Capture)
    --provenance    Directory to write JSON provenance log
    --dry_run       Print results without writing files

Output TSV columns:
    species             Binomial (Genus_species)
    sra_run             Best SRA run ID (SRR/DRR/ERR)
    strategy            WGS or Targeted-Capture
    bases               Total bases in the run
    size_mb             Estimated file size
    gap_markers         Markers this species is missing from Plan B
    destination         data/planA/ subdirectory for this download
    assembly_tool       Recommended assembly tool (GetOrganelle for both strategies)

Tool version requirements:
    Python >= 3.8
    Entrez Direct (esearch, efetch) — must be on PATH

Downstream commands (after running this script):
    # Download SRA runs listed in sra_out
    bash $PHYLOSKILLS_ROOT/scripts/data/download_sra.sh \\
        -l data/planA/sra_gaps.tsv \\
        -o data/planA/sra_raw \\
        -P results/planA/provenance

    # Assemble each downloaded run with GetOrganelle
    bash $PHYLOSKILLS_ROOT/scripts/assembly/assemble_plastome_getorganelle.sh \\
        -i data/planA/sra_raw \\
        -o data/planA/plastomes \\
        -P results/planA/provenance

    # Extract CDS from assembled plastomes → merge into Plan B
    python $PHYLOSKILLS_ROOT/scripts/assembly/extract_cds.py \\
        --input    data/planA/plastomes \\
        --output   data/planA/cds \\
        --markers  matK,rbcL,trnL,psbA-trnH,rpoB,rpoC1,atpB,ndhF,ycf1,ycf2 \\
        --provenance results/planA/provenance

    # Merge extracted CDS into Plan B genbank directories
    python $PHYLOSKILLS_ROOT/scripts/data/merge_sra_cds_to_planB.py \\
        --cds_dir   data/planA/cds \\
        --planB_dir data/planB/genbank_dedup \\
        --markers   matK,rbcL,trnL,psbA-trnH,rpoB,rpoC1,atpB,ndhF,ycf1,ycf2
"""

import argparse
import csv
import io
import os
import subprocess
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
        description="Find gap species (in SRA but not in Plan B) for gap-filling."
    )
    p.add_argument("--group",       required=True,
                   help='Taxonomic group, e.g. "Zingiberaceae"')
    p.add_argument("--genbank_dir", required=True,
                   help="Root dir of Plan B genbank downloads (one subdir per marker)")
    p.add_argument("--markers",     required=True,
                   help="Comma-separated marker names")
    p.add_argument("--sra_out",     required=True,
                   help="Output TSV of gap SRA runs")
    p.add_argument("--min_bases",   type=int, default=100_000_000,
                   help="Min bases for a run to qualify (default: 100M = ~1× genome coverage)")
    p.add_argument("--strategy",    default="WGS,Targeted-Capture",
                   help="Comma-separated SRA library strategies (default: WGS,Targeted-Capture)")
    p.add_argument("--provenance",  default=None,
                   help="Directory for JSON provenance log")
    p.add_argument("--dry_run",     action="store_true",
                   help="Print results without writing output files")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Plan B species scanning
# ---------------------------------------------------------------------------

def get_planB_species(genbank_dir, markers):
    """
    Return a dict: species_binomial → set(markers present).
    Scans genbank_dir/<marker>/ for FASTA files named Genus_species_*.fasta.
    """
    species_markers = defaultdict(set)
    for marker in markers:
        marker_dir = os.path.join(genbank_dir, marker)
        if not os.path.isdir(marker_dir):
            continue
        for fname in os.listdir(marker_dir):
            if not fname.endswith(".fasta"):
                continue
            parts = fname.split("_")
            if len(parts) < 3:
                continue
            binomial = f"{parts[0]}_{parts[1]}"
            species_markers[binomial].add(marker)
    return species_markers


# ---------------------------------------------------------------------------
# SRA querying via Entrez Direct
# ---------------------------------------------------------------------------

def run_cmd(cmd, description=""):
    """Run a shell command and return stdout as string. Exits on error."""
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR running: {description or cmd}", file=sys.stderr)
        print(result.stderr[:500], file=sys.stderr)
        sys.exit(1)
    return result.stdout.strip()


def fetch_sra_runinfo(group, strategies):
    """
    Fetch SRA runinfo CSV for the given taxonomic group and library strategies.
    Returns list of dicts (one per SRA run).
    """
    strategy_query = " OR ".join(f'"{s}"[Strategy]' for s in strategies)
    query = f'"{group}"[Organism] AND ("GENOMIC"[Source]) AND ({strategy_query})'

    print(f"  SRA query: {query}")
    cmd = (
        f'esearch -db sra -query \'{query}\' '
        f'| efetch -format runinfo 2>/dev/null'
    )
    raw = run_cmd(cmd, "esearch+efetch SRA runinfo")
    if not raw:
        print("  WARNING: No SRA results returned.", file=sys.stderr)
        return []

    lines = raw.splitlines()
    if not lines:
        return []

    # Parse CSV
    reader = csv.DictReader(io.StringIO(raw))
    runs = []
    for row in reader:
        # Filter: skip empty rows and runs with missing key fields
        if not row.get("Run") or not row.get("ScientificName"):
            continue
        # Parse bases as int
        try:
            bases = int(row.get("bases", 0))
        except ValueError:
            bases = 0
        row["_bases_int"] = bases
        try:
            size_mb = int(row.get("size_MB", 0))
        except ValueError:
            size_mb = 0
        row["_size_mb_int"] = size_mb
        runs.append(row)

    print(f"  SRA runs fetched: {len(runs)}")
    return runs


def build_sra_species_map(runs, strategies, min_bases):
    """
    For each species, collect all qualifying runs by strategy.
    Returns: species_binomial → {strategy → [run_dict sorted by bases desc]}
    Normalizes species names: replaces spaces with underscores.
    """
    species_map = defaultdict(lambda: defaultdict(list))

    for run in runs:
        sci_name = run.get("ScientificName", "").strip()
        if not sci_name or sci_name.lower() in ("", "nan"):
            continue
        # Normalize to Genus_species (handle variety/subspecies names)
        parts = sci_name.split()
        if len(parts) < 2:
            continue
        binomial = f"{parts[0]}_{parts[1]}"

        strategy = run.get("LibraryStrategy", "").strip()
        if strategy not in strategies:
            continue

        bases = run["_bases_int"]
        if bases < min_bases:
            continue

        species_map[binomial][strategy].append(run)

    return species_map


def select_best_run(strat_dict, priority_order):
    """
    Select the best single SRA run for a gap species.
    Follows priority_order: first strategy that has qualifying runs wins.
    Within that strategy, picks the run with the most bases.
    Returns (run_dict, strategy_str) or (None, None).
    """
    for strategy in priority_order:
        runs = strat_dict.get(strategy, [])
        if runs:
            return max(runs, key=lambda r: r["_bases_int"]), strategy
    return None, None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    start    = time.time()
    args     = parse_args()
    markers  = [m.strip() for m in args.markers.split(",")]
    strategies = [s.strip() for s in args.strategy.split(",")]
    # Priority order: WGS first (works better with GetOrganelle), then remaining strategies
    priority = (["WGS"] if "WGS" in strategies else []) + [s for s in strategies if s != "WGS"]

    print("find_sra_gaps.py")
    print(f"  Group:         {args.group}")
    print(f"  Plan B dir:    {args.genbank_dir}")
    print(f"  Markers:       {', '.join(markers)}")
    print(f"  Strategies:    {', '.join(strategies)} (priority: {' > '.join(priority)})")
    print(f"  Min bases:     {args.min_bases:,}")
    print(f"  Output:        {args.sra_out}")
    print()

    # 1. Get Plan B species and their marker coverage
    print("Scanning Plan B dataset for existing species...")
    planB_species = get_planB_species(args.genbank_dir, markers)
    print(f"  Species in Plan B: {len(planB_species)}")
    print()

    # 2. Fetch SRA runinfo
    print("Querying NCBI SRA...")
    runs = fetch_sra_runinfo(args.group, strategies)
    if not runs:
        print("No SRA results. Exiting.")
        return 0
    print()

    # 3. Build SRA species map (qualifying runs only)
    print(f"Building SRA species map (min_bases={args.min_bases:,})...")
    sra_species_map = build_sra_species_map(runs, strategies, args.min_bases)
    print(f"  SRA species with qualifying runs: {len(sra_species_map)}")
    print()

    # 4. Find gap species: in SRA but NOT (or only partially) in Plan B
    print("Identifying gap species (absent from Plan B)...")
    gap_species = {}
    for sra_sp, strat_dict in sra_species_map.items():
        if sra_sp not in planB_species:
            # Completely absent from Plan B
            gap_species[sra_sp] = {
                "missing_markers": markers,    # all markers missing
                "sra_strategies":  strat_dict,
            }

    print(f"  Gap species (in SRA, absent from Plan B): {len(gap_species)}")
    print()

    if not gap_species:
        print("No gap species found. Plan B already covers all SRA species.")
        return 0

    # 5. Select best SRA run per gap species (accumulate by_strategy during selection)
    print("Selecting best SRA run per gap species...")
    output_rows = []
    by_strategy = defaultdict(int)
    for species, info in sorted(gap_species.items()):
        best_run, chosen_strategy = select_best_run(info["sra_strategies"], priority)
        if best_run is None:
            continue
        by_strategy[chosen_strategy] += 1
        output_rows.append({
            "species":       species,
            "sra_run":       best_run.get("Run", ""),
            "strategy":      chosen_strategy,
            "bases":         best_run["_bases_int"],
            "size_mb":       best_run["_size_mb_int"],
            "gap_markers":   ",".join(info["missing_markers"]),
            "destination":   "data/planA/sra_raw",
            "assembly_tool": "GetOrganelle",
        })

    print(f"  Gap runs selected: {len(output_rows)}")
    for strat, cnt in sorted(by_strategy.items()):
        print(f"    {strat}: {cnt} runs")
    print()

    # 6. Print preview
    print("Preview (first 20 gap species):")
    print(f"  {'Species':<35} {'Run':<15} {'Strategy':<20} {'Bases (M)':>10}")
    print(f"  {'-'*35} {'-'*15} {'-'*20} {'-'*10}")
    for row in output_rows[:20]:
        print(f"  {row['species']:<35} {row['sra_run']:<15} {row['strategy']:<20} "
              f"{row['bases']//1_000_000:>9}M")
    if len(output_rows) > 20:
        print(f"  ... and {len(output_rows) - 20} more (see {args.sra_out})")
    print()

    # 7. Write output TSV
    if not args.dry_run:
        os.makedirs(os.path.dirname(args.sra_out) or ".", exist_ok=True)
        fieldnames = ["species", "sra_run", "strategy", "bases", "size_mb",
                      "gap_markers", "destination", "assembly_tool"]
        with open(args.sra_out, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(output_rows)
        print(f"  Gap runs written: {args.sra_out}")

        # Also write a plain run-ID list for download_sra.sh -l flag
        run_list_path = args.sra_out.replace(".tsv", "_runlist.txt")
        with open(run_list_path, "w") as fh:
            for row in output_rows:
                fh.write(row["sra_run"] + "\n")
        print(f"  Run ID list:      {run_list_path}")
    else:
        print("  [dry_run] No files written.")

    # 8. Provenance
    if args.provenance and not args.dry_run:
        write_provenance(
            script="find_sra_gaps",
            tool="python",
            version=sys.version.split()[0],
            provenance_dir=args.provenance,
            parameters={
                "group":       args.group,
                "genbank_dir": args.genbank_dir,
                "markers":     markers,
                "strategies":  strategies,
                "min_bases":   args.min_bases,
                "sra_out":     args.sra_out,
                "planB_species":  len(planB_species),
                "sra_species":    len(sra_species_map),
                "gap_species":    len(gap_species),
                "runs_selected":  len(output_rows),
                "by_strategy":    dict(by_strategy),
            },
            input_files={"genbank_dir": args.genbank_dir},
            output_files={"sra_out": args.sra_out},
            runtime_seconds=int(time.time() - start),
            exit_code=0,
        )

    print()
    print("Next steps:")
    print(f"  1. Review gap list: {args.sra_out}")
    print(f"  2. Download SRA runs:")
    print(f"       bash $PHYLOSKILLS_ROOT/scripts/data/download_sra.sh \\")
    print(f"            -l {args.sra_out.replace('.tsv','_runlist.txt')} \\")
    print(f"            -o data/planA/sra_raw")
    print(f"  3. Assemble plastomes:")
    print(f"       bash $PHYLOSKILLS_ROOT/scripts/assembly/assemble_plastome_getorganelle.sh \\")
    print(f"            -i data/planA/sra_raw -o data/planA/plastomes")
    print(f"  4. Extract CDS and merge into Plan B genbank_dedup/")
    print()
    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
