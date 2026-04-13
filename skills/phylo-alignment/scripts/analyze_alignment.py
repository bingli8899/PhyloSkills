#!/usr/bin/env python3
"""
analyze_alignment.py — Alignment quality statistics and diagnostics

Computes per-marker and per-taxon statistics for aligned FASTA files.
Flags outlier sequences, high-gap columns, and low parsimony-informative sites.
Outputs a summary TSV and prints a human-readable report.

Usage:
    python analyze_alignment.py \
        --input <aligned.fasta or directory> \
        [--min_length <int>] \
        [--max_gap_pct <float>] \
        [--min_pis_pct <float>] \
        [--output <report.tsv>]

Arguments:
    --input        Aligned FASTA file OR directory with *_aligned.fasta files
    --min_length   Minimum sequence length to flag as short (default: 100 bp)
    --max_gap_pct  Maximum allowed gap % per sequence before flagging (default: 50.0)
    --min_pis_pct  Minimum parsimony-informative sites % to flag (default: 5.0)
    --output       Output TSV file path (default: alignment_diagnostics.tsv)

Tool version requirements:
    Biopython >= 1.83  — pip install biopython
    (Latest: 1.85, 2025)

Output columns in TSV:
    marker, n_seqs, alignment_length, mean_gap_pct, pis_count, pis_pct,
    n_outlier_seqs, flags, outlier_taxa
"""

import argparse
import sys
import os
from pathlib import Path
from collections import Counter

# ── Version check ────────────────────────────────────────────────────────────
try:
    from Bio import SeqIO, AlignIO
    from Bio.Align import MultipleSeqAlignment
    import Bio
    print(f"# Biopython version: {Bio.__version__}  (latest known: 1.85)")
except ImportError:
    print("ERROR: Biopython not found. Install: pip install biopython", file=sys.stderr)
    sys.exit(1)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Alignment quality diagnostics"
    )
    parser.add_argument("--input", required=True, help="FASTA file or directory")
    parser.add_argument("--min_length", type=int, default=100)
    parser.add_argument("--max_gap_pct", type=float, default=50.0)
    parser.add_argument("--min_pis_pct", type=float, default=5.0)
    parser.add_argument("--output", default="alignment_diagnostics.tsv")
    return parser.parse_args()


def find_fasta_files(input_path: str) -> list[Path]:
    p = Path(input_path)
    if p.is_file():
        return [p]
    elif p.is_dir():
        files = []
        for pat in ["*_aligned.fasta", "*_trimmed.fasta", "*.fasta", "*.fa"]:
            files.extend(sorted(p.glob(pat)))
        # Deduplicate (prefer _aligned over raw)
        seen = set()
        deduped = []
        for f in files:
            if f.stem not in seen:
                seen.add(f.stem)
                deduped.append(f)
        return deduped
    else:
        print(f"ERROR: Input not found: {input_path}", file=sys.stderr)
        sys.exit(1)


def analyze_alignment(fasta_path: Path, args) -> dict:
    """Compute statistics for a single aligned FASTA."""
    marker = fasta_path.stem.replace("_aligned", "").replace("_trimmed", "")

    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    if len(records) < 2:
        return {
            "marker": marker, "n_seqs": len(records),
            "alignment_length": 0, "mean_gap_pct": 0,
            "pis_count": 0, "pis_pct": 0,
            "n_outlier_seqs": 0,
            "flags": "TOO_FEW_SEQUENCES",
            "outlier_taxa": "",
        }

    seqs = {r.id: str(r.seq).upper() for r in records}
    aln_len = max(len(s) for s in seqs.values())
    n_seqs = len(seqs)

    # Per-sequence gap analysis
    outlier_taxa = []
    gap_pcts = []
    for taxon, seq in seqs.items():
        gap_count = seq.count("-") + seq.count("?")
        eff_len = len(seq)
        gap_pct = gap_count / eff_len * 100 if eff_len > 0 else 0
        gap_pcts.append(gap_pct)

        # Flag outliers
        seq_len = len(seq.replace("-", "").replace("?", "").replace("N", ""))
        if gap_pct > args.max_gap_pct:
            outlier_taxa.append(f"{taxon}(gap={gap_pct:.0f}%)")
        elif seq_len < args.min_length:
            outlier_taxa.append(f"{taxon}(len={seq_len}bp)")

    mean_gap_pct = sum(gap_pcts) / len(gap_pcts) if gap_pcts else 0

    # Parsimony informative sites
    seqs_list = list(seqs.values())
    pis = 0
    for i in range(aln_len):
        col = [s[i] for s in seqs_list if i < len(s)]
        col = [c for c in col if c not in ("-", "?", "N")]
        counts = Counter(col)
        n_variable_chars = sum(1 for c, n in counts.items() if n >= 2)
        if n_variable_chars >= 2:
            pis += 1

    pis_pct = pis / aln_len * 100 if aln_len > 0 else 0

    # Aggregate flags
    flags = []
    if mean_gap_pct > 30:
        flags.append(f"HIGH_GAP({mean_gap_pct:.0f}%)")
    if pis_pct < args.min_pis_pct:
        flags.append(f"LOW_PIS({pis_pct:.1f}%)")
    if outlier_taxa:
        flags.append(f"OUTLIER_SEQS({len(outlier_taxa)})")

    return {
        "marker": marker,
        "n_seqs": n_seqs,
        "alignment_length": aln_len,
        "mean_gap_pct": round(mean_gap_pct, 1),
        "pis_count": pis,
        "pis_pct": round(pis_pct, 1),
        "n_outlier_seqs": len(outlier_taxa),
        "flags": ",".join(flags) if flags else "PASS",
        "outlier_taxa": ";".join(outlier_taxa),
    }


def main():
    args = parse_args()
    print(f"# Date: {__import__('datetime').date.today()}")
    print(f"# max_gap_pct threshold: {args.max_gap_pct}%")
    print(f"# min_pis_pct threshold: {args.min_pis_pct}%")
    print(f"# min_length threshold:  {args.min_length} bp")
    print()

    files = find_fasta_files(args.input)
    print(f"Found {len(files)} alignment file(s) to analyze.")
    print()

    results = []
    for f in files:
        result = analyze_alignment(f, args)
        results.append(result)

    # Print human-readable report
    print("=" * 70)
    print(f"{'Marker':<25} {'Seqs':>5} {'Len':>7} {'Gap%':>6} {'PIS':>5} {'PIS%':>6}  Flags")
    print("-" * 70)
    n_pass = 0
    n_warn = 0
    for r in results:
        status = "PASS" if r["flags"] == "PASS" else r["flags"]
        print(
            f"{r['marker']:<25} {r['n_seqs']:>5} {r['alignment_length']:>7} "
            f"{r['mean_gap_pct']:>6.1f} {r['pis_count']:>5} {r['pis_pct']:>6.1f}  {status}"
        )
        if r["flags"] == "PASS":
            n_pass += 1
        else:
            n_warn += 1
            if r["outlier_taxa"]:
                for taxon in r["outlier_taxa"].split(";"):
                    if taxon:
                        print(f"  → Outlier: {taxon}")
    print("=" * 70)
    print(f"PASS: {n_pass}  |  Flagged: {n_warn}")
    print()

    # Recommend actions
    flagged = [r for r in results if r["flags"] != "PASS"]
    if flagged:
        print("Recommended actions:")
        for r in flagged:
            print(f"  {r['marker']}: {r['flags']}")
            if "HIGH_GAP" in r["flags"]:
                print("    → Re-trim with trimal -automated1, or remove partial sequences")
            if "LOW_PIS" in r["flags"]:
                print("    → Consider excluding marker or adding more divergent taxa")
            if "OUTLIER_SEQS" in r["flags"]:
                print(f"    → Investigate or remove: {r['outlier_taxa']}")
        print()

    # Write TSV
    output_path = Path(args.output)
    with open(output_path, "w") as fh:
        cols = ["marker", "n_seqs", "alignment_length", "mean_gap_pct",
                "pis_count", "pis_pct", "n_outlier_seqs", "flags", "outlier_taxa"]
        fh.write("\t".join(cols) + "\n")
        for r in results:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    print(f"Full report written to: {output_path}")
    print("Next step: phylo-model-selection (if all markers PASS or reviewed)")


if __name__ == "__main__":
    main()
