#!/usr/bin/env python3
"""
snp_density_scan.py — Sliding-window SNP density scan for alignment outlier detection

Divides an aligned FASTA into windows (default: alignment_length / 100) and computes
per-sequence SNP density in each window relative to the column consensus. Flags sequences
whose SNP density in any window exceeds a Z-score threshold — a signature of whole-genome
contamination, chimeric sequences, or misaligned regions.

Usage:
    python snp_density_scan.py -i alignment.fasta [options]

Required:
    -i / --input        Aligned FASTA (single marker or concatenated supermatrix)

Optional:
    -w / --window       Window size in bp (default: alignment_length / 100)
    -z / --zscore       Z-score threshold for flagging (default: 3.0)
    -o / --output       Output TSV file (default: snp_density_scan.tsv)
    -f / --flagged      Output file listing flagged sequences only (default: snp_density_flagged.tsv)
    -m / --min_cov      Minimum fraction of non-gap sites per window to include in stats
                        (default: 0.5 — skip windows where >50% of sites are gaps)
    -v / --verbose      Print per-window stats for flagged sequences

Output columns (snp_density_scan.tsv):
    taxon | window_start | window_end | snp_count | window_non_gap_sites |
    snp_density | window_mean | window_sd | zscore | FLAGGED

Output columns (snp_density_flagged.tsv):
    taxon | n_flagged_windows | max_zscore | max_zscore_window | snp_density_at_max |
    total_windows_tested | flag_rate_pct | recommendation

Notes:
    - SNP density = (sites where seq differs from consensus) / (non-gap sites in window)
    - Consensus at each column = most frequent non-gap base; ties broken arbitrarily
    - Windows with < min_cov non-gap coverage are excluded from per-window Z-score stats
    - A sequence flagged in many windows suggests systematic data quality issues
    - Run AFTER trimAl and BEFORE (or during) tree inference
"""

import argparse
import statistics
import sys
from collections import Counter


# ── Argument parsing ──────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description="Sliding-window SNP density scan for alignment outlier detection"
    )
    p.add_argument("-i", "--input", required=True, help="Aligned FASTA file")
    p.add_argument(
        "-w", "--window", type=int, default=None,
        help="Window size in bp (default: alignment_length / 100)"
    )
    p.add_argument(
        "-z", "--zscore", type=float, default=3.0,
        help="Z-score threshold for flagging a window (default: 3.0)"
    )
    p.add_argument(
        "-o", "--output", default="snp_density_scan.tsv",
        help="Output TSV with per-window per-sequence stats"
    )
    p.add_argument(
        "-f", "--flagged", default="snp_density_flagged.tsv",
        help="Output TSV listing only flagged sequences"
    )
    p.add_argument(
        "-m", "--min_cov", type=float, default=0.5,
        help="Min fraction of non-gap sites per window to include in stats (default: 0.5)"
    )
    p.add_argument(
        "-v", "--verbose", action="store_true",
        help="Print per-window details for flagged sequences to stdout"
    )
    return p.parse_args()


# ── FASTA reader ──────────────────────────────────────────────────────────────
def read_fasta(path):
    seqs = {}
    header = None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                header = line[1:].split()[0]
                seqs[header] = []
            elif header is not None:
                seqs[header].append(line.upper())
    return {h: "".join(parts) for h, parts in seqs.items()}


# ── Consensus column ──────────────────────────────────────────────────────────
def column_consensus(col_chars):
    """Most frequent non-gap, non-N base in column. Returns '-' if none."""
    counts = Counter(c for c in col_chars if c not in ("-", "?", "N"))
    if not counts:
        return "-"
    return counts.most_common(1)[0][0]


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()

    print(f"Reading alignment: {args.input}")
    seqs = read_fasta(args.input)
    if not seqs:
        print("ERROR: No sequences found in input file.", file=sys.stderr)
        sys.exit(1)

    taxa = list(seqs.keys())
    aln_len = max(len(s) for s in seqs.values())
    n_taxa = len(taxa)

    # Pad sequences shorter than alignment length
    seqs = {h: s.ljust(aln_len, "-") for h, s in seqs.items()}

    # Determine window size
    window_size = args.window if args.window else max(10, aln_len // 100)
    print(f"Alignment length: {aln_len} bp | Taxa: {n_taxa} | Window size: {window_size} bp")
    print(f"Z-score threshold: {args.zscore} | Min coverage: {args.min_cov*100:.0f}%")
    print()

    # Build column consensus array (fast)
    print("Computing column consensus...")
    consensus = []
    for i in range(aln_len):
        col = [seqs[t][i] for t in taxa]
        consensus.append(column_consensus(col))

    # ── Sliding window scan ───────────────────────────────────────────────────
    windows = []
    start = 0
    while start < aln_len:
        end = min(start + window_size, aln_len)
        windows.append((start, end))
        start = end

    print(f"Scanning {len(windows)} windows...")

    # Per-window, per-taxon SNP density
    # all_rows: list of (taxon, win_start, win_end, snp_count, non_gap, density)
    all_rows = []

    for win_start, win_end in windows:
        win_consensus = consensus[win_start:win_end]
        win_len = win_end - win_start

        densities = {}
        for taxon in taxa:
            seq_win = seqs[taxon][win_start:win_end]
            non_gap = sum(1 for c in seq_win if c not in ("-", "?", "N"))
            if non_gap < args.min_cov * win_len:
                continue  # Too gappy — skip this taxon for this window
            snp_count = sum(
                1 for c, cons in zip(seq_win, win_consensus)
                if c not in ("-", "?", "N") and cons not in ("-", "?", "N") and c != cons
            )
            density = snp_count / non_gap if non_gap > 0 else 0.0
            densities[taxon] = (snp_count, non_gap, density)

        if len(densities) < 3:
            continue  # Not enough taxa to compute meaningful Z-scores

        dens_values = [d for _, _, d in densities.values()]
        win_mean = statistics.mean(dens_values)
        win_sd = statistics.stdev(dens_values) if len(dens_values) > 1 else 0.0

        for taxon, (snp_count, non_gap, density) in densities.items():
            zscore = (density - win_mean) / win_sd if win_sd > 0 else 0.0
            all_rows.append((taxon, win_start + 1, win_end, snp_count, non_gap,
                             density, win_mean, win_sd, zscore))

    # ── Write full output ─────────────────────────────────────────────────────
    print(f"Writing full results to: {args.output}")
    with open(args.output, "w") as out:
        out.write("taxon\twindow_start\twindow_end\tsnp_count\tnon_gap_sites\t"
                  "snp_density\twindow_mean\twindow_sd\tzscore\tFLAGGED\n")
        for row in all_rows:
            taxon, ws, we, snp, ng, dens, wmean, wsd, z = row
            flagged = "YES" if z >= args.zscore else ""
            out.write(
                f"{taxon}\t{ws}\t{we}\t{snp}\t{ng}\t"
                f"{dens:.4f}\t{wmean:.4f}\t{wsd:.4f}\t{z:.2f}\t{flagged}\n"
            )

    # ── Summarize flagged sequences ───────────────────────────────────────────
    flagged_summary = {}  # taxon -> list of (win_start, win_end, density, zscore)
    for row in all_rows:
        taxon, ws, we, snp, ng, dens, wmean, wsd, z = row
        if z >= args.zscore:
            if taxon not in flagged_summary:
                flagged_summary[taxon] = []
            flagged_summary[taxon].append((ws, we, dens, z))

    # Total windows tested per taxon
    windows_tested = {}
    for row in all_rows:
        taxon = row[0]
        windows_tested[taxon] = windows_tested.get(taxon, 0) + 1

    print(f"\nWriting flagged sequences to: {args.flagged}")
    with open(args.flagged, "w") as out:
        out.write("taxon\tn_flagged_windows\tmax_zscore\tmax_zscore_window\t"
                  "snp_density_at_max\ttotal_windows_tested\tflag_rate_pct\trecommendation\n")
        for taxon in sorted(flagged_summary, key=lambda t: -max(z for _, _, _, z in flagged_summary[t])):
            flags = flagged_summary[taxon]
            n_flag = len(flags)
            max_z_entry = max(flags, key=lambda x: x[3])
            max_z = max_z_entry[3]
            max_win = f"{max_z_entry[0]}-{max_z_entry[1]}"
            max_dens = max_z_entry[2]
            total_tested = windows_tested.get(taxon, 0)
            flag_rate = n_flag / total_tested * 100 if total_tested > 0 else 0.0

            if flag_rate >= 20:
                rec = "INVESTIGATE: high flag rate — possible whole-genome/chimeric sequence"
            elif n_flag >= 3:
                rec = "INVESTIGATE: multiple flagged windows — possible partial contamination"
            else:
                rec = "REVIEW: single/few flagged windows — may be genuine rate variation"

            out.write(f"{taxon}\t{n_flag}\t{max_z:.2f}\t{max_win}\t"
                      f"{max_dens:.4f}\t{total_tested}\t{flag_rate:.1f}\t{rec}\n")

            if args.verbose:
                print(f"  {taxon}: {n_flag} flagged windows (max Z={max_z:.2f} at {max_win})")

    # ── Console summary ───────────────────────────────────────────────────────
    n_flagged_taxa = len(flagged_summary)
    print(f"\n{'='*60}")
    print(f"SUMMARY")
    print(f"{'='*60}")
    print(f"Alignment: {aln_len} bp | {n_taxa} taxa | {len(windows)} windows ({window_size} bp each)")
    print(f"Z-score threshold: {args.zscore}")
    print(f"Flagged taxa: {n_flagged_taxa} / {n_taxa}")

    if flagged_summary:
        print(f"\nTop flagged sequences (by max Z-score):")
        print(f"  {'Taxon':<50} {'FlagWin':>7} {'MaxZ':>7} {'FlagRate%':>10} {'Recommendation'}")
        print(f"  {'-'*50} {'-'*7} {'-'*7} {'-'*10} {'-'*40}")
        for taxon in sorted(flagged_summary, key=lambda t: -max(z for _, _, _, z in flagged_summary[t]))[:20]:
            flags = flagged_summary[taxon]
            max_z = max(z for _, _, _, z in flags)
            total_tested = windows_tested.get(taxon, 0)
            flag_rate = len(flags) / total_tested * 100 if total_tested > 0 else 0
            if flag_rate >= 20:
                rec = "INVESTIGATE (high rate)"
            elif len(flags) >= 3:
                rec = "INVESTIGATE (multiple)"
            else:
                rec = "REVIEW"
            print(f"  {taxon:<50} {len(flags):>7} {max_z:>7.2f} {flag_rate:>10.1f}   {rec}")
    else:
        print("No flagged sequences — all taxa within Z-score threshold.")

    print(f"\nFull results: {args.output}")
    print(f"Flagged summary: {args.flagged}")


if __name__ == "__main__":
    main()
