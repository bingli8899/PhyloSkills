#!/usr/bin/env python3
"""
remove_N_fasta.py — Filter or mask FASTA sequences with excessive ambiguous bases

Removes sequences exceeding a threshold of N/ambiguous characters, or optionally
replaces ambiguous regions with gaps for downstream alignment. Useful after
plastome assembly or consensus calling to remove low-quality sequences.

Usage:
    python remove_N_fasta.py \
        --input <input.fasta> \
        --output <output.fasta> \
        [--max_n_pct <float>] \
        [--max_gap_pct <float>] \
        [--min_length <int>] \
        [--report <report.tsv>]

Arguments:
    --input       Input FASTA file
    --output      Output FASTA file (filtered)
    --max_n_pct   Maximum percent of N/ambiguous bases allowed (default: 20.0)
    --max_gap_pct Maximum percent of gap characters (-) allowed (default: 50.0)
    --min_length  Minimum ungapped sequence length in bp (default: 100)
    --report      Optional: write per-sequence statistics to TSV

Ambiguous bases counted: N, R, Y, S, W, K, M, B, D, H, V (IUPAC)

Tool version requirements:
    Python >= 3.8  (no external dependencies)
"""

import argparse
import sys
from pathlib import Path

AMBIGUOUS = set("NRYSWKMBDHV")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter FASTA sequences with excessive N/ambiguous characters"
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--max_n_pct", type=float, default=20.0)
    parser.add_argument("--max_gap_pct", type=float, default=50.0)
    parser.add_argument("--min_length", type=int, default=100)
    parser.add_argument("--report", default="")
    return parser.parse_args()


def parse_fasta(path: str):
    """Yield (header, sequence) tuples from a FASTA file."""
    header = None
    seq_parts = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            elif line:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


def main():
    args = parse_args()
    print(f"# remove_N_fasta.py")
    print(f"# Input:       {args.input}")
    print(f"# max_n_pct:   {args.max_n_pct}%")
    print(f"# max_gap_pct: {args.max_gap_pct}%")
    print(f"# min_length:  {args.min_length} bp")

    report_rows = []
    n_kept = 0
    n_removed = 0

    with open(args.output, "w") as fout:
        for header, seq in parse_fasta(args.input):
            seq_upper = seq.upper()
            total_len = len(seq_upper)
            gap_count = seq_upper.count("-")
            ambig_count = sum(1 for c in seq_upper if c in AMBIGUOUS)
            ungapped_len = total_len - gap_count

            gap_pct = gap_count / total_len * 100 if total_len > 0 else 0
            n_pct = ambig_count / total_len * 100 if total_len > 0 else 0

            # Determine reason for removal
            reasons = []
            if n_pct > args.max_n_pct:
                reasons.append(f"N_pct={n_pct:.1f}%>{args.max_n_pct}%")
            if gap_pct > args.max_gap_pct:
                reasons.append(f"gap_pct={gap_pct:.1f}%>{args.max_gap_pct}%")
            if ungapped_len < args.min_length:
                reasons.append(f"ungapped_len={ungapped_len}<{args.min_length}")

            seq_id = header.split()[0]
            status = "removed" if reasons else "kept"
            report_rows.append({
                "id": seq_id,
                "total_len": total_len,
                "ungapped_len": ungapped_len,
                "gap_pct": round(gap_pct, 1),
                "n_pct": round(n_pct, 1),
                "status": status,
                "reason": ";".join(reasons),
            })

            if reasons:
                print(f"  REMOVE {seq_id}: {'; '.join(reasons)}", file=sys.stderr)
                n_removed += 1
            else:
                fout.write(f">{header}\n")
                # Write sequence in 60-char lines
                for i in range(0, len(seq), 60):
                    fout.write(seq[i:i+60] + "\n")
                n_kept += 1

    print(f"# Kept:    {n_kept}")
    print(f"# Removed: {n_removed}")
    print(f"# Output:  {args.output}")

    if args.report:
        with open(args.report, "w") as fh:
            fh.write("id\ttotal_len\tungapped_len\tgap_pct\tn_pct\tstatus\treason\n")
            for r in report_rows:
                fh.write(f"{r['id']}\t{r['total_len']}\t{r['ungapped_len']}\t"
                         f"{r['gap_pct']}\t{r['n_pct']}\t{r['status']}\t{r['reason']}\n")
        print(f"# Report:  {args.report}")


if __name__ == "__main__":
    main()
