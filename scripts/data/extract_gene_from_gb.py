#!/usr/bin/env python3
# =============================================================================
# extract_gene_from_gb.py — Extract a single marker gene from a GenBank file
# =============================================================================
# Usage:
#   python3 extract_gene_from_gb.py \
#       --gb <file.gb> --gene <marker> --out <output.fasta> --header <">header">
#
# Arguments:
#   --gb      Path to GenBank (.gb) file
#   --gene    Marker name: matK, rbcL, rpoB, rpoC1, atpB, ndhF, ycf1, ycf2,
#             trnL, or psbA-trnH
#   --out     Output FASTA path (will be overwritten if it exists)
#   --header  FASTA header line (without leading ">"; the ">" is added here)
#
# Exit codes:
#   0 — success (marker extracted and written)
#   1 — failure (gene not found, file missing, or other error)
#
# Intended to be called from download_genbank.sh after a size-check reveals
# that efetch returned a complete plastome rather than a gene-length record.
# =============================================================================

import argparse
import sys

# ---------------------------------------------------------------------------
# Gene alias table — maps caller-facing marker names to GenBank /gene qualifier
# values used in CDS, gene, and tRNA features.
# ---------------------------------------------------------------------------
GENE_ALIASES = {
    "matK":  ["matK", "matk"],
    "rbcL":  ["rbcL", "rbcl", "rbcL1"],
    "rpoB":  ["rpoB", "rpob"],
    "rpoC1": ["rpoC1", "rpoc1"],
    "atpB":  ["atpB", "atpb"],
    "ndhF":  ["ndhF", "ndhf"],
    "ycf1":  ["ycf1"],
    "ycf2":  ["ycf2"],
    "trnL":  ["trnL-UAA", "trnL"],   # tRNA gene with intron
}

# Feature types to search (in priority order — checked in iteration order)
FEATURE_TYPES = {"CDS", "gene", "tRNA", "misc_RNA"}


def extract_psba_trnh_spacer(rec):
    """
    Extract the psbA–trnH intergenic spacer.
    Finds the end of psbA and start of trnH-GUG, returns the sequence between.
    Both genes are on the minus strand in most plastomes; we take coordinate
    extremes rather than strand-aware positions.
    Returns str or None.
    """
    psba_end   = None
    trnh_start = None

    for feat in rec.features:
        gene_name = feat.qualifiers.get("gene", [""])[0]
        if gene_name in ("psbA",) and feat.type in ("CDS", "gene"):
            positions = [int(p) for p in feat.location]
            psba_end = (max(positions)
                        if psba_end is None
                        else max(psba_end, max(positions)))
        if gene_name in ("trnH-GUG", "trnH") and feat.type in ("tRNA", "gene"):
            positions = [int(p) for p in feat.location]
            trnh_start = (min(positions)
                          if trnh_start is None
                          else min(trnh_start, min(positions)))

    if psba_end is not None and trnh_start is not None:
        lo = min(trnh_start, psba_end)
        hi = max(trnh_start, psba_end)
        spacer = rec.seq[lo:hi]
        if len(spacer) > 50:
            return str(spacer)
    return None


def extract_gene(rec, marker):
    """
    Extract marker sequence from a SeqRecord.
    Returns str or None.
    """
    if marker == "psbA-trnH":
        return extract_psba_trnh_spacer(rec)

    aliases = GENE_ALIASES.get(marker, [marker])

    for feature in rec.features:
        if feature.type not in FEATURE_TYPES:
            continue
        gene_name = feature.qualifiers.get("gene", [""])[0]
        if any(a.lower() == gene_name.lower() for a in aliases):
            try:
                seq = feature.extract(rec.seq)
                if len(seq) > 50:
                    return str(seq)
            except Exception:
                continue
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Extract a single marker gene from a GenBank annotation file."
    )
    parser.add_argument("--gb",     required=True, help="Path to GenBank file")
    parser.add_argument("--gene",   required=True, help="Marker/gene name (e.g. matK, rbcL)")
    parser.add_argument("--out",    required=True, help="Output FASTA path")
    parser.add_argument("--header", required=True,
                        help="FASTA header (without '>'; the '>' will be prepended)")
    args = parser.parse_args()

    # Import Biopython — give a clear error if absent
    try:
        from Bio import SeqIO
    except ImportError:
        print("ERROR: Biopython is not installed. "
              "Run: pip install biopython", file=sys.stderr)
        sys.exit(1)

    # Parse the GenBank file
    try:
        records = list(SeqIO.parse(args.gb, "genbank"))
    except Exception as exc:
        print(f"ERROR: Could not parse GenBank file '{args.gb}': {exc}", file=sys.stderr)
        sys.exit(1)

    if not records:
        print(f"ERROR: No records found in '{args.gb}'", file=sys.stderr)
        sys.exit(1)

    # Try each record (most GB files have exactly one, but plastome files
    # occasionally have two records when they are split at the IR junction)
    seq_str = None
    for rec in records:
        seq_str = extract_gene(rec, args.gene)
        if seq_str:
            break

    if not seq_str:
        print(f"ERROR: Could not find '{args.gene}' feature in '{args.gb}'",
              file=sys.stderr)
        sys.exit(1)

    # Write FASTA output
    header = args.header.lstrip(">")   # defensive — strip leading > if caller included it
    try:
        with open(args.out, "w") as fh:
            fh.write(f">{header}\n")
            # Wrap sequence at 60 characters (standard FASTA)
            for i in range(0, len(seq_str), 60):
                fh.write(seq_str[i:i + 60] + "\n")
    except OSError as exc:
        print(f"ERROR: Could not write output file '{args.out}': {exc}", file=sys.stderr)
        sys.exit(1)

    print(f"Extracted {len(seq_str)} bp for '{args.gene}' → {args.out}")
    sys.exit(0)


if __name__ == "__main__":
    main()
