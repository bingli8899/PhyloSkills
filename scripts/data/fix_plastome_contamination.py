#!/usr/bin/env python3
"""
fix_plastome_contamination.py — For each complete-plastome FASTA that was
incorrectly stored in a marker directory, download the GenBank annotation,
extract the correct marker region, and replace the contaminated file.

The problem: download_genbank.sh searches for markers with
  "marker"[All Fields]
which also matches complete plastome records. efetch -format fasta returns the
entire plastome (100-170 kb) instead of just the ~1-2 kb gene.

This script fixes the already-downloaded contaminated files by:
  1. Identifying all FASTA files in genbank_dedup/<marker>/ that are over-length
  2. Downloading the GenBank annotation (.gb) for each unique accession
  3. Extracting the gene/region from the annotation
  4. Replacing the contaminated FASTA with the correctly extracted sequence
  5. Deleting the file if extraction fails (missing data is preferable to wrong data)

Usage:
    python fix_plastome_contamination.py \\
        --genbank_dedup  data/planB/genbank_dedup \\
        --markers        matK,rbcL,trnL,psbA-trnH,rpoB,rpoC1,atpB,ndhF,ycf1,ycf2 \\
        [--dry_run]      # preview without modifying files

Tool requirements:
    Biopython >= 1.79
    efetch (Entrez Direct)
    Python >= 3.8
"""

import argparse
import os
import subprocess
import sys
import time
import tempfile
from collections import defaultdict
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("ERROR: Biopython required. pip install biopython", file=sys.stderr)
    sys.exit(1)


# ---------------------------------------------------------------------------
# Marker size bounds — sequences larger than max are whole-genome contaminants
# ---------------------------------------------------------------------------
# Any sequence > 20 kb in a marker directory is a complete plastome contaminant.
# We do NOT use tighter per-marker thresholds here because sequences slightly
# over the CDS length are often valid amplicons:
#   - matK amplicons often include the trnK intron (~2500-2800 bp)
#   - psbA-trnH spacer can reach ~1200 bp in Zingiberaceae
# Only whole plastomes (100-170 kb) need to be replaced.
CONTAMINATION_THRESHOLD = 20_000   # bp

# Keep MARKER_SIZE_MAX only for context / documentation; not used for detection.
MARKER_SIZE_MAX = {
    "matK":      2_500,
    "rbcL":      1_800,
    "trnL":      2_500,
    "psbA-trnH": 1_000,
    "rpoB":      4_000,
    "rpoC1":     2_500,
    "atpB":      2_000,
    "ndhF":      3_000,
    "ycf1":     10_000,
    "ycf2":     10_000,
}
DEFAULT_MAX = CONTAMINATION_THRESHOLD

# Gene name aliases for feature lookup in GenBank annotations
GENE_ALIASES = {
    "matK":      ["matK", "matk"],
    "rbcL":      ["rbcL", "rbcl", "rbcL1"],
    "rpoB":      ["rpoB", "rpob"],
    "rpoC1":     ["rpoC1", "rpoc1"],
    "atpB":      ["atpB", "atpb"],
    "ndhF":      ["ndhF", "ndhf"],
    "ycf1":      ["ycf1"],
    "ycf2":      ["ycf2"],
    "trnL":      ["trnL-UAA", "trnL", "trnl-uaa", "trnl"],
}


def parse_args():
    p = argparse.ArgumentParser(
        description="Replace complete-plastome contamination in marker directories."
    )
    p.add_argument("--genbank_dedup", required=True,
                   help="Root of genbank_dedup (one subdir per marker)")
    p.add_argument("--markers", required=True,
                   help="Comma-separated marker names")
    p.add_argument("--dry_run", action="store_true",
                   help="Preview actions without modifying files")
    return p.parse_args()


def find_contaminated(genbank_dedup, markers):
    """
    Returns dict: accession → {marker: fasta_path, ...}
    for every FASTA file that exceeds the marker's expected size maximum.
    """
    by_accession = defaultdict(dict)

    for marker in markers:
        marker_dir = os.path.join(genbank_dedup, marker)
        if not os.path.isdir(marker_dir):
            continue
        for fasta_path in sorted(Path(marker_dir).glob("*.fasta")):
            try:
                rec = next(SeqIO.parse(str(fasta_path), "fasta"))
            except Exception:
                continue
            seq_len = len(rec.seq)
            if seq_len > CONTAMINATION_THRESHOLD:
                # Get accession from FASTA header (first word after >)
                acc = rec.id.split()[0]
                by_accession[acc][marker] = str(fasta_path)

    return by_accession


def download_gb(accession, tmpdir):
    """
    Download GenBank flat file for accession using efetch.
    Returns path to .gb file or None on failure.
    """
    gb_path = os.path.join(tmpdir, f"{accession.replace('.', '_')}.gb")
    if os.path.exists(gb_path) and os.path.getsize(gb_path) > 0:
        return gb_path   # already downloaded

    try:
        result = subprocess.run(
            ["efetch", "-db", "nuccore", "-id", accession, "-format", "gb"],
            capture_output=True, text=True, timeout=120
        )
        if result.returncode != 0 or not result.stdout.strip():
            return None
        with open(gb_path, "w") as fh:
            fh.write(result.stdout)
        return gb_path
    except Exception as e:
        print(f"    efetch error for {accession}: {e}", file=sys.stderr)
        return None


def extract_marker_from_gb(gb_path, marker, accession):
    """
    Extract the marker sequence from a GenBank annotation file.
    Returns (sequence_string, feature_type) or (None, None).
    """
    try:
        rec = next(SeqIO.parse(gb_path, "genbank"))
    except Exception as e:
        print(f"    Cannot parse {gb_path}: {e}", file=sys.stderr)
        return None, None

    # psbA-trnH: intergenic spacer between psbA and trnH-GUG
    if marker == "psbA-trnH":
        psba_end = None
        trnh_start = None
        for feat in rec.features:
            gene_name = feat.qualifiers.get("gene", [""])[0]
            if gene_name in ("psbA",) and feat.type in ("CDS", "gene"):
                positions = sorted(int(p) for p in feat.location)
                end = max(positions)
                psba_end = end if psba_end is None else max(psba_end, end)
            if gene_name in ("trnH-GUG", "trnH") and feat.type in ("tRNA", "gene"):
                positions = sorted(int(p) for p in feat.location)
                start = min(positions)
                trnh_start = start if trnh_start is None else min(trnh_start, start)
        if psba_end is not None and trnh_start is not None:
            lo = min(trnh_start, psba_end)
            hi = max(trnh_start, psba_end)
            spacer = rec.seq[lo:hi]
            if len(spacer) > 50:
                return str(spacer), "intergenic_spacer"
        return None, None

    # All other markers: find gene by name
    aliases = GENE_ALIASES.get(marker, [marker])
    for feat in rec.features:
        gene_name = feat.qualifiers.get("gene", [""])[0]
        if any(a.lower() == gene_name.lower() for a in aliases):
            if feat.type in ("CDS", "gene", "tRNA"):
                try:
                    seq = feat.extract(rec.seq)
                    if len(seq) > 50:
                        return str(seq), feat.type
                except Exception:
                    continue

    return None, None


def build_fasta_header(fasta_path, marker, acc):
    """
    Reconstruct the expected FASTA header from the filename.
    File naming convention: Genus_species_ACCESSION_marker.fasta
    """
    base = os.path.basename(fasta_path)
    # Strip the trailing _marker.fasta to get Genus_species_ACCESSION
    prefix = base[: -(len(marker) + len(".fasta") + 1)]  # remove _marker.fasta
    return f">{prefix}_{marker} source=GenBank_annotation"


def main():
    args = parse_args()
    markers = [m.strip() for m in args.markers.split(",")]

    print("fix_plastome_contamination.py")
    print(f"  genbank_dedup: {args.genbank_dedup}")
    print(f"  markers:       {', '.join(markers)}")
    if args.dry_run:
        print("  [DRY RUN] — no files will be modified")
    print()

    # Step 1: Find all contaminated files
    print("Scanning for over-length sequences...")
    by_acc = find_contaminated(args.genbank_dedup, markers)
    total_files = sum(len(v) for v in by_acc.values())
    print(f"  {len(by_acc)} unique accessions; {total_files} files to fix\n")

    if not by_acc:
        print("Nothing to fix.")
        return 0

    n_fixed = 0
    n_deleted = 0
    n_failed = 0
    skipped_acc = []

    with tempfile.TemporaryDirectory(prefix="fix_plastome_") as tmpdir:
        for i, (acc, marker_files) in enumerate(sorted(by_acc.items()), 1):
            print(f"[{i}/{len(by_acc)}] {acc}  ({len(marker_files)} markers)")

            if args.dry_run:
                for m, path in sorted(marker_files.items()):
                    print(f"    would fix: {os.path.basename(path)}")
                continue

            # Download GenBank annotation (once per accession)
            gb_path = download_gb(acc, tmpdir)
            if gb_path is None:
                print(f"    FAILED to download GenBank for {acc} — deleting contaminated files")
                for m, path in marker_files.items():
                    os.remove(path)
                    n_deleted += 1
                skipped_acc.append(acc)
                # Rate limit
                time.sleep(0.4)
                continue

            # Extract each marker from the annotation
            for marker, fasta_path in sorted(marker_files.items()):
                seq, feat_type = extract_marker_from_gb(gb_path, marker, acc)

                if seq is None:
                    print(f"    {marker:<14} NOT FOUND in annotation — deleting")
                    os.remove(fasta_path)
                    n_deleted += 1
                else:
                    header = build_fasta_header(fasta_path, marker, acc)
                    with open(fasta_path, "w") as fh:
                        fh.write(f"{header} len={len(seq)} feat={feat_type}\n{seq}\n")
                    print(f"    {marker:<14} OK  {len(seq)} bp  ({feat_type})")
                    n_fixed += 1

            # Rate limit: NCBI allows ~3 req/sec without API key, ~10 with
            time.sleep(0.4)

    print()
    print(f"Fixed:   {n_fixed}")
    print(f"Deleted: {n_deleted}  (gene not in annotation or download failed)")
    if skipped_acc:
        print(f"Failed downloads ({len(skipped_acc)}): {', '.join(skipped_acc[:5])}" +
              (" ..." if len(skipped_acc) > 5 else ""))
    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
