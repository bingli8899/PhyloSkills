#!/usr/bin/env python3
"""
revise_hybpiper_sequences.py — Post-process HybPiper retrieved sequences

HybPiper v2 retrieve_sequences outputs one FASTA per locus with headers in the
format: ">SampleName". This script:
  1. Renames headers to Genus_species_samplename format using a metadata table
  2. Removes sequences shorter than a minimum length threshold
  3. Optionally removes samples with too many missing loci across all locus files
  4. Reports per-locus and per-sample recovery statistics

Usage:
    python revise_hybpiper_sequences.py \
        --input <retrieved_sequences_dir> \
        --output <outdir> \
        --metadata <metadata.tsv> \
        [--min_length <int>] \
        [--max_missing_pct <float>] \
        [--id_col <col_name>] \
        [--genus_col <col_name>] \
        [--species_col <col_name>]

Arguments:
    --input          Directory with per-locus FASTA files from hybpiper retrieve_sequences
    --output         Output directory for revised FASTA files
    --metadata       TSV file mapping sample IDs to taxonomic names
                     Must contain columns for sample ID, genus, and species
    --min_length     Minimum sequence length in bp to retain (default: 100)
    --max_missing_pct  Exclude samples missing more than this % of loci (default: 80.0)
    --id_col         Column name for sample ID in metadata (default: sample_id)
    --genus_col      Column name for genus (default: genus)
    --species_col    Column name for species (default: species)

Tool version requirements:
    Python >= 3.8  (no external dependencies)
"""

import argparse
import sys
import os
from pathlib import Path
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Post-process HybPiper retrieved sequences"
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--min_length", type=int, default=100)
    parser.add_argument("--max_missing_pct", type=float, default=80.0)
    parser.add_argument("--id_col", default="sample_id")
    parser.add_argument("--genus_col", default="genus")
    parser.add_argument("--species_col", default="species")
    return parser.parse_args()


def load_metadata(meta_file: str, id_col: str, genus_col: str,
                  species_col: str) -> dict[str, str]:
    """Returns {sample_id: 'Genus_species'} mapping."""
    mapping = {}
    with open(meta_file) as fh:
        header = None
        for line in fh:
            line = line.rstrip()
            if line.startswith("#"):
                continue
            if header is None:
                header = line.split("\t")
                # Validate required columns
                for col in (id_col, genus_col, species_col):
                    if col not in header:
                        print(f"ERROR: Column '{col}' not found in metadata file.",
                              file=sys.stderr)
                        print(f"  Available columns: {header}", file=sys.stderr)
                        sys.exit(1)
                continue
            parts = line.split("\t")
            row = dict(zip(header, parts))
            sample_id = row.get(id_col, "").strip()
            genus = row.get(genus_col, "Unknown").strip()
            species = row.get(species_col, "sp").strip()
            if sample_id:
                taxon = f"{genus}_{species}"
                # Sanitize
                taxon = "".join(c if c.isalnum() or c == "_" else "_" for c in taxon)
                mapping[sample_id] = taxon
    return mapping


def parse_fasta(path: str):
    header, seq_parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].split()[0]
                seq_parts = []
            elif line:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


def main():
    args = parse_args()
    print(f"# revise_hybpiper_sequences.py")
    print(f"# Input dir:     {args.input}")
    print(f"# Metadata:      {args.metadata}")
    print(f"# min_length:    {args.min_length} bp")
    print(f"# max_missing:   {args.max_missing_pct}%")
    print()

    metadata = load_metadata(
        args.metadata, args.id_col, args.genus_col, args.species_col
    )
    print(f"Loaded {len(metadata)} sample-to-taxon mappings")

    indir = Path(args.input)
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    fasta_files = sorted(
        list(indir.glob("*.fasta")) + list(indir.glob("*.fa"))
    )
    if not fasta_files:
        print(f"ERROR: No FASTA files found in {indir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(fasta_files)} locus FASTA files")
    print()

    # First pass: collect all sample IDs across all loci
    # locus → {taxon: seq}
    locus_data: dict[str, dict[str, str]] = {}
    all_samples: set[str] = set()

    for fasta in fasta_files:
        locus = fasta.stem
        locus_data[locus] = {}
        for sample_id, seq in parse_fasta(str(fasta)):
            all_samples.add(sample_id)
            ungapped_len = len(seq.replace("-", "").replace("N", ""))
            if ungapped_len < args.min_length:
                print(f"  REMOVE {sample_id} from {locus}: "
                      f"length {ungapped_len} bp < {args.min_length}", file=sys.stderr)
                continue
            # Rename using metadata
            if sample_id in metadata:
                taxon = f"{metadata[sample_id]}_{sample_id}"
            else:
                print(f"  WARNING: No metadata for '{sample_id}' — using original ID",
                      file=sys.stderr)
                taxon = sample_id
            locus_data[locus][taxon] = seq

    n_loci = len(locus_data)

    # Second pass: identify samples with excessive missing data
    sample_locus_count: dict[str, int] = defaultdict(int)
    for locus, seqs in locus_data.items():
        for taxon in seqs:
            sample_locus_count[taxon] += 1

    excluded_samples = set()
    for taxon, count in sample_locus_count.items():
        missing_pct = (1 - count / n_loci) * 100
        if missing_pct > args.max_missing_pct:
            excluded_samples.add(taxon)
            print(f"  EXCLUDE {taxon}: {missing_pct:.0f}% loci missing "
                  f"({count}/{n_loci} recovered)")

    print(f"\nExcluded samples (>{args.max_missing_pct}% missing): "
          f"{len(excluded_samples)}")
    print()

    # Write revised FASTA files
    n_written = 0
    for locus, seqs in sorted(locus_data.items()):
        out_seqs = {t: s for t, s in seqs.items() if t not in excluded_samples}
        if not out_seqs:
            print(f"  SKIP {locus}: no sequences after filtering")
            continue

        out_path = outdir / f"{locus}.fasta"
        with open(out_path, "w") as fout:
            for taxon, seq in sorted(out_seqs.items()):
                fout.write(f">{taxon}\n")
                for i in range(0, len(seq), 60):
                    fout.write(seq[i:i+60] + "\n")
        print(f"  {locus}: {len(out_seqs)} sequences → {out_path.name}")
        n_written += 1

    print(f"\nWrote {n_written} locus files to {outdir}/")
    print("Next step: phylo-alignment — run align_markers.sh on output directory")


if __name__ == "__main__":
    main()
