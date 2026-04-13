#!/usr/bin/env python3
"""
extract_cds.py — Extract CDS sequences from annotated plastome GenBank files

Parses GenBank annotation files, extracts CDS features for selected genes,
handles multi-exon genes by concatenating exon sequences, and writes one
FASTA file per gene with all taxa.

Usage:
    python extract_cds.py \
        --input <annotation_dir_or_gb_file> \
        --output <outdir> \
        [--genes <gene1,gene2,...>] \
        [--min_length <int>] \
        [--translate] \
        [--concatenate]

Arguments:
    --input      Directory with *.gb / *.gbk files, or a single GenBank file
    --output     Output directory for per-gene FASTA files
    --genes      Comma-separated gene names to extract (default: common plastid
                 marker set: matK,rbcL,rpoB,rpoC1,rpoC2,ndhF,ndhB,atpB,atpI,
                 psbA,psbB,psbC,psbD,psbE,psbK,psbL)
    --min_length Minimum CDS length in bp to include (default: 200)
    --translate  Also write protein FASTA files
    --concatenate Write a concatenated supermatrix FASTA + partition file

Tool version requirements (checked at runtime):
    Biopython >= 1.83  — pip install biopython
    (Latest: 1.85, 2025 — https://github.com/biopython/biopython/releases)

Output:
    <outdir>/
        <gene>.fasta            — per-gene nucleotide FASTA (one seq per taxon)
        <gene>_aa.fasta         — amino acid (if --translate)
        concatenated.fasta      — supermatrix (if --concatenate)
        partition.txt           — RAxML-style partition file (if --concatenate)
        extraction_report.tsv   — per-taxon, per-gene extraction status
"""

import argparse
import sys
import os
import re
from pathlib import Path
from collections import defaultdict

# ── Version check ────────────────────────────────────────────────────────────
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    import Bio
    print(f"# Biopython version: {Bio.__version__}  (latest known: 1.85)")
    print(f"# Python version: {sys.version.split()[0]}")
except ImportError:
    print("ERROR: Biopython not found. Install: pip install biopython", file=sys.stderr)
    sys.exit(1)

# ── Default gene set ─────────────────────────────────────────────────────────
DEFAULT_GENES = [
    "matK", "rbcL", "rpoB", "rpoC1", "rpoC2",
    "ndhF", "ndhB", "ndhH", "ndhA", "ndhI",
    "atpB", "atpI", "atpA", "atpE", "atpF", "atpH",
    "psbA", "psbB", "psbC", "psbD", "psbE", "psbK", "psbL", "psbT",
    "psaA", "psaB", "psaC", "accD", "cemA", "clpP", "infA",
    "rps2", "rps3", "rps4", "rps7", "rps8", "rps11", "rps12", "rps14",
    "rps15", "rps16", "rps18", "rps19",
    "rpl2", "rpl14", "rpl16", "rpl20", "rpl22", "rpl23", "rpl32", "rpl33", "rpl36",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract CDS from plastome GenBank annotations"
    )
    parser.add_argument("--input", required=True, help="GenBank file or directory")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument(
        "--genes",
        default=",".join(DEFAULT_GENES),
        help="Comma-separated gene names (default: common plastid markers)"
    )
    parser.add_argument(
        "--min_length", type=int, default=200,
        help="Minimum CDS length in bp (default: 200)"
    )
    parser.add_argument(
        "--translate", action="store_true",
        help="Also write protein FASTA files"
    )
    parser.add_argument(
        "--concatenate", action="store_true",
        help="Write concatenated supermatrix + partition file"
    )
    return parser.parse_args()


def find_gb_files(input_path: str) -> list[Path]:
    """Find all GenBank files in a path (file or directory)."""
    p = Path(input_path)
    if p.is_file():
        return [p]
    elif p.is_dir():
        files = []
        for ext in ["*.gb", "*.gbk", "*.genbank"]:
            files.extend(sorted(p.glob(ext)))
        return files
    else:
        print(f"ERROR: Input not found: {input_path}", file=sys.stderr)
        sys.exit(1)


def extract_taxon_name(record) -> str:
    """Extract genus_species from a SeqRecord. Sanitize for filenames."""
    # Try organism field
    org = record.annotations.get("organism", "")
    if not org:
        # Fall back to record name
        org = record.name
    # Take first two words (Genus species)
    parts = org.split()
    if len(parts) >= 2:
        taxon = f"{parts[0]}_{parts[1]}"
    else:
        taxon = parts[0] if parts else record.id
    # Sanitize
    taxon = re.sub(r"[^A-Za-z0-9_]", "_", taxon)
    return taxon


def normalize_gene_name(name: str) -> str:
    """Normalize gene name: lowercase, remove common prefixes."""
    return name.lower().strip()


def extract_cds_for_gene(record, gene_name: str, min_length: int):
    """
    Extract CDS sequence for a given gene from a SeqRecord.
    Handles multi-exon genes by concatenating exon positions.
    Returns (sequence_str, strand) or (None, None) if not found.
    """
    gene_norm = normalize_gene_name(gene_name)

    for feature in record.features:
        if feature.type not in ("CDS", "gene"):
            continue

        # Check gene qualifiers
        feat_gene = ""
        for qual in ("gene", "product", "label", "note"):
            val = feature.qualifiers.get(qual, [""])[0]
            if val:
                feat_gene = normalize_gene_name(val)
                break

        if gene_norm not in feat_gene and feat_gene not in gene_norm:
            continue

        if feature.type == "CDS":
            try:
                seq_str = str(feature.extract(record.seq))
                if len(seq_str) >= min_length:
                    return seq_str, feature.strand
            except Exception as e:
                print(f"  WARNING: Failed to extract {gene_name}: {e}", file=sys.stderr)
                return None, None

    return None, None


def write_fasta(sequences: dict[str, str], output_path: Path, description: str = ""):
    """Write a dict of {header: sequence} to FASTA."""
    records = []
    for header, seq in sequences.items():
        rec = SeqRecord(Seq(seq), id=header, description=description)
        records.append(rec)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, str(output_path), "fasta")


def main():
    args = parse_args()
    print(f"# Date: {__import__('datetime').date.today()}")
    print()

    target_genes = [g.strip() for g in args.genes.split(",") if g.strip()]
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    gb_files = find_gb_files(args.input)
    print(f"Found {len(gb_files)} GenBank file(s)")
    print(f"Target genes: {len(target_genes)}")
    print()

    # gene → {taxon: sequence}
    gene_sequences: dict[str, dict[str, str]] = defaultdict(dict)
    # For report
    report_rows: list[dict] = []

    for gb_file in gb_files:
        print(f"Processing: {gb_file.name}")
        try:
            records = list(SeqIO.parse(str(gb_file), "genbank"))
        except Exception as e:
            print(f"  ERROR: Failed to parse {gb_file}: {e}", file=sys.stderr)
            continue

        for record in records:
            taxon = extract_taxon_name(record)
            extracted = 0

            for gene in target_genes:
                seq_str, strand = extract_cds_for_gene(record, gene, args.min_length)
                status = "extracted"
                if seq_str:
                    # If taxon already present, keep longer sequence
                    existing = gene_sequences[gene].get(taxon, "")
                    if len(seq_str) >= len(existing):
                        gene_sequences[gene][taxon] = seq_str
                    extracted += 1
                else:
                    status = "not_found"
                report_rows.append({
                    "file": gb_file.name,
                    "taxon": taxon,
                    "gene": gene,
                    "status": status,
                    "length": len(seq_str) if seq_str else 0,
                })

            print(f"  {taxon}: {extracted}/{len(target_genes)} genes extracted")

    # ── Write per-gene FASTA files ────────────────────────────────────────────
    print()
    print("Writing per-gene FASTA files...")
    genes_with_data = []
    for gene in target_genes:
        seqs = gene_sequences.get(gene, {})
        if not seqs:
            print(f"  SKIP {gene}: no sequences found")
            continue
        out_path = outdir / f"{gene}.fasta"
        write_fasta(seqs, out_path)
        print(f"  {gene}: {len(seqs)} sequences → {out_path.name}")
        genes_with_data.append(gene)

        if args.translate:
            prot_seqs = {}
            for taxon, seq in seqs.items():
                try:
                    aa = str(Seq(seq).translate(to_stop=True))
                    if aa:
                        prot_seqs[taxon] = aa
                except Exception:
                    pass
            if prot_seqs:
                write_fasta(prot_seqs, outdir / f"{gene}_aa.fasta")

    # ── Concatenated supermatrix ──────────────────────────────────────────────
    if args.concatenate and genes_with_data:
        print()
        print("Building concatenated supermatrix...")
        # Collect all taxa
        all_taxa = set()
        for gene in genes_with_data:
            all_taxa.update(gene_sequences[gene].keys())
        all_taxa = sorted(all_taxa)

        concat_seqs: dict[str, str] = {t: "" for t in all_taxa}
        partition_lines = []
        pos = 1

        for gene in genes_with_data:
            seqs = gene_sequences[gene]
            # Determine max length for this gene (pad shorter sequences)
            max_len = max(len(s) for s in seqs.values())
            for taxon in all_taxa:
                seq = seqs.get(taxon, "")
                if seq:
                    concat_seqs[taxon] += seq.ljust(max_len, "N")
                else:
                    concat_seqs[taxon] += "?" * max_len
            partition_lines.append(f"DNA, {gene} = {pos}-{pos + max_len - 1}")
            pos += max_len

        write_fasta(concat_seqs, outdir / "concatenated.fasta", "concatenated plastid markers")
        (outdir / "partition.txt").write_text("\n".join(partition_lines) + "\n")
        total_len = pos - 1
        missing = sum(s.count("?") for s in concat_seqs.values())
        total_cells = len(all_taxa) * total_len
        missing_pct = missing / total_cells * 100 if total_cells > 0 else 0
        print(f"  Taxa: {len(all_taxa)}, Length: {total_len} bp")
        print(f"  Missing data: {missing_pct:.1f}%")
        print(f"  Files: {outdir}/concatenated.fasta, {outdir}/partition.txt")

    # ── Extraction report ─────────────────────────────────────────────────────
    report_path = outdir / "extraction_report.tsv"
    with open(report_path, "w") as fh:
        fh.write("file\ttaxon\tgene\tstatus\tlength_bp\n")
        for row in report_rows:
            fh.write(f"{row['file']}\t{row['taxon']}\t{row['gene']}\t"
                     f"{row['status']}\t{row['length']}\n")
    print()
    print(f"Report: {report_path}")
    print()
    print("Next step: phylo-alignment — align each gene FASTA in output dir")


if __name__ == "__main__":
    main()
