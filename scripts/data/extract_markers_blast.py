#!/usr/bin/env python3
"""
extract_markers_blast.py — Extract plastid marker sequences from an assembled
plastome via BLAST, without requiring genome annotation.

For each marker, the script selects a representative reference sequence from
the Plan B genbank_dedup directory, BLASTs it against the assembled plastome,
and extracts the matching region. Handles both coding (matK, rbcL, rpoB, etc.)
and non-coding markers (trnL, psbA-trnH) with the same approach.

Usage:
    python extract_markers_blast.py \\
        --plastome  data/planA/plastomes/Amomum_biphyllum/*.fasta \\
        --ref_dir   data/planB/genbank_dedup \\
        --markers   matK,rbcL,trnL,psbA-trnH,rpoB,rpoC1,atpB,ndhF,ycf1,ycf2 \\
        --species   Amomum_biphyllum \\
        --run_id    SRR12824540 \\
        --output    data/planA/cds \\
        [--ref_plastome data/planA/reference/Zingiber_officinale_NC_037455.1.fasta]
        [--min_pct  70] \\
        [--min_cov  50] \\
        [--evalue   1e-10] \\
        [--threads  4]

Arguments:
    --plastome      Assembled plastome FASTA (single file or glob pattern)
    --ref_dir       Root directory of genbank_dedup (one subdir per marker)
    --markers       Comma-separated marker names
    --species       Binomial species name (Genus_species) for output filenames
    --run_id        SRA run ID for output filenames
    --output        Output root directory; marker subdirs created inside
    --ref_plastome  Optional: use Zingiber officinale (or another complete
                    plastome) as reference for all markers instead of
                    per-marker sequences from genbank_dedup
    --min_pct       Minimum % identity to accept a BLAST hit (default: 70)
    --min_cov       Minimum % of query covered by hit (default: 50)
    --evalue        BLAST e-value cutoff (default: 1e-10)
    --threads       Threads for BLAST (default: 4)

Output:
    <output>/<marker>/
        <species>_<run_id>_<marker>.fasta  — extracted sequence (one per marker)

    Markers with no BLAST hit are skipped (missing data — acceptable in
    partitioned ML with 20–40% missing data).

Tool version requirements:
    BLAST+ >= 2.12  (makeblastdb, blastn)
    Python >= 3.8   (standard library + Biopython)
    Biopython >= 1.79
"""

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils"))
from provenance_py import write_provenance

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("ERROR: Biopython required. Install: pip install biopython", file=sys.stderr)
    sys.exit(1)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="BLAST-based extraction of plastid markers from assembled plastomes."
    )
    p.add_argument("--plastome",     required=True,
                   help="Assembled plastome FASTA file (one sequence)")
    p.add_argument("--ref_dir",      required=True,
                   help="Root of genbank_dedup (one subdir per marker)")
    p.add_argument("--markers",      required=True,
                   help="Comma-separated marker names")
    p.add_argument("--species",      required=True,
                   help="Genus_species for output filenames")
    p.add_argument("--run_id",       required=True,
                   help="SRA run ID for output filenames")
    p.add_argument("--output",       required=True,
                   help="Output root directory")
    p.add_argument("--ref_plastome", default=None,
                   help="Optional reference plastome FASTA (use instead of per-marker refs)")
    p.add_argument("--min_pct",      type=float, default=70.0,
                   help="Minimum BLAST %% identity (default: 70)")
    p.add_argument("--min_cov",      type=float, default=50.0,
                   help="Minimum %% of query covered (default: 50)")
    p.add_argument("--evalue",       default="1e-10",
                   help="BLAST e-value cutoff (default: 1e-10)")
    p.add_argument("--threads",      type=int, default=4,
                   help="BLAST threads (default: 4)")
    p.add_argument("--min_assembly", type=int, default=50000,
                   help="Skip assemblies shorter than this (bp). Default: 50000")
    p.add_argument("--provenance",   default=None,
                   help="Directory for JSON provenance log")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Expected marker size bounds (bp) — used to exclude full-plastome contaminants
# download_genbank.sh occasionally stores the full plastome sequence in a
# marker directory when efetch returns the parent record. Any reference
# outside these bounds is excluded from BLAST reference selection.
# ---------------------------------------------------------------------------

MARKER_SIZE_BOUNDS = {
    "matK":      (700,   2_500),   # CDS ~1.5 kb; trnK intron amplicons up to ~2.4 kb
    "rbcL":      (400,   1_800),   # CDS ~1.4 kb
    "trnL":      (150,   2_500),   # intron + spacer; variable length
    "psbA-trnH": (100,   1_000),   # intergenic spacer
    "rpoB":      (600,   4_000),   # CDS ~3.2 kb; partial amplicons smaller
    "rpoC1":     (400,   2_500),   # CDS ~1.5 kb
    "atpB":      (400,   2_000),   # CDS ~1.5 kb
    "ndhF":      (500,   3_000),   # CDS ~2.3 kb
    "ycf1":      (500,  10_000),   # CDS varies; partial ycf1 commonly ~700 bp–5.5 kb
    "ycf2":      (500,  10_000),   # CDS ~6.8 kb; partial amplicons smaller
}

DEFAULT_BOUNDS = (200, 20_000)

# Per-marker minimum query-coverage overrides.
# ycf1 and ycf2 are very large genes (5–7 kb) that frequently assemble only
# partially from genome-skimming reads; even 300–600 bp of ycf1/ycf2 is
# phylogenetically informative. The global min_cov default (50%) would reject
# these short but valid hits, so we relax coverage for these two markers.
MARKER_MIN_COV = {
    "ycf1": 5.0,
    "ycf2": 5.0,
}


# ---------------------------------------------------------------------------
# Reference sequence selection
# ---------------------------------------------------------------------------

def pick_reference(marker, ref_dir):
    """
    Pick one reference sequence for BLAST from genbank_dedup/<marker>/.

    Filters out sequences outside MARKER_SIZE_BOUNDS (catches full-plastome
    records stored in marker directories by download_genbank.sh). Among
    qualifying candidates, prefers accessions without 'NC_' prefix (those
    are barcoding amplicons with clean gene boundaries) over NC_ plastome
    records, then picks the median-length candidate as the most representative.
    Returns the file path or None if no qualifying candidate found.
    """
    marker_dir = os.path.join(ref_dir, marker)
    if not os.path.isdir(marker_dir):
        return None

    fasta_files = [f for f in os.listdir(marker_dir) if f.endswith(".fasta")]
    if not fasta_files:
        return None

    lo, hi = MARKER_SIZE_BOUNDS.get(marker, DEFAULT_BOUNDS)

    # Score each candidate: measure length, filter by bounds
    candidates = []
    for fname in fasta_files:
        fpath = os.path.join(marker_dir, fname)
        try:
            rec = next(SeqIO.parse(fpath, "fasta"))
            seq_len = len(rec.seq)
            if lo <= seq_len <= hi:
                candidates.append((seq_len, fpath))
        except Exception:
            continue

    if not candidates:
        return None

    # Prefer non-NC_ (barcoding amplicons) for cleaner gene boundaries;
    # fall back to NC_ if only plastome-derived sequences pass the filter.
    non_nc = [(l, p) for l, p in candidates if "_NC_" not in os.path.basename(p)]
    pool = non_nc if non_nc else candidates

    # Pick the median-length sequence (avoids both shortest partial seqs
    # and any borderline long sequences)
    pool.sort(key=lambda x: x[0])
    _, best_file = pool[len(pool) // 2]
    return best_file


def extract_marker_from_ref_plastome(marker, ref_plastome_path):
    """
    Extract a marker sequence from the Zingiber reference plastome GenBank
    annotation. Handles CDS genes, tRNA genes (trnL-UAA), and the
    psbA-trnH intergenic spacer.

    Looks for a .gb file alongside the FASTA (same basename, .gb extension).
    Returns the sequence as a string, or None if extraction fails.
    """
    gb_path = ref_plastome_path.replace(".fasta", ".gb")
    if not os.path.exists(gb_path):
        return None

    try:
        rec = next(SeqIO.parse(gb_path, "genbank"))
    except Exception:
        return None

    # --- Marker-to-gene-name mapping (matches GenBank qualifier 'gene') ---
    gene_aliases = {
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

    # --- psbA-trnH: extract the intergenic spacer ---
    if marker == "psbA-trnH":
        psba_end   = None
        trnh_start = None
        for feat in rec.features:
            gene_name = feat.qualifiers.get("gene", [""])[0]
            if gene_name in ("psbA",) and feat.type in ("CDS", "gene"):
                loc = feat.location
                # Both psbA and trnH are on the minus strand in most plastomes;
                # collect all positions and take extremes.
                positions = [int(p) for p in loc]
                psba_end = max(positions) if psba_end is None else max(psba_end, max(positions))
            if gene_name in ("trnH-GUG", "trnH") and feat.type in ("tRNA", "gene"):
                loc = feat.location
                positions = [int(p) for p in loc]
                trnh_start = min(positions) if trnh_start is None else min(trnh_start, min(positions))

        if psba_end is not None and trnh_start is not None:
            lo = min(trnh_start, psba_end)
            hi = max(trnh_start, psba_end)
            spacer = rec.seq[lo:hi]
            if len(spacer) > 50:
                return str(spacer)
        return None

    # --- All other markers: find by gene name across CDS, gene, and tRNA features ---
    aliases = gene_aliases.get(marker, [marker])
    for feature in rec.features:
        gene_name = feature.qualifiers.get("gene", [""])[0]
        if any(a.lower() == gene_name.lower() for a in aliases):
            try:
                seq = feature.extract(rec.seq)
                if len(seq) > 50:
                    return str(seq)
            except Exception:
                continue
    return None


# ---------------------------------------------------------------------------
# BLAST execution and hit extraction
# ---------------------------------------------------------------------------

def blast_marker(query_fasta, subject_fasta, evalue, min_pct, min_cov, threads, tmpdir):
    """
    BLAST query_fasta against subject_fasta.
    Returns (sstart, send, sstrand, pct_identity, qcov) for best hit,
    or None if no qualifying hit found.
    """
    db_path = os.path.join(tmpdir, "subject_db")

    # Build BLAST database
    subprocess.run(
        ["makeblastdb", "-in", subject_fasta, "-dbtype", "nucl", "-out", db_path],
        check=True, capture_output=True
    )

    # Run BLAST: outfmt 6 with custom columns.
    # Notes:
    # - -max_hsps is intentionally omitted — it causes empty output for some
    #   markers (e.g. matK) even when a valid hit exists.
    # - -word_size 7 enables sensitive mode required for markers with ~70–75%
    #   identity across genera (e.g. matK). Default word_size 11 misses these.
    blast_out = os.path.join(tmpdir, "blast_hits.tsv")
    cols = "qseqid sseqid pident length qlen slen sstart send sstrand evalue"
    subprocess.run(
        [
            "blastn", "-query", query_fasta, "-db", db_path,
            "-outfmt", f"6 {cols}",
            "-evalue", evalue,
            "-num_threads", str(threads),
            "-max_target_seqs", "5",
            "-word_size", "7",
            "-out", blast_out,
        ],
        check=True, capture_output=True
    )

    if not os.path.exists(blast_out) or os.path.getsize(blast_out) == 0:
        return None

    best = None
    best_pct = 0.0
    with open(blast_out) as fh:
        for line in fh:
            row = line.strip().split("\t")
            if len(row) < 10:
                continue
            pct    = float(row[2])
            length = int(row[3])
            qlen   = int(row[4])
            sstart = int(row[6])
            send   = int(row[7])
            strand = row[8]   # "plus" or "minus"
            qcov   = 100.0 * length / qlen

            if pct < min_pct or qcov < min_cov:
                continue
            if pct > best_pct:
                best_pct = pct
                sseqid = row[1]   # subject sequence ID — needed for multi-scaffold assemblies
                best = (sseqid, sstart, send, strand, pct, qcov)

    return best


def extract_subsequence(fasta_path, sseqid, sstart, send, strand):
    """
    Extract a subsequence from fasta_path.  For multi-scaffold assemblies the
    BLAST hit may be on any contig, so we match by sseqid (the subject sequence
    ID reported by BLAST, which is the FASTA header word after '>').
    Falls back to the first record if sseqid is not found (single-sequence files).
    BLAST coordinates are 1-based inclusive.
    Reverse-complements if strand == 'minus'.
    """
    target_rec = None
    first_rec = None
    for rec in SeqIO.parse(fasta_path, "fasta"):
        if first_rec is None:
            first_rec = rec
        if rec.id == sseqid or rec.name == sseqid:
            target_rec = rec
            break
    rec = target_rec if target_rec is not None else first_rec
    if rec is None:
        return Seq("")
    seq = rec.seq

    lo = min(sstart, send) - 1   # convert to 0-based
    hi = max(sstart, send)
    subseq = seq[lo:hi]

    if strand == "minus":
        subseq = subseq.reverse_complement()

    return subseq


# ---------------------------------------------------------------------------
# Main extraction loop
# ---------------------------------------------------------------------------

def main():
    import time
    start = time.time()
    args  = parse_args()
    markers = [m.strip() for m in args.markers.split(",")]

    print(f"extract_markers_blast.py")
    print(f"  Plastome:  {args.plastome}")
    print(f"  Species:   {args.species}")
    print(f"  Run ID:    {args.run_id}")
    print(f"  Markers:   {', '.join(markers)}")
    print(f"  Output:    {args.output}")
    print()

    if not os.path.exists(args.plastome):
        print(f"ERROR: plastome FASTA not found: {args.plastome}", file=sys.stderr)
        sys.exit(1)

    # Assembly length check — skip fragmented assemblies that are too short
    # to reliably contain most markers.
    assembly_len = sum(len(r.seq) for r in SeqIO.parse(args.plastome, "fasta"))
    print(f"  Assembly length: {assembly_len:,} bp")
    if assembly_len < args.min_assembly:
        print(f"  Assembly too short ({assembly_len:,} bp < {args.min_assembly:,} bp threshold) — skipping all markers.")
        return 0

    n_extracted = 0
    n_failed    = 0
    results     = {}   # marker → "extracted" | "no_hit" | "no_ref"

    with tempfile.TemporaryDirectory(prefix="blast_extract_") as tmpdir:
        for marker in markers:
            print(f"  [{marker}]", end=" ", flush=True)

            # --- Pick reference sequence ---
            ref_fasta = None
            if args.ref_plastome:
                # Try to extract from reference plastome annotation first
                ref_seq = extract_marker_from_ref_plastome(marker, args.ref_plastome)
                if ref_seq:
                    ref_fasta = os.path.join(tmpdir, f"ref_{marker}.fasta")
                    with open(ref_fasta, "w") as fh:
                        fh.write(f">{marker}_reference\n{ref_seq}\n")

            if ref_fasta is None:
                ref_file = pick_reference(marker, args.ref_dir)
                if ref_file is None:
                    print("no reference → skip")
                    results[marker] = "no_ref"
                    n_failed += 1
                    continue
                ref_fasta = ref_file

            # --- BLAST ---
            effective_min_cov = MARKER_MIN_COV.get(marker, args.min_cov)
            try:
                hit = blast_marker(
                    ref_fasta, args.plastome,
                    args.evalue, args.min_pct, effective_min_cov, args.threads, tmpdir
                )
            except subprocess.CalledProcessError as e:
                print(f"BLAST error → skip ({e})")
                results[marker] = "blast_error"
                n_failed += 1
                continue

            if hit is None:
                print(f"no hit (pct≥{args.min_pct}%% cov≥{effective_min_cov}%%) → skip")
                results[marker] = "no_hit"
                n_failed += 1
                continue

            sseqid, sstart, send, strand, pct, qcov = hit
            hit_len = abs(send - sstart) + 1

            # --- Extract subsequence ---
            subseq = extract_subsequence(args.plastome, sseqid, sstart, send, strand)

            # --- Write output ---
            out_marker_dir = os.path.join(args.output, marker)
            os.makedirs(out_marker_dir, exist_ok=True)
            out_fname = f"{args.species}_{args.run_id}_{marker}.fasta"
            out_path  = os.path.join(out_marker_dir, out_fname)

            header = f">{args.species}_{args.run_id}_{marker} extracted_by=BLAST pct={pct:.1f} cov={qcov:.1f} strand={strand} len={hit_len} scaffold={sseqid}"
            with open(out_path, "w") as fh:
                fh.write(f"{header}\n{subseq}\n")

            print(f"OK ({hit_len} bp, {pct:.1f}%% id, {qcov:.1f}%% qcov, {strand})")
            results[marker] = "extracted"
            n_extracted += 1

    print()
    print(f"Summary: {n_extracted} markers extracted, {n_failed} skipped")
    print(f"Results: {results}")

    if args.provenance:
        write_provenance(
            script="extract_markers_blast",
            tool="blastn",
            version="2.16.0+",
            provenance_dir=args.provenance,
            parameters={
                "plastome":   args.plastome,
                "species":    args.species,
                "run_id":     args.run_id,
                "markers":    markers,
                "min_pct":    args.min_pct,
                "min_cov":    args.min_cov,
                "evalue":     args.evalue,
                "n_extracted": n_extracted,
                "n_failed":   n_failed,
                "results":    results,
            },
            input_files={"plastome": args.plastome, "ref_dir": args.ref_dir},
            output_files={"output": args.output},
            runtime_seconds=int(time.time() - start),
            exit_code=0,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
