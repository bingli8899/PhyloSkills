#!/usr/bin/env bash
# =============================================================================
# assemble_and_extract.sh — Per-sample wrapper: assemble plastome from SRA
# reads, extract 10 Plan B markers, merge into genbank_dedup.
#
# Called by download_sra.sh in streaming mode:
#   bash assemble_and_extract.sh <acc> <R1.fastq> <R2.fastq>
#
# Required environment variables (set by caller or exported before running):
#   PHYLOSKILLS_ROOT   — root of the PhyloSkills repo
#   PROJECT_ROOT       — root of the Zingiberaceae project (where data/ lives)
#   SRA_GAPS_TSV       — path to data/planA/sra_gaps_filtered.tsv
#   GETORG_THREADS     — CPU threads for GetOrganelle (default: 16)
#   BLAST_THREADS      — CPU threads for BLAST (default: 4)
#   REF_PLASTOME       — path to reference plastome FASTA (Zingiber officinale)
#
# Output contract (for download_sra.sh delete_raw):
#   Prints the GetOrganelle output directory as the LAST line of stdout.
#   This path is passed to delete_raw() to confirm assembly succeeded before
#   deleting raw reads.
#
# Exit codes:
#   0  — assembly + extraction succeeded (even if some markers not extracted)
#   1  — GetOrganelle failed (raw reads are NOT deleted by download_sra.sh)
#
# Per-sample outputs:
#   data/planA/plastomes/<species>/      — GetOrganelle assembly directory
#   data/planA/cds/<marker>/             — extracted marker FASTA files
#   data/planB/genbank_dedup/<marker>/   — merged into Plan B (after extraction)
#
# Log:
#   logs/sra_assembly_<acc>.log          — full per-sample log
# =============================================================================

set -euo pipefail

# ── Arguments ─────────────────────────────────────────────────────────────────
ACC="${1:-}"
READ1="${2:-}"
READ2="${3:-}"

if [[ -z "$ACC" || -z "$READ1" || -z "$READ2" ]]; then
    echo "Usage: assemble_and_extract.sh <acc> <R1.fastq> <R2.fastq>" >&2
    exit 1
fi

# ── Environment defaults ───────────────────────────────────────────────────────
PHYLOSKILLS_ROOT="${PHYLOSKILLS_ROOT:-/nobackup2/bingl/PhyloSkills}"
PROJECT_ROOT="${PROJECT_ROOT:-/nobackup2/bingl/PhyloSkills/Zingiberaceae}"
SRA_GAPS_TSV="${SRA_GAPS_TSV:-$PROJECT_ROOT/data/planA/sra_gaps_filtered.tsv}"
GETORG_THREADS="${GETORG_THREADS:-16}"
BLAST_THREADS="${BLAST_THREADS:-4}"
REF_PLASTOME="${REF_PLASTOME:-$PROJECT_ROOT/data/planA/reference/Zingiber_officinale_NC_037455.1.fasta}"
MARKERS="matK,rbcL,trnL,psbA-trnH,rpoB,rpoC1,atpB,ndhF,ycf1,ycf2"

export PATH="$HOME/miniconda3/bin:$PATH"

LOG_DIR="$PROJECT_ROOT/logs"
mkdir -p "$LOG_DIR"
LOG="$LOG_DIR/sra_assembly_${ACC}.log"

# Redirect all output to log and also stderr to terminal for monitoring
exec > >(tee -a "$LOG") 2>&1

echo "============================================================"
echo "assemble_and_extract.sh — $(date '+%Y-%m-%d %H:%M:%S')"
echo "  Accession: $ACC"
echo "  R1:        $READ1"
echo "  R2:        $READ2"
echo "============================================================"

# ── Step 1: Resolve species name from TSV ─────────────────────────────────────
SPECIES=$(awk -F'\t' -v acc="$ACC" 'NR>1 && $2==acc {print $1; exit}' "$SRA_GAPS_TSV")
if [[ -z "$SPECIES" ]]; then
    echo "ERROR: accession $ACC not found in $SRA_GAPS_TSV" >&2
    exit 1
fi
echo "Species: $SPECIES"

# ── Step 2: GetOrganelle assembly ─────────────────────────────────────────────
PLASTOME_OUTDIR="$PROJECT_ROOT/data/planA/plastomes/$SPECIES"
mkdir -p "$PLASTOME_OUTDIR"

GETORG_OUT="$PLASTOME_OUTDIR/$ACC"

if [[ -d "$GETORG_OUT" ]]; then
    echo "GetOrganelle output already exists: $GETORG_OUT — skipping assembly"
else
    echo ""
    echo "Running GetOrganelle..."
    REF_FLAG=""
    [[ -f "$REF_PLASTOME" ]] && REF_FLAG="-R 10 -w 105 --continue"

    get_organelle_from_reads.py \
        -1 "$READ1" \
        -2 "$READ2" \
        -F embplant_pt \
        -o "$GETORG_OUT" \
        -t "$GETORG_THREADS" \
        -k 21,45,65,85,105 \
        --overwrite \
        2>&1

    echo "GetOrganelle finished."
fi

# ── Step 3: Find assembled plastome FASTA ─────────────────────────────────────
# GetOrganelle outputs *.complete.graph1.1.path_sequence.fasta when successful
# Fall back to any .path_sequence.fasta (partial assembly also usable for BLAST)
PLASTOME_FASTA=""
if [[ -d "$GETORG_OUT" ]]; then
    PLASTOME_FASTA=$(find "$GETORG_OUT" -name "*.path_sequence.fasta" \
                         -not -name "*.selected_graph*" \
                         2>/dev/null | sort | head -1 || true)
fi

if [[ -z "$PLASTOME_FASTA" || ! -f "$PLASTOME_FASTA" ]]; then
    echo "ERROR: No plastome assembly found for $ACC ($SPECIES)." >&2
    echo "  GetOrganelle may have failed — check: $GETORG_OUT/get_org.log.txt" >&2
    # Still print GETORG_OUT as output path (delete_raw will check it's non-empty)
    echo "$GETORG_OUT"
    exit 1
fi

echo "Assembly: $PLASTOME_FASTA"
PLASTOME_LEN=$(awk '/^>/{next}{len+=length($0)}END{print len}' "$PLASTOME_FASTA")
echo "Assembly length: ${PLASTOME_LEN} bp"

# Warn if assembly is too short (likely incomplete)
if [[ "$PLASTOME_LEN" -lt 100000 ]]; then
    echo "WARNING: Assembly is short (${PLASTOME_LEN} bp < 100,000 bp). " \
         "Markers may not be recoverable. Proceeding anyway."
fi

# ── Step 4: Extract markers via BLAST ─────────────────────────────────────────
echo ""
echo "Extracting markers via BLAST..."
CDS_DIR="$PROJECT_ROOT/data/planA/cds"
python "$PHYLOSKILLS_ROOT/scripts/data/extract_markers_blast.py" \
    --plastome    "$PLASTOME_FASTA" \
    --ref_dir     "$PROJECT_ROOT/data/planB/genbank_dedup" \
    --ref_plastome "$REF_PLASTOME" \
    --markers     "$MARKERS" \
    --species     "$SPECIES" \
    --run_id      "$ACC" \
    --output      "$CDS_DIR" \
    --threads     "$BLAST_THREADS" \
    --provenance  "$PROJECT_ROOT/results/planA/provenance"

# ── Step 5: Merge extracted sequences into Plan B ──────────────────────────────
echo ""
echo "Merging extracted sequences into Plan B genbank_dedup..."
python "$PHYLOSKILLS_ROOT/scripts/data/merge_sra_cds_to_planB.py" \
    --cds_dir   "$CDS_DIR" \
    --planB_dir "$PROJECT_ROOT/data/planB/genbank_dedup" \
    --markers   "$MARKERS" \
    --provenance "$PROJECT_ROOT/results/planA/provenance"

echo ""
echo "============================================================"
echo "Completed: $SPECIES ($ACC) — $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"

# IMPORTANT: Print assembly output directory as final line of stdout.
# download_sra.sh streaming mode reads this to know where assembly landed
# before deleting raw reads.
echo "$GETORG_OUT"
