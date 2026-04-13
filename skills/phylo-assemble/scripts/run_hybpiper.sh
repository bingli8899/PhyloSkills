#!/usr/bin/env bash
# =============================================================================
# run_hybpiper.sh — Target enrichment (HybSeq) assembly with HybPiper v2
# =============================================================================
# Runs the full HybPiper v2 pipeline:
#   1. hybpiper assemble  — per-sample target assembly
#   2. hybpiper stats     — recovery statistics
#   3. hybpiper retrieve_sequences — collect sequences across samples
#
# Usage:
#   bash run_hybpiper.sh \
#     -r <target_refs.fasta> \
#     -s <sample_list.txt> \
#     -d <reads_dir> \
#     -o <outdir> \
#     [-t <threads>] [-m <min_length>] [-T <seq_type>]
#
# Arguments:
#   -r  Target reference FASTA (DNA or amino acid sequences for target loci)
#   -s  Sample list: one sample name per line (must match read file prefixes)
#   -d  Directory containing paired reads:
#         <reads_dir>/<sample>_R1.fastq[.gz] and <sample>_R2.fastq[.gz]
#   -o  Output directory
#   -t  CPU threads per sample (default: auto-detect via nproc)
#   -m  Minimum percent target length for sequence retrieval (default: 10)
#   -T  Sequence type for retrieval: dna | aa | intron | supercontig (default: dna)
#   -e  Run hybpiper paralog_retriever after retrieve_sequences (flag)
#
# Tool version requirements (checked at runtime):
#   HybPiper >= 2.3   — https://github.com/mossmatters/HybPiper
#     Latest release: v2.3.4 (2025)
#     Install: conda install -c bioconda hybpiper
#   NOTE: HybPiper v2 uses different commands from v1.
#         v1 used: reads_first.py, retrieve_sequences.py
#         v2 uses: hybpiper assemble, hybpiper retrieve_sequences
#
# QC thresholds (SKILL.md):
#   - Target recovery >= 50% of loci per sample considered acceptable
#   - Samples with < 20% recovery should be investigated or excluded
#   - Check hybpiper stats output for per-locus and per-sample recovery
#
# =============================================================================

set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────
THREADS=$(nproc 2>/dev/null || echo 8)
MIN_LENGTH=10
SEQ_TYPE="dna"
RUN_PARALOGS=false

# ── Argument parsing ──────────────────────────────────────────────────────────
usage() {
    grep '^#' "$0" | grep -v '#!/' | sed 's/^# \{0,1\}//'
    exit 1
}

while getopts "r:s:d:o:t:m:T:eh" opt; do
    case $opt in
        r) TARGET_REFS="$OPTARG" ;;
        s) SAMPLE_LIST="$OPTARG" ;;
        d) READS_DIR="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        m) MIN_LENGTH="$OPTARG" ;;
        T) SEQ_TYPE="$OPTARG" ;;
        e) RUN_PARALOGS=true ;;
        h|*) usage ;;
    esac
done

[[ -z "${TARGET_REFS:-}" || -z "${SAMPLE_LIST:-}" || \
   -z "${READS_DIR:-}" || -z "${OUTDIR:-}" ]] && {
    echo "ERROR: -r, -s, -d, and -o are required." >&2; usage; }

[[ ! -f "$TARGET_REFS" ]] && { echo "ERROR: Target refs not found: ${TARGET_REFS}" >&2; exit 1; }
[[ ! -f "$SAMPLE_LIST" ]] && { echo "ERROR: Sample list not found: ${SAMPLE_LIST}" >&2; exit 1; }
[[ ! -d "$READS_DIR" ]]   && { echo "ERROR: Reads dir not found: ${READS_DIR}" >&2; exit 1; }

# ── Version check ─────────────────────────────────────────────────────────────
command -v hybpiper &>/dev/null || {
    echo "ERROR: HybPiper not found. Install:" >&2
    echo "  conda install -c bioconda hybpiper" >&2
    exit 1
}

HP_VERSION=$(hybpiper --version 2>&1 | head -1 || echo "unknown")
echo "# HybPiper version: ${HP_VERSION}  (latest known: 2.3.4)"
echo "# Target refs: ${TARGET_REFS}"
echo "# Threads per sample: ${THREADS}"
echo "# Seq type: ${SEQ_TYPE}"
echo "# Date: $(date +%Y-%m-%d)"
echo ""

# Warn if using an old HybPiper installation
if echo "$HP_VERSION" | grep -qE '^1\.' 2>/dev/null; then
    echo "ERROR: HybPiper v1.x detected. This script requires HybPiper v2.x." >&2
    echo "  v2 changed commands: 'reads_first.py' → 'hybpiper assemble'" >&2
    echo "  Update: conda install -c bioconda 'hybpiper>=2'" >&2
    exit 1
fi

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# ── Step 1 — Per-sample assembly ──────────────────────────────────────────────
echo "=== Step 1: Per-sample target assembly ==="
echo ""

FAILED_SAMPLES=()

while read -r sample; do
    [[ -z "$sample" || "$sample" =~ ^# ]] && continue

    echo "--- Processing: ${sample} ---"

    # Locate read files (accept .fastq or .fastq.gz, _R1/_1 suffix variants)
    R1=""
    R2=""
    for suffix in "_R1.fastq.gz" "_R1.fastq" "_1.fastq.gz" "_1.fastq"; do
        [[ -f "${READS_DIR}/${sample}${suffix}" ]] && R1="${READS_DIR}/${sample}${suffix}" && break
    done
    for suffix in "_R2.fastq.gz" "_R2.fastq" "_2.fastq.gz" "_2.fastq"; do
        [[ -f "${READS_DIR}/${sample}${suffix}" ]] && R2="${READS_DIR}/${sample}${suffix}" && break
    done

    if [[ -z "$R1" || -z "$R2" ]]; then
        echo "  WARNING: Read files not found for ${sample} in ${READS_DIR}" >&2
        FAILED_SAMPLES+=("$sample")
        continue
    fi

    hybpiper assemble \
        --targetfile_dna "$TARGET_REFS" \
        --readfiles "$R1" "$R2" \
        --prefix "$sample" \
        --cpu "$THREADS" \
        --run_intronerate \
        2>&1 | tee "${sample}_hybpiper.log"

    echo "  Done: ${sample}"
    echo ""

done < "$SAMPLE_LIST"

if [[ "${#FAILED_SAMPLES[@]}" -gt 0 ]]; then
    echo "WARNING: The following samples had missing reads and were skipped:"
    printf '  %s\n' "${FAILED_SAMPLES[@]}"
    echo ""
fi

# ── Step 2 — Recovery statistics ─────────────────────────────────────────────
echo "=== Step 2: Computing recovery statistics ==="

hybpiper stats \
    --targetfile_dna "$TARGET_REFS" \
    --stats_filename hybpiper_stats \
    --sequence_type "$SEQ_TYPE" \
    $(grep -v '^#' "$SAMPLE_LIST" | grep -v '^$') \
    2>&1 | tee hybpiper_stats_run.log

echo ""
echo "  Stats written to: ${OUTDIR}/hybpiper_stats.tsv"

# Quick recovery summary
if [[ -f "hybpiper_stats.tsv" ]]; then
    echo ""
    echo "  Recovery summary (% loci recovered per sample):"
    awk -F'\t' 'NR>1 {
        if($3>0) pct=sprintf("%.0f%%", $2/$3*100); else pct="N/A";
        printf "    %-30s %s\n", $1, pct
    }' hybpiper_stats.tsv
fi

# ── Step 3 — Retrieve sequences ───────────────────────────────────────────────
echo ""
echo "=== Step 3: Retrieving sequences across samples ==="

hybpiper retrieve_sequences \
    --targetfile_dna "$TARGET_REFS" \
    --sample_names "$SAMPLE_LIST" \
    --fasta_dir retrieved_sequences \
    --sequence_type "$SEQ_TYPE" \
    --min_length_percentage "$MIN_LENGTH" \
    2>&1 | tee hybpiper_retrieve.log

echo ""
N_RETRIEVED=$(find retrieved_sequences/ -name "*.fasta" 2>/dev/null | wc -l || echo 0)
echo "  Retrieved ${N_RETRIEVED} locus FASTA files to: ${OUTDIR}/retrieved_sequences/"

# ── Step 4 — Paralog retriever (optional) ────────────────────────────────────
if [[ "$RUN_PARALOGS" == true ]]; then
    echo ""
    echo "=== Step 4: Paralog retriever ==="
    hybpiper paralog_retriever \
        --targetfile_dna "$TARGET_REFS" \
        --sample_names "$SAMPLE_LIST" \
        2>&1 | tee hybpiper_paralogs.log
    echo "  Paralog output: ${OUTDIR}/paralogs/"
    echo ""
    echo "  NOTE: Multi-copy loci detected — use ASTRAL-Pro3 for coalescent analysis"
fi

# ── QC gate ───────────────────────────────────────────────────────────────────
echo ""
echo "=== QC Summary ==="
echo "  Stats file:       ${OUTDIR}/hybpiper_stats.tsv"
echo "  Retrieved loci:   ${OUTDIR}/retrieved_sequences/"
echo "  Next step:        phylo-alignment (align each locus in retrieved_sequences/)"
echo ""
echo "  QC thresholds:"
echo "    - >= 50% loci recovered per sample: acceptable"
echo "    - < 20% loci recovered: investigate or exclude sample"
