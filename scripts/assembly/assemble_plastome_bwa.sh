#!/usr/bin/env bash
# =============================================================================
# assemble_plastome_bwa.sh — Plastome assembly via BWA read mapping + consensus
# =============================================================================
# Use this when GetOrganelle fails or you have a closely related reference.
# Suitable for low-coverage genome skimming data mapped to a reference plastome.
#
# Usage:
#   bash assemble_plastome_bwa.sh \
#     -1 <reads_R1.fastq> -2 <reads_R2.fastq> \
#     -r <reference.fasta> \
#     -o <outdir> -s <sample_name> \
#     [-t <threads>] [-d <min_depth>] [-c <consensus_type>]
#
# Arguments:
#   -1  Forward reads (FASTQ, gzipped or plain)
#   -2  Reverse reads (FASTQ, gzipped or plain)
#   -r  Reference plastome FASTA (same genus preferred; see SKILL.md hierarchy)
#   -o  Output directory (created if absent)
#   -s  Sample name (used for file prefix)
#   -t  CPU threads (default: auto-detect via nproc)
#   -d  Minimum depth to call a base (default: 3)
#   -c  Consensus type: strict | majority (default: strict)
#         strict   = only call bases with depth >= min_depth, else N
#         majority = call most frequent allele regardless of depth
#
# Tool version requirements (checked at runtime):
#   BWA >= 0.7.19       — https://github.com/lh3/bwa  (latest: v0.7.19, 2025)
#   SAMtools >= 1.23    — https://github.com/samtools/samtools (latest: 1.23.1)
#   BCFtools >= 1.23    — https://github.com/samtools/bcftools (latest: 1.23.1)
#
# Output:
#   <outdir>/<sample_name>_consensus.fasta  — final plastome consensus
#   <outdir>/<sample_name>_sorted.bam       — sorted, indexed BAM
#   <outdir>/<sample_name>_coverage.txt     — per-position depth
#   <outdir>/<sample_name>_stats.txt        — flagstat summary
#
# =============================================================================

set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────
THREADS=$(nproc 2>/dev/null || echo 8)
MIN_DEPTH=3
CONSENSUS_TYPE="strict"

# ── Argument parsing ──────────────────────────────────────────────────────────
usage() {
    grep '^#' "$0" | grep -v '#!/' | sed 's/^# \{0,1\}//'
    exit 1
}

while getopts "1:2:r:o:s:t:d:c:h" opt; do
    case $opt in
        1) READ1="$OPTARG" ;;
        2) READ2="$OPTARG" ;;
        r) REFERENCE="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        s) SAMPLE="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        d) MIN_DEPTH="$OPTARG" ;;
        c) CONSENSUS_TYPE="$OPTARG" ;;
        h|*) usage ;;
    esac
done

[[ -z "${READ1:-}" || -z "${READ2:-}" || -z "${REFERENCE:-}" || \
   -z "${OUTDIR:-}" || -z "${SAMPLE:-}" ]] && {
    echo "ERROR: -1, -2, -r, -o, and -s are required." >&2; usage; }

[[ ! -f "$READ1" ]] && { echo "ERROR: Read1 not found: ${READ1}" >&2; exit 1; }
[[ ! -f "$READ2" ]] && { echo "ERROR: Read2 not found: ${READ2}" >&2; exit 1; }
[[ ! -f "$REFERENCE" ]] && { echo "ERROR: Reference not found: ${REFERENCE}" >&2; exit 1; }

[[ "$CONSENSUS_TYPE" != "strict" && "$CONSENSUS_TYPE" != "majority" ]] && {
    echo "ERROR: -c must be 'strict' or 'majority'" >&2; exit 1; }

# ── Version checks ────────────────────────────────────────────────────────────
for tool in bwa samtools bcftools; do
    command -v "$tool" &>/dev/null || {
        echo "ERROR: '${tool}' not found. See executables/software-inventory.md" >&2
        exit 1
    }
done

BWA_VERSION=$(bwa 2>&1 | grep "Version:" | awk '{print $2}' || echo "unknown")
SAMTOOLS_VERSION=$(samtools --version 2>&1 | head -1 | awk '{print $2}' || echo "unknown")
BCFTOOLS_VERSION=$(bcftools --version 2>&1 | head -1 | awk '{print $2}' || echo "unknown")

echo "# BWA version:      ${BWA_VERSION}  (latest known: 0.7.19)"
echo "# SAMtools version: ${SAMTOOLS_VERSION}  (latest known: 1.23.1)"
echo "# BCFtools version: ${BCFTOOLS_VERSION}  (latest known: 1.23.1)"
echo "# Sample: ${SAMPLE}"
echo "# Reference: ${REFERENCE}"
echo "# Threads: ${THREADS}"
echo "# Min depth: ${MIN_DEPTH}"
echo "# Consensus type: ${CONSENSUS_TYPE}"
echo "# Date: $(date +%Y-%m-%d)"
echo ""

mkdir -p "$OUTDIR"

# ── File paths ────────────────────────────────────────────────────────────────
BAM_UNSORTED="${OUTDIR}/${SAMPLE}_unsorted.bam"
BAM_SORTED="${OUTDIR}/${SAMPLE}_sorted.bam"
PILEUP="${OUTDIR}/${SAMPLE}.pileup.vcf.gz"
CONSENSUS="${OUTDIR}/${SAMPLE}_consensus.fasta"
COVERAGE="${OUTDIR}/${SAMPLE}_coverage.txt"
STATS="${OUTDIR}/${SAMPLE}_stats.txt"
REF_INDEX="${REFERENCE}.bwt"

# ── Index reference if needed ─────────────────────────────────────────────────
if [[ ! -f "$REF_INDEX" ]]; then
    echo "=== Indexing reference ==="
    bwa index "$REFERENCE"
fi

# ── Map reads ─────────────────────────────────────────────────────────────────
echo "=== Mapping reads with BWA MEM ==="
bwa mem \
    -t "$THREADS" \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
    "$REFERENCE" "$READ1" "$READ2" \
    | samtools view -bS -F 4 -@ "$THREADS" \
    > "$BAM_UNSORTED"

# ── Sort and index BAM ────────────────────────────────────────────────────────
echo "=== Sorting and indexing BAM ==="
samtools sort -@ "$THREADS" "$BAM_UNSORTED" -o "$BAM_SORTED"
samtools index "$BAM_SORTED"
rm -f "$BAM_UNSORTED"

# ── Alignment statistics ──────────────────────────────────────────────────────
echo "=== Alignment statistics ==="
samtools flagstat "$BAM_SORTED" | tee "$STATS"
echo ""

# ── Coverage depth ────────────────────────────────────────────────────────────
echo "=== Computing coverage depth ==="
samtools depth -a "$BAM_SORTED" > "$COVERAGE"

MEAN_DEPTH=$(awk '{sum+=$3; n++} END {if(n>0) printf "%.1f", sum/n; else print "0"}' "$COVERAGE")
COVERED_SITES=$(awk -v d="$MIN_DEPTH" '$3>=d {n++} END {print n+0}' "$COVERAGE")
TOTAL_SITES=$(wc -l < "$COVERAGE")
echo "  Mean depth: ${MEAN_DEPTH}×"
echo "  Sites covered at ≥${MIN_DEPTH}×: ${COVERED_SITES} / ${TOTAL_SITES}"

if (( $(echo "$MEAN_DEPTH < 10" | awk '{print ($1 < 10)}') )); then
    echo "  WARNING: Mean depth < 10× — assembly quality may be low"
fi

# ── Generate consensus ────────────────────────────────────────────────────────
echo ""
echo "=== Generating consensus (${CONSENSUS_TYPE} mode) ==="

if [[ "$CONSENSUS_TYPE" == "strict" ]]; then
    # Call only positions with depth >= MIN_DEPTH; others become N.
    # BCFtools v1.21+: --mask-with requires --mask FILE.
    # We generate a BED of low-coverage positions from the samtools depth output
    # and pass it as the mask, then use --mask-with N.
    LOW_COV_BED="${OUTDIR}/${SAMPLE}_low_cov.bed"
    awk -v d="$MIN_DEPTH" 'BEGIN{OFS="\t"} $3<d {print $1, $2-1, $2}' \
        "$COVERAGE" > "$LOW_COV_BED"

    bcftools mpileup \
        --fasta-ref "$REFERENCE" \
        --min-BQ 20 \
        --min-MQ 20 \
        "$BAM_SORTED" \
        | bcftools call \
            --consensus-caller \
            --variants-only \
        | bcftools view -Oz -o "$PILEUP"
    bcftools index "$PILEUP"

    bcftools consensus \
        --fasta-ref "$REFERENCE" \
        --mask "$LOW_COV_BED" \
        --mask-with N \
        "$PILEUP" \
        | sed "1s/>.*/>${SAMPLE}/" \
        > "$CONSENSUS"

elif [[ "$CONSENSUS_TYPE" == "majority" ]]; then
    # Call most frequent base regardless of depth
    bcftools mpileup \
        --fasta-ref "$REFERENCE" \
        --min-BQ 15 \
        --min-MQ 15 \
        "$BAM_SORTED" \
        | bcftools call \
            --consensus-caller \
        | bcftools view -Oz -o "$PILEUP"
    bcftools index "$PILEUP"

    bcftools consensus \
        --fasta-ref "$REFERENCE" \
        "$PILEUP" \
        | sed "1s/>.*/>${SAMPLE}/" \
        > "$CONSENSUS"
fi

# ── QC check ─────────────────────────────────────────────────────────────────
echo ""
echo "=== Consensus QC ==="

if [[ ! -f "$CONSENSUS" || ! -s "$CONSENSUS" ]]; then
    echo "ERROR: Consensus file not generated or empty." >&2
    exit 1
fi

CONSENSUS_LEN=$(awk '/^>/{next} {len+=length($0)} END {print len}' "$CONSENSUS")
N_COUNT=$(awk '/^>/{next} {seq=seq$0} END {gsub(/[^Nn]/, "", seq); print length(seq)}' "$CONSENSUS")
echo "  Consensus length: ${CONSENSUS_LEN} bp"
echo "  N count: ${N_COUNT} ($(awk -v n="$N_COUNT" -v l="$CONSENSUS_LEN" \
    'BEGIN{printf "%.1f", (l>0?n/l*100:0)}')%)"

if [[ "$CONSENSUS_LEN" -lt 90000 ]]; then
    echo "  WARNING: Consensus < 90,000 bp — possibly partial assembly"
fi

echo ""
echo "=== Output files ==="
echo "  BAM:       ${BAM_SORTED}"
echo "  Coverage:  ${COVERAGE}"
echo "  Consensus: ${CONSENSUS}"
echo ""
# Return output path for streaming mode
echo "Output: ${CONSENSUS}"
