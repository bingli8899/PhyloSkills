#!/usr/bin/env bash
# =============================================================================
# assemble_plastome_getorganelle.sh — Plastome assembly with GetOrganelle
# =============================================================================
# Usage:
#   bash assemble_plastome_getorganelle.sh \
#     -1 <reads_R1.fastq> -2 <reads_R2.fastq> \
#     -o <outdir> -s <sample_name> \
#     [-t <threads>] [-k <kmer_list>] [-F <plant_type>]
#
# Arguments:
#   -1  Forward reads (FASTQ, gzipped or plain)
#   -2  Reverse reads (FASTQ, gzipped or plain)
#   -o  Output base directory (sample subdirectory created inside)
#   -s  Sample name (used for subdirectory and output prefix)
#   -t  CPU threads (default: 8; check available CPUs first with nproc)
#   -k  Comma-separated k-mer sizes (default: 21,45,65,85,105)
#   -F  GetOrganelle target type (default: embplant_pt)
#         embplant_pt  = land plant plastome
#         embplant_mt  = land plant mitochondrial genome
#         animal_mt    = animal mitochondrial genome
#         fungus_mt    = fungal mitochondrial genome
#         other_pt     = other plastome
#   -R  Reference FASTA for -w seed reads (optional, improves assembly)
#   -P  Max reads to use (default: GetOrganelle default ~75000)
#
# Tool version requirements (checked at runtime):
#   GetOrganelle >= 1.7.7  — https://github.com/Kinggerm/GetOrganelle
#     Latest release: 1.7.7.1 (2024)
#     Install: pip install getorganelle
#     or: conda install -c bioconda getorganelle
#   Also requires: SPAdes, Bowtie2, BLAST+ (installed with GetOrganelle)
#
# QC thresholds (SKILL.md):
#   - Plastome completeness >= 90% of expected length
#   - If assembly yields multiple circular paths, inspect graph in Bandage
#   - Mean coverage depth >= 10× (check GetOrganelle *.graph1.fastg or log)
#
# Output:
#   <outdir>/<sample_name>/
#     *.complete.graph1.1.path_sequence.fasta  — primary circular assembly
#     *.complete.graph1.selected_graph.gfa     — assembly graph (view in Bandage)
#     get_org.log.txt                          — full log
#
# =============================================================================

set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────
THREADS=$(nproc 2>/dev/null || echo 8)
KMER="21,45,65,85,105"
PLANT_TYPE="embplant_pt"
REFERENCE=""
MAX_READS=""

# ── Argument parsing ──────────────────────────────────────────────────────────
usage() {
    grep '^#' "$0" | grep -v '#!/' | sed 's/^# \{0,1\}//'
    exit 1
}

while getopts "1:2:o:s:t:k:F:R:P:h" opt; do
    case $opt in
        1) READ1="$OPTARG" ;;
        2) READ2="$OPTARG" ;;
        o) OUTBASE="$OPTARG" ;;
        s) SAMPLE="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        k) KMER="$OPTARG" ;;
        F) PLANT_TYPE="$OPTARG" ;;
        R) REFERENCE="$OPTARG" ;;
        P) MAX_READS="$OPTARG" ;;
        h|*) usage ;;
    esac
done

[[ -z "${READ1:-}" || -z "${READ2:-}" || -z "${OUTBASE:-}" || -z "${SAMPLE:-}" ]] && {
    echo "ERROR: -1, -2, -o, and -s are required." >&2; usage; }
[[ ! -f "$READ1" ]] && { echo "ERROR: Read1 not found: ${READ1}" >&2; exit 1; }
[[ ! -f "$READ2" ]] && { echo "ERROR: Read2 not found: ${READ2}" >&2; exit 1; }

# ── Version check ─────────────────────────────────────────────────────────────
command -v get_organelle_from_reads.py &>/dev/null || {
    echo "ERROR: GetOrganelle not found. Install:" >&2
    echo "  pip install getorganelle" >&2
    echo "  or: conda install -c bioconda getorganelle" >&2
    exit 1
}

GO_VERSION=$(get_organelle_from_reads.py --version 2>&1 | head -1 || echo "unknown")
echo "# GetOrganelle version: ${GO_VERSION}"
echo "# Sample: ${SAMPLE}"
echo "# Target type: ${PLANT_TYPE}"
echo "# Threads: ${THREADS}"
echo "# Date: $(date +%Y-%m-%d)"
echo ""

# ── Warn if version is outdated ───────────────────────────────────────────────
# Latest known release: 1.7.7.1 (checked 2026-04-13)
# If this script reports an older version, update GetOrganelle before proceeding.
echo "# NOTE: Latest known GetOrganelle release is 1.7.7.1 (2026-04-13)."
echo "#       If version above is older, consider updating: pip install -U getorganelle"
echo ""

# ── Setup ─────────────────────────────────────────────────────────────────────
OUTDIR="${OUTBASE}/${SAMPLE}"
mkdir -p "$OUTDIR"

# ── Build command ─────────────────────────────────────────────────────────────
CMD="get_organelle_from_reads.py"
CMD+=" -1 ${READ1}"
CMD+=" -2 ${READ2}"
CMD+=" -F ${PLANT_TYPE}"
CMD+=" -o ${OUTDIR}"
CMD+=" -t ${THREADS}"
CMD+=" -k ${KMER}"
CMD+=" --overwrite"

[[ -n "$REFERENCE" ]] && CMD+=" -s ${REFERENCE}"
[[ -n "$MAX_READS" ]] && CMD+=" -P ${MAX_READS}"

echo "=== Running GetOrganelle ==="
echo "Command: ${CMD}"
echo ""

eval "$CMD"

# ── Verify output ─────────────────────────────────────────────────────────────
echo ""
echo "=== Assembly QC ==="

FASTA_OUT=$(find "${OUTDIR}" -name "*.path_sequence.fasta" 2>/dev/null | head -1 || true)
if [[ -z "$FASTA_OUT" ]]; then
    echo "ERROR: Assembly failed — no path_sequence.fasta found in ${OUTDIR}" >&2
    echo "  Check ${OUTDIR}/get_org.log.txt for details"
    exit 1
fi

echo "  Primary assembly: ${FASTA_OUT}"

# Check assembly length (rough plastome size: 100–170 kb for land plants)
SEQ_LEN=$(awk '/^>/{next} {len+=length($0)} END {print len}' "$FASTA_OUT" || echo "0")
echo "  Assembly length: ${SEQ_LEN} bp"

if [[ "$SEQ_LEN" -lt 90000 ]]; then
    echo "  WARNING: Assembly length < 90,000 bp — may be incomplete"
    echo "           Expected: 100,000–170,000 bp for land plant plastome"
elif [[ "$SEQ_LEN" -gt 180000 ]]; then
    echo "  WARNING: Assembly length > 180,000 bp — possible chimera or contamination"
else
    echo "  QC PASS: Assembly length within expected range for land plant plastome"
fi

# Check for multiple paths (indicates unresolved repeats)
N_PATHS=$(find "${OUTDIR}" -name "*.path_sequence.fasta" 2>/dev/null | wc -l)
if [[ "$N_PATHS" -gt 1 ]]; then
    echo "  WARNING: ${N_PATHS} assembly paths found (unresolved inverted repeat)"
    echo "           Inspect assembly graph: ${OUTDIR}/*.gfa in Bandage"
fi

# ── Return output path for streaming mode integration ────────────────────────
echo ""
echo "Output: ${FASTA_OUT}"
