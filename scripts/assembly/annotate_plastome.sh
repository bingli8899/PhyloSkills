#!/usr/bin/env bash
# =============================================================================
# annotate_plastome.sh — Annotate plastome assemblies with PLANN or chloe
# =============================================================================
# Annotates one or more plastome FASTA files using PLANN (primary) or
# chloe (alternative). Both tools are auto-detected if in PATH or executables/.
#
# Usage:
#   bash annotate_plastome.sh \
#     -i <input_dir_or_fasta> \
#     -o <outdir> \
#     [-T <tool>] \
#     [-p <plann_path>] \
#     [-c <chloe_path>] \
#     [-g <genome_size_min>]
#
# Arguments:
#   -i  Input: single FASTA file OR directory containing *.fasta files
#   -o  Output directory for annotation results
#   -T  Tool to use: plann | chloe | auto (default: auto — tries plann first)
#   -p  Path to PLANN executable or script (if not in PATH)
#         GitHub: https://github.com/ian-small/PLANN
#   -c  Path to chloe executable (if not in PATH)
#         GitHub: https://github.com/ian-small/chloe
#         chloe is a Julia tool; run as: julia chloe.jl or via compiled binary
#   -g  Minimum assembly length to attempt annotation (default: 50000 bp)
#
# Tool version requirements:
#   PLANN   — https://github.com/ian-small/PLANN  (Python, pip install)
#   chloe   — https://github.com/ian-small/chloe  (Julia; requires Julia >= 1.9)
#   At least one must be available.
#
# Output per assembly:
#   <outdir>/<sample>.gb    — GenBank format annotation
#   <outdir>/<sample>.gff3  — GFF3 annotation (if available from tool)
#   <outdir>/<sample>.sff   — SFF format (chloe only)
#
# =============================================================================

set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────
TOOL="auto"
PLANN_BIN=""
CHLOE_BIN=""
MIN_LEN=50000

# ── Argument parsing ──────────────────────────────────────────────────────────
usage() {
    grep '^#' "$0" | grep -v '#!/' | sed 's/^# \{0,1\}//'
    exit 1
}

while getopts "i:o:T:p:c:g:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        T) TOOL="$OPTARG" ;;
        p) PLANN_BIN="$OPTARG" ;;
        c) CHLOE_BIN="$OPTARG" ;;
        g) MIN_LEN="$OPTARG" ;;
        h|*) usage ;;
    esac
done

[[ -z "${INPUT:-}" || -z "${OUTDIR:-}" ]] && {
    echo "ERROR: -i and -o are required." >&2; usage; }

# ── Tool detection ────────────────────────────────────────────────────────────
detect_plann() {
    if [[ -n "$PLANN_BIN" && -x "$PLANN_BIN" ]]; then
        echo "$PLANN_BIN"; return 0
    fi
    local found
    found=$(command -v plann 2>/dev/null || \
            command -v PLANN.py 2>/dev/null || \
            find executables/plann/ -name "*.py" -o -name "plann" 2>/dev/null | head -1 || \
            echo "")
    echo "$found"
}

detect_chloe() {
    if [[ -n "$CHLOE_BIN" && -x "$CHLOE_BIN" ]]; then
        echo "$CHLOE_BIN"; return 0
    fi
    local found
    found=$(command -v chloe 2>/dev/null || \
            find executables/chloe/ -name "chloe" 2>/dev/null | head -1 || \
            echo "")
    echo "$found"
}

PLANN_PATH=$(detect_plann)
CHLOE_PATH=$(detect_chloe)

# Resolve tool choice
if [[ "$TOOL" == "auto" ]]; then
    if [[ -n "$PLANN_PATH" ]]; then
        TOOL="plann"
    elif [[ -n "$CHLOE_PATH" ]]; then
        TOOL="chloe"
    else
        echo "ERROR: Neither PLANN nor chloe found." >&2
        echo "  Install PLANN:  pip install plann" >&2
        echo "  Install chloe:  see https://github.com/ian-small/chloe" >&2
        echo "  Or specify paths with -p/-c" >&2
        exit 1
    fi
fi

if [[ "$TOOL" == "plann" && -z "$PLANN_PATH" ]]; then
    echo "ERROR: PLANN not found. Install: pip install plann" >&2
    echo "  Or specify path with -p <path>" >&2
    exit 1
fi

if [[ "$TOOL" == "chloe" && -z "$CHLOE_PATH" ]]; then
    echo "ERROR: chloe not found." >&2
    echo "  See: https://github.com/ian-small/chloe" >&2
    echo "  Or specify path with -c <path>" >&2
    exit 1
fi

echo "# Annotation tool: ${TOOL}"
echo "# PLANN path: ${PLANN_PATH:-not used}"
echo "# chloe path: ${CHLOE_PATH:-not used}"
echo "# Date: $(date +%Y-%m-%d)"
echo ""

mkdir -p "$OUTDIR"

# ── Build file list ───────────────────────────────────────────────────────────
if [[ -f "$INPUT" ]]; then
    FASTA_FILES=("$INPUT")
elif [[ -d "$INPUT" ]]; then
    mapfile -t FASTA_FILES < <(find "$INPUT" -maxdepth 2 -name "*.fasta" -o -name "*.fa" | sort)
else
    echo "ERROR: Input not found: ${INPUT}" >&2; exit 1
fi

echo "Found ${#FASTA_FILES[@]} FASTA file(s) to annotate."
echo ""

# ── Annotation functions ──────────────────────────────────────────────────────
annotate_with_plann() {
    local fasta="$1"
    local base
    base=$(basename "$fasta" .fasta)
    base=$(basename "$base" .fa)

    local seq_len
    seq_len=$(awk '/^>/{next} {len+=length($0)} END {print len+0}' "$fasta")
    if [[ "$seq_len" -lt "$MIN_LEN" ]]; then
        echo "  SKIP ${base}: sequence length ${seq_len} bp < minimum ${MIN_LEN} bp"
        return
    fi

    echo "  Annotating: ${base} (${seq_len} bp)"
    "$PLANN_PATH" \
        --input "$fasta" \
        --output "${OUTDIR}/${base}" \
        2>&1 | tee "${OUTDIR}/${base}_plann.log" || {
        echo "  WARNING: PLANN failed for ${base}" >&2
        return 1
    }
    echo "  Done: ${OUTDIR}/${base}.gb"
}

annotate_with_chloe() {
    local fasta="$1"
    local base
    base=$(basename "$fasta" .fasta)
    base=$(basename "$base" .fa)

    local seq_len
    seq_len=$(awk '/^>/{next} {len+=length($0)} END {print len+0}' "$fasta")
    if [[ "$seq_len" -lt "$MIN_LEN" ]]; then
        echo "  SKIP ${base}: sequence length ${seq_len} bp < minimum ${MIN_LEN} bp"
        return
    fi

    echo "  Annotating: ${base} (${seq_len} bp)"

    # chloe can be invoked as binary or julia script
    if echo "$CHLOE_PATH" | grep -q "\.jl$"; then
        julia "$CHLOE_PATH" annotate \
            --output "${OUTDIR}/${base}.sff" \
            "$fasta" \
            2>&1 | tee "${OUTDIR}/${base}_chloe.log" || {
            echo "  WARNING: chloe failed for ${base}" >&2; return 1; }
    else
        "$CHLOE_PATH" annotate \
            --output "${OUTDIR}/${base}.sff" \
            "$fasta" \
            2>&1 | tee "${OUTDIR}/${base}_chloe.log" || {
            echo "  WARNING: chloe failed for ${base}" >&2; return 1; }
    fi
    echo "  Done: ${OUTDIR}/${base}.sff"
}

# ── Main loop ─────────────────────────────────────────────────────────────────
SUCCESS=0
FAILED=0

for fasta in "${FASTA_FILES[@]}"; do
    if [[ "$TOOL" == "plann" ]]; then
        annotate_with_plann "$fasta" && (( SUCCESS++ )) || (( FAILED++ ))
    else
        annotate_with_chloe "$fasta" && (( SUCCESS++ )) || (( FAILED++ ))
    fi
done

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "=== Annotation complete ==="
echo "  Succeeded: ${SUCCESS}"
echo "  Failed:    ${FAILED}"
echo "  Output dir: ${OUTDIR}/"
echo ""
echo "Next step: extract_cds.py — extract coding sequences from annotations"
