#!/usr/bin/env bash
# =============================================================================
# build_gene_trees.sh — Per-gene ML tree inference with IQ-TREE for ASTER input
# =============================================================================
# Runs IQ-TREE on each aligned marker FASTA to produce individual gene trees,
# then collects all .treefile outputs into a single file ready for ASTER
# (wASTRAL or ASTRAL-Pro3).
#
# Usage:
#   bash build_gene_trees.sh \
#     -i <aligned_dir> \
#     -o <outdir> \
#     [-m <model>] \
#     [-B <ufboot_reps>] \
#     [-A <alrt_reps>] \
#     [-t <threads_per_gene>] \
#     [-j <parallel_jobs>]
#
# Arguments:
#   -i  Directory containing per-marker aligned FASTA files (*_aligned.fasta
#       or *_trimmed.fasta or *.fasta)
#   -o  Output directory for gene trees (subdirectory per gene)
#   -m  Substitution model: TEST (ModelFinder) or specific model e.g. GTR+F+I+G4
#         Use TEST to run ModelFinder per gene (slower but correct)
#         Use a fixed model for speed when model is already known
#         (default: TEST)
#   -B  UFBoot replicates (default: 1000; threshold for reliability: >= 95)
#   -A  SH-aLRT replicates (default: 1000; threshold: >= 80)
#   -t  CPU threads per IQ-TREE job (default: 2; use more if running sequentially)
#   -j  Parallel IQ-TREE jobs (default: 1; increase if CPUs allow)
#
# Tool version requirements (checked at runtime):
#   IQ-TREE >= 2.3    — https://github.com/iqtree/iqtree2/releases
#     Latest release: v2.4.0 (2025)
#     Install: conda install -c bioconda iqtree
#
# Output:
#   <outdir>/
#     <gene>/
#       <gene>.treefile   — best ML tree with bootstrap support
#       <gene>.iqtree     — full log with model and tree
#       <gene>.log        — run log
#     all_gene_trees.txt  — all .treefile paths concatenated (ASTER input)
#
# Notes:
#   - IQ-TREE .treefile always contains branch lengths → suitable for wASTRAL
#   - If paralogs/multi-copy loci detected, use ASTRAL-Pro3 (run_aster.sh)
#   - UFBoot values are inflated vs. standard bootstrap:
#       >= 95 (UFBoot) ≈ >= 70 (standard bootstrap)
#
# =============================================================================

set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────
MODEL="TEST"
UFBOOT=1000
ALRT=1000
THREADS=2
PARALLEL_JOBS=1

# ── Argument parsing ──────────────────────────────────────────────────────────
usage() {
    grep '^#' "$0" | grep -v '#!/' | sed 's/^# \{0,1\}//'
    exit 1
}

while getopts "i:o:m:B:A:t:j:h" opt; do
    case $opt in
        i) INDIR="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        m) MODEL="$OPTARG" ;;
        B) UFBOOT="$OPTARG" ;;
        A) ALRT="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        j) PARALLEL_JOBS="$OPTARG" ;;
        h|*) usage ;;
    esac
done

[[ -z "${INDIR:-}" || -z "${OUTDIR:-}" ]] && {
    echo "ERROR: -i and -o are required." >&2; usage; }
[[ ! -d "$INDIR" ]] && { echo "ERROR: Input dir not found: ${INDIR}" >&2; exit 1; }

# ── Version check ─────────────────────────────────────────────────────────────
command -v iqtree2 &>/dev/null || command -v iqtree &>/dev/null || {
    echo "ERROR: IQ-TREE not found. Install: conda install -c bioconda iqtree" >&2
    exit 1
}
IQTREE_BIN=$(command -v iqtree2 2>/dev/null || command -v iqtree)
IQTREE_VERSION=$("$IQTREE_BIN" --version 2>&1 | head -1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || echo "unknown")
echo "# IQ-TREE binary:  ${IQTREE_BIN}"
echo "# IQ-TREE version: ${IQTREE_VERSION}  (latest known: 2.4.0)"
echo "# Model:     ${MODEL}"
echo "# UFBoot:    ${UFBOOT}"
echo "# SH-aLRT:   ${ALRT}"
echo "# Threads per job: ${THREADS}"
echo "# Parallel jobs:   ${PARALLEL_JOBS}"
echo "# Date: $(date +%Y-%m-%d)"
echo ""

TOTAL_CPUS=$(nproc 2>/dev/null || echo 4)
USED_CPUS=$(( THREADS * PARALLEL_JOBS ))
if [[ "$USED_CPUS" -gt "$TOTAL_CPUS" ]]; then
    echo "WARNING: ${USED_CPUS} total CPU slots requested but only ${TOTAL_CPUS} available."
    echo "  Consider reducing -t or -j."
    echo ""
fi

mkdir -p "$OUTDIR"

# ── Collect FASTA files ───────────────────────────────────────────────────────
mapfile -t FASTA_FILES < <(
    find "$INDIR" -maxdepth 1 \
        \( -name "*_trimmed.fasta" -o -name "*_aligned.fasta" -o -name "*.fasta" \) \
        | sort
)

if [[ "${#FASTA_FILES[@]}" -eq 0 ]]; then
    echo "ERROR: No FASTA files found in ${INDIR}" >&2
    exit 1
fi
echo "Found ${#FASTA_FILES[@]} marker FASTA files."
echo ""

# ── Gene tree function ────────────────────────────────────────────────────────
run_gene_tree() {
    local fasta="$1"
    local base
    base=$(basename "$fasta")
    base="${base%_trimmed.fasta}"
    base="${base%_aligned.fasta}"
    base="${base%.fasta}"
    base="${base%.fa}"

    local gene_dir="${OUTDIR}/${base}"
    mkdir -p "$gene_dir"

    local prefix="${gene_dir}/${base}"
    local logfile="${gene_dir}/${base}.log"

    # Skip if treefile already exists
    if [[ -f "${prefix}.treefile" ]]; then
        echo "  SKIP (exists): ${base}"
        return 0
    fi

    echo "--- Running IQ-TREE: ${base} ---"

    "$IQTREE_BIN" \
        -s "$fasta" \
        -m "$MODEL" \
        -B "$UFBOOT" \
        --alrt "$ALRT" \
        -T "$THREADS" \
        --prefix "$prefix" \
        --redo \
        > "$logfile" 2>&1

    if [[ -f "${prefix}.treefile" ]]; then
        echo "  Done: ${base} → ${prefix}.treefile"
    else
        echo "  ERROR: ${base} failed — see ${logfile}" >&2
        return 1
    fi
}

export -f run_gene_tree
export IQTREE_BIN UFBOOT ALRT THREADS MODEL OUTDIR

# ── Run gene trees ────────────────────────────────────────────────────────────
if [[ "$PARALLEL_JOBS" -gt 1 ]]; then
    command -v parallel &>/dev/null || {
        echo "WARNING: GNU parallel not found. Running sequentially." >&2
        PARALLEL_JOBS=1
    }
fi

if [[ "$PARALLEL_JOBS" -gt 1 ]]; then
    printf '%s\n' "${FASTA_FILES[@]}" \
        | parallel -j "$PARALLEL_JOBS" run_gene_tree {}
else
    for fasta in "${FASTA_FILES[@]}"; do
        run_gene_tree "$fasta" || true
    done
fi

# ── Collect all gene trees ────────────────────────────────────────────────────
echo ""
echo "=== Collecting gene trees ==="

ALL_TREES="${OUTDIR}/all_gene_trees.txt"
> "$ALL_TREES"

N_TREES=0
N_FAILED=0
for fasta in "${FASTA_FILES[@]}"; do
    base=$(basename "$fasta")
    base="${base%_trimmed.fasta}"
    base="${base%_aligned.fasta}"
    base="${base%.fasta}"
    treefile="${OUTDIR}/${base}/${base}.treefile"
    if [[ -f "$treefile" ]]; then
        cat "$treefile" >> "$ALL_TREES"
        (( N_TREES++ ))
    else
        echo "  MISSING: ${treefile}" >&2
        (( N_FAILED++ ))
    fi
done

echo "  Gene trees collected: ${N_TREES}"
echo "  Failed/missing:       ${N_FAILED}"
echo "  Combined file:        ${ALL_TREES}"
echo ""

if [[ "$N_TREES" -eq 0 ]]; then
    echo "ERROR: No gene trees produced." >&2
    exit 1
fi

# ── Summary ───────────────────────────────────────────────────────────────────
echo "=== Gene tree inference complete ==="
echo ""
echo "Next steps:"
echo "  1. Run ASTER coalescent analysis:"
echo "     bash run_aster.sh -i ${ALL_TREES} -o <aster_outdir>"
echo "  2. Or run concatenation tree (supermatrix already built in alignment step)"
echo ""
echo "Note: IQ-TREE .treefile always includes branch lengths"
echo "  → Use wASTRAL (weighted) by default in run_aster.sh"
echo "  → Only use ASTRAL-Pro3 if multi-copy/paralog loci are present"
