#!/usr/bin/env bash
# =============================================================================
# run_aster.sh — Coalescent species tree with ASTER (wASTRAL or ASTRAL-Pro3)
# =============================================================================
# Auto-detects the appropriate ASTER tool based on dataset type:
#   - Multi-copy/paralog genes  → ASTRAL-Pro3
#   - IQ-TREE gene trees (default, branch lengths present) → wASTRAL
#   - No branch lengths         → standard ASTRAL
#
# Usage:
#   bash run_aster.sh \
#     -i <all_gene_trees.txt> \
#     -o <outdir> \
#     [-T <tool>] \
#     [-t <threads>] \
#     [-u <weight_mode>] \
#     [-c <concat_tree>] \
#     [-a <aster_bin_dir>]
#
# Arguments:
#   -i  Input gene trees file (one Newick tree per line)
#   -o  Output directory
#   -T  Tool override: wastral | astral-pro3 | astral | auto (default: auto)
#         auto checks for branch lengths in the gene tree file:
#           branch lengths present → wastral
#           no branch lengths      → astral
#         Use -T astral-pro3 explicitly if multi-copy loci are present
#   -t  CPU threads (default: auto-detect via nproc)
#   -u  wASTRAL weight mode (default: 2 = use both branch lengths + local posterior)
#         0 = ASTRAL mode (no weights)
#         1 = branch lengths only
#         2 = branch lengths + local posterior (recommended for IQ-TREE output)
#   -c  Concatenation tree FASTA or Newick file for concordance factor computation
#         If provided, also runs IQ-TREE concordance factors after ASTER
#   -a  Directory containing ASTER binaries (if not in PATH or executables/)
#
# Tool version requirements (checked at runtime):
#   ASTER v1.23  — https://github.com/chaoszhang/ASTER
#     Latest release: v1.23 (2025)
#     Provides: wastral, astral-pro3, astral (as separate binaries in bin/)
#     Install: download release zip, extract, add bin/ to PATH
#              or place in executables/aster/
#   IQ-TREE >= 2.3 (for concordance factors, optional)
#
# Output:
#   <outdir>/
#     species_tree_wastral.nwk     — wASTRAL species tree
#     species_tree_astral_pro3.nwk — ASTRAL-Pro3 species tree
#     concordance/                  — concordance factor output (if -c provided)
#
# =============================================================================

set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────
TOOL="auto"
THREADS=$(nproc 2>/dev/null || echo 8)
WEIGHT_MODE=2
CONCAT_TREE=""
ASTER_BIN_DIR=""

# ── Argument parsing ──────────────────────────────────────────────────────────
usage() {
    grep '^#' "$0" | grep -v '#!/' | sed 's/^# \{0,1\}//'
    exit 1
}

while getopts "i:o:T:t:u:c:a:h" opt; do
    case $opt in
        i) GENE_TREES="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        T) TOOL="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        u) WEIGHT_MODE="$OPTARG" ;;
        c) CONCAT_TREE="$OPTARG" ;;
        a) ASTER_BIN_DIR="$OPTARG" ;;
        h|*) usage ;;
    esac
done

[[ -z "${GENE_TREES:-}" || -z "${OUTDIR:-}" ]] && {
    echo "ERROR: -i and -o are required." >&2; usage; }
[[ ! -f "$GENE_TREES" ]] && {
    echo "ERROR: Gene trees file not found: ${GENE_TREES}" >&2; exit 1; }

# ── Tool detection ────────────────────────────────────────────────────────────
find_aster_bin() {
    local name="$1"
    # 1. Explicit bin dir
    [[ -n "$ASTER_BIN_DIR" && -x "${ASTER_BIN_DIR}/${name}" ]] && {
        echo "${ASTER_BIN_DIR}/${name}"; return 0; }
    # 2. PATH
    command -v "$name" &>/dev/null && { echo "$(command -v "$name")"; return 0; }
    # 3. executables/aster/bin/
    local exe_path
    exe_path=$(find executables/aster/ -name "$name" -type f 2>/dev/null | head -1 || true)
    [[ -n "$exe_path" ]] && { echo "$exe_path"; return 0; }
    echo ""
}

WASTRAL_BIN=$(find_aster_bin "wastral")
ASTRAL_PRO3_BIN=$(find_aster_bin "astral-pro3")
ASTRAL_BIN=$(find_aster_bin "astral")

echo "# ASTER tool scan:"
echo "#   wastral:      ${WASTRAL_BIN:-not found}"
echo "#   astral-pro3:  ${ASTRAL_PRO3_BIN:-not found}"
echo "#   astral:       ${ASTRAL_BIN:-not found}"
echo "# Latest known ASTER release: v1.23 (2025)"
echo "#   Install: https://github.com/chaoszhang/ASTER/releases"
echo ""

# ── Auto-detect tool from gene tree file ─────────────────────────────────────
detect_tool_from_trees() {
    # Check if trees have branch lengths (colon followed by number)
    if grep -qE ':[0-9]' "$GENE_TREES" 2>/dev/null; then
        echo "wastral"
    else
        echo "astral"
    fi
}

if [[ "$TOOL" == "auto" ]]; then
    TOOL=$(detect_tool_from_trees)
    echo "Auto-detected tool: ${TOOL}"
    echo "  (based on $(grep -c '' "$GENE_TREES") gene trees in ${GENE_TREES})"
fi

# Validate tool availability
case "$TOOL" in
    wastral)
        [[ -z "$WASTRAL_BIN" ]] && {
            echo "ERROR: wastral not found. Download ASTER release:" >&2
            echo "  https://github.com/chaoszhang/ASTER/releases" >&2
            echo "  Extract and add bin/ to PATH, or use -a <path>" >&2
            exit 1
        }
        ;;
    astral-pro3)
        [[ -z "$ASTRAL_PRO3_BIN" ]] && {
            echo "ERROR: astral-pro3 not found." >&2
            echo "  Download ASTER release: https://github.com/chaoszhang/ASTER/releases" >&2
            exit 1
        }
        ;;
    astral)
        [[ -z "$ASTRAL_BIN" ]] && {
            echo "ERROR: astral not found." >&2
            exit 1
        }
        ;;
    *)
        echo "ERROR: Unknown tool: ${TOOL}. Use: wastral | astral-pro3 | astral | auto" >&2
        exit 1
        ;;
esac

N_TREES=$(grep -c '' "$GENE_TREES" || echo 0)
echo "# Tool selected: ${TOOL}"
echo "# Gene trees: ${N_TREES}"
echo "# Threads: ${THREADS}"
echo "# Date: $(date +%Y-%m-%d)"
echo ""

mkdir -p "$OUTDIR"

# ── Run ASTER ─────────────────────────────────────────────────────────────────
case "$TOOL" in

    wastral)
        SPECIES_TREE="${OUTDIR}/species_tree_wastral.nwk"
        echo "=== Running wASTRAL (weighted coalescent) ==="
        echo "  Weight mode: -u ${WEIGHT_MODE}"
        echo "  (0=no weights, 1=branch lengths, 2=branch lengths + local posterior)"
        echo ""
        "$WASTRAL_BIN" \
            --input "$GENE_TREES" \
            --output "$SPECIES_TREE" \
            -u "$WEIGHT_MODE" \
            --thread "$THREADS" \
            2>&1 | tee "${OUTDIR}/wastral.log"
        ;;

    astral-pro3)
        SPECIES_TREE="${OUTDIR}/species_tree_astral_pro3.nwk"
        echo "=== Running ASTRAL-Pro3 (multi-copy/paralog-aware) ==="
        echo ""
        "$ASTRAL_PRO3_BIN" \
            --input "$GENE_TREES" \
            --output "$SPECIES_TREE" \
            --thread "$THREADS" \
            2>&1 | tee "${OUTDIR}/astral_pro3.log"
        ;;

    astral)
        SPECIES_TREE="${OUTDIR}/species_tree_astral.nwk"
        echo "=== Running ASTRAL (standard, no branch length weights) ==="
        echo ""
        "$ASTRAL_BIN" \
            --input "$GENE_TREES" \
            --output "$SPECIES_TREE" \
            --thread "$THREADS" \
            2>&1 | tee "${OUTDIR}/astral.log"
        ;;
esac

# ── Verify output ─────────────────────────────────────────────────────────────
echo ""
echo "=== ASTER output check ==="
if [[ -f "$SPECIES_TREE" && -s "$SPECIES_TREE" ]]; then
    echo "  Species tree: ${SPECIES_TREE}"
    echo "  Content preview: $(head -c 200 "$SPECIES_TREE")"
else
    echo "ERROR: Species tree not generated." >&2
    exit 1
fi

# ── Concordance factors (optional) ───────────────────────────────────────────
if [[ -n "$CONCAT_TREE" ]]; then
    echo ""
    echo "=== Computing concordance factors ==="

    IQTREE_BIN=$(command -v iqtree2 2>/dev/null || command -v iqtree || true)
    if [[ -z "$IQTREE_BIN" ]]; then
        echo "WARNING: IQ-TREE not found — skipping concordance factors." >&2
    else
        IQTREE_VERSION=$("$IQTREE_BIN" --version 2>&1 | head -1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || echo "unknown")
        echo "# IQ-TREE version: ${IQTREE_VERSION}  (latest known: 2.4.0)"

        CF_DIR="${OUTDIR}/concordance"
        mkdir -p "$CF_DIR"

        "$IQTREE_BIN" \
            -t "$SPECIES_TREE" \
            --gcf "$GENE_TREES" \
            -p "$CONCAT_TREE" \
            --scf 100 \
            --prefix "${CF_DIR}/concordance" \
            -T "$THREADS" \
            2>&1 | tee "${CF_DIR}/concordance.log"

        echo ""
        echo "  Concordance factors: ${CF_DIR}/"
        echo "  gCF = gene concordance factor (% gene trees supporting branch)"
        echo "  sCF = site concordance factor (% informative sites supporting branch)"
    fi
fi

# ── Interpretation guide ─────────────────────────────────────────────────────
echo ""
echo "=== Interpreting ASTER output ==="
echo "  Local posterior support values on branches: range 0–1"
echo "  Well-supported: >= 0.95"
echo "  Compare this species tree to the concatenation tree:"
echo "    - Concordant topology → ILS unlikely to affect conclusions"
echo "    - Discordant topology → ILS present; species tree is more appropriate"
echo ""
echo "Output: ${SPECIES_TREE}"
