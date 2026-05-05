#!/usr/bin/env bash
# =============================================================================
# align_markers.sh — Per-marker alignment with MAFFT + optional trimAl
# =============================================================================
# Aligns each marker FASTA independently, optionally trims, then optionally
# concatenates into a supermatrix with partition file using AMAS.
#
# Usage:
#   bash align_markers.sh \
#     -i <input_dir> \
#     -o <outdir> \
#     [-s <strategy>] \
#     [-T] \
#     [-t <threads>] \
#     [-c] \
#     [-m <trim_mode>]
#
# Arguments:
#   -i  Input directory containing per-marker FASTA files (*.fasta)
#   -o  Output directory for aligned files
#   -s  MAFFT strategy: auto | linsi | ginsi | einsi | fftnsi | fftns2
#         auto    = MAFFT chooses internally (safe default)
#         linsi   = iterative local alignment; best for ≤200 seqs, moderate divergence
#         ginsi   = global/local pairwise; high divergence, gappy datasets
#         einsi   = E-INS-i; multiple conserved domains with large unalignable regions
#         fftnsi  = FFT-based; fast for 200–10,000 sequences
#         fftns2  = fastest; intraspecific/population level data
#         (default: auto)
#   -T  Run trimAl after alignment (flag; default: off)
#   -m  trimAl mode: automated1 | gappyout (default: automated1)
#   -t  CPU threads (default: auto-detect via nproc)
#   -c  Concatenate all trimmed/aligned markers with AMAS (flag; default: off)
#   -A  Path to AMAS.py (default: auto-detect in PATH; use when AMAS is not in PATH)
#   -d  Data type for AMAS: dna | aa (default: dna)
#   -f  Input format for AMAS: fasta | phylip | nexus (default: fasta)
#
# Tool version requirements (checked at runtime):
#   MAFFT >= 7.5      — https://mafft.cbrc.jp/alignment/software/
#     Latest: 7.526 (2024)
#     Install: conda install -c bioconda mafft
#   trimAl >= 1.5     — https://github.com/inab/trimal (latest: v1.5.1, 2024)
#     Install: conda install -c bioconda trimal
#   AMAS (for -c)     — https://github.com/marekborowiec/AMAS
#     Install: pip install amas  OR  conda install -c bioconda amas
#
# Output:
#   <outdir>/
#     <marker>_aligned.fasta         — per-marker alignment
#     <marker>_trimmed.fasta         — after trimAl (if -T)
#     concatenated.fasta             — supermatrix (if -c)
#     partition.txt                  — RAxML-style partition file (if -c)
#     alignment_stats.tsv            — length, gap%, parsimony-informative sites
#
# =============================================================================

set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────
STRATEGY="auto"
TRIM=false
TRIM_MODE="automated1"
THREADS=$(nproc 2>/dev/null || echo 4)
CONCATENATE=false
DATA_TYPE="dna"
FORMAT="fasta"
AMAS_PATH=""

# ── Argument parsing ──────────────────────────────────────────────────────────
usage() {
    grep '^#' "$0" | grep -v '#!/' | sed 's/^# \{0,1\}//'
    exit 1
}

PROVENANCE_DIR=""

while getopts "i:o:s:Tm:t:cA:d:f:P:h" opt; do
    case $opt in
        i) INDIR="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        s) STRATEGY="$OPTARG" ;;
        T) TRIM=true ;;
        m) TRIM_MODE="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        c) CONCATENATE=true ;;
        A) AMAS_PATH="$OPTARG" ;;
        d) DATA_TYPE="$OPTARG" ;;
        f) FORMAT="$OPTARG" ;;
        P) PROVENANCE_DIR="$OPTARG" ;;
        h|*) usage ;;
    esac
done

[[ -z "${INDIR:-}" || -z "${OUTDIR:-}" ]] && {
    echo "ERROR: -i and -o are required." >&2; usage; }
[[ ! -d "$INDIR" ]] && { echo "ERROR: Input dir not found: ${INDIR}" >&2; exit 1; }

# ── Version checks ────────────────────────────────────────────────────────────
command -v mafft &>/dev/null || {
    echo "ERROR: MAFFT not found. Install: conda install -c bioconda mafft" >&2; exit 1; }

MAFFT_VERSION=$(mafft --version 2>&1 | head -1 | awk '{print $1}' || echo "unknown")
echo "# MAFFT version: ${MAFFT_VERSION}  (latest known: 7.526)"
echo "# Available CPUs: $(nproc 2>/dev/null || echo unknown)"
echo "# Using threads: ${THREADS}"
echo "# Strategy: ${STRATEGY}"
echo "# Trimming: ${TRIM} (mode: ${TRIM_MODE})"
echo "# Concatenate: ${CONCATENATE}"
echo "# Date: $(date +%Y-%m-%d)"
echo ""

if [[ "$TRIM" == true ]]; then
    command -v trimal &>/dev/null || {
        echo "ERROR: trimAl not found. Install: conda install -c bioconda trimal" >&2; exit 1; }
    TRIMAL_VERSION=$(trimal --version 2>&1 | head -1 || echo "unknown")
    echo "# trimAl version: ${TRIMAL_VERSION}  (latest known: 1.5.1)"
fi

if [[ "$CONCATENATE" == true ]]; then
    if [[ -n "$AMAS_PATH" ]]; then
        # Expand ~ in path
        AMAS_PATH="${AMAS_PATH/#\~/$HOME}"
        [[ ! -f "$AMAS_PATH" ]] && {
            echo "ERROR: AMAS not found at: ${AMAS_PATH}" >&2; exit 1; }
    else
        command -v AMAS.py &>/dev/null || command -v amas &>/dev/null || {
            echo "ERROR: AMAS not found. Install: pip install amas  OR  use -A <path/to/AMAS.py>" >&2; exit 1; }
    fi
fi

mkdir -p "$OUTDIR"

# ── Stats file ────────────────────────────────────────────────────────────────
STATS_FILE="${OUTDIR}/alignment_stats.tsv"
echo -e "marker\tn_seqs\taligned_len\tgap_pct\ttrimmed_len\tparsimony_informative_sites\tpis_pct\tnotes" \
    > "$STATS_FILE"

# ── MAFFT strategy flag ───────────────────────────────────────────────────────
mafft_strategy_flag() {
    local s="$1" n_seqs="$2"
    case "$s" in
        linsi)   echo "--localpair --maxiterate 1000" ;;
        ginsi)   echo "--globalpair --maxiterate 1000" ;;
        einsi)   echo "--ep 0 --genafpair --maxiterate 1000" ;;
        fftnsi)  echo "--retree 2 --maxiterate 2" ;;
        fftns2)  echo "--retree 2 --maxiterate 0" ;;
        auto|*)  echo "--auto" ;;
    esac
}

# ── Compute alignment stats ───────────────────────────────────────────────────
compute_stats() {
    local fasta="$1"
    python3 - "$fasta" <<'PYEOF'
import sys, re
from collections import Counter

fasta = sys.argv[1]
seqs = {}
header = ""
with open(fasta) as f:
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            header = line[1:].split()[0]
            seqs[header] = ""
        else:
            seqs[header] += line.upper()

if not seqs:
    print("0\t0\t0\t0\t0")
    sys.exit(0)

n_seqs = len(seqs)
lengths = [len(s) for s in seqs.values()]
aln_len = max(lengths) if lengths else 0

# Gap percentage across all positions
total_chars = sum(lengths)
gap_count = sum(s.count("-") + s.count("?") + s.count("N") for s in seqs.values())
gap_pct = round(gap_count / total_chars * 100, 1) if total_chars > 0 else 0

# Parsimony informative sites
pis = 0
seqs_list = list(seqs.values())
for i in range(aln_len):
    col = [s[i] for s in seqs_list if i < len(s)]
    col = [c for c in col if c not in ("-", "?", "N")]
    counts = Counter(col)
    n_variable = sum(1 for c, n in counts.items() if n >= 2)
    if n_variable >= 2:
        pis += 1

pis_pct = round(pis / aln_len * 100, 1) if aln_len > 0 else 0
print(f"{n_seqs}\t{aln_len}\t{gap_pct}\t{pis}\t{pis_pct}")
PYEOF
}

# ── Main alignment loop ───────────────────────────────────────────────────────
ALIGNED_FILES=()
N_SUCCESS=0
N_SKIP=0

for fasta in "${INDIR}"/*.fasta "${INDIR}"/*.fa; do
    [[ ! -f "$fasta" ]] && continue

    base=$(basename "$fasta")
    base="${base%.fasta}"
    base="${base%.fa}"

    # Count sequences
    n_seqs=$(grep -c '^>' "$fasta" 2>/dev/null || echo 0)
    if [[ "$n_seqs" -lt 2 ]]; then
        echo "SKIP ${base}: only ${n_seqs} sequence(s) — need >= 2 to align"
        (( N_SKIP++ ))
        continue
    fi

    echo "=== Aligning: ${base} (${n_seqs} sequences) ==="

    aligned="${OUTDIR}/${base}_aligned.fasta"
    STRAT_FLAGS=$(mafft_strategy_flag "$STRATEGY" "$n_seqs")

    mafft $STRAT_FLAGS \
        --thread "$THREADS" \
        --quiet \
        "$fasta" \
        > "$aligned"

    echo "  Aligned: ${aligned}"

    # Trim if requested
    FINAL_ALIGNED="$aligned"
    trimmed_len="NA"
    if [[ "$TRIM" == true ]]; then
        trimmed="${OUTDIR}/${base}_trimmed.fasta"
        if [[ "$TRIM_MODE" == "automated1" ]]; then
            trimal -in "$aligned" -out "$trimmed" -automated1
        else
            trimal -in "$aligned" -out "$trimmed" -gappyout
        fi
        trimmed_len=$(awk '/^>/{next} {len+=length($0)} END {print len+0}' "$trimmed")
        echo "  Trimmed: ${trimmed} (${trimmed_len} bp)"
        FINAL_ALIGNED="$trimmed"
    fi

    # Stats
    stats=$(compute_stats "$FINAL_ALIGNED")
    read -r n_seqs_aln aln_len gap_pct pis pis_pct <<< "$stats"

    # Warnings
    notes=""
    if (( $(echo "$gap_pct > 30" | awk '{print ($1>30)}') )) 2>/dev/null; then
        notes="HIGH_GAP"
        echo "  WARNING: Gap% = ${gap_pct}% (threshold: 30%)"
    fi
    if (( $(echo "$pis_pct < 5" | awk '{print ($1<5)}') )) 2>/dev/null; then
        [[ -n "$notes" ]] && notes="${notes},LOW_PIS" || notes="LOW_PIS"
        echo "  WARNING: Parsimony-informative sites = ${pis_pct}% (threshold: 5%)"
    fi

    echo -e "${base}\t${n_seqs_aln}\t${aln_len}\t${gap_pct}\t${trimmed_len}\t${pis}\t${pis_pct}\t${notes}" \
        >> "$STATS_FILE"

    ALIGNED_FILES+=("$FINAL_ALIGNED")
    (( N_SUCCESS++ ))
    echo ""
done

echo "=== Alignment complete: ${N_SUCCESS} markers aligned, ${N_SKIP} skipped ==="
echo "Stats: ${STATS_FILE}"
echo ""

# ── Concatenation ─────────────────────────────────────────────────────────────
if [[ "$CONCATENATE" == true && "${#ALIGNED_FILES[@]}" -gt 0 ]]; then
    echo "=== Concatenating ${#ALIGNED_FILES[@]} markers with AMAS ==="

    AMAS_CMD=""
    if [[ -n "$AMAS_PATH" ]]; then
        AMAS_CMD="python3 ${AMAS_PATH}"
    else
        command -v AMAS.py &>/dev/null && AMAS_CMD="AMAS.py" || AMAS_CMD="amas"
    fi

    $AMAS_CMD concat \
        -i "${ALIGNED_FILES[@]}" \
        -f fasta \
        -d "$DATA_TYPE" \
        -t "${OUTDIR}/concatenated.fasta" \
        -p "${OUTDIR}/partition.txt" \
        -u fasta \
        -y raxml

    echo "  Supermatrix: ${OUTDIR}/concatenated.fasta"
    echo "  Partition:   ${OUTDIR}/partition.txt"

    # Quick missing data estimate
    python3 - "${OUTDIR}/concatenated.fasta" <<'PYEOF'
import sys
seqs = {}
header = ""
with open(sys.argv[1]) as f:
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            header = line[1:].split()[0]; seqs[header] = ""
        else:
            seqs[header] += line.upper()
if seqs:
    total = sum(len(s) for s in seqs.values())
    missing = sum(s.count("?") + s.count("N") for s in seqs.values())
    print(f"  Total taxa: {len(seqs)}, Length: {max(len(s) for s in seqs.values())} bp")
    print(f"  Missing data: {missing/total*100:.1f}%")
PYEOF
fi

echo ""
echo "Next step: phylo-model-selection"

# ── Provenance ────────────────────────────────────────────────────────────────
if [[ -n "${PROVENANCE_DIR:-}" ]]; then
    MAFFT_VER=$(mafft --version 2>&1 | head -1 | awk '{print $1}')
    TRIMAL_VER=""
    [[ "$TRIM" == true ]] && TRIMAL_VER=$(trimal --version 2>&1 | head -1 || echo "unknown")
    END_T=$(date +%s)
    mkdir -p "$PROVENANCE_DIR"
    PROV_FILE="${PROVENANCE_DIR}/align_markers_$(date +%Y-%m-%d).json"
    python3 - <<PYEOF
import json
data = {
    "script": "align_markers",
    "tool": "mafft",
    "version": "${MAFFT_VER}",
    "date": "$(date +%Y-%m-%dT%H:%M:%S)",
    "parameters": {
        "strategy": "${STRATEGY}",
        "threads": ${THREADS},
        "trim": ${TRIM},
        "trim_mode": "${TRIM_MODE}",
        "concatenate": ${CONCATENATE},
        "data_type": "${DATA_TYPE}"
    },
    "input_files": {"input_dir": "${INDIR}"},
    "output_files": {
        "alignment_dir": "${OUTDIR}",
        "stats": "${OUTDIR}/alignment_stats.tsv"
    },
    "trimal_version": "${TRIMAL_VER}",
    "runtime_seconds": $(( $(date +%s) - _PROV_START )),
    "exit_code": 0,
    "working_dir": "$(pwd)"
}
import pathlib
pathlib.Path("${PROV_FILE}").parent.mkdir(parents=True, exist_ok=True)
with open("${PROV_FILE}", "w") as f:
    json.dump(data, f, indent=2)
print(f"Provenance written: ${PROV_FILE}")
PYEOF
fi
