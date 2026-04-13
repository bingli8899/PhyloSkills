#!/usr/bin/env bash
# =============================================================================
# download_sra.sh — Storage-aware SRA download (bulk or streaming mode)
# =============================================================================
# Usage:
#   bash download_sra.sh -l <sra_list.txt> -o <outdir> [-a <assemble_script>] \
#                        [-t <data_type>] [-m <mode>]
#
# Arguments:
#   -l  Text file with one SRA accession per line (SRR/ERR/DRR IDs)
#   -o  Output directory for raw reads (created if absent)
#   -a  Optional: path to assembly script to call in streaming mode
#         If provided, streaming mode calls this script per sample after download
#   -t  Data type hint: wgs | genome_skim | hybseq | transcriptome (default: wgs)
#   -m  Override mode: bulk | streaming | auto (default: auto)
#   -c  Number of CPU threads for fasterq-dump (default: 4)
#   -e  Buffer size for fasterq-dump (default: 250MB)
#
# Tool version requirements (checked at runtime):
#   SRA Toolkit >= 3.1  — https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit
#     Provides: prefetch, fasterq-dump, vdb-dump
#   Entrez Direct >= 21.0 (optional, for size estimation)
#
# Storage logic:
#   - Estimates total dataset size from SRA metadata before any download
#   - If available_storage > estimated_size × 1.5  → bulk mode
#   - Otherwise                                    → streaming mode
#   - Streaming mode: download → assemble → verify → delete raw → next sample
#   - Override with -m bulk|streaming to skip auto-detection
#
# =============================================================================

set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────
DATA_TYPE="wgs"
MODE="auto"
THREADS=4
BUFFER="250MB"
ASSEMBLE_SCRIPT=""

# ── Argument parsing ──────────────────────────────────────────────────────────
usage() {
    grep '^#' "$0" | grep -v '#!/' | sed 's/^# \{0,1\}//'
    exit 1
}

while getopts "l:o:a:t:m:c:e:h" opt; do
    case $opt in
        l) SRA_LIST="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        a) ASSEMBLE_SCRIPT="$OPTARG" ;;
        t) DATA_TYPE="$OPTARG" ;;
        m) MODE="$OPTARG" ;;
        c) THREADS="$OPTARG" ;;
        e) BUFFER="$OPTARG" ;;
        h|*) usage ;;
    esac
done

[[ -z "${SRA_LIST:-}" || -z "${OUTDIR:-}" ]] && {
    echo "ERROR: -l and -o are required." >&2; usage; }
[[ ! -f "$SRA_LIST" ]] && {
    echo "ERROR: SRA list file not found: ${SRA_LIST}" >&2; exit 1; }

# ── Version checks ────────────────────────────────────────────────────────────
check_tool() {
    command -v "$1" &>/dev/null || {
        echo "ERROR: '$1' not found. Install SRA Toolkit:" >&2
        echo "  https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit" >&2
        exit 1
    }
}
check_tool prefetch
check_tool fasterq-dump
check_tool vdb-dump

SRATOOLS_VERSION=$(fasterq-dump --version 2>&1 | head -1 || echo "unknown")
echo "# SRA Toolkit version: ${SRATOOLS_VERSION}"
echo "# SRA list: ${SRA_LIST}"
echo "# Output dir: ${OUTDIR}"
echo "# Date: $(date +%Y-%m-%d)"
echo ""

mkdir -p "${OUTDIR}"

# ── Storage estimation ────────────────────────────────────────────────────────
estimate_sra_size_mb() {
    local accession="$1"
    # Try vdb-dump first (fast, no network beyond metadata)
    local size_mb
    size_mb=$(vdb-dump --info "$accession" 2>/dev/null \
        | grep -i "file size" | grep -oE '[0-9]+' | head -1 || echo "0")
    # Convert bytes to MB if needed (vdb-dump reports bytes)
    if [[ "$size_mb" -gt 100000 ]]; then
        size_mb=$(( size_mb / 1048576 ))
    fi
    echo "$size_mb"
}

available_storage_mb() {
    # df output: use the directory where data will be stored
    df -m "${OUTDIR}" | awk 'NR==2 {print $4}'
}

echo "=== Storage Assessment ==="
AVAIL_MB=$(available_storage_mb)
echo "  Available storage: ${AVAIL_MB} MB ($(( AVAIL_MB / 1024 )) GB)"

TOTAL_EST_MB=0
declare -A SRA_SIZES

while read -r acc; do
    [[ -z "$acc" || "$acc" =~ ^# ]] && continue
    size_mb=$(estimate_sra_size_mb "$acc")
    SRA_SIZES["$acc"]=$size_mb
    TOTAL_EST_MB=$(( TOTAL_EST_MB + size_mb ))
    echo "  ${acc}: ~${size_mb} MB"
done < "$SRA_LIST"

SAFE_THRESHOLD=$(( TOTAL_EST_MB * 3 / 2 ))  # × 1.5
echo ""
echo "  Estimated total raw data: ${TOTAL_EST_MB} MB"
echo "  Safety threshold (×1.5): ${SAFE_THRESHOLD} MB"

if [[ "$MODE" == "auto" ]]; then
    if [[ "$AVAIL_MB" -gt "$SAFE_THRESHOLD" ]]; then
        MODE="bulk"
    else
        MODE="streaming"
    fi
fi
echo "  Selected mode: ${MODE}"
echo ""

# ── Log file ──────────────────────────────────────────────────────────────────
LOG_FILE="${OUTDIR}/download_log_$(date +%Y-%m-%d).tsv"
echo -e "accession\tstatus\traw_deleted\tnotes" > "$LOG_FILE"

# ── Download function (single sample) ─────────────────────────────────────────
download_sample() {
    local acc="$1"
    echo "--- Downloading: ${acc} ---"

    local raw_sra_dir="${OUTDIR}/${acc}"
    local read1="${OUTDIR}/${acc}_1.fastq"
    local read2="${OUTDIR}/${acc}_2.fastq"

    # prefetch: download .sra file
    prefetch "$acc" -O "${OUTDIR}/" --max-size 100G || {
        echo "  ERROR: prefetch failed for ${acc}" >&2
        echo -e "${acc}\tFAILED_PREFETCH\tno\tprefetch error" >> "$LOG_FILE"
        return 1
    }

    # fasterq-dump: extract to FASTQ
    fasterq-dump "${raw_sra_dir}/" \
        -O "${OUTDIR}/" \
        --split-files \
        --threads "$THREADS" \
        --bufsize "$BUFFER" || {
        echo "  ERROR: fasterq-dump failed for ${acc}" >&2
        echo -e "${acc}\tFAILED_FASTERQDUMP\tno\tfasterq-dump error" >> "$LOG_FILE"
        return 1
    }

    echo "  Downloaded: ${read1}, ${read2}"
}

# ── Delete raw reads (streaming mode only) ────────────────────────────────────
delete_raw() {
    local acc="$1"
    local assembly_path="$2"

    # Verify assembly output exists and is non-empty before deleting
    if [[ -n "$assembly_path" && -e "$assembly_path" && -s "$assembly_path" ]]; then
        rm -rf "${OUTDIR}/${acc}/"
        rm -f "${OUTDIR}/${acc}_1.fastq" "${OUTDIR}/${acc}_2.fastq"
        echo "  Raw reads deleted: ${acc}"
        echo -e "${acc}\tCOMPLETE\tyes\t" >> "$LOG_FILE"
    elif [[ -z "$assembly_path" ]]; then
        # No assembly script — streaming mode without assembly step
        # Do not delete; log for manual review
        echo "  WARNING: No assembly path specified — raw reads retained: ${acc}" >&2
        echo -e "${acc}\tDOWNLOAD_ONLY\tno\tno assembly script" >> "$LOG_FILE"
    else
        echo "  WARNING: Assembly output not found at: ${assembly_path}" >&2
        echo "  Raw reads RETAINED for debugging: ${acc}" >&2
        echo -e "${acc}\tASSEMBLY_FAILED\tno\toutput not found: ${assembly_path}" >> "$LOG_FILE"
    fi
}

# ── Bulk mode ─────────────────────────────────────────────────────────────────
run_bulk() {
    echo "=== BULK MODE: downloading all samples first ==="
    while read -r acc; do
        [[ -z "$acc" || "$acc" =~ ^# ]] && continue
        download_sample "$acc" && \
            echo -e "${acc}\tDOWNLOAD_COMPLETE\tno\t" >> "$LOG_FILE"
    done < "$SRA_LIST"
    echo ""
    echo "All downloads complete. Proceed to phylo-assemble for assembly."
}

# ── Streaming mode ────────────────────────────────────────────────────────────
run_streaming() {
    echo "=== STREAMING MODE: download → assemble → delete raw → next ==="
    while read -r acc; do
        [[ -z "$acc" || "$acc" =~ ^# ]] && continue

        download_sample "$acc" || continue

        if [[ -n "$ASSEMBLE_SCRIPT" && -x "$ASSEMBLE_SCRIPT" ]]; then
            echo "  Running assembly: ${ASSEMBLE_SCRIPT} for ${acc}"
            assembly_output=""
            # Call assembly script; it should print output path to stdout on last line
            assembly_output=$(bash "$ASSEMBLE_SCRIPT" "$acc" \
                "${OUTDIR}/${acc}_1.fastq" "${OUTDIR}/${acc}_2.fastq" \
                2>&1 | tee /dev/stderr | tail -1 || true)
            delete_raw "$acc" "$assembly_output"
        else
            echo "  No assembly script provided; raw reads retained: ${acc}"
            echo -e "${acc}\tDOWNLOAD_ONLY\tno\tno assembly script" >> "$LOG_FILE"
        fi

        echo "  Storage after ${acc}: $(available_storage_mb) MB available"
        echo ""
    done < "$SRA_LIST"
}

# ── Run ───────────────────────────────────────────────────────────────────────
if [[ "$MODE" == "bulk" ]]; then
    run_bulk
else
    run_streaming
fi

echo "=== Download session complete ==="
echo "Log: ${LOG_FILE}"
echo ""
echo "Summary:"
awk -F'\t' 'NR>1 {counts[$2]++} END {for(s in counts) print "  " s ": " counts[s]}' "$LOG_FILE"
