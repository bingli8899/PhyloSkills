#!/usr/bin/env bash
# =============================================================================
# download_genbank.sh — Download sequences from GenBank via Entrez Direct
# =============================================================================
# Usage:
#   bash download_genbank.sh -g <group> -m <marker1,marker2,...> -o <outdir> \
#                            [-n <max_seqs>] [-f <format>]
#
# Arguments:
#   -g  Taxonomic group name (e.g. "Zingiberaceae", "Zingiber", "Poales")
#   -m  Comma-separated list of markers (e.g. "matK,rbcL,ITS,psbA-trnH")
#   -o  Output directory (will be created if absent)
#   -n  Maximum sequences per marker to download (default: 500)
#   -f  Output format: fasta or gb (default: fasta)
#   -s  Survey only — print coverage matrix, do not download (flag)
#
# Tool version requirements (checked at runtime):
#   Entrez Direct (edirect) >= 21.0  — https://www.ncbi.nlm.nih.gov/books/NBK179288/
#     Install: sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
#
# Output naming convention:  <Genus>_<species>_<accession>_<marker>.fasta
#   Spaces in organism names are replaced with underscores.
#   If species is not parseable, "sp" is used as placeholder.
#
# =============================================================================

set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────
MAX_SEQS=500
FORMAT="fasta"
SURVEY_ONLY=false

# ── Marker size thresholds (bp) ───────────────────────────────────────────────
# Any downloaded FASTA longer than this limit is assumed to be a complete
# plastome record rather than a single-gene amplicon.  The script will then
# re-fetch the GenBank annotation and extract the gene with extract_gene_from_gb.py.
declare -A MARKER_SIZE_MAX=(
    [matK]=2500
    [rbcL]=1800
    [trnL]=2500
    [psbA-trnH]=1000
    [rpoB]=4000
    [rpoC1]=2500
    [atpB]=2000
    [ndhF]=3000
    [ycf1]=10000
    [ycf2]=10000
)

# Resolve the repo root so extract_gene_from_gb.py can be found regardless of
# the caller's working directory.
PHYLOSKILLS_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

# Return the size threshold for a given marker; fall back to 5000 if unknown.
get_marker_size_max() {
    local m="$1"
    echo "${MARKER_SIZE_MAX[$m]:-5000}"
}

# ── Argument parsing ──────────────────────────────────────────────────────────
usage() {
    grep '^#' "$0" | grep -v '#!/' | sed 's/^# \{0,1\}//'
    exit 1
}

while getopts "g:m:o:n:f:sh" opt; do
    case $opt in
        g) GROUP="$OPTARG" ;;
        m) MARKERS_RAW="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        n) MAX_SEQS="$OPTARG" ;;
        f) FORMAT="$OPTARG" ;;
        s) SURVEY_ONLY=true ;;
        h|*) usage ;;
    esac
done

[[ -z "${GROUP:-}" || -z "${MARKERS_RAW:-}" || -z "${OUTDIR:-}" ]] && {
    echo "ERROR: -g, -m, and -o are required." >&2; usage; }

# ── Version check ─────────────────────────────────────────────────────────────
check_tool() {
    command -v "$1" &>/dev/null || {
        echo "ERROR: '$1' not found. Install Entrez Direct:" >&2
        echo "  sh -c \"\$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)\"" >&2
        exit 1
    }
}
check_tool esearch
check_tool efetch
check_tool xtract

# Record version (edirect does not have a single --version flag; use einfo as proxy)
EDIRECT_VERSION=$(einfo --version 2>/dev/null | head -1 || echo "unknown")
echo "# Entrez Direct version: ${EDIRECT_VERSION}"
echo "# Group: ${GROUP}"
echo "# Markers: ${MARKERS_RAW}"
echo "# Date: $(date +%Y-%m-%d)"
echo ""

# ── Parse markers ─────────────────────────────────────────────────────────────
IFS=',' read -ra MARKERS <<< "$MARKERS_RAW"

mkdir -p "${OUTDIR}"

# ── Survey function ───────────────────────────────────────────────────────────
survey_marker() {
    local marker="$1"
    local count
    count=$(esearch -db nuccore \
        -query "\"${GROUP}\"[Organism] AND (\"${marker}\"[All Fields]) \
                AND biomol_genomic[PROP]" \
        | xtract -pattern ENTREZ_DIRECT -element Count 2>/dev/null || echo "0")
    echo "  ${marker}: ${count} records"
}

echo "=== Coverage Survey: ${GROUP} ==="
for marker in "${MARKERS[@]}"; do
    survey_marker "$marker"
done
echo ""

[[ "$SURVEY_ONLY" == true ]] && { echo "Survey complete (no download)."; exit 0; }

# ── Download function ─────────────────────────────────────────────────────────
download_marker() {
    local marker="$1"
    local marker_dir="${OUTDIR}/${marker}"
    mkdir -p "${marker_dir}"

    echo "=== Downloading: ${marker} ==="

    # Fetch doc summaries to get accession + organism
    local tmpfile
    tmpfile=$(mktemp /tmp/genbank_meta_XXXXXX.tsv)

    esearch -db nuccore \
        -query "\"${GROUP}\"[Organism] AND (\"${marker}\"[All Fields]) \
                AND biomol_genomic[PROP]" \
        | efetch -format docsum \
        | xtract -pattern DocumentSummary \
            -element AccessionVersion Organism \
        | head -n "$MAX_SEQS" > "$tmpfile"

    local n_seqs
    n_seqs=$(wc -l < "$tmpfile")
    echo "  Found ${n_seqs} records (capped at ${MAX_SEQS})"

    if [[ "$n_seqs" -eq 0 ]]; then
        echo "  WARNING: No sequences found for ${marker}" >&2
        rm -f "$tmpfile"
        return
    fi

    # Download each sequence individually with standardized naming
    while IFS=$'\t' read -r accession organism; do
        # Build filename: Genus_species_accession_marker.fasta
        local genus species filename
        genus=$(echo "$organism" | awk '{print $1}')
        species=$(echo "$organism" | awk '{print $2}')
        # Sanitize: remove special characters
        genus=$(echo "$genus" | tr -cd '[:alnum:]_')
        species=$(echo "$species" | tr -cd '[:alnum:]_')
        species=${species:-sp}
        accession_clean=$(echo "$accession" | tr -cd '[:alnum:]._')
        filename="${marker_dir}/${genus}_${species}_${accession_clean}_${marker}.fasta"

        if [[ -f "$filename" ]]; then
            echo "  Skipping (exists): $(basename "$filename")"
            continue
        fi

        efetch -db nuccore -id "$accession" -format fasta 2>/dev/null \
            > "$filename" || {
            echo "  WARNING: Failed to download ${accession}" >&2
            rm -f "$filename"
            continue
        }

        # Verify non-empty
        if [[ ! -s "$filename" ]]; then
            echo "  WARNING: Empty file for ${accession}, removing" >&2
            rm -f "$filename"
            continue
        fi

        # ── Size check: reject complete plastome sequences ────────────────────
        # "marker"[All Fields] also matches full plastome records; efetch -format
        # fasta returns the entire sequence (~100-170 kb).  If the file exceeds
        # the per-marker threshold, re-fetch the GenBank annotation and extract
        # just the gene feature.
        local file_len size_max
        file_len=$(awk '/^>/{next}{len+=length($0)}END{print len+0}' "$filename")
        size_max=$(get_marker_size_max "$marker")

        if [[ "$file_len" -gt "$size_max" ]]; then
            echo "  Oversized (${file_len} bp > ${size_max} bp max for ${marker}) — attempting GB extraction..."
            local gb_file="${filename%.fasta}.gb"
            efetch -db nuccore -id "$accession" -format gb 2>/dev/null > "$gb_file"
            if "$HOME/miniconda3/bin/python3" \
                    "${PHYLOSKILLS_ROOT}/scripts/data/extract_gene_from_gb.py" \
                    --gb "$gb_file" \
                    --gene "$marker" \
                    --out "$filename" \
                    --header "${genus}_${species}_${accession_clean}_${marker} source=GenBank_annotation"; then
                echo "  Extracted ${marker} from complete plastome ${accession}"
            else
                echo "  WARNING: Could not extract ${marker} from ${accession} — skipping" >&2
                rm -f "$filename"
            fi
            rm -f "$gb_file"
        else
            echo "  Downloaded: $(basename "$filename")"
        fi

    done < "$tmpfile"

    rm -f "$tmpfile"

    local final_count
    final_count=$(find "${marker_dir}" -name "*.fasta" | wc -l)
    echo "  Marker ${marker}: ${final_count} sequences saved to ${marker_dir}/"
    echo ""
}

# ── Main download loop ────────────────────────────────────────────────────────
for marker in "${MARKERS[@]}"; do
    download_marker "$marker"
done

# ── Summary ───────────────────────────────────────────────────────────────────
echo "=== Download complete ==="
echo "Output directory: ${OUTDIR}"
for marker in "${MARKERS[@]}"; do
    n=$(find "${OUTDIR}/${marker}" -name "*.fasta" 2>/dev/null | wc -l)
    echo "  ${marker}: ${n} sequences"
done
echo ""
echo "Next step: Check coverage matrix, then proceed to phylo-alignment or phylo-assemble."
