#!/usr/bin/env bash
# =============================================================================
# extract_taxa_from_fasta.sh — Extract specific taxa from a FASTA file
# =============================================================================
# Filters a FASTA file to keep only sequences matching a list of taxon names.
# Matching is done against the sequence ID (portion of header before first space).
# Supports exact match and prefix/partial match modes.
#
# Usage:
#   bash extract_taxa_from_fasta.sh \
#     -i <input.fasta> \
#     -l <taxa_list.txt> \
#     -o <output.fasta> \
#     [-m <match_mode>] \
#     [-v]
#
# Arguments:
#   -i  Input FASTA file
#   -l  Text file with one taxon name per line (lines starting with # ignored)
#         Names are matched against the sequence ID (before first space in header)
#   -o  Output FASTA file
#   -m  Match mode: exact | prefix | contains (default: exact)
#         exact    = full ID must match entry in list
#         prefix   = ID must start with entry in list
#         contains = entry in list is a substring of ID
#   -v  Invert selection — exclude listed taxa (keep everything else)
#
# Tool version requirements:
#   Python >= 3.8  (no external dependencies beyond standard library)
#   AWK (for --stats only)
#
# =============================================================================

set -euo pipefail

MATCH_MODE="exact"
INVERT=false

usage() {
    grep '^#' "$0" | grep -v '#!/' | sed 's/^# \{0,1\}//'
    exit 1
}

while getopts "i:l:o:m:vh" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        l) TAXA_LIST="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        m) MATCH_MODE="$OPTARG" ;;
        v) INVERT=true ;;
        h|*) usage ;;
    esac
done

[[ -z "${INPUT:-}" || -z "${TAXA_LIST:-}" || -z "${OUTPUT:-}" ]] && {
    echo "ERROR: -i, -l, and -o are required." >&2; usage; }
[[ ! -f "$INPUT" ]]     && { echo "ERROR: Input not found: ${INPUT}" >&2; exit 1; }
[[ ! -f "$TAXA_LIST" ]] && { echo "ERROR: Taxa list not found: ${TAXA_LIST}" >&2; exit 1; }

echo "# extract_taxa_from_fasta.sh"
echo "# Input:      ${INPUT}"
echo "# Taxa list:  ${TAXA_LIST}"
echo "# Match mode: ${MATCH_MODE}"
echo "# Invert:     ${INVERT}"

# Use Python for robust FASTA parsing (avoids multiline/wrapping issues)
python3 - "$INPUT" "$TAXA_LIST" "$OUTPUT" "$MATCH_MODE" "$INVERT" <<'PYEOF'
import sys

fasta_path, list_path, out_path, match_mode, invert_str = sys.argv[1:6]
invert = invert_str == "True"

# Load taxa list
taxa = set()
with open(list_path) as fh:
    for line in fh:
        line = line.strip()
        if line and not line.startswith("#"):
            taxa.add(line)

print(f"# Loaded {len(taxa)} taxon names from list")

def matches(seq_id: str) -> bool:
    if match_mode == "exact":
        return seq_id in taxa
    elif match_mode == "prefix":
        return any(seq_id.startswith(t) for t in taxa)
    elif match_mode == "contains":
        return any(t in seq_id for t in taxa)
    return seq_id in taxa

n_kept = 0
n_skipped = 0
current_keep = False

with open(fasta_path) as fin, open(out_path, "w") as fout:
    seq_buffer = []
    current_header = None

    def flush():
        nonlocal n_kept, n_skipped
        if current_header is None:
            return
        seq_id = current_header.split()[0]
        match = matches(seq_id)
        keep = match ^ invert  # XOR: invert flips the keep decision
        if keep:
            fout.write(f">{current_header}\n")
            for chunk in seq_buffer:
                fout.write(chunk + "\n")
            n_kept += 1
        else:
            n_skipped += 1

    for line in fin:
        line = line.rstrip()
        if line.startswith(">"):
            flush()
            current_header = line[1:]
            seq_buffer = []
        elif line:
            seq_buffer.append(line)
    flush()

print(f"# Kept:    {n_kept}")
print(f"# Skipped: {n_skipped}")
print(f"# Output:  {out_path}")
PYEOF
