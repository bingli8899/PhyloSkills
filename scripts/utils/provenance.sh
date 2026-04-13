#!/usr/bin/env bash
# =============================================================================
# provenance.sh — Shared provenance JSON writer for pipeline scripts
# =============================================================================
# Source this file in any pipeline script, then call write_provenance() at
# completion. Every script that calls this produces a JSON log in results/provenance/
# that methods_gen.py reads to auto-populate the Methods section.
#
# Usage (in your script):
#   source scripts/utils/provenance.sh
#   # ... run your analysis ...
#   write_provenance \
#     --script "align_markers" \
#     --tool "mafft" \
#     --version "$(mafft --version 2>&1 | head -1)" \
#     --params '{"strategy":"linsi","threads":8,"trim":true}' \
#     --inputs '{"input_dir":"data/cds/"}' \
#     --outputs '{"alignment_dir":"data/aligned/","stats":"alignment_stats.tsv"}' \
#     --provenance_dir "results/provenance" \
#     --exit_code 0 \
#     --start_time "$START_TIME"
# =============================================================================

# Capture script start time at source time
_PROV_START_EPOCH=$(date +%s)

write_provenance() {
    local script="" tool="" version="" params="{}" inputs="{}" outputs="{}"
    local provenance_dir="results/provenance" exit_code=0 start_time=""

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --script)       script="$2"; shift 2 ;;
            --tool)         tool="$2"; shift 2 ;;
            --version)      version="$2"; shift 2 ;;
            --params)       params="$2"; shift 2 ;;
            --inputs)       inputs="$2"; shift 2 ;;
            --outputs)      outputs="$2"; shift 2 ;;
            --provenance_dir) provenance_dir="$2"; shift 2 ;;
            --exit_code)    exit_code="$2"; shift 2 ;;
            --start_time)   start_time="$2"; shift 2 ;;
            *) shift ;;
        esac
    done

    [[ -z "$script" ]] && { echo "write_provenance: --script is required" >&2; return 1; }

    local end_epoch
    end_epoch=$(date +%s)
    local start_epoch="${start_time:-$_PROV_START_EPOCH}"
    local runtime=$(( end_epoch - start_epoch ))
    local datestamp
    datestamp=$(date +%Y-%m-%dT%H:%M:%S)
    local date_short
    date_short=$(date +%Y-%m-%d)

    mkdir -p "$provenance_dir"
    local outfile="${provenance_dir}/${script}_${date_short}.json"

    # Build checksums for input files (if inputs is a JSON object with file paths)
    # Simple approach: extract path values and compute md5
    local input_checksums="{}"
    if command -v python3 &>/dev/null; then
        input_checksums=$(python3 - "$inputs" <<'PYEOF'
import sys, json, hashlib, os
try:
    inp = json.loads(sys.argv[1])
    checksums = {}
    for k, v in inp.items():
        if isinstance(v, str) and os.path.exists(v):
            h = hashlib.md5()
            if os.path.isfile(v):
                with open(v, 'rb') as f:
                    for chunk in iter(lambda: f.read(65536), b''):
                        h.update(chunk)
                checksums[k] = {"path": v, "md5": h.hexdigest()}
            else:
                checksums[k] = {"path": v, "md5": "directory"}
        else:
            checksums[k] = {"path": str(v)}
    print(json.dumps(checksums))
except Exception as e:
    print("{}", file=sys.stderr)
    print("{}")
PYEOF
        )
    fi

    # Write JSON
    python3 - <<PYEOF
import json
from datetime import datetime

data = {
    "script": "$script",
    "tool": "$tool",
    "version": "$version",
    "date": "$datestamp",
    "parameters": $params,
    "input_files": $input_checksums,
    "output_files": $outputs,
    "runtime_seconds": $runtime,
    "exit_code": $exit_code,
    "host": "$(hostname 2>/dev/null || echo unknown)",
    "working_dir": "$(pwd)",
}
with open("$outfile", "w") as f:
    json.dump(data, f, indent=2)
print(f"Provenance written: $outfile")
PYEOF
}

# Convenience: call at script exit with the collected exit code
# Usage: trap 'finalize_provenance $? script tool version params inputs outputs prov_dir' EXIT
finalize_provenance() {
    local exit_code="$1"
    shift
    write_provenance "$@" --exit_code "$exit_code"
}
