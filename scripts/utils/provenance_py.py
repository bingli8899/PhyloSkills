"""
provenance_py.py — Shared provenance writer for Python pipeline scripts

Import this module in any Python pipeline script and call write_provenance()
after analysis completes. Produces a JSON log readable by methods_gen.py.

Usage:
    from scripts.utils.provenance_py import write_provenance
    import time

    start = time.time()
    # ... run analysis ...
    write_provenance(
        script="align_markers",
        tool="mafft",
        version="7.526",
        provenance_dir="results/provenance",
        parameters={"strategy": "linsi", "threads": 8},
        input_files={"input_dir": "data/cds/"},
        output_files={"alignment_dir": "data/aligned/"},
        runtime_seconds=int(time.time() - start),
        exit_code=0,
    )
"""

import json
import hashlib
import os
import time
from datetime import datetime
from pathlib import Path


def checksum_path(path: str) -> str:
    """Return MD5 checksum of a file, or 'directory' for directories."""
    p = Path(path)
    if not p.exists():
        return "not_found"
    if p.is_dir():
        return "directory"
    h = hashlib.md5()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def write_provenance(
    script: str,
    tool: str,
    version: str,
    provenance_dir: str = "results/provenance",
    parameters: dict = None,
    input_files: dict = None,
    output_files: dict = None,
    runtime_seconds: int = 0,
    exit_code: int = 0,
) -> str:
    """
    Write a provenance JSON file to results/provenance/.

    Returns the path to the written file.
    """
    if parameters is None:
        parameters = {}
    if input_files is None:
        input_files = {}
    if output_files is None:
        output_files = {}

    # Compute checksums for input files that exist
    checksummed_inputs = {}
    for k, v in input_files.items():
        if isinstance(v, str) and (Path(v).exists()):
            checksummed_inputs[k] = {"path": v, "md5": checksum_path(v)}
        else:
            checksummed_inputs[k] = {"path": str(v)}

    date_short = datetime.now().strftime("%Y-%m-%d")
    datestamp = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")

    data = {
        "script": script,
        "tool": tool,
        "version": version,
        "date": datestamp,
        "parameters": parameters,
        "input_files": checksummed_inputs,
        "output_files": output_files,
        "runtime_seconds": runtime_seconds,
        "exit_code": exit_code,
        "host": os.uname().nodename if hasattr(os, "uname") else "unknown",
        "working_dir": str(Path.cwd()),
    }

    out_dir = Path(provenance_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{script}_{date_short}.json"

    with open(out_path, "w") as f:
        json.dump(data, f, indent=2)

    print(f"Provenance written: {out_path}")
    return str(out_path)
