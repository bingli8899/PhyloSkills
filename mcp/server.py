"""
PhyloSkills MCP Server
======================
Thin wrapper that exposes pipeline scripts as MCP tools.
All domain logic lives in scripts/. This server only validates inputs and calls scripts.

Usage:
    # Install dependencies
    pip install fastmcp

    # Run locally (stdio transport — for Claude Desktop / Claude Code)
    python mcp/server.py

    # Run as HTTP server (for remote / server deployment)
    fastmcp run mcp/server.py --transport streamable-http --port 8765

Configuration:
    Set PHYLOSKILLS_ROOT environment variable to override the repo root path.
    Defaults to the directory containing this file's parent (repo root).

Tool naming convention:  <domain>_<action>
    data_download_genbank
    data_download_sra
    assembly_run_getorganelle
    assembly_run_bwa
    assembly_run_hybpiper
    assembly_annotate_plastome
    assembly_extract_cds
    alignment_run
    alignment_analyze
    inference_build_gene_trees
    inference_run_aster
    manuscript_generate_methods
    manuscript_generate_captions
    utils_change_headers
    utils_remove_n_sequences
    utils_extract_taxa
"""

import os
import sys
import subprocess
import json
from pathlib import Path

try:
    from fastmcp import FastMCP
except ImportError:
    print("ERROR: fastmcp not installed. Run: pip install fastmcp", file=sys.stderr)
    sys.exit(1)

# ── Root path resolution ──────────────────────────────────────────────────────
REPO_ROOT = Path(os.environ.get("PHYLOSKILLS_ROOT", Path(__file__).parent.parent)).resolve()
SCRIPTS = REPO_ROOT / "scripts"

mcp = FastMCP(
    name="PhyloSkills",
    instructions=(
        "Phylogenetic analysis pipeline. Tools call deterministic scripts in scripts/. "
        "All tools accept a provenance_dir argument (default: results/provenance/) to log "
        "tool versions and parameters for auto-generating the Methods section. "
        "Read skills/<module>/SKILL.md for decision logic before calling tools. "
        "For plant markers, read skills/alignment/references/locus-guide.md before alignment."
    ),
)


# ── Shared helper ─────────────────────────────────────────────────────────────
def run_script(cmd: list[str], cwd: str | None = None) -> dict:
    """
    Run a subprocess command and return a structured result.
    Never raises — always returns stdout/stderr/exit_code.
    """
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=cwd or str(REPO_ROOT),
        )
        return {
            "exit_code": result.returncode,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "success": result.returncode == 0,
            "command": " ".join(cmd),
        }
    except FileNotFoundError as e:
        return {
            "exit_code": 1,
            "stdout": "",
            "stderr": f"Script not found: {e}",
            "success": False,
            "command": " ".join(cmd),
        }


def script_path(*parts: str) -> str:
    return str(SCRIPTS.joinpath(*parts))


# ═══════════════════════════════════════════════════════════════════════════════
# DATA ACQUISITION TOOLS
# ═══════════════════════════════════════════════════════════════════════════════

@mcp.tool()
def data_download_genbank(
    group: str,
    markers: str,
    output_dir: str,
    max_seqs: int = 500,
    survey_only: bool = False,
    provenance_dir: str = "results/provenance",
) -> dict:
    """
    Survey GenBank and download sequences for a taxonomic group and marker list.

    Calls scripts/data/download_genbank.sh. Enforces Genus_species_accession_marker.fasta
    naming. Use survey_only=True to check coverage before downloading.

    Args:
        group: Taxonomic group name (e.g. "Zingiberaceae", "Zingiber")
        markers: Comma-separated marker list (e.g. "matK,rbcL,ITS")
        output_dir: Output directory for downloaded sequences
        max_seqs: Maximum sequences per marker (default 500)
        survey_only: If True, only print coverage statistics; do not download
        provenance_dir: Directory for provenance JSON output
    """
    cmd = [
        "bash", script_path("data", "download_genbank.sh"),
        "-g", group,
        "-m", markers,
        "-o", output_dir,
        "-n", str(max_seqs),
    ]
    if survey_only:
        cmd.append("-s")
    return run_script(cmd)


@mcp.tool()
def data_download_sra(
    sra_list: str,
    output_dir: str,
    mode: str = "auto",
    threads: int = 4,
    assembly_script: str = "",
    provenance_dir: str = "results/provenance",
) -> dict:
    """
    Storage-aware SRA download with bulk or streaming mode.

    Calls scripts/data/download_sra.sh. Auto mode checks available storage vs.
    estimated dataset size and selects bulk (download all) or streaming
    (download→assemble→delete raw→next) accordingly.

    Args:
        sra_list: Path to text file with one SRA accession per line
        output_dir: Output directory for raw reads
        mode: "auto" | "bulk" | "streaming"
        threads: CPU threads for fasterq-dump
        assembly_script: Optional path to assembly script for streaming mode
        provenance_dir: Directory for provenance JSON output
    """
    cmd = [
        "bash", script_path("data", "download_sra.sh"),
        "-l", sra_list,
        "-o", output_dir,
        "-m", mode,
        "-c", str(threads),
    ]
    if assembly_script:
        cmd += ["-a", assembly_script]
    return run_script(cmd)


# ═══════════════════════════════════════════════════════════════════════════════
# ASSEMBLY TOOLS
# ═══════════════════════════════════════════════════════════════════════════════

@mcp.tool()
def assembly_run_getorganelle(
    read1: str,
    read2: str,
    output_dir: str,
    sample_name: str,
    plant_type: str = "embplant_pt",
    threads: int = 0,
    reference: str = "",
    provenance_dir: str = "results/provenance",
) -> dict:
    """
    Assemble plastome de novo with GetOrganelle (v1.7.7.1).

    Use for land plant plastome assembly from genome skimming or WGS reads.
    plant_type options: embplant_pt (plastome), embplant_mt (mitochondria),
    animal_mt, fungus_mt.

    Args:
        read1: Forward reads FASTQ path
        read2: Reverse reads FASTQ path
        output_dir: Base output directory (sample subdir created inside)
        sample_name: Sample name for subdirectory and output prefix
        plant_type: GetOrganelle -F flag (default: embplant_pt)
        threads: CPU threads (0 = auto-detect)
        reference: Optional reference FASTA for seed reads
        provenance_dir: Directory for provenance JSON output
    """
    cmd = [
        "bash", script_path("assembly", "assemble_plastome_getorganelle.sh"),
        "-1", read1, "-2", read2,
        "-o", output_dir,
        "-s", sample_name,
        "-F", plant_type,
    ]
    if threads > 0:
        cmd += ["-t", str(threads)]
    if reference:
        cmd += ["-R", reference]
    return run_script(cmd)


@mcp.tool()
def assembly_run_bwa(
    read1: str,
    read2: str,
    reference: str,
    output_dir: str,
    sample_name: str,
    consensus_type: str = "strict",
    min_depth: int = 3,
    threads: int = 0,
    provenance_dir: str = "results/provenance",
) -> dict:
    """
    Reference-guided plastome assembly: BWA + SAMtools + BCFtools consensus.

    Use when GetOrganelle fails or a closely related reference is available.
    consensus_type "strict" calls N at positions below min_depth (recommended);
    "majority" calls the most frequent allele regardless of depth.

    Args:
        read1: Forward reads FASTQ path
        read2: Reverse reads FASTQ path
        reference: Reference plastome FASTA
        output_dir: Output directory
        sample_name: Sample name prefix for output files
        consensus_type: "strict" (default) or "majority"
        min_depth: Minimum depth to call a base (default 3)
        threads: CPU threads (0 = auto-detect)
        provenance_dir: Directory for provenance JSON output
    """
    cmd = [
        "bash", script_path("assembly", "assemble_plastome_bwa.sh"),
        "-1", read1, "-2", read2,
        "-r", reference,
        "-o", output_dir,
        "-s", sample_name,
        "-c", consensus_type,
        "-d", str(min_depth),
    ]
    if threads > 0:
        cmd += ["-t", str(threads)]
    return run_script(cmd)


@mcp.tool()
def assembly_run_hybpiper(
    target_refs: str,
    sample_list: str,
    reads_dir: str,
    output_dir: str,
    threads: int = 0,
    seq_type: str = "dna",
    min_length: int = 10,
    run_paralogs: bool = False,
    provenance_dir: str = "results/provenance",
) -> dict:
    """
    Target enrichment assembly with HybPiper v2 (assemble + stats + retrieve).

    Use for Angiosperms353, custom HybSeq bait kits, or any target enrichment data.
    HybPiper v2 syntax — will error if v1.x is installed.
    Set run_paralogs=True to also run hybpiper paralog_retriever (use ASTRAL-Pro3
    for coalescent analysis if paralogs detected).

    Args:
        target_refs: Target reference FASTA (Angiosperms353_targetSequences.fasta or custom)
        sample_list: Text file with one sample name per line
        reads_dir: Directory containing paired reads (<sample>_R1.fastq.gz, etc.)
        output_dir: Output directory
        threads: CPU threads per sample (0 = auto-detect)
        seq_type: Sequence type for retrieval: "dna" | "aa" | "intron" | "supercontig"
        min_length: Minimum % target length for retrieval (default 10)
        run_paralogs: Also run paralog_retriever (use if multi-copy loci expected)
        provenance_dir: Directory for provenance JSON output
    """
    cmd = [
        "bash", script_path("assembly", "run_hybpiper.sh"),
        "-r", target_refs,
        "-s", sample_list,
        "-d", reads_dir,
        "-o", output_dir,
        "-m", str(min_length),
        "-T", seq_type,
    ]
    if threads > 0:
        cmd += ["-t", str(threads)]
    if run_paralogs:
        cmd.append("-e")
    return run_script(cmd)


@mcp.tool()
def assembly_extract_cds(
    input_path: str,
    output_dir: str,
    genes: str = "",
    min_length: int = 200,
    translate: bool = False,
    concatenate: bool = False,
    provenance_dir: str = "results/provenance",
) -> dict:
    """
    Extract CDS sequences from GenBank-annotated plastome files.

    Handles multi-exon genes by concatenating exons in order. By default extracts
    all 50 standard plastid protein-coding genes. Pass genes="" for the full default
    set or provide a comma-separated custom gene list.

    Args:
        input_path: GenBank .gb file or directory of .gb files
        output_dir: Output directory for per-gene FASTA files
        genes: Comma-separated gene names (empty = full default set)
        min_length: Minimum CDS length in bp (default 200)
        translate: Also write protein FASTA files
        concatenate: Write concatenated supermatrix + partition file
        provenance_dir: Directory for provenance JSON output
    """
    cmd = [
        "python3", script_path("assembly", "extract_cds.py"),
        "--input", input_path,
        "--output", output_dir,
        "--min_length", str(min_length),
    ]
    if genes:
        cmd += ["--genes", genes]
    if translate:
        cmd.append("--translate")
    if concatenate:
        cmd.append("--concatenate")
    return run_script(cmd)


# ═══════════════════════════════════════════════════════════════════════════════
# ALIGNMENT TOOLS
# ═══════════════════════════════════════════════════════════════════════════════

@mcp.tool()
def alignment_run(
    input_dir: str,
    output_dir: str,
    strategy: str = "auto",
    trim: bool = False,
    trim_mode: str = "automated1",
    concatenate: bool = False,
    threads: int = 0,
    provenance_dir: str = "results/provenance",
) -> dict:
    """
    Per-marker alignment with MAFFT (v7.526) + optional trimAl + optional concatenation.

    IMPORTANT: Read skills/alignment/references/locus-guide.md before choosing strategy.
    strategy options:
      auto    — safe default; MAFFT chooses internally
      linsi   — best for ≤200 seqs, moderate divergence (rbcL, matK, rpoB, Angiosperms353 same-family)
      einsi   — multiple conserved domains, large indels (trnL-F, ITS, psbA-trnH, A353 cross-order)
      ginsi   — high divergence, gappy datasets
      fftnsi  — fast; 200–10,000 sequences
      fftns2  — fastest; intraspecific/population data

    Args:
        input_dir: Directory with per-marker FASTA files
        output_dir: Output directory for aligned files
        strategy: MAFFT strategy (see above)
        trim: Run trimAl after alignment
        trim_mode: "automated1" (default) or "gappyout"
        concatenate: Concatenate all markers with AMAS into supermatrix
        threads: CPU threads (0 = auto-detect)
        provenance_dir: Directory for provenance JSON (used by methods_gen.py)
    """
    cmd = [
        "bash", script_path("alignment", "align_markers.sh"),
        "-i", input_dir,
        "-o", output_dir,
        "-s", strategy,
        "-P", provenance_dir,
    ]
    if trim:
        cmd += ["-T", "-m", trim_mode]
    if concatenate:
        cmd.append("-c")
    if threads > 0:
        cmd += ["-t", str(threads)]
    return run_script(cmd)


@mcp.tool()
def alignment_analyze(
    input_path: str,
    output_tsv: str = "alignment_diagnostics.tsv",
    max_gap_pct: float = 50.0,
    min_pis_pct: float = 5.0,
    min_length: int = 100,
) -> dict:
    """
    Compute alignment quality statistics and flag outlier sequences.

    Reports gap%, parsimony-informative sites (PIS%), and per-sequence outliers.
    Flags to route to the debug skill: HIGH_GAP, LOW_PIS, OUTLIER_SEQS.

    Args:
        input_path: Aligned FASTA file or directory of aligned FASTA files
        output_tsv: Output TSV file path for full report
        max_gap_pct: Maximum gap% per sequence before flagging (default 50.0)
        min_pis_pct: Minimum PIS% before flagging low-variation marker (default 5.0)
        min_length: Minimum ungapped sequence length in bp (default 100)
    """
    cmd = [
        "python3", script_path("alignment", "analyze_alignment.py"),
        "--input", input_path,
        "--output", output_tsv,
        "--max_gap_pct", str(max_gap_pct),
        "--min_pis_pct", str(min_pis_pct),
        "--min_length", str(min_length),
    ]
    return run_script(cmd)


# ═══════════════════════════════════════════════════════════════════════════════
# TREE INFERENCE TOOLS
# ═══════════════════════════════════════════════════════════════════════════════

@mcp.tool()
def inference_build_gene_trees(
    input_dir: str,
    output_dir: str,
    model: str = "TEST",
    ufboot: int = 1000,
    alrt: int = 1000,
    threads: int = 2,
    parallel_jobs: int = 1,
    provenance_dir: str = "results/provenance",
) -> dict:
    """
    Per-gene ML tree inference with IQ-TREE 2 (v2.4.0) for ASTER input.

    Runs IQ-TREE on each aligned marker FASTA in input_dir, then collects all
    .treefile outputs into all_gene_trees.txt. IQ-TREE output always includes
    branch lengths → use inference_run_aster with tool="wastral" (default).

    UFBoot threshold: ≥ 95 (NOT ≥ 70 — UFBoot is inflated).
    SH-aLRT threshold: ≥ 80. A node is well-supported when BOTH thresholds met.

    Args:
        input_dir: Directory with aligned marker FASTA files
        output_dir: Output directory for gene trees (subdir per gene)
        model: "TEST" (ModelFinder per gene, slower but correct) or specific model
        ufboot: Ultrafast bootstrap replicates (default 1000)
        alrt: SH-aLRT replicates (default 1000)
        threads: CPU threads per IQ-TREE job (default 2)
        parallel_jobs: Number of IQ-TREE jobs to run in parallel (default 1)
        provenance_dir: Directory for provenance JSON (used by methods_gen.py)
    """
    cmd = [
        "bash", script_path("inference", "build_gene_trees.sh"),
        "-i", input_dir,
        "-o", output_dir,
        "-m", model,
        "-B", str(ufboot),
        "-A", str(alrt),
        "-t", str(threads),
        "-j", str(parallel_jobs),
        "-P", provenance_dir,
    ]
    return run_script(cmd)


@mcp.tool()
def inference_run_aster(
    gene_trees: str,
    output_dir: str,
    tool: str = "auto",
    threads: int = 0,
    weight_mode: int = 2,
    concat_alignment: str = "",
    aster_bin_dir: str = "",
    provenance_dir: str = "results/provenance",
) -> dict:
    """
    Coalescent species tree inference with ASTER (v1.23).

    Auto-detection logic:
      - If gene trees have branch lengths (IQ-TREE output always does) → wASTRAL
      - If paralogs/multi-copy genes detected → override with tool="astral-pro3"
      - If no branch lengths → standard ASTRAL

    Local posterior support threshold: ≥ 0.95 for well-supported nodes.

    Args:
        gene_trees: Path to all_gene_trees.txt (one Newick per line)
        output_dir: Output directory for species tree
        tool: "auto" | "wastral" | "astral-pro3" | "astral"
        threads: CPU threads (0 = auto-detect)
        weight_mode: wASTRAL weighting (2 = branch lengths + local posterior, default)
        concat_alignment: Alignment file for concordance factor computation (optional)
        aster_bin_dir: Directory containing ASTER binaries (if not in PATH)
        provenance_dir: Directory for provenance JSON (used by methods_gen.py)
    """
    cmd = [
        "bash", script_path("inference", "run_aster.sh"),
        "-i", gene_trees,
        "-o", output_dir,
        "-T", tool,
        "-u", str(weight_mode),
        "-P", provenance_dir,
    ]
    if threads > 0:
        cmd += ["-t", str(threads)]
    if concat_alignment:
        cmd += ["-c", concat_alignment]
    if aster_bin_dir:
        cmd += ["-a", aster_bin_dir]
    return run_script(cmd)


# ═══════════════════════════════════════════════════════════════════════════════
# MANUSCRIPT TOOLS
# ═══════════════════════════════════════════════════════════════════════════════

@mcp.tool()
def manuscript_generate_methods(
    provenance_dir: str = "results/provenance",
    output: str = "",
    journal: str = "",
    sections: str = "",
) -> dict:
    """
    Auto-generate Methods section draft from provenance JSON logs.

    Reads all JSON files in provenance_dir, matches each to the appropriate
    paragraph template in skills/manuscript/references/methods-templates.md,
    substitutes actual parameter values, and writes a Markdown draft.
    All unfilled placeholders are marked [FILL IN: KEY] for manual completion.

    Args:
        provenance_dir: Directory containing pipeline provenance JSON files
        output: Output Markdown file path (default: reports/methods_draft_<date>.md)
        journal: Target journal name (for citation style note)
        sections: Comma-separated sections to include (empty = all)
    """
    from datetime import date as _date
    if not output:
        output = f"reports/methods_draft_{_date.today()}.md"

    templates = str(REPO_ROOT / "skills" / "manuscript" / "references" / "methods-templates.md")
    cmd = [
        "python3", script_path("manuscript", "methods_gen.py"),
        "--provenance", provenance_dir,
        "--templates", templates,
        "--output", output,
    ]
    if journal:
        cmd += ["--journal", journal]
    if sections:
        cmd += ["--section", sections]
    return run_script(cmd)


@mcp.tool()
def manuscript_generate_captions(
    trees_dir: str = "results/trees",
    provenance_dir: str = "results/provenance",
    output: str = "",
    journal: str = "",
) -> dict:
    """
    Generate figure captions for phylogenetic tree figures.

    Detects tree types (ML, ASTER, MrBayes, BEAST2) from file extensions and
    provenance logs, then produces ready-to-edit captions with [FILL IN:] for
    taxon-specific sentences (outgroup, biological context).

    Args:
        trees_dir: Directory containing tree files (.treefile, .nwk, .mcc.tree)
        provenance_dir: Directory containing pipeline provenance JSON files
        output: Output Markdown file path (default: reports/figure_captions_<date>.md)
        journal: Target journal (affects caption formatting notes)
    """
    from datetime import date as _date
    if not output:
        output = f"reports/figure_captions_{_date.today()}.md"

    cmd = [
        "python3", script_path("manuscript", "fig_caption_gen.py"),
        "--provenance", provenance_dir,
        "--trees", trees_dir,
        "--output", output,
    ]
    if journal:
        cmd += ["--journal", journal]
    return run_script(cmd)


# ═══════════════════════════════════════════════════════════════════════════════
# UTILITY TOOLS
# ═══════════════════════════════════════════════════════════════════════════════

@mcp.tool()
def utils_change_headers(
    input_fasta: str,
    mapping_file: str,
    output_fasta: str,
    strict: bool = False,
) -> dict:
    """
    Rename FASTA sequence headers from a two-column TSV mapping file.

    Column 1 = old name (matched against sequence ID before first space).
    Column 2 = new name. Use to standardize to Genus_species_accession_marker format.

    Args:
        input_fasta: Input FASTA file
        mapping_file: TSV with old_name (col 0) and new_name (col 1)
        output_fasta: Output FASTA file with renamed headers
        strict: Exit with error if any sequence has no mapping (default: warn only)
    """
    cmd = [
        "python3", script_path("utils", "change_header_name.py"),
        "--input", input_fasta,
        "--map", mapping_file,
        "--output", output_fasta,
    ]
    if strict:
        cmd.append("--strict")
    return run_script(cmd)


@mcp.tool()
def utils_remove_n_sequences(
    input_fasta: str,
    output_fasta: str,
    max_n_pct: float = 20.0,
    max_gap_pct: float = 50.0,
    min_length: int = 100,
    report: str = "",
) -> dict:
    """
    Remove sequences with excessive N/ambiguous bases or gaps.

    Useful after plastome assembly or consensus calling to remove low-quality
    sequences before alignment.

    Args:
        input_fasta: Input FASTA file
        output_fasta: Output FASTA file (filtered)
        max_n_pct: Maximum % N/ambiguous bases allowed (default 20.0)
        max_gap_pct: Maximum % gap characters allowed (default 50.0)
        min_length: Minimum ungapped sequence length in bp (default 100)
        report: Optional TSV file for per-sequence statistics
    """
    cmd = [
        "python3", script_path("utils", "remove_N_fasta.py"),
        "--input", input_fasta,
        "--output", output_fasta,
        "--max_n_pct", str(max_n_pct),
        "--max_gap_pct", str(max_gap_pct),
        "--min_length", str(min_length),
    ]
    if report:
        cmd += ["--report", report]
    return run_script(cmd)


@mcp.tool()
def utils_extract_taxa(
    input_fasta: str,
    taxa_list: str,
    output_fasta: str,
    match_mode: str = "exact",
    invert: bool = False,
) -> dict:
    """
    Extract (or exclude) specific taxa from a FASTA file.

    Args:
        input_fasta: Input FASTA file
        taxa_list: Text file with one taxon name per line
        output_fasta: Output FASTA file
        match_mode: "exact" | "prefix" | "contains" (default: exact)
        invert: If True, EXCLUDE listed taxa (keep everything else)
    """
    cmd = [
        "bash", script_path("utils", "extract_taxa_from_fasta.sh"),
        "-i", input_fasta,
        "-l", taxa_list,
        "-o", output_fasta,
        "-m", match_mode,
    ]
    if invert:
        cmd.append("-v")
    return run_script(cmd)


# ── Entry point ───────────────────────────────────────────────────────────────
if __name__ == "__main__":
    mcp.run()
