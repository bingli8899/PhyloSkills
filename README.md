# PhyloSkills

AI-assisted phylogenetic research toolkit — domain-expert skill library, deterministic
pipeline scripts, and an MCP server that wires them together.

## What this is

PhyloSkills separates **knowledge** from **execution**:

- `skills/<domain>/` — SKILL.md files that teach an AI how to reason about each
  pipeline stage (marker selection, alignment strategy, model choice, etc.)
- `scripts/<domain>/` — bash/Python scripts that execute each step deterministically;
  callable from the CLI without any AI involvement
- `mcp/server.py` — thin FastMCP wrapper; exposes scripts as AI-callable tools

Every script writes a provenance JSON log to `results/provenance/`. The
`scripts/manuscript/methods_gen.py` script reads those logs and auto-generates a
Methods section draft with `[FILL IN: ...]` placeholders for taxon-specific text.

## Repository layout

```
skills/
  pipeline/          # Routing hub — start here
  data-acquisition/  # GenBank / SRA download strategy
  assembly/          # Plastome assembly, HybPiper, CDS extraction
  alignment/         # MAFFT strategy, locus-guide.md, trimAl
  model-selection/   # ModelTest-NG, partition schemes
  tree-inference/    # IQ-TREE, ASTER, MrBayes, BEAST2
  visualization/     # ggtree, FigTree
  manuscript/        # Methods templates, journal formats
  environment/       # Tool installation, version tracking
  debug/             # Error diagnosis

scripts/
  data/              # download_genbank.sh, download_sra.sh
  assembly/          # assemble_plastome_*.sh, extract_cds.sh, run_hybpiper.sh
  alignment/         # align_markers.sh
  inference/         # build_gene_trees.sh, run_aster.sh
  manuscript/        # methods_gen.py, fig_caption_gen.py
  utils/             # change_headers.sh, provenance.sh, provenance_py.py

mcp/
  server.py          # FastMCP server — 15 tools wrapping scripts/

project-template/    # Copy this to start a new research project
  CLAUDE.md          # AI instructions for the project
  .gitignore         # Ignores raw data; keeps provenance logs
  data/accessions/   # GenBank/SRA accession lists go here
  results/provenance/
  reports/
    pipeline-status.md
  executables/
    software-inventory.md
```

## Quick start

### Running the MCP server

```bash
# Install dependencies
pip install -e .

# Run server (stdio transport — for Claude Desktop / claude CLI)
python -m mcp.server

# Or set repo root explicitly if running from outside the repo:
PHYLOSKILLS_ROOT=/path/to/PhyloSkills python -m mcp.server
```

Add to `~/.claude/claude_desktop_config.json`:
```json
{
  "mcpServers": {
    "phyloskills": {
      "command": "python",
      "args": ["-m", "mcp.server"],
      "cwd": "/path/to/PhyloSkills",
      "env": {"PHYLOSKILLS_ROOT": "/path/to/PhyloSkills"}
    }
  }
}
```

### Running scripts directly (no MCP)

```bash
# Align a set of CDS markers
bash scripts/alignment/align_markers.sh \
    -i data/cds/ -o data/aligned/ -s auto -P results/provenance

# Build gene trees
bash scripts/inference/build_gene_trees.sh \
    -i data/aligned/ -o results/trees/ -t 8 -P results/provenance

# Run wASTRAL coalescent
bash scripts/inference/run_aster.sh \
    -i results/trees/ -o results/trees/ -P results/provenance

# Generate Methods draft
python scripts/manuscript/methods_gen.py \
    --provenance results/provenance/ \
    --templates skills/manuscript/references/methods-templates.md \
    --journal "Systematic Botany" \
    --output reports/methods_draft_$(date +%Y-%m-%d).md
```

### Starting a new project

```bash
cp -r project-template/ ~/projects/my-taxon-study
cd ~/projects/my-taxon-study
# Edit CLAUDE.md and reports/pipeline-status.md with project details
git init && git add . && git commit -m "Initial project scaffold"
```

## Provenance

Every script accepts `-P <dir>` (bash) or `--provenance <dir>` (Python) to write
a JSON log:

```json
{
  "script": "align_markers",
  "tool": "mafft",
  "version": "7.526",
  "date": "2026-04-13T14:23:01",
  "parameters": {"strategy": "linsi", "threads": 8},
  "input_files": {"input_dir": {"path": "data/cds/", "md5": "directory"}},
  "output_files": {"alignment_dir": "data/aligned/"},
  "runtime_seconds": 142,
  "exit_code": 0
}
```

`methods_gen.py` reads these logs and substitutes values into prose templates.

## Skills overview

| Skill | Purpose |
|-------|---------|
| `pipeline` | Routing hub — which step to run next |
| `data-acquisition` | GenBank batch download, SRA streaming, accession management |
| `assembly` | GetOrganelle, BWA-plastome, HybPiper v2, CDS extraction |
| `alignment` | MAFFT strategy per marker, locus guide, trimAl, AMAS |
| `model-selection` | ModelTest-NG, partitioned models, AIC/BIC |
| `tree-inference` | IQ-TREE UFBoot/SH-aLRT, wASTRAL/ASTRAL-Pro3, MrBayes, BEAST2 |
| `visualization` | ggtree R package, FigTree, figure export |
| `manuscript` | Methods templates, figure captions, journal formats |
| `environment` | Tool install, version audit, conda environments |
| `debug` | Error diagnosis and recovery |
