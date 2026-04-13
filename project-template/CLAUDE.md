# Project AI Instructions

This is a phylogenetic research project using the PhyloSkills pipeline.
AI assistants working in this directory should follow these guidelines.

## Project Overview

[FILL IN: one sentence describing the taxon group and research question]

**Example:** "This project infers a dated phylogeny of [TAXON] using [MARKERS] to resolve
[PROBLEM], with [N] taxa sampled from [GEOGRAPHIC SCOPE]."

## PhyloSkills Location

The PhyloSkills skill library and pipeline scripts are at:
```
PHYLOSKILLS_ROOT=~/path/to/PhyloSkills
```
Scripts live at `$PHYLOSKILLS_ROOT/scripts/<domain>/`.
Skills (domain expertise) live at `$PHYLOSKILLS_ROOT/skills/<domain>/SKILL.md`.

When helping with this project, read the relevant SKILL.md before writing commands.
Start with `skills/pipeline/SKILL.md` for the routing overview.

## Directory Layout

```
data/
  accessions/        # GenBank/SRA accession lists (.txt)
  cds/               # Extracted CDS FASTA files, one per marker
  aligned/           # Post-alignment FASTA files
  angiosperms353/    # A353 target file and extracted loci

results/
  provenance/        # Pipeline provenance JSON logs (auto-generated)
  trees/             # Tree files (.treefile, .contree, *_wastral.nwk, *.mcc.tree)
  beast2/            # BEAST2 XML and output
  mrbayes/           # MrBayes nexus and output
  getorganelle/      # GetOrganelle assembly output

reports/
  methods_draft_*.md      # Auto-generated from provenance
  figure_captions_*.md    # Auto-generated from tree files
  pipeline-status.md      # Hand-maintained progress tracker
```

## Provenance Convention

Every pipeline script writes a JSON log to `results/provenance/` automatically
when called with `-P results/provenance`. Do not skip this flag.
The `methods_gen.py` script reads these logs to generate the Methods section.

## Active Task

[FILL IN: what is currently being done]

## Notes for AI

- Read `executables/software-inventory.md` before citing tool versions.
- Use `reports/pipeline-status.md` to understand current progress.
- All scripts are at `$PHYLOSKILLS_ROOT/scripts/<domain>/<script>.sh` or `.py`.
- Never hardcode absolute paths; use relative paths from project root.
- When generating code, follow the patterns in existing scripts.
