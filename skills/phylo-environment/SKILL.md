---
name: phylo-environment
description: Use when starting a phylogenetic project, before running any analysis module, or when a module reports a missing or incompatible tool. Checks, installs, and versions all third-party software required for the pipeline. Use when the researcher needs to set up tools, verify an existing environment, or when any downstream module cannot find a required executable.
---

# Managing the Phylogenetic Software Environment

## Overview

Verify, install, and document all pipeline software before analysis begins. No module runs until its required tools are confirmed present and version-recorded.

## When to Use

- Project start — before `phylo-research-design` or any other module
- Any module reports a missing or unrecognized executable
- Researcher switches machines or environments
- A tool needs upgrading mid-project

## Setup

### Executables folder

All tools install to `executables/<tool-name>/` at the project root. Always resolve tool paths from here first; fall back to system PATH with a warning.

### Two modes — ask the researcher which to use

**Mode 1 — Human-managed**
Researcher installs tools themselves. For each missing tool, provide:
- Official download URL and recommended install method (conda, brew, binary release)
- Expected path in `executables/`
- Verification command

Wait for researcher confirmation before proceeding.

**Mode 2 — AI-assisted**
For each missing tool:
1. Identify correct binary for the researcher's OS
2. Download from official source into `executables/<tool-name>/`
3. Verify with a test call (e.g., `--version` or `--help`)
4. Record in `executables/software-inventory.md`

## Required Tools by Module

| Module | Tools |
|--------|-------|
| `phylo-data-acquisition` | NCBI Entrez Direct (`edirect`), SRA Toolkit (`prefetch`, `fasterq-dump`) |
| `phylo-assemble` | GetOrganelle, NOVOPlasty, HybPiper, BWA, SAMtools, BLAST+, Trinity *(load only what the chosen assembly strategy needs)* |
| `phylo-alignment` | MAFFT *(primary)*; MUSCLE *(alternative)* |
| `phylo-model-selection` | IQ-TREE 2 *(includes ModelFinder)*; ModelTest-NG *(alternative)* |
| `phylo-tree-inference` | IQ-TREE 2, RAxML-NG *(alternative)*, MrBayes *(Bayesian)*, BEAST2 *(divergence time)* |
| `phylo-visualization` | R ≥ 4.0, `ape`, `phytools`, `ggtree`, `ggplot2` |

Only check tools needed for the planned modules — do not install the full list upfront unless requested.

## Version Recording

Run `--version` (or equivalent) for every confirmed tool. Record output verbatim — never guess version numbers.

Every downstream module report must include this table (copy from inventory):

```markdown
| Tool | Version | Source | Install date |
|------|---------|--------|--------------|
| MAFFT | 7.520 | https://mafft.cbrc.jp | 2026-04-13 |
| IQ-TREE | 2.3.6 | github.com/iqtree/iqtree2/releases | 2026-04-13 |
```

## Software Inventory File

Maintain `executables/software-inventory.md` as a running record throughout the project:

```markdown
# Software Inventory

| Tool | Version | Source URL | Install date | Platform |
|------|---------|-----------|-------------|---------|
| ... | ... | ... | ... | ... |
```

Update this file every time a tool is added or upgraded. It is the project's reproducibility record.

## Hard Gate

**Do not allow any module to proceed if its required tools are absent or unverified.**

If a tool cannot be installed in either mode, report the blocker clearly:
- What is missing
- Why it is needed
- What the researcher must do to resolve it

## Common Mistakes

| Mistake | Fix |
|---------|-----|
| Recording version from memory or docs | Always run `--version` at install time |
| Installing to system PATH only | Also symlink or copy to `executables/` so inventory stays accurate |
| Checking all tools upfront | Only verify tools for planned modules — avoids unnecessary installs |
| Skipping version record mid-project when upgrading | Update `software-inventory.md` immediately after any upgrade |
| Proceeding when a tool test call fails | Stop and debug — a silent failure here causes cryptic errors downstream |
