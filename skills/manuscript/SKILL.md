---
name: manuscript
description: Use after tree inference and visualization are complete. Drafts the Methods section by reading provenance JSON logs produced by each analysis script, filling journal-specific paragraph templates with actual parameter values. Also generates figure captions from tree files and run metadata. Use when the researcher needs to write or finalize the Methods section, generate figure captions, or format references for a target journal.
---

# Manuscript Writing

## Overview

The Methods section writes itself from the analysis logs. Every script in this pipeline writes a
provenance JSON file (`results/provenance/<script>_YYYY-MM-DD.json`) recording tool versions,
parameters, input checksums, and runtime. `scripts/manuscript/methods_gen.py` reads these logs,
selects the correct paragraph template from `references/methods-templates.md`, fills in the actual
values, and outputs a draft Methods section ready for editing.

## Step 1 — Verify provenance logs exist

```bash
ls results/provenance/
# Expected files for a complete run:
# download_genbank_YYYY-MM-DD.json
# assemble_plastome_YYYY-MM-DD.json   (or run_hybpiper_YYYY-MM-DD.json)
# align_markers_YYYY-MM-DD.json
# build_gene_trees_YYYY-MM-DD.json
# run_aster_YYYY-MM-DD.json
# (optionally) beast_run_YYYY-MM-DD.json
```

If any provenance file is missing, the corresponding Methods paragraph will use placeholder
values — flag these to the researcher for manual completion.

## Step 2 — Choose target journal

Read `skills/manuscript/references/journal-formats.md` for:
- Word limits for Methods section
- Citation style (numbered vs. author-year)
- Figure format requirements (dimensions, resolution, format)
- Supplementary data policy (are raw sequences / trees required as supplementary?)

Common target journals for plant phylogenetics:

| Journal | Impact | Citation style | Methods limit |
|---------|--------|---------------|--------------|
| *Taxon* | ~3.0 | Author-year | None specified |
| *Systematic Botany* | ~2.5 | Author-year | None specified |
| *American Journal of Botany* | ~3.5 | Author-year | None specified |
| *Molecular Phylogenetics and Evolution* | ~4.0 | Numbered | None specified |
| *Systematic Biology* | ~14 | Numbered | Concise |
| *PLOS ONE* | ~3.5 | Numbered | None |
| PhytoKeys | ~1.0 | Author-year | None specified | 
| PhytoTaxa | ~1.0 | Author-year | None specified | 

Ask researcher to confirm target journal before generating output.

## Step 3 — Generate Methods draft

```bash
python scripts/manuscript/methods_gen.py \
  --provenance results/provenance/ \
  --templates skills/manuscript/references/methods-templates.md \
  --journal "Systematic Botany" \
  --output reports/methods_draft_YYYY-MM-DD.md
```

The script will:
1. Scan all JSON files in `results/provenance/`
2. Match each JSON to the appropriate template paragraph
3. Substitute actual values (tool versions, parameters, taxon counts, etc.)
4. Assemble paragraphs in pipeline order
5. Flag any missing provenance with `[FILL IN: ...]` placeholders

## Step 4 — Generate figure captions

This step is optional but recommended. It ensures that figure captions are consistent with the Methods section and contain all necessary details about tree inference and support values. The AI agent should ask the human researcher to see if to skip this step, but if not, run: 
```bash
python scripts/manuscript/fig_caption_gen.py \
  --trees results/trees/ \
  --provenance results/provenance/ \
  --output reports/figure_captions_YYYY-MM-DD.md
```

Generates captions for:
- ML tree (IQ-TREE): lists bootstrap method, replicates, model, partitioning
- ASTER species tree: states tool, number of gene trees, support metric
- BEAST2 chronogram (if present): calibration nodes, clock model, chain length

## Step 5 — Review and finalize

The generated Methods section is a draft, not a final product:
- All `[FILL IN: ...]` placeholders must be resolved
- Researcher must verify all parameter values match their intentions
- Add biological context sentences (why this group, why these markers) — the script cannot generate these
- Check journal-specific formatting (subheadings, paragraph order, abbreviations)

## QC Checklist

- [ ] All tools cited with version numbers (never "IQ-TREE was used")
- [ ] Bootstrap replicates stated explicitly
- [ ] Model selection criterion stated (BIC vs. AIC vs. AICc)
- [ ] Alignment strategy per marker documented
- [ ] Convergence diagnostics reported (ASDSF, ESS) for Bayesian runs
- [ ] BEAST2 calibration priors cited with source references
- [ ] Data availability statement: GenBank accessions, TreeBASE study ID, or Dryad DOI
- [ ] Software versions match `executables/software-inventory.md`
- [ ] All figures at journal-required resolution and format

## Report

**Mandatory:** Every log file generated during this module must be listed with its exact path in the report so the researcher can monitor background processes and audit what ran.

Write to `reports/[planX/]manuscript_YYYY-MM-DD.md`:

```markdown
# Manuscript Report
Date: YYYY-MM-DD
Journal: [target journal]

## Methods Draft
[paste or link generated methods section]

## Figure Captions
[paste or link generated captions]

## Placeholders Remaining
[list any [FILL IN] items not yet resolved]

## Data Availability
- GenBank accessions: [list or file path]
- Trees deposited: [TreeBASE / Dryad / GitHub]
- Analysis scripts: [GitHub repo URL]
```

## Log Files Generated
[List every log file created during this module with its exact path]
[Examples:]
[  reports/methods_draft_YYYY-MM-DD.md      (generated Methods section)]
[  reports/figure_captions_YYYY-MM-DD.md   (generated figure captions)]
[  scripts/manuscript/methods_gen.log       (script stdout, if redirected)]

## Common Mistakes

| Mistake | Fix |
|---------|-----|
| Not stating tool version | Read from provenance JSON; never write "IQ-TREE was used" without version |
| Reporting UFBoot as "bootstrap" without qualification | Write "ultrafast bootstrap (UFBoot)" and cite Hoang et al. 2018 |
| Omitting model selection criterion | Always state "BIC criterion" or "AICc criterion" |
| Not reporting convergence for Bayesian runs | ASDSF, ESS, and trace plot inspection must appear in Methods |
| Using placeholder accessions | Sequences must be deposited before submission; placeholder accessions are a common rejection reason |
| Generic "sequences were aligned with MAFFT" | State strategy flag (--linsi, --einsi), version, and per-marker justification |
| Forgetting data availability statement | Required by most journals; TreeBASE or Dryad deposition often mandatory |
