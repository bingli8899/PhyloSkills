---
name: phylo-research-design
description: Use when starting a new phylogenetic study, when the biological question is undefined or needs refinement, or when a researcher needs guidance on taxon sampling strategy, marker selection, or outgroup choice. Use before phylo-data-acquisition. Requires literature review before any design recommendation is made.
---

# Designing a Phylogenetic Research Plan

## Overview

Turn a vague biological question into a concrete, literature-grounded research design. Always conduct a literature review first — never propose taxa, markers, or outgroups before understanding the current state of the field.

## Process

### Step 1 — Assess researcher level

Ask one question to gauge experience: have they run a phylogenetic analysis before?  
Calibrate all subsequent explanations accordingly:
- **Novice:** explain concepts briefly before asking each decision question
- **Intermediate:** offer options with trade-offs, minimal background
- **Expert:** present options directly, skip definitions

### Step 2 — Literature review (always first)

Before asking the researcher anything about their design, search the literature:

**Search sources** (use all that are accessible):
- PubMed / NCBI — peer-reviewed molecular systematics papers
- bioRxiv — recent preprints for the group
- Google Scholar — broader coverage, reviews, book chapters
- Web of Science — citation-rich searches for landmark papers

**Search terms to adapt per group:**
```
"<group name>" AND (phylogeny OR phylogenetic OR phylogenom* OR plastome
  OR "genome skimming" OR "target enrichment" OR transcriptome)
"<group name>" AND (rbcL OR matK OR trnL OR ITS OR psbA OR nrITS)
"<group name>" AND (systematics OR taxonomy OR "species delimitation"
  OR "divergence time" OR biogeography)
```

**What to extract from the literature:**
- Current accepted taxonomy and major generic/tribal boundaries
- Which markers have been used and how well they resolved relationships
- Known problem areas (non-monophyletic genera, cryptic species, unsampled lineages)
- Published datasets that may be reusable (TreeBASE, Dryad, GenBank submissions)
- Calibration points used in previous divergence time studies (if relevant)

Present a concise literature review summary to the researcher before proceeding. Flag conflicts or gaps in current knowledge.

### Step 3 — Define the biological question (question-first)

Ask the researcher to state their biological question. If unclear, offer prompts:
- Are these taxa monophyletic / distinct?
- What are the relationships among genera / families?
- When did this group diversify?
- Where did this group originate (biogeography)?
- How did a trait evolve across the tree?

The question drives every downstream decision — do not proceed without a clear statement.

### Step 4 — Research design

Work backward from the question to each design element. Ask one at a time:

**Taxonomic scope and level**
- What is the ingroup? (any level: populations → orders)
- How many species/accessions are needed to answer the question?
- Are there known gaps in taxonomic coverage that would weaken the conclusion?

**Marker / genomic data strategy**
- What markers have worked for this group (from literature review)?
- Is whole-plastome or whole-genome data available? (more markers = better)
- For divergence time: are clock-like markers available?
- Prefer multi-locus or genomic data over single markers wherever possible

**Outgroup selection**
- Identify 2–3 outgroup taxa from a well-established sister group
- Confirm outgroup placement is supported in the literature
- Avoid outgroups too distant (long-branch attraction) or too close (inadequate rooting)

### Step 5 — Confirm design with researcher

Present the complete design as a summary and get explicit approval before writing the report. If the researcher requests changes, revise and re-confirm.

## QC Gate

Do not proceed to `phylo-data-acquisition` until:
- [ ] Biological question is clearly and specifically stated
- [ ] Target taxonomic scope and level are defined
- [ ] Marker/data strategy is chosen and grounded in literature
- [ ] Outgroup is identified and justified
- [ ] Researcher has approved the design

## Report

Write to `reports/research-design_YYYY-MM-DD.md` upon researcher approval.

```markdown
# Research Design Report
Date: YYYY-MM-DD

## Biological Question
[Exact statement of the research question]

## Literature Review Summary
[Key papers, current taxonomic consensus, known problem areas,
 markers previously used and their performance, available datasets]

## Research Design
- **Ingroup:** [taxa, taxonomic level, approximate species count]
- **Outgroup:** [taxa and justification]
- **Markers / data type:** [genes or genomic strategy, with rationale]
- **Taxonomic scope:** [population / species / genus / family / order]

## Known Gaps and Risks
[Unsampled lineages, conflicting taxonomies, data availability concerns]

## Next Module
phylo-data-acquisition
```

Update the report if the design changes after initial approval.

## Common Mistakes

| Mistake | Fix |
|---------|-----|
| Proposing markers before reviewing the literature | Always review first — prior studies reveal what works for the group |
| Asking multiple design questions at once | One question at a time; the researcher may not know all answers immediately |
| Choosing outgroups by familiarity rather than phylogenetic position | Confirm outgroup placement in a recent backbone phylogeny |
| Treating the question as optional | Every design decision must be traceable back to the biological question |
| Skipping the literature review for "well-known" groups | Even well-studied groups have recent revisions — always search |
