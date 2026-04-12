# PhyloSkills Design Spec
**Date:** 2026-04-12  
**Status:** Approved  
**Author:** Collaborative brainstorm (human + AI)

---

## Overview

PhyloSkills is a hub-and-spoke skill suite that coaches an AI agent to assist a phylogeneticist through the full research pipeline — from forming a biological question to producing a publication-quality tree. Skills adapt to the researcher's experience level and are independently invokable or chainable as a complete pipeline.

**Target user:** AI agent assisting a phylogeneticist with limited bioinformatics experience  
**Organism scope:** Any (population to order-level; organism-agnostic)  
**Tool philosophy:** Recommend common open-source defaults; guide the decision, not just the command

---

## Architecture

```
phylogenetic-analysis (hub)
├── phylo-research-design
├── phylo-data-acquisition
├── phylo-alignment
├── phylo-model-selection
├── phylo-tree-inference
├── phylo-visualization
└── phylo-debug
```

The hub reads existing checkpoint reports to determine pipeline state and routes the researcher to the correct module. Each module ends with a finalized report written to the `reports/` folder.

---

## Cross-Cutting: Report System

Every module writes reports to `reports/` at completion (and updated continuously during the module).

### Report conventions
- **Filename format:** `reports/<module-name>_YYYY-MM-DD.<ext>`
- **Format:** Markdown (`.md`) by default; LaTeX (`.tex`) where specified
- **Readable by:** Both humans and AI agents
- **Purpose (dual):**
  1. Human audit trail — what the AI did, what decisions were made, why
  2. AI checkpoint — agent reads reports on session start to determine resume point

### Report minimum sections
- What was done (actions taken, commands run)
- Decisions made and justification
- Tool outputs / key results
- QC outcomes
- Current status (complete / blocked / needs review)
- Recommended next step

---

## Module Specifications

### 1. `phylo-research-design`

**Trigger:** Start of a new phylogenetic project or when research question is undefined

**Process:**
1. Assess researcher's experience level; calibrate explanation depth accordingly
2. **Literature review first** — search online for recent manuscripts and papers (PubMed, Google Scholar, bioRxiv, Web of Science) before proposing anything
3. Present a concise but comprehensive literature review summary to the researcher
4. Ask the biological question the researcher wants to answer (question-first approach)
5. Work backward from the question to: target taxa, appropriate markers/genes, outgroup strategy, taxonomic level (population / species / genus / family / order)
6. Confirm research design with researcher before proceeding

**QC gate:** Is the research question clearly defined and answerable with available public data?

**Report:** `reports/research-design_YYYY-MM-DD.md`
- Literature review summary (key papers, current state of knowledge)
- Biological question statement
- Target taxa and taxonomic level
- Selected markers / genomic regions
- Outgroup strategy and justification
- Recommended next module: `phylo-data-acquisition`

---

### 2. `phylo-data-acquisition`

**Trigger:** Research design is finalized; need to download sequence data

**Process:**
1. Read `phylo-research-design` report to understand organism group, markers, taxonomic scope
2. Adaptively recommend databases based on organism group + markers + taxonomic level (not a fixed list — reasoning driven by context)
3. Guide search query construction, taxon filters, sequence length filters, date filters
4. Help organize downloaded files with consistent naming conventions
5. Track every accession downloaded: ID, database, date, sequence metadata

**QC gate:** Sequence count check, format validation, file naming consistency, no duplicate accessions

**Reports (two outputs):**
- `reports/data-acquisition_YYYY-MM-DD.tex` — LaTeX format  
  - Table of every accession: accession ID, taxon name, database source, download date, sequence length, gene/marker
  - Suitable for supplementary materials in a manuscript
- `reports/data-acquisition_YYYY-MM-DD.md` — Narrative format  
  - Search strategy used for each database
  - Filters applied and why
  - Sequences included vs. excluded with reasoning
  - Data quality notes
  - Recommended next module: `phylo-alignment`

---

### 3. `phylo-alignment`

**Trigger:** Sequences downloaded and organized; ready for alignment

**Process:**
1. Read data-acquisition report to understand dataset size, sequence type (DNA / protein / rRNA), number of markers
2. Recommend alignment tool based on dataset characteristics:
   - Primary: MAFFT (most datasets), MUSCLE (alternative)
   - Large datasets: MAFFT --auto or --linsi
3. Guide parameter selection
4. Review alignment output quality (length, gap %, poorly aligned regions)
5. If multi-locus: guide concatenation or coalescent approach decision

**QC gate:** Alignment length reasonable, gap % within acceptable range, no obviously misaligned sequences  
**On QC failure:** Route to `phylo-debug`

**Report:** `reports/alignment_YYYY-MM-DD.md`
- Tool and version used
- Parameters applied
- Alignment statistics (number of sequences, alignment length, gap %)
- QC outcome
- Any sequences trimmed or removed with reasoning
- Recommended next module: `phylo-model-selection`

---

### 4. `phylo-model-selection`

**Trigger:** Alignment complete; need to select substitution model before tree inference

**Process:**
1. Read alignment report
2. Explain substitution models at the researcher's calibrated level
3. Run model testing:
   - Primary: IQ-TREE built-in ModelFinder (`-m TEST`)
   - Alternative: ModelTest-NG (standalone)
4. Guide interpretation of AIC / BIC / AICc results
5. For multi-locus data: guide partitioning scheme selection
6. Output selected model(s) ready for tree inference

**QC gate:** Model is appropriate for data type and sequence length

**Report:** `reports/model-selection_YYYY-MM-DD.md`
- Selected model(s) per partition
- AIC/BIC scores summary
- Justification for chosen model
- Partitioning scheme (if applicable)
- Model string ready to copy into tree inference command
- Recommended next module: `phylo-tree-inference`

---

### 5. `phylo-tree-inference`

**Trigger:** Alignment and model selection complete

**Process:**
1. Read alignment + model-selection reports
2. Guide ML vs Bayesian decision based on research question and available compute:
   - ML primary: IQ-TREE (most cases)
   - ML alternative: RAxML-NG
   - Bayesian: MrBayes (topology focus), BEAST2 (divergence times / phylodynamics)
3. Configure run: substitution model, bootstrap replicates or posterior sampling, partitions
4. Assess support values: bootstrap thresholds (≥70 UFBoot, ≥95 standard), posterior probabilities (≥0.95)
5. Convergence check for Bayesian runs (ESS, PSRF via Tracer)
6. Sanity check tree topology against known biology

**QC gate:** Support values adequate, topology biologically plausible, convergence confirmed (Bayesian)  
**On QC failure:** Route to `phylo-debug`

**Report:** `reports/tree-inference_YYYY-MM-DD.md`
- Tool, version, and exact command used
- Model(s) applied
- Bootstrap/posterior support summary
- Convergence diagnostics (Bayesian only)
- Notable topology features
- Output file paths (tree files)
- Recommended next module: `phylo-visualization`

---

### 6. `phylo-visualization`

**Trigger:** Tree inference complete; ready to visualize and annotate

**Process:**
1. Read tree-inference report; locate output tree files
2. **R-based exclusively** — guide R script writing using:
   - `ape` — tree reading, manipulation, basic plotting
   - `phytools` — advanced visualization, contmap, phylomorphospace
   - `ggtree` — grammar-of-graphics tree plotting (primary for publication figures)
   - `ggplot2` — annotations, themes, multi-panel figures
3. Guide annotation: support values, scale bar, tip labels, clade highlighting, outgroup rooting
4. Produce publication-quality figures (SVG / PDF output)

**QC gate:** Tree correctly rooted, support values displayed accurately, tip labels match original taxon names  
**On QC failure:** Route to `phylo-debug`

**Report:** `reports/visualization_YYYY-MM-DD.md`
- R session info and package versions
- R script used (inline or path reference)
- Description of each figure produced
- Output file paths
- Recommended next step: manuscript / supplementary materials

---

### 7. `phylo-debug`

**Trigger:** Any module QC gate failure OR researcher-reported problem at any stage

**Process:**
1. Identify which module triggered the debug call and read its report
2. Match reported symptom to known failure categories:
   - **Bad alignment quality** — high gap %, misaligned regions, wrong sequences
   - **Bad sequence data** — contamination, chimeras, wrong gene, truncated sequences
   - **Download failures** — network errors, accession not found, database access issues
   - **File name mismatches** — inconsistent naming between alignment, model, and tree files
   - **Non-convergence** — Bayesian ESS too low, split frequencies too high
   - **Unexpected topology** — known clades not recovered, outgroup placement wrong
   - **Software errors** — version conflicts, missing dependencies, malformed input files
3. For each diagnosis: explain cause, provide fix, re-enter pipeline at correct stage
4. Update the original module's report with debug outcome

**Proactive mode:** Each module's QC gate automatically calls `phylo-debug` on failure  
**Reactive mode:** Researcher can invoke directly when something goes wrong

**Report:** `reports/debug_YYYY-MM-DD.md`
- Originating module and symptom
- Diagnosis
- Fix applied
- Outcome
- Pipeline re-entry point

---

### 8. `phylogenetic-analysis` (Hub)

**Trigger:** Start of any phylogenetic project or when researcher needs guidance on where to begin

**Process:**
1. Check `reports/` folder for existing checkpoint reports
2. Assess pipeline state from reports — determine resume point
3. Assess researcher's experience level
4. Route to appropriate module
5. Track overall pipeline state across the session

**Decision logic:**
- No reports exist → start at `phylo-research-design`
- Research design report exists, no data report → `phylo-data-acquisition`
- Data report exists, no alignment → `phylo-alignment`
- Alignment done, no model → `phylo-model-selection`
- Model done, no tree → `phylo-tree-inference`
- Tree done, no visualization → `phylo-visualization`
- Any report shows blocked/failed status → `phylo-debug`

---

## Case Study: Plant Systematics

### Scenario
A graduate student asks: "Are these two genera in family Apiaceae (carrots/parsley family) actually distinct, or should one be synonymized with the other?"

### Dataset
- **Organism:** Flowering plants (Apiaceae), genus-level question
- **Markers:** Plastid (*matK*, *rbcL*) + nuclear (ITS2) — tri-locus dataset
- **Taxonomic scope:** ~30 ingroup species, 3 outgroup species

### Pipeline walkthrough
1. **Research design:** Literature review of Apiaceae systematics, prior treatments of target genera, established marker utility for family
2. **Data acquisition:** NCBI GenBank for all three markers; BOLD as secondary for ITS2; LaTeX accession table + MD search narrative
3. **Alignment:** MAFFT `--auto` per marker; concatenation for combined analysis
4. **Model selection:** IQ-TREE ModelFinder with partitioned scheme (one model per marker)
5. **Tree inference:** IQ-TREE ML + 1000 ultrafast bootstrap replicates
6. **Visualization:** `ggtree` + `ggplot2` — annotated tree with bootstrap support, clade boxes, scale bar
7. **Debug scenario:** File naming mismatch — sequences downloaded from NCBI used accession IDs as names, BOLD used species names; alignment file has mixed naming causing downstream errors → `phylo-debug` diagnoses, standardizes names, re-runs alignment

### Expected outputs
- `reports/research-design_2026-04-12.md`
- `reports/data-acquisition_2026-04-12.tex` + `reports/data-acquisition_2026-04-12.md`
- `reports/alignment_2026-04-12.md`
- `reports/model-selection_2026-04-12.md`
- `reports/tree-inference_2026-04-12.md`
- `reports/visualization_2026-04-12.md`
- `reports/debug_2026-04-12.md`

---

## Directory Structure

```
PhyloSkills/
├── README.md
├── docs/
│   └── superpowers/
│       └── specs/
│           └── 2026-04-12-phyloskills-design.md   ← this file
├── reports/                                         ← generated per project run
│   ├── research-design_YYYY-MM-DD.md
│   ├── data-acquisition_YYYY-MM-DD.tex
│   ├── data-acquisition_YYYY-MM-DD.md
│   ├── alignment_YYYY-MM-DD.md
│   ├── model-selection_YYYY-MM-DD.md
│   ├── tree-inference_YYYY-MM-DD.md
│   ├── visualization_YYYY-MM-DD.md
│   └── debug_YYYY-MM-DD.md
└── skills/
    ├── phylogenetic-analysis/
    │   └── SKILL.md
    ├── phylo-research-design/
    │   └── SKILL.md
    ├── phylo-data-acquisition/
    │   └── SKILL.md
    ├── phylo-alignment/
    │   └── SKILL.md
    ├── phylo-model-selection/
    │   └── SKILL.md
    ├── phylo-tree-inference/
    │   └── SKILL.md
    ├── phylo-visualization/
    │   └── SKILL.md
    └── phylo-debug/
        └── SKILL.md
```

---

## Open Questions / Future Scope

- Should `reports/` be per-project (subfolder by project name/date) or flat? Currently flat — revisit if multiple concurrent projects needed.
- BEAST2 / divergence time estimation is mentioned as an alternative in `phylo-tree-inference` but not fully specced — could become its own module later.
- Population-level analyses (haplotype networks, PopART, DnaSP) not covered in v1 — future `phylo-population` module.
