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
├── phylo-assemble          ← NEW: assembly/extraction when raw reads available
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

**Trigger:** Research design is finalized; need to survey and download sequence data

**Process:**
1. Read `phylo-research-design` report to understand organism group, markers, taxonomic scope
2. **Data landscape survey first** — before downloading anything, search both GenBank and SRA to understand what data types are available for the group:
   - GenBank: assembled sequences, gene records, plastomes, whole genomes
   - SRA: raw read datasets (genome skimming, target enrichment/HybSeq, transcriptomes, WGS)
3. **Search strategy** — use structured keyword combinations:
   - Primary: `"<group name>" AND (phylogeny OR phylogenetic OR phylogenom* OR plastome OR "genome skimming" OR "target enrichment" OR transcriptome)`
   - Marker-specific: `(rbcL OR matK OR trnL OR ITS OR psbA OR nrITS)` combined with group name
   - Adapt terms to organism group (e.g., add plastid gene names specific to the lineage)
4. **Sampling matrix analysis** — after the survey, AI builds a taxa × markers/data-type availability matrix:
   - For each taxon found, record which markers/data types are available (e.g., GenBank *matK* only, SRA genome skimming, both, none)
   - Identify coverage gaps: taxa with rich data vs. taxa with only 1–2 markers
   - Compute coverage statistics: total taxa available, per-marker taxon counts, data type breakdown
5. **Generate sampling plans** — AI proposes 2–4 explicit plans representing different taxon/marker trade-offs, for example:
   - **Plan A (broad taxon sampling):** maximum taxa, minimum marker requirement (e.g., 80 taxa × 3 markers — include any taxon with ≥3 markers)
   - **Plan B (balanced):** moderate taxa, moderate markers (e.g., 50 taxa × 6 markers)
   - **Plan C (marker-rich):** fewer taxa, maximum markers/genomic data (e.g., 20 taxa × 10+ markers or whole plastome)
   - **Plan D (all data):** include everything available regardless of per-taxon gaps; use missing-data-tolerant analysis
   - Plans scale to what the actual data landscape supports — AI generates realistic plans based on the survey, not hypothetical ones
   - For each plan: list exact taxon count, marker/data-type count, estimated missing data %, and recommended downstream approach
6. **Human selects plan(s)** — researcher chooses one plan OR multiple plans (all selected plans are downloaded and run as parallel analyses through the pipeline)
   - If multiple plans selected: each plan gets its own subdirectory and report series (e.g., `reports/planA/`, `reports/planB/`)
7. **Data type decision** — applied per plan:
   - Assembled gene sequences → download FASTA from GenBank (route → `phylo-alignment`)
   - SRA genome skimming / WGS data → plastome or whole-genome approach (route → `phylo-assemble`)
   - SRA target enrichment (HybSeq/Angiosperms353) data → target capture assembly (route → `phylo-assemble`)
   - Mixed → combined strategy within the plan
   - **Principle: more markers = better; whole plastome/genome preferred over single markers when data permits**
8. Guide search query construction, taxon filters, sequence length filters, date ranges
9. Download confirmed data per plan; enforce consistent file naming conventions
10. Track every accession: ID, database, data type, download date, sequence metadata, plan assignment

**QC gate:** Sequence/SRA run count check, format validation, file naming consistency, no duplicate accessions  
**On QC failure:** Route to `phylo-debug`

**Reports (two outputs, written once per run; plan-specific data goes into plan subdirs):**
- `reports/data-acquisition_YYYY-MM-DD.tex` — LaTeX format
  - Table of every accession/SRA run across all plans: ID, taxon name, database source, data type, download date, sequence length or read count, gene/marker or library strategy, plan assignment
  - Suitable for supplementary materials in a manuscript
- `reports/data-acquisition_YYYY-MM-DD.md` — Narrative format
  - Data landscape survey summary (what was found, what data types exist for the group)
  - Sampling matrix: taxa × markers/data-type availability
  - All plans proposed, with taxon count, marker count, missing data %, and trade-off reasoning
  - Plan(s) selected by researcher and rationale
  - Search queries used for each database
  - Filters applied and why
  - Sequences/runs included vs. excluded with reasoning
  - Data quality notes
  - Recommended next module per plan: `phylo-assemble` (raw reads) or `phylo-alignment` (assembled sequences)

**Multi-plan directory convention (when researcher selects multiple plans):**
```
reports/
  data-acquisition_YYYY-MM-DD.tex   ← combined accession table across all plans
  data-acquisition_YYYY-MM-DD.md    ← survey narrative + plan comparison
  planA/
    assembly_YYYY-MM-DD.md
    alignment_YYYY-MM-DD.md
    model-selection_YYYY-MM-DD.md
    tree-inference_YYYY-MM-DD.md
    visualization_YYYY-MM-DD.md
  planB/
    alignment_YYYY-MM-DD.md
    ...
```

---

### 3. `phylo-assemble`

**Trigger:** `phylo-data-acquisition` determines raw reads are available in SRA (genome skimming, WGS, target enrichment, transcriptomes)

**Process:**
1. Read data-acquisition report to understand data type, SRA run IDs, and reference availability
2. Determine assembly strategy based on data type:
   - **Genome skimming / WGS → plastome assembly:**
     - Primary: GetOrganelle (chloroplast/mitochondrial genome assembly from WGS)
     - Alternative: NOVOPlasty (seed-and-extend assembler)
     - Output: assembled plastome FASTA → extract individual genes with custom scripts or mfannot/PGA
   - **Genome skimming / WGS → nuclear marker extraction:**
     - Map reads to reference gene sequences using BWA or Bowtie2 + SAMtools; call consensus
     - Or use target-bait references (Angiosperms353 probes) as mapping targets
   - **Target enrichment (HybSeq / Angiosperms353 / custom baits) → marker assembly:**
     - Primary: HybPiper (purpose-built for target enrichment data)
     - Output: assembled target gene FASTA files per sample
   - **Transcriptome → gene extraction:**
     - Assemble with Trinity; extract target genes via BLAST against reference
3. Assess reference availability:
   - Published whole-genome or plastome reference for the group → use as mapping reference
   - No close reference → use nearest published relative; flag to researcher
   - Multiple references available → recommend the most complete and recently published
4. Download SRA reads (prefetch + fasterq-dump)
5. Run appropriate assembly pipeline; review assembly statistics (coverage, completeness, length)
6. Collect all assembled sequences into organized FASTA files ready for alignment

**QC gate:** Assembly completeness check, minimum coverage threshold, gene recovery rate (for target enrichment), no contamination flags  
**On QC failure:** Route to `phylo-debug`

**Report:** `reports/assembly_YYYY-MM-DD.md`
- SRA runs processed
- Assembly method and reference used, with justification
- Assembly statistics per sample (coverage, length, completeness %)
- Genes/markers successfully recovered
- Samples that failed assembly and why
- Output file paths
- Recommended next module: `phylo-alignment`

---

### 4. `phylo-alignment`

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
   - Bayesian topology: MrBayes
   - Divergence time / phylodynamics: BEAST2
3. Configure run: substitution model, bootstrap replicates or posterior sampling, partitions
4. **BEAST2 divergence time gate (HARD STOP):** If researcher selects BEAST2 for time-calibrated analysis:
   - AI proposes candidate fossil/secondary calibration points from literature
   - **Pipeline PAUSES — human agent must review and confirm calibration points before proceeding**
   - Only after explicit human approval does the AI configure the BEAST2 XML and run
5. Assess support values: bootstrap thresholds (≥70 UFBoot, ≥95 standard), posterior probabilities (≥0.95)
6. Convergence check for Bayesian/BEAST2 runs (ESS ≥200, PSRF ≈1.0 via Tracer)
7. Sanity check tree topology against known biology

**QC gate:** Support values adequate, topology biologically plausible, convergence confirmed (Bayesian/BEAST2)  
**On QC failure:** Route to `phylo-debug`

**Report:** `reports/tree-inference_YYYY-MM-DD.md`
- Tool, version, and exact command used
- Model(s) applied
- Bootstrap/posterior support summary
- Convergence diagnostics (Bayesian/BEAST2 only)
- Calibration points used (BEAST2 only) — with human-approval record
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
- Data report exists, data type = raw reads (SRA) → `phylo-assemble`
- Data report exists, data type = assembled sequences → `phylo-alignment`
- Assembly report exists, no alignment → `phylo-alignment`
- Alignment done, no model → `phylo-model-selection`
- Model done, no tree → `phylo-tree-inference`
- Tree inference report shows BEAST2 pending human calibration approval → PAUSE, await human
- Tree done, no visualization → `phylo-visualization`
- Any report shows blocked/failed status → `phylo-debug`

---

## Case Study: Plant Systematics (Zingiberaceae)

### Scenario
A graduate student asks: "What are the phylogenetic relationships within Zingiberaceae (the ginger family), and can we resolve generic boundaries using the best available molecular data?"

### Dataset strategy (determined by `phylo-data-acquisition` survey)
- **Organism:** Zingiberaceae, family-level with genus-level resolution
- **Taxonomic scope:** Representative species across all major genera; outgroups from Musaceae/Cannaceae
- **Data survey outcome:** SRA contains genome skimming datasets for ~25 species across key genera; GenBank has assembled *matK*, *rbcL*, *trnL-F*, and ITS for ~90 species → two plans proposed:
  - **Plan A (broad):** 90 taxa × 4 markers (GenBank only; any taxon with ≥2 markers included)
  - **Plan B (marker-rich):** 25 taxa × full plastome via SRA assembly + nuclear ITS (GetOrganelle + GenBank supplement)
  - **Researcher selects both** → parallel analyses run under `reports/planA/` and `reports/planB/`

### Pipeline walkthrough
1. **Research design:** Literature review of Zingiberaceae systematics, landmark papers (e.g., Kress et al.), current generic concepts, known problem taxa
2. **Data acquisition:** SRA + GenBank survey using structured search terms; data strategy decision (assembly vs. direct download); LaTeX accession/run table + MD narrative
3. **Assembly:** GetOrganelle on genome skimming SRA runs → assembled plastomes; gene extraction via annotation; supplement with GenBank sequences for missing taxa
4. **Alignment:** MAFFT `--auto` per marker or whole-plastome alignment; partitioned concatenated matrix
5. **Model selection:** IQ-TREE ModelFinder with partitioned scheme per marker/region
6. **Tree inference:** IQ-TREE ML + 1000 ultrafast bootstrap; optionally BEAST2 for divergence times (calibration points from Zingiberales fossil record — **human approval required before run**)
7. **Visualization:** `ggtree` + `ggplot2` — annotated tree with bootstrap support, clade boxes, divergence time bars if BEAST2 used
8. **Debug scenario:** File naming mismatch — SRA-assembled sequences use run accession IDs (e.g., `SRR12345_plastome`) while GenBank sequences use species names; mixed naming breaks alignment → `phylo-debug` diagnoses, standardizes to `Genus_species_accession` format, re-runs

### Expected outputs
- `reports/research-design_2026-04-12.md`
- `reports/data-acquisition_2026-04-12.tex` + `reports/data-acquisition_2026-04-12.md`
- `reports/assembly_2026-04-12.md`
- `reports/alignment_2026-04-12.md`
- `reports/model-selection_2026-04-12.md`
- `reports/tree-inference_2026-04-12.md` *(includes BEAST2 calibration approval record if used)*
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
│   ├── assembly_YYYY-MM-DD.md
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
    ├── phylo-assemble/
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
- Population-level analyses (haplotype networks, PopART, DnaSP) not covered in v1 — future `phylo-population` module.
- `phylo-assemble` currently covers GetOrganelle, NOVOPlasty, HybPiper, and Trinity; additional assemblers (e.g., Captus, IOGA) could be added as the skill matures.
