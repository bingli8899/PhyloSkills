---
name: phylo-data-acquisition
description: Use after phylo-research-design is complete and a biological question, target taxa, and marker strategy are defined. Surveys GenBank and SRA for available data, builds a taxon-by-marker coverage matrix, proposes sampling plans with explicit trade-offs, and downloads approved data. Use when the researcher needs to find, evaluate, and obtain molecular sequence data for phylogenetic analysis.
---

# Acquiring Phylogenetic Sequence Data

## Overview

Survey available data before downloading anything. Build a coverage matrix, propose sampling plans that maximize both taxa and markers, get researcher approval, then download. Never download without an approved plan.

## Process

### Step 1 — Read the research design report

Load `reports/research-design_YYYY-MM-DD.md`. Extract:
- Target group name (at any taxonomic level: population → order)
- Required markers or genomic data strategy
- Taxonomic scope and approximate taxon count goal
- Outgroup taxa

### Step 2 — Data landscape survey

Search GenBank and SRA before touching any download command.

**GenBank search** (assembled sequences, gene records, plastomes, whole genomes):
```bash
# Entrez Direct — adapt terms per group and markers
esearch -db nuccore -query '"<group>" AND ("<marker1>" OR "<marker2>")' \
  | efetch -format docsum | xtract -pattern DocumentSummary \
    -element AccessionVersion Organism Title Length
```

**SRA search** (raw reads — genome skimming, WGS, HybSeq, transcriptomes):
```bash
esearch -db sra -query '"<group>" AND ("genome skimming" OR "WGS" \
  OR "target enrichment" OR "transcriptome")' \
  | efetch -format docsum | xtract -pattern DocumentSummary \
    -element Run@acc ScientificName LibraryStrategy LibrarySource
```

**Keyword templates** — adapt per group:
```
"<group name>" AND (phylogeny OR phylogenetic OR phylogenom* OR plastome
  OR "genome skimming" OR "target enrichment" OR transcriptome)
"<group name>" AND (rbcL OR matK OR trnL OR ITS OR psbA OR nrITS)
```

Supplement with database-specific searches as appropriate:
- **BOLD** — barcoding loci (COI, ITS, rbcL) with taxonomic filters
- **TreeBASE / Dryad** — published alignment matrices for reuse
- Any specialist database relevant to the organism group

### Step 3 — Build the coverage matrix

Compile a taxa × data-type availability table from the survey results:

| Taxon | rbcL | matK | ITS | Plastome (SRA) | HybSeq (SRA) |
|-------|------|------|-----|----------------|--------------|
| Sp. A | ✓ | ✓ | ✓ | — | ✓ |
| Sp. B | ✓ | — | ✓ | ✓ | — |
| ... | | | | | |

Compute summary statistics: total taxa found, per-marker taxon counts, data type breakdown, estimated missing data % at different cut-offs.

### Step 4 — Generate sampling plans

**Core principle: maximize both taxon count and marker count.** Plans represent natural trade-off thresholds in the data — not arbitrary points.

Identify thresholds where the trade-off changes sharply (e.g., dropping from 8 to 3 markers adds 40 more taxa — that is a natural plan boundary). Generate as many plans as meaningfully distinct thresholds exist.

For each plan present:
- Taxon count and which taxa are included/excluded
- Marker or data-type count
- Estimated missing data %
- Downstream route: `phylo-assemble` (raw SRA reads) or `phylo-alignment` (assembled sequences)
- Key trade-off: what is gained and lost vs. adjacent plans

Label plans descriptively from the data (e.g., *breadth-optimized*, *depth-optimized*, *balanced*, *all-inclusive*) — do not pre-assign labels.

### Step 5 — Human approves plan(s)

Present all plans to the researcher. They may select one or multiple.

**If multiple plans selected:** each runs as an independent parallel analysis. Create plan subdirectories:
```
reports/planA/
reports/planB/
```
All downstream module reports for each plan go into its subdirectory.

**Do not download anything until at least one plan is approved.**

### Step 6 — Estimate storage and choose download mode

Before downloading any SRA data, compare estimated dataset size against available server storage.

**Check available storage:**
```bash
df -h .        # available space in current working directory
df -h data/    # or wherever data will be stored
```

**Estimate SRA dataset sizes from metadata (before downloading):**
```bash
# Get file size for each SRA run from metadata
esearch -db sra -query "<SRR_ID>" | efetch -format runinfo \
  | cut -d',' -f7,10,11   # size_MB, LibraryStrategy, ScientificName

# Or use SRA toolkit
vdb-dump --info <SRR_ID> | grep -i "size"
```

**Decision — apply per plan:**

| Condition | Mode | Description |
|-----------|------|-------------|
| Available storage > estimated raw data × 1.5 | **Bulk mode** | Download all → assemble all → keep assemblies |
| Available storage ≤ estimated raw data × 1.5 | **Streaming mode** | Download one sample → assemble → keep assembly → delete raw → next |

The 1.5× buffer accounts for intermediate files created during assembly (GetOrganelle, HybPiper, etc. produce large temporary directories).

**Bulk mode** (storage is sufficient):
```bash
# Download all SRA runs first
for SRR_ID in $(cat sra_list.txt); do
  prefetch $SRR_ID -O data/raw/
  fasterq-dump data/raw/$SRR_ID/ -O data/raw/ --split-files
done
# Then hand off to phylo-assemble
```

**Streaming mode** (storage is limited — download → assemble → delete raw, one sample at a time):
```bash
for SRR_ID in $(cat sra_list.txt); do
  echo "Processing $SRR_ID..."

  # 1. Download
  prefetch $SRR_ID -O data/raw/
  fasterq-dump data/raw/$SRR_ID/ -O data/raw/ --split-files

  # 2. Assemble immediately (adapt command to assembly strategy)
  get_organelle_from_reads.py \
    -1 data/raw/${SRR_ID}_1.fastq -2 data/raw/${SRR_ID}_2.fastq \
    -F embplant_pt -o data/assembled/${SRR_ID}/ -t 8

  # 3. Verify assembly output exists before deleting raw data
  if [ -f "data/assembled/${SRR_ID}/*.fasta" ]; then
    rm -rf data/raw/$SRR_ID/
    rm -f data/raw/${SRR_ID}*.fastq
    echo "$SRR_ID: assembly complete, raw reads removed"
  else
    echo "$SRR_ID: assembly FAILED — raw reads retained for debugging"
  fi

  # 4. Log storage after each sample
  df -h . | tail -1
done
```

**Critical:** In streaming mode, only delete raw reads after confirming the assembly output file exists and is non-empty. Never delete before verifying. Log each deletion in the report.

**Assembled sequences from GenBank** (not affected by storage mode — these are small):
```bash
efetch -db nuccore -id <accession_list> -format fasta > sequences.fasta
```

**File naming convention — enforce strictly:**
```
<Genus>_<species>_<accession>_<marker>.fasta
# e.g. Zingiber_officinale_MN123456_matK.fasta
#      Curcuma_longa_SRR9876543_WGS.fastq
```

Mixed-source downloads are a common source of naming mismatches — standardize at download time, not later.

### Step 7 — QC check

- Sequence count matches expected plan count
- All files named consistently
- No duplicate accessions across sources
- FASTA headers parseable (no special characters breaking downstream tools)
- SRA downloads complete (file size > 0, paired files matched)

On any failure → route to `phylo-debug`.

## Reports

### `reports/data-acquisition_YYYY-MM-DD.md` — narrative

```markdown
# Data Acquisition Report
Date: YYYY-MM-DD

## Data Landscape Survey
[Summary of what was found in GenBank and SRA: total records,
 data types available, notable gaps or richly sampled lineages]

## Coverage Matrix
[Taxa × markers/data-type table]

## Sampling Plans Proposed
[Each plan: taxon count, marker count, missing data %, route, trade-off]

## Plan(s) Selected
[Which plan(s) the researcher chose and why]

## Search Queries Used
[Exact queries run against each database]

## Inclusion / Exclusion Decisions
[Sequences included or excluded with reasoning]

## Storage Assessment
- Available storage at start: [X GB]
- Estimated raw SRA data size: [X GB]
- Download mode selected: bulk / streaming — justification

## Streaming Mode Log (if applicable)
| SRR ID | Assembly output | Raw reads deleted | Notes |
|--------|----------------|-----------------|-------|

## Data Quality Notes
[Suspect sequences, misidentified taxa, truncated records flagged]

## Software Versions
[Entrez Direct version, SRA Toolkit version]

## Next Module
phylo-assemble (if SRA raw reads) / phylo-alignment (if assembled)
```

### `reports/data-acquisition_YYYY-MM-DD.tex` — LaTeX accession table

```latex
\begin{table}[h]
\caption{Sequence data used in this study}
\begin{tabular}{lllllll}
\hline
Taxon & Accession/Run & Database & Data type & Marker/Library & Length/Reads & Download date \\
\hline
Zingiber officinale & MN123456 & GenBank & Assembled & matK & 873 bp & 2026-04-13 \\
Curcuma longa       & SRR9876543 & SRA   & WGS       & —    & 4.2M reads & 2026-04-13 \\
\hline
\end{tabular}
\end{table}
```

One row per accession/run. Spans all plans — include a Plan column if multiple plans were run.

## Scripts

Pre-built scripts for this module are in `skills/phylo-data-acquisition/scripts/`. Load when needed:

| Script | Purpose |
|--------|---------|
| `download_genbank.sh` | Survey GenBank and download sequences per marker; enforces `Genus_species_accession_marker.fasta` naming |
| `download_sra.sh` | Storage-aware SRA download; auto-selects bulk vs. streaming mode; calls assembly script in streaming mode |

Usage examples:
```bash
# Survey only (no download)
bash skills/phylo-data-acquisition/scripts/download_genbank.sh \
  -g "Zingiberaceae" -m "matK,rbcL,ITS" -o data/genbank -s

# Download GenBank sequences (≤300 per marker)
bash skills/phylo-data-acquisition/scripts/download_genbank.sh \
  -g "Zingiberaceae" -m "matK,rbcL,ITS,psbA" -o data/genbank -n 300

# SRA download (auto mode — detects bulk vs. streaming)
bash skills/phylo-data-acquisition/scripts/download_sra.sh \
  -l sra_list.txt -o data/raw
```

## Common Mistakes

| Mistake | Fix |
|---------|-----|
| Downloading before surveying | Survey first — data landscape shapes which plan is even possible |
| Hardcoding a fixed number of taxa or markers | Let the coverage matrix reveal natural thresholds; don't pre-decide plan shapes |
| Mixing naming conventions from different databases | Standardize to `Genus_species_accession_marker` at download time |
| Treating missing data as a binary pass/fail | Estimate % per plan — some missing data is acceptable and expected |
| Downloading SRA reads without checking library strategy | Confirm WGS/genome skimming vs. amplicon vs. RNA-seq before routing to phylo-assemble |
| Forgetting to record Entrez Direct and SRA Toolkit versions | Run `edirect -version` and `fasterq-dump --version`; log both in report |
| Skipping storage estimation before bulk SRA download | Large WGS datasets can exceed hundreds of GB; always estimate first |
| Deleting raw reads before verifying assembly output | Check file exists and is non-empty before `rm`; a failed assembly with deleted reads cannot be recovered |
| Using bulk mode when storage is borderline | Apply the 1.5× buffer conservatively; streaming mode is safer and produces the same result |
