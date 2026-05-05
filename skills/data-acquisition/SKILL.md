---
name: data-acquisition
description: Use after research-design is complete and a biological question, target taxa, and marker strategy are defined. Surveys GenBank and SRA for available data, builds a taxon-by-marker coverage matrix, proposes sampling plans with explicit trade-offs, and downloads approved data. Use when the researcher needs to find, evaluate, and obtain molecular sequence data for phylogenetic analysis.
---

# Acquiring Phylogenetic Sequence Data

## Overview

> **General-purpose skill — do not hard-code taxon names or marker sets.**
> This skill applies to any taxonomic group and any marker strategy. All taxon names, marker lists, accession numbers, and path examples in this document use generic placeholders such as `<TaxonName>`, `<marker1>`, and `<plan>`. Replace these with actual values for each study. When modifying this skill, keep all examples generic so the skill remains reusable across projects.

Survey available data before downloading anything. Build a coverage matrix, propose sampling plans that maximize both taxa and markers, get researcher approval, then download. Never download without an approved plan.

---

## CRITICAL tips: Read Before Any Download

### Marker search always returns complete plastomes — verify before accepting sequences
This section is applied to when **plastid markers** will be used and downloaded. 
If the plan is to download SRA raw reads or complete plastomes, this section does not apply. 

When searching GenBank for individual plastid markers (matK, rbcL, atpB, etc.), the query
`"marker"[All Fields]` also matches **complete chloroplast genome records** because those
records contain the gene. `efetch -format fasta` returns the **entire sequence** of whatever
record is found — so a search for "atpB" will return complete 150 kb plastomes alongside
genuine 1.5 kb atpB amplicons.

**You must check sequence lengths after any GenBank marker download.** A file > 20 kb in a
marker directory is a complete plastome, not a marker sequence. This will corrupt alignment.

`download_genbank.sh` handles this automatically with a post-download size check (see §Scripts).
If you see oversized files in existing data, run `fix_plastome_contamination.py` to fix them
(see §Known Issues).


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

This is an example only applied to plastid markers. The actual matrix should include all data types, depending on specific research plan decided by human, relevant to the research design (e.g., nuclear genes, transcriptomes, etc.) and all markers specified in the design. 
| Taxon | rbcL | matK | trnL | psbA-trnH | Plastome (SRA) |
|-------|------|------|------|-----------|----------------|
| Sp. A | ✓ | ✓ | ✓ | ✓ | — |
| Sp. B | ✓ | — | ✓ | — | ✓ |
| ... | | | | | |
The above is an example. 

Compute summary statistics: total taxa found, per-marker taxon counts, data type breakdown, estimated missing data % at different cut-offs.

This matrix could be flexible. For example, if the research design specifies a core set of 5 markers but also allows optional additional markers, the matrix could have columns for "core marker count" and "total marker count" to reflect this. Or if the design allows for either plastid markers or transcriptomes, the matrix could have separate columns for "plastid marker count" and "transcriptome availability". The key is to capture the relevant dimensions of data availability that will inform sampling plans.
If the research design focuses solely on full plastome, then we do not need to show each individual marker in the matrix. Instead, we can have a column for "Plastome (SRA)" that indicates whether a complete plastome assembly is available from SRA data, which would be the preferred data type for phylogenetic analysis in that case. 
If the research design includes both plastid markers and nuclear genes, the matrix should have separate columns for each data type, and the sampling plans can be based on different combinations of these data types. The matrix should be designed to reflect the specific requirements and flexibility of the research design, while providing a clear overview of data availability across taxa and markers/data types. 
If the research design includes a large number of nuclear markers (>100 markers), it may be helpful to summarize the marker availability in terms of "core marker count" (number of required markers available) and "total marker count" (number of all markers available) to simplify the matrix and focus on the key trade-offs for sampling plans. 

### Step 4 — Generate sampling plans

**Core principle: maximize both taxon count and marker count.** Plans represent natural trade-off thresholds in the data — not arbitrary points.

Identify thresholds where the trade-off changes sharply (e.g., dropping from 8 to 3 markers adds 40 more taxa — that is a natural plan boundary). Generate as many plans as meaningfully distinct thresholds exist.

For each plan present:
- Taxon count and which taxa are included/excluded
- Marker or data-type count
- Estimated missing data %
- Downstream route: `assembly` (raw SRA reads) or `alignment` (assembled sequences)
- Key trade-off: what is gained and lost vs. adjacent plans

Label plans descriptively from the data (e.g., *breadth-optimized*, *depth-optimized*, *balanced*, *all-inclusive*) — do not pre-assign labels.

### Step 5 — Human approves plan(s)

Present all plans to the researcher. They may select one or multiple.
Human approval is required. Do not proceed to download data until at least one plan is approved. 

**If multiple plans selected:** each runs as an independent parallel analysis. Create plan subdirectories:
```
reports/planA/
reports/planB/
```
All downstream module reports for each plan go into its subdirectory.

**Do not download anything until at least one plan is approved.**

### Step 6 — Download GenBank sequences

Use `download_genbank.sh` for all GenBank marker downloads. This script handles the
complete-plastome contamination problem automatically (see §Known Issues #1).

```bash
# Survey only (no download — run this first)
bash scripts/data/download_genbank.sh \
  -g "<TaxonName>" \
  -m "<marker1>,<marker2>,..." \
  -o data/<plan>/genbank -s

# Download (≤2000 per marker)
bash scripts/data/download_genbank.sh \
  -g "<TaxonName>" \
  -m "<marker1>,<marker2>,..." \
  -o data/<plan>/genbank -n 2000
```

After download, immediately run `select_best_accession.py` to deduplicate to one sequence
per species per marker (see §Scripts).

### Step 6b — Estimate storage and choose SRA download mode

Before downloading any SRA data, compare estimated dataset size against available server storage.

**Check available storage:**
```bash
df -h .        # available space in current working directory
df -h data/    # or wherever data will be stored
```

**Estimate SRA dataset sizes from metadata (before downloading):**
```bash
esearch -db sra -query "<SRR_ID>" | efetch -format runinfo \
  | cut -d',' -f7,10,11   # size_MB, LibraryStrategy, ScientificName

vdb-dump --info <SRR_ID> | grep -i "size"
```

**Decision — apply per plan:**

| Condition | Mode | Description |
|-----------|------|-------------|
| Available storage > estimated raw data × 1.5 | **Bulk mode** | Download all → assemble all → keep assemblies |
| Available storage ≤ estimated raw data × 1.5 | **Streaming mode** | Download one sample → assemble → keep assembly → delete raw → next |

The 1.5× buffer accounts for intermediate files created during assembly (GetOrganelle, HybPiper, etc.).

**Streaming mode** (storage is limited — download → assemble → delete raw, one at a time):
```bash
bash scripts/data/download_sra.sh \
  -l data/<plan>/sra_runlist.txt \
  -o data/<plan>/plastomes \
  -a scripts/data/assemble_and_extract.sh
```

`download_sra.sh` selects bulk vs. streaming automatically. In streaming mode it calls
`assemble_and_extract.sh` per sample (see §Scripts).

**Critical:** In streaming mode, only delete raw reads after confirming the assembly output
file exists and is non-empty. Never delete before verifying.

**File naming convention — enforce strictly:**
```
<Genus>_<species>_<accession>_<marker>.fasta
# e.g. Zingiber_officinale_MN123456_matK.fasta
#      Curcuma_longa_SRR9876543_WGS.fastq
```

### Step 6c — Fill SRA gaps (species in SRA absent from GenBank markers)

After the GenBank download and deduplication, identify species that are in SRA but have no
sequences in the GenBank marker dataset. These are "gap" species that can be filled by
assembling their plastomes from SRA reads.

sra meta data is important to be saved for future usage. The SRA metadata TSV should contain at least the following columns:
- Run@acc (SRA run accession)
- ScientificName (species name)
- BioProject (project accession) 

Below script assumes that planB is the GenBank marker dataset and planA is the SRA dataset. Adjust paths if different. Plan name is arbitrary — just be consistent across scripts. 
```bash
python scripts/data/find_sra_gaps.py \
    --genbank_dedup data/planB/genbank_dedup \
    --sra_tsv       data/planA/sra_metadata.tsv \
    --markers       matK,rbcL,trnL,psbA-trnH,rpoB,rpoC1,atpB,ndhF,ycf1,ycf2 \
    --output        data/planA/sra_gaps.tsv
```

Review the output TSV before proceeding. Exclude:
- Species with extremely large SRA runs (> 10 GB, e.g. high-coverage WGS not needed for plastome assembly)
- Unidentified/undescribed species you choose to exclude

Then download + assemble + extract + merge with:
```bash
bash scripts/data/download_sra.sh \
    -l data/planA/sra_gaps_runlist_filtered.txt \
    -o data/planA/plastomes \
    -a scripts/data/assemble_and_extract.sh
```
Again, the above plan names are arbitrary — just be consistent across scripts. The key is that the gap species are assembled and extracted in planA, then merged into planB's genbank_dedup/ by `merge_sra_cds_to_planB.py` (see §Scripts). 

`assemble_and_extract.sh` calls `extract_markers_blast.py` and `merge_sra_cds_to_planB.py`
automatically after each successful assembly.

### Step 7 — QC check

QC check is mandatory after every download and before alignment. Do not skip. This is the last chance to catch problems before they propagate to alignment and tree inference. **QC check reports will be reported to human for approval before proceeding to alignment.**  

After all downloads are complete:
- Sequence count per marker matches expected plan count
- All files named consistently (`Genus_species_accession_marker.fasta`)
- No duplicate accessions across sources
- FASTA headers parseable (no special characters breaking downstream tools)
- **No complete plastome contamination** — run the detection command below
- SRA downloads complete (file size > 0, paired files matched)

**Contamination audit (run after every GenBank download):**
```bash
# Report any FASTA file > 20 kb in marker directories — these are complete plastomes
for f in data/planB/genbank_dedup/*/*.fasta; do
    size=$(awk '/^>/{next}{len+=length($0)}END{print len+0}' "$f")
    [[ "$size" -gt 20000 ]] && echo "$size bp  $f"
done

# Grep FASTA headers for complete genome language
grep -rl "complete genome" data/planB/genbank_dedup/
```
Again file name should be arbitrary — just be consistent across scripts. The key is that the contamination check is run after every GenBank download and before alignment, and any contaminated files are identified and fixed before proceeding. 

If contamination is found, run `fix_plastome_contamination.py` (see §Known Issues #1).

**Sequence header normalization:**

FASTA headers from GenBank often contain full titles like:
```
>PP542015.1 Aframomum alboviolaceum chloroplast, complete genome
```
These must be normalized to the pipeline convention before alignment:
```
>Aframomum_alboviolaceum_PP542015.1
```
The file name carries the marker (`_atpB.fasta`) if maker is used, so the header does not need it. However, in general, using the species + accession form is sufficient and consistent across sources.

Headers already generated by `download_genbank.sh`, `extract_markers_blast.py`, and
`fix_plastome_contamination.py` follow this convention. For any file that does not, use:
```bash
python scripts/utils/change_header_name.py --mapping header_map.tsv --fasta file.fasta
```

On any failure → route to `debug`.

---

## Known Issues and Fixes

This section documents every data-quality problem encountered in the Zingiberaceae pipeline
and the exact script to fix it. Run diagnostics before alignment.

---

### Issue #1 — Complete plastome sequences in marker directories

**Symptom:** Files > 20 kb in `data/planB/genbank_dedup/<marker>/`. FASTA header contains
"complete genome" or "chloroplast genome". These sequences are 100–170 kb full plastomes,
not the intended 1–4 kb marker amplicons.

**Root cause:** `esearch … "marker"[All Fields]` matches any GenBank record annotated with
that gene, including complete plastomes. `efetch -format fasta` returns the full record.

**Prevention (automatic):** `download_genbank.sh` now checks sequence length after each
download and re-fetches `.gb` + extracts the gene if oversized. No manual action needed for
new downloads.

**Fix for existing data — batch:**
```bash
# Scans ALL marker directories, downloads .gb for each contaminated accession,
# extracts the correct gene, replaces the file. Safe to re-run (idempotent).
$HOME/miniconda3/bin/python3 scripts/data/fix_plastome_contamination.py \
    --genbank_dedup data/planB/genbank_dedup \
    --markers       matK,rbcL,trnL,psbA-trnH,rpoB,rpoC1,atpB,ndhF,ycf1,ycf2

# Dry run first to see scope:
$HOME/miniconda3/bin/python3 scripts/data/fix_plastome_contamination.py \
    --genbank_dedup data/planB/genbank_dedup \
    --markers       matK,rbcL,trnL,psbA-trnH,rpoB,rpoC1,atpB,ndhF,ycf1,ycf2 \
    --dry_run
```

**Fix for a single file:**
```bash
efetch -db nuccore -id "PP542015.1" -format gb > /tmp/PP542015.1.gb
$HOME/miniconda3/bin/python3 scripts/data/extract_gene_from_gb.py \
    --gb     /tmp/PP542015.1.gb \
    --gene   atpB \
    --out    data/planB/genbank_dedup/atpB/Aframomum_alboviolaceum_PP542015.1_atpB.fasta \
    --header "Aframomum_alboviolaceum_PP542015.1"
```

---

### Issue #2 — Duplicate accessions: multiple records per species per marker

**Symptom:** Many FASTA files for the same species in a marker directory (e.g. 5 different
matK sequences for *Alpinia galanga*). The alignment step needs exactly one sequence per
species to avoid redundant taxa inflating tree support.

**Root cause:** GenBank has multiple submissions for the same species and marker. The download
cap (`-n 2000`) downloads all of them.

**Fix:**
```bash
# Deduplicates to one best accession per species per marker.
# Criteria: longest sequence, then most recent accession, then NC_ preferred.
# Writes kept files to data/planB/genbank_dedup/ and a selection report.
$HOME/miniconda3/bin/python3 scripts/data/select_best_accession.py \
    --input      data/planB/genbank \
    --output     data/planB/genbank_dedup \
    --markers    matK,rbcL,trnL,psbA-trnH,rpoB,rpoC1,atpB,ndhF,ycf1,ycf2 \
    --report     reports/planB/accession_selection.tsv \
    --provenance results/planB/provenance
```

---

### Issue #3 — SRA gap species: in SRA but absent from GenBank marker data

**Symptom:** Species list from SRA metadata contains taxa with no sequences in any
`genbank_dedup/<marker>/` directory. These would be completely missing from the tree.

**Diagnosis:**
```bash
$HOME/miniconda3/bin/python3 scripts/data/find_sra_gaps.py \
    --genbank_dedup data/planB/genbank_dedup \
    --sra_tsv       data/planA/sra_metadata.tsv \
    --markers       matK,rbcL,trnL,psbA-trnH,rpoB,rpoC1,atpB,ndhF,ycf1,ycf2 \
    --output        data/planA/sra_gaps.tsv
```

**Fix:** Download, assemble, extract, and merge for gap species (see Step 6c above).
After extraction via `assemble_and_extract.sh`, the sequences land in `data/planA/cds/`
and are merged to `data/planB/genbank_dedup/` by `merge_sra_cds_to_planB.py`.

---

### Issue #4 — SRA assembly too fragmented to extract markers

**Symptom:** GetOrganelle produces a scaffold assembly < 50 kb for a sample. BLAST extraction
finds few or no markers.

**Root cause:** Low sequencing depth, contamination, or poor library quality. GetOrganelle
cannot complete the circular plastome with insufficient coverage.

**Diagnosis:**
```bash
# Check assembly lengths for all SRA assemblies
find data/planA/plastomes -name "*.path_sequence.fasta" ! -name "*selected_graph*" | \
  while read f; do
    len=$(awk '/^>/{next}{len+=length($0)}END{print len}' "$f")
    echo "$len  $f"
  done | sort -n | head -20
```

**Fix:** Skip these samples — do not extract or merge. `extract_markers_blast.py` enforces
`--min_assembly 50000` by default and will print "Assembly too short … skipping all markers."
No action needed if using the standard scripts. The sample is simply absent from the tree
(acceptable missing data).

---

### Issue #5 — BLAST extraction fails for specific markers (matK, ycf1, ycf2)

**Symptom:** `extract_markers_blast.py` reports "no hit" for matK even on a complete
(>150 kb) assembly; or reports "no hit" for ycf1/ycf2 when the gene is clearly present.

**Root causes and fixes — already applied in the current script:**

| Problem | Wrong setting | Correct setting |
|---------|--------------|-----------------|
| matK not found at ~75% identity across genera | `-word_size 11` (default) | `-word_size 7` in blastn call |
| matK empty output with max_hsps | `-max_hsps 1` (causes empty output bug) | Omit `-max_hsps` flag entirely |
| ycf1/ycf2 fails min_cov check | `min_cov=50%` global | `min_cov=5%` for ycf1/ycf2 via `MARKER_MIN_COV` dict |

These are **already fixed** in `scripts/data/extract_markers_blast.py`. Do not revert them.

**If you add new BLAST calls elsewhere**, remember:
- Always use `-word_size 7` for cross-genus plant plastid marker searches
- Never add `-max_hsps 1` to blastn when extracting plastid markers
- Apply per-marker min coverage thresholds for large genes (ycf1, ycf2)

---

### Issue #6 — Wrong files merged from a previous extraction run

**Symptom:** SRA-extracted marker files in `data/planB/genbank_dedup/` have a header but
empty or incorrect sequence (e.g. from a run with wrong BLAST reference).

**Diagnosis:**
```bash
# Find zero-length or near-empty FASTA files
for marker in matK rbcL trnL psbA-trnH rpoB rpoC1 atpB ndhF ycf1 ycf2; do
    dir=data/planB/genbank_dedup/$marker
    for f in $dir/*_SRR* $dir/*_ERR* $dir/*_DRR* 2>/dev/null; do
        [[ -f "$f" ]] || continue
        len=$(awk '/^>/{next}{len+=length($0)}END{print len+0}' "$f")
        [[ "$len" -eq 0 ]] && echo "EMPTY: $f"
        [[ "$len" -gt 0 && "$len" -lt 100 ]] && echo "SUSPECT ($len bp): $f"
    done
done
```

**Fix:**
```bash
# 1. Remove wrong SRA files from genbank_dedup
for marker in matK rbcL trnL psbA-trnH rpoB rpoC1 atpB ndhF ycf1 ycf2; do
    rm -f data/planB/genbank_dedup/$marker/*_SRR*.fasta \
          data/planB/genbank_dedup/$marker/*_ERR*.fasta \
          data/planB/genbank_dedup/$marker/*_DRR*.fasta
done

# 2. Clear the planA/cds output directory to force re-extraction
rm -rf data/planA/cds/

# 3. Re-run extraction for all SRA assemblies
bash scripts/extract_all_sra_markers.sh
```

---

### Issue #7 — Incorrect PLANN annotation on scaffold assemblies

**Symptom:** PLANN annotates only 40–50 genes on a scaffold assembly (expected > 100 for a
complete plastome). The `.tbl` feature table is nearly empty.

**Root cause:** PLANN requires a **complete, circular plastome**. Scaffold assemblies from
GetOrganelle (produced from low-coverage genome-skimming) are fragmented; PLANN cannot
align reference gene coordinates to discontiguous scaffolds.

**Fix:** Do not use PLANN on scaffold assemblies. Use BLAST-based extraction instead:
```bash
$HOME/miniconda3/bin/python3 scripts/data/extract_markers_blast.py \
    --plastome     <assembly.fasta> \
    --ref_dir      data/planB/genbank_dedup \
    --ref_plastome data/planA/reference/Zingiber_officinale_NC_037455.1.fasta \
    --markers      matK,rbcL,trnL,psbA-trnH,rpoB,rpoC1,atpB,ndhF,ycf1,ycf2 \
    --species      <Genus_species> \
    --run_id       <SRR_ID> \
    --output       data/planA/cds \
    --threads      4
```

PLANN is only appropriate when GetOrganelle produces a **complete circular assembly**
(file named `*.complete.graph1.1.path_sequence.fasta`, length ≥ 145 kb for Zingiberaceae).

---

## Scripts Reference

All scripts are in `scripts/data/` unless noted. Use `$HOME/miniconda3/bin/python3` for
all Python scripts (has Biopython). BASH scripts call this Python automatically.

### Download scripts

| Script | Purpose | Key arguments |
|--------|---------|---------------|
| `download_genbank.sh` | Survey + download GenBank marker sequences; size-checks each download; auto-extracts gene from complete plastomes via `extract_gene_from_gb.py` | `-g <group>` `-m <markers>` `-o <outdir>` `-n <max>` `-s` (survey-only) |
| `download_sra.sh` | Storage-aware SRA download; auto-selects bulk vs. streaming; calls assembly script in streaming mode | `-l <runlist>` `-o <outdir>` `-a <assembly_script>` |

### Deduplication and gap scripts

| Script | Purpose | Key arguments |
|--------|---------|---------------|
| `select_best_accession.py` | Deduplicate to one best sequence per species per marker; writes kept files to `genbank_dedup/` | `--input` `--output` `--markers` `--report` `--provenance` |
| `find_sra_gaps.py` | Find species in SRA metadata that have no GenBank marker sequences; outputs gap species + best SRA run per species | `--genbank_dedup` `--sra_tsv` `--markers` `--output` |

### Extraction and fix scripts

| Script | Purpose | Key arguments |
|--------|---------|---------------|
| `extract_markers_blast.py` | BLAST-based marker extraction from an assembled plastome FASTA; handles all 10 Plan B markers; skips assemblies < 50 kb | `--plastome` `--ref_dir` `--ref_plastome` `--markers` `--species` `--run_id` `--output` `--threads` `--min_assembly` |
| `extract_gene_from_gb.py` | Extract one marker gene from a GenBank annotation file; used by `download_genbank.sh` and for manual one-off fixes | `--gb` `--gene` `--out` `--header` |
| `fix_plastome_contamination.py` | **Batch repair:** scans all marker dirs for complete-plastome files (>20 kb), downloads `.gb` for each, extracts correct genes, replaces files in-place | `--genbank_dedup` `--markers` `[--dry_run]` |
| `merge_sra_cds_to_planB.py` | Copy SRA-extracted marker sequences from `planA/cds/` into `planB/genbank_dedup/`; idempotent (skips existing files) | `--cds_dir` `--planB_dir` `--markers` `[--dry_run]` `--provenance` |

### Per-sample assembly wrapper

| Script | Purpose | Called by |
|--------|---------|-----------|
| `assemble_and_extract.sh` | Per-SRA-run wrapper: GetOrganelle → `extract_markers_blast.py` → `merge_sra_cds_to_planB.py`; prints assembly dir as last line of stdout for `download_sra.sh` delete-raw logic | `download_sra.sh` (streaming mode) |

### Batch re-extraction (project-specific)

| Script | Location | Purpose |
|--------|----------|---------|
| `extract_all_sra_markers.sh` | `Zingiberaceae/scripts/` | Re-run BLAST extraction for all 39 SRA assemblies; picks best assembly per species; merges into genbank_dedup |

---

## Reports

**Mandatory:** Every log file generated during this module must be listed with its exact path in the report so the researcher can monitor background processes and audit what ran.

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

## Log Files Generated
[List every log file created during this module with its exact path so the researcher
 can monitor background downloads and diagnose failures]
[Examples:]
[  data/planB/genbank_download_2026-04-19.log  (download_genbank.sh stdout/stderr)]
[  data/planA/sra_download_2026-04-19.log      (download_sra.sh per-sample log)]
[  data/planA/plastomes/SRR<ID>/get_organelle.log  (GetOrganelle per-sample log)]

## Next Module
assembly (if SRA raw reads) / alignment (if assembled)
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

---

## Common Mistakes

| Mistake | Consequence | Fix |
|---------|-------------|-----|
| Using `"marker"[All Fields]` without post-download size check | Complete plastomes (~150 kb) in marker directories corrupt alignment | Run `fix_plastome_contamination.py`; or use `download_genbank.sh` which checks automatically |
| Downloading before surveying | Data landscape shapes which plan is even possible | Survey first with `-s` flag |
| Hardcoding a fixed number of taxa or markers | Misses natural data thresholds | Let coverage matrix reveal thresholds |
| Mixing naming conventions from different databases | Header parsing failures in alignment | Standardize to `Genus_species_accession_marker` at download time |
| Treating missing data as a binary pass/fail | Plans are over- or under-filtered | Estimate % per plan — 20–40% missing is acceptable for partitioned ML |
| Downloading SRA reads without checking library strategy | RNA-seq or amplicon reads sent to wrong assembly workflow | Confirm WGS/genome-skimming vs. amplicon vs. RNA-seq before routing |
| Forgetting to record tool versions | Methods section incomplete | Run `edirect -version` and `fasterq-dump --version`; log in report |
| Skipping storage estimation before bulk SRA download | Runs out of disk mid-download | Always estimate first; use streaming mode if borderline |
| Deleting raw reads before verifying assembly output | Unrecoverable data loss | Check file exists and is non-empty before `rm` |
| Using `get_organelle` on fragmented assemblies then running PLANN | PLANN returns only 40–50/130 genes from scaffold assemblies | Use BLAST extraction (`extract_markers_blast.py`) for scaffold assemblies |
| Adding `-max_hsps 1` to blastn | Causes empty output for matK and other markers | Omit `-max_hsps`; already fixed in `extract_markers_blast.py` |
| Using default blastn `word_size 11` for cross-genus searches | Misses matK (~75% identity); no hits reported | Use `-word_size 7`; already set in `extract_markers_blast.py` |
| Applying global `min_cov=50%` to ycf1/ycf2 | Large genes fail coverage check even when partially recovered | Use `MARKER_MIN_COV = {"ycf1": 5.0, "ycf2": 5.0}`; already set in script |
