---
name: phylo-assemble
description: Use when phylo-data-acquisition has identified raw SRA reads (genome skimming, WGS, target enrichment/HybSeq, or transcriptomes) that need assembly or marker extraction before alignment. Runs between phylo-data-acquisition and phylo-alignment. Use when assembled FASTA sequences are not directly available and must be generated from raw reads using reference-guided or de novo approaches.
---

# Assembling Markers and Genomes from Raw Reads

## Overview

Convert raw SRA reads into assembled sequences ready for alignment. Strategy depends entirely on the data type — choose the right tool for the data, not the most familiar one.

## Strategy Selection

| Data type (from SRA survey) | Goal | Primary tool | Alternative |
|-----------------------------|------|-------------|-------------|
| Genome skimming / WGS | Full plastome | GetOrganelle | NOVOPlasty |
| Genome skimming / WGS | Nuclear markers from reference | BWA + SAMtools consensus | Bowtie2 + SAMtools |
| Target enrichment (HybSeq, Angiosperms353, custom baits) | Target gene recovery | HybPiper | Captus |
| Transcriptome (RNA-seq) | Gene extraction | Trinity → BLAST+ extraction | — |
| Mixed (multiple data types in same plan) | Per-sample strategy | Apply per row | — |

Read `reports/[planX/]data-acquisition_YYYY-MM-DD.md` to confirm data type per sample before choosing strategy.

## Step 1 — Select and obtain a reference

Reference quality determines assembly quality. Assess availability in this order:

1. **Published whole plastome / genome for the same genus** — use if available; most reliable
2. **Published reference for the same family** — acceptable; flag divergence to researcher
3. **Published reference for the same order** — use with caution; expect lower recovery rates
4. **No close reference** — for plastome: use GetOrganelle with a distant seed; for target enrichment: use probe sequences directly as reference; alert researcher

For GetOrganelle: reference databases are built-in; specify the correct database flag (`-F embplant_pt` for plant plastome, `-F embplant_mt` for mitochondria, `-F animal_mt`, etc.).

For HybPiper: the target file is the probe/bait sequence set (e.g., Angiosperms353 `mega353.fasta`, or a custom bait file from the literature).

## Step 2 — Run assembly per sample

Apply the selected strategy to every sample in the approved plan. Process samples in batch where possible.

### Plastome assembly (GetOrganelle)

```bash
get_organelle_from_reads.py \
  -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  -F embplant_pt \
  -o output/sample_plastome/ \
  -t 8
```

Key flags: `-R` (rounds, increase if assembly fails), `-k` (kmer sizes), `-w` (word size for seed).  
Output: circular plastome FASTA + assembly graph. Inspect graph in Bandage if assembly is fragmented.

**After plastome assembly — extract genes:**
```bash
# Annotate with PGA or GeSeq (web), then extract CDS/gene regions
# Or use mfannot for a command-line option
python GetOrganelleAnnotation/GetAnnotation.py ...
```

### Plastome assembly (NOVOPlasty — alternative)

Seed-and-extend; requires a short seed sequence (e.g., rbcL from a close relative):
```bash
perl NOVOPlasty.pl -c config.txt
# config.txt: set Seed_input, Genome_range, Read_length, Insert_size
```

### Nuclear marker extraction (reference-guided)

```bash
# Index reference
bwa index reference_markers.fasta

# Map and extract consensus per sample
bwa mem reference_markers.fasta sample_R1.fastq.gz sample_R2.fastq.gz \
  | samtools sort -o sample.bam
samtools index sample.bam
samtools consensus -f fasta sample.bam > sample_markers_consensus.fasta
```

Minimum mapping depth for a usable consensus: ≥5×; flag samples below 10× as low-confidence.

### Target enrichment assembly (HybPiper)

```bash
# Assemble all samples
hybpiper assemble -t_dna mega353.fasta \
  -r sample_R1.fastq.gz sample_R2.fastq.gz \
  --prefix sample --cpu 8

# Retrieve sequences across all samples
hybpiper retrieve_sequences dna -t_dna mega353.fasta \
  --sample_names namelist.txt
```

Review recovery rates with:
```bash
hybpiper stats -t_dna mega353.fasta gene --sample_names namelist.txt
hybpiper recovery_heatmap seq_lengths.tsv
```

Acceptable gene recovery: ≥50% of targets per sample. Flag samples below this threshold.

### Transcriptome gene extraction (Trinity + BLAST+)

```bash
# Assemble transcriptome
Trinity --seqType fq --left sample_R1.fastq.gz --right sample_R2.fastq.gz \
  --max_memory 50G --CPU 8 --output trinity_out/

# Extract target genes by BLAST
makeblastdb -in reference_genes.fasta -dbtype nucl
blastn -query trinity_out/Trinity.fasta -db reference_genes.fasta \
  -outfmt 6 -evalue 1e-10 -out blast_results.txt
# Filter top hits and extract sequences
```

## Step 3 — Collect and organize outputs

After assembly, consolidate per-marker FASTA files across all samples:

```
data/assembled/
  matK_all_samples.fasta        # one file per marker, all samples
  rbcL_all_samples.fasta
  plastome_complete/
    Genus_species_accession.fasta
```

Apply the naming convention from `phylo-data-acquisition`:
```
>Genus_species_accession_marker
```

Rename FASTA headers immediately after assembly — do not carry through tool-generated headers.

## Step 4 — QC

| Check | Threshold | Action on failure |
|-------|-----------|-------------------|
| Plastome completeness | ≥90% of expected length | Retry with adjusted `-R`/`-k`; flag to researcher |
| Mean read depth (reference-guided) | ≥10× recommended; ≥5× minimum | Flag low-coverage samples; consider excluding |
| HybPiper gene recovery rate | ≥50% targets per sample | Flag low-recovery samples |
| Assembly graph (GetOrganelle) | Single circular path | Inspect in Bandage if fragmented |
| FASTA header consistency | All headers match naming convention | Rename before proceeding |
| No empty output files | All expected files >0 bytes | Debug with `phylo-debug` |

On any failure → route to `phylo-debug` before proceeding.

## Report

Write to `reports/[planX/]assembly_YYYY-MM-DD.md`:

```markdown
# Assembly Report
Date: YYYY-MM-DD
Plan: [planA / planB / ...]

## Assembly Strategy
- Data type: [genome skimming / WGS / HybSeq / transcriptome]
- Tool: [GetOrganelle / HybPiper / BWA+SAMtools / Trinity]
- Reference used: [name, accession, source, divergence notes]

## Per-Sample Results
| Sample | Accession/Run | Completeness / Recovery | Mean depth | Notes |
|--------|--------------|------------------------|------------|-------|

## Markers / Regions Recovered
[List of genes/regions successfully assembled across samples]

## Failed / Low-Quality Samples
[Samples excluded or flagged, with reason]

## Output Files
[Paths to consolidated per-marker FASTA files]

## Software Versions
| Tool | Version | Source | Install date |
|------|---------|--------|-------------|

## Next Module
phylo-alignment
```

## Scripts

Pre-built scripts for this module are in `skills/phylo-assemble/scripts/`. Load when needed:

| Script | Purpose |
|--------|---------|
| `assemble_plastome_getorganelle.sh` | De novo plastome assembly with GetOrganelle; QC checks assembly length |
| `assemble_plastome_bwa.sh` | Reference-guided plastome assembly: BWA → SAMtools → BCFtools consensus; strict or majority mode |
| `run_hybpiper.sh` | Full HybPiper v2 pipeline: assemble → stats → retrieve_sequences → optional paralog_retriever |
| `annotate_plastome.sh` | Plastome annotation with PLANN or chloe (auto-detected); configurable paths |
| `extract_cds.py` | Extract CDS from GenBank annotations; handles multi-exon genes; outputs per-gene FASTA + partition |

Utility scripts (in `skills/utils/`):

| Script | Purpose |
|--------|---------|
| `change_header_name.py` | Rename FASTA headers from a two-column TSV mapping |
| `remove_N_fasta.py` | Remove sequences with excessive N or gap characters |
| `revise_hybpiper_sequences.py` | Post-process HybPiper output: rename headers, filter short seqs, flag low-recovery samples |

Usage examples:
```bash
# GetOrganelle plastome assembly (land plant)
bash skills/phylo-assemble/scripts/assemble_plastome_getorganelle.sh \
  -1 SRR123_1.fastq -2 SRR123_2.fastq -o assemblies -s Zingiber_officinale

# BWA-based assembly (reference-guided)
bash skills/phylo-assemble/scripts/assemble_plastome_bwa.sh \
  -1 SRR123_1.fastq -2 SRR123_2.fastq \
  -r reference_plastome.fasta -o assemblies -s Zingiber_officinale -c strict

# HybPiper pipeline
bash skills/phylo-assemble/scripts/run_hybpiper.sh \
  -r Angiosperms353_targetfile.fasta -s sample_list.txt \
  -d data/raw -o hybpiper_output

# Extract CDS from annotations
python skills/phylo-assemble/scripts/extract_cds.py \
  --input annotations/ --output data/cds/ --concatenate
```

## Common Mistakes

| Mistake | Fix |
|---------|-----|
| Using the wrong `-F` database flag in GetOrganelle | Check organism kingdom: `embplant_pt` for land plant plastome, `animal_mt` for animal mitochondria |
| Keeping tool-generated FASTA headers | Rename to `Genus_species_accession_marker` immediately after assembly |
| Accepting fragmented plastome assemblies without inspection | Open the assembly graph in Bandage; fragmentation often means low coverage or contamination |
| Running HybPiper without checking recovery heatmap | Low per-sample recovery is invisible until you visualize it; always run `recovery_heatmap` |
| Mixing assembled and directly-downloaded sequences with inconsistent headers | Standardize headers before merging into per-marker files |
| Skipping depth check for reference-guided assemblies | Consensuses from <5× depth are unreliable; flag or exclude |
