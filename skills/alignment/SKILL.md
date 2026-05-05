---
name: alignment
description: Use after assembled or downloaded sequences are organized and ready for alignment. Selects alignment strategy based on dataset size and sequence type, runs per-marker alignments, optionally trims, and guides the concatenation vs. coalescent decision for multi-locus datasets. Use when raw FASTA files need to be aligned before model selection and tree inference.
---

# Aligning Phylogenetic Sequences

## Overview

> **General-purpose skill — do not hard-code taxon names or marker sets.**
> This skill applies to any taxonomic group and any marker strategy. All species names, accession numbers, and path examples in this document use generic placeholders. Replace these with actual values for each study. When modifying this skill, keep all examples generic so the skill remains reusable across projects.

The alignment pipeline has five phases, each with a required report to the researcher:

```
Phase 1 → Pre-alignment QC     (report length distributions → human approves any filtering)
Phase 2 → MAFFT alignment      (per-marker, strategy from locus guide)
Phase 3 → Post-MAFFT QC        (remove >75% gap sequences BEFORE trimAl)
Phase 4 → trimAl + AMAS        (trim columns, then concatenate)
Phase 5 → Concatenated QC      (report per-taxon and per-marker missing data → human decides threshold)
```

**Do not skip any phase. Do not remove sequences or markers without human approval
(except Phase 3 automatic >75% gap filter, which must still be reported).**

---

## Phase 1 — Pre-Alignment Quality Check

Run before MAFFT. Every finding must be reported; any proposed filtering requires human approval.

### 1a. Check sequence headers
Verify consistent formatting across all files. Headers must be unique and follow the project
naming convention (`Genus_species_ACCESSION` or equivalent). If headers are inconsistent,
normalize before alignment and document the normalization rule used.

### 1b. Check sequence length distribution per marker

For each marker, compute:
- Count, minimum, maximum, median, mean length
- Number and percentage of sequences below 50% of the median length

```python
# Quick check per marker
from statistics import median
lengths = [len(seq) for seq in sequences]
med = median(lengths)
n_short = sum(1 for l in lengths if l < 0.5 * med)
```

Report this table for every marker:

| Marker | Count | Min | Max | Median | Mean | Seqs < 50% median | Notes |
|--------|-------|-----|-----|--------|------|-------------------|-------|

### 1c. Propose filtering to the researcher — HUMAN APPROVAL REQUIRED

If any marker has sequences substantially shorter than the average for that marker
(suggested flag: < 50% of the median length), **do not remove automatically**.
Instead:

1. Present the researcher with a table showing the flagged sequences:
   - Taxon name, accession, length, % of median
2. Explain the potential consequence: trimAl will over-trim the alignment if partial
   sequences create sparse columns
3. **Wait for the researcher to approve or reject removal before proceeding**

Example report to present:
```
Proposed pre-alignment filtering for ycf1:
  42 sequences below 1947 bp (50% of median 3894 bp):
  - Alpinia_cf._DRR502928: 510 bp (13% of median)
  - Amomum_biphyllum_SRR12824540: 510 bp (13% of median)
  ... [full list]
  Reason: these partial sequences will cause trimAl to remove most alignment columns.
  Recommend: remove sequences < 1500 bp (practical threshold).
  Await researcher decision before proceeding.
```

Only filter after explicit researcher approval. Document every removed sequence
(taxon, marker, length, reason) in the report.

---

## Phase 2 — Select MAFFT Strategy and Align

**Always read `skills/alignment/references/locus-guide.md` before choosing a strategy.**
It contains per-marker recommendations for all standard plastid markers, ITS/ETS, and
Angiosperms353. The table below is a summary; the locus guide has details on problematic
regions and codon-aware alignment.

| Dataset | Sequence type | Recommended strategy | Notes |
|---------|--------------|---------------------|-------|
| ≤200 sequences, moderate divergence | DNA / protein | `--linsi` | Most accurate; iterative local alignment |
| ≤200 sequences, high divergence | DNA / protein | `--ginsi` or `--einsi` | Global/local pairwise; handles gaps better |
| 200–10,000 sequences | DNA | `--fftnsi` or `--auto` | FFT-based; fast, good quality |
| >10,000 sequences | DNA | `--auto` or `--parttree` | Speed-accuracy trade-off |
| Highly similar (intraspecific) | DNA | `--fftns2` or `--auto` | Population-level |
| Protein sequences | Amino acid | `--auto --amino` | Always flag `--amino` explicitly |
| Unknown / mixed | Any | `--auto` | Safe default |

```bash
# ≤200 taxa, moderate divergence
mafft --localpair --maxiterate 1000 --thread N input.fasta > aligned.fasta

# >200 taxa
mafft --retree 2 --maxiterate 2 --thread N input.fasta > aligned.fasta

# Intergenic spacers (trnL, psbA-trnH) — MANDATORY einsi
mafft --ep 0 --genafpair --maxiterate 1000 --thread N input.fasta > aligned.fasta
```

Check available CPUs (`nproc`) before setting `--thread`.

Align each marker independently. Never concatenate unaligned sequences.

**Alternative — MUSCLE** (if MAFFT produces poor results on a specific marker):
```bash
muscle -align input.fasta -output aligned.fasta -threads N
```

---

## Phase 3 — Post-MAFFT QC: Remove High-Gap Sequences BEFORE trimAl

**This step must be done before trimAl.** Sequences with large proportions of gaps in the
MAFFT alignment create sparse columns that cause trimAl to over-aggressively remove
alignment positions, discarding real phylogenetic signal.

### 3a. Check per-sequence gap percentage in each aligned file

For every `*_aligned.fasta`, compute the gap percentage per sequence:
```python
gap_pct = (seq.count("-") + seq.count("?") + seq.count("N")) / len(seq) * 100
```

### 3b. Remove sequences with >75% gaps

Sequences with >75% gaps in the MAFFT-aligned file cover <25% of the alignment and
contribute minimal phylogenetic signal while biasing trimAl. Remove them automatically
**from the aligned file** before running trimAl.

```python
# Filter high-gap sequences from aligned FASTA
kept = {id_: seq for id_, seq in seqs.items()
        if (seq.count("-") + seq.count("?")) / len(seq) * 100 <= 75.0}
```

### 3c. Report all removals

Report every removed sequence: taxon, marker, gap%, non-gap bp retained.
Use `scripts/alignment/analyze_alignment.py` to generate this report:

```bash
python scripts/alignment/analyze_alignment.py \
  --input data/aligned/ \
  --max_gap_pct 75.0 \
  --output alignment_diagnostics.tsv
```

Report format (per marker):
```
matK: removed N sequences (>75% gaps)
  Taxon_name_ACC: 93% gaps (108 non-gap bp)
  ...
```

**Note:** This 75% threshold is automatic — it does not require human approval because
sequences exceeding it are nearly unusable for phylogenetic inference. However, the
removals must still be fully documented in the alignment report.

---

## Phase 4 — trimAl Column Trimming and AMAS Concatenation

Run trimAl on the gap-filtered aligned files (output of Phase 3).

**trimAl** (primary):
```bash
# Automated mode — good for phylogenomics
trimal -in filtered_aligned.fasta -out trimmed.fasta -automated1

# For intergenic spacers (trnL, psbA-trnH) — conservative column trimming
trimal -in filtered_aligned.fasta -out trimmed.fasta -gappyout
```

Do not trim if the alignment is already clean (no ragged ends, low gap columns).
Unnecessary trimming removes phylogenetic signal.

**Gblocks** (alternative — more conservative, produces discrete blocks):
```bash
Gblocks aligned.fasta -t=d -b4=5 -b5=h
```

**AMAS concatenation:**
```bash
python AMAS.py concat \
  -i marker1_trimmed.fasta marker2_trimmed.fasta ... \
  -f fasta -d dna \
  -t concatenated.fasta \
  -p partition.txt \
  -u fasta -y raxml
```

Record alignment length before and after trimming for each marker. If trimming removes
>30% of sites, investigate — this may indicate residual partial sequences or alignment
strategy mismatch.

---

## Phase 5 — Multi-Locus Decision and Concatenated Matrix QC

### 5a. Concatenation vs. coalescent decision

**Concatenation (supermatrix)** — recommended when:
- Markers are largely congruent (plastid markers are genetically linked)
- Dataset has incomplete taxon coverage per gene (missing data handled better)
- Faster and simpler; IQ-TREE handles partitioned models natively

**Coalescent (gene-tree summary — ASTRAL)** — recommended when:
- Rapid radiation expected (high ILS)
- Gene trees show substantial discordance
- Note: plastid markers are physically linked and do not represent independent loci —
  coalescent is primarily appropriate for nuclear markers

When uncertain: default to concatenation; revisit after seeing gene trees.

### 5b. Per-marker statistics report (required)

After trimming, report for every marker:

| Marker | Sequences | Aligned length | Gap % | Trimmed length | PIS | PIS% | Notes |
|--------|-----------|---------------|-------|---------------|-----|------|-------|

QC thresholds:
- Aligned length: within expected range for the marker
- Gap % after trimming: < 30%
- PIS%: > 5% (flag if lower — may indicate insufficient variation)

### 5c. Concatenated matrix missing data — REPORT TO HUMAN, AWAIT DECISION

After AMAS, compute and report two tables:

**Table 1 — Per-marker coverage in the concatenated matrix:**
For each marker partition, report:
- Number of taxa with data for that marker
- % of all taxa covered
- % of cells that are missing (gap or N) across the marker columns

| Marker | Taxa with data | Coverage % | Missing data % |
|--------|---------------|------------|----------------|

**Table 2 — Per-taxon missing data in the concatenated matrix:**
For each taxon (row), report:
- Total non-gap positions across all marker partitions
- % missing across the full supermatrix

Report the distribution (e.g., how many taxa have >80%, >60%, >40% missing), and list
the worst 20 taxa.

**Then pause and present both tables to the researcher.**

The researcher will specify:
- A **sequence threshold**: remove taxa with > X% missing data (e.g., 80%)
- A **marker threshold**: remove markers with < Y% taxon coverage (e.g., 10%)

The AI agent must not decide these thresholds independently. Present the data and wait.

### 5d. Apply researcher's thresholds and re-concatenate

After the researcher specifies thresholds:
1. Remove sequences exceeding the sequence threshold from all trimmed marker files
2. Remove marker files below the marker coverage threshold
3. Re-run AMAS concatenation
4. Re-report the same tables so the researcher can confirm the result

Document in the report:
- The thresholds chosen (and who chose them)
- Every taxon removed and from which markers
- Every marker removed

---

## Phase 5.5 — Post-Tree SNP Density Scan (Long Branch QC)

Run after tree inference if any taxon shows an anomalously long branch (long branch attraction),
or as a routine QC step before submitting to tree inference on datasets that include GenBank
sequences (which may be complete plastomes or other non-target-length records).

### When to run

- After tree inference reveals a taxon with a branch length ≥5× the median for its clade
- Whenever GenBank downloads are used without strict length filtering (whole plastomes are
  frequently downloaded when searching for individual genes because the gene is annotated on them)
- As a standard QC pass on any concatenated supermatrix before publication
- Before running this analysis, the AI agent should check with the human for approval to proceed. The AI agent should present the proposed plan for k-mer window size (default: `100` bp). After running this analysis, the AI agent should present the results to the human and ask for approval before removing any flagged windows. 
When removing flagged windows, the AI Agent should only remove flagged windows instead of the whole sequence. Then, missing percentage of the alignment should be re-calculated and reported to the human. The AI agent should also report the number of flagged windows removed and the number of flagged windows remaining for each sequence. Everything should be clearly documented in the report. 

### What it detects

The script divides the alignment into windows (default: `alignment_length / 100` bp) and
computes per-sequence SNP density in each window relative to the column consensus. It then
flags sequences whose SNP density in any window exceeds a Z-score threshold (default: 3.0).

| Pattern | Likely cause |
|---------|-------------|
| Very high Z-score (>8) in many windows | Whole-plastome sequence misaligned to individual gene |
| Multiple flagged windows in one gene partition | Chimeric sequence or IR-boundary artifact |
| Single flagged window at alignment edge | Alignment edge effect; usually not a problem |
| All outgroup taxa flagged in one partition | Outgroup divergence (not a data error — review manually) |

### Command

```bash
python scripts/alignment/snp_density_scan.py \
  -i data/aligned/concatenated_final.fasta \
  -o data/aligned/qc/snp_density_scan.tsv \
  -f data/aligned/qc/snp_density_flagged.tsv \
  -v
```

Options:
- `-w` — override window size (bp); default = alignment_length / 100
- `-z` — Z-score threshold (default 3.0; raise to 4–5 to reduce false positives from divergent outgroups)
- `-m` — minimum non-gap coverage per window to include in stats (default 0.5)

### Interpreting the output

The flagged summary file (`snp_density_flagged.tsv`) lists:

| Column | Meaning |
|--------|---------|
| `n_flagged_windows` | Number of windows exceeding Z-score threshold |
| `max_zscore` | Highest Z-score across all windows |
| `max_zscore_window` | Alignment positions of the worst window |
| `flag_rate_pct` | % of all tested windows that were flagged |
| `recommendation` | INVESTIGATE (multiple flags) / REVIEW (1–2 flags) |

### Action thresholds — REPORT TO HUMAN

Present the flagged summary to the researcher. Do not remove sequences automatically.

| Criterion | Action |
|-----------|--------|
| `max_zscore > 8` AND `n_flagged_windows > 3` | Strong evidence of data error — recommend removal or re-extraction |
| `n_flagged_windows ≥ 3` AND `flag_rate ≥ 10%` | Investigate source record (check if whole plastome was used) |
| Outgroup taxa flagged in a consistent partition | Likely genuine divergence — do not remove without researcher approval |
| Single flagged window (`n_flagged_windows = 1`) | Usually benign — report but do not recommend removal |

The researcher decides whether to remove or retain flagged sequences after reviewing the findings.

---

## Step 6 — Organize Outputs

```
data/aligned/
  <marker>_aligned.fasta       — MAFFT output (before gap filtering)
  <marker>_filtered.fasta      — after Phase 3 high-gap removal (if any removed)
  <marker>_trimmed.fasta       — after trimAl
  concatenated.fasta           — supermatrix
  partition.txt                — RAxML-style partition file
  alignment_stats.tsv          — per-marker statistics
  removed_high_gap.tsv         — Phase 3 removals log
  alignment_diagnostics.tsv    — full per-sequence diagnostics
```

Partition file format (RAxML-style, also accepted by IQ-TREE):
```
DNA, matK = 1-873
DNA, rbcL = 874-1423
DNA, ITS  = 1424-2035
```

---

## Report Template

**Mandatory:** Every log file generated during this module must be listed with its exact path in the report so the researcher can monitor background processes and audit what ran.

Write to `reports/[planX/]alignment_YYYY-MM-DD.md`:

```markdown
# Alignment Report
Date: YYYY-MM-DD
Plan: [planA / planB / ...]

## Phase 1 — Pre-Alignment QC

### Sequence Length Distributions
| Marker | Count | Min | Max | Median | Seqs < 50% median |
|--------|-------|-----|-----|--------|-------------------|

### Pre-Alignment Filtering (human-approved)
[List every removed sequence: taxon, marker, length, reason]
[State the threshold used and confirm human approved it]

## Phase 2 — Alignment Strategy
- Tool: MAFFT [version]
- Strategy per marker: [marker: strategy, ...]
- Justification: [why each strategy was chosen]

## Phase 3 — Post-MAFFT High-Gap Removal
[List every removed sequence: taxon, marker, gap%, non-gap bp]
[Threshold used: >75% gaps — automatic, no human approval required]

## Phase 4 — Trimming
- trimAl mode per marker: [marker: mode, ...]

## Phase 5 — Per-Marker Statistics (post-trimming)
| Marker | Sequences | Aligned length | Gap % | Trimmed length | PIS | PIS% | QC |
|--------|-----------|---------------|-------|---------------|-----|------|----|

## Concatenated Matrix — Missing Data Report

### Table 1: Per-Marker Coverage
| Marker | Taxa with data | Coverage % | Missing % |
|--------|---------------|------------|-----------|

### Table 2: Per-Taxon Missing Data (worst 20)
| Taxon | Non-gap bp | Missing % |
|-------|-----------|-----------|

[Distribution summary: N taxa > 80%, N taxa 60-80%, N taxa < 60% missing]

### Researcher Decision
- Sequence threshold: [X%] — applied / not applied
- Marker threshold: [Y%] — applied / not applied
- [List any taxa or markers removed at this stage]

## Final Matrix Summary
- Total length: X bp
- Total taxa: N
- Overall missing data: X%
- Partition file: [path]

## All Sequences Removed (complete log)
[Phase 1, Phase 3, Phase 5 removals consolidated]

## Software Versions
| Tool | Version | Source |
|------|---------|--------|

## Log Files Generated
[List every log file created during this module with its exact path so the researcher
 can monitor background alignments and review per-marker details]
[Examples:]
[  data/aligned/align_markers_2026-04-19.log        (align_markers.sh stdout)]
[  data/aligned/qc/snp_density_scan.tsv             (full per-window SNP density table)]
[  data/aligned/qc/snp_density_flagged.tsv          (flagged sequences summary)]
[  data/aligned/alignment_diagnostics.tsv           (analyze_alignment.py output)]

## Next Module
model-selection
```

---

## Scripts

Pre-built scripts for this module are in `scripts/alignment/`:

| Script | Purpose |
|--------|---------|
| `align_markers.sh` | MAFFT alignment per marker with strategy selection; optional trimAl; optional AMAS concatenation. Use `-A <path/to/AMAS.py>` when AMAS is not in PATH |
| `analyze_alignment.py` | Per-marker and per-sequence diagnostics: gap%, parsimony-informative sites, outlier detection |
| `snp_density_scan.py` | Sliding-window SNP density scan on aligned FASTA; flags sequences with anomalously high SNP density per window — detects whole-plastome misalignment, chimeric sequences, and long-branch-attraction candidates |

Usage examples:
```bash
# Align one marker with einsi, trim, no concatenation
bash scripts/alignment/align_markers.sh \
  -i data/for_alignment_trnL/ -o data/aligned/ -s einsi -T -m gappyout -t 8

# Analyze alignment quality and flag outliers
python scripts/alignment/analyze_alignment.py \
  --input data/aligned/ \
  --max_gap_pct 75.0 \
  --output data/aligned/alignment_diagnostics.tsv

# Sliding-window SNP density scan (run on concatenated supermatrix or any aligned FASTA)
python scripts/alignment/snp_density_scan.py \
  -i data/aligned/concatenated_final.fasta \
  -o data/aligned/qc/snp_density_scan.tsv \
  -f data/aligned/qc/snp_density_flagged.tsv \
  -v
```

---

## Common Mistakes

| Mistake | Fix |
|---------|-----|
| Running trimAl before removing high-gap sequences | Always remove >75% gap sequences from the MAFFT output BEFORE trimAl — large gaps bias trimAl and cause it to discard real signal |
| Removing sequences without human approval | Pre-alignment length filtering (Phase 1) and post-AMAS taxon/marker filtering (Phase 5) both require human approval. Only the >75% gap filter (Phase 3) is automatic |
| Using `--auto` for all datasets | High-divergence datasets need `--ginsi`/`--einsi`; `--auto` may misalign intergenic spacers |
| Aligning all markers together in one file | Align each marker separately, then concatenate post-alignment |
| Trimming clean alignments | Only trim when ragged ends or high gap columns are present |
| Silently removing sequences | Every removed sequence must appear in the report with taxon, marker, and reason |
| Forgetting to report missing data to the researcher | Steps 5c–5d are mandatory — the researcher decides thresholds, not the AI agent |
| Forgetting to record alignment length before and after trimming | Both lengths go in the report; trimming that removes >30% of sites needs justification |
