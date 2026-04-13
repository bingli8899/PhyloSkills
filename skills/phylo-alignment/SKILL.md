---
name: phylo-alignment
description: Use after assembled or downloaded sequences are organized and ready for alignment. Selects alignment strategy based on dataset size and sequence type, runs per-marker alignments, optionally trims, and guides the concatenation vs. coalescent decision for multi-locus datasets. Use when raw FASTA files need to be aligned before model selection and tree inference.
---

# Aligning Phylogenetic Sequences

## Overview

Align each marker separately, assess quality, trim if needed, then build the combined matrix. Tool and parameter choice depends on dataset size and sequence type — one setting does not fit all.

## Step 1 — Read upstream reports

Load `reports/[planX/]assembly_YYYY-MM-DD.md` or `reports/data-acquisition_YYYY-MM-DD.md`. Extract:
- Number of sequences per marker
- Sequence type (DNA / protein / rRNA)
- Expected sequence length range per marker
- Number of markers (single vs. multi-locus)

## Step 2 — Select MAFFT strategy

MAFFT is the primary tool. Strategy selection is non-obvious — use this table:

| Dataset | Sequence type | Recommended strategy | Notes |
|---------|--------------|---------------------|-------|
| ≤200 sequences, moderate divergence | DNA / protein | `--linsi` | Most accurate; iterative local alignment |
| ≤200 sequences, high divergence | DNA / protein | `--ginsi` or `--einsi` | Global/local pairwise; slower but handles gaps better |
| 200–10,000 sequences | DNA | `--fftnsi` or `--auto` | FFT-based; fast, good quality |
| >10,000 sequences | DNA | `--auto` or `--parttree` | Speed-accuracy trade-off; `--parttree` for very large |
| Highly similar (intraspecific / population) | DNA | `--fftns2` or `--auto` | Over-alignment not a concern at this scale |
| Protein sequences | Amino acid | `--auto --amino` | Always flag `--amino` explicitly |
| Unknown / mixed | Any | `--auto` | Safe default; MAFFT chooses internally |

```bash
# Typical phylogenetic dataset (≤200 taxa, gene-level)
mafft --linsi --thread 4 input.fasta > aligned.fasta

# Large dataset
mafft --auto --thread 8 input.fasta > aligned.fasta

# Protein
mafft --auto --amino --thread 4 input.fasta > aligned.fasta
```
Note that the server could have different number of processors to be used. Need to check the number of processors first and then adjust `--thread` accordingly.

**Alternative — MUSCLE** (use if MAFFT produces poor results on a specific marker):
```bash
muscle -align input.fasta -output aligned.fasta -threads 4
```

Align each marker independently. Never concatenate unaligned sequences.

## Step 3 — Inspect alignment quality

Open in AliView or inspect programmatically before any trimming:
- Scan for obviously misaligned sequences (single sequences out of phase with the rest)
- Check for sequences that are mostly gaps (likely wrong gene or chimeric)
- Confirm expected alignment length (e.g., matK ~900 bp, rbcL ~550 bp, ITS ~600 bp)

Flag problematic sequences — remove and note in report rather than silently trimming around them.

## Step 4 — Trim alignment (optional but recommended)

Trimming removes poorly aligned and gap-rich columns. Use when:
- Alignment has ragged ends
- Gap % across columns is high (>50% in many columns)
- Dataset mixes sequences of different lengths (e.g., partial records)

**trimAl** (primary):
```bash
# Automated mode — good for phylogenomics
trimal -in aligned.fasta -out trimmed.fasta -automated1

# Conservative gap-based trimming
trimal -in aligned.fasta -out trimmed.fasta -gappyout

# Check what was removed
trimal -in aligned.fasta -out trimmed.fasta -automated1 -htmlout trim_report.html
```

**Gblocks** (alternative — more conservative, produces blocks):
```bash
Gblocks aligned.fasta -t=d -b4=5 -b5=h
```

Do not trim if the alignment is already clean — unnecessary trimming removes phylogenetic signal.

## Step 5 — Multi-locus decision: concatenation vs. coalescent

If the dataset has multiple markers, ask the researcher:

**Concatenation (supermatrix)** — recommended when:
- Markers are largely congruent (no strong gene-tree discordance expected)
- Dataset has incomplete taxon coverage per gene (missing data handled better)
- Faster and simpler; IQ-TREE handles partitioned models natively

```bash
# Concatenate with AMAS
python AMAS.py concat -i marker1.fasta marker2.fasta marker3.fasta \
  -f fasta -d dna -o concatenated.fasta --part-format raxml
```

**Coalescent (gene-tree summary — ASTRAL)** — recommended when:
- Rapid radiation expected (high incomplete lineage sorting / ILS)
- Gene trees show substantial discordance
- All markers have good taxon coverage

Coalescent requires individual gene trees first (run `phylo-tree-inference` per marker, then summarize with ASTRAL). Note this in the report and coordinate with the researcher.

**When uncertain:** default to concatenation, note the caveat, revisit after seeing the gene trees.

## Step 6 — QC gate

| Check | Threshold | Action on failure |
|-------|-----------|-------------------|
| Alignment length | Within expected range for marker | Inspect for misaligned seqs; route to `phylo-debug` |
| Mean gap % per column | <30% after trimming | Re-trim or remove partial sequences |
| Sequences removed during QC | Document every removal | Note taxon, reason, and marker in report |
| Concatenated matrix completeness | Note overall missing data % | Inform researcher; >50% missing data warrants discussion |
| Parsimony-informative sites | >10% of alignment length | Flag if very low — may indicate insufficient variation |

On any failure → route to `phylo-debug`.

## Step 7 — Organize outputs

```
data/aligned/
  matK_aligned.fasta
  rbcL_aligned.fasta
  ITS_aligned.fasta
  concatenated.fasta          # supermatrix (if concatenation chosen)
  partition.txt               # partition file for IQ-TREE / RAxML
```

Partition file format (RAxML-style, also accepted by IQ-TREE):
```
DNA, matK = 1-873
DNA, rbcL = 874-1423
DNA, ITS  = 1424-2035
```

## Report

Write to `reports/[planX/]alignment_YYYY-MM-DD.md`:

```markdown
# Alignment Report
Date: YYYY-MM-DD
Plan: [planA / planB / ...]

## Alignment Strategy
- Tool: MAFFT [version] / MUSCLE [version]
- Strategy per marker: [e.g., matK --linsi, rbcL --auto]
- Trimming: [trimAl automated1 / Gblocks / none] with justification

## Per-Marker Statistics
| Marker | Sequences | Aligned length | Gap % | Trimmed length | Removed seqs |
|--------|-----------|---------------|-------|---------------|-------------|

## Multi-Locus Approach
- [Concatenation / Coalescent] — justification

## Concatenated Matrix (if applicable)
- Total length: X bp
- Total taxa: N
- Overall missing data: X%
- Partition file: [path]

## Sequences Removed
[Taxon, marker, reason for each removal]

## Software Versions
| Tool | Version | Source | Install date |
|------|---------|--------|-------------|

## Next Module
phylo-model-selection
```

## Scripts

Pre-built scripts for this module are in `skills/phylo-alignment/scripts/`. Load when needed:

| Script | Purpose |
|--------|---------|
| `align_markers.sh` | MAFFT alignment per marker with strategy selection; optional trimAl; optional AMAS concatenation |
| `analyze_alignment.py` | Per-marker and per-sequence diagnostics: gap%, parsimony-informative sites, outlier detection |

Usage examples:
```bash
# Align all markers in a directory, trim with trimAl, then concatenate
bash skills/phylo-alignment/scripts/align_markers.sh \
  -i data/cds/ -o data/aligned/ -s linsi -T -c

# Analyze alignment quality and flag problems
python skills/phylo-alignment/scripts/analyze_alignment.py \
  --input data/aligned/ --output alignment_diagnostics.tsv
```

## Common Mistakes

| Mistake | Fix |
|---------|-----|
| Using `--auto` for all datasets without considering divergence | High-divergence datasets need `--ginsi`/`--einsi`; `--auto` may misalign them |
| Aligning all markers together in one file | Align each marker separately, then concatenate post-alignment |
| Trimming clean alignments | Only trim when ragged ends or high gap columns are present |
| Silently removing sequences without documenting | Every removed sequence must appear in the report with a reason |
| Assuming concatenation is always appropriate | Check for rapid radiations or known ILS — coalescent may be warranted |
| Forgetting to record alignment length before and after trimming | Both lengths go in the report; trimming that removes >30% of sites needs justification |
