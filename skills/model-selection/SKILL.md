---
name: model-selection
description: Use after alignment is complete and before tree inference. Selects the best-fit substitution model(s) for each marker or partition using ModelFinder or ModelTest-NG. Use when the researcher needs to determine which substitution model to apply in IQ-TREE, RAxML, MrBayes, or BEAST2, including partitioned schemes for multi-locus datasets.
---

# Selecting Substitution Models for Phylogenetic Analysis

## Overview

Run model testing on the aligned matrix and extract the best-fit model string(s) ready for tree inference. IQ-TREE's ModelFinder is the primary approach — it tests models and selects the tree in one run if needed.

## Step 1 — Read the alignment report

Load `reports/[planX/]alignment_YYYY-MM-DD.md`. Extract:
- Aligned FASTA path(s) and partition file path (if multi-locus)
- Sequence type (DNA / protein)
- Whether the analysis is unpartitioned or partitioned
- Downstream tool (IQ-TREE, RAxML-NG, MrBayes, BEAST2) — determines output format needed

## Step 2 — Choose criterion: BIC, AIC, or AICc

| Criterion | When to use |
|-----------|-------------|
| **BIC** (default) | Most phylogenetic analyses; penalizes complexity more than AIC; preferred for large alignments |
| **AICc** | Short alignments (<500 bp) or small taxon sets (<50 taxa); corrects AIC for small samples |
| **AIC** | Acceptable; less preferred in phylogenetics — tends to select over-parameterized models |

ModelFinder uses BIC by default. Override with `--crit AIC` or `--crit AICc` if needed.

## Step 3 — Run model selection

### Unpartitioned (single marker)

```bash
# Model test only — no tree search
iqtree2 -s aligned.fasta -m TEST --crit BIC -T AUTO --prefix model_test

# Output: model_test.iqtree — look for "Best-fit model:" line
grep "Best-fit model" model_test.iqtree
```

### Partitioned (multi-locus)

```bash
# Test models per partition
iqtree2 -s concatenated.fasta -p partition.txt -m TEST --crit BIC \
  -T AUTO --prefix partitioned_model_test

# With partition merging — finds best grouping of similar partitions
iqtree2 -s concatenated.fasta -p partition.txt -m TEST --crit BIC \
  --merge rclusterf -T AUTO --prefix merged_model_test
```

`--merge rclusterf` is recommended for large partitioned datasets — it can reduce model complexity by merging partitions that share the same best-fit model, improving computational efficiency without sacrificing fit.

### Alternative — ModelTest-NG (standalone)

```bash
modeltest-ng -i aligned.fasta -d nt -o modeltest_output \
  --force -p 4 --crit BIC
# Results in modeltest_output.log and modeltest_output.out
```

Use when IQ-TREE is not available or for independent verification.

## Step 4 — Interpret results

**Key output from IQ-TREE (`*.iqtree` file):**
```
Best-fit model: GTR+F+I+G4 chosen according to BIC

Model of substitution: GTR+F+I+G4
  Rate matrix:        [6 substitution rates]
  State frequencies:  Empirical (from data)
  Proportion of invariable sites: 0.12
  Gamma shape parameter: 0.87
```

**Common model suffixes:**

| Suffix | Meaning |
|--------|---------|
| `+G4` | Gamma rate variation, 4 categories |
| `+I` | Proportion of invariable sites |
| `+I+G4` | Both invariable sites and gamma |
| `+F` | Empirical base frequencies |
| `+FC` | Counted base frequencies |
| `+R4` | FreeRate model, 4 categories (more flexible than +G) |

**Common DNA base models:**

| Model | Complexity | Typical use |
|-------|-----------|-------------|
| JC | Simplest | Rarely selected in real data |
| K2P / K80 | 2 rate classes | Very similar sequences |
| HKY | Transitions/transversions + freq | Moderate divergence |
| TrN, TVM, TIM | Intermediate | Varies |
| GTR | Most parameter-rich | High divergence; most common in phylogenomics |

**Protein models:** LG, WAG, JTT are most common. ModelFinder tests these automatically with `--amino`.

**Explain to a novice researcher (calibrate to level):**
- More complex models (GTR) fit better but use more parameters
- BIC penalizes unnecessary parameters — the selected model is the best balance
- GTR+G4 or GTR+I+G4 being selected is normal and expected for most datasets

## Step 5 — Extract model string for tree inference

From IQ-TREE output, copy the exact model string. For partitioned analyses, the updated partition file (`.best_model.nex` or `-p` output) contains per-partition models — use this file directly in tree inference.

```bash
# The merged/updated partition file is the key output
ls *.best_model.nex     # IQ-TREE writes this automatically
```

For MrBayes or BEAST2: translate the IQ-TREE model to the tool's syntax:
- `GTR+G4` → MrBayes: `lset nst=6 rates=gamma ngammacat=4`
- `HKY+G4` → MrBayes: `lset nst=2 rates=gamma ngammacat=4`
- BEAST2: set model in BEAUti XML; note the IQ-TREE model in the report for reference

## QC Gate

- [ ] Model selected for every partition / marker
- [ ] Criterion (BIC/AICc/AIC) documented and justified
- [ ] Model string verified as valid input for the downstream tree inference tool
- [ ] Partition merging decision documented (merged or not, with reason)

## Report

Write to `reports/[planX/]model-selection_YYYY-MM-DD.md`:

```markdown
# Model Selection Report
Date: YYYY-MM-DD
Plan: [planA / planB / ...]

## Configuration
- Tool: IQ-TREE ModelFinder [version] / ModelTest-NG [version]
- Criterion: BIC / AICc / AIC — justification
- Partitioned: yes / no

## Selected Models
| Partition / Marker | Best-fit model | BIC score | Notes |
|-------------------|---------------|-----------|-------|
| matK | GTR+F+I+G4 | 12345.6 | |
| rbcL | HKY+F+G4 | 8901.2 | |
| Merged (matK+rbcL) | GTR+F+I+G4 | 20100.3 | Merged: similar models |

## Partition Merging
[Merged or not, which partitions were merged, impact on model count]

## Model Strings for Tree Inference
[Exact strings / updated partition file path ready for next module]

## Notes for Researcher
[Plain-language summary of what the model means, calibrated to level]

## Software Versions
| Tool | Version | Source | Install date |
|------|---------|--------|-------------|

## Next Module
tree-inference
```

## Common Mistakes

| Mistake | Fix |
|---------|-----|
| Using AIC instead of BIC by default | BIC is standard in phylogenetics; use `--crit BIC` or verify the default |
| Copying a model string that doesn't match the inference tool's syntax | Translate IQ-TREE model names to MrBayes/BEAST2 syntax before handing off |
| Ignoring partition merging for large datasets | `--merge rclusterf` can halve computation time with no accuracy cost |
| Assuming GTR+G is always correct without testing | Always run ModelFinder — short markers or similar sequences often select simpler models |
| Forgetting to save the `.best_model.nex` file | This file is the direct input to IQ-TREE tree inference; losing it means re-running model selection |
