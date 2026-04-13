# Support Value Reference — Tree Inference

Quick reference for interpreting branch support values across methods.
Read this before reporting support thresholds in the Methods section.

---

## Support Value Thresholds by Method

| Method | Flag | Well-supported threshold | Note |
|--------|------|--------------------------|------|
| UFBoot (IQ-TREE ultrafast bootstrap) | `-B` | **≥ 95** | Inflated relative to standard bootstrap; ≥95 ≈ ≥70 standard |
| Standard bootstrap | `-b` | **≥ 70** | Felsenstein 1985 convention; ≥70 = reasonably well-supported |
| SH-aLRT (likelihood ratio test) | `--alrt` | **≥ 80** | Fast; run alongside UFBoot for independent confirmation |
| Bayesian posterior probability (MrBayes, BEAST2) | — | **≥ 0.95** | Highest reliability when MCMC has converged |
| ASTER local posterior | — | **≥ 0.95** | Branch-specific support from gene tree concordance |
| Transfer bootstrap expectation (TBE) | — | **≥ 0.70** | Used for very large datasets; more stable than FBP |

**Combined threshold (IQ-TREE recommended):**  
A node is considered well-supported when **UFBoot ≥ 95 AND SH-aLRT ≥ 80**.  
Either alone is less reliable than both together.

---

## Common Mistake: UFBoot ≠ Standard Bootstrap

UFBoot values are systematically inflated. Quoting a UFBoot value of 75 as "moderate support" is incorrect.

| UFBoot value | Approximate standard bootstrap equivalent |
|-------------|------------------------------------------|
| 95 | ~70 (threshold) |
| 85 | ~50–60 (weak) |
| 75 | ~40–50 (unreliable) |

**Methods section language:**  
> "Branch support was assessed using 1,000 ultrafast bootstrap replicates (Hoang et al. 2018) and the SH-aLRT test (Guindon et al. 2010). Nodes with UFBoot ≥ 95 and SH-aLRT ≥ 80 were considered well-supported."

---

## ASTER / Coalescent Support

ASTER (wASTRAL / ASTRAL-Pro3) places **local posterior probabilities** on branches.

| Value | Interpretation |
|-------|---------------|
| ≥ 0.95 | Well-supported |
| 0.70–0.95 | Moderate support; note in text |
| < 0.70 | Poorly supported; likely ILS or insufficient gene trees |

Local posterior values are distinct from bootstrap values — do not mix them in figures.

---

## Concordance Factors (IQ-TREE `--gcf / --scf`)

Run after species tree is available to quantify actual gene-tree support:

| Factor | What it measures | Interpretation |
|--------|-----------------|----------------|
| gCF (gene concordance factor) | % of informative gene trees supporting the branch | gCF ≥ 50% → majority of genes agree |
| sCF (site concordance factor) | % of informative alignment sites supporting the branch | Useful when gene tree numbers are low |

**gCF vs. UFBoot:**  
High UFBoot + low gCF = concatenation is overconfident; likely ILS present.  
Always compute gCF for multi-locus datasets spanning multiple families or orders.

```bash
# After ASTER species tree
iqtree2 -t output/aster/species_tree_wastral.nwk \
  --gcf output/gene_trees/all_gene_trees.txt \
  -p data/aligned/concatenated.fasta \
  --scf 100 --prefix concordance -T AUTO
```

---

## Bayesian Convergence Diagnostics

### MrBayes

| Diagnostic | Threshold | Where to check |
|-----------|-----------|---------------|
| ASDSF (average std deviation of split frequencies) | **< 0.01** | MrBayes console output |
| PSRF (potential scale reduction factor) | **≈ 1.0** (< 1.01 ideal) | MrBayes `.p` file |
| ESS (effective sample size) | **> 200** for all parameters | Tracer (`.p` file) |
| Burnin | 25% default; extend if trace plots show long burnin | Visual inspection in Tracer |

### BEAST2

| Diagnostic | Threshold | Where to check |
|-----------|-----------|---------------|
| ESS (all parameters) | **> 200** | Tracer (`.log` file) |
| TreeHeight ESS | **> 200** | Tracer |
| clockRate ESS | **> 200** | Tracer |
| Trace plots | Stationary (no trend after burnin) | Tracer visual |

**If ESS < 200:** multiply chain length by ×2–×5. Do not report results from unconverged runs.

---

## RAxML-NG Support

RAxML-NG uses standard bootstrap by default (same as `-b` in IQ-TREE):
- Threshold: **≥ 70** for standard bootstrap

For fast bootstrap in RAxML-NG use `--bs-trees autoMRE` (converges when bootstrap support stabilizes).

---

## Reporting Support in Figures

**Best practice:**
1. Show UFBoot (or BS) as numbers on branches for key nodes only
2. For figures with many nodes: use symbols (●= ≥95, ○= ≥70)
3. Never show every bootstrap value on every branch — clutters the figure
4. For combined UFBoot/SH-aLRT: `95/82` format on branches
5. For combined concatenation + ASTER: use separate trees or flag discordant nodes

---

## References

- Felsenstein J. 1985. Confidence limits on phylogenies: an approach using the bootstrap. *Evolution* 39:783–791.
- Guindon S et al. 2010. New algorithms and methods for estimating ML phylogenies. *Syst Biol* 59:307–321.
- Hoang DT et al. 2018. UFBoot2: improving the ultrafast bootstrap approximation. *Mol Biol Evol* 35:518–522.
- Zhang C et al. 2023. ASTRAL-Pro 2: ultrafast species tree reconstruction from multi-copy gene family trees. *Bioinformatics* 39.
- Minh BQ et al. 2020. IQ-TREE 2: new models and methods for phylogenetic inference. *Mol Biol Evol* 37:1530–1534.
