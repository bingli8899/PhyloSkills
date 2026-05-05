---
name: visualization
description: Use after tree inference is complete. Produces annotated, publication-quality phylogenetic tree figures using R (ggtree, ape, phytools, ggplot2). Use when the researcher needs to plot a phylogenetic tree with support values, scale bars, clade labels, tip annotations, or divergence time axes. All visualization is R-based.
---

# Visualizing Phylogenetic Trees in R

## Overview

Read the inferred tree, root it correctly, annotate support values and clades, and export publication-quality figures. `ggtree` is the primary tool — it extends ggplot2 grammar to trees and handles most annotation needs cleanly.

## Package roles

| Package | Role |
|---------|------|
| `ape` | Tree reading, rooting, manipulation, basic plotting |
| `treeio` | Reading annotated trees (IQ-TREE, BEAST2, RAxML) |
| `ggtree` | Publication figures, layered annotation (primary) |
| `ggplot2` | Themes, multi-panel layouts, color scales |
| `phytools` | Trait mapping, phenograms, specialized visualizations |

Install if needed:
```r
install.packages(c("ape", "ggplot2"))
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("ggtree", "treeio"))
install.packages("phytools")
```

## Step 1 — Read the tree inference report

Load `reports/[planX/]tree-inference_YYYY-MM-DD.md`. Extract:
- Tree file path and format (`.treefile` from IQ-TREE, `.con.tre` from MrBayes, `.mcc.tree` from BEAST2)
- Support value type (UFBoot, standard bootstrap, posterior probability)
- Outgroup taxon(a) for rooting
- Whether this is a time-calibrated tree (BEAST2)

## Step 2 — Read and root the tree

```r
library(ape)
library(treeio)
library(ggtree)
library(ggplot2)

# IQ-TREE — preserves bootstrap on nodes
tree <- read.iqtree("output/tree_inference.treefile")

# MrBayes consensus tree
tree <- read.mrbayes("output/tree.con.tre")

# BEAST2 MCC tree (includes node ages and HPD intervals)
tree <- read.beast("output/beast_mcc.tree")

# Plain Newick (any tool)
tree <- read.tree("output/tree.nwk")

# Root on outgroup — always do this before plotting
tree_rooted <- root(as.phylo(tree), outgroup = "Outgroup_species",
                    resolve.root = TRUE)
```

Verify rooting visually before annotating — an incorrectly rooted tree invalidates the figure.

## Step 3 — Basic annotated tree (ggtree)

```r
# Core plot
p <- ggtree(tree_rooted) +
  geom_tiplab(size = 3, offset = 0.001) +
  theme_tree2()   # adds time/substitution axis

p
```

## Step 4 — Add support values

Support values sit on internal nodes. Show only values above a threshold to avoid clutter.

```r
# UFBoot or standard bootstrap from IQ-TREE
# Node labels are in $node.label after read.iqtree
p <- p + geom_nodelab(
  aes(label = ifelse(as.numeric(label) >= 70, label, "")),
  size = 2.5, hjust = 1.2, vjust = -0.3
)

# Posterior probabilities from MrBayes (often stored as prob)
p <- p + geom_nodelab(
  aes(label = ifelse(as.numeric(prob) >= 0.95,
                     round(as.numeric(prob), 2), "")),
  size = 2.5
)
```

Show both UFBoot and SH-aLRT if both were computed:
```r
# Format as "UFBoot/SH-aLRT" on node
p <- p + geom_nodelab(
  aes(label = paste0(UFboot, "/", SHaLRT)),
  size = 2, hjust = 1.2
)
```

## Step 5 — Scale bar and rooting indicator

```r
# Scale bar (substitutions per site for ML; time for BEAST2)
p <- p + geom_treescale(x = 0, y = -1, fontsize = 3)

# Or use theme_tree2() axis (shows branch lengths as axis)
p + theme_tree2() + xlab("Substitutions per site")
```

## Step 6 — Highlight clades

```r
# First identify node numbers
ggtree(tree_rooted) + geom_text(aes(label = node), size = 2)

# Then highlight by node number
p <- p +
  geom_hilight(node = 45, fill = "steelblue", alpha = 0.3) +
  geom_hilight(node = 62, fill = "tomato", alpha = 0.3)

# Add clade labels
p <- p +
  geom_cladelabel(node = 45, label = "Clade A", color = "steelblue",
                  fontsize = 3.5, offset = 0.01) +
  geom_cladelabel(node = 62, label = "Clade B", color = "tomato",
                  fontsize = 3.5, offset = 0.01)
```

## Step 7 — BEAST2 time-calibrated tree

```r
library(treeio)
beast_tree <- read.beast("beast_mcc.tree")

# Time axis (x-axis = time, root at oldest)
p_time <- ggtree(beast_tree, mrsd = "2026-01-01") +  # set most-recent sampling date if needed
  theme_tree2() +
  geom_tiplab(size = 2.5) +
  geom_range("height_0.95_HPD", color = "steelblue",
             size = 1.5, alpha = 0.5) +    # 95% HPD bars on nodes
  geom_nodelab(aes(label = round(height, 1)), size = 2, hjust = 1.3)

# Add geological time scale (optional — requires deeptime package)
# install.packages("deeptime")
# library(deeptime)
# p_time + coord_geo(...)
```

## Step 8 — Multi-panel and final styling

```r
library(ggplot2)

# Publication theme
p_final <- p +
  theme_tree2() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 8),
    plot.margin = margin(10, 20, 10, 10)
  )

# Multi-panel with patchwork or cowplot
# install.packages("patchwork")
library(patchwork)
p_final_left + p_final_right + plot_layout(ncol = 2)
```

## Step 9 — Export

```r
# PDF — preferred for publication (vector, scalable)
ggsave("figures/tree_figure.pdf", plot = p_final,
       width = 180, height = 240, units = "mm", dpi = 300)

# SVG — vector, editable in Illustrator/Inkscape
ggsave("figures/tree_figure.svg", plot = p_final,
       width = 180, height = 240, units = "mm")

# PNG — raster, for presentations
ggsave("figures/tree_figure.png", plot = p_final,
       width = 180, height = 240, units = "mm", dpi = 600)
```

Standard journal figure widths: single column ≈ 85 mm, double column ≈ 170–180 mm.

## QC Gate

| Check | Action on failure |
|-------|-------------------|
| Outgroup roots tree correctly | Re-root with `root()`; inspect for long-branch issues |
| All tip labels present and correctly spelled | Cross-check against original FASTA headers |
| Support values match tree inference report summary | Verify node label parsing — check `tree$node.label` directly |
| No overlapping tip labels | Reduce `size`, increase figure height, or use `geom_tiplab(align=TRUE)` |
| Scale bar unit matches tree type (substitutions vs. time) | ML tree → substitutions/site; BEAST2 → time units |
| HPD bars visible and correctly placed (BEAST2) | Check `geom_range` column name matches `.beast` object slots |

On any failure → route to `debug`.

## Report

**Mandatory:** Every log file generated during this module must be listed with its exact path in the report so the researcher can monitor background processes and audit what ran.

Write to `reports/[planX/]visualization_YYYY-MM-DD.md`:

```markdown
# Visualization Report
Date: YYYY-MM-DD
Plan: [planA / planB / ...]

## R Script
[Inline script or path to saved .R file]

## Figures Produced
| File | Description | Dimensions | Format |
|------|-------------|-----------|--------|

## Annotation Decisions
[Support value threshold displayed, clade highlights, tip label formatting]

## R Session Info
[Output of sessionInfo() — paste verbatim]

## Software Versions
| Package | Version |
|---------|---------|
| R | x.x.x |
| ggtree | x.x.x |
| ape | x.x.x |
| phytools | x.x.x |
| ggplot2 | x.x.x |
```

Always include `sessionInfo()` output — it captures all package versions in one call.

## Log Files Generated
[List every log file created during this module with its exact path]
[Examples:]
[  figures/tree_figure.pdf     (exported figure)]
[  figures/tree_figure.svg     (vector figure)]
[  scripts/visualization_<date>.R   (R script used to produce figures)]
[If R warnings were captured to a file, list that path here]

## Common Mistakes

| Mistake | Fix |
|---------|-----|
| Plotting before verifying rooting | Always root first; an unrooted tree silently produces a wrong figure |
| Showing all bootstrap values including low ones | Filter with `ifelse(as.numeric(label) >= 70, label, "")` |
| Using PNG for manuscript submission | Use PDF or SVG — journals require vector format for line art |
| Hardcoding node numbers for highlights without checking | Node numbers change when tree is re-read or re-rooted; always verify with `geom_text(aes(label=node))` |
| Skipping `sessionInfo()` in the report | Package versions are needed for reproducibility; one line, no excuse to skip |
| Using `read.tree()` for IQ-TREE files with bootstrap values | Use `treeio::read.iqtree()` — `read.tree()` drops node annotations |
