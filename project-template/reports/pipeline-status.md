# Pipeline Status

Track progress through the phylogenetic pipeline.
Update this file as each step completes.

**Project:** [FILL IN: taxon group and study question]
**Last updated:** [FILL IN: date]

---

## Progress

| Step | Script | Status | Date | Notes |
|------|--------|--------|------|-------|
| 1. Data acquisition | `download_genbank.sh` / `download_sra.sh` | pending | | |
| 2. Plastome assembly | `assemble_plastome_getorganelle.sh` | pending | | |
| 3. CDS extraction | `extract_cds.sh` | pending | | |
| 4. Alignment | `align_markers.sh` | pending | | |
| 5. Gene tree inference | `build_gene_trees.sh` | pending | | |
| 6. Coalescent species tree | `run_aster.sh` | pending | | |
| 7. Bayesian inference | (MrBayes / BEAST2) | pending | | |
| 8. Visualization | `ggtree_plot.R` | pending | | |
| 9. Methods draft | `methods_gen.py` | pending | | |
| 10. Figure captions | `fig_caption_gen.py` | pending | | |

Status values: `pending` / `running` / `done` / `failed` / `skipped`

---

## Data Summary

- **Taxa included:** [N]
- **Taxa excluded:** [list with reasons]
- **Markers used:** [list]
- **GenBank accessions:** `data/accessions/genbank_accessions.txt`
- **SRA accessions:** `data/accessions/sra_accessions.txt` (if applicable)

---

## Issues / Decisions

Record any data quality issues, excluded taxa, parameter choices, and why.

| Date | Issue | Decision |
|------|-------|----------|
| | | |

---

## Provenance Logs

Auto-generated logs in `results/provenance/`:

```
ls -lh results/provenance/
```

Run methods draft:
```bash
python $PHYLOSKILLS_ROOT/scripts/manuscript/methods_gen.py \
    --provenance results/provenance/ \
    --templates $PHYLOSKILLS_ROOT/skills/manuscript/references/methods-templates.md \
    --journal "Systematic Botany" \
    --output reports/methods_draft_$(date +%Y-%m-%d).md
```
