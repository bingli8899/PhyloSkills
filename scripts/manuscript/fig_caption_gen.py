#!/usr/bin/env python3
"""
fig_caption_gen.py — Generate figure captions from tree files and provenance logs

Reads provenance JSON files and tree metadata to produce ready-to-edit figure
captions for phylogenetic trees. Outputs captions formatted for the target journal.

Usage:
    python scripts/manuscript/fig_caption_gen.py \
        --provenance results/provenance/ \
        --trees results/trees/ \
        --output reports/figure_captions_YYYY-MM-DD.md \
        [--journal "Systematic Botany"]

Arguments:
    --provenance  Directory with provenance JSON files
    --trees       Directory containing tree files (.treefile, .nwk, .mcc.tree)
    --output      Output Markdown file for captions
    --journal     Target journal (affects caption style)

Tool version requirements:
    Python >= 3.8  (no external dependencies)
"""

import argparse
import json
import re
from pathlib import Path
from datetime import date
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate figure captions from tree files and provenance logs"
    )
    parser.add_argument("--provenance", required=True)
    parser.add_argument("--trees", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--journal", default="")
    return parser.parse_args()


def load_provenance_by_script(prov_dir: str) -> dict[str, dict]:
    """Load provenance files keyed by script base name."""
    result = {}
    pdir = Path(prov_dir)
    if not pdir.exists():
        return result
    for jf in sorted(pdir.glob("*.json")):
        try:
            data = json.load(open(jf))
        except Exception:
            continue
        script = re.sub(r"_\d{4}-\d{2}-\d{2}$", "", data.get("script", jf.stem))
        result[script] = data
    return result


def count_taxa_in_tree(tree_path: Path) -> int:
    """Count number of leaf taxa in a Newick tree file."""
    try:
        content = tree_path.read_text()
        # Count leaf names: strings before ':', before ',', before ')' that aren't numbers
        leaves = re.findall(r"([A-Za-z][A-Za-z0-9_]+)(?=[:,()\s])", content)
        # Filter internal node labels (usually numbers or short strings)
        taxa = [l for l in leaves if len(l) > 3 and not l.replace(".", "").isdigit()]
        return len(set(taxa))
    except Exception:
        return 0


def find_tree_files(trees_dir: str) -> dict[str, list[Path]]:
    """Find tree files by type."""
    tdir = Path(trees_dir)
    trees = defaultdict(list)
    if not tdir.exists():
        return trees
    for f in sorted(tdir.rglob("*")):
        if f.name.endswith(".treefile"):
            trees["ml"].append(f)
        elif f.name.endswith(".contree"):
            trees["consensus"].append(f)
        elif f.name.endswith("_wastral.nwk") or "wastral" in f.name:
            trees["wastral"].append(f)
        elif f.name.endswith("_astral_pro3.nwk") or "astral_pro" in f.name:
            trees["astral_pro3"].append(f)
        elif f.name.endswith(".mcc.tree") or "mcc" in f.name:
            trees["beast2"].append(f)
        elif f.name.endswith(".con.tre") or f.name.endswith("_mrbayes.nwk"):
            trees["mrbayes"].append(f)
    return trees


def caption_ml_tree(prov: dict, tree_path: Path | None) -> str:
    """Generate ML tree figure caption."""
    params = prov.get("parameters", {})
    iqtree_ver = prov.get("version", "[FILL IN: iqtree_version]")
    model = params.get("model", "[FILL IN: model_string]")
    ufboot = params.get("ufboot", params.get("bootstrap_reps", "1000"))
    alrt = params.get("alrt", params.get("alrt_reps", "1000"))
    n_taxa = count_taxa_in_tree(tree_path) if tree_path else "[FILL IN: n_taxa]"
    partition_info = params.get("partition_info", "")

    partition_text = f" under a partitioned model ({partition_info})" if partition_info else ""
    n_taxa_text = str(n_taxa) if isinstance(n_taxa, int) and n_taxa > 0 else "[FILL IN: n_taxa]"

    return (
        f"Maximum likelihood phylogram of {n_taxa_text} [FILL IN: taxon_group] inferred "
        f"using IQ-TREE {iqtree_ver} (Minh et al. 2020) with the {model} substitution model"
        f"{partition_text}. "
        f"Branch support was assessed using {ufboot} ultrafast bootstrap replicates (UFBoot; "
        f"Hoang et al. 2018) and the SH-aLRT test ({alrt} replicates; Guindon et al. 2010). "
        f"Values above/below branches represent UFBoot/SH-aLRT support; filled circles indicate "
        f"UFBoot ≥ 95 and SH-aLRT ≥ 80 (well-supported nodes). "
        f"[FILL IN: outgroup_rooting_sentence] "
        f"Scale bar represents [FILL IN: scale_bar_units]."
    )


def caption_aster_tree(prov: dict, tree_path: Path | None) -> str:
    """Generate ASTER coalescent species tree caption."""
    params = prov.get("parameters", {})
    aster_ver = prov.get("version", "[FILL IN: aster_version]")
    tool = params.get("tool", prov.get("tool", "wASTRAL"))
    n_gene_trees = params.get("n_gene_trees", "[FILL IN: n_gene_trees]")
    n_taxa = count_taxa_in_tree(tree_path) if tree_path else "[FILL IN: n_taxa]"
    n_taxa_text = str(n_taxa) if isinstance(n_taxa, int) and n_taxa > 0 else "[FILL IN: n_taxa]"

    tool_display = {
        "wastral": "wASTRAL",
        "astral-pro3": "ASTRAL-Pro3",
        "astral": "ASTRAL",
    }.get(str(tool).lower(), str(tool))

    citation = {
        "wastral": "Zhang et al. 2023",
        "astral-pro3": "Zhang et al. 2023",
        "astral": "Zhang et al. 2018",
    }.get(str(tool).lower(), "Zhang et al. 2023")

    return (
        f"Coalescent species tree of {n_taxa_text} [FILL IN: taxon_group] inferred from "
        f"{n_gene_trees} gene trees using {tool_display} {aster_ver} ({citation}). "
        f"Numbers on branches represent local posterior probability support values; "
        f"values ≥ 0.95 indicate well-supported nodes. "
        f"[FILL IN: comparison_with_concatenation_tree_sentence] "
        f"[FILL IN: outgroup_rooting_sentence]"
    )


def caption_beast2_tree(prov: dict, tree_path: Path | None) -> str:
    """Generate BEAST2 chronogram caption."""
    params = prov.get("parameters", {})
    beast_ver = prov.get("version", "[FILL IN: beast2_version]")
    clock_model = params.get("clock_model", "[FILL IN: clock_model]")
    n_calibrations = params.get("n_calibrations", "[FILL IN: n_calibrations]")
    n_taxa = count_taxa_in_tree(tree_path) if tree_path else "[FILL IN: n_taxa]"
    n_taxa_text = str(n_taxa) if isinstance(n_taxa, int) and n_taxa > 0 else "[FILL IN: n_taxa]"

    return (
        f"Time-calibrated maximum clade credibility (MCC) chronogram of {n_taxa_text} "
        f"[FILL IN: taxon_group] inferred using BEAST2 {beast_ver} (Bouckaert et al. 2019) "
        f"with a {clock_model} clock model. "
        f"{n_calibrations} calibration constraints were applied (indicated by filled triangles "
        f"at nodes; see Table [FILL IN: calibration_table_number] for prior distributions and sources). "
        f"Node bars represent 95% highest posterior density (HPD) intervals for divergence times. "
        f"Posterior probability support values are shown above branches; "
        f"values ≥ 0.95 are considered well-supported. "
        f"Time scale in million years (Ma)."
    )


def caption_mrbayes_tree(prov: dict, tree_path: Path | None) -> str:
    """Generate MrBayes tree caption."""
    params = prov.get("parameters", {})
    mb_ver = prov.get("version", "[FILL IN: mrbayes_version]")
    model = params.get("model_string", "[FILL IN: model_string]")
    ngen = params.get("ngen", "[FILL IN: ngen]")
    n_taxa = count_taxa_in_tree(tree_path) if tree_path else "[FILL IN: n_taxa]"
    n_taxa_text = str(n_taxa) if isinstance(n_taxa, int) and n_taxa > 0 else "[FILL IN: n_taxa]"

    return (
        f"Bayesian majority-rule consensus phylogram of {n_taxa_text} [FILL IN: taxon_group] "
        f"inferred using MrBayes {mb_ver} (Ronquist et al. 2012) with the {model} substitution model. "
        f"Two independent runs of four MCMC chains were executed for {ngen} generations. "
        f"Posterior probability values above branches indicate node support; "
        f"values ≥ 0.95 are considered well-supported. "
        f"[FILL IN: outgroup_rooting_sentence] "
        f"Scale bar represents substitutions per site."
    )


def main():
    args = parse_args()
    print(f"# fig_caption_gen.py")
    print(f"# Date: {date.today()}")
    print(f"# Journal: {args.journal or '(not specified)'}")
    print()

    prov = load_provenance_by_script(args.provenance)
    trees = find_tree_files(args.trees)

    output_lines = [
        f"# Figure Captions Draft",
        f"",
        f"**Generated:** {date.today()}  ",
        f"**Journal:** {args.journal or '[specify]'}  ",
        f"",
        "> AI-generated draft. Fill all `[FILL IN: ...]` placeholders before submission.",
        "> Add figure numbers once final figure order is determined.",
        "",
        "---",
        "",
    ]

    fig_num = 1

    # ML tree
    if "build_gene_trees" in prov or trees.get("ml"):
        ml_prov = prov.get("build_gene_trees", {})
        ml_tree = trees["ml"][0] if trees.get("ml") else None
        output_lines.append(f"**Figure {fig_num}.** {caption_ml_tree(ml_prov, ml_tree)}")
        output_lines.append("")
        fig_num += 1

    # ASTER species tree
    if "run_aster" in prov or trees.get("wastral") or trees.get("astral_pro3"):
        aster_prov = prov.get("run_aster", {})
        aster_tree = (trees.get("wastral") or trees.get("astral_pro3") or [None])[0]
        output_lines.append(f"**Figure {fig_num}.** {caption_aster_tree(aster_prov, aster_tree)}")
        output_lines.append("")
        fig_num += 1

    # MrBayes tree
    if trees.get("mrbayes"):
        mb_prov = prov.get("mrbayes", {})
        mb_tree = trees["mrbayes"][0]
        output_lines.append(f"**Figure {fig_num}.** {caption_mrbayes_tree(mb_prov, mb_tree)}")
        output_lines.append("")
        fig_num += 1

    # BEAST2 chronogram
    if "beast2" in prov or trees.get("beast2"):
        b2_prov = prov.get("beast2", {})
        b2_tree = trees["beast2"][0] if trees.get("beast2") else None
        output_lines.append(f"**Figure {fig_num}.** {caption_beast2_tree(b2_prov, b2_tree)}")
        output_lines.append("")
        fig_num += 1

    if fig_num == 1:
        output_lines.append(
            "⚠️  No tree files found in `results/trees/` and no matching provenance logs. "
            "Run pipeline scripts and ensure they write provenance JSON files."
        )

    # Count placeholders
    full_text = "\n".join(output_lines)
    n_placeholders = len(re.findall(r"\[FILL IN:", full_text))
    output_lines.extend([
        "---",
        "",
        f"**Placeholders remaining:** {n_placeholders}",
    ])

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(output_lines))
    print(f"Captions written to: {out_path}")
    print(f"Figures generated: {fig_num - 1}")
    print(f"Placeholders remaining: {n_placeholders}")


if __name__ == "__main__":
    main()
