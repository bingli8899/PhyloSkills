#!/usr/bin/env python3
"""
methods_gen.py — Auto-generate Methods section from provenance JSON logs

Reads provenance JSON files written by pipeline scripts (results/provenance/*.json),
matches each to the appropriate paragraph template in methods-templates.md, substitutes
actual parameter values, and assembles a complete Methods draft.

Usage:
    python scripts/manuscript/methods_gen.py \
        --provenance results/provenance/ \
        --templates skills/manuscript/references/methods-templates.md \
        --journal "Systematic Botany" \
        --output reports/methods_draft_YYYY-MM-DD.md \
        [--section data-acquisition,alignment,inference]

Arguments:
    --provenance  Directory containing pipeline provenance JSON files
    --templates   Path to methods-templates.md
    --journal     Target journal name (for citation style note in output)
    --output      Output Markdown file path
    --section     Comma-separated list of sections to generate (default: all)
                  Options: data-acquisition, assembly, alignment,
                           model-selection, ml-tree, coalescent, bayesian,
                           divergence-dating, visualization

Tool version requirements:
    Python >= 3.8  (no external dependencies)

Provenance JSON schema:
    Each JSON file must contain at minimum:
      {"script": "<name>", "tool": "<name>", "version": "<ver>",
       "date": "<ISO8601>", "parameters": {}, "input_files": {},
       "output_files": {}, "runtime_seconds": <int>, "exit_code": <int>}
"""

import argparse
import json
import sys
import re
from pathlib import Path
from datetime import date


SECTION_ORDER = [
    "data-acquisition",
    "assembly",
    "alignment",
    "model-selection",
    "ml-tree",
    "coalescent",
    "bayesian",
    "divergence-dating",
    "visualization",
]

# Maps script names (from provenance JSON "script" field) to section labels
SCRIPT_TO_SECTION = {
    "download_genbank": "data-acquisition",
    "download_sra": "data-acquisition",
    "assemble_plastome_getorganelle": "assembly",
    "assemble_plastome_bwa": "assembly",
    "run_hybpiper": "assembly",
    "annotate_plastome": "assembly",
    "extract_cds": "assembly",
    "align_markers": "alignment",
    "build_gene_trees": "ml-tree",
    "run_aster": "coalescent",
    "mrbayes": "bayesian",
    "beast2": "divergence-dating",
    "ggtree_plot": "visualization",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate Methods section from provenance JSON logs"
    )
    parser.add_argument("--provenance", required=True,
                        help="Directory with provenance JSON files")
    parser.add_argument("--templates", required=True,
                        help="Path to methods-templates.md")
    parser.add_argument("--journal", default="",
                        help="Target journal name")
    parser.add_argument("--output", required=True,
                        help="Output Markdown file")
    parser.add_argument("--section", default="",
                        help="Comma-separated sections to include (default: all)")
    return parser.parse_args()


def load_provenance(provenance_dir: str) -> dict[str, list[dict]]:
    """Load all JSON provenance files, grouped by section."""
    pdir = Path(provenance_dir)
    if not pdir.exists():
        print(f"WARNING: Provenance directory not found: {pdir}", file=sys.stderr)
        return {}

    section_data: dict[str, list[dict]] = {}
    for jf in sorted(pdir.glob("*.json")):
        try:
            with open(jf) as f:
                data = json.load(f)
        except json.JSONDecodeError as e:
            print(f"WARNING: Could not parse {jf}: {e}", file=sys.stderr)
            continue

        script_name = data.get("script", jf.stem)
        # Strip date suffix if present (e.g., "align_markers_2026-04-13" → "align_markers")
        script_base = re.sub(r"_\d{4}-\d{2}-\d{2}$", "", script_name)
        section = SCRIPT_TO_SECTION.get(script_base, "unknown")

        if section not in section_data:
            section_data[section] = []
        section_data[section].append(data)

    return section_data


def extract_template(templates_path: str, section_header: str) -> str:
    """Extract a template block for a given section from templates.md."""
    try:
        content = Path(templates_path).read_text()
    except FileNotFoundError:
        return f"[TEMPLATE FILE NOT FOUND: {templates_path}]"

    # Find section by header (## SECTION: <name>)
    pattern = rf"## SECTION: {re.escape(section_header)}.*?(?=\n## SECTION:|\Z)"
    match = re.search(pattern, content, re.DOTALL | re.IGNORECASE)
    if match:
        return match.group(0).strip()
    return f"[TEMPLATE NOT FOUND for section: {section_header}]"


def fill_template(template: str, data: dict) -> str:
    """
    Substitute {KEY} placeholders in template with values from provenance data.
    Unfilled placeholders become [FILL IN: KEY].
    """
    params = data.get("parameters", {})
    tool_versions = {
        data.get("tool", "").lower(): data.get("version", ""),
    }
    # Flatten all keys for substitution
    fill_map = {}
    fill_map.update(params)
    fill_map.update(tool_versions)
    # Add common top-level fields
    for k in ("version", "date", "tool", "script"):
        fill_map[k] = data.get(k, "")

    def replacer(m):
        key = m.group(1)
        if key.lower() in {k.lower() for k in fill_map}:
            # Case-insensitive lookup
            for k, v in fill_map.items():
                if k.lower() == key.lower():
                    return str(v) if v else f"[FILL IN: {key}]"
        return f"[FILL IN: {key}]"

    return re.sub(r"\{(\w+)\}", replacer, template)


def generate_section(section_label: str, provenance_list: list[dict],
                     templates_path: str) -> str:
    """Generate one Methods section paragraph."""
    # Map section label to template section header
    section_to_template = {
        "data-acquisition": "Sequence Data Acquisition",
        "assembly": "Plastome Assembly",
        "alignment": "Sequence Alignment",
        "model-selection": "Model Selection",
        "ml-tree": "Maximum Likelihood Tree Inference",
        "coalescent": "Coalescent Species Tree (ASTER)",
        "bayesian": "Bayesian Tree Inference (MrBayes)",
        "divergence-dating": "Divergence Time Estimation (BEAST2)",
        "visualization": "Visualization",
    }

    template_header = section_to_template.get(section_label, section_label)
    template_block = extract_template(templates_path, template_header)

    # Use the first (most recent) provenance file for this section
    # Multiple files in same section → append as separate paragraphs
    paragraphs = []
    for prov in provenance_list:
        # Extract the first ```...``` code block as the template paragraph
        code_blocks = re.findall(r"```\n(.*?)\n```", template_block, re.DOTALL)
        if code_blocks:
            # Use first code block that looks like prose (not bash)
            prose_blocks = [b for b in code_blocks if not b.strip().startswith("#")]
            if prose_blocks:
                para = fill_template(prose_blocks[0], prov)
                paragraphs.append(para)
            else:
                paragraphs.append(f"[TEMPLATE BLOCK FOUND BUT NO PROSE TEMPLATE for {section_label}]")
        else:
            paragraphs.append(f"[NO TEMPLATE CODE BLOCK FOUND for {section_label}]")

    return "\n\n".join(paragraphs)


def count_placeholders(text: str) -> int:
    return len(re.findall(r"\[FILL IN:", text))


def main():
    args = parse_args()
    print(f"# methods_gen.py")
    print(f"# Provenance dir: {args.provenance}")
    print(f"# Templates:      {args.templates}")
    print(f"# Journal:        {args.journal or '(not specified)'}")
    print(f"# Output:         {args.output}")
    print(f"# Date:           {date.today()}")
    print()

    # Determine which sections to generate
    if args.section:
        sections_to_run = [s.strip() for s in args.section.split(",")]
    else:
        sections_to_run = SECTION_ORDER

    # Load provenance
    section_data = load_provenance(args.provenance)
    print(f"Loaded provenance sections: {list(section_data.keys())}")

    # Build output
    output_lines = [
        f"# Methods Section Draft",
        f"",
        f"**Generated:** {date.today()}  ",
        f"**Target journal:** {args.journal or '[specify journal]'}  ",
        f"**Source:** `{args.provenance}`  ",
        f"",
        "> This is an AI-generated draft. All `[FILL IN: ...]` placeholders must be",
        "> resolved before submission. Verify all parameter values against your analysis",
        "> logs. Add biological context sentences (why this group, why these markers).",
        "",
        "---",
        "",
        "## Materials and Methods",
        "",
    ]

    total_placeholders = 0

    for section in sections_to_run:
        if section not in section_data:
            output_lines.append(f"### {section.replace('-', ' ').title()}")
            output_lines.append("")
            output_lines.append(f"> ⚠️  No provenance data found for section `{section}`.")
            output_lines.append(f"> Run the corresponding pipeline script with `--provenance results/provenance/`.")
            output_lines.append("")
            continue

        prov_list = section_data[section]
        # Sort by date if available
        prov_list.sort(key=lambda x: x.get("date", ""), reverse=True)

        section_text = generate_section(section, prov_list, args.templates)
        n_ph = count_placeholders(section_text)
        total_placeholders += n_ph

        section_title = section.replace("-", " ").title()
        output_lines.append(f"### {section_title}")
        output_lines.append("")
        output_lines.append(section_text)
        output_lines.append("")

        if n_ph > 0:
            output_lines.append(f"> ⚠️  {n_ph} placeholder(s) need manual completion in this section.")
            output_lines.append("")

    # Summary
    output_lines.extend([
        "---",
        "",
        "## Software Citations Required",
        "",
        "See `skills/manuscript/references/methods-templates.md` for full citation list.",
        "Verify all tool versions against `executables/software-inventory.md`.",
        "",
        "---",
        "",
        f"**Total placeholders remaining:** {total_placeholders}",
        "",
    ])

    # Write output
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(output_lines))
    print(f"\nMethods draft written to: {out_path}")
    print(f"Placeholders remaining: {total_placeholders}")

    if total_placeholders > 0:
        print("Run grep '[FILL IN:' on the output file to find all placeholders.")


if __name__ == "__main__":
    main()
