# Software Inventory

Record installed tool versions here after running `--version` checks.
This file is the source of truth for versions cited in the Methods section.
`methods_gen.py` falls back to provenance JSON logs, but this file is for manual
overrides and tools not covered by provenance.

## How to populate

Run each tool's version command and paste the output below.
Update whenever you upgrade a tool.

```bash
mafft --version 2>&1
iqtree2 --version 2>&1 | head -2
trimal --version 2>&1
python3 -c "import Bio; print(Bio.__version__)"
java -jar ~/bin/ASTER/bin/wastral.jar --version 2>&1 || echo "ASTER: check manually"
mb --version 2>&1 | head -3
beast -version 2>&1 | head -3
```

## Installed Versions

| Tool | Version | Install date | Notes |
|------|---------|--------------|-------|
| MAFFT | | | |
| IQ-TREE | | | |
| trimAl | | | |
| AMAS | | | |
| ASTER / wASTRAL | | | |
| ASTRAL-Pro3 | | | |
| MrBayes | | | |
| BEAST2 | | | |
| GetOrganelle | | | |
| HybPiper | | | |
| BWA | | | |
| SAMtools | | | |
| Python | | | |
| Biopython | | | |
| R | | | |
| ggtree | | | |

## Conda / Container

If using conda or a container, record the environment name and freeze file:

```bash
# To freeze current conda env:
conda env export > environment.yml
```

Environment name: [FILL IN]
Freeze file: `environment.yml` (not committed if large)
