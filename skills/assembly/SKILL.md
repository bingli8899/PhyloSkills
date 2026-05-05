---
name: assembly
description: Use when data-acquisition has identified raw SRA reads (genome skimming, WGS, target enrichment/HybSeq, or transcriptomes) that need assembly or marker extraction before alignment. Runs between data-acquisition and alignment. Use when assembled FASTA sequences are not directly available and must be generated from raw reads using reference-guided or de novo approaches.
---

# Assembling Markers and Genomes from Raw Reads

## Overview

Convert raw SRA reads into assembled sequences ready for alignment. The pipeline has six
stages that apply to all data types:

```
Stage 0 → Test assembly (2–5 samples) + QC + human approval  ← MANDATORY before full run
Stage 1 → Read QC and cleaning      (fastp — mandatory before any assembly)
Stage 2 → Select reference
Stage 3 → Assembly                  (GetOrganelle / HybPiper / BWA / Trinity)
Stage 4 → Scaffold incomplete assemblies  (RagTag — plastome only, when needed)
Stage 5 → Standardize orientation   (Dnaapler — plastome only, all samples)
```

**Do not skip Stage 0.** Run a small pilot before committing to the full dataset.
Catching parameter errors, naming issues, or tool failures on 2–5 samples costs minutes;
catching them halfway through 100 samples costs days.

**Do not skip Stage 1.** Adapter contamination and low-quality reads degrade assembly
continuity and introduce false variants. All downstream tools assume clean input.

---

## Stage 0 — Test Assembly, QC, and Human Approval

**This stage is mandatory before the full-dataset run.** Select 2–5 representative
samples and run the complete pipeline (Stages 1–5) on them. Inspect all outputs,
run QC, and report results to the researcher. **Do not proceed to the full run
without explicit human approval.**

### Sample selection

Choose samples that represent the range of your dataset:
- Include 1–2 smallest runs (fastest to test)
- Include 1 run from a divergent clade if present
- Exclude known-problematic samples (very large runs, flagged accessions) from the pilot

### Run the pilot

Run the complete pipeline script on the selected samples, using the same parameters
and environment variables that will be used for the full run:

```bash
for RUN in <run1> <run2>; do
    prefetch "$RUN" -O <raw_dir>/
    fasterq-dump <raw_dir>/$RUN/ -O <raw_dir>/ --split-files --threads $THREADS
    bash <assembly_script> "$RUN" <raw_dir>/${RUN}_1.fastq <raw_dir>/${RUN}_2.fastq
done
```

### QC checks after pilot

Run all of the following and include every result in the report to the researcher:

**1. Assembly length and type**
```bash
for f in <output_dir>/*.fasta; do
    name=$(basename "$f" .fasta)
    len=$(awk '/^>/{next}{l+=length($0)}END{print l+0}' "$f")
    atype=$(grep -o "assembly_type=[^ ]*" "$f" | head -1)
    echo "$name  ${len} bp  $atype"
done
```

Expected range for complete plastomes: 120,000–220,000 bp. Flag anything outside this range.

**2. FASTA header format**
```bash
grep "^>" <output_dir>/*.fasta
```
All headers must follow the project naming convention (`>Genus_species_RunID`).
Report any that do not match.

**3. Missing data / N content**
```bash
python3 - << 'EOF'
import glob, re
for f in sorted(glob.glob("<output_dir>/*.fasta")):
    seq = "".join(l.strip() for l in open(f) if not l.startswith(">"))
    n_pct = seq.upper().count("N") / max(len(seq), 1) * 100
    gap_pct = seq.count("-") / max(len(seq), 1) * 100
    print(f"{f.split('/')[-1]:<50}  N={n_pct:.1f}%  gaps={gap_pct:.1f}%")
EOF
```
Flag assemblies with > 10% N (heavily scaffolded or failed closure).

**4. Orientation check (psbA position)**
```bash
for f in <output_dir>/*.fasta; do
    name=$(basename "$f" .fasta)
    pos=$(blastn -query <psba_anchor_nt.fasta> -subject "$f" \
                 -outfmt "6 sstart send" -num_alignments 1 2>/dev/null | head -1)
    echo "$name  psbA at: ${pos:-NOT_FOUND}"
done
```
All plastomes should show psbA within the first 500 bp. `NOT_FOUND` = orientation failed;
flag for manual inspection.

**5. fastp read QC pass rate**
```bash
python3 - << 'EOF'
import json, glob
for f in glob.glob("<fastp_qc_dir>/*_fastp.json"):
    d = json.load(open(f))
    in_r  = d["summary"]["before_filtering"]["total_reads"]
    out_r = d["summary"]["after_filtering"]["total_reads"]
    q30   = d["summary"]["after_filtering"]["q30_rate"] * 100
    pct   = out_r / in_r * 100
    flag  = "WARN" if pct < 70 or q30 < 60 else "OK"
    print(f"{flag}  {f.split('/')[-1].replace('_fastp.json',''):<40}  {pct:.1f}% pass  Q30={q30:.1f}%")
EOF
```

**6. Provenance / log check**

Confirm per-sample provenance JSON exists and is complete:
```bash
ls <provenance_dir>/assembly_*.json | wc -l  # should equal number of test samples
cat <provenance_dir>/assembly_<RUN>.json      # spot-check one
```

### Decision table

| QC outcome | Action |
|-----------|--------|
| All checks PASS | Report to researcher → await approval → proceed to full run |
| Length out of range for 1 sample | Investigate that sample; re-run with adjusted parameters if needed; re-QC |
| Length out of range for all samples | Stop; pipeline parameters are wrong; fix and re-run pilot |
| Wrong header format | Fix naming convention in script; re-run pilot |
| N content > 10% in most samples | Scaffolding or assembly is incomplete; tune parameters |
| Orientation failed | Fix anchor file or Dnaapler parameters; re-run Stage 5 only |
| fastp < 70% pass rate | Investigate raw data quality; may be expected for old or degraded libraries |

### Report template for human approval

Present results in this format before requesting approval:

```
## Test Assembly QC Report
Samples tested: <N> of <total>
Date: YYYY-MM-DD

| Sample | Run | Length (bp) | Assembly type | N% | psbA pos | fastp pass% | Flag |
|--------|-----|------------|--------------|-----|----------|-------------|------|

Issues found:
- [list any flags or anomalies]

Recommendation: [PROCEED / DO NOT PROCEED — reason]
```

**Do not start the full run until the researcher explicitly approves.**

---

---

## Strategy Selection

| Data type | Goal | Primary tool | Alternative |
|-----------|------|-------------|-------------|
| Genome skimming / WGS | Full plastome | GetOrganelle | NOVOPlasty |
| Genome skimming / WGS | Nuclear markers from reference | BWA + SAMtools consensus | Bowtie2 + SAMtools |
| Target enrichment (HybSeq, Angiosperms353, custom baits) | Target gene recovery | HybPiper | Captus |
| Transcriptome (RNA-seq) | Gene extraction | Trinity → BLAST+ extraction | — |
| Mixed (multiple data types in same plan) | Per-sample strategy | Apply per row | — |

Read `reports/[planX/]data-acquisition_YYYY-MM-DD.md` to confirm data type per sample
before choosing strategy.

---

## Stage 1 — Read QC and Cleaning (mandatory)

Run **fastp** on every sample before any assembly step. fastp performs adapter trimming,
quality filtering, polyG tail removal (Illumina NextSeq/NovaSeq artifact), and per-sample
QC reports in one pass.

### Paired-end WGS / genome skimming

```bash
fastp \
  -i "${INPUT_DIR}/${forward}" \
  -o "${CLEANED_DIR}/${name}_cleaned_1.fastq.gz" \
  -I "${INPUT_DIR}/${reverse}" \
  -O "${CLEANED_DIR}/${name}_cleaned_2.fastq.gz" \
  --detect_adapter_for_pe \
  --trim_poly_g \
  --length_required 100 \
  --thread "$THREADS" \
  -h "${QC_DIR}/${name}_fastp.html" \
  -j "${QC_DIR}/${name}_fastp.json"
```

### Paired-end target enrichment (HybSeq / Angiosperms353)

Same command; adapter detection is automatic. Lower `--length_required` to 50 for
short amplicon libraries if needed.

### Single-end reads

```bash
fastp \
  -i "${INPUT_DIR}/${sample}.fastq.gz" \
  -o "${CLEANED_DIR}/${name}_cleaned.fastq.gz" \
  --trim_poly_g \
  --length_required 100 \
  --thread "$THREADS" \
  -h "${QC_DIR}/${name}_fastp.html" \
  -j "${QC_DIR}/${name}_fastp.json"
```

### Long reads (PacBio / Oxford Nanopore)

fastp does not support long reads. Use **NanoFilt** (ONT) or **HiFiAdapterFilt** (PacBio):

```bash
# ONT — filter by quality and length
gunzip -c raw.fastq.gz | NanoFilt -q 8 -l 500 | gzip > cleaned.fastq.gz

# PacBio HiFi — remove adapter-only reads
bash HiFiAdapterFilt.sh -p sample -t "$THREADS"
```

### QC review before proceeding

Check the fastp JSON/HTML reports:

| Metric | Flag if |
|--------|---------|
| Reads passing filter | < 70% of input reads |
| Mean quality after filtering | < Q20 |
| Adapter content before filtering | > 30% (suggests unusual library prep) |
| Insert size mode | < 100 bp (very short fragments — assembly may fail) |
| % bases ≥ Q30 | < 60% |

If any sample fails the filter-rate check (< 70% reads passing), investigate before
proceeding — low yield may indicate a contaminated or degraded library.

---

## Stage 2 — Select and Obtain a Reference

Reference quality determines assembly quality. Assess availability in this order:

1. **Published assembly for the same species** — best; use directly
2. **Published assembly for a congener** — reliable for plastome and most nuclear markers
3. **Published assembly for the same family** — acceptable; flag divergence
4. **Published assembly for the same order** — use with caution; expect lower recovery
5. **No close reference** — for plastome: use GetOrganelle with a distant seed and
   built-in database; alert researcher to expected scaffold output

For **GetOrganelle**: built-in databases cover most plant plastomes; specify the correct
`-F` flag (see Stage 3).

For **HybPiper**: the target file is the probe/bait set
(e.g., `mega353.fasta` for Angiosperms353, or a custom bait file).

---

## Stage 3 — Assembly

Use **cleaned reads** from Stage 1 as input to all assemblers.

### Plastome assembly (GetOrganelle) — decision tree

For whole-plastome analyses, only **complete circular assemblies** proceed to orientation
and alignment. Scaffolded outputs require researcher review before use. Follow this
decision tree for every sample:

```
Step 1 — Standard run
  get_organelle_from_reads.py ... -k 21,45,65,85,105,127
  → complete circular  ──→ Stage 5 (orient) ── DONE (exit 0)
  → scaffold / failed  ──→ Step 2

Step 2 — Retry A: smaller word size (-w 65 -R 20)
  get_organelle_from_reads.py ... -k 21,45,65,85,105,127 -w 65 -R 20
  → complete circular  ──→ Stage 5 (orient) ── DONE (exit 0)
  → scaffold / failed  ──→ Step 3

Step 3 — Retry B: larger word size (-w 127 -R 20)
  get_organelle_from_reads.py ... -k 21,45,65,85,105,127 -w 127 -R 20
  → complete circular  ──→ Stage 5 (orient) ── DONE (exit 0)
  → total bases ≥ 100,000 bp  ──→ save scaffold → REPORT TO RESEARCHER (exit 2)
  → total bases <  100,000 bp  ──→ REPORT AS FAILED (exit 1)
```

**Only complete circular assemblies are oriented and included automatically.**
Long scaffolds (≥ 100,000 bp) are saved but require researcher review and manual
reference-guided assembly before alignment. The AI agent does not perform this step.
Failed assemblies (< 100,000 bp after all retries) are excluded pending researcher decision.

**Standard run:**
```bash
get_organelle_from_reads.py \
  -1 "${CLEANED_DIR}/${name}_cleaned_1.fastq.gz" \
  -2 "${CLEANED_DIR}/${name}_cleaned_2.fastq.gz" \
  -F embplant_pt \
  -o "${ASSEMBLY_DIR}/${name}" \
  -t "$THREADS" \
  -k 21,45,65,85,105,127
```

**Retry A (smaller word size — increases seed sensitivity for low coverage):**
```bash
get_organelle_from_reads.py ... -k 21,45,65,85,105,127 -w 65 -R 20
```

**Retry B (larger word size — reduces spurious seeds in complex/contaminated libraries):**
```bash
get_organelle_from_reads.py ... -k 21,45,65,85,105,127 -w 127 -R 20
```

**Stage 3B — BWA reference-guided fallback (if all GetOrganelle attempts fail):**

When all three GetOrganelle attempts fail to produce a complete circular plastome,
fall back to reference-guided consensus assembly using BWA + BCFtools:

```bash
bash scripts/assembly/assemble_plastome_bwa.sh \
  -1 cleaned_R1.fastq.gz -2 cleaned_R2.fastq.gz \
  -r reference_plastome_oriented.fasta \
  -o bwa_outdir/ -s Sample_RunID \
  -t "$THREADS" -d 3 -c strict
```

BWA output is labeled `assembly_type=bwa_reference_guided` in the FASTA header,
distinguishing it from de novo GetOrganelle assemblies. These sequences:
- Are reference-dependent (any regions absent from the reference will be N-masked)
- May have lower accuracy at divergent regions
- Should be noted separately in the assembly report

**Cleaned reads: always deleted at script exit**, regardless of assembly outcome.
Use a shell `trap` to guarantee deletion even on early exit:
```bash
cleanup() { rm -rf "$CLEANED_DIR"; }
trap cleanup EXIT
```
This prevents large gzip-compressed cleaned read files from accumulating when running
many samples in sequence.

**Key `-F` flags by organism:**

| Organism group | Flag |
|---------------|------|
| Land plant plastome | `embplant_pt` |
| Land plant mitochondria | `embplant_mt` |
| Animal mitochondria | `animal_mt` |
| Fungal mitochondria | `fungus_mt` |

Inspect output in Bandage if the assembly graph is complex:
```bash
bandage image "${ASSEMBLY_DIR}/${name}/${name}.fastg" "${name}_graph.png"
```

### Plastome assembly (NOVOPlasty — alternative)

Seed-and-extend; requires a short seed sequence (~500 bp from a closely related species):

```bash
perl NOVOPlasty.pl -c config.txt
# config.txt key settings:
#   Project name = <sample_name>
#   Seed input   = seed_rbcL.fasta      # ~500 bp from close relative
#   Genome range = 120000-220000        # expected plastome size range
#   K-mer        = 33
#   Read length  = <avgLength>
#   Insert size  = <insertSize>
#   Forward reads = cleaned_1.fastq.gz
#   Reverse reads = cleaned_2.fastq.gz
```

### Nuclear marker extraction (reference-guided)

```bash
# Index reference
bwa index reference_markers.fasta

# Map and extract consensus per sample
bwa mem -t "$THREADS" reference_markers.fasta \
    "${CLEANED_DIR}/${name}_cleaned_1.fastq.gz" \
    "${CLEANED_DIR}/${name}_cleaned_2.fastq.gz" \
  | samtools sort -@ "$THREADS" -o "${name}.bam"
samtools index "${name}.bam"
samtools consensus -f fasta "${name}.bam" > "${name}_markers_consensus.fasta"
```

Minimum mapping depth for a usable consensus: ≥5×; flag samples below 10×.

### Target enrichment assembly (HybPiper)

```bash
# Assemble all samples
hybpiper assemble -t_dna target_file.fasta \
  -r "${CLEANED_DIR}/${name}_cleaned_1.fastq.gz" \
     "${CLEANED_DIR}/${name}_cleaned_2.fastq.gz" \
  --prefix "$name" --cpu "$THREADS"

# Retrieve sequences across all samples
hybpiper retrieve_sequences dna -t_dna target_file.fasta \
  --sample_names namelist.txt
```

Review recovery rates:
```bash
hybpiper stats -t_dna target_file.fasta gene --sample_names namelist.txt
hybpiper recovery_heatmap seq_lengths.tsv
```

Acceptable gene recovery: ≥50% of targets per sample. Flag samples below this threshold.

### Transcriptome gene extraction (Trinity + BLAST+)

```bash
Trinity --seqType fq \
  --left  "${CLEANED_DIR}/${name}_cleaned_1.fastq.gz" \
  --right "${CLEANED_DIR}/${name}_cleaned_2.fastq.gz" \
  --max_memory 50G --CPU "$THREADS" --output "${name}_trinity/"

makeblastdb -in reference_genes.fasta -dbtype nucl
blastn -query "${name}_trinity/Trinity.fasta" -db reference_genes.fasta \
  -outfmt 6 -evalue 1e-10 -out "${name}_blast.txt"
```

---

## Stage 4 — Scaffold Incomplete Plastome Assemblies (Researcher-Only)

**This stage is NOT performed automatically by the AI agent.** It is a manual step
performed by the researcher after reviewing the Stage 3 failure report.

The AI pipeline saves long scaffolds (≥ 100,000 bp total) to a dedicated directory and
reports them to the researcher. The researcher then decides whether to:
1. Perform reference-guided assembly manually using RagTag + nucmer (see below)
2. Use the raw scaffold directly (with gaps) for alignment
3. Exclude the sample

The AI agent reports the scaffold path and length; the researcher makes the final decision
and performs any manual assembly steps.

**When the researcher requests scaffolding**, use the tools below.

This stage applies **only to plastome assemblies** where GetOrganelle produced scaffold
contigs (multiple sequences) rather than a single circular assembly after all automated
retries (Stage 3 decision tree). Skip this stage if GetOrganelle produced a complete
circular plastome — those proceed directly to Stage 5.

### When to scaffold vs. accept scaffolds as-is

| Assembly output | Action |
|----------------|--------|
| Single sequence ≥ 90% expected length, complete circular | **Accept** — proceed to Stage 5 |
| 2–5 large contigs (total ≥ 80% expected length) | **Scaffold** with RagTag |
| Many small contigs (total < 80% expected length) | **Attempt scaffolding** — if still < 50% expected, flag for exclusion |
| No output / failed assembly | Retry with adjusted `-R`, `-k`, `-w`; then try NOVOPlasty; if still fails, exclude |

### Reference-guided scaffolding (RagTag + nucmer)

RagTag uses nucmer (MUMmer4) alignment to order and orient contigs against a reference
plastome. Use the closest available reference (same genus preferred).

```bash
# 1. Scaffold contigs against reference
ragtag.py scaffold \
  reference_plastome.fasta \
  "${ASSEMBLY_DIR}/${name}_contigs.fasta" \
  -o "${SCAFFOLD_DIR}/${name}" \
  -t "$THREADS" \
  --aligner nucmer

# Output: ${SCAFFOLD_DIR}/${name}/ragtag.scaffold.fasta
# Unplaced contigs: ${SCAFFOLD_DIR}/${name}/ragtag.scaffold.agp
```

RagTag inserts 100 Ns at each gap between scaffolded contigs. The resulting sequence
contains the correct order and orientation even if gaps remain.

**Inspect the AGP file** to understand scaffold quality:
```bash
# Count gaps and placed contigs
grep -v "^#" "${SCAFFOLD_DIR}/${name}/ragtag.scaffold.agp" | \
  awk '{print $5}' | sort | uniq -c
# W = placed contig; N = gap; U = gap (size unknown)
```

### Check for inverted repeat (IR) integrity

Plastomes have two copies of an inverted repeat (IRa and IRb) flanking the large
single-copy (LSC) and small single-copy (SSC) regions. After scaffolding, verify
both IR copies are present:

```bash
# BLAST the assembly against itself — both IRs should appear as hits of ~25–30 kb
makeblastdb -in "${SCAFFOLD_DIR}/${name}/ragtag.scaffold.fasta" -dbtype nucl
blastn -query "${SCAFFOLD_DIR}/${name}/ragtag.scaffold.fasta" \
       -db    "${SCAFFOLD_DIR}/${name}/ragtag.scaffold.fasta" \
       -outfmt 6 -perc_identity 95 -qcov_hsp_perc 10 \
       -out "${name}_selfblast.txt"
# Expect: 2 hits of ~25,000–30,000 bp on opposite strands (the two IR copies)
```

If only one IR is present or IR sizes differ substantially from the reference, the
scaffold is partially complete — include but flag as `scaffold_partial` in the assembly report.

### Alternative: Mauve Contig Mover (MCM)

MCM is an alternative to RagTag when contigs are highly divergent from the reference
(e.g., structural rearrangements are suspected). Requires the Mauve GUI or command-line
`progressiveMauve`:

```bash
# Run as command-line (batch mode)
java -Xmx4000m -cp Mauve.jar org.gel.mauve.contigs.ContigOrderer \
  -output "${name}_mauve/" \
  -ref reference_plastome.fasta \
  -draft "${ASSEMBLY_DIR}/${name}_contigs.fasta"
```

MCM produces a reordered FASTA that can be used directly in Stage 5.

---

## Stage 5 — Standardize Plastome Orientation

This stage applies **only to plastome assemblies** (complete or scaffolded). All plastomes
in a dataset must be:

1. **Oriented to the same strand** — plastomes can assemble in two orientations
   due to recombination across the inverted repeats; both represent the same molecule
2. **Linearized at the same gene** — circular sequences need a consistent start point
   for pairwise alignment

Skip this stage for non-plastome data (nuclear markers, transcriptomes, target enrichment).

### Why this matters

Without standardization, MAFFT or any other aligner will produce nonsensical alignments
because sequences that are biologically identical will appear completely different at the
sequence level (one is the reverse complement of the other, or starts at a different position).

### Recommended tool: Dnaapler (2024)

Dnaapler reorients circular assemblies to start at a user-specified anchor gene.
Use `dnaapler custom` to reorient all plastomes to start at **psbA**
(a conserved gene at the LSC/IRb boundary, present in nearly all land plant plastomes):

```bash
# Reorient a single assembly
dnaapler custom \
  -i "${SCAFFOLD_DIR}/${name}/ragtag.scaffold.fasta" \
  -o "${ORIENTED_DIR}/${name}" \
  -p "$name" \
  -c psba_anchor.fasta \
  -t "$THREADS"
# psba_anchor_aa.fasta = psbA **protein** (amino acid) sequence — dnaapler custom requires AA, not nucleotide
# Translate nucleotide psbA CDS: python3 -c "from Bio.Seq import Seq; s=open('psba_nt.fasta').read().split('\n',1)[1].replace('\n',''); print('>psbA_protein\n'+str(Seq(s).translate(to_stop=True)))" > psba_anchor_aa.fasta
```

**Batch reorientation of all plastomes:**

```bash
mkdir -p "$ORIENTED_DIR"
for fasta in "${PLASTOME_DIR}"/*.fasta; do
    name=$(basename "$fasta" .fasta)
    dnaapler custom \
      -i "$fasta" \
      -o "${ORIENTED_DIR}/${name}_tmp" \
      -p "$name" \
      -c psba_anchor.fasta \
      -t "$THREADS" 2>/dev/null
    # Move final FASTA to oriented dir
    mv "${ORIENTED_DIR}/${name}_tmp/${name}_reoriented.fasta" \
       "${ORIENTED_DIR}/${name}.fasta" 2>/dev/null || \
    cp "$fasta" "${ORIENTED_DIR}/${name}.fasta"  # fallback: keep original if Dnaapler fails
done
```

**Extract psbA anchor from reference** (run once per project):

```bash
# Using BLAST to extract psbA from the reference plastome
makeblastdb -in reference_plastome.fasta -dbtype nucl
blastn -query psbA_seed.fasta -db reference_plastome.fasta \
  -outfmt "6 sseqid sstart send" -num_alignments 1 | \
awk '{print $1, $2, $3}' | \
while read id start end; do
  samtools faidx reference_plastome.fasta "${id}:${start}-${end}" > psba_anchor.fasta
done
```

### Alternative: nucmer + custom rotation script

If Dnaapler is unavailable, use nucmer to find the psbA position and rotate with Python:

```bash
# Find psbA position in each assembly
nucmer --mum psba_anchor.fasta assembly.fasta -p "${name}_nucmer"
show-coords -r -T "${name}_nucmer.delta" | tail -n +5 | head -1
# Column 4 (sstart) gives the position to rotate to

# Python rotation (Biopython)
python3 - << 'EOF'
from Bio import SeqIO
import sys
record = next(SeqIO.parse("assembly.fasta", "fasta"))
rotate_pos = int(sys.argv[1])   # position from nucmer
rotated = record[rotate_pos:] + record[:rotate_pos]
rotated.id = record.id
SeqIO.write(rotated, "assembly_rotated.fasta", "fasta")
EOF
```

### Verify orientation after standardization

```bash
# Quick check: all plastomes should start with psbA region
# BLAST psbA anchor against all oriented plastomes; expect hit at position ~1
for f in "${ORIENTED_DIR}"/*.fasta; do
    name=$(basename "$f" .fasta)
    pos=$(blastn -query psba_anchor.fasta -subject "$f" \
                 -outfmt "6 sstart" -num_alignments 1 2>/dev/null | head -1)
    echo "$name  psbA at: ${pos:-NOT_FOUND}"
done
```

All plastomes should show psbA at position 1–100. Any showing `NOT_FOUND` need manual
inspection — psbA may be absent (rare) or the anchor sequence diverged too much.

### Include reference plastomes in orientation step

**Reference plastomes downloaded from GenBank must also be reoriented** — they are
deposited at arbitrary starting positions. Run all plastomes (SRA-assembled + GenBank)
through Stage 5 before alignment.

---

## Stage 6 — Collect and Organize Outputs

After assembly, scaffolding, and orientation, consolidate outputs:

```
data/assembled/
  plastomes/
    oriented/
      Species_name_Accession.fasta      # final oriented plastome, all sources
    raw/
      Species_name_RunID/               # GetOrganelle working dirs (delete after QC)
  markers/
    matK_all_samples.fasta              # for marker-based analyses
    rbcL_all_samples.fasta
```

Apply the naming convention from `data-acquisition`:
```
>Genus_species_Accession_or_RunID
```

Rename FASTA headers immediately — do not carry tool-generated headers into alignment.

---

## Stage 7 — QC

### Assembly QC

| Check | Threshold | Action on failure |
|-------|-----------|-------------------|
| Plastome total length | ≥ 80% of reference length | Flag; scaffold if not already done |
| Assembly type | Complete circular preferred | Note scaffold in report; acceptable for alignment |
| Both IR copies present | Sizes match reference ± 10% | Flag as scaffold_partial |
| Dnaapler orientation | psbA at position 1–100 | Manual check; re-run with adjusted anchor |
| Mean read depth (reference-guided) | ≥ 10× recommended; ≥ 5× minimum | Flag; consider excluding |
| HybPiper gene recovery | ≥ 50% targets per sample | Flag low-recovery samples |
| FASTA header consistency | All match naming convention | Rename before proceeding |
| No empty output files | All expected files > 0 bytes | Debug |

### Read QC check (after fastp)

Review fastp JSON reports programmatically:

```bash
python3 - << 'EOF'
import json, glob, os
for f in glob.glob("qc/*_fastp.json"):
    d = json.load(open(f))
    name = os.path.basename(f).replace("_fastp.json", "")
    in_reads  = d["summary"]["before_filtering"]["total_reads"]
    out_reads = d["summary"]["after_filtering"]["total_reads"]
    pct_pass  = out_reads / in_reads * 100
    q30       = d["summary"]["after_filtering"]["q30_rate"] * 100
    flag = "WARN" if pct_pass < 70 or q30 < 60 else "OK"
    print(f"{flag}  {name:<40}  {pct_pass:.1f}% pass  Q30={q30:.1f}%")
EOF
```

On any failure → route to `debug` before proceeding.

---

## Report

**Mandatory:** Every log file generated during this module must be listed with its exact path in the report so the researcher can monitor background processes and audit what ran.

Write to `reports/[planX/]assembly_YYYY-MM-DD.md`:

```markdown
# Assembly Report
Date: YYYY-MM-DD
Plan: [planA / planB / ...]

## Read QC Summary (fastp)
| Sample | Input reads | Pass rate | Q30 rate | Flag |
|--------|------------|----------|---------|------|

## Assembly Strategy
- Data type: [genome skimming / WGS / HybSeq / transcriptome]
- Assembly tool: [GetOrganelle / HybPiper / BWA+SAMtools / Trinity]
- Reference used: [name, accession, source, divergence notes]
- Scaffolding: [RagTag / none needed]
- Orientation tool: [Dnaapler / nucmer+rotation / none needed]
- Anchor gene for orientation: [psbA / other]

## Per-Sample Results
| Sample | Run/Accession | Assembly length | Assembly type | IR check | Orientation | QC flag |
|--------|--------------|----------------|--------------|----------|-------------|---------|

## Failed / Low-Quality Samples
[Samples excluded or flagged, with reason]

## Output Files
[Paths to oriented plastomes / per-marker FASTA files]

## Software Versions
| Tool | Version |
|------|---------|
| fastp | |
| GetOrganelle | |
| RagTag | |
| Dnaapler | |
| MUMmer / nucmer | |

## Log Files Generated
[List every log file created during this module with its exact path so the researcher
 can monitor background assemblies and review per-sample details]
[Examples:]
[  qc/<sample>_fastp.json         (fastp per-sample QC log)]
[  qc/<sample>_fastp.html         (fastp HTML report)]
[  assemblies/<sample>/get_organelle.log   (GetOrganelle stdout)]
[  hybpiper_output/<sample>.log            (HybPiper per-sample log)]

## Next Module
alignment
```

---

## Scripts

Pre-built scripts in `scripts/assembly/`:

| Script | Purpose |
|--------|---------|
| `assemble_plastome_getorganelle.sh` | GetOrganelle plastome assembly with QC length check |
| `assemble_plastome_bwa.sh` | Reference-guided assembly: BWA → SAMtools → BCFtools consensus |
| `run_hybpiper.sh` | Full HybPiper v2 pipeline: assemble → stats → retrieve_sequences |
| `annotate_plastome.sh` | Plastome annotation with PLANN or chloe (auto-detected) |
| `extract_cds.py` | Extract CDS from GenBank `.gb` annotations; handles multi-exon genes |

Utility scripts in `skills/utils/`:

| Script | Purpose |
|--------|---------|
| `change_header_name.py` | Rename FASTA headers from a two-column TSV mapping |
| `remove_N_fasta.py` | Remove sequences with excessive N or gap content |
| `revise_hybpiper_sequences.py` | Post-process HybPiper output: rename, filter short seqs |

Usage examples:

```bash
# Plastome assembly pipeline (Stage 3)
bash scripts/assembly/assemble_plastome_getorganelle.sh \
  -1 sample_cleaned_1.fastq.gz -2 sample_cleaned_2.fastq.gz \
  -o assemblies -s Species_name -t 8

# Reference-guided scaffolding (Stage 4)
ragtag.py scaffold reference_plastome.fasta contigs.fasta \
  -o scaffolds/Species_name -t 8 --aligner nucmer

# Orientation standardization (Stage 5)
dnaapler custom -i scaffold.fasta -o oriented/Species_name \
  -p Species_name -c psba_anchor.fasta -t 8

# HybPiper pipeline (Stage 3, target enrichment)
bash scripts/assembly/run_hybpiper.sh \
  -r target_file.fasta -s sample_list.txt \
  -d data/cleaned -o hybpiper_output

# Extract CDS from annotation
python scripts/assembly/extract_cds.py \
  --input annotations/ --output data/cds/ --concatenate
```

---

## Annotation and Marker Extraction after Plastome Assembly

For plastome-based analyses where CDS extraction is needed (not whole-plastome alignment):

### Option A — BLAST-based extraction (works on complete and scaffold assemblies)

Use `scripts/data/extract_markers_blast.py`. Requires a reference GenBank file for the
target gene coordinates and a per-marker reference FASTA directory.

Key BLAST settings that must not be changed:
- `-word_size 7` for cross-genus plant plastid searches (~70–75% identity)
- Omit `-max_hsps` — adding `-max_hsps 1` causes empty output for some markers
- Use per-marker `min_cov` thresholds for large genes (ycf1, ycf2): set to ≤ 5%

Skip assemblies shorter than 50,000 bp — too fragmented for reliable extraction.

### Option B — PLANN annotation (complete assemblies only)

PLANN requires a **complete, circular plastome**. It annotates scaffold assemblies but
recovers far fewer genes. Use PLANN only for:
- Complete circular assemblies (≥ 90% of expected length for the target taxon)
- Preparing annotated GenBank files for NCBI submission

Call the Perl script directly to avoid a `FindBin` symlink bug:
```bash
perl /path/to/plann/plann.pl \
  -reference reference.gb \
  -fasta assembly.fasta \
  -out output_prefix \
  -organism "Genus species"
```

PLANN outputs `.tbl` + `.fsa`. Convert to GenBank format with `tbl2asn` if downstream
tools require `.gb`.

---

## Common Mistakes

| Mistake | Fix |
|---------|-----|
| Skipping read QC | Adapter contamination causes assembly artifacts; always run fastp first |
| Running fastp with default `--length_required 15` on WGS | Set `--length_required 100` for genome skimming; 50 for HybSeq |
| Using `--trim_poly_g` without `--detect_adapter_for_pe` | Both flags are needed for Illumina PE data |
| Keeping cleaned reads after assembly | Use `trap cleanup EXIT` to delete cleaned reads at script exit; large runs produce 10–30 GB of cleaned reads that accumulate silently |
| Not labeling BWA fallback assemblies | BWA consensus must carry `assembly_type=bwa_reference_guided` in the FASTA header; mixing with de novo assemblies without labeling creates a silent quality difference |
| Running BWA fallback on circular GetOrganelle assemblies | BWA is only for samples that failed GetOrganelle after all retries; circular assemblies should proceed directly to orientation |
| Using wrong `-F` flag in GetOrganelle | `embplant_pt` for plant plastome; `animal_mt` for animal mito |
| Accepting scaffold assemblies without orientation check | Scaffolds can be on opposite strand; always run Stage 5 |
| Skipping plastome orientation standardization | Aligned plastomes appear completely different if half are reverse-complemented |
| Using psbA from a distantly related reference as anchor | Anchor divergence > 20% causes Dnaapler misidentification; use same-genus psbA |
| Not standardizing GenBank-downloaded plastomes | GenBank plastomes start at arbitrary positions; include them in Stage 5 |
| Mixing oriented and unoriented plastomes in alignment | Silent error — half the alignment will be garbage; verify all are oriented |
| Running PLANN on scaffold assemblies | PLANN only works on complete circular plastomes; use BLAST extraction for scaffolds |
| Using `-max_hsps 1` in blastn | Causes empty BLAST output for some markers even when a hit exists; omit this flag |
| Using default blastn word_size (11) for cross-genus searches | Misses divergent markers (~75% identity); use `-word_size 7` |
| Keeping tool-generated FASTA headers | Rename to `Genus_species_accession` immediately after assembly |
| Skipping Bandage inspection on complex assembly graphs | Fragmentation or repeat tangles are invisible without visual inspection |
| Deleting raw reads before verifying final FASTA exists | Unrecoverable; check file exists and length > 0 before deleting |
