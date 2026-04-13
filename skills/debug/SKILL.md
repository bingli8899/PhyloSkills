---
name: debug
description: Use when any pipeline module fails its QC gate, when a tool produces unexpected output, or when the researcher reports a problem at any stage. Diagnoses common failures in phylogenetic pipelines: bad sequence quality, alignment problems, assembly failures, download errors, file naming mismatches, non-convergence, unexpected tree topology, and software errors. Always identifies re-entry point after fix.
---

# Debugging Phylogenetic Pipelines

## Overview

Match the symptom to a category, diagnose the root cause, apply the fix, verify it resolved the issue, then re-enter the pipeline at the correct module. Never proceed past a failure without documenting both the cause and the fix.

## Invocation modes

- **Proactive:** called automatically by a module QC gate failure — the failed module passes its report and the specific check that failed
- **Reactive:** researcher reports a problem — ask which module they were in and what error or unexpected output they saw

In both cases: read the originating module's report first before diagnosing.

---

## Symptom Categories

### 1. Data / Sequence Quality

**Symptom:** BLAST of sequences against expected gene returns no hit, low identity, or a completely different gene  
**Diagnosis:** Misannotated GenBank record, wrong gene in the accession, or chimeric sequence  
**Fix:**
```bash
# Verify sequence against NCBI nt
blastn -query suspect_sequence.fasta -db nt -remote \
  -outfmt "6 qseqid sseqid pident length evalue stitle" \
  -max_target_seqs 5 -out blast_check.txt
```
Remove sequences that do not match the expected gene or are chimeric. Document in the data-acquisition report.  
**Re-entry:** `data-acquisition` (replace sequence) or `alignment` (re-align without it)

---

**Symptom:** Many sequences are very short relative to others for the same marker  
**Diagnosis:** Partial records downloaded from GenBank; common for some markers in under-sequenced groups  
**Fix:** Set a minimum length filter in Entrez search:
```bash
esearch -db nuccore -query '"Genus" AND "matK" AND 500:900[SLEN]' ...
```
Or remove post-download based on length. Flag removed accessions in report.  
**Re-entry:** `data-acquisition`

---

**Symptom:** Alignment has one or a few sequences dramatically misaligned (long rows of gaps or shifted positions)  
**Diagnosis:** Wrong gene, reverse complement issue, or very divergent sequence  
**Fix:**
```bash
# Check orientation
makeblastdb -in reference.fasta -dbtype nucl
blastn -query suspect.fasta -db reference.fasta -outfmt 6
# If hit is on minus strand, reverse complement the sequence
# awk or seqkit: seqkit seq --reverse --complement suspect.fasta
```
Remove or reorient the sequence, then re-run alignment.  
**Re-entry:** `alignment`

---

### 2. File Naming and Format

**Symptom:** Tool fails with "sequence not found", "taxon mismatch", or mismatched counts between files  
**Diagnosis:** Inconsistent naming between FASTA headers and partition/tree files; common when mixing NCBI accessions with species names  
**Fix:**
```bash
# Inspect headers
grep ">" sequences.fasta | head -20

# Standardize to Genus_species_accession_marker with seqkit or awk
seqkit replace -p "^(\\S+).*" -r "{kv}" -k rename_map.tsv sequences.fasta \
  > sequences_renamed.fasta

# Verify all names are unique
grep ">" sequences_renamed.fasta | sort | uniq -d
```
Rebuild any files (partition, taxon list) that reference the old names after renaming.  
**Re-entry:** `alignment` (re-concatenate), `model-selection`, or `tree-inference` depending on where the mismatch was caught

---

**Symptom:** FASTA headers contain spaces, parentheses, colons, or other special characters causing tool crashes  
**Diagnosis:** Special characters are illegal in Newick format and cause silent or cryptic failures in IQ-TREE, RAxML, MrBayes  
**Fix:**
```bash
# Replace spaces and special characters
sed 's/[ ():|,;]/_/g' sequences.fasta > sequences_clean.fasta
```
**Re-entry:** `alignment`

---

### 3. Assembly Failures

**Symptom:** GetOrganelle produces no output or very fragmented assembly  
**Diagnosis (check in order):**
1. Wrong `-F` database flag for the organism
2. Coverage too low (genome skimming reads <5× plastome depth)
3. Seed sequence not matching the reads
4. Contamination dominating the assembly

**Fix:**
```bash
# Verify database flag: embplant_pt, embplant_mt, animal_mt, fungi_mt
# Increase rounds and adjust word size
get_organelle_from_reads.py -1 R1.fastq.gz -2 R2.fastq.gz \
  -F embplant_pt -R 15 -k 21,45,65,85,105 -o output/ -t 8

# Inspect coverage in assembly graph with Bandage
# If fragmented: try NOVOPlasty as alternative
```
**Re-entry:** `assembly`

---

**Symptom:** HybPiper recovery rate <50% for many samples  
**Diagnosis:** Low read depth, poor read quality, wrong target file, or highly divergent sequences  
**Fix:**
```bash
# Check read quality first
fastqc sample_R1.fastq.gz

# Verify target file matches organism group
# If using Angiosperms353 on a non-angiosperm: switch to appropriate bait set

# Check per-gene recovery stats
hybpiper stats -t_dna target.fasta gene --sample_names namelist.txt

# Low-coverage samples: consider excluding or supplementing from GenBank
```
**Re-entry:** `assembly` (re-run with adjusted settings) or `data-acquisition` (supplement with GenBank sequences)

---

### 4. Download Failures

**Symptom:** `efetch` or `prefetch` returns error or empty output  
**Diagnosis:** Network issue, NCBI rate limiting, accession retired or updated, disk space  
**Fix:**
```bash
# Check accession status directly
esearch -db nuccore -query "KJ123456" | efetch -format acc
# If empty: accession may have been superseded — search for updated version

# Add retry with delay for rate limiting
for acc in $(cat accessions.txt); do
  efetch -db nuccore -id $acc -format fasta >> output.fasta
  sleep 0.5   # respect NCBI rate limits: max 3 requests/sec without API key
done

# With API key (recommended for bulk downloads)
export NCBI_API_KEY="your_key_here"
```
**Re-entry:** `data-acquisition`

---

**Symptom:** SRA download completes but FASTQ files are 0 bytes or truncated  
**Diagnosis:** Incomplete prefetch, disk space exhausted, or corrupt SRA cache  
**Fix:**
```bash
# Check disk space
df -h .

# Re-download with validation
prefetch --verify yes SRR123456
fasterq-dump SRR123456 --split-files -O data/raw/

# Verify output
ls -lh data/raw/SRR123456*
```
**Re-entry:** `data-acquisition`

---

### 5. Alignment Quality

**Symptom:** Alignment has >50% gap columns across most sequences after trimming  
**Diagnosis:** Highly divergent sequences mixed with conserved ones, or short partial sequences inflating gaps  
**Fix:** Remove sequences shorter than 50% of the median alignment length, then re-align:
```bash
seqkit seq --min-len $(echo "median_len * 0.5" | bc) aligned.fasta \
  > filtered.fasta
mafft --linsi filtered.fasta > realigned.fasta
```
**Re-entry:** `alignment`

---

**Symptom:** Parsimony-informative sites <5% of alignment length  
**Diagnosis:** Marker is too conserved for the taxonomic level, or too few taxa  
**Fix:** Discuss with researcher — consider adding a more variable marker or expanding taxon sampling. This may require going back to `data-acquisition`.  
**Re-entry:** `research-design` (reconsider marker strategy) or `data-acquisition` (add marker)

---

### 6. Tree Inference Failures

**Symptom:** MrBayes ASDSF does not drop below 0.01 after expected run length  
**Diagnosis (check in order):**
1. Runs stuck in different tree topologies (bimodal posterior)
2. Model overly complex for the data
3. Long-branch taxa dominating topology
4. Insufficient generations

**Fix:**
```bash
# Check in Tracer: are the two runs converging?
# If bimodal: inspect for problematic long-branch taxa and consider removing them
# If model complex: simplify (e.g., unpartitioned run)
# Extend run — add to existing run file rather than starting over:
# In MrBayes: 'mcmc' continues from checkpoint
```
**Re-entry:** `tree-inference`

---

**Symptom:** BEAST2 ESS values <200 for key parameters  
**Diagnosis:** Chain too short, poor mixing, over-parameterized model, or problematic calibrations  
**Fix:** Extend chain (×2–×5 current length). If still failing after extension:
- Simplify the clock model (strict → relaxed or vice versa)
- Check calibration priors are not too tight
- Remove taxa with extreme branch lengths
- Run with fewer partitions

**Re-entry:** `tree-inference`

---

**Symptom:** Tree topology strongly contradicts well-established clades from the literature  
**Diagnosis (check in order):**
1. Outgroup misplacement → long-branch attraction
2. Alignment artifact → misaligned region pulling sequences together
3. Contamination in one or more sequences
4. Genuine novel signal (less common)

**Fix:**
```bash
# Identify suspect long branches
# Plot branch lengths: in R
library(ape)
tree <- read.tree("tree.nwk")
dotTree(tree)   # or barplot(tree$edge.length)

# Remove the longest-branch outlier and re-run
# Check alignment of sequences forming unexpected clade
```
If contamination suspected → BLAST sequences against NCBI nt (see category 1).  
**Re-entry:** `alignment` or `tree-inference`

---

### 7. Software Errors

**Symptom:** Executable not found or `command not found`  
**Fix:** Route immediately to `environment`. Do not attempt workarounds.  
**Re-entry:** `environment`, then back to originating module

---

**Symptom:** IQ-TREE or RAxML crashes with out-of-memory error  
**Fix:**
```bash
# Reduce threads to free RAM
iqtree2 -s concatenated.fasta -p partition.nex -B 1000 -T 4

# Or reduce bootstrap replicates for initial exploratory run
iqtree2 -s concatenated.fasta -p partition.nex -B 100 -T 4
```
**Re-entry:** `tree-inference`

---

**Symptom:** R package not available or version conflict in `visualization`  
**Fix:**
```r
# Check installed version
packageVersion("ggtree")

# Reinstall via Bioconductor
BiocManager::install("ggtree", force = TRUE)

# If conflict with another package, use renv for project isolation
# install.packages("renv"); renv::init()
```
**Re-entry:** `visualization`

---

## Report

Write to `reports/debug_YYYY-MM-DD.md` (append if file exists — multiple debug events in one project):

```markdown
# Debug Report
Date: YYYY-MM-DD
Plan: [planA / planB / ...]

## Event [N]
- **Originating module:** [module name]
- **Trigger:** QC gate failure / researcher-reported
- **Symptom:** [exact error message or QC check that failed]
- **Diagnosis:** [root cause identified]
- **Fix applied:** [commands run or changes made]
- **Outcome:** [resolved / partially resolved / escalated to researcher]
- **Pipeline re-entry point:** [module name]
- **Affected files:** [any files renamed, removed, or replaced]
```

Append a new Event block for each debug session — do not overwrite prior events. The full debug history is part of the project record.
