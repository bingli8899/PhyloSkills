# Methods Paragraph Templates

Each section below is a reusable paragraph stub for `scripts/manuscript/methods_gen.py`.
Placeholders use `{KEY}` syntax. The script substitutes values from provenance JSON logs.
Researcher adds biological context sentences (taxon group, research question) manually.

---

## SECTION: Taxon Sampling

```
A total of {N_TAXA} {TAXON_GROUP} accessions were included in this study, representing {N_GENERA} genera and {N_FAMILIES} families. {OUTGROUP_SENTENCE} Voucher information and GenBank accession numbers are provided in {SUPPLEMENTARY_TABLE}.
```

**Provenance keys:** `n_taxa`, `n_genera`, `n_families`, `outgroup_taxa`  
**Manual additions required:** voucher information, supplementary table reference

---

## SECTION: Sequence Data Acquisition

### GenBank

```
Sequence data were retrieved from GenBank (https://www.ncbi.nlm.nih.gov/nuccore/) using NCBI Entrez Direct {EDIRECT_VERSION} (Wheeler et al. 2003). Searches were conducted for {MARKER_LIST} using the query terms "{SEARCH_QUERY}". A total of {N_SEQUENCES} sequences were downloaded on {DOWNLOAD_DATE}. Sequences were renamed following the convention Genus_species_accession_marker to standardize downstream processing.
```

**Provenance keys:** `edirect_version`, `marker_list`, `search_query`, `n_sequences`, `download_date`

### SRA Raw Reads

```
Raw sequencing reads for {N_SAMPLES} samples were obtained from the NCBI Sequence Read Archive (SRA) using SRA Toolkit {SRATOOLS_VERSION} (Leinonen et al. 2011). Accession numbers are listed in {SUPPLEMENTARY_TABLE}. {DOWNLOAD_MODE_SENTENCE}
```

**Download mode sentences:**
- Bulk: "All reads were downloaded prior to assembly."
- Streaming: "Reads were downloaded and assembled sequentially to conserve storage, with raw reads deleted after successful assembly verification."

**Provenance keys:** `sratools_version`, `n_samples`, `download_mode`, `download_date`

---

## SECTION: Plastome Assembly

### GetOrganelle

```
Plastid genomes were assembled de novo from raw reads using GetOrganelle {GETORGANELLE_VERSION} (Jin et al. 2020) with the database flag {PLANT_TYPE_FLAG} and k-mer sizes {KMER_LIST}. Assembly graphs were inspected in Bandage (Wick et al. 2015) for fragmented assemblies. Assemblies shorter than {MIN_LENGTH_BP} bp were excluded from downstream analyses.
```

**Provenance keys:** `getorganelle_version`, `plant_type_flag`, `kmer_list`, `min_length_bp`

### BWA Reference-Guided Assembly

```
Plastid genome sequences were assembled by mapping raw reads to a reference plastome ({REFERENCE_TAXON}; {REFERENCE_ACCESSION}) using BWA {BWA_VERSION} (Li & Durbin 2009). Mapped reads were sorted and indexed with SAMtools {SAMTOOLS_VERSION} (Danecek et al. 2021), and consensus sequences were called using BCFtools {BCFTOOLS_VERSION} with {CONSENSUS_TYPE} consensus mode and a minimum depth of {MIN_DEPTH}×. Samples with mean mapping depth below 10× were excluded.
```

**Provenance keys:** `bwa_version`, `samtools_version`, `bcftools_version`, `reference_taxon`, `reference_accession`, `consensus_type`, `min_depth`

### HybPiper (Target Enrichment)

```
Target locus recovery from {BAIT_KIT} enrichment data was performed using HybPiper {HYBPIPER_VERSION} (Johnson et al. 2016; McLay et al. 2021) with the {TARGET_FILE} target file containing {N_TARGET_LOCI} loci. Per-sample recovery statistics were assessed using `hybpiper stats`, and sequences were retrieved with `hybpiper retrieve_sequences`. Samples recovering fewer than {MIN_RECOVERY_PCT}% of target loci were excluded from phylogenetic analyses.
```

**Provenance keys:** `hybpiper_version`, `bait_kit`, `target_file`, `n_target_loci`, `min_recovery_pct`

---

## SECTION: Plastome Annotation and CDS Extraction

```
Plastid genome assemblies were annotated using {ANNOTATION_TOOL} {ANNOTATION_VERSION} ({ANNOTATION_CITATION}). Coding sequences (CDS) for {N_GENES} genes were extracted from GenBank-format annotation files using a custom Python script (extract_cds.py; available at {REPO_URL}). Multi-exon genes were assembled by concatenating exon sequences in order.
```

**Provenance keys:** `annotation_tool`, `annotation_version`, `n_genes`

---

## SECTION: Sequence Alignment

### Standard (per-marker)

```
Sequences for each locus were aligned independently using MAFFT {MAFFT_VERSION} (Katoh & Standley 2013). {STRATEGY_SENTENCE} Alignment quality was assessed using a custom diagnostic script (analyze_alignment.py), and sequences with gap content exceeding {MAX_GAP_PCT}% or length below {MIN_LENGTH_BP} bp were excluded. {TRIMMING_SENTENCE}
```

**Strategy sentences (choose based on provenance `strategy` key):**
- linsi: "The `--localpair --maxiterate 1000` algorithm was used for datasets of ≤200 sequences."
- einsi: "The `--ep 0 --genafpair --maxiterate 1000` algorithm was used to accommodate large insertions between conserved regions."
- fftnsi: "The FFT-NS-2 heuristic algorithm (`--retree 2 --maxiterate 2`) was used for datasets exceeding 200 sequences."
- auto: "MAFFT's automatic strategy selection was applied."

**Trimming sentences:**
- trimAl automated1: "Poorly aligned positions were removed using trimAl {TRIMAL_VERSION} (Capella-Gutiérrez et al. 2009) with the `-automated1` mode."
- gappyout: "Gap-rich columns were removed using trimAl {TRIMAL_VERSION} with the `-gappyout` mode."
- none: "No trimming was applied as alignment quality was assessed as sufficient."

**Provenance keys:** `mafft_version`, `strategy`, `max_gap_pct`, `min_length_bp`, `trimal_version`, `trim_mode`

### Concatenation (supermatrix)

```
Individual locus alignments were concatenated into a supermatrix using AMAS {AMAS_VERSION} (Borowiec 2016). The final matrix comprised {N_TAXA} taxa and {TOTAL_LENGTH} bp across {N_PARTITIONS} partitions, with {MISSING_DATA_PCT}% missing data. Partition boundaries are provided in {SUPPLEMENTARY_FILE}.
```

**Provenance keys:** `amas_version`, `n_taxa`, `total_length`, `n_partitions`, `missing_data_pct`

---

## SECTION: Model Selection

```
The best-fit substitution model for each partition was determined using ModelFinder (Kalyaanamoorthy et al. 2017) implemented in IQ-TREE {IQTREE_VERSION} (Minh et al. 2020) with the {MODEL_CRITERION} criterion. {PARTITION_MERGE_SENTENCE}
```

**Partition merge sentences:**
- With merging: "Partitions were merged using the `--merge rclusterf` algorithm to reduce model complexity."
- Without: "Each partition was assigned its best-fit model independently."

**Provenance keys:** `iqtree_version`, `model_criterion`, `partition_merge`

---

## SECTION: Maximum Likelihood Tree Inference

```
Maximum likelihood phylogenetic analysis was conducted using IQ-TREE {IQTREE_VERSION} (Minh et al. 2020) with the best-fit model {MODEL_STRING} {PARTITION_INFO}. Branch support was evaluated using {BOOTSTRAP_REPS} ultrafast bootstrap replicates (UFBoot; Hoang et al. 2018) and the SH-aLRT test (Guindon et al. 2010) with {ALRT_REPS} replicates. Nodes with UFBoot ≥ 95 and SH-aLRT ≥ 80 were considered well-supported.
```

**Provenance keys:** `iqtree_version`, `model_string`, `partition_info`, `bootstrap_reps`, `alrt_reps`

---

## SECTION: Coalescent Species Tree (ASTER)

### wASTRAL

```
Individual gene trees for {N_GENE_TREES} loci were inferred with IQ-TREE {IQTREE_VERSION} as above. A coalescent species tree was then estimated from the gene trees using wASTRAL {ASTER_VERSION} (Zhang et al. 2023), which weights branches by local posterior probability and branch length information (weighting mode {WEIGHT_MODE}). Local posterior support values ≥ 0.95 were considered well-supported.
```

### ASTRAL-Pro3

```
Gene trees for {N_GENE_TREES} loci (including multi-copy loci) were inferred with IQ-TREE {IQTREE_VERSION}. A coalescent species tree accounting for gene duplication and loss was estimated using ASTRAL-Pro3 {ASTER_VERSION} (Zhang et al. 2023).
```

**Provenance keys:** `n_gene_trees`, `iqtree_version`, `aster_version`, `aster_tool`, `weight_mode`

### Concordance Factors

```
Gene concordance factors (gCF) and site concordance factors (sCF) were computed for all branches of the species tree using IQ-TREE {IQTREE_VERSION} (Minh et al. 2020) to quantify the proportion of gene trees and alignment sites, respectively, that support each branch.
```

---

## SECTION: Bayesian Tree Inference (MrBayes)

```
Bayesian phylogenetic inference was conducted in MrBayes {MRBAYES_VERSION} (Ronquist et al. 2012) using the {MODEL_STRING} model {PARTITION_INFO}. Two independent runs of four Markov chain Monte Carlo (MCMC) chains were run for {NGEN} generations, sampling every {SAMPLEFREQ} generations. Convergence was assessed by the average standard deviation of split frequencies (ASDSF < 0.01) and by examining trace plots and effective sample sizes (ESS > 200 for all parameters) in Tracer {TRACER_VERSION} (Rambaut et al. 2018). The first 25% of samples were discarded as burnin. Nodes with posterior probability ≥ 0.95 were considered well-supported.
```

**Provenance keys:** `mrbayes_version`, `model_string`, `ngen`, `samplefreq`, `asdsf`, `tracer_version`

---

## SECTION: Divergence Time Estimation (BEAST2)

```
Divergence times were estimated in BEAST2 {BEAST2_VERSION} (Bouckaert et al. 2019) using a {CLOCK_MODEL} clock model and a {TREE_PRIOR} tree prior. The substitution model was set to {MODEL_STRING}. {N_CALIBRATIONS} fossil or secondary calibration constraints were applied (Table {CALIBRATION_TABLE}). Markov chain Monte Carlo analysis was run for {CHAIN_LENGTH} generations, sampling every {LOG_EVERY} steps. Convergence was assessed in Tracer {TRACER_VERSION} (ESS > 200 for all parameters). After discarding the first 25% as burnin, a maximum clade credibility (MCC) tree was summarized with TreeAnnotator {TREEANNOTATOR_VERSION} using mean node heights.
```

**Provenance keys:** `beast2_version`, `clock_model`, `tree_prior`, `model_string`, `n_calibrations`, `chain_length`, `log_every`, `tracer_version`, `treeannotator_version`

---

## SECTION: Visualization

```
Phylogenetic trees were visualized and annotated using ggtree {GGTREE_VERSION} (Yu et al. 2017) and ape {APE_VERSION} (Paradis & Schliep 2019) in R {R_VERSION} (R Core Team {R_YEAR}). {SUPPORT_DISPLAY_SENTENCE} Figures were exported at {DPI} dpi in {FORMAT} format at publication dimensions ({WIDTH} mm wide).
```

**Support display sentences:**
- UFBoot: "Ultrafast bootstrap support values (UFBoot ≥ 95) are shown above branches."
- PP: "Posterior probabilities (≥ 0.95) are displayed above branches."
- Both: "Ultrafast bootstrap / posterior probability values are shown above/below branches, respectively."

**Provenance keys:** `ggtree_version`, `ape_version`, `r_version`, `r_year`, `dpi`, `format`, `width`

---

## Software Citations

Include all citations below in the References section. Use this as a checklist.

| Tool | Citation |
|------|---------|
| IQ-TREE 2 | Minh BQ et al. 2020. *Mol Biol Evol* 37:1530–1534 |
| ModelFinder | Kalyaanamoorthy S et al. 2017. *Nat Methods* 14:587–589 |
| UFBoot2 | Hoang DT et al. 2018. *Mol Biol Evol* 35:518–522 |
| SH-aLRT | Guindon S et al. 2010. *Syst Biol* 59:307–321 |
| MAFFT | Katoh K & Standley DM. 2013. *Mol Biol Evol* 30:772–780 |
| trimAl | Capella-Gutiérrez S et al. 2009. *Bioinformatics* 25:1972–1973 |
| AMAS | Borowiec ML. 2016. *PeerJ* 4:e1660 |
| ASTRAL-Pro3 / wASTRAL | Zhang C et al. 2023. *Bioinformatics* 39:btad536 |
| MrBayes | Ronquist F et al. 2012. *Syst Biol* 61:539–542 |
| BEAST2 | Bouckaert R et al. 2019. *PLOS Comput Biol* 15:e1006650 |
| Tracer | Rambaut A et al. 2018. *Syst Biol* 67:901–904 |
| GetOrganelle | Jin JJ et al. 2020. *Genome Biol* 21:241 |
| HybPiper | Johnson MG et al. 2016. *Appl Plant Sci* 4:1600016 |
| BWA | Li H & Durbin R. 2009. *Bioinformatics* 25:1754–1760 |
| SAMtools | Danecek P et al. 2021. *Gigascience* 10:giab008 |
| BCFtools | Danecek P et al. 2021. *Gigascience* 10:giab008 |
| ggtree | Yu G et al. 2017. *Mol Biol Evol* 34:1exxx |
| ape | Paradis E & Schliep K. 2019. *Bioinformatics* 35:526–528 |
| Entrez Direct | Wheeler DL et al. 2003. *Nucleic Acids Res* 31:28–33 |
| Angiosperms353 | Johnson MG et al. 2019. *Syst Biol* 68:594–606 |
| RAxML-NG | Kozlov AM et al. 2019. *Bioinformatics* 35:4453–4455 |
