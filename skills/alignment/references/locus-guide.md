# Locus Guide — Plant Phylogenetic Markers

Read this before choosing MAFFT strategy and trimming parameters in the `alignment` skill.
Organized by marker class: plastid single-copy → plastid whole-genome → nuclear → Angiosperms353.

---

## 1. Standard Plastid Single-Copy Markers

These are PCR-amplifiable loci used for targeted Sanger or amplicon sequencing.
Alignment strategy is always per-marker (never concatenate before aligning).

### rbcL — ribulose-1,5-bisphosphate carboxylase large subunit

| Property | Value |
|----------|-------|
| Genome region | Plastid, large single-copy (LSC) |
| Typical amplified length | 550–650 bp (partial); 1,428 bp (full CDS) |
| Taxonomic resolution | Order–family (poor at species level) |
| MAFFT strategy | `--linsi` (≤200 seqs); `--auto` (>200) |
| Trimming | Usually minimal; avoid aggressive trimming |
| Codon-aware alignment | Recommended for divergent orders; use `--adjustdirectionaccurately` if orientation varies |
| Known issues | Conservative — low variation at species/genus level; use as backbone only |
| Universal primers | rbcLa-F / rbcLa-R (Kress & Erickson 2007); amplifies ~550 bp |
| Alignment expected length | 550–560 bp (partial), 1,400–1,430 bp (full) |

### matK — maturase K

| Property | Value |
|----------|-------|
| Genome region | Plastid, within trnK intron, LSC |
| Typical amplified length | 800–900 bp |
| Taxonomic resolution | Genus–species; better than rbcL at lower levels |
| MAFFT strategy | `--linsi` (≤200 seqs); `--fftnsi` (200–1,000) |
| Trimming | `trimAl -automated1`; matK has variable indel regions |
| Codon-aware alignment | Not recommended — intron context complicates CDS boundary |
| Known issues | Primer difficulty in some families (Orchidaceae); alignment has hypervariable regions in some clades |
| Universal primers | matK-390F / matK-1326R (Cuenoud et al. 2002); 3F-KIM / 1R-KIM for difficult groups |
| Alignment expected length | 800–930 bp after trimming |

### trnL-F — trnL intron + trnL-trnF intergenic spacer

| Property | Value |
|----------|-------|
| Genome region | Plastid LSC, between trnL and trnF |
| Typical amplified length | 400–1,100 bp (highly variable) |
| Taxonomic resolution | Family–genus; good for angiosperms |
| MAFFT strategy | `--einsi` (multiple conserved domains separated by large indels); do NOT use `--linsi` |
| Trimming | `trimAl -gappyout`; the spacer region has length polymorphism |
| Codon-aware alignment | No — non-coding |
| Known issues | Large indels in spacer region; highly variable length between taxa; requires `--einsi` or alignment quality degrades significantly |
| Universal primers | c / f (Taberlet et al. 1991); most widely used non-coding plastid marker |
| Alignment expected length | 400–1,100 bp depending on lineage (variable) |

### psbA-trnH — psbA-trnH intergenic spacer

| Property | Value |
|----------|-------|
| Genome region | Plastid LSC, between psbA and trnH genes |
| Typical amplified length | 200–800 bp (extremely variable) |
| Taxonomic resolution | Species–genus; one of the most variable plastid regions |
| MAFFT strategy | `--einsi`; high length variation requires local alignment |
| Trimming | Aggressive trimming needed; many regions are lineage-specific |
| Codon-aware alignment | No — non-coding spacer |
| Known issues | Inverted repeat proximity causes inconsistent amplification; very short in some families (<200 bp); large length variation makes alignment across orders unreliable |
| Universal primers | psbA3-f / trnHf (Tate & Simpson 2003) |
| Alignment expected length | Highly variable; exclude from cross-order analysis |

### ITS — nuclear ribosomal Internal Transcribed Spacer (ITS1 + 5.8S + ITS2)

| Property | Value |
|----------|-------|
| Genome region | Nuclear 18S–26S ribosomal array |
| Typical amplified length | ITS1: 180–300 bp; ITS2: 180–300 bp; Full (ITS1+5.8S+ITS2): 500–750 bp |
| Taxonomic resolution | Species–genus; occasionally family |
| MAFFT strategy | `--einsi` (ITS1/ITS2 have conserved ends and variable core); `--linsi` for closely related taxa |
| Trimming | Trim only ragged ends; core ITS region is informative |
| Codon-aware alignment | No — ribosomal RNA gene (secondary structure-aware alignment preferred for ITS2 specifically) |
| Known issues | Concerted evolution not always complete → multiple divergent copies; paralogs may exist; pseudogenes common in polyploids; intragenomic variation inflates branch lengths |
| Universal primers | ITS4 / ITS5 (White et al. 1990); most widely cited plant marker pair |
| Alignment expected length | 550–750 bp (full); 180–320 bp (ITS1 or ITS2 alone) |
| Special note | Separate ITS1 and ITS2 for secondary structure alignment; use `mfold` or `4SALE` for RNA folding |

### ETS — nuclear ribosomal External Transcribed Spacer

| Property | Value |
|----------|-------|
| Genome region | Nuclear, 5' end of 18S rDNA |
| Typical amplified length | 400–700 bp (partial) |
| Taxonomic resolution | Species–genus; often used with ITS for combined nuclear marker |
| MAFFT strategy | `--einsi` |
| Trimming | `trimAl -automated1` |
| Known issues | Less standardized than ITS; primer binding varies by family; fewer database records |

### rpoB, rpoC1, rpoC2 — RNA polymerase beta subunits

| Property | Value |
|----------|-------|
| Genome region | Plastid LSC |
| Typical amplified length | rpoB: 800–1,000 bp; rpoC1: 500–700 bp; rpoC2: 400–600 bp |
| Taxonomic resolution | Family–order; useful for deep-level phylogenetics |
| MAFFT strategy | `--linsi` (CDS; generally conserved) |
| Trimming | Minimal needed |
| Codon-aware alignment | Yes; code within CDS regions |

### ndhF — NADH dehydrogenase F subunit

| Property | Value |
|----------|-------|
| Genome region | Plastid small single-copy (SSC) |
| Typical amplified length | 700–2,000 bp depending on primers |
| Taxonomic resolution | Family–order |
| MAFFT strategy | `--linsi`; CDS, moderate divergence |
| Trimming | `trimAl -automated1` |
| Known issues | Sometimes combined with rpl32-trnL spacer for more variation |

---

## 2. Whole-Plastome (Genome Skimming) Strategy

When using genome skimming data assembled with GetOrganelle or BWA:

| Approach | Description |
|----------|-------------|
| Full plastome tree | Concatenate LSC + SSC + IR into single matrix; use partition by region |
| CDS extraction | Use `scripts/assembly/extract_cds.py` to pull individual genes from GenBank annotation |
| Recommended partition | LSC / SSC / IR1 / IR2 as minimum; or by individual gene |

**Alignment strategy for whole-plastome CDS:**
1. Extract individual CDS with `extract_cds.py`
2. Align each gene independently (`--linsi` for conserved genes; `--einsi` for spacers)
3. Concatenate with AMAS to build supermatrix
4. Use partitioned model in IQ-TREE (BIC criterion for model selection)

**Key plastome structural features:**
- Inverted repeat (IR): ~25 kb; identical or near-identical copies flanking SSC
- IR genes (rpl2, rpl23, ycf2, ycf15, rps7, rps12, ndhB) will be duplicated — keep only one copy
- Common IR boundary genes: rps19 (LSC/IRA), ndhF (SSC/IRB)
- Exclude IR from substitution rate analyses (rate suppression due to copy correction)

---

## 3. Angiosperms353 Target Enrichment Markers

**Source:** Johnson et al. (2019) *Systematic Biology* 68(4):594–606  
**DOI:** 10.1093/sysbio/syy086

### Overview

| Property | Value |
|----------|-------|
| Total target loci | 353 nuclear protein-coding genes |
| Total target CDS length | ~260,802 bp across all loci |
| Average target length | ~740 bp per locus |
| Taxonomic scope | All angiosperms (flowering plants) |
| Resolution | Species → all angiosperms |
| Probe kit | Available from Arbor Biosciences (myBaits); Kew Genomics probe set |
| Reference sequences | `Angiosperms353_targetSequences.fasta` (GitHub: mossmatters/Angiosperms353) |

### Gene selection criteria

- Single-copy in most angiosperms (low-copy nuclear; LCN)
- At most 30% divergence within 95% of angiosperm sequences at time of design
- Protein-coding (enables codon-aware alignment and cross-taxa comparison)
- 143 loci have organelle-related function (but are nuclear-encoded)
- Bias toward well-represented taxa in 1KP transcriptome data

### Assembly

Use HybPiper v2 with `Angiosperms353_targetSequences.fasta` as the target file:
```bash
bash scripts/assembly/run_hybpiper.sh \
  -r Angiosperms353_targetSequences.fasta \
  -s sample_list.txt -d data/raw -o hybpiper_output
```

### Alignment strategy per Angiosperms353 locus

| Condition | Strategy |
|-----------|----------|
| ≤200 samples, single family | `--linsi` |
| ≤200 samples, order or above | `--linsi` or `--einsi` if indels present |
| >200 samples | `--auto` or `--fftnsi` |
| Any locus with indels >5% of columns | Use `trimAl -automated1` |
| Codon-aware | Recommended: align at protein level with `--amino`, back-translate to nucleotide for divergent orders |

**Protein-level alignment (recommended for divergent clades):**
```bash
# 1. Translate CDS to protein
python scripts/utils/extract_taxa_from_fasta.sh  # ensure only CDS in file
# 2. Align protein
mafft --auto --amino --thread 8 locus_aa.fasta > locus_aa_aligned.fasta
# 3. Back-translate using PAL2NAL or in Biopython
pal2nal.pl locus_aa_aligned.fasta locus_nt.fasta -output fasta > locus_codon_aligned.fasta
```

### Known problematic Angiosperms353 loci

| Issue | Affected loci | Action |
|-------|--------------|--------|
| Low recovery in some families | Varies; check HybPiper stats per locus | Exclude if <20% taxon recovery |
| Multiple paralog copies | ~57 loci excluded from original set for >30% divergence; some may still have copies | Run `hybpiper paralog_retriever`; use ASTRAL-Pro3 if paralogs present |
| Very short target (<400 bp) | Some loci in lower-copy gene families | Flag in alignment stats |
| High indel rate | Loci with rapid evolution between orders | Use `--einsi`; check PIS% |

### Recovery expectations

| Taxonomic scope | Expected recovery |
|-----------------|------------------|
| Same genus | >90% of loci typically |
| Same family | 70–90% |
| Same order | 50–80% |
| Cross-order | 30–70% (highly variable) |
| Monocots vs. eudicots | Can drop to <50% for some loci |

Samples recovering <20% of loci should be investigated before inclusion. Check:
1. Library quality (low reads, poor QC)
2. Correct taxonomic identity (wrong species?)
3. Genome duplication masking hybridization

### Concatenation vs. coalescent for Angiosperms353

- **Use coalescent (wASTRAL)** for any dataset spanning multiple families or orders
- **Use concatenation** only for within-genus or within-family studies where ILS is low
- **Always run both** for cross-order studies; compare topologies
- Concordance factors (IQ-TREE `--gcf`) are recommended to quantify discordance

---

## 4. Choosing Markers by Research Question

| Research question | Recommended markers |
|-------------------|---------------------|
| Species-level barcoding | ITS + rbcL (CBOL standard); add matK if needed |
| Genus-level phylogeny | matK + rbcL + ITS ± trnL-F |
| Family-level phylogeny | Angiosperms353 (preferred) OR matK + rbcL + rpoC1 + ndhF |
| Order-level / deep | Angiosperms353; whole-plastome CDS; multi-gene nuclear |
| Intraspecific / population | trnS-trnG, rpl32-trnL spacers; plastid microsatellites; RADseq |
| Divergence time | Add rpoB, atpB for additional clock signal; whole-plastome |

---

## 5. Read this SKILL.md before calling align_markers.sh

The `alignment` SKILL.md and `scripts/alignment/align_markers.sh` require you to choose a MAFFT strategy. Use this table:

| Marker | Recommended flag |
|--------|----------------|
| rbcL, matK, rpoB, rpoC1, ndhF, atpB | `--linsi` (≤200 seqs) |
| trnL-F, psbA-trnH, ITS, ETS | `--einsi` (variable-length regions) |
| Angiosperms353 single-copy CDS (same family) | `--linsi` |
| Angiosperms353 CDS (cross-order) | `--einsi` or protein-level then back-translate |
| Whole-plastome CDS extracted per gene | `--linsi` per gene |
| Any dataset >200 sequences | `--auto` |
