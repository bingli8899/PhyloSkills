# Sequence Database Field Guide

Reference for the `data-acquisition` skill and `scripts/data/download_genbank.sh`.
Describes each database, what it contains, and how to query it effectively.

---

## 1. NCBI GenBank / Nucleotide

**URL:** https://www.ncbi.nlm.nih.gov/nuccore/  
**Entrez Direct docs:** https://www.ncbi.nlm.nih.gov/books/NBK179288/

### What it contains

- Assembled gene sequences from Sanger, amplicon, and genomic sequencing
- Complete and partial plastomes (cp genomes)
- Whole genome shotgun (WGS) contigs
- RefSeq curated sequences

### Effective search patterns

```bash
# Standard markers for a family
esearch -db nuccore \
  -query '"Zingiberaceae"[Organism] AND ("matK"[All Fields] OR "maturase K"[All Fields]) AND biomol_genomic[PROP]'

# Whole plastomes only (RefSeq or complete sequences)
esearch -db nuccore \
  -query '"Zingiberaceae"[Organism] AND "complete genome"[Title] AND chloroplast[Filter]'

# Multiple markers at once (OR logic)
esearch -db nuccore \
  -query '"Zingiberaceae"[Organism] AND ("matK" OR "rbcL" OR "ITS" OR "psbA") AND biomol_genomic[PROP]'
```

### Useful filters

| Filter | Entrez syntax |
|--------|--------------|
| Genomic DNA only | `biomol_genomic[PROP]` |
| Plastid sequences | `chloroplast[Filter]` or `plastid[Filter]` |
| Length range | `500:2000[SLEN]` |
| Recent submissions | `2020:2026[PDAT]` |
| Exclude partial sequences | Add `NOT partial[Title]` (imperfect — use length filter too) |

### Output formats

```bash
efetch -format fasta        # FASTA sequence
efetch -format gb           # GenBank flat file (for annotation)
efetch -format docsum       # Document summary (for metadata)
efetch -format acc          # Accession list only
```

---

## 2. NCBI SRA (Sequence Read Archive)

**URL:** https://www.ncbi.nlm.nih.gov/sra/  
**Used for:** Raw reads from genome skimming, WGS, HybSeq, RNA-seq

### Library strategy keywords

| Strategy | SRA term | Route to |
|----------|----------|----------|
| Genome skimming / WGS | `WGS` | `assembly` (GetOrganelle or BWA) |
| Target enrichment (HybSeq) | `Targeted-Capture` | `assembly` (HybPiper) |
| RNA-seq / transcriptome | `RNA-Seq` | `assembly` (Trinity + BLAST) |
| Amplicon / barcode | `AMPLICON` | `alignment` directly (assembled) |

### Search queries

```bash
# Genome skimming for a group
esearch -db sra \
  -query '"Zingiberaceae"[Organism] AND ("genome skimming"[Strategy] OR "WGS"[Strategy])'

# Target enrichment (Angiosperms353 or other HybSeq)
esearch -db sra \
  -query '"Zingiberaceae"[Organism] AND "target enrichment"[Strategy]'

# Estimate sizes before downloading
esearch -db sra -query "SRR12345678" | efetch -format runinfo \
  | cut -d',' -f1,7,10   # Run, size_MB, LibraryStrategy
```

### Download (see download_sra.sh)

```bash
# Bulk mode (storage sufficient)
bash scripts/data/download_sra.sh -l sra_accessions.txt -o data/raw/ -m bulk

# Streaming mode (storage limited — download→assemble→delete)
bash scripts/data/download_sra.sh -l sra_accessions.txt -o data/raw/ -m streaming \
  -a scripts/assembly/assemble_plastome_getorganelle.sh
```

---

## 3. BOLD Systems (Barcode of Life Data System)

**URL:** https://www.boldsystems.org  
**Contains:** COI (animals), rbcL + matK (plants), ITS (fungi) barcoding records

### When to use

- Fast species-level identification and barcode comparison
- Cross-checking NCBI sequences against BOLD voucher specimens
- Filling gaps in rbcL / matK coverage for angiosperms

### Query approach

BOLD has a public API:
```bash
# Download all rbcL records for a genus (replace spaces with +)
curl "https://www.boldsystems.org/index.php/API_Public/sequence?taxon=Zingiber&marker=rbcL" \
  > bold_zingiber_rbcL.fasta

# matK
curl "https://www.boldsystems.org/index.php/API_Public/sequence?taxon=Zingiberaceae&marker=matK" \
  > bold_zingiberaceae_matK.fasta
```

### Caveats

- BOLD sequences are not deposited in GenBank by default — treat as supplemental
- Sequence quality varies; always check length and N content
- Headers are in BOLD format; rename to `Genus_species_BOLDID_marker` format before merging
- May have taxonomic misidentifications at variety/subspecies level

---

## 4. TreeBASE

**URL:** https://treebase.org  
**Contains:** Published alignment matrices and tree files from peer-reviewed phylogenetic studies

### When to use

- Reuse published alignments (avoids re-aligning well-curated data)
- Benchmark alignment quality against published work
- Supplement with newly sequenced taxa

### Search

Search at https://treebase.org/treebase-web/search/studySearch.html  
Or use the R package `treebase`:
```r
library(treebase)
results <- search_treebase("Zingiberaceae", by="taxon")
```

---

## 5. GBIF (Global Biodiversity Information Facility)

**URL:** https://www.gbif.org  
**Use for:** Taxonomic validation, not sequences. Verify accepted species names before searching GenBank.

```bash
# Check accepted name for a taxon via GBIF API
curl "https://api.gbif.org/v1/species?name=Zingiber+officinale" | python3 -m json.tool
```

---

## 6. Choosing Between Databases for Common Tasks

| Task | Primary source | Supplement |
|------|---------------|-----------|
| Standard marker phylogeny (rbcL, matK, ITS) | GenBank nuccore | BOLD for barcoding loci |
| Whole plastome | GenBank (chloroplast filter) | SRA (genome skimming) |
| Angiosperms353 / HybSeq | SRA (Targeted-Capture) | — |
| Reuse published alignment | TreeBASE | GenBank for additions |
| Taxonomic name validation | GBIF | NCBI taxonomy |

---

## 7. Checklist Before Downloading

- [ ] Survey all databases before downloading anything (see `data-acquisition` SKILL.md Step 2)
- [ ] Confirm library strategy for SRA records (WGS vs. RNA-seq vs. amplicon)
- [ ] Estimate total SRA dataset size with `vdb-dump --info` or `esearch runinfo`
- [ ] Check storage availability: `df -h .`
- [ ] Verify naming convention will be: `Genus_species_accession_marker.fasta`
- [ ] Record Entrez Direct and SRA Toolkit versions in report
