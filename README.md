# Small RNA-Seq Analysis Pipeline

## Overview

This repository contains a **Small RNA-Seq analysis pipeline** designed for the identification, quantification, and differential expression analysis of small RNAs, including **miRNAs, snoRNAs, snRNAs, rRNAs, and tRNAs** from Illumina sequencing data.

The pipeline is optimized for **short read lengths (18–35 nt)** and uses **Bowtie (v1)** for alignment, **featureCounts** for quantification, and **edgeR** for downstream statistical analysis.

---

## Pipeline Workflow

1. Environment setup (conda + mamba)
2. Reference genome and annotation preparation
3. Quality control and adapter trimming
4. Read length filtering
5. Genome alignment (Bowtie)
6. Read quantification (featureCounts)
7. Data integration and differential expression analysis in R

---

## Requirements

### Operating System

* Linux or WSL (Windows Subsystem for Linux)

### Conda Environment

```bash
conda create -n rnaseq_env python=3.9
conda activate rnaseq_env
conda install -n base -c conda-forge mamba
```

### Software Dependencies

Installed via **mamba**:

```bash
mamba install -c conda-forge -c bioconda \
    fastqc \
    fastq-screen \
    multiqc \
    trim-galore \
    cutadapt \
    bowtie \
    samtools \
    subread \
    sra-tools \
    adapterremoval \
    mirtrace \
    seqkit
```

---

## Reference Preparation

### Genome

Human reference genome (GRCh38, Ensembl release 113):

```bash
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

### Annotation

```bash
wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
gunzip Homo_sapiens.GRCh38.113.gtf.gz
```

Filter annotation for small RNA biotypes:

```bash
grep -P 'gene_biotype "(miRNA|snRNA|snoRNA|rRNA|tRNA)"' \
Homo_sapiens.GRCh38.113.gtf > small_RNA.gtf
```

> **Note:** Ensembl does not reliably annotate piRNAs. For piRNA analysis, use external databases such as **piRBase**.

### Bowtie Index

```bash
bowtie-build Homo_sapiens.GRCh38.dna.primary_assembly.fa bowtie_index_hg38
```

---

## Step 1: Quality Control & Adapter Trimming

### Option A: AdapterRemoval

```bash
AdapterRemoval \
  --file1 SRR10905119.fastq.gz \
  --output1 SRR10905119.trimmed.fastq.gz \
  --trimns \
  --trimqualities \
  --minlength 18
```

### Option B: Trim Galore

```bash
trim_galore --fastqc \
    --quality 20 \
    --length 18 \
    --stringency 3 \
    --cores 8 \
    --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    SRR10905119.fastq.gz
```

### Quality Assessment

```bash
fastqc SRR10905119.trimmed.fastq.gz
mirtrace qc --species hsa SRR10905119.trimmed.fastq.gz
```

---

## Step 2: Length Filtering

Filter reads ≤35 nt:

```bash
seqkit seq -M 35 SRR10905119.trimmed.fastq.gz | gzip > SRR10905119.trimmed.maxlen35.fastq.gz
```

---

## Step 3: Genome Alignment (Bowtie)

```bash
bowtie -v 1 --best --strata -l 18 \
    -x bowtie_index_hg38 \
    -p 20 \
    -q SRR10905119.trimmed.maxlen35.fastq.gz \
    -S SRR10905119.trimmed.maxlen35.sam
```

Convert and process alignments:

```bash
samtools view -bS SRR10905119.trimmed.maxlen35.sam > SRR10905119.bam
samtools sort SRR10905119.bam -o SRR10905119.sorted.bam
samtools index SRR10905119.sorted.bam
samtools flagstat SRR10905119.sorted.bam > SRR10905119.flagstat
```

---

## Step 4: Quantification

```bash
featureCounts -T 20 \
  -t miRNA \
  -g gene_id \
  -s 0 \
  -a small_RNA.gtf \
  -o SRR10905119_featureCounts_gene.txt \
  SRR10905119.sorted.bam
```

> If miRNA counts are zero, try `-t exon` or `-s 1` depending on library preparation.

---

## Step 5: Data Integration & Differential Expression (R)

### R Packages

```r
library(edgeR)
library(data.table)
library(rtracklayer)
```

### Import Counts

```r
setwd("/mnt/d/miRNA-trial.2")
files <- list.files(pattern = "_featureCounts_gene\\.txt$", full.names = TRUE)
dgl <- readDGE(files, columns = c(1,7), skip = 1)
```

### Annotation Mapping

```r
gtf <- import("small_RNA.gtf")
anno <- unique(as.data.frame(gtf)[, c("gene_id","gene_name")])
```

### Downstream Analysis

Proceed with standard **edgeR** workflow:

* Filtering low counts
* Normalization (TMM)
* Dispersion estimation
* Differential expression testing

---

## Output Files

* `*.trimmed.fastq.gz` – adapter-trimmed reads
* `*.sorted.bam` – aligned reads
* `*_featureCounts_gene.txt` – raw counts
* `counts_raw.tsv` – processed count table


---

## License

This pipeline is provided for academic and research use.
