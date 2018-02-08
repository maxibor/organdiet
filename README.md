# OrganDiet

<img src="./img/logo.png" width="75">

**Currently in development. For now, you can only run it on a Linux based machine**

## Introduction

**OrganDiet** is a [Nextflow](https://www.nextflow.io/) pipeline to infer a human diet based on [shotgun metagenomics](https://en.wikipedia.org/wiki/Metagenomics#Shotgun_metagenomics) data.

# Dependancies

- [Conda](https://conda.io/miniconda.html)  
- [BASTA](https://github.com/timkahlke/BASTA) and installed databases.


# Installation

```
git clone https://github.com/maxibor/organdiet.git
cd organdiet
conda env create -f environment.yml
source activate organdiet
```

### Set up BASTA databases
- Install taxonomy database: `./bin/basta taxonomy`
- Install *prot* database:   `./bin/basta download prot`

### Set up `nr` database for [Diamond](https://github.com/bbuchfink/diamond)
```
mkdir nr_diamond_db
cd nr_diamond_db
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
gunzip nr.gz
mv nr nr.fa
diamond makedb --in nr.fa -d nr
```

### Download the Bowtie2 index for the host genome
From [illumina **iGenomes**](https://support.illumina.com/sequencing/sequencing_software/igenome.html)

```
mkdir hs_genome
cd hs_genome
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/Ensembl/GRCh37/Homo_sapiens_Ensembl_GRCh37.tar.gz
tar -xvzf Homo_sapiens_Ensembl_GRCh37.tar.gz
```

### Download the organellome database and build Bowtie2 index
From NCBI [Refseq organelles genomes](https://www.ncbi.nlm.nih.gov/genome/organelle/)


# Usage

```
nextflow run organdiet.nf --reads '/path/to/reads/*_R{1,2}.fastq.gz' --ctrl '/path/to/negative/control/reads/*_R{1,2}.fastq.gz' --btindex '/path/to/bowtie2/index/organellome' --hgindex 'hs_genome/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome' --nrdb 'nr_diamond_db/nr' -with-report run_report.html -with-timeline timeline_report.html
```
