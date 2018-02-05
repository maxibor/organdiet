# OrganDiet

<img src="./img/logo.png" width="75">

**Currently in development. Until now, can only run on a Linux based machine**

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


# Usage

`nextflow run organdiet.nf`
