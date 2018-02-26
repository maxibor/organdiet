
<img src="./img/logo.png" width="150">

## Introduction

**OrganDiet** is a [Nextflow](https://www.nextflow.io/) pipeline to infer a human diet based on [shotgun metagenomics](https://en.wikipedia.org/wiki/Metagenomics#Shotgun_metagenomics) data.

**Currently in development. For now, you can only run it on a Linux based machine**

# Dependancies

- [Conda](https://conda.io/miniconda.html)  

# Installation

### 1. Set up conda environments
```
wget https://github.com/maxibor/organdiet/archive/0.2.zip
unzip 0.2.zip
cd organdiet-0.2
conda env create -f envs/organdiet.yml
conda env create -f envs/py36.yml
source activate organdiet
```

### 2. Download the Bowtie2 index for the host genome
From [illumina **iGenomes**](https://support.illumina.com/sequencing/sequencing_software/igenome.html)

```
mkdir hs_genome
cd hs_genome
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/Ensembl/GRCh37/Homo_sapiens_Ensembl_GRCh37.tar.gz
tar -xvzf Homo_sapiens_Ensembl_GRCh37.tar.gz
cd ..
```

### 3. Download the organellome database and build Bowtie2 index
From NCBI [Refseq organelles genomes](https://www.ncbi.nlm.nih.gov/genome/organelle/)
```
./bin/download_organellome_db.sh
bowtie2-build organellome_db/organellome.fa organellome_db/organellome
```


### 3. Case 1: You plan on using the nr database

#### 3.1 Set up BASTA databases
- Install taxonomy database: `./bin/basta taxonomy`
- Install *prot* database:   `./bin/basta download prot`

#### 3.2 Set up `nr` database for [Diamond](https://github.com/bbuchfink/diamond)
```
mkdir nr_diamond_db
cd nr_diamond_db
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
gunzip nr.gz
mv nr nr.fa
diamond makedb --in nr.fa -d nr
cd ..
```

### 3. Case 2: You plan on using the nt database

#### 3.1 Download and extract the centrifuge database

```
mkdir nt_db
cd nt_db
wget http://som1.ific.uv.es/nt/nt.cf.7z
7z e nt.cf.7z
cd ..
```

#### 3.2 Set up the recentrifuge taxonomy database
```
source activate py36
./bin/recentrifuge/retaxdump.py
source deactivate
```




# Run Example

```
nextflow run maxibor/organdiet --reads '*_R{1,2}.fastq.gz' -with-report run_report.html -with-dag flowchart.png
```

# Get help

```
nextflow run maxibor/organdiet --help
```

# An example workflow for this pipeline

![](./img/flowchart.png)
