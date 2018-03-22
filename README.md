<img src="./img/logo.png" width="150">

## Introduction

**OrganDiet** is a [Nextflow](https://www.nextflow.io/) pipeline to infer a human diet based on [shotgun metagenomics](https://en.wikipedia.org/wiki/Metagenomics#Shotgun_metagenomics) data.

**Currently in development. For now, you can only run it on a Linux based machine**

# Dependancies

-   [Conda](https://conda.io/miniconda.html)  

# Quick start

Assuming you already have all databases and the conda environment installed

    conda activate organdiet
    nextflow run maxibor/organdiet --reads '*_R{1,2}.fastq.gz' -with-report run_report.html -with-dag flowchart.png

# Installation

### 1. Set up conda environments

    wget https://github.com/maxibor/organdiet/archive/v0.2.2.zip
    unzip v0.2.2.zip
    cd organdiet-0.2.2
    conda env create -f envs/organdiet.yml
    source activate organdiet

### 2. Set up Taxonomy database

-   Install taxonomy database: `./bin/basta taxonomy -o ./taxonomy`

### 3. Download the Bowtie2 index for the host genome

From [illumina **iGenomes**](https://support.illumina.com/sequencing/sequencing_software/igenome.html)

    mkdir hs_genome
    cd hs_genome
    wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/Ensembl/GRCh37/Homo_sapiens_Ensembl_GRCh37.tar.gz
    tar -xvzf Homo_sapiens_Ensembl_GRCh37.tar.gz
    cd ..

### 4. Download the organellome database and build Bowtie2 index

From NCBI [Refseq organelles genomes](https://www.ncbi.nlm.nih.gov/genome/organelle/)

    ./bin/download_organellome_db.sh
    bowtie2-build organellome_db/organellome.fa organellome_db/organellome

### 4. nt/nr database set up: two solutions

#### Case 1: You plan on using the nr database

##### 4.1.1 Set up `nr` database for [Diamond](https://github.com/bbuchfink/diamond)

    mkdir nr_diamond_db
    cd nr_diamond_db
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
    gunzip nr.gz
    mv nr nr.fa
    diamond makedb --in nr.fa -d nr
    cd ..

#### 4.1.2 Set up TaxID mapping database

-   Install _prot_ database:   


    ./bin/basta download prot -d ./taxonomy

#### Case 2: You plan on using the nt database

##### 4.1.1 Download and extract the centrifuge database

    mkdir nt_db
    cd nt_db
    wget http://som1.ific.uv.es/nt/nt.cf.7z
    7z e nt.cf.7z
    cd ..

#### 4.1.2 Set up krona mapping database

-   Install _krona_ database:  


    ktUpdateTaxonomy.sh ./taxonomy

# Get help

    nextflow run maxibor/organdiet --help

# An example workflow for this pipeline

![](./img/flowchart.png)

# Credits

The OrganDiet pipeline uses many tools listed below:

-   [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
-   [AdapterRemoval](https://github.com/MikkelSchubert/adapterremoval)
-   [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
-   [Samtools](http://www.htslib.org/)
-   [Diamond](https://github.com/bbuchfink/diamond)
-   [Centrifuge](https://ccb.jhu.edu/software/centrifuge/manual.shtml)
-   [BASTA](https://github.com/timkahlke/BASTA)
-   [Krona](https://github.com/marbl/Krona/wiki)
-   [Nextflow](https://www.nextflow.io/)

The author of OrganDiet also got some inspiration and help from the following awesome developers:

-   [Paolo Di Tomasso](https://twitter.com/paoloditommaso)
-   [Phil Ewels](https://twitter.com/tallphil)
-   [Tim Kahlke](https://twitter.com/AdvancedTwigTec)
