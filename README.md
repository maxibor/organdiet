# OrganDiet

<img src="./img/logo.png" width="200">

**Currently in development. Until now, can only run on a Linux based machine**

## Introduction

**OrganDiet** is a [Nextflow](https://www.nextflow.io/) pipeline to infer a human diet based on [shotgun metagenomics](https://en.wikipedia.org/wiki/Metagenomics#Shotgun_metagenomics) data.

# Dependancies

- [Conda](https://conda.io/miniconda.html)  
- [BASTA](https://github.com/timkahlke/BASTA) with installed `prot` database. 

# Installation

```
git clone https://github.com/maxibor/organdiet.git
cd organdiet
conda env create -f environment.yml
source activate organdiet
```

# Usage

`nextflow run organdiet.nf`
