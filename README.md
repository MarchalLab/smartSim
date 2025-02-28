# smartSim:  Simulation of splice aware single cell Smart-seq3 data

## Prerequisites

In order to use this package: you need folowing software:
* conda (https://www.anaconda.com/)

We provide a conda environment, which can be created with:
```bash
conda env create -f smartSim_env.yml
```
## Input

* reference gtf file
* reference fasta file
* aligned and barcoded Smart-seq3 data similar to the data you want to simulate (alternatively you can use [example_data/aligned_chr1.bam](example_data/aligned_chr1.bam))

## Usage
We provide a tutorial to simulate reads using the example data in example_data: [tutorial.R](tutorial.R)

Before running the tutorial make sure you unzip [example_data/chr1.fa.gz](example_data/chr1.fa.gz) and [example_data/chr1.gtf.gz](example_data/chr1.gtf.gz)

