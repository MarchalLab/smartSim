# smartSim:  Simulation of splice aware single cell Smart-seq3 data

## Prerequisites

In order to use this package: you need folowing software:
* conda (https://www.anaconda.com/)

We provide a conda environment containing all required dependencies, which can be created and activated with:
```bash
conda env create -f smartSim_env.yml
conda activate smartSim
```

> **Note:** smartSim is available for Linux (linux-64) and macOS (osx-64).  
> It is currently not compatible with Windows due to its dependency on HTSeq and Bioconda. 
> We advise Windows users to install and run smartSim within a Windows Subsystem for Linux (WSL) environment.  
> Mac users with the osx-arm64 architecture can emulate using the osx-64 subdirectory: hereto they can create the conda environment using : 
> ```bash
> CONDA_SUBDIR=osx-64 conda env create -f smartSim_env.yml
> ```

## Installation

Make sure the smartSim conda environment is active.
Next, smartSim can be installed directly from GitHub using:
```r
library(devtools)
install_github("MarchalLab/smartSim", dependencies = TRUE, repos = BiocManager::repositories())
```

Alternatively, it can be installed from source using:
```r
install.packages("path/to/smartSim, repos = NULL, type="source")
```

## Required input

* reference gtf file
* reference fasta file
* aligned and barcoded Smart-seq3 data similar to the data you want to simulate (alternatively you can use [example_data/aligned](example_data/aligned))

## Usage
We provide a tutorial to simulate reads using provided example data [here](https://github.com/pages/MarchalLab/smartSim/tutorial.html).

## Test environment
The R package was developed and tested using:

* R version : 4.4.1 (platform: x86\_64-conda-linux-gnu)
* all dependencies are provided in the corresponding conda environment: [smartSim_env.yml](smartSim_env.yml).
* test environment: All tests were run on a system with: 
  * Intel(R) Xeon(R) CPU E5-2698 v3 @ 2.30GHz
  * 256 GB RAM
  * Ubuntu 22.04.4 LTS 
* session info can be found at: [docs/sessionInfo.txt](docs/sessionInfo.txt)

