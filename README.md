# smartSim:  Simulation of splice aware single cell Smart-seq3 data

## Prerequisites

In order to use this package: you need folowing software:
* conda (https://www.anaconda.com/)

We provide a conda environment containing all required dependencies, which can be created and activated with:
```bash
conda env create -f smartSim_env.yml
conda activate smartSim
```

> **Note:** smartSim is currently only available for Linux (linux-64).
> This limitation is due to its dependency on the HTSeq package, which is not compatible with Windows or macOS.
> To ensure full functionality for users on Windows or macOS, we recommend running smartSim within a Linux-based virtual machine (VM), such as an Ubuntu VM (e.g. [WSL2](https://learn.microsoft.com/en-us/windows/wsl/about) or [multipass](https://canonical.com/multipass)).

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
* aligned and barcoded Smart-seq3 data similar to the data you want to simulate (alternatively you can use the provided example dataset as described in the [tutorial](https://marchallab.github.io/smartSim/tutorial.html))

## Usage
We provide a tutorial to simulate reads using provided example data [here](https://marchallab.github.io/smartSim/tutorial.html).

## Test environment
The R package was developed and tested using:

* R version : 4.4.1 (platform: x86\_64-conda-linux-gnu)
* all dependencies are provided in the corresponding conda environment: [smartSim_env.yml](smartSim_env.yml).
* test environment: All tests were run on a system with: 
  * Intel(R) Xeon(R) CPU E5-2698 v3 @ 2.30GHz
  * 256 GB RAM
  * Ubuntu 22.04.4 LTS 
* session info can be found at: [docs/sessionInfo.txt](docs/sessionInfo.txt)


