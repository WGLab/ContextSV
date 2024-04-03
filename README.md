
<img src="https://github.com/WGLab/ContextSV/assets/14855676/2749b420-90ac-46ef-9311-0dbd0d086e72" width=15% height=15%>

# ContextSV

An alignment-based, generalized structural variant caller for long-read
sequencing/mapping data 

[![build
tests](https://github.com/WGLab/ContextSV/actions/workflows/build-tests.yml/badge.svg)](https://github.com/WGLab/ContextSV/actions/workflows/build-tests.yml)

Class documentation available at [https://wglab.openbioinformatics.org/ContextSV](https://wglab.openbioinformatics.org/ContextSV)

contextSV takes as input a long read alignments file (BAM), a 
reference genome file (FASTA), a VCF file with high-quality SNPs 
 (e.g. via GATK, Deepvariant, Nanocaller), and (optionally) database
 VCF files for each chromosomes with SNP population frequencies (specifics
 provided [below](link)).

## Software requirements
Please refer to the conda
[environment.yml](link)
file for all required packages.

# Installation using Anaconda (Linux)
First, install [Anaconda](https://www.anaconda.com/).

Next, create a new environment. This installation has been tested with Python 3.11:

```
conda create -n contextsv
conda activate contextsv
```

contextSV can then be installed using the following command:

```
conda install -c bioconda -c wglab contextsv=1.0.0
```

# Building from source
First install [Anaconda](https://www.anaconda.com/). Then follow the instructions below to install LongReadSum and its dependencies:

```
git clone https://github.com/WGLab/ContextSV
cd ContextSV
conda env create -f environment.yml
make

```



## Running
Activate the conda environment and then run with arguments:
```
conda activate contextsv
python contextsv [arguments]
```

# General Usage

Specifying input files:

```
ContextSV: A tool for integrative structural variant detection.

options:
  -h, --help            show this help message and exit
  -r REGION, --region REGION
                        region to analyze (e.g. chr1, chr1:1000-2000). If not provided, the entire genome will be analyzed
  -o OUTPUT, --output OUTPUT
                        path to the output directory
  -v, --version         print the version number and exit
  -d, --debug           debug mode (verbose logging)
  -t THREADS, --threads THREADS
                        number of threads to use
  -sr SHORT_READ, --short-read SHORT_READ
                        path to the short read alignment BAM file
  -lr LONG_READ, --long-read LONG_READ
                        path to the long read alignment BAM file
  -g REFERENCE, --reference REFERENCE
                        path to the reference genome FASTA file
  -s SNPS, --snps SNPS  path to the SNPs VCF file
  --pfb PFB             path to the text file listing chromosome SNP population allele frequency VCF filepaths (see docs for format)
  --hmm HMM             path to the PennCNV HMM file
```

# Revision history
For release history, please visit [here](https://github.com/WGLab/LongReadSum/releases). 

# Getting help
Please refer to the [LongReadSum issue pages](https://github.com/WGLab/LongReadSum/issues) for posting your issues. We will also respond your questions quickly. Your comments are criticl to improve our tool and will benefit other users.
