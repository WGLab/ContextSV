
<img src="https://github.com/WGLab/ContextSV/assets/14855676/2749b420-90ac-46ef-9311-0dbd0d086e72" width=15% height=15%>

# ContextSV

An alignment-based, generalized structural variant caller for long-read
sequencing/mapping data 

[![build
tests](https://github.com/WGLab/ContextSV/actions/workflows/build-tests.yml/badge.svg)](https://github.com/WGLab/ContextSV/actions/workflows/build-tests.yml)

Class documentation available at [https://wglab.openbioinformatics.org/ContextSV](https://wglab.openbioinformatics.org/ContextSV)

ContextSV takes as input a long read alignments file (BAM), a 
reference genome file (FASTA), a VCF file with high-quality SNPs 
 (e.g. via GATK, Deepvariant, Nanocaller), and (optionally) [gnomAD](https://gnomad.broadinstitute.org/downloads) database
 VCF files with SNP population frequencies for each chromosome.

## Installation (Linux)
### Using Anaconda (tecommended)
First, install [Anaconda](https://www.anaconda.com/).

Next, create a new environment. This installation has been tested with Python 3.11:

```
conda create -n contextsv python=3.11
conda activate contextsv
```

ContextSV can then be installed using the following command:

```
conda install -c bioconda -c wglab contextsv=1.0.0
```

### Building from source (for testing/development)
First install [Anaconda](https://www.anaconda.com/). Then follow the instructions below to install LongReadSum and its dependencies:

```
git clone https://github.com/WGLab/ContextSV
cd ContextSV
conda env create -f environment.yml
make
```

### Downloading gnomAD data with SNP population frequencies
Population allele frequency
information is used for SNP copy number predictions in this tool (see
[PennCNV](http://www.genome.org/cgi/doi/10.1101/gr.6861907) for specifics). We
recommend downloading this data from the Genome Aggregation Database (gnomAD).

Download links for genome VCF files are located here (last updated April 3,
2024):

**gnomAD v4.0.0 (GRCh38)**: https://gnomad.broadinstitute.org/downloads#4

**gnomAD v2.1.1 (GRCh37)**: https://gnomad.broadinstitute.org/downloads#2


Example download:
```
download_dir="~/data/gnomad/v4.0.0/"

chr_list=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")

for chr in "${chr_list[@]}"; do
    echo "Downloading chromosome ${chr}..."
    
    # V4 (hg38)
    wget "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.chr${chr}.vcf.bgz" -P "${download_dir}"
    echo "Done"
done
```

Finally, create a text file that specifies the chromosome and its corresponding
gnomAD filepath:





# Running
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
