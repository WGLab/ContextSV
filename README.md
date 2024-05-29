
[![build
tests](https://github.com/WGLab/ContextSV/actions/workflows/build-tests.yml/badge.svg)](https://github.com/WGLab/ContextSV/actions/workflows/build-tests.yml)

![contextsv_small_15p](https://github.com/WGLab/ContextSV/assets/14855676/79d70c76-a34a-472e-a14c-e49489ae0f09)

# ContextSV
_This is a work in progress, software is under development and not ready for official release._

An alignment-based, generalized structural variant caller for long-read
sequencing/mapping data.

ContextSV takes as input a long read alignments file (BAM), a 
corresponding reference genome file (FASTA), a VCF file with high-quality SNPs 
 (e.g. via GATK, Deepvariant, [NanoCaller](https://github.com/WGLab/NanoCaller)), and [gnomAD](https://gnomad.broadinstitute.org/downloads) database
 VCF files with SNP population frequencies for each chromosome.

Class documentation is available at [https://wglab.openbioinformatics.org/ContextSV](https://wglab.openbioinformatics.org/ContextSV)

## Installation (Linux)
### Using Anaconda (recommended)
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

## Downloading gnomAD SNP population frequencies
SNP population allele frequency
information is used for copy number predictions in this tool (see
[PennCNV](http://www.genome.org/cgi/doi/10.1101/gr.6861907) for specifics). We
recommend downloading this data from the Genome Aggregation Database (gnomAD).

Download links for genome VCF files are located here (last updated April 3,
2024):

 - **gnomAD v4.0.0 (GRCh38)**: https://gnomad.broadinstitute.org/downloads#4

 - **gnomAD v2.1.1 (GRCh37)**: https://gnomad.broadinstitute.org/downloads#2


### Example download
```
download_dir="~/data/gnomad/v4.0.0/"

chr_list=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")

for chr in "${chr_list[@]}"; do
    echo "Downloading chromosome ${chr}..."
    wget "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.chr${chr}.vcf.bgz" -P "${download_dir}"
done
```

Finally, create a text file that specifies the chromosome and its corresponding
gnomAD filepath. This file will be passed in as an argument:

**gnomadv4_filepaths.txt**
```
1=~/data/gnomad/v4.0.0/gnomad.genomes.v4.0.sites.chr1.vcf.bgz
2=~/data/gnomad/v4.0.0/gnomad.genomes.v4.0.sites.chr2.vcf.bgz
3=~/data/gnomad/v4.0.0/gnomad.genomes.v4.0.sites.chr3.vcf.bgz
...
X=~/data/gnomad/v4.0.0/gnomad.genomes.v4.0.sites.chrX.vcf.bgz
Y=~/data/gnomad/v4.0.0/gnomad.genomes.v4.0.sites.chrY.vcf.bgz
```

## Calling structural variants
### Example full script generating a merged VCF of structural variants
```
# Activate the environment
conda activate contextsv

# Set the input reference genome
ref_file="~/data/GRCh38.fa"

# Set the input alignment file (e.g. from minimap2)
long_read_bam="~/data/HG002.GRCh38.bam"

# Set the input SNPs file (e.g. from NanoCaller)
snps_file="~/data/variant_calls.snps.vcf.gz"

# Set the SNP population frequencies filepath
pfb_file="~/data/gnomadv4_filepaths.txt"

# Set the output directory
output_dir=~/data/contextSV_output

# Specify the number of threads (system-specific)
thread_count=40

# Run SV calling (~3-4 hours for whole-genome, 40 cores)
python contextsv --threads $thread_count -o $output_dir -lr $long_read_bam --snps $snps_file --reference $ref_file --pfb $pfb_file

# The output VCF filepath is located here:
output_vcf=$output_dir/output.vcf

# Merge SVs (~3-4 hours for whole-genome, 40 cores)
python contextsv --merge $output_vcf

# The final merged VCF filepath is located here:
merged_vcf=$output_dir/output.merged.vcf
```

## Input arguments

```
python contextsv --help

ContextSV: A tool for integrative structural variant detection.

options:
  -h, --help            show this help message and exit
  -lr LONG_READ, --long-read LONG_READ
                        path to the long read alignment BAM file
  -g REFERENCE, --reference REFERENCE
                        path to the reference genome FASTA file
  -s SNPS, --snps SNPS  path to the SNPs VCF file
  --pfb PFB             path to the file with SNP population frequency VCF filepaths (see docs for format)
  -o OUTPUT, --output OUTPUT
                        path to the output directory
  -r REGION, --region REGION
                        region to analyze (e.g. chr1, chr1:1000-2000). If not provided, the entire genome will be analyzed
  -t THREADS, --threads THREADS
                        number of threads to use
  --hmm HMM             path to the PennCNV HMM file
  --window-size WINDOW_SIZE
                        window size for calculating log2 ratios for CNV predictions (default: 10 kb)
  -d, --debug           debug mode (verbose logging)
  -v, --version         print the version number and exit
```

## Revision history
For release history, please visit [here](https://github.com/WGLab/ContextSV/releases). 

## Getting help
Please refer to the [contextSV issue pages](https://github.com/WGLab/ContextSV/issues) for posting your issues. We will also respond your questions quickly. Your comments are critical to improve our tool and will benefit other users.
