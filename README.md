[![build
tests](https://github.com/WGLab/ContextSV/actions/workflows/build-tests.yml/badge.svg)](https://github.com/WGLab/ContextSV/actions/workflows/build-tests.yml)

# ContextSV

<p>
<img src="https://github.com/user-attachments/assets/03603ad1-df9d-438d-911c-81af0cf612e3" alt="ContextSV" align="left" style="width:100px;"/>
A long-read, whole-genome structural variant (SV) caller. It takes as input long read alignments (BAM), the 
corresponding reference genome (FASTA), a VCF with high-quality SNPs 
 (e.g. via GATK, Deepvariant, <a href="https://github.com/WGLab/NanoCaller">NanoCaller</a>, and <a href="https://gnomad.broadinstitute.org/downloads">gnomAD</a> database
 VCF files with SNP population frequencies for each chromosome.
Class documentation is available at <a href="https://wglab.openbioinformatics.org/ContextSV">https://wglab.openbioinformatics.org/ContextSV</a>
</p>

## Installation
### Anaconda
First, install [Anaconda](https://www.anaconda.com/).

Next, create a new environment. This installation has been tested with Python 3.10, Linux 64-bit.

```bash
conda create -n contextsv python=3.10
conda activate contextsv
```

ContextSV and its dependencies can then be installed using the following command:

```bash
conda install -c wglab -c conda-forge -c bioconda contextsv

# Or using mamba (faster dependency resolution):
mamba install -c wglab contextsv
```

After installation, you should have access to the following commands in your terminal:

- `contextsv`: the main SV caller
- `contextsv-cnv-plot`: utility to generate CNV plots from ContextSV JSON output
- `contextscore`: [ContextScore](https://github.com/WGLab/ContextScore) utility for post-filtering of low-confidence SV calls

Example usage:

```bash
# SV calling example:
contextsv \
  --bam sample.bam \
  --ref hg38.fa \
  --outdir output/ \
  --threads 4 \
  --snp snps.vcf \
  --eth nfe \
  --pfb gnomadv4_filepaths.txt \
  --assembly-gaps hg38-gaps.bed \   # optional: assembly gaps file
  --save-cnv                        # optional: save CNV calls in JSON

# SV post-filtering example:
contextscore \
  --input input.vcf \
  --output scored.vcf \
  --sample-coverage 30 \
  --buildver hg38 \
  --threshold 0.2 \
  --annovar /path/to/annovar \
  --annovar-db /path/to/humandb


# CNV plotting example:
contextsv-cnv-plot ./output/CNVCalls.json chr3 --formats html,svg --output-dir ./CNV_Plots
```

### Docker
First, install [Docker](https://docs.docker.com/engine/install/).
Pull the latest image from Docker hub, which contains the latest release and its dependencies.

```bash
docker pull genomicslab/contextsv
```

Example usage:

```bash
# SV calling:
docker run --rm genomicslab/contextsv --help

# SV post-filtering:
docker run --rm \
  -v /path/to/data:/mnt \
  genomicslab/contextsv \
  contextscore \
  --help

# CNV plotting:
docker run --rm \
  -v /path/to/data:/mnt \
  genomicslab/contextsv \
  contextsv-cnv-plot \
  --help
```


## Building from source (for testing/development)
ContextSV requires HTSLib as a dependency that can be installed using  [Anaconda](https://www.anaconda.com/). Create an environment
containing HTSLib: 

```bash
conda create -n htsenv -c bioconda -c conda-forge htslib
conda activate htsenv
```

Then follow the instructions below to build ContextSV:

```bash
git clone https://github.com/WGLab/ContextSV
cd ContextSV
make
```

ContextSV can then be run:
```bash
./build/contextsv --help

Options:
  -b, --bam <bam_file>          Long-read BAM file (required)
  -r, --ref <ref_file>          Reference genome FASTA file (required)
  -s, --snp <vcf_file>          Long-read SNP VCF file (required)
  -o, --outdir <output_dir>     Output directory (required)
  -t, --threads <thread_count>  Number of threads, chromosome-level parallelization (default: 1)
  -h, --hmm <hmm_file>          HMM parameter file for copy number predictions (included in the repository)
  -e, --eth <eth_file>          Ethnicity as used in gnomAD (e.g. "asj" for Ashkenazi Jewish, "nfe" for Non-Finnish European, etc.)
  -p, --pfb <pfb_file>          File containing per-chromosome population allele frequency filepaths as described in this documentation
     --assembly-gaps <gaps_file> Assembly gaps file in BED format available from UCSC Genome Browser (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz for GRCh38)
     --save-cnv                 Save CNV data in JSON for downstream plotting with contextsv-cnv-plot
     --debug                    Debug mode with verbose logging
     --version                  Print version and exit
  -h, --help                    Print usage and exit
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


### Script for downloading gnomAD VCFs
```bash
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
```bash
1=~/data/gnomad/v4.0.0/gnomad.genomes.v4.0.sites.chr1.vcf.bgz
2=~/data/gnomad/v4.0.0/gnomad.genomes.v4.0.sites.chr2.vcf.bgz
3=~/data/gnomad/v4.0.0/gnomad.genomes.v4.0.sites.chr3.vcf.bgz
...
X=~/data/gnomad/v4.0.0/gnomad.genomes.v4.0.sites.chrX.vcf.bgz
Y=~/data/gnomad/v4.0.0/gnomad.genomes.v4.0.sites.chrY.vcf.bgz
```

## Revision history
For release history, please visit [here](https://github.com/WGLab/ContextSV/releases). 

## Getting help
Please refer to the [contextSV issue pages](https://github.com/WGLab/ContextSV/issues) for posting your issues. We will also respond your questions quickly. Your comments are critical to improve our tool and will benefit other users.
