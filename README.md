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

### Building from source (for testing/development)
ContextSV requires HTSLib as a dependency that can be installed using  [Anaconda](https://www.anaconda.com/). Create an environment
containing HTSLib: 

```
conda create -n htsenv -c bioconda -c conda-forge htslib
conda activate htsenv
```

Then follow the instructions below to build ContextSV:

```
git clone https://github.com/WGLab/ContextSV
cd ContextSV
make
```

ContextSV can then be run:
```
./build/contextsv --help

Usage: ./build/contextsv [options]
Options:
  -b, --bam <bam_file>          Long-read BAM file (required)
  -r, --ref <ref_file>          Reference genome FASTA file (required)
  -s, --snp <vcf_file>          SNPs VCF file (required)
  -o, --outdir <output_dir>     Output directory (required)
  -c, --chr <chromosome>        Chromosome
  -r, --region <region>         Region (start-end)
  -t, --threads <thread_count>  Number of threads
  -h, --hmm <hmm_file>          HMM file
  -n, --sample-size <size>      Sample size for HMM predictions
     --min-cnv <min_length>     Minimum CNV length
     --eps <epsilon>             DBSCAN epsilon
     --min-pts-pct <min_pts_pct> Percentage of mean chr. coverage to use for DBSCAN minimum points
  -e, --eth <eth_file>          ETH file
  -p, --pfb <pfb_file>          PFB file
     --save-cnv                 Save CNV data
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

## Revision history
For release history, please visit [here](https://github.com/WGLab/ContextSV/releases). 

## Getting help
Please refer to the [contextSV issue pages](https://github.com/WGLab/ContextSV/issues) for posting your issues. We will also respond your questions quickly. Your comments are critical to improve our tool and will benefit other users.
