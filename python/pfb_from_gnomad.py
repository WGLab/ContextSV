"""Run benchmarking on the ContextSV output VCF file."""
import os
import sys
import logging as log
import pandas as pd

from .utils import get_info_field_column, get_info_field_value


def save_snp_af_to_tsv(input_vcf, output_tsv):
    """
    Get the SNP Filtered Allele Frequency (FAF) from the input VCF file and save
    to a TSV file.
    This is the population B-allele frequency (PFB) of the SNP.
    Any allele frequency less than 0.01 (1%) is considered a very rare variant
    and is ignored in this analysis.
    (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9160216/)
    """
    # Read the input VCF file.
    vcf_df = pd.read_csv(input_vcf, sep="\t", comment="#", header=None)

    # Get the column index of the INFO field.

    # Return the SNP FAF.
    return snp_faf


# Run the program.
if __name__ == "__main__":

    # Get the GNOMAD VCF file of benchmarking variants.
    benchmark_vcf = sys.argv[1]
    
    # Get the output TSV file.
    output_tsv = sys.argv[2]

    # Run the program.
    run_cnv_benchmarking(benchmark_vcf, output_tsv)
    