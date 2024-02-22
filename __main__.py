"""
__main__.py: Run the program.
"""

__version__ = "0.0.1"

import os
import sys
import argparse
import logging as log
from io import StringIO

from lib import contextsv
from python import cnv_plots

# Set up logging.
log.basicConfig(
    level=log.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        log.StreamHandler(sys.stdout)
    ]
)

# Define a class for redirecting c++ stdout to python logging info.
class LogRedirect(StringIO):
    """Redirect c++ stdout to python logging info."""
    def write(self, buf):
        super(LogRedirect, self).write(buf)
        log.info(buf)

# Redirect c++ stdout to python logging info.
sys.stdout = LogRedirect()

def main():
    """Entry point and user interface for the program."""

    # Grab the command line arguments using argparse.
    parser = argparse.ArgumentParser(
        description="ContextSV: A tool for integrative structural variant detection."
    )

    # Add common arguments.
    parser.add_argument(
        "-r", "--region",
        help="The region to analyze. If not provided, the entire genome will be analyzed.",
        required=False,
    )

    parser.add_argument(
        "-o", "--output",
        help="The path to the output directory.",
        required=True
    )

    parser.add_argument(
        "-v", "--version",
        help="Print the version number and exit.",
        action="version",
        version=f"contextSV version {__version__}"
    )

    # Verbose mode.
    parser.add_argument(
        "-d", "--debug",
        help="Enable debug mode (verbose logging).",
        action="store_true",
        default=False
    )

    # Thread count.
    parser.add_argument(
        "-t", "--threads",
        help="The number of threads to use.",
        required=False,
        default=1,
        type=int
    )

    # CNV data path.
    parser.add_argument(
        "--cnv",
        help="The path to the CNV data in TSV format.",
        required=False
    )

    # Mode 1: SV detection mode.
    # Short read alignment file (BAM), reference genome, and short read SNPs file.
    parser.add_argument(
        "-sr", "--short-read",
        help="The path to the short read alignment file (BAM).",
        required=False
    )
    parser.add_argument(
        "-lr", "--long-read",
        help="The path to the long read alignment file (BAM).",
        required=False
    )

    parser.add_argument(
        "-g", "--reference",
        help="The path to the reference genome.",
        required=False
    )
    parser.add_argument(
        "-s", "--snps",
        help="The path to the SNPs file.",
        required=False
    )

    # Text file with VCF filepaths of SNP population allele frequencies for each
    # chromosome from a database such as gnomAD (e.g. 1=chr1.vcf.gz\n2=chr2.vcf.gz\n...).
    parser.add_argument(
        "--pfb",
        help="Path to the text file listing VCF population allele frequency filepaths for each chromosome.",
        required=False
    )

    # HMM file path.
    parser.add_argument(
        "--hmm",
        help="The path to the HMM file.",
        required=False
    )

    # Pass in chromosome mean coverage data for speediness.
    parser.add_argument(
        "--chr-cov",
        help="Chromosome mean coverage values passed in as a comma-separated list (e.g. chr1:100,chr2:200,chr3:300).",
        required=False
    )

    # Turn off CIGAR string SV detection (split-read only)
    parser.add_argument(
        "--disable-cigar",
        help="Turn off CIGAR string SV detection.",
        required=False,
        action="store_true",
        default=False
    )

    # Turn off SNP-based CNV predictions for SV classification.
    parser.add_argument(
        "--disable-snp-cnv",
        help="Turn off SNP-based CNV predictions.",
        required=False,
        action="store_true",
        default=False
    )

    # Mode 2: CNV plots mode. If only the VCF file is provided, then the program
    # will run in this mode.
    parser.add_argument(
        "--vcf",
        help="The path to the VCF file of SV calls.",
        required=False
    )
    
    # Get the command line arguments.
    args = parser.parse_args()

    # Determine the selected program mode (SV detection or CNV plots).
    if (args.vcf is not None):
        # Run SV analysis mode.
        
        # Ensure that the CNV data path is provided.
        arg_error = False
        if (args.cnv is None):
            log.error("Please provide the CNV data path.")
            arg_error = True

        # Exit if there are any errors.
        if (arg_error):
            sys.exit(1)

        # Create the output directory if it doesn't exist.
        if (not os.path.exists(args.output)):
            os.makedirs(args.output)

        # Set the data paths from user input for downstream analysis.
        vcf_path = args.vcf
        cnv_data_path = args.cnv
        output_dir = args.output

        log.info("Analyzing SV calls in VCF file %s...", args.vcf)

    else:
        # Run SV detection mode.

        # Ensure BAM, reference, and SNPs files are provided.
        arg_error = False
        if (args.long_read is None):
            log.error("Please provide the long read alignment file (BAM).")
            arg_error = True

        if (args.reference is None):
            log.error("Please provide the reference genome.")
            arg_error = True

        # Short read alignment file is optional. Use the long read alignment
        # file if it is not provided.
        if (args.short_read is None):
            log.warning("Short read alignment file not provided. Using long read alignment file in its place.")
            args.short_read = args.long_read

        # SNPs file is required unless SNP-based CNV predictions are disabled.
        if (args.snps is None and not args.disable_snp_cnv):
            log.error("Please provide the SNPs file.")
            arg_error = True

        # # PFB file is required unless SNP-based CNV predictions are disabled.
        # if (args.pfb is None and not args.disable_snp_cnv):
        #     log.error("Please provide the PFB file.")
        #     arg_error = True

        # Exit if there are any errors.
        if (arg_error):
            # Exit with error code 1.
            sys.exit(1)

        # Set all None values to empty strings.
        for key, value in vars(args).items():
            if value is None:
                setattr(args, key, "")
            else:
                log.info("Setting %s to %s", key, value)
        
        # Loop and set all None values to empty strings.
        for key, value in vars(args).items():
            if value is None:
                setattr(args, key, "")

        # Set input parameters.
        input_data = contextsv.InputData()
        input_data.setVerbose(args.debug)
        input_data.setShortReadBam(args.short_read)
        input_data.setLongReadBam(args.long_read)
        input_data.setRefGenome(args.reference)
        input_data.setSNPFilepath(args.snps)
        input_data.setRegion(args.region)
        input_data.setThreadCount(args.threads)
        input_data.setMeanChromosomeCoverage(args.chr_cov)
        input_data.setAlleleFreqFilepaths(args.pfb)
        input_data.setHMMFilepath(args.hmm)
        input_data.setOutputDir(args.output)
        input_data.setDisableCIGAR(args.disable_cigar)
        input_data.setDisableSNPCNV(args.disable_snp_cnv)
        input_data.setCNVFilepath(args.cnv)

        # Run the analysis.
        contextsv.run(input_data)

        # Determine the data paths for downstream analysis.
        vcf_path = os.path.join(args.output, "sv_calls.vcf")
        output_dir = args.output
        region = args.region

        if (args.cnv == ""):
            cnv_data_path = os.path.join(args.output, "cnv_data.tsv")
        else:
            cnv_data_path = args.cnv

    # Generate python-based CNV plots if SNP-based CNV predictions are enabled.
    if (not args.disable_snp_cnv):
        log.info("Generating CNV plots...")
        cnv_plots.run(vcf_path, cnv_data_path, output_dir, region)

    log.info("Complete. Thank you for using contextSV!")

if __name__ == '__main__':

    # Run the program.
    main()
