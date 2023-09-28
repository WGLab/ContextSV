"""
__main__.py: Run the program.
"""

import os
import sys
import argparse
import logging as log
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

def main():

    # Grab the command line arguments using argparse.
    parser = argparse.ArgumentParser(
        description="ContextSV: A tool for integrative structural variant detection."
    )

    # Add common arguments.
    parser.add_argument(
        "-r", "--region",
        help="The region to analyze.",
        required=True
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
        version="%(prog)s 0.0.1"
    )

    # Two modes: SV detection and CNV plots. Create a subparser for each mode.
    subparsers = parser.add_subparsers(
        title="Program mode",
        description="The program mode to run.",
        dest="mode"
    )

    # Mode 1: SV detection mode.
    sv_parser = subparsers.add_parser(
        "sv",
        help="Run SV detection."
    )
    sv_parser.add_argument(
        "-b", "--bam",
        help="The path to the BAM file.",
        required=True
    )
    sv_parser.add_argument(
        "-g", "--reference",
        help="The path to the reference genome.",
        required=True
    )
    sv_parser.add_argument(
        "-s", "--snps",
        help="The path to the SNPs file.",
        required=True
    )

    # PFB file of population allele frequencies.
    sv_parser.add_argument(
        "-p", "--pfb",
        help="The path to the PFB file of population allele frequencies.",
        required=False
    )

    # Pass in chromosome mean coverage data for speediness.
    sv_parser.add_argument(
        "-c", "--chr-cov",
        help="Chromosome mean coverage values passed in as a comma-separated list (e.g. chr1:100,chr2:200,chr3:300).",
        required=False
    )

    # Mode 2: CNV plots mode.
    cnv_parser = subparsers.add_parser(
        "plot_cnv",
        help="Run CNV plots."
    )
    cnv_parser.add_argument(
        "-c", "--cnv",
        help="The path to the CNV data in TSV format.",
        required=False
    )
    cnv_parser.add_argument(
        "-v", "--vcf",
        help="The path to the VCF file.",
        required=False
    )
    
    # Get the command line arguments.
    args = parser.parse_args()

    # Remove commas or spaces from the region.
    region = args.region
    region = region.replace(",", "")
    region = region.replace(" ", "")

    # Determine the selected program mode (SV detection or CNV plots).
    program_mode = args.mode
    mode_str = "SV detection" if program_mode == "sv" else "CNV plots"
    print("Program mode: {}".format(mode_str))

    if (program_mode == "sv"):
        
        # Run SV detection.
        log.info("Running SV detection.")
        log.info("BAM filepath: %s", args.bam)
        log.info("Reference filepath: %s", args.reference)
        log.info("SNPs filepath: %s", args.snps)
        log.info("Output directory: %s", args.output)
        log.info("Region: %s", args.region)
        
        if args.pfb is None:
            args.pfb = ""
        if args.chr_cov is None:
            args.chr_cov = ""
            
        log.info("Chromosome mean coverage: %s", args.chr_cov)
        log.info("PFB filepath: %s", args.pfb)
    
        contextsv.run(
            args.bam,
            args.reference,
            args.snps,
            args.output,
            args.region,
            args.chr_cov,
            args.pfb
        )

        vcf_path = os.path.join(args.output, "sv_calls.vcf")
        cnv_data_path = os.path.join(args.output, "cnv_data.tsv")
        output_dir = args.output
    else:
        vcf_path = args.vcf
        cnv_data_path = args.cnv
        output_dir = args.output

    # Run the python-based analysis.
    log.info("Running python-based analysis...")
    cnv_plots.run(vcf_path, cnv_data_path, output_dir, region)
    log.info("Done.")


if __name__ == '__main__':

    # Run the program.
    main()
