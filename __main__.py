"""
__main__.py: Run the program.
"""

import os
import argparse
from lib import contextsv
from python import cnv_plots

def main():

    # Grab the command line arguments using argparse.
    parser = argparse.ArgumentParser(
        description="ContextSV: A tool for integrative structural variant detection."
    )
    parser.add_argument(
        "-b", "--bam",
        help="The path to the BAM file.",
        required=True
    )
    parser.add_argument(
        "-g", "--reference",
        help="The path to the reference genome.",
        required=True
    )
    parser.add_argument(
        "-s", "--snps",
        help="The path to the SNPs file.",
        required=True
    )
    parser.add_argument(
        "-o", "--output",
        help="The path to the output file.",
        required=True
    )
    parser.add_argument(
        "-r", "--region",
        help="The region to analyze.",
        required=True
    )
    
    # Get the command line arguments.
    args = parser.parse_args()

    # Run the program.
    print("Running contextsv with the following arguments:")
    print("BAM: {}".format(args.bam))
    print("Reference: {}".format(args.reference))
    print("SNPs: {}".format(args.snps))
    print("Output: {}".format(args.output))
    print("Region: {}".format(args.region))
    
    contextsv.run(
        args.bam,
        args.reference,
        args.snps,
        args.output,
        args.region
    )

    # Run the python-based analysis.
    print("Running python-based analysis.")
    vcf_path = os.path.join(args.output, "sv_calls.vcf")
    cnv_data_path = os.path.join(args.output, "cnv_data.tsv")
    region = args.region
    output_dir = args.output
    print("VCF: {}".format(vcf_path))
    print("CNV Data: {}".format(cnv_data_path))
    print("Region: {}".format(region))
    
    cnv_plots.run(vcf_path, cnv_data_path, output_dir, region)
    #cnv_plots.run(vcf_path, cnv_data_path, args.output, args.region)


if __name__ == '__main__':

    # Run the program.
    main()
