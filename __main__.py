"""
__main__.py: Run the program.
"""

import os
import argparse
from lib import contextsv

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


if __name__ == '__main__':

    # Run the program.
    main()
