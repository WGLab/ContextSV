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
        "-lr", "--long-read",
        help="path to the long read alignment BAM file",
        required=True
    )

    parser.add_argument(
        "-g", "--reference",
        help="path to the reference genome FASTA file",
        required=False
    )

    parser.add_argument(
        "-s", "--snps",
        help="path to the SNPs VCF file",
        required=False
    )

    parser.add_argument(
        '-c', '--chr',
        help="chromosome to analyze (e.g. 1, 2, 3, ..., X, Y)",
        required=False,
        default="",
        type=str
    )

    parser.add_argument(
        "-r", "--region",
        help="region to analyze (e.g. 1:1000-2000)",
        required=False,
        default="",
        type=str
    )

    # Specify the ethnicity of the sample for obtaining population allele
    # frequencies from a database such as gnomAD. If not provided, the allele
    # frequencies will be obtained for all populations.
    parser.add_argument(
        "-e", "--ethnicity",
        help="ethnicity of the sample (e.g. afr, amr, eas, fin, nfe, oth, sas, asj)",
        required=False
    )

    # Text file with VCF filepaths of SNP population allele frequencies for each
    # chromosome from a database such as gnomAD (e.g. 1=chr1.vcf.gz\n2=chr2.vcf.gz\n...).
    parser.add_argument(
        "--pfb",
        help="path to the file with SNP population frequency VCF filepaths (see docs for format)",
        required=False
    )

    parser.add_argument(
        "-o", "--output",
        help="path to the output directory",
        required=True
    )

    # Thread count.
    parser.add_argument(
        "-t", "--threads",
        help="number of threads to use",
        required=False,
        default=1,
        type=int
    )

    # HMM file path.
    parser.add_argument(
        "--hmm",
        help="path to the PennCNV HMM file",
        required=False
    )

    # Window size for calculating log2 ratios for CNV predictions.
    parser.add_argument(
        "--window-size",
        help="window size for calculating log2 ratios for CNV predictions (default: 2500 bp)",
        required=False,
        type=int,
        default=2500
    )

    # Minimum SV length for copy number variation (CNV) predictions.
    parser.add_argument(
        "--min-cnv-length",
        help="minimum SV length for CNV predictions (default: 1000 bp)",
        required=False,
        type=int,
        default=1000
    )

    # Verbose mode.
    parser.add_argument(
        "-d", "--debug",
        help="debug mode (verbose logging)",
        action="store_true",
        default=False
    )

    parser.add_argument(
        "-v", "--version",
        help="print the version number and exit",
        action="version",
        version=f"contextSV version {__version__}"
    )

    # Extend SNP-based CNV predictions to regions surrounding SVs (+/- 1/2 SV
    # length) and save CNV data to TSV. This will be useful for plotting CNV
    # data around SVs, but takes longer to run.
    parser.add_argument(
        "--save-cnv",
        required=False,
        action="store_true",
        default=False,
        help=argparse.SUPPRESS
    )

    # Mode 1: SV detection mode.
    # Short read alignment file (BAM), reference genome, and short read SNPs file.
    parser.add_argument(
        "-sr", "--short-read",
        required=False,
        help=argparse.SUPPRESS
    )

    # Chromosome mean coverage values passed in as a comma-separated list (e.g. chr1:100,chr2:200,chr3:300)
    parser.add_argument(
        "--chr-cov",
        required=False,
        help=argparse.SUPPRESS
    )
    
    # ----------------------------------------------------------------------- #

    # Run the program.

    # Get the command line arguments.
    args = parser.parse_args()

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

    # SNPs file is required
    if (args.snps is None):
        log.error("Please provide the SNPs file.")
        arg_error = True

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

    # Set input parameters
    input_data = contextsv.InputData()
    input_data.setVerbose(args.debug)
    input_data.setShortReadBam(args.short_read)
    input_data.setLongReadBam(args.long_read)
    input_data.setRefGenome(args.reference)
    input_data.setSNPFilepath(args.snps)
    input_data.setEthnicity(args.ethnicity)
    input_data.setThreadCount(args.threads)
    input_data.setChromosome(args.chr)
    input_data.setRegion(args.region)
    input_data.setMeanChromosomeCoverage(args.chr_cov)
    input_data.setAlleleFreqFilepaths(args.pfb)
    input_data.setHMMFilepath(args.hmm)
    input_data.setOutputDir(args.output)
    input_data.saveCNVData(args.save_cnv)
    input_data.setWindowSize(args.window_size)
    input_data.setMinCNVLength(args.min_cnv_length)

    # Run the analysis
    contextsv.run(input_data)

    # Determine the data paths for downstream analysis.
    vcf_path = os.path.join(args.output, "output.vcf")
    output_dir = args.output
    # cnv_data_path = os.path.join(args.output, "cnv_data.tsv")

    # Generate python-based CNV plots if SNP-based CNV predictions are enabled
    if (args.save_cnv):
        log.info("Generating CNV plots...")

        # Find all TSV files in the output directory
        for file in os.listdir(output_dir):
            if file.endswith(".tsv"):
                cnv_data_path = os.path.join(output_dir, file)

                # Set the HTML output path by changing the file extension
                output_html = os.path.splitext(cnv_data_path)[0] + ".html"

                # Generate the CNV plot for the current TSV file
                cnv_plots.run(cnv_data_path, output_html)
        # cnv_plots.run(vcf_path, cnv_data_path, output_dir, region)

    log.info("Complete. File saved to %s\nThank you for using ContextSV!", vcf_path)

if __name__ == '__main__':

    # Check if the user specified the --merge flag.
    if "--merge" in sys.argv:
        # Ensure the user provided the correct number of arguments (last 2 are
        # optional).
        if len(sys.argv) < 2:
            log.error("Usage: python __main__.py --merge <input_vcf> (optional: <epsilon> <suffix>)")
            sys.exit(1)

        # The second argument is the input VCF file.
        input_vcf = sys.argv[2]

        # The third argument is the epsilon value for the DBSCAN clustering. If
        # empty, set to 34.
        epsilon = sys.argv[3] if len(sys.argv) >= 4 else 34

        # The fourth argument is the suffix for the output file. If empty, set
        # to ".merged"
        suffix = sys.argv[4] if len(sys.argv) >= 5 else ".merged"

        # Run the SV merger.
        from python import sv_merger
        sv_merger.sv_merger(input_vcf, mode='dbscan', eps=int(epsilon), suffix=suffix)

    # Run the program.
    main()
