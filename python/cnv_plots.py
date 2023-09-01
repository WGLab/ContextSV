"""
cnv_plots.py: Plot the copy number variants and their LRR, BAF values.
"""

import os
import sys
import plotly
import pandas as pd

def run(vcf_path, cnv_data_path, output_path, region):
    # Open the VCf file and get the CNV data.
    vcf_file = open(vcf_path, "r")
    cnv_data = pd.read_csv(cnv_data_path, sep="\t")

    # Get the chromosome and start and end positions.
    chromosome, pos = region.split(":")
    start, end = pos.split("-")
    start = int(start)
    end = int(end)

    # Get the CNV data for the region.
    cnv_data = cnv_data[(cnv_data["chromosome"] == chromosome) & (cnv_data["start"] >= start) & (cnv_data["end"] <= end)]
    