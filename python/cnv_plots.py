"""
cnv_plots.py: Plot the copy number variants and their log2_ratio, BAF values.
"""

import os
import sys
import plotly
import pandas as pd

def run(vcf_path, cnv_data_path, output_path):
    # Get the CNV data.
    vcf_file = open(vcf_path, "r")
    cnv_data = pd.read_csv(cnv_data_path, sep="\t")

    # Set up the CNV type string list for the 6 CNV states.
    cnv_types = ["NAN", "DEL", "DEL", "NEUT", "NEUT", "DUP", "DUP"]

    # Get the log2_ratio and BAF values.
    l2r_values = cnv_data["log2_ratio"]
    baf_values = cnv_data["b_allele_freq"]

    # # Get the CNV types as a list of integers (non-64 bit integers) using numpy.
    # cnv_states = cnv_data["cnv_state"].to_numpy()
    # cnv_states = cnv_states.astype(int)

    # Get the CNV colors.
    cnv_colors = []
    cnv_states = cnv_data["cnv_state"]
    for cnv_state in cnv_states:
        if cnv_state in [1, 2]:  # DEL
            cnv_colors.append("red")
        elif cnv_state in [3, 4]:  # NEUT
            cnv_colors.append("black")
        elif cnv_state in [5, 6]:  # DUP
            cnv_colors.append("blue")

    # Get the CNV names.
    cnv_names = []
    for index, row in cnv_data.iterrows():
        cnv_state = row["cnv_state"]
        cnv_names.append(f"{cnv_types[cnv_state]} {row['chromosome']}:{row['position']}")

    # Get the CNV hover text.
    cnv_hover_text = []
    for index, row in cnv_data.iterrows():
        cnv_hover_text.append(f"TYPE: {cnv_types[row['cnv_state']]}<br>CHR: {row['chromosome']}<br>POS: {row['position']}<br>L2R: {row['log2_ratio']}<br>BAF: {row['b_allele_freq']}")

    # Plot the CNVs.
    cnv_trace = plotly.graph_objs.Scatter(
        x = cnv_data["position"],
        y = cnv_data["log2_ratio"],
        mode = "markers+lines",
        name = "CNVs",
        text = cnv_hover_text,
        hoverinfo = "text",
        marker = dict(
            color = cnv_colors,
            size = 10,
        )
    )

    # Read the VCF file.
    vcf_lines = vcf_file.readlines()

    # Get all CNVs from the VCF file (SVTYPE = DEL, DUP) and store the start and
    # end positions.
    cnv_start_positions = []
    cnv_end_positions = []
    for line in vcf_lines:
        if line.startswith("#"):
            continue
        line = line.strip()
        line_parts = line.split("\t")
        if line_parts[4] in ["DEL", "DUP"]:
            cnv_start_positions.append(int(line_parts[1]))
            cnv_end_positions.append(int(line_parts[7].split(";")[0].split("=")[1]))
        

    # Define the layout.
    layout = plotly.graph_objs.Layout(
        title = f"CNVs for region {cnv_data['chromosome'][0]}:{cnv_data['position'][0]}-{cnv_data['position'][len(cnv_data) - 1]}",
        xaxis = dict(
            title = "Chromosome Position"
        ),
        yaxis = dict(
            title = r"Log<sub>2</sub> Ratio "
        )
    )

    # Create the figure.
    fig = plotly.graph_objs.Figure(
        data = [cnv_trace],
        layout = layout
    )

    # Save the figure.
    filepath = os.path.join(output_path, "cnv_plot.html")
    plotly.offline.plot(fig, filename=filepath, auto_open=False)

    # Close the VCF file.
    vcf_file.close()

    print("Saved CNV plot to {}".format(filepath))
