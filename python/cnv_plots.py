"""
cnv_plots.py: Plot the copy number variants and their log2_ratio, BAF values.
"""

import os
import sys
import logging as log
import plotly
import pandas as pd

# Set up logging.
log.basicConfig(
    level=log.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        log.StreamHandler(sys.stdout)
    ]
)

def get_info_field_value(info_field, field_name):
    """
    Get the value of a field in the INFO field of a VCF file.

    Args:
        info_field (str): The INFO field.
        field_name (str): The name of the field to get the value of.

    Returns:
        str: The value of the field.
    """

    # Split the INFO field into its parts.
    info_field_parts = info_field.split(";")

    # Get the field value.
    field_value = ""
    for info_field_part in info_field_parts:
        if info_field_part.startswith("{}=".format(field_name)):
            field_value = info_field_part.split("=")[1]
            break

    # Return the field value.
    return field_value

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

    # Get the CNV positions.
    cnv_types = []
    cnv_chromosomes = []
    cnv_start_positions = []
    cnv_end_positions = []
    for line in vcf_lines:
        if line.startswith("#"):
            continue
        line = line.strip()
        line_parts = line.split("\t")
        
        # Get the INFO field.
        info_field = line_parts[7]

        # Get the SVTYPE field value.
        svtype = get_info_field_value(info_field, "SVTYPE")

        # Skip if the SVTYPE is not DEL or DUP.
        if svtype != "DEL" and svtype != "DUP":
            continue

        # Get the start and end positions.
        start_position = int(line_parts[1])
        end_position = int(get_info_field_value(info_field, "END"))

        # Get the chromosome.
        chromosome = line_parts[0]

        # Add the SV type and positions to the lists.
        if svtype == "DEL" or svtype == "DUP":
            cnv_types.append(svtype)
            cnv_chromosomes.append(chromosome)
            cnv_start_positions.append(start_position)
            cnv_end_positions.append(end_position)

    # # Create the figure.
    # fig = plotly.graph_objs.Figure(
    #     data = [cnv_trace]
    # )


    # # Create the figure.
    # fig = plotly.graph_objs.Figure(
    #     data = [cnv_trace]
    # )

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

    # Create a shaded rectangle for each CNV, layering them below the CNV trace
    # and labeling them with the CNV type.
    for i in range(len(cnv_chromosomes)):
        fig.add_vrect(
            x0 = cnv_start_positions[i],
            x1 = cnv_end_positions[i],
            fillcolor = "Black",
            layer = "below",
            line_width = 0,
            opacity = 0.1,
            annotation_text = cnv_types[i],
            annotation_position = "top left",
            annotation_font_size = 20,
            annotation_font_color = "black"
        )

        # Add vertical lines at the start and end positions of the CNV.
        fig.add_vline(
            x = cnv_start_positions[i],
            line_width = 1,
            line_color = "black",
            layer = "below"
        )

        fig.add_vline(
            x = cnv_end_positions[i],
            line_width = 1,
            line_color = "black",
            layer = "below"
        )

        log.info("Added CNV rectangle for {}:{}-{}, SVLEN={}".format(
            cnv_chromosomes[i],
            cnv_start_positions[i],
            cnv_end_positions[i],
            cnv_end_positions[i] - cnv_start_positions[i]
        ))

    # Save the figure.
    filepath = os.path.join(output_path, "cnv_plot.html")
    plotly.offline.plot(fig, filename=filepath, auto_open=False)

    # Close the VCF file.
    vcf_file.close()

    print("Saved CNV plot to {}".format(filepath))
