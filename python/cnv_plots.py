# Plot the copy number variants and their log2_ratio, BAF values.

import os
import sys
import logging as log
import plotly
from plotly.subplots import make_subplots
import pandas as pd

from .utils import parse_region

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

def run(vcf_file, cnv_data_file, output_path, region):
    """
    Saves a plot of the CNVs and their log2 ratio and B-allele frequency
    values.
    
    Args:
        vcf_path (str): The path to the VCF file.
        cnv_data_path (str): The path to the CNV data file.
        output_path (str): The path to the output directory.

    Returns:
        None
    """
    # Set the maximum number of CNVs to plot.
    max_cnvs = 10

    # Set the plot range for each CNV to 3 times its length on each side.
    plot_range = 3

    # Parse the region.
    chromosome, start_position, end_position = parse_region(region)

    # Set up the CNV type string list for the 6 CNV states.
    cnv_types = ["NAN", "DEL", "DEL", "NEUT", "NEUT", "DUP", "DUP"]

    # Read the VCF file.
    # vcf_file = open(vcf_path, "r")

    # Filter the VCF file to the region using pandas.
    vcf_data = pd.read_csv(vcf_file, sep="\t", comment="#", header=None)
    if start_position is not None and end_position is not None:
        vcf_data = vcf_data[(vcf_data[0] == chromosome) & (vcf_data[1] >= start_position) & (vcf_data[1] <= end_position)]
    else:
        vcf_data = vcf_data[(vcf_data[0] == chromosome)]

    # Filter the CNV data to the region using pandas.
    cnv_data = pd.read_csv(cnv_data_file, sep="\t")
    if start_position is not None and end_position is not None:
        cnv_data = cnv_data[(cnv_data["chromosome"] == chromosome) & (cnv_data["position"] >= start_position) & (cnv_data["position"] <= end_position)]
    else:
        cnv_data = cnv_data[(cnv_data["chromosome"] == chromosome)]

    # Create an output html file where we will append the CNV plots.
    if start_position is not None and end_position is not None:
        html_filename = "cnv_plots_{}_{}_{}.html".format(chromosome, start_position, end_position)
    else:
        html_filename = "cnv_plots_{}.html".format(chromosome)
    
    output_html_file = open(os.path.join(output_path, html_filename), "w")

    #output_html_file = open(os.path.join(output_path, "cnv_plots.html"), "w")

    # Loop through the VCF data and plot each CNV (DEL or DUP) along with log2
    # ratio and BAF values for the SNPs in the CNV.
    cnv_count = 0
    for index, row in vcf_data.iterrows():
        # Get the INFO field.
        info_field = row[7]

        # Get the SVTYPE field value.
        svtype = get_info_field_value(info_field, "SVTYPE")

        # Analyze the CNV if it is a DEL or DUP.
        if svtype == "DEL" or svtype == "DUP":

            # Get the read depth (DP) value.
            read_depth = int(get_info_field_value(info_field, "DP"))

            # Get the start and end positions.
            start_position = int(row[1])
            end_position = int(get_info_field_value(info_field, "END"))

            # Get the chromosome.
            chromosome = row[0]

            # Get the length of the CNV.
            cnv_length = end_position - start_position + 1

            # Get the plot range.
            plot_start_position = start_position - (plot_range * cnv_length)
            plot_end_position = end_position + (plot_range * cnv_length)

            # Get the CNV state, log2 ratio, and BAF values for all SNPs in the
            # plot range.
            sv_data = cnv_data[(cnv_data["position"] >= plot_start_position) & (cnv_data["position"] <= plot_end_position)]
            
            # Get the marker colors for the state sequence.
            marker_colors = []
            for state in sv_data["cnv_state"]:
                if state in [1, 2]:
                    marker_colors.append("red")
                elif state in [3, 4]:
                    marker_colors.append("black")
                elif state in [5, 6]:
                    marker_colors.append("blue")

            # Get the hover text for the state sequence markers.
            hover_text = []
            for index, row in sv_data.iterrows():
                hover_text.append(f"TYPE: {cnv_types[row['cnv_state']]}<br>CHR: {row['chromosome']}<br>POS: {row['position']}<br>L2R: {row['log2_ratio']}<br>BAF: {row['b_allele_freq']}")

            # Create the log2 ratio trace.
            log2_ratio_trace = plotly.graph_objs.Scatter(
                x = sv_data["position"],
                y = sv_data["log2_ratio"],
                mode = "markers+lines",
                name = r'Log<sub>2</sub> Ratio',
                text = hover_text,
                hoverinfo = "text",
                marker = dict(
                    color = marker_colors,
                    size = 10,
                ),
                line = dict(
                    color = "gray",
                    width = 2
                )
            )

            # Create the B-allele frequency trace.
            baf_trace = plotly.graph_objs.Scatter(
                x = sv_data["position"],
                y = sv_data["b_allele_freq"],
                mode = "markers+lines",
                name = "B-Allele Frequency",
                text = hover_text,
                marker = dict(
                    color = marker_colors,
                    size = 10,
                ),
                line = dict(
                    color = "gray",
                    width = 2
                ),
            )

            # Create a subplot for the CNV plot and the BAF plot.
            fig = make_subplots(
                rows=2,
                cols=1,
                shared_xaxes=True,
                vertical_spacing=0.05,
                subplot_titles=(r"SNP Log<sub>2</sub> Ratio", "SNP B-Allele Frequency")
            )

            # Add the traces to the figure.
            fig.append_trace(log2_ratio_trace, 1, 1)
            fig.append_trace(baf_trace, 2, 1)

            # Set the x-axis title.
            fig.update_xaxes(
                title_text = "Chromosome Position",
                row = 2,
                col = 1
            )

            # Set the y-axis titles.
            fig.update_yaxes(
                title_text = r"Log<sub>2</sub> Ratio",
                row = 1,
                col = 1
            )

            fig.update_yaxes(
                title_text = "B-Allele Frequency",
                row = 2,
                col = 1
            )

            # Set the figure title.
            fig.update_layout(
                title_text = f"CNV {svtype} {chromosome}:{start_position}-{end_position} (Read Depth: {read_depth})"
            )

            # Create a shaded rectangle for the CNV, layering it below the CNV
            # trace and labeling it with the CNV type.
            fig.add_vrect(
                x0 = start_position,
                x1 = end_position,
                fillcolor = "Black",
                layer = "below",
                line_width = 0,
                opacity = 0.1,
                annotation_text = svtype,
                annotation_position = "top left",
                annotation_font_size = 20,
                annotation_font_color = "black"
            )

            # Add vertical lines at the start and end positions of the CNV.
            fig.add_vline(
                x = start_position,
                line_width = 2,
                line_color = "black",
                layer = "below"
            )

            fig.add_vline(
                x = end_position,
                line_width = 2,
                line_color = "black",
                layer = "below"
            )

            # Hide the legends.
            fig.update_layout(
                showlegend = False
            )

            # Add the figure to the output html file.
            output_html_file.write(fig.to_html(full_html=False, include_plotlyjs="cdn"))

            log.info("Plotted CNV {} {}:{}-{}.".format(svtype, chromosome, start_position, end_position))

            # # Save the figure.
            # filepath = os.path.join(output_path, f"cnv_{svtype}_{chromosome}_{start_position}_{end_position}.html")
            # plotly.offline.plot(fig, filename=filepath, auto_open=False)

            # log.info("Saved CNV plot to {}".format(filepath))

            # Increment the CNV count.
            cnv_count += 1

            # Break if the maximum number of CNVs has been reached.
            if cnv_count == max_cnvs:
                break

    # Close the output html file.
    output_html_file.close()
    
    log.info("Saved CNV plots to {}".format(os.path.join(output_path, html_filename)))
