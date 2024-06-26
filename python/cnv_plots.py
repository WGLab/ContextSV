"""Plot the copy number variants and their log2_ratio, BAF values."""

import os
import sys
import logging as log
import plotly
from plotly.subplots import make_subplots
import pandas as pd

from .utils import parse_region, get_info_field_column, get_info_field_value

# Set up logging.
log.basicConfig(
    level=log.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        log.StreamHandler(sys.stdout)
    ]
)

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
    snp_cnvs_only = False

    # Set the maximum number of CNVs to plot.
    max_cnvs = 10

    # Parse the region.
    chromosome, start_position, end_position = parse_region(region)

    if start_position is not None and end_position is not None:
        log.info("Plotting CNVs in region %s:%d-%d.", chromosome, start_position, end_position)
    else:
        log.info("Plotting CNVs in region %s.", chromosome)

    # Set up the CNV type string list for the 6 CNV states.
    cnv_types = ["NAN", "DEL", "DEL", "NEUT", "NEUT", "DUP", "DUP"]

    # Filter the VCF file to the region using pandas, and make the chromosome
    # column a string.
    log.info("Loading VCF data from %s", vcf_file)
    try:
        vcf_data = pd.read_csv(vcf_file, sep="\t", comment="#", header=None, dtype={0: str})
    except pd.errors.EmptyDataError:
        log.info("No variants found in %s.", vcf_file)
        return
    
    if start_position is not None and end_position is not None:
        vcf_data = vcf_data[(vcf_data[0] == chromosome) & (vcf_data[1] >= start_position) & (vcf_data[1] <= end_position)]
    else:
        vcf_data = vcf_data[(vcf_data[0] == chromosome)]

    log.info("Loaded %d variants from %s", len(vcf_data), vcf_file)

    # Filter the CNV data to the region using pandas, and make the chromosome
    # column a string.
    log.info("Loading CNV data from %s", cnv_data_file)
    cnv_data = pd.read_csv(cnv_data_file, sep="\t", header=0, dtype={"chromosome": str})
    if start_position is not None and end_position is not None:
        cnv_data = cnv_data[(cnv_data["chromosome"] == chromosome) & (cnv_data["position"] >= start_position) & (cnv_data["position"] <= end_position)]
    else:
        cnv_data = cnv_data[(cnv_data["chromosome"] == chromosome)]

    # Filter the CNV data to only include SNPs.
    cnv_data = cnv_data[cnv_data["snp"] == 1]

    log.info("Loaded %d SNPs from %s", len(cnv_data), cnv_data_file)

    # Create an output html file where we will append the CNV plots.
    if start_position is not None and end_position is not None:
        html_filename = f"cnv_plots_{chromosome}_{start_position}_{end_position}.html"
    else:
        html_filename = f"cnv_plots_{chromosome}.html"
    
    # Create the output directory if it doesn't exist.
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Create the output html file.
    output_html_filepath = os.path.join(output_path, html_filename)
    if os.path.exists(output_html_filepath):
        os.remove(output_html_filepath)

    with open(output_html_filepath, "w", encoding="utf-8") as output_html_file:

        # Get the INFO field column by searching for the first column that
        # contains SVTYPE=.
        info_column = get_info_field_column(vcf_data)

        # Loop through the VCF data and plot each CNV (DEL or DUP) along with log2
        # ratio and BAF values for the SNPs in the CNV.
        cnv_count = 0
        for _, sv_data in vcf_data.iterrows():

            # Get the INFO field
            info_field = sv_data[info_column]

            # Get the SVTYPE field value.
            svtype = get_info_field_value(info_field, "SVTYPE")
            if svtype == "INS" and get_info_field_value(info_field, "REPTYPE") == "DUP":
                svtype = "DUP"

            # Get the ALN field value (alignment type used to call the SV).
            aln = get_info_field_value(info_field, "ALN")

            # # Skip the SV if SNP CNV data was not used to call it.
            if snp_cnvs_only and "SNPCNV" not in aln:
                continue

            log.info("Found CNV %s %s:%d-%d, LEN=%d", svtype, sv_data[0], sv_data[1], sv_data[1] + int(get_info_field_value(info_field, "SVLEN")) - 1, int(get_info_field_value(info_field, "SVLEN")))

            # Analyze the CNV if it is a DEL or DUP (=INS with INFO/REPTYE=DUP)
            if svtype in ("DEL", "DUP"):

                # Get the read support for the CNV.
                read_support = int(get_info_field_value(info_field, "SUPPORT"))

                # # Skip the CNV if the support is < 2.
                # if read_support < 2:
                #     continue

                # Get the start position.
                start_position = int(sv_data[1])

                # Get the SV length.
                cnv_length = int(get_info_field_value(info_field, "SVLEN"))
                log.info("CNV length: %d", cnv_length)

                # Continue if the CNV length is < 100kb.
                if abs(cnv_length) < 100000:
                    continue

                # Use absolute value of CNV length (deletions are negative).
                cnv_length = abs(cnv_length)

                # Get the end position using the start position and SV length.
                end_position = start_position + cnv_length - 1

                # Get the chromosome.
                chromosome = sv_data[0]

                # Get the plot range as a multiple of the CNV length.
                plot_start_position = start_position - (cnv_length/2)
                plot_end_position = end_position + (cnv_length/2)

                # Get the CNV state, log2 ratio, and BAF values for all SNPs in the
                # plot range.
                log.info("Getting SNPs in CNV %s %s:%d-%d.", svtype, chromosome, plot_start_position, plot_end_position)
                # sv_data = cnv_data[(cnv_data["position"] >=
                # plot_start_position) & (cnv_data["position"] <=
                # plot_end_position)]
                sv_data = cnv_data[(cnv_data["position"] >= start_position) & (cnv_data["position"] <= end_position)]
                
                # If there are no SNPs in the plot range, skip the CNV.
                if len(sv_data) == 0:
                    log.info("No SNPs found in CNV %s %s:%d-%d.", svtype, chromosome, start_position, end_position)
                    continue
                else:
                    log.info("Found %d SNPs in CNV %s %s:%d-%d.", len(sv_data), svtype, chromosome, start_position, end_position)

                # Get the marker colors for the state sequence.
                marker_colors = []
                for state in sv_data["cnv_state"]:
                    if state in [1, 2]:
                        marker_colors.append("red")
                    elif state in [3, 4]:
                        marker_colors.append("black")
                    elif state in [5, 6]:
                        marker_colors.append("blue")

                # Now get the values before and after the CNV.
                sv_data_before = cnv_data[(cnv_data["position"] >= plot_start_position) & (cnv_data["position"] < start_position)]
                sv_data_after = cnv_data[(cnv_data["position"] > end_position) & (cnv_data["position"] <= plot_end_position)]

                # Set the marker colors for the SNPs before and after the CNV to
                # gray.
                marker_colors_before = ["gray"] * len(sv_data_before)
                marker_colors_after = ["gray"] * len(sv_data_after)

                # Concatenate the SNP data before, during, and after the CNV.
                sv_data = pd.concat([sv_data_before, sv_data, sv_data_after])

                # Concatenate the marker colors before, during, and after the
                # CNV.
                marker_colors = marker_colors_before + marker_colors + marker_colors_after

                # Get the hover text for the state sequence markers.
                hover_text = []
                for _, row in sv_data.iterrows():
                    hover_text.append(f"SNP: {row['snp']}<br>TYPE: {cnv_types[row['cnv_state']]}<br>CHR: {row['chromosome']}<br>POS: {row['position']}<br>L2R: {row['log2_ratio']}<br>BAF: {row['b_allele_freq']}<br>PFB: {row['population_freq']}")
                    # hover_text.append(f"TYPE: {cnv_types[row['cnv_state']]}<br>CHR: {row['chromosome']}<br>POS: {row['position']}<br>L2R: {row['log2_ratio']}<br>BAF: {row['b_allele_freq']}<br>PFB: {row['population_freq']}")

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
                        color = "black",
                        width = 0
                    ),
                    showlegend = False
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
                        color = "black",
                        width = 0
                    ),
                    showlegend = False
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

                # Set the Y-axis range for the log2 ratio plot.
                fig.update_yaxes(
                    range = [-1.2, 1.2],
                    row = 1,
                    col = 1
                )

                # Set the Y-axis range for the BAF plot.
                fig.update_yaxes(
                    range = [-0.2, 1.2],
                    row = 2,
                    col = 1
                )

                # Set the figure title.
                fig.update_layout(
                    title_text = f"{svtype} (SUPPORT={read_support}, LEN={cnv_length}bp) at {chromosome}:{start_position}-{end_position} [ALN={aln}]",
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

                # Append the figure to the output html file.
                output_html_file.write(fig.to_html(full_html=False, include_plotlyjs="cdn"))
                log.info("Plotted CNV %s %s:%d-%d.", svtype, chromosome, start_position, end_position)

                # Increment the CNV count.
                cnv_count += 1

                # Break if the maximum number of CNVs has been reached.
                if cnv_count == max_cnvs:
                    break
    
    log.info("Saved CNV plots to %s.", output_html_filepath)
