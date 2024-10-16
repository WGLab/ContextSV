"""Plot the copy number variants and their log2_ratio, BAF values."""

import os
import sys
import logging as log
import plotly
from plotly.subplots import make_subplots
import pandas as pd

try:
    from .utils import parse_region, get_info_field_column, get_info_field_value
except ImportError:
    from utils import parse_region, get_info_field_column, get_info_field_value

MIN_CNV_LENGTH = 10000

# Set up logging.
log.basicConfig(
    level=log.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        log.StreamHandler(sys.stdout)
    ]
)

def parse_region(region):
    """
    Parses the region string to get the chromosome, start position, and end
    position.
    
    Args:
        region (str): The region string in the format "chr:start-end".
    
    Returns:
        tuple: A tuple containing the chromosome, start position, and end
        position.
    """

    # Split the region string by ":" and "-".
    region_parts = region.split(":")
    chromosome = region_parts[0]
    region_parts = region_parts[1].split("-")
    start_position = int(region_parts[0])
    end_position = int(region_parts[1])

    return chromosome, start_position, end_position

def run(cnv_data_file, output_html):
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

    # Filter the CNV data to the region using pandas, and make the chromosome
    # column a string.
    log.info("Loading CNV data from %s", cnv_data_file)

    # Read the first 3 lines of the file to get metadata.
    # Metadata is formatted as follows:
    # "SVTYPE="
    # "POS="
    # "HMM_LOGLH="
    metadata = {}
    metadata_row_count = 3
    with open(cnv_data_file, "r", encoding="utf-8") as f:
        # Read the first 3 lines of the file.
        for _ in range(metadata_row_count):
            line = f.readline().strip()
            if '=' in line:
                key, value = line.split("=")
                log.info("Metadata: %s=%s", key, value)
                value = value.strip()
                metadata[key] = value

    sv_type = metadata["SVTYPE"]
    position = metadata["POS"]
    chromosome, start_position, end_position = parse_region(position)
    hmm_loglh = float(metadata["HMM_LOGLH"])

    # Extract information from the metadata.
    log.info("SV type: %s, chromosome: %s, start position: %d, end position: %d, HMM log likelihood: %f", sv_type, chromosome, start_position, end_position, hmm_loglh)

    # Read the CNV data from the file.
    sv_data = pd.read_csv(cnv_data_file, sep="\t", header=metadata_row_count, dtype={"chromosome": str})
    if len(sv_data) == 0:
        log.info("No predictions found in %s", cnv_data_file)
        return
    else:
        log.info("Found %d predictions in %s", len(sv_data), cnv_data_file)

    # Create an output html file where we will append the CNV plots.
    if start_position is not None and end_position is not None:
        html_filename = f"cnv_plots_{chromosome}_{start_position}_{end_position}.html"
    else:
        html_filename = f"cnv_plots_{chromosome}.html"

    # Create the output html file.
    if os.path.exists(output_html):
        os.remove(output_html)

    with open(output_html, "w", encoding="utf-8") as output_html_file:

        # Use absolute value of CNV length (deletions are negative).
        # cnv_length = abs(cnv_length)
        cnv_length = end_position - start_position + 1

        # Return if the CNV length is less than the minimum CNV length.
        if cnv_length < MIN_CNV_LENGTH:
            log.info("Skipping CNV %s:%d-%d due to length < %d.", chromosome, start_position, end_position, MIN_CNV_LENGTH)
            return

        # Get the plot range as the minimum and maximum positions in the CNV
        # data.
        plot_start_position = sv_data["position"].min()
        plot_end_position = sv_data["position"].max()

        # Get the CNV state, log2 ratio, and BAF values for all SNPs in the
        # plot range.
        log.info("Getting SNPs in CNV %s:%d-%d.", chromosome, plot_start_position, plot_end_position)

        # If there are no SNPs in the plot range, skip the CNV.
        if len(sv_data) == 0:
            log.info("No SNPs found in CNV %s:%d-%d.", chromosome, start_position, end_position)
            # continue
        else:
            log.info("Found %d SNPs in CNV %s:%d-%d.", len(sv_data), chromosome, start_position, end_position)

        # Get the marker colors for the state sequence.
        marker_colors = []
        for state in sv_data["cnv_state"]:
            if state in [1, 2]:
                marker_colors.append("red")
            elif state in [3, 4]:
                marker_colors.append("black")
            elif state in [5, 6]:
                marker_colors.append("blue")

        # [TEST] Set the marker colors for the SNPs before and after the CNV to
        # gray (no state prediction).
        # for i in range(len(sv_data)):
        #     if sv_data["position"].iloc[i] < start_position or sv_data["position"].iloc[i] > end_position:
        #         marker_colors[i] = "gray"

        
        # Use row['snp'] to get whether SNP or not (0=not SNP, 1=SNP).
        marker_symbols = ["circle" if snp == 1 else "circle-open" for snp in sv_data["snp"]]
        # marker_symbols = marker_symbols_before + marker_symbols + marker_symbols_after

        # Concatenate the SNP data before, during, and after the CNV.
        # sv_data = pd.concat([sv_data_before, sv_data, sv_data_after])

        # Set all -1 B-allele frequency values to 0.
        sv_data.loc[sv_data["b_allele_freq"] == -1, "b_allele_freq"] = 0

        # Get the hover text for the state sequence markers.
        hover_text = []
        for _, row in sv_data.iterrows():
            hover_text.append(f"SNP: {row['snp']}<br>TYPE: {'NA'}<br>CHR: {row['chromosome']}<br>POS: {row['position']}<br>L2R: {row['log2_ratio']}<br>BAF: {row['b_allele_freq']}<br>PFB: {row['population_freq']}<br>STATE: {row['cnv_state']}")
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
                symbol = marker_symbols
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
                symbol = marker_symbols
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
            range = [-2.0, 2.0],
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
        # fig.update_layout(
        #     title_text = f"{svtype} (SUPPORT={read_support}, LEN={cnv_length}bp) at {chromosome}:{start_position}-{end_position} [ALN={aln}]",
        # )

        # Create a shaded rectangle for the CNV, layering it below the CNV
        # trace and labeling it with the CNV type.
        fig.add_vrect(
            x0 = start_position,
            x1 = end_position,
            fillcolor = "Black",
            layer = "below",
            line_width = 0,
            opacity = 0.1,
            annotation_text = '',
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
        log.info("Plotted CNV %s %s:%d-%d.", 'SVType', chromosome, start_position, end_position)

        # Increment the CNV count.
        # cnv_count += 1

        # # Break if the maximum number of CNVs has been reached.
        # if cnv_count == max_cnvs:
        #     break
    
    log.info("Saved CNV plots to %s.", output_html)

if __name__ == "__main__":
    cnv_data_file = sys.argv[1]
    output_path = sys.argv[2]

    run(cnv_data_file, output_path)
