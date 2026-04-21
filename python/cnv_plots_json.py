"""
pip install plotly kaleido==0.2.1
"""
import os
import argparse
import json
import numpy as np
import plotly
from plotly.subplots import make_subplots

min_sv_length = 50000 # Minimum SV length in base pairs

# Set up argument parser
parser = argparse.ArgumentParser(description='Generate CNV plots from JSON data.')
parser.add_argument('json_file', type=str, help='Path to the JSON file containing SV data')
parser.add_argument('chromosome', type=str, help='Chromosome to filter the SVs by (e.g., "chr3")', nargs='?', default=None)
parser.add_argument('--formats', type=str, default='html,svg', help='Comma-separated output formats (e.g., html,svg,pdf,png)')
parser.add_argument('--width', type=int, default=1800, help='Figure width in pixels for static exports')
parser.add_argument('--height', type=int, default=1200, help='Figure height in pixels for static exports')
parser.add_argument('--scale', type=float, default=2.0, help='Scale factor for raster exports (png,jpg,webp)')
args = parser.parse_args()

output_formats = [fmt.strip().lower() for fmt in args.formats.split(',') if fmt.strip()]

repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
output_dir = os.path.join(repo_root, 'linktoscripts', 'CNV_Plots')
os.makedirs(output_dir, exist_ok=True)

# Load your JSON data
with open(args.json_file) as f:
    sv_data = json.load(f)

# State marker colors
# https://community.plotly.com/t/plotly-colours-list/11730/6
state_colors_dict = {
    '1': 'darkred',
    '2': 'red',
    '3': 'gray',
    '4': 'green',
    '5': 'blue',
    '6': 'darkblue',
}

sv_type_dict = {
    'DEL': 'Deletion',
    'DUP': 'Duplication',
    'INV': 'Inversion'
}

# Loop through each SV (assuming your JSON contains multiple SVs)
skip_count_chrom = 0
skip_count_length = 0
save_count = 0
for sv in sv_data:

    # If a chromosome is specified, filter the SVs by that chromosome
    if args.chromosome and sv['chromosome'] != args.chromosome:
        skip_count_chrom += 1
        continue

    # Filter out SVs that are smaller than the minimum length
    if np.abs(sv['size']) < min_sv_length:
        skip_count_length += 1
        continue

    # Extract data for plotting
    positions_before = sv['before_sv']['positions']
    b_allele_freq_before = sv['before_sv']['b_allele_freq']
    positions_after = sv['after_sv']['positions']
    b_allele_freq_after = sv['after_sv']['b_allele_freq']

    # Create a subplot for the CNV plot and the BAF plot.
    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.05,
        subplot_titles=(r"SNP Log<sub>2</sub> Ratio", "SNP B-Allele Frequency")
        )

    # Get the chromosome, start, end, and sv_type from the SV data
    chromosome = sv['chromosome']
    start = sv['start']
    end = sv['end']
    sv_type = sv['sv_type']
    likelihood = sv['likelihood']
    sv_length = sv['size']

    # Plot the data for 'before_sv', 'sv', and 'after_sv'
    for section in ["before_sv", "sv", "after_sv"]:
        positions = sv[section]['positions']
        b_allele_freq = sv[section]['b_allele_freq']
        population_freq = sv[section]['population_freq']
        log2_ratio = sv[section]['log2_ratio']
        is_snp = sv[section]['is_snp']

        # Set all b-allele frequencies to NaN if not SNPs
        b_allele_freq = [freq if is_snp_val else float('nan') for freq, is_snp_val in zip(b_allele_freq, is_snp)]

        if section == "sv":
            # is_snp = sv[section]['is_snp']
            states = sv[section]['states']
            state_colors = [state_colors_dict[str(state)] for state in states]
            marker_symbols = ['circle' if is_snp_val else 'circle-open' for is_snp_val in is_snp]

            # Set the hover text
            hover_text = []
            for i, position in enumerate(positions):
                # Add hover text for each point
                hover_text.append(
                    f"Position: {position}<br>"
                    f"State: {states[i]}<br>"
                    f"Log2 Ratio: {log2_ratio[i]}<br>"
                    f"SNP: {is_snp[i]}<br>"
                    f"BAF: {b_allele_freq[i]}<br>"
                    f"Population Frequency: {population_freq[i]}<br>"
                )
        else:
            # is_snp = sv[section]['is_snp']
            state_colors = ['black'] * len(positions)
            # marker_symbols = ['circle-open'] * len(positions)
            marker_symbols = ['circle' if is_snp_val else 'circle-open' for is_snp_val in is_snp]
            hover_text = []
            for i, position in enumerate(positions):
                # Add hover text for each point
                hover_text.append(
                    f"Position: {position}<br>"
                    f"Log2 Ratio: {log2_ratio[i]}<br>"
                    f"BAF: {b_allele_freq[i]}<br>"
                    f"Population Frequency: {population_freq[i]}<br>"
                )

        # Create the log2 trace
        log2_trace = plotly.graph_objs.Scatter(
            x=positions,
            y=log2_ratio,
            mode='markers+lines',
            name=r'Log<sub>2</sub> Ratio',
            text=hover_text,
            hoverinfo='text',
            marker=dict(
                color=state_colors,
                size=5,
                symbol=marker_symbols,
            ),
            line=dict(
                color='black',
                width=0
            ),
            showlegend=False
        )

        # Create the BAF trace
        baf_trace = plotly.graph_objs.Scatter(
            x=positions,
            y=b_allele_freq,
            mode='markers+lines',
            name='B-Allele Frequency',
            text=hover_text,
            hoverinfo='text',
            marker=dict(
                color=state_colors,
                size=5,
                symbol=marker_symbols,
            ),
            line=dict(
                color='black',
                width=0
            ),
            showlegend=False
        )

        if section == "sv":
            # Create a shaded rectangle for the CNV, layering it below the CNV
            # trace and labeling it with the CNV type.
            fig.add_vrect(
                x0 = start,
                x1 = end,
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
                x = start,
                line_width = 2,
                line_color = "black",
                layer = "below"
            )

            fig.add_vline(
                x = end,
                line_width = 2,
                line_color = "black",
                layer = "below"
            )

        # Add traces to the figure
        fig.append_trace(log2_trace, row=1, col=1)
        fig.append_trace(baf_trace, row=2, col=1)
    
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

    # Set the title of the plot.
    fig.update_layout(
        title_text = f"{sv_type_dict[sv_type]} at {chromosome}:{start}-{end} ({sv_length} bp)",
        title_x = 0.5,
        showlegend = False,
        template = 'simple_white',
        font = dict(family='Arial', size=20, color='black'),
        width = args.width,
        height = args.height,
        margin = dict(l=100, r=30, t=120, b=90)
    )

    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True, ticks='outside')
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True, ticks='outside')
    #     height = 800,
    #     width = 800
    # )
    # Save plots into a dedicated repository output directory.
    svlen_kb = sv_length // 1000
    base_name = f"SV_{chromosome}_{start}_{end}_{sv_type}_{svlen_kb}kb"

    if 'html' in output_formats:
        html_path = os.path.join(output_dir, f"{base_name}.html")
        fig.write_html(html_path)
        print(f"Plot saved as {html_path}")

    static_formats = {'svg', 'pdf', 'png', 'jpg', 'jpeg', 'webp', 'eps'}
    requested_static_formats = [fmt for fmt in output_formats if fmt in static_formats]

    if requested_static_formats:
        try:
            for fmt in requested_static_formats:
                out_path = os.path.join(output_dir, f"{base_name}.{fmt}")
                fig.write_image(out_path, format=fmt, width=args.width, height=args.height, scale=args.scale)
                print(f"Plot saved as {out_path}")
        except ValueError as e:
            print("Static image export requires Kaleido. Install with: pip install -U kaleido")
            raise e

    save_count += 1

print(f"Finished processing {save_count} SVs. Skipped {skip_count_chrom} SVs due to chromosome filter and {skip_count_length} SVs due to length filter.")

