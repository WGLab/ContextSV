"""
plot_distributions.py: Plot the distributions of SV sizes in the input VCF file
and save the plot as a PNG file.

Usage:
    plot_distributions.py <input> <output>

Arguments:
    <input>     Path to the input VCF file.
    <output>    Path to the output PNG file.

Output:
    A PNG file with the SV size distributions.

Example:
    python plot_distributions.py input.vcf output.png
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# import plotly
import plotly.graph_objects as go

def generate_sv_size_plot(input_vcf, output_png, plot_title="SV Caller"):
    # Read VCF file into a pandas DataFrame
    vcf_df = pd.read_csv(input_vcf, sep='\t', comment='#', header=None, \
                         names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'], \
                            dtype={'CHROM': str, 'POS': np.int64, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, \
                                   'FILTER': str, 'INFO': str, 'FORMAT': str, 'SAMPLE': str})

    # Initialize dictionaries to store SV sizes for each type of SV
    sv_sizes = {}

    # Iterate over each record in the VCF file
    print("SV CALLER: ", plot_title)
    for _, record in vcf_df.iterrows():

        # Get the POS
        pos = record['POS']

        # Get the SV data by splitting semi-colon separated INFO field and
        # extracting SVTYPE and SVLEN
        info_fields = record['INFO'].split(';')
        sv_type = None
        sv_len = None  # INFO/SVLEN
        sv_span = None  # INFO/END - POS
        alignment = "NA"
        for field in info_fields:
            if field.startswith('SVTYPE='):
                sv_type = field.split('=')[1]                
            elif field.startswith('SVLEN='):
                sv_len = abs(int(field.split('=')[1]))
            elif field.startswith('END='):
                sv_span = int(field.split('=')[1]) - pos
            elif field.startswith('ALN='):
                alignment = field.split('=')[1]

        # Continue if SV type is BND (no SV size)
        if sv_type == "BND":
            continue
        # If the SV caller is DELLY, then we use the second SV size for non-INS
        # (they don't have SVLEN) and the first SV size for INS
        sv_size = None
        if plot_title == "DELLY" and sv_type != "INS":
            sv_size = sv_span
        else:
            sv_size = sv_len

        # If the plot title is GIAB, then we need to convert INS to DUP if
        # INFO/SVTYPE is INS and INFO/REPTYPE is DUP
        if plot_title == "GIAB" and sv_type == "INS":
            if 'REPTYPE=DUP' in record['INFO']:
                sv_type = "DUP"

        # Add the SV type if it's not in the dictionary
        if sv_type not in sv_sizes:
            sv_sizes[sv_type] = []

        # Add the SV size to the dictionary
        sv_sizes[sv_type].append(sv_size)

    # Create a tiled plot where each tile shows the SV size distribution for a
    # different SV type
    sv_type_count = len(sv_sizes)
    fig, axes = plt.subplots(sv_type_count, 1, figsize=(10, 5 * sv_type_count))
    print(f'Number of SV types: {sv_type_count}')

    # Create a dictionary of SV types and their corresponding colors.
    sv_colors = {'DEL': 'crimson', 'DUP': 'royalblue', 'INV': 'orange', 'INS': 'limegreen'}

    # Create a dictionary of SV types and their corresponding labels
    sv_labels = {'DEL': 'Deletion', 'DUP': 'Duplication', 'INV': 'Inversion', 'INS': 'Insertion'}

    # Get the list of SV types and sort them in the order of the labels
    sv_types = sorted(sv_sizes.keys(), key=lambda x: sv_labels[x])

    # Print the number of SVs for each type, starting with the label
    print("SV Caller: ", plot_title)
    print("Total number of SVs: ", len(vcf_df))

    print('Number of SVs for each type:')
    total_sv_count = 0
    for sv_type in sv_types:
        print(f'{sv_labels[sv_type]}: {len(sv_sizes[sv_type])}')
        total_sv_count += len(sv_sizes[sv_type])

    print(f'Total number of SVs (sum): {total_sv_count}')

    # Print the number of SVs for each type with size > 50kb
    print('Number of SVs for each type with size > 50kb:')
    for sv_type in sv_types:
        print(f'{sv_labels[sv_type]}: {len([x for x in sv_sizes[sv_type] if abs(x) > 50000])}')

    # Summary statistics
    all_sv_sizes = []
    for sv_type in sv_types:
        all_sv_sizes.extend(sv_sizes[sv_type])
    print('Summary statistics:')
    print(f'Minimum SV size: {min(all_sv_sizes)}')
    print(f'Maximum SV size: {max(all_sv_sizes)}')
    print(f'Mean SV size: {np.mean(all_sv_sizes)}')
    print(f'Median SV size: {np.median(all_sv_sizes)}')
    print(f'Standard deviation of SV sizes: {np.std(all_sv_sizes)}')
    print(f'Number of SVs >10kb: {len([x for x in all_sv_sizes if abs(x) > 10000])}')
    print(f'Number of SVs >50kb: {len([x for x in all_sv_sizes if abs(x) > 50000])}')
    print(f'Number of SVs >100kb: {len([x for x in all_sv_sizes if abs(x) > 100000])}')

    # Plot the SV size distributions
    size_scale = 1000 # Convert SV sizes from bp to kb. Use abs() to handle negative deletion sizes
    for i, sv_type in enumerate(sv_types):
        sizes = np.array(sv_sizes[sv_type])
        axes[i].hist(np.abs(sizes) / size_scale, bins=100, color=sv_colors[sv_type], alpha=0.7, label=sv_labels[sv_type], edgecolor='black')
        axes[i].set_xlabel('SV size (kb)')
        axes[i].set_ylabel('Frequency (log scale)')
        axes[i].set_title(f'{plot_title}: {sv_labels[sv_type]}')

        # Use a log scale for the y-axis
        axes[i].set_yscale('log')

        # # In the same axis, plot a known duplication if within the range of the plot
        if sv_type == 'DUP':
            print("TEST: Found DUP")
            cnv_size = 776237 / size_scale
            x_min, x_max = axes[i].get_xlim()
            if cnv_size > x_min and cnv_size < x_max:
                axes[i].axvline(x=cnv_size, color='black', linestyle='--')
            else:
                # Print the values
                print(f'CNV size: {cnv_size}, x_min: {x_min}, x_max: {x_max}')

        # Refresh the plot
        plt.draw()

    # Save the plot as a PNG file
    plt.tight_layout()
    plt.savefig(output_png)

    # Plot an additional plot with suffix _full.png that includes all SV types
    # (using plotly to avoid overlapping histograms)
    fig = go.Figure()
    for sv_type in sv_types:
        sizes = np.array(np.abs(sv_sizes[sv_type]))
        fig.add_trace(go.Histogram(x=sizes / size_scale, name=sv_labels[sv_type], marker_color=sv_colors[sv_type], xbins=dict(size=100)))

    # Update the layout to group the bars side by side
    fig.update_layout(
        barmode='group',
        title=f'{plot_title}: All SV types',
        xaxis_title='SV size (kb)',
        yaxis_title='Frequency (log scale)',
        yaxis_type='log',
        bargap=0.2,
    )

    # Move the legend to the top right inside the plot
    fig.update_layout(legend=dict(
        orientation='v',
        yanchor='top',
        y=0.75,
        xanchor='right',
        x=0.75,
    ))

    # Set a larger font size for all text in the plot
    fig.update_layout(font=dict(size=32))
    
    # # Save the plot as a high-resolution PNG file for using in posters
    fig.write_image(output_png.replace('.png', '_full.png'), width=1200, height=800)


if __name__ == '__main__':
    # Get the input and output file paths from the command line arguments
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    plot_title = sys.argv[3]

    print(f'Input file: {input_file}')
    print(f'Output file: {output_file}')

    # Generate the SV size plot
    generate_sv_size_plot(input_file, output_file, plot_title=plot_title)
