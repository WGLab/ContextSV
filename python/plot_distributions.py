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


def generate_sv_size_plot(input_vcf, output_png, plot_title="SV Caller"):
    # Read VCF file into a pandas DataFrame
    vcf_df = pd.read_csv(input_vcf, sep='\t', comment='#', header=None, \
                         names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'], \
                            dtype={'CHROM': str, 'POS': np.int64, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, \
                                   'FILTER': str, 'INFO': str, 'FORMAT': str, 'SAMPLE': str})

    # Initialize dictionaries to store SV sizes for each type of SV
    sv_sizes = {}

    # Iterate over each record in the VCF file
    for _, record in vcf_df.iterrows():

        # Check if INFO/SVLEN is present. These are not present in BND and imprecise calls
        if 'SVLEN' not in record['INFO']:
            continue

        # Get the SV type
        sv_type = record['INFO'].split('SVTYPE=')[1].split(';')[0]

        # Get the SV size
        sv_size = int(record['INFO'].split('SVLEN=')[1].split(';')[0])

        # Print the SV type and size
        # print(f'SV type: {sv_type}, SV size: {sv_size}')

        # If the SV type is not in the dictionary, add it
        if sv_type not in sv_sizes:
            sv_sizes[sv_type] = []

        # Add the SV size to the dictionary
        sv_sizes[sv_type].append(sv_size)

    # Create a tiled plot where each tile shows the SV size distribution for a
    # different SV type
    sv_type_count = len(sv_sizes)
    fig, axes = plt.subplots(sv_type_count, 1, figsize=(10, 5 * sv_type_count))

    # Create a dictionary of SV types and their corresponding colors. Use light
    # colors to make the plot more readable, such as 'skyblue' for deletions,
    # 'lightgreen' for duplications, 'lightcoral' for inversions, and 'violet'
    # for insertions
    sv_colors = {'DEL': 'skyblue', 'DUP': 'lightgreen', 'INV': 'lightcoral', 'INS': 'violet'}

    # Create a dictionary of SV types and their corresponding labels
    sv_labels = {'DEL': 'Deletion', 'DUP': 'Duplication', 'INV': 'Inversion', 'INS': 'Insertion'}

    # Get the list of SV types and sort them in the order of the labels
    sv_types = sorted(sv_sizes.keys(), key=lambda x: sv_labels[x])

    # Print the number of SVs for each type, starting with the label
    print("SV Caller: ", plot_title)
    print('Number of SVs for each type:')
    for sv_type in sv_types:
        print(f'{sv_labels[sv_type]}: {len(sv_sizes[sv_type])}')


    # Plot the SV size distributions
    size_scale = 1000 # Convert SV sizes from bp to kb. Use abs() to handle negative deletion sizes
    for i, sv_type in enumerate(sv_types):
        sizes = np.array(sv_sizes[sv_type])
        axes[i].hist(np.abs(sizes) / size_scale, bins=100, color=sv_colors[sv_type], alpha=0.7, label=sv_labels[sv_type], edgecolor='black')
        axes[i].set_xlabel('SV size (kb)')
        axes[i].set_ylabel('Frequency')
        axes[i].set_title(f'{plot_title}: {sv_labels[sv_type]}')

        # Use a log scale for the y-axis
        axes[i].set_yscale('log')

    # Save the plot as a PNG file
    plt.tight_layout()
    plt.savefig(output_png)


if __name__ == '__main__':
    # Get the input and output file paths from the command line arguments
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    plot_title = sys.argv[3]

    print(f'Input file: {input_file}')
    print(f'Output file: {output_file}')

    # Generate the SV size plot
    generate_sv_size_plot(input_file, output_file, plot_title=plot_title)
