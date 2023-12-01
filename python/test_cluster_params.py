import os
import sys

import matplotlib.pyplot as plt

def get_precision_recall(file_path, sv_type='DEL'):
    """Parse text file containing epsilon, precision, and recall values."""
    epsilon_values = []
    precision_values = []
    recall_values = []
    fp_counts = []
    fn_counts = []

    with open(file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()

        epsilon = None
        sv_section_found = False
        for i, line in enumerate(lines):

            if "#EPSILON=" in line:
                epsilon = float(line.split('=')[1])
                epsilon_values.append(epsilon)

            # SV type sections
            elif "Running truvari" in line:
                if sv_type in line:
                    sv_section_found = True

            # Get the number of FPs
            elif "FP" in line and sv_section_found:
                # Get the value after the ':'
                fp = line.split(':')[1]

                # Clean up the string
                fp = fp.replace('\n', '')
                fp = fp.replace(',', '')
                fp = int(fp)
                fp_counts.append(fp)

            # Get the number of FNs
            elif "FN" in line and sv_section_found:
                # Get the value after the ':'
                fn = line.split(':')[1]

                # Clean up the string
                fn = fn.replace('\n', '')
                fn = fn.replace(',', '')
                fn = int(fn)
                fn_counts.append(fn)

            elif "precision" in line and sv_section_found:
                # Get the value after the ':'
                p = line.split(':')[1]

                # Clean up the string
                p = p.replace('\n', '')
                p = p.replace(',', '')
                p = float(p)
                precision_values.append(p)

            elif "recall" in line and sv_section_found:
                # Get the value after the ':'
                r = line.split(':')[1]

                # Clean up the string
                r = r.replace('\n', '')
                r = r.replace(',', '')
                r = float(r)
                recall_values.append(r)

                # Reset epsilon and sv_section_found
                epsilon = None
                sv_section_found = False

    # Get the maximum recall value, and then the maximum precision value at that
    # recall value
    max_recall = max(recall_values)
    max_precision = None
    max_index = None  # Index of the maximum recall and corresponding precision
    for i, recall in enumerate(recall_values):
        if recall == max_recall:
            if max_precision is None:
                max_precision = precision_values[i]
                max_index = i
            elif precision_values[i] > max_precision:
                max_precision = precision_values[i]
                max_index = i

    # Print the maximum precision and recall values
    print(f'SV Type: {sv_type}')
    print(f'Maximum Recall: {max_recall}')
    print(f'Maximum Precision at Maximum Recall: {max_precision}')

    # Print the FP and FN counts at the maximum recall and corresponding
    # precision
    print(f'FP Count at Maximum Recall: {fp_counts[max_index]}')
    print(f'FN Count at Maximum Recall: {fn_counts[max_index]}')

    return epsilon_values, precision_values, recall_values


def get_f1_scores(file_path, sv_type='DEL'):
    """Parse text file containing epsilon and F1 scores."""
    epsilon_values = []
    f1_values = []

    with open(file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()

        epsilon = None
        sv_section_found = False
        for i, line in enumerate(lines):

            if "#EPSILON=" in line:
                epsilon = float(line.split('=')[1])
                epsilon_values.append(epsilon)

            # SV type sections
            elif "Running truvari" in line:
                if sv_type in line:
                    sv_section_found = True

            elif "f1" in line and sv_section_found:
                # Get the value after the ':'
                f1 = line.split(':')[1]

                # Clean up the string
                f1 = f1.replace('\n', '')
                f1 = f1.replace(',', '')
                f1 = float(f1)
                f1_values.append(f1)
                
                # Reset epsilon and sv_section_found
                epsilon = None
                sv_section_found = False

    return epsilon_values, f1_values


def plot_precision_recall(epsilon, precision, recall, title="Precision and Recall vs. Epsilon", parameter_name='Epsilon'):
    """Plot precision and recall values."""
    # Create figure
    plt.figure()

    # Plot precision and recall vs. epsilon on same plot but different axes
    ax1 = plt.gca()
    ax2 = ax1.twinx()

    # Plot precision vs. epsilon on ax1
    ax1.plot(epsilon, precision, label='Precision', color='black')

    # Plot recall vs. epsilon on ax2
    ax2.plot(epsilon, recall, label='Recall', color='blue')

    # # Show ticks for all epsilon values
    # ax1.set_xticks(epsilon)

    # # Make X-ticks vertical
    # plt.xticks(rotation=90)

    # # Double the figure width
    # plt.gcf().set_size_inches(18.5, 10.5)

    # Add axis labels
    ax1.set_xlabel(parameter_name, color='black')
    ax1.set_ylabel('Precision', color='black')
    ax2.set_xlabel('Epsilon', color='black')
    ax2.set_ylabel('Recall', color='blue')

    # Set tick colors
    ax1.tick_params(axis='y', colors='black')
    ax2.tick_params(axis='y', colors='blue')

    # Add title
    plt.title(title)

    return plt


def plot_f1(epsilon, f1_scores, title="F1 vs. Epsilon", parameter_name='Epsilon'):
    """Plot F1 values."""
    # Create figure
    plt.figure()

    # Plot F1 vs. epsilon
    plt.plot(epsilon, f1_scores, label='F1')

    # # Show ticks for all epsilon values
    # plt.xticks(epsilon)

    # # Make X-ticks vertical
    # plt.xticks(rotation=90)

    # # Double the figure width
    # plt.gcf().set_size_inches(18.5, 10.5)

    # Add axis labels
    plt.xlabel(parameter_name)
    plt.ylabel('F1')

    # Add title
    plt.title(title)

    # Return figure
    return plt


if __name__ == '__main__':
    # Take in benchmark file path as command line argument
    file_path = sys.argv[1]

    print(f'Input file path: {file_path}')

    # Take in cluster type as command line argument
    cluster_type = sys.argv[2]
    if cluster_type not in ['dbscan', 'agglo']:
        print(f"Invalid cluster type: {cluster_type}")
        sys.exit(1)

    # Get the cluster type string
    cluster_string = 'DBSCAN' if cluster_type == 'dbscan' else 'Agglomerative'

    # Determine the parameter to test
    parameter_name = 'Epsilon' if cluster_type == 'dbscan' else 'Distance Threshold'

    # Create the plot title
    plot_title = cluster_string + ' Cluster + Merge'

    # Create directory to store plots
    output_dir = cluster_type + '_tests'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print(f"Plotting precision and recall vs. {parameter_name}...")

    # Plot precision and recall values
    eps, prec, rec = get_precision_recall(file_path, sv_type='DEL')
    fig = plot_precision_recall(eps, prec, rec, title=plot_title + ' (DEL)', parameter_name=parameter_name)
    fig.savefig(output_dir + '/Precision_Recall_DEL.png')

    eps, prec, rec = get_precision_recall(file_path, sv_type='DUP')
    fig = plot_precision_recall(eps, prec, rec, title=plot_title + ' (DUP)', parameter_name=parameter_name)
    fig.savefig(output_dir + '/Precision_Recall_DUP.png')

    # Plot F1 scores
    eps, f1 = get_f1_scores(file_path, sv_type='DEL')
    fig = plot_f1(eps, f1, title=plot_title + ' (DEL)', parameter_name=parameter_name)
    fig.savefig(output_dir + '/F1_DEL.png')

    eps, f1 = get_f1_scores(file_path, sv_type='DUP')
    fig = plot_f1(eps, f1, title=plot_title + ' (DUP)', parameter_name=parameter_name)
    fig.savefig(output_dir + '/F1_DUP.png')