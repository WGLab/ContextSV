"""
sv_merger.py
Use DBSCAN to merge SVs with the same breakpoint.
Mode can be 'dbscan', 'gmm', or 'agglomerative'.
https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html

Usage: python sv_merger.py <VCF file path>
Output: <VCF file path>.merged.vcf
"""

import os
import numpy as np
import pandas as pd

import logging
logging.basicConfig(level=logging.INFO)

import matplotlib.pyplot as plt  # For plotting merge behavior

# DBSCAN clustering algorithm
from sklearn.cluster import DBSCAN

# Agglomerative clustering algorithm
from sklearn.cluster import AgglomerativeClustering


def plot_dbscan(breakpoints, chosen_breakpoints, filename='dbscan_clustering.png'):
    """
    Plot the DBSCAN clustering behavior for SV breakpoints.
    """
    # logging.info the filename
    logging.info(f"Plotting DBSCAN clustering behavior to {filename}...")

    # logging.info all breakpoints
    logging.info(f"Breakpoints:")
    for i in range(breakpoints.shape[0]):
        logging.info(f"Row {i+1} - Breakpoints: {breakpoints[i, :]}")

    # Remove the chosen breakpoints from the breakpoints array
    breakpoints = np.delete(breakpoints, np.where(breakpoints == chosen_breakpoints), axis=0)

    # Plot the SV breakpoints as individual lines in each row, and the chosen
    # SV breakpoint as a red line at the top
    # Create a new figure
    plt.close()
    plt.clf()
    plt.cla()
    plt.figure(figsize=(10, 10))
    for i in range(breakpoints.shape[0]):
        row = i+1
        plt.plot(breakpoints[i, :], [row, row], 'b-')

    plt.plot(chosen_breakpoints, [0, 0], 'r-')
    logging.info(f"Chosen breakpoints: {chosen_breakpoints}")

    # Set plot labels
    plt.title('DBSCAN Clustering Behavior')
    plt.xlabel('Breakpoint Position')
    plt.ylabel('SVs')
    plt.legend()

    # Save the plot
    plt.savefig(filename)


def sv_merger(vcf_file_path, mode='dbscan', eps=100, min_samples=2, suffix='.merged'):
    """
    Use DBSCAN to merge SVs with the same breakpoint.
    Mode can be 'dbscan', 'gmm', or 'agglomerative'.
    https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html
    """

    logging.info(f"Merging SVs in {vcf_file_path} using {mode} with eps={eps} and min_samples={min_samples}...")

    # # Read VCF file into a pandas DataFrame
    # vcf_df = pd.read_csv(vcf_file_path, sep='\t', comment='#', header=None, \
    #                      names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'], \
    #                         dtype={'CHROM': str, 'POS': np.int64, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, \
    #                                'FILTER': str, 'INFO': str, 'FORMAT': str, 'SAMPLE': str})
    
    # Read VCF file into a pandas DataFrame, using only CHROM, POS, and INFO
    # columns
    logging.info("Reading VCF file into a pandas DataFrame...")
    vcf_df = pd.read_csv(vcf_file_path, sep='\t', comment='#', header=None, usecols=[0, 1, 7], \
                            names=['CHROM', 'POS', 'INFO'], \
                                dtype={'CHROM': str, 'POS': np.int64, 'INFO': str})
    logging.info("VCF file read into a pandas DataFrame.")

    # Print total number of records
    logging.info(f"Total number of records: {vcf_df.shape[0]}")

    # Store a list of record indices that will form the merged VCF file
    merged_indices = []
    merge_records = []
    
    # Create a set with each chromosome in the VCF file
    chromosomes = set(vcf_df['CHROM'].values)

    # Number of clustering plots to generate
    max_plots = 10
    num_plots = 0

    # # Open a new VCF file for writing
    # merged_vcf = os.path.splitext(vcf_file_path)[0] + suffix + '.vcf'
    # with open(merged_vcf, 'w', encoding='utf-8') as merged_vcf_file:

        # # Write the VCF header to the merged VCF file
        # with open(vcf_file_path, 'r', encoding='utf-8') as vcf_file:
        #     for line in vcf_file:
        #         if line.startswith('#'):
        #             merged_vcf_file.write(line)
        #         else:
        #             break

    # Iterate over each chromosome
    first_del_record = None
    first_ins_dup_record = None
    records_processed = 0
    for chromosome in chromosomes:
        logging.info(f"Clustering SVs on chromosome {chromosome}...")

        # Get the chromosome deletion, insertion, and duplication breakpoints
        chr_del_df = vcf_df[(vcf_df['CHROM'] == chromosome) & (vcf_df['INFO'].str.contains('SVTYPE=DEL'))]
        chr_ins_dup_df = vcf_df[(vcf_df['CHROM'] == chromosome) & ((vcf_df['INFO'].str.contains('SVTYPE=INS')) | (vcf_df['INFO'].str.contains('SVTYPE=DUP')))]

        # Get the deletion start and end positions
        chr_del_start = chr_del_df['POS'].values
        chr_del_end = chr_del_df['INFO'].str.extract(r'END=(\d+)', expand=False).astype(np.int32)

        # Get the insertion and duplication start and end positions
        chr_ins_dup_start = chr_ins_dup_df['POS'].values
        chr_ins_dup_len = chr_ins_dup_df['INFO'].str.extract(r'SVLEN=(-?\d+)', expand=False).astype(np.int32)
        chr_ins_dup_end = chr_ins_dup_start + chr_ins_dup_len - 1

        # Format the deletion breakpoints
        chr_del_breakpoints = np.column_stack((chr_del_start, chr_del_end))
        logging.info("Number of deletion breakpoints: " + str(len(chr_del_breakpoints)))

        # Format the insertion and duplication breakpoints
        chr_ins_dup_breakpoints = np.column_stack((chr_ins_dup_start, chr_ins_dup_end))
        logging.info("Number of insertion and duplication breakpoints: " + str(len(chr_ins_dup_breakpoints)))

        # Get the SV depth and clipped base evidence scores for deletions
        chr_del_depth_scores = chr_del_df['INFO'].str.extract(r'DP=(\d+)', expand=False).astype(np.int32)
        chr_del_clipped_bases = chr_del_df['INFO'].str.extract(r'CLIPDP=(\d+)', expand=False).astype(np.int32)
        chr_del_depth_scores = chr_del_depth_scores + chr_del_clipped_bases

        # Get the SV depth and clipped base evidence scores for insertions
        # and duplications
        chr_ins_dup_depth_scores = chr_ins_dup_df['INFO'].str.extract(r'DP=(\d+)', expand=False).astype(np.int32)
        chr_ins_dup_clipped_bases = chr_ins_dup_df['INFO'].str.extract(r'CLIPDP=(\d+)', expand=False).astype(np.int32)
        chr_ins_dup_depth_scores = chr_ins_dup_depth_scores + chr_ins_dup_clipped_bases

        # Cluster SV breakpoints using the specified mode
        if mode == 'dbscan':
            # Cluster SV breakpoints using DBSCAN
            dbscan = DBSCAN(eps=eps, min_samples=min_samples)

            # Cluster deletion breakpoints
            if len(chr_del_breakpoints) > 0:
                logging.info(f"Clustering deletion breakpoints using DBSCAN with eps={eps} and min_samples={min_samples}...")
                del_labels = dbscan.fit_predict(chr_del_breakpoints)
                logging.info(f"Deletion label counts: {len(np.unique(del_labels))}")
            else:
                del_labels = []

            # Cluster insertion and duplication breakpoints together since
            # insertions are a subset of duplications
            if len(chr_ins_dup_breakpoints) > 0:
                logging.info(f"Clustering insertion and duplication breakpoints using DBSCAN with eps={eps} and min_samples={min_samples}...")
                ins_dup_labels = dbscan.fit_predict(chr_ins_dup_breakpoints)
                logging.info(f"Insertion and duplication label counts: {len(np.unique(ins_dup_labels))}")
            else:
                ins_dup_labels = []

        elif mode == 'agglomerative':
            # Cluster SV breakpoints using agglomerative clustering
            logging.info(f"Clustering deletion breakpoints using agglomerative clustering with eps={eps}...")
            agglomerative = AgglomerativeClustering(n_clusters=None, distance_threshold=eps, compute_full_tree=True)

            # Cluster deletion breakpoints
            del_labels = agglomerative.fit_predict(chr_del_breakpoints)
            logging.info(f"Deletion label counts: {len(np.unique(del_labels))}")

            # Cluster insertion breakpoints
            logging.info(f"Clustering insertion and duplication breakpoints using agglomerative clustering with eps={eps}...")
            ins_labels = agglomerative.fit_predict(chr_ins_dup_breakpoints)
            logging.info(f"Insertion label counts: {len(np.unique(ins_labels))}")

        # Get the unique deletion labels for the chromosome
        unique_del_labels = np.unique(del_labels)

        # Get the unique insertion and duplication labels for the chromosome
        unique_ins_dup_labels = np.unique(ins_dup_labels)

        # Merge deletions with the same label
        del_count = 0
        for label in unique_del_labels:

            # Skip label -1 (outliers)
            if label == -1:
                continue

            # Get the indices of SVs with the same label
            idx = del_labels == label

            # logging.info the number of SVs with the same label
            #logging.info(f"Chromosome {chromosome} - {len(chr_del_breakpoints[idx])} deletions with label {label}")

            # Get the SV depth scores with the same label
            depth_scores = chr_del_depth_scores[idx]

            # Get the index of the SV with the highest depth score
            max_depth_score_idx = np.argmax(depth_scores)

            # Plot the DBSCAN clustering behavior if there are 10 < X < 20 SVs with the same label
            plot_enabled = False
            if plot_enabled:
                if len(chr_del_breakpoints[idx]) > 10 and len(chr_del_breakpoints[idx]) < 20 and num_plots < max_plots:

                    # Increment the number of plots
                    num_plots += 1

                    # Convert the max depth score index (index within labels) to the index within the original deletion DataFrame
                    chosen_idx = np.where(idx)[0][max_depth_score_idx]
                    chosen_breakpoints = chr_del_breakpoints[chosen_idx]
                    plot_dbscan(chr_del_breakpoints[idx], chosen_breakpoints, filename=f"dbscan_clustering_{num_plots}.png")

                    # TEST: Return if the number of plots is reached
                    if num_plots == max_plots:
                        return

            # Get the index of the SV with the highest depth score (index within
            # full VCF file), need to convert DEL label -> DEL index -> VCF
            # index
            best_del_vcf_idx = chr_del_df.index[idx][max_depth_score_idx]

            # Add the index of the SV with the highest depth score to the list
            # of indices that will form the merged VCF file
            merged_indices.append(best_del_vcf_idx)

            # Append the chosen record to the list of records that will form the
            # merged VCF file
            merge_records.append(chr_del_df.iloc[idx, :].iloc[max_depth_score_idx, :])

            # chosen_idx = np.where(idx)[0][max_depth_score_idx]
            # Get the VCF record with the highest depth score
            vcf_record = chr_del_df.iloc[idx, :].iloc[max_depth_score_idx, :]

            # Update the first deletion record
            if first_del_record is None:
                first_del_record = vcf_record

            # Write the full VCF record to the merged VCF file
            # merged_vcf_file.write(f"{vcf_record['CHROM']}\t{vcf_record['POS']}\t{vcf_record['ID']}\t{vcf_record['REF']}\t{vcf_record['ALT']}\t{vcf_record['QUAL']}\t{vcf_record['FILTER']}\t{vcf_record['INFO']}\t{vcf_record['FORMAT']}\t{vcf_record['SAMPLE']}\n")

        # Merge insertions and duplications with the same label
        ins_count = 0
        dup_count = 0
        ins_dup_count = 0
        for label in unique_ins_dup_labels:

            # Skip label -1 (outliers)
            if label == -1:
                continue

            # Get the indices of SVs with the same label
            idx = ins_dup_labels == label

            # Get the SV depth scores with the same label
            depth_scores = chr_ins_dup_depth_scores[idx]

            # Get the index of the SV with the highest depth score
            max_depth_score_idx = np.argmax(depth_scores)

            # # Get the VCF record with the highest depth score (This also determines if the SV is an insertion or duplication)
            # vcf_record = chr_ins_dup_df.iloc[idx, :].iloc[max_depth_score_idx, :]

            # Get the VCF record with the highest depth score
            max_ins_dup_record = chr_ins_dup_df.iloc[idx, :].iloc[max_depth_score_idx, :]

            # Update the first insertion or duplication record
            if first_ins_dup_record is None:
                first_ins_dup_record = max_ins_dup_record

            # Get the index for the original VCF file
            max_ins_dup_index = chr_ins_dup_df.index[idx][max_depth_score_idx]

            # Print the two records to see if they are the same
            # logging.info(f"max_ins_dup_record: {max_ins_dup_record}")
            # logging.info(f"max_ins_dup_index record: {vcf_df.iloc[max_ins_dup_index, :]}")

            # Get the index of the SV with the highest depth score
            max_depth_score_idx = np.argmax(depth_scores)

            # Get the index in the original VCF file
            best_ins_dup_vcf_idx = chr_ins_dup_df.index[idx][max_depth_score_idx]

            # Add the index of the SV with the highest depth score to the list
            # of indices that will form the merged VCF file
            merged_indices.append(best_ins_dup_vcf_idx)

            # Append the chosen record to the list of records that will form the
            # merged VCF file
            merge_records.append(chr_ins_dup_df.iloc[idx, :].iloc[max_depth_score_idx, :])

            # Update SV type count
            # if 'SVTYPE=INS' in vcf_record['INFO']:
            #     ins_count += 1
            # elif 'SVTYPE=DUP' in vcf_record['INFO']:
            #     dup_count += 1

            # Write the full VCF record to the merged VCF file
            # merged_vcf_file.write(f"{vcf_record['CHROM']}\t{vcf_record['POS']}\t{vcf_record['ID']}\t{vcf_record['REF']}\t{vcf_record['ALT']}\t{vcf_record['QUAL']}\t{vcf_record['FILTER']}\t{vcf_record['INFO']}\t{vcf_record['FORMAT']}\t{vcf_record['SAMPLE']}\n")

        logging.info(f"Chromosome {chromosome} - {del_count} deletions, {ins_dup_count} insertions, and duplications merged.")

        records_processed += chr_del_breakpoints.shape[0] + chr_ins_dup_breakpoints.shape[0]

        # logging.info(f"Chromosome {chromosome} - {del_count} deletions,
        # {ins_count} insertions, and {dup_count} duplications merged.")
        
    # Free up memory
    del chr_del_df
    del chr_ins_dup_df
    del chr_del_start
    del chr_del_end
    del chr_ins_dup_start
    del chr_ins_dup_len
    del chr_ins_dup_end
    del chr_del_breakpoints
    del chr_ins_dup_breakpoints
    del chr_del_depth_scores
    del chr_ins_dup_depth_scores
    del del_labels
    del ins_dup_labels
    del unique_del_labels
    del unique_ins_dup_labels

    # logging.info("Saved merged VCF file to " + merged_vcf)
    logging.info(f"Processed {records_processed} records of {vcf_df.shape[0]} total records.")

    # Open a new VCF file for writing
    logging.info("Writing merged VCF file...")
    merged_vcf = os.path.splitext(vcf_file_path)[0] + suffix + '.vcf'

    logging.info(f"Writing {len(merged_indices)} records to merged VCF file...")

    # # Read the specific rows from the original VCF file
    # vcf_df = pd.read_csv(vcf_file_path, sep='\t', comment='#', header=None, \
    #                      names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'], \
    #                         dtype={'CHROM': str, 'POS': np.int64, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, \
    #                                'FILTER': str, 'INFO': str, 'FORMAT': str, 'SAMPLE': str}, \
    #                         skiprows=lambda x: x not in merged_indices)

    # Read the entire VCF file into a DataFrame
    logging.info("Reading the entire VCF file into a DataFrame...")
    vcf_df = pd.read_csv(vcf_file_path, sep='\t', comment='#', header=None, \
                         names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'], \
                            dtype={'CHROM': str, 'POS': np.int64, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, \
                                   'FILTER': str, 'INFO': str, 'FORMAT': str, 'SAMPLE': str})

    with open(merged_vcf, 'w', encoding='utf-8') as merged_vcf_file:

        # Write the VCF header to the merged VCF file
        with open(vcf_file_path, 'r', encoding='utf-8') as vcf_file:
            for line in vcf_file:
                if line.startswith('#'):
                    merged_vcf_file.write(line)
                else:
                    break

        # Print number of records written
        logging.info(f"Writing {len(merge_records)} records to merged VCF file...")

        # Write the specific rows from the original VCF file to the merged VCF
        # file
        merge_count = 0
        for merged_record in merge_records:
            # Find the matching record in the original VCF file
            matching_record = vcf_df[(vcf_df['CHROM'] == merged_record['CHROM']) & (vcf_df['POS'] == merged_record['POS']) & (vcf_df['INFO'] == merged_record['INFO'])]

            # Write the matching record to the merged VCF file
            merge_count += 1
            merged_vcf_file.write(f"{matching_record['CHROM'].values[0]}\t{matching_record['POS'].values[0]}\t{matching_record['ID'].values[0]}\t{matching_record['REF'].values[0]}\t{matching_record['ALT'].values[0]}\t{matching_record['QUAL'].values[0]}\t{matching_record['FILTER'].values[0]}\t{matching_record['INFO'].values[0]}\t{matching_record['FORMAT'].values[0]}\t{matching_record['SAMPLE'].values[0]}\n")
            
            # Log every 1000 records written
            logging.info(f"Wrote {merge_count} records to merged VCF file...")

        # # Write the merged VCF records to the merged VCF file
        # printed_first_del_record = False
        # printed_first_ins_dup_record = False
        # for i in range(vcf_df.shape[0]):
        #     vcf_record = vcf_df.iloc[i, :]

        #     # Compare the first deletion or insertion/duplication record to the first record in the merged VCF file
        #     if first_del_record is not None and first_ins_dup_record is not None:
        #         if vcf_record['INFO'].find('SVTYPE=DEL') != -1 and not printed_first_del_record:
        #             # Print the first deletion record
        #             logging.info(f"First deletion record: {first_del_record}")
        #             logging.info(f"First deletion merged record: {vcf_record}")
        #             printed_first_del_record = True

        #         elif (vcf_record['INFO'].find('SVTYPE=INS') != -1 or vcf_record['INFO'].find('SVTYPE=DUP') != -1) and not printed_first_ins_dup_record:
        #             # Print the first insertion or duplication record
        #             logging.info(f"First insertion or duplication record: {first_ins_dup_record}")
        #             logging.info(f"First insertion or duplication merged record: {vcf_record}")
        #             printed_first_ins_dup_record = True

        #     merged_vcf_file.write(f"{vcf_record['CHROM']}\t{vcf_record['POS']}\t{vcf_record['ID']}\t{vcf_record['REF']}\t{vcf_record['ALT']}\t{vcf_record['QUAL']}\t{vcf_record['FILTER']}\t{vcf_record['INFO']}\t{vcf_record['FORMAT']}\t{vcf_record['SAMPLE']}\n")

    logging.info("Merged VCF file written to " + merged_vcf)

    return merged_vcf

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        logging.info(f"Usage: {sys.argv[0]} <VCF file path>")
        sys.exit(1)

    # Get the VCF file path from the command line
    vcf_file_path = sys.argv[1]

    # Get the epsilon value from the command line
    if len(sys.argv) > 2:
        eps = int(sys.argv[2])
    else:
        eps = 30

    # DBSCAN 
    sv_merger(sys.argv[1], mode='dbscan', eps=eps, suffix='.merged_eps' + str(eps))
