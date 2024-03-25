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
    # merge_records = []
    # Store a dataframe of records that will form the merged VCF file
    merged_records = pd.DataFrame(columns=['CHROM', 'POS', 'INFO'])
    
    # Create a set with each chromosome in the VCF file
    chromosomes = set(vcf_df['CHROM'].values)

    # Number of clustering plots to generate
    max_plots = 10
    num_plots = 0

    # Iterate over each chromosome
    records_processed = 0
    current_chromosome = 0
    chromosome_count = len(chromosomes)
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

            # Get the SV depth scores with the same label
            depth_scores = chr_del_depth_scores[idx]

            # Get the index of the SV with the highest depth score
            max_depth_score_idx = np.argmax(depth_scores)

            # Get the VCF record with the highest depth score
            max_del_record = chr_del_df.iloc[idx, :].iloc[max_depth_score_idx, :]

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

            # Append the chosen record to the list of records that will form the
            # merged VCF file
            # merge_records.append(chr_del_df.iloc[idx, :].iloc[max_depth_score_idx, :])

            # Append the chosen record to the dataframe of records that will
            # form the merged VCF file
            merged_records.loc[merged_records.shape[0]] = max_del_record

        # Merge insertions and duplications with the same label
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

            # Get the VCF record with the highest depth score
            max_ins_dup_record = chr_ins_dup_df.iloc[idx, :].iloc[max_depth_score_idx, :]

            # Append the chosen record to the list of records that will form the
            # merged VCF file
            # merge_records.append(chr_ins_dup_df.iloc[idx, :].iloc[max_depth_score_idx, :])

            # Append the chosen record to the dataframe of records that will
            # form the merged VCF file
            merged_records.loc[merged_records.shape[0]] = max_ins_dup_record

        logging.info(f"Chromosome {chromosome} - {del_count} deletions, {ins_dup_count} insertions, and duplications merged.")
        
        current_chromosome += 1
        logging.info(f"Processed {current_chromosome} of {chromosome_count} chromosomes.")

        records_processed += chr_del_breakpoints.shape[0] + chr_ins_dup_breakpoints.shape[0]

        # logging.info(f"Chromosome {chromosome} - {del_count} deletions,
        # {ins_count} insertions, and {dup_count} duplications merged.")
        
    # logging.info("Saved merged VCF file to " + merged_vcf)
    logging.info(f"Processed {records_processed} records of {vcf_df.shape[0]} total records.")

    # Free up memory
    del vcf_df
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

    # Open a new VCF file for writing
    logging.info("Writing merged VCF file...")
    merged_vcf = os.path.splitext(vcf_file_path)[0] + suffix + '.vcf'

    logging.info(f"Writing {merged_records.shape[0]} records to merged VCF file...")

    merge_count = 0
    with open(merged_vcf, 'w', encoding='utf-8') as merged_vcf_file:

        # Write the VCF header to the merged VCF file
        with open(vcf_file_path, 'r', encoding='utf-8') as vcf_file:
            for line in vcf_file:
                if line.startswith('#'):
                    merged_vcf_file.write(line)
                else:
                    break

        # Read the next 1000 records from the original VCF file
        logging.info("Reading a chunk of 1000 records from the original VCF file...")
        for chunk in pd.read_csv(vcf_file_path, sep='\t', comment='#', header=None, \
                                    names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'], \
                                    dtype={'CHROM': str, 'POS': np.int64, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, \
                                             'FILTER': str, 'INFO': str, 'FORMAT': str, 'SAMPLE': str}, \
                                             chunksize=1000):
            
            # Get the matching records from the chunk by merging on CHROM, POS,
            # and INFO, but only keep the records from the chunk since they
            # contain the full VCF record
            matching_records = pd.merge(chunk, merged_records, on=['CHROM', 'POS', 'INFO'], how='inner')
            matching_records = matching_records.drop_duplicates(subset=['CHROM', 'POS', 'INFO'])  # Drop duplicate records

            # Remove the matching records from the merged records dataframe
            merged_records = merged_records[~merged_records.isin(matching_records)].dropna()

            # Write the matching records to the merged VCF file
            for _, matching_record in matching_records.iterrows():
                merge_count += 1
                merged_vcf_file.write(f"{matching_record['CHROM']}\t{matching_record['POS']}\t{matching_record['ID']}\t{matching_record['REF']}\t{matching_record['ALT']}\t{matching_record['QUAL']}\t{matching_record['FILTER']}\t{matching_record['INFO']}\t{matching_record['FORMAT']}\t{matching_record['SAMPLE']}\n")

            logging.info(f"Wrote {merge_count} of {merged_records.shape[0]} records to merged VCF file...")

    logging.info(f"Processed {merge_count} records of {merged_records.shape[0]} total records.")

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

    # Get the suffix from the command line
    suffix = '.merged_eps' + str(eps)
    if len(sys.argv) > 3:
        suffix += sys.argv[3]

    # DBSCAN 
    sv_merger(sys.argv[1], mode='dbscan', eps=eps, suffix=suffix)
    