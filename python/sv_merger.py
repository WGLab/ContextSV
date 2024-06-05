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
# from sklearn.cluster import DBSCAN

# OPTICS clustering algorithm
# from sklearn.cluster import OPTICS

# HDBSCAN clustering algorithm
from sklearn.cluster import HDBSCAN


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


def update_support(record, cluster_size):
    """
    Set the SUPPORT field in the INFO column of a VCF record to the cluster size.
    """
    # Get the INFO column
    info = record['INFO']

    # Parse the INFO columns
    info_fields = info.split(';')

    # Loop and update the SUPPORT field, while creating a new INFO string
    updated_info = ''
    for field in info_fields:
        if field.startswith('SUPPORT='):
            # Set the SUPPORT field to the cluster size
            updated_info += f'SUPPORT={cluster_size};'
        else:
            updated_info += field + ';'  # Append the field to the updated INFO
            
    # Update the INFO column
    record['INFO'] = updated_info

    return record

def cluster_breakpoints(vcf_df, sv_type, min_samples=2):
    """
    Cluster SV breakpoints using HDBSCAN.
    """
    # Set up the output DataFrame
    merged_records = pd.DataFrame(columns=['INDEX', 'CHROM', 'POS', 'INFO'])

    # Format the SV breakpoints
    breakpoints = None
    if sv_type == 'DEL':
        sv_start = vcf_df['POS'].values
        sv_end = vcf_df['INFO'].str.extract(r'END=(\d+)', expand=False).astype(np.int32)

        # Format the deletion breakpoints
        breakpoints = np.column_stack((sv_start, sv_end))

    elif sv_type == 'INS' or sv_type == 'DUP':
        sv_start = vcf_df['POS'].values
        sv_len = vcf_df['INFO'].str.extract(r'SVLEN=(-?\d+)', expand=False).astype(np.int32)
        sv_end = sv_start + sv_len - 1

        # Format the insertion and duplication breakpoints
        breakpoints = np.column_stack((sv_start, sv_end))
    else:
        logging.error(f"Invalid SV type: {sv_type}")
        return
    
    # Get the combined SV read and clipped base support
    sv_support = vcf_df['INFO'].str.extract(r'SUPPORT=(\d+)', expand=False).astype(np.int32)
    sv_clipped_base_support = vcf_df['INFO'].str.extract(r'CLIPSUP=(\d+)', expand=False).astype(np.int32)
    sv_support = sv_support + sv_clipped_base_support

    # Get the HMM likelihood scores
    hmm_scores = vcf_df['INFO'].str.extract(r'HMM=(\d+)', expand=False).astype(np.float32)

    # Cluster SV breakpoints using HDBSCAN
    cluster_labels = []
    dbscan = HDBSCAN(min_cluster_size=min_samples)
    if len(breakpoints) > 0:
        logging.info(f"Clustering deletion breakpoints using HDSCAN with minimum cluster size={min_samples}...")
        cluster_labels = dbscan.fit_predict(breakpoints)

        logging.info(f"Label counts: {len(np.unique(cluster_labels))}")

    # Merge SVs with the same label
    unique_labels = np.unique(cluster_labels)
    for label in unique_labels:

        # Skip label -1 (outliers)
        if label == -1:
            # Print the positions if any are within a certain range
            pos_min = 180915940
            pos_max = 180950356
            idx = cluster_labels == label
            pos_values = breakpoints[idx][:, 0]
            if (np.any(pos_values >= pos_min) and np.any(pos_values <= pos_max)):
                # Print all within range
                pos_within_range = pos_values[(pos_values >= pos_min) & (pos_values <= pos_max)]
                logging.info(f"Outlier deletion positions: {pos_within_range}")
                # logging.info(f"Outlier deletion positions: {pos_values}")

            continue

        # Get the indices of SVs with the same label
        idx = cluster_labels == label

        # Use the SV with the lowest log likelihood score, if available
        # (values are not all the same)
        max_score_idx = None
        cluster_hmm_scores = hmm_scores[idx]
        cluster_depth_scores = sv_support[idx]
        if len(np.unique(cluster_hmm_scores)) > 1:
            max_score_idx = np.argmin(cluster_hmm_scores)

        # Use the SV with the highest read support if the log likelihood
        # scores are all the same
        elif len(np.unique(cluster_depth_scores)) > 1:
            max_score_idx = np.argmax(cluster_depth_scores)

        # Use the first SV in the cluster if the depth scores are all the
        # same
        else:
            max_score_idx = 0

        # Get the VCF record with the highest depth score
        max_del_record = vcf_df.iloc[idx, :].iloc[max_score_idx, :]

        # Get the number of SVs in this cluster
        cluster_size = np.sum(idx)
        # logging.info("DEL Cluster size: %s", cluster_size)

        # Update the SUPPORT field in the INFO column
        max_del_record = update_support(max_del_record, cluster_size)

        # Get all position values in the cluster
        pos_values = breakpoints[idx][:, 0]

        # If the POS value is a certain value, plot the support
        pos_min = 180915940
        pos_max = 180950356
        if (np.any(pos_values >= pos_min) and np.any(pos_values <= pos_max)) or cluster_size > 1000:
            logging.info(f"Cluster size: {cluster_size}")
            logging.info(f"Pos values: {pos_values}")

        # Append the chosen record to the dataframe of records that will
        # form the merged VCF file
        merged_records.loc[merged_records.shape[0]] = max_del_record

    return merged_records


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

    # Add a column at the beginning with the index
    vcf_df.insert(0, 'INDEX', range(0, len(vcf_df)))
    logging.info("VCF file read into a pandas DataFrame.")

    # Print total number of records
    logging.info(f"Total number of records: {vcf_df.shape[0]}")

    # Store a dataframe of records that will form the merged VCF file
    merged_records = pd.DataFrame(columns=['INDEX', 'CHROM', 'POS', 'INFO'])

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

        # Cluster deletions
        logging.info(f"Clustering deletions on chromosome {chromosome}...")
        chr_del_df = vcf_df[(vcf_df['CHROM'] == chromosome) & (vcf_df['INFO'].str.contains('SVTYPE=DEL'))]
        del_records = cluster_breakpoints(chr_del_df, 'DEL', min_samples=min_samples)

        # Cluster insertions and duplications
        logging.info(f"Clustering insertions and duplications on chromosome {chromosome}...")
        chr_ins_dup_df = vcf_df[(vcf_df['CHROM'] == chromosome) & ((vcf_df['INFO'].str.contains('SVTYPE=INS')) | (vcf_df['INFO'].str.contains('SVTYPE=DUP')))]
        ins_dup_records = cluster_breakpoints(chr_ins_dup_df, 'INS', min_samples=min_samples)

        # Summarize the number of deletions and insertions/duplications
        del_count = del_records.shape[0]
        ins_dup_count = ins_dup_records.shape[0]
        records_processed += del_count + ins_dup_count
        logging.info(f"Chromosome {chromosome} - {del_count} deletions, {ins_dup_count} insertions, and duplications merged.")

        # Append the deletion and insertion/duplication records to the merged
        # records DataFrame
        merged_records = pd.concat([merged_records, del_records, ins_dup_records], ignore_index=True)

        # logging.info(f"Chromosome {chromosome} - {del_count} deletions, {ins_dup_count} insertions, and duplications merged.")
       
        current_chromosome += 1
        logging.info(f"Processed {current_chromosome} of {chromosome_count} chromosomes.")
        
    logging.info(f"Processed {records_processed} records of {vcf_df.shape[0]} total records.")


    # Free up memory
    del vcf_df
    del chr_del_df
    del chr_ins_dup_df

    # Open a new VCF file for writing
    logging.info("Writing merged VCF file...")
    merged_vcf = os.path.splitext(vcf_file_path)[0] + suffix + '.vcf'

    logging.info(f"Writing {merged_records.shape[0]} records to merged VCF file...")

    merge_count = 0
    index_start = 0
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
            
            # Add an INDEX column to the chunk
            chunk.insert(0, 'INDEX', range(index_start, index_start + chunk.shape[0]))
            index_start += chunk.shape[0]

            # Merge on INDEX, and use all information from the original VCF file
            # (chunk) but update the INFO field with the merged INFO field.
            # This is done by dropping the INFO column from the chunk so that
            # the INFO column from the merged_records dataframe is used.
            matching_records = pd.merge(chunk.drop(columns=['INFO']), merged_records[['INDEX', 'INFO']], on=['INDEX'], how='inner')
            matching_records = matching_records.drop_duplicates(subset=['INDEX'])  # Drop duplicate records
            matching_records = matching_records.drop(columns=['INDEX'])  # Drop the INDEX column

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

    # Check if the file exists
    if not os.path.exists(vcf_file_path):
        logging.error(f"Error: {vcf_file_path} not found.")
        sys.exit(1)

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
    