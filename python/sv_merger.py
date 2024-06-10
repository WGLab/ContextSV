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
            # Get the current SUPPORT field
            previous_support = int(field.split('=')[1])

            # Add the cluster size to the SUPPORT field
            updated_info += f'SUPPORT={previous_support + cluster_size};'
            # updated_info += f'SUPPORT={cluster_size};'
        else:
            updated_info += field + ';'  # Append the field to the updated INFO
            
    # Update the INFO column
    record['INFO'] = updated_info

    return record

def cluster_breakpoints(vcf_df, sv_type, cluster_size_min):
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

    elif sv_type == 'INS/DUP':
        sv_start = vcf_df['POS'].values
        sv_len = vcf_df['INFO'].str.extract(r'SVLEN=(-?\d+)', expand=False).astype(np.int32)
        sv_end = sv_start + sv_len - 1

        # Format the insertion and duplication breakpoints
        breakpoints = np.column_stack((sv_start, sv_end))
    else:
        logging.error("Invalid SV type: %s", sv_type)
        return
    
    # Get the combined SV read and clipped base support
    sv_support = vcf_df['INFO'].str.extract(r'SUPPORT=(\d+)', expand=False).astype(np.int32)
    sv_clipped_base_support = vcf_df['INFO'].str.extract(r'CLIPSUP=(\d+)', expand=False).astype(np.int32)
    sv_support = sv_support + sv_clipped_base_support

    # Get the HMM likelihood scores
    hmm_scores = vcf_df['INFO'].str.extract(r'HMM=(\d+)', expand=False).astype(np.float32)

    # Set all 0 values to NaN
    hmm_scores[hmm_scores == 0] = np.nan

    # Cluster SV breakpoints using HDBSCAN
    cluster_labels = []
    dbscan = HDBSCAN(min_cluster_size=cluster_size_min, min_samples=3)
    if len(breakpoints) > 0:
        logging.info("Clustering %d SV breakpoints...", len(breakpoints))
        cluster_labels = dbscan.fit_predict(breakpoints)

        logging.info("Label counts: %d", len(np.unique(cluster_labels)))

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
        max_hmm_score = None
        max_support = None
        if len(np.unique(cluster_hmm_scores)) > 1:
            # Find the index of the maximum score, excluding NaN values
            max_score_idx = np.nanargmax(cluster_hmm_scores)
            # max_score_idx = np.argmax(cluster_hmm_scores)
            max_hmm_score = cluster_hmm_scores[max_score_idx]

        # Use the SV with the highest read support if the log likelihood
        # scores are all the same
        elif len(np.unique(cluster_depth_scores)) > 1:
            max_score_idx = np.argmax(cluster_depth_scores)
            max_support = cluster_depth_scores[max_score_idx]

        # Use the first SV in the cluster if the depth scores are all the
        # same
        else:
            max_score_idx = 0

        # Get the VCF record with the highest depth score
        max_record = vcf_df.iloc[idx, :].iloc[max_score_idx, :]

        # Get the number of SVs in this cluster
        cluster_size = np.sum(idx)
        # logging.info("DEL Cluster size: %s", cluster_size)

        # Update the SUPPORT field in the INFO column
        max_record = update_support(max_record, cluster_size)

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
        merged_records.loc[merged_records.shape[0]] = max_record

    return merged_records

def sv_merger(vcf_file_path, cluster_size_min=3, suffix='.merged'):
    """
    Use DBSCAN to merge SVs with the same breakpoint.
    Mode can be 'dbscan', 'gmm', or 'agglomerative'.
    https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html
    """

    logging.info("Merging SVs in %s using HDBSCAN with minimum cluster size=%d...", vcf_file_path, cluster_size_min)
    
    # Read VCF file into a pandas DataFrame, using only CHROM, POS, and INFO
    # columns
    logging.info("Reading VCF file into a pandas DataFrame...")
    vcf_df = pd.read_csv(vcf_file_path, sep='\t', comment='#', header=None, usecols=[0, 1, 7], \
                            names=['CHROM', 'POS', 'INFO'], \
                                dtype={'CHROM': str, 'POS': np.int64, 'INFO': str})

    # Add a column at the beginning with the index
    vcf_df.insert(0, 'INDEX', range(0, len(vcf_df)))
    logging.info("Reading complete.")

    # Print total number of records
    logging.info("Total number of records: %d", vcf_df.shape[0])

    # Store a dataframe of records that will form the merged VCF file
    merged_records = pd.DataFrame(columns=['INDEX', 'CHROM', 'POS', 'INFO'])

    # Create a set with each chromosome in the VCF file
    chromosomes = set(vcf_df['CHROM'].values)

    # Iterate over each chromosome
    records_processed = 0
    current_chromosome = 0
    chromosome_count = len(chromosomes)
    for chromosome in chromosomes:

        # Cluster deletions
        logging.info("Clustering deletions on chromosome %s...", chromosome)
        chr_del_df = vcf_df[(vcf_df['CHROM'] == chromosome) & (vcf_df['INFO'].str.contains('SVTYPE=DEL'))]
        del_records = cluster_breakpoints(chr_del_df, 'DEL', cluster_size_min)

        # Cluster insertions and duplications
        logging.info("Clustering insertions and duplications on chromosome %s...", chromosome)
        chr_ins_dup_df = vcf_df[(vcf_df['CHROM'] == chromosome) & ((vcf_df['INFO'].str.contains('SVTYPE=INS')) | (vcf_df['INFO'].str.contains('SVTYPE=DUP')))]
        ins_dup_records = cluster_breakpoints(chr_ins_dup_df, 'INS/DUP', cluster_size_min)

        # Summarize the number of deletions and insertions/duplications
        del_count = del_records.shape[0]
        ins_dup_count = ins_dup_records.shape[0]
        records_processed += del_count + ins_dup_count
        logging.info("Chromosome %s - %d deletions, %d insertions, and duplications merged.", chromosome, del_count, ins_dup_count)

        # Append the deletion and insertion/duplication records to the merged
        # records DataFrame
        merged_records = pd.concat([merged_records, del_records, ins_dup_records], ignore_index=True)

        current_chromosome += 1
        logging.info("Processed %d of %d chromosomes.", current_chromosome, chromosome_count)
        
    logging.info("Processed %d records of %d total records.", records_processed, vcf_df.shape[0])


    # Free up memory
    del vcf_df
    del chr_del_df
    del chr_ins_dup_df

    # Open a new VCF file for writing
    logging.info("Writing merged VCF file...")
    merged_vcf = os.path.splitext(vcf_file_path)[0] + suffix + '.vcf'

    logging.info("Writing %d records to merged VCF file...", merged_records.shape[0])

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

            logging.info("Wrote %d of %d records to merged VCF file...", merge_count, merged_records.shape[0])

    logging.info("Merged VCF file written to %s", merged_vcf)

    return merged_vcf

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        logging.info("Usage: %s <VCF file path>", sys.argv[0])
        sys.exit(1)

    # Get the VCF file path from the command line
    vcf_file_path = sys.argv[1]

    # Check if the file exists
    if not os.path.exists(vcf_file_path):
        logging.error("Error: %s not found.", vcf_file_path)
        sys.exit(1)

    # Get the minimum cluster size from the command line
    if len(sys.argv) > 2:
        cluster_size_min = int(sys.argv[2])
    else:
        cluster_size_min = 2

    # Get the suffix from the command line
    suffix = '.merged'
    if len(sys.argv) > 3:
        suffix += sys.argv[3]

    # DBSCAN 
    sv_merger(vcf_file_path, cluster_size_min=cluster_size_min, suffix=suffix)
    