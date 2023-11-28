import os
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN


def merge_svs_by_breakpoint(vcf_file_path, eps=1000, min_samples=2, suffix='.merged'):
    """
    Use DBSCAN to merge SVs with the same breakpoint.
    https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html
    """
    # Read VCF file into a pandas DataFrame
    vcf_df = pd.read_csv(vcf_file_path, sep='\t', comment='#', header=None, \
                         names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'], \
                            dtype={'CHROM': str, 'POS': np.int64, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, \
                                   'FILTER': str, 'INFO': str, 'FORMAT': str, 'SAMPLE': str})
    
    # vcf_df = pd.read_csv(vcf_file_path, sep='\t', comment='#', header=None, \
    #                      names=[vcf_headers], dtype={'CHROM': str, 'POS': np.int32})

    # Set up arrays for storing the indices for true SV breakpoints from the
    # VCF file
    del_idx = vcf_df['INFO'].str.contains('SVTYPE=DEL')

    # Create a set with each chromosome in the VCF file
    chromosomes = set(vcf_df['CHROM'].values)

    # Create a dictionary with each chromosome as a key and a list of SV
    # breakpoints as the value. This will be used to cluster SVs by chromosome.
    # chr_del_breakpoints = {}
    # chr_ins_breakpoints = {}

    # Create arrays containing the read depth + clipped base evidence for each
    # SV (INFO/DP + INFO/CLIPDP)
    sv_read_depth = vcf_df['INFO'].str.extract(r'DP=(\d+)', expand=False).astype(np.int32)
    sv_clipped_bases = vcf_df['INFO'].str.extract(r'CLIPDP=(\d+)', expand=False).astype(np.int32)
    sv_depth_scores = sv_read_depth + sv_clipped_bases

    # # Create dictionaries with each chromosome as a key and a list of SV depth scores
    # chr_del_depth_scores = {}
    # chr_ins_depth_scores = {}

    # Open a new VCF file for writing
    merged_vcf = os.path.splitext(vcf_file_path)[0] + suffix + '.vcf'
    with open(merged_vcf, 'w', encoding='utf-8') as merged_vcf_file:

        # Write the VCF header to the merged VCF file
        with open(vcf_file_path, 'r', encoding='utf-8') as vcf_file:
            for line in vcf_file:
                if line.startswith('#'):
                    merged_vcf_file.write(line)
                else:
                    break

        # Iterate over each chromosome
        for chromosome in chromosomes:
            # Get the chromosome deletion and insertion records
            chr_del_df = vcf_df[(vcf_df['CHROM'] == chromosome) & (vcf_df['INFO'].str.contains('SVTYPE=DEL'))]
            chr_ins_df = vcf_df[(vcf_df['CHROM'] == chromosome) & (vcf_df['INFO'].str.contains('SVTYPE=INS'))]

            # # Further filter insertions to duplications only (INFO/REPTYPE=DUP)
            # chr_ins_df = chr_ins_df[chr_ins_df['INFO'].str.contains('REPTYPE=DUP')]

            # Get the deletion start and end positions
            chr_del_start = chr_del_df['POS'].values
            chr_del_end = chr_del_df['INFO'].str.extract(r'END=(\d+)', expand=False).astype(np.int32)

            # Get the insertion start and end positions (end = start + length)
            chr_ins_start = chr_ins_df['POS'].values
            chr_ins_len = chr_ins_df['INFO'].str.extract(r'SVLEN=(-?\d+)', expand=False).astype(np.int32)
            chr_ins_end = chr_ins_start + chr_ins_len - 1

            # Format the deletion breakpoints
            chr_del_breakpoints = np.column_stack((chr_del_start, chr_del_end))

            # Format the insertion breakpoints
            chr_ins_breakpoints = np.column_stack((chr_ins_start, chr_ins_end))

            # Get the SV depth and clipped base evidence scores for deletions
            chr_del_depth_scores = chr_del_df['INFO'].str.extract(r'DP=(\d+)', expand=False).astype(np.int32)
            chr_del_clipped_bases = chr_del_df['INFO'].str.extract(r'CLIPDP=(\d+)', expand=False).astype(np.int32)
            chr_del_depth_scores = chr_del_depth_scores + chr_del_clipped_bases

            # Get the SV depth and clipped base evidence scores for insertions
            chr_ins_depth_scores = chr_ins_df['INFO'].str.extract(r'DP=(\d+)', expand=False).astype(np.int32)
            chr_ins_clipped_bases = chr_ins_df['INFO'].str.extract(r'CLIPDP=(\d+)', expand=False).astype(np.int32)
            chr_ins_depth_scores = chr_ins_depth_scores + chr_ins_clipped_bases

            # Cluster SV breakpoints using DBSCAN
            dbscan = DBSCAN(eps=eps, min_samples=min_samples)
            
            # Cluster deletion breakpoints
            del_labels = dbscan.fit_predict(chr_del_breakpoints)

            # Cluster insertion breakpoints
            ins_labels = dbscan.fit_predict(chr_ins_breakpoints)

            # Get the unique deletion labels for the chromosome
            unique_del_labels = np.unique(del_labels)

            # Get the unique insertion labels for the chromosome
            unique_ins_labels = np.unique(ins_labels)

            # Merge deletions with the same label
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
                vcf_record = chr_del_df.iloc[idx, :].iloc[max_depth_score_idx, :]

                # Write the full VCF record to the merged VCF file
                merged_vcf_file.write(f"{vcf_record['CHROM']}\t{vcf_record['POS']}\t{vcf_record['ID']}\t{vcf_record['REF']}\t{vcf_record['ALT']}\t{vcf_record['QUAL']}\t{vcf_record['FILTER']}\t{vcf_record['INFO']}\t{vcf_record['FORMAT']}\t{vcf_record['SAMPLE']}\n")

            # Merge insertions with the same label
            for label in unique_ins_labels:

                # Skip label -1 (outliers)
                if label == -1:
                    continue

                # Get the indices of SVs with the same label
                idx = ins_labels == label

                # Get the SV depth scores with the same label
                depth_scores = chr_ins_depth_scores[idx]

                # Get the index of the SV with the highest depth score
                max_depth_score_idx = np.argmax(depth_scores)

                # Get the VCF record with the highest depth score
                vcf_record = chr_ins_df.iloc[idx, :].iloc[max_depth_score_idx, :]

                # Write the full VCF record to the merged VCF file
                merged_vcf_file.write(f"{vcf_record['CHROM']}\t{vcf_record['POS']}\t{vcf_record['ID']}\t{vcf_record['REF']}\t{vcf_record['ALT']}\t{vcf_record['QUAL']}\t{vcf_record['FILTER']}\t{vcf_record['INFO']}\t{vcf_record['FORMAT']}\t{vcf_record['SAMPLE']}\n")


    return merged_vcf


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <VCF file path>")
        sys.exit(1)

    # Testing approach: Find optimal eps value for DBSCAN, then test different
    # min_samples values

    # # Test different eps values from 1000 to 100 by 10
    # for eps in range(1000, 90, -10):
    #     merge_svs_by_breakpoint(sys.argv[1], eps=eps, min_samples=2, suffix=f'.merged_eps{eps}_min2')

    # Test different eps values from 1 to 99 by 1
    for eps in range(1, 100):
        merge_svs_by_breakpoint(sys.argv[1], eps=eps, min_samples=2, suffix=f'.merged_eps{eps}_min2')

    # # Test different min_samples values from 2 to 10
    # for min_samples in range(2, 11):
    #     merge_svs_by_breakpoint(sys.argv[1], eps=100, min_samples=min_samples, suffix=f'.merged_eps100_min{min_samples}')
