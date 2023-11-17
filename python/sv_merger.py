import os
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN


def merge_svs_by_breakpoint(vcf_file_path, eps=1000, min_samples=2):
    """
    Use DBSCAN to merge SVs with the same breakpoint.
    https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html
    """
    # Read VCF file into a pandas DataFrame
    vcf_headers = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']
    vcf_df = pd.read_csv(vcf_file_path, sep='\t', comment='#', header=None, \
                         names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'], \
                            dtype={'CHROM': str, 'POS': np.int64, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, \
                                   'FILTER': str, 'INFO': str, 'FORMAT': str, 'SAMPLE': str})
    
    # vcf_df = pd.read_csv(vcf_file_path, sep='\t', comment='#', header=None, \
    #                      names=[vcf_headers], dtype={'CHROM': str, 'POS': np.int32})

    # Extract SV breakpoints from VCF DataFrame (CHROM, POS, INFO/END)
    sv_breakpoints = vcf_df[['CHROM', 'POS']].values
    sv_end = vcf_df['INFO'].str.extract(r'END=(\d+)', expand=False).astype(np.int32)
    sv_breakpoints = np.column_stack((sv_breakpoints, sv_end))

    # Cluster SV breakpoints using DBSCAN
    dbscan = DBSCAN(eps=eps, min_samples=min_samples)
    sv_labels = dbscan.fit_predict(sv_breakpoints)

    # Merge SVs with the same label
    merged_svs = {}
    for i, label in enumerate(sv_labels):
        if label not in merged_svs:
            merged_svs[label] = {'CHROM': sv_breakpoints[i, 0], 'POS': sv_breakpoints[i, 1]}
        else:
            merged_svs[label]['POS'] = (merged_svs[label]['POS'] + sv_breakpoints[i, 1]) // 2

    # Convert merged SVs dictionary to a pandas DataFrame
    merged_svs_df = pd.DataFrame.from_dict(merged_svs, orient='index')

    # Sort merged SVs DataFrame by chromosome and position
    merged_svs_df.sort_values(by=['CHROM', 'POS'], inplace=True)

    # Create a new VCF file with the same header as the original VCF file
    vcf_header = ''
    with open(vcf_file_path, 'r', encoding='utf-8') as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                vcf_header += line
            else:
                break

    # Format the merged VCf file path
    merged_vcf = os.path.splitext(vcf_file_path)[0] + '.merged.vcf'

    # Open the new VCF file for writing
    with open(merged_vcf, 'w', encoding='utf-8') as merged_vcf_file:
        # Write the header to the new VCF file
        merged_vcf_file.write(vcf_header)

        # Write the merged SVs to the new VCF file
        for _, row in merged_svs_df.iterrows():
            # Get each column value from the merged SVs DataFrame
            chrom = row['CHROM']
            pos = row['POS']

            merged_vcf_file.write(f"{row['CHROM']}\t{row['POS']}\t.\tN\t<SV>\t.\tPASS\t.\tGT\t1/1\n")

    return merged_svs_df


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <VCF file path>")
        sys.exit(1)

    merge_svs_by_breakpoint(sys.argv[1])
