"""
extract_features.py: Extract features from the input VCF file.

Usage:
    extract_features.py <input>

Arguments:
    <input>     Path to the input VCF file.

Output:
    A dataframe with a column for each feature.
"""

import os
import sys
import logging
import numpy as np
import pandas as pd


def read_vcf(filepath):
    """Read in the VCF file."""
    vcf_df = pd.read_csv(filepath, sep='\t', comment='#', header=None, usecols=[0, 1, 7], \
                         names=['CHROM', 'POS', 'INFO'], \
                            dtype={'CHROM': str, 'POS': np.int64, 'INFO': str})
    return vcf_df

def extract_features(input_vcf):
    """Extract the features from the VCF file's data."""
    # Read in the VCF file.
    vcf_df = read_vcf(input_vcf)

    # Extract the read and clipped base support from the INFO column.
    read_support = vcf_df['INFO'].str.extract(r'SUPPORT=(\d+)', expand=False).astype(np.int32)

    # Check if any read depths are missing.
    if read_support.isnull().values.any():
        logging.error('Read support is missing.')
        sys.exit(1)

    clipped_bases = vcf_df['INFO'].str.extract(r'CLIPSUP=(\d+)', expand=False).astype(np.int32)

    # Check if any clipped bases are missing.
    if clipped_bases.isnull().values.any():
        logging.error('Clipped bases is missing.')
        sys.exit(1)

    # Get the array of chromosome names.
    chrom = vcf_df['CHROM']

    # Keep only chromosomes 1-22, X, Y, and MT.
    chrom = chrom[chrom.isin(['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY', 'chrMT'])]

    # If empty, try without the 'chr' prefix.
    if chrom.empty:
        chrom = vcf_df['CHROM']
        chrom = chrom[chrom.isin([str(i) for i in range(1, 23)] + ['X', 'Y', 'MT'])]

    # Remove the 'chr' prefix.
    chrom = chrom.str.replace('chr', '')

    # Convert the chromosome names to integers.
    chrom = chrom.replace('X', '23')
    chrom = chrom.replace('Y', '24')
    chrom = chrom.replace('MT', '25')
    chrom = chrom.astype(int)

    # Check if any chromosome names are missing.
    if chrom.isnull().values.any():
        logging.error('Chromosome name is missing.')
        sys.exit(1)
    else:
        # Print space-separated chromosome names.
        logging.info('Chromosomes: ' + ' '.join(chrom.unique().astype(str)))

    # Get the start and end positions.
    start = vcf_df['POS']

    # Check if any start positions are missing.
    if start.isnull().values.any():
        logging.error('Start position is missing.')
        sys.exit(1)

    # Get the SV length from the INFO column.
    sv_length = vcf_df['INFO'].str.extract(r'SVLEN=(-?\d+)', expand=False).astype(np.int32)

    # Check if any SV lengths are missing.
    if sv_length.isnull().values.any():
        logging.error('SV length is missing.')
        sys.exit(1)

    # Get the SV type from the INFO column.
    sv_type = vcf_df['INFO'].str.extract(r'SVTYPE=(\w+)', expand=False)

    # If INFO/REPTYPE=DUP, then the SV type is a duplication.
    sv_type[vcf_df['INFO'].str.contains('REPTYPE=DUP')] = 'DUP'

    # Convert the SV type to integers.
    sv_type = sv_type.replace('DEL', '0')
    sv_type = sv_type.replace('DUP', '1')
    sv_type = sv_type.replace('INV', '2')
    sv_type = sv_type.replace('INS', '3')
    sv_type = sv_type.replace('BND', '4')
    sv_type = sv_type.astype(np.int32)

    # Check if any SV types are missing.
    if sv_type.isnull().values.any():
        logging.error('SV type is missing.')
        sys.exit(1)

    # Loop through the columns and check if any values are missing for all of
    # the feature arrays.
    for col in [chrom, start, sv_length, sv_type, read_support, clipped_bases]:
        if col.isnull().values.all():
            logging.error('All values are missing for a feature.')
            logging.error(col)
            sys.exit(1)

    # Print the first 4 rows of the features.
    logging.info('Features:')
    logging.info(pd.DataFrame({'chrom': chrom.head(4), 'start': start.head(4), 'sv_length': sv_length.head(4), \
                               'sv_type': sv_type.head(4), 'read_support': read_support.head(4), \
                               'clipped_bases': clipped_bases.head(4)}))
    
    # Check that all features have the same length.
    if not all(len(col) == len(chrom) for col in [start, sv_length, sv_type, read_support, clipped_bases]):
        logging.error('Features do not have the same length.')

        # Print the length of each feature.
        logging.error('Chromosomes: ' + str(len(chrom)))
        logging.error('Start positions: ' + str(len(start)))
        logging.error('SV lengths: ' + str(len(sv_length)))
        logging.error('SV types: ' + str(len(sv_type)))
        logging.error('Read support: ' + str(len(read_support)))
        logging.error('Clipped bases: ' + str(len(clipped_bases)))
        
        sys.exit(1)

    # Create a dataframe of the features.
    features = pd.DataFrame({'chrom': chrom, 'start': start, 'sv_length': sv_length, 'sv_type': sv_type, \
                             'read_support': read_support, 'clipped_bases': clipped_bases})

    # Check if any features are missing.
    if features.isnull().values.any():
        logging.error('Features are missing.')

        # Get the rows with missing features.
        missing_features = features[features.isnull().any(axis=1)]

        # Print the rows with missing features.
        logging.error(missing_features)
        sys.exit(1)

    # Return the features.
    return features
