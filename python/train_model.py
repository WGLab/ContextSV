"""
train_model.py: Train the binary classification model.

This script trains the binary classification model using the true positive and
false positive data. The true positive data is obtained from a benchmarking
dataset. The false positive data is obtained from running the caller on data
that is known to be negative for SVs. This data can be obtained by running the
caller on a normal sample with known SVs accounted for in the reference genome.

For example for HG002, the true positive data is obtained from the Genome in a
Bottle benchmarking dataset, and the false positive data is obtained from
running the caller on the HG002 normal sample and extracting the SV calls that
are not in the benchmarking dataset. This can be repeated for other samples such
as HG001 and HG005 as long as the known SVs are accounted for.

In the HG002 SV v0.6 dataset, there are low-confidence regions which
are excluded from the true positive data. Thus, we must include true SVs from
other publicly available normal samples with information from complex regions,
such as those aligned to CHM13. 

The model is trained using logistic regression. The features are the LRR and
BAF values. The labels are 1 for true positives and 0 for false positives.

The model is saved to the output directory as a pickle file.

Usage:
    python train_model.py <true_positives_filepath> <false_positives_filepath>
    <output_directory>
    
    true_positives_filepath: Path to the VCF of true positive SV calls obtained
        from a benchmarking dataset.
    false_positives_filepath: Path to the VCF of false positive SV calls
        obtained from running the caller on data that is known to be negative
        for SVs. This data can be obtained by running the caller on a normal
        sample with known SVs accounted for in the reference genome.

    output_directory: Path to the output directory.

Output:
    model.pkl: The binary classification model.

Example:
    python train_model.py data/sv_scoring_dataset/true_positives.vcf
    sv_scoring_dataset/false_positives.vcf data/sv_scoring_dataset/model
"""

import os
import sys
import joblib
import numpy as np
import logging
from sklearn.linear_model import LogisticRegression
import pandas as pd

# Set up the logger.
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# # Set up logging to print to standard output.
# console = logging.StreamHandler()
# console.setLevel(logging.INFO)
# formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
# console.setFormatter(formatter)
# logging.getLogger('').addHandler(console)


def read_vcf(filepath):
    """Read in the VCF file."""
    vcf_df = pd.read_csv(filepath, sep='\t', comment='#', header=None, \
                         names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'], \
                            dtype={'CHROM': str, 'POS': np.int64, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, \
                                   'FILTER': str, 'INFO': str, 'FORMAT': str, 'SAMPLE': str})
    return vcf_df

def extract_features(vcf_df):
    """Extract the features from the VCF file's data."""
    # Extract the read depth and clipped bases from the INFO column.
    read_depth = vcf_df['INFO'].str.extract(r'DP=(\d+)', expand=False).astype(np.int32)

    # Check if any read depths are missing.
    if read_depth.isnull().values.any():
        logging.error('Read depth is missing.')
        sys.exit(1)

    clipped_bases = vcf_df['INFO'].str.extract(r'CLIPDP=(\d+)', expand=False).astype(np.int32)

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
    chrom = chrom.astype(np.int32)

    # Check if any chromosome names are missing.
    if chrom.isnull().values.any():
        logging.error('Chromosome name is missing.')
        sys.exit(1)

    # Get the start and end positions.
    start = vcf_df['POS']

    # Check if any start positions are missing.
    if start.isnull().values.any():
        logging.error('Start position is missing.')
        sys.exit(1)

    end = vcf_df['INFO'].str.extract(r'END=(\d+)', expand=False).astype(np.int32)

    # Check if any end positions are missing.
    if end.isnull().values.any():
        logging.error('End position is missing.')
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
    sv_type = sv_type.replace('CNV', '5')
    sv_type = sv_type.astype(np.int32)

    # Check if any SV types are missing.
    if sv_type.isnull().values.any():
        logging.error('SV type is missing.')
        sys.exit(1)

    # Loop through the columns and check if any values are missing for all of
    # the feature arrays.
    for col in [chrom, start, end, sv_length, sv_type, read_depth, clipped_bases]:
        if col.isnull().values.all():
            logging.error('All values are missing for a feature.')
            logging.error(col)
            sys.exit(1)

    # Create a dataframe of the features.
    features = pd.DataFrame({'chrom': chrom, 'start': start, 'end': end, 'sv_length': sv_length, 'sv_type': sv_type, \
                             'read_depth': read_depth, 'clipped_bases': clipped_bases})

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

def train(true_positives_filepath, false_positives_filepath):
    """Train the binary classification model."""
    # Read in the true positive and false positive VCF files.
    tp_vcf = read_vcf(true_positives_filepath)
    fp_vcf = read_vcf(false_positives_filepath)

    # Extract the features from the VCF files.
    logging.info('Extracting features from the true positive VCF file.')
    tp_data = extract_features(tp_vcf)

    # Check if any features are missing.
    if tp_data.isnull().values.any():
        logging.error('Features are missing.')

        # Get the rows with missing features.
        missing_features = tp_data[tp_data.isnull().any(axis=1)]

        # Print the rows with missing features.
        logging.error(missing_features)
        sys.exit(1)

    logging.info('Extracting features from the false positive VCF file.')
    fp_data = extract_features(fp_vcf)

    # Check if any features are missing.
    if fp_data.isnull().values.any():
        logging.error('Features are missing.')

        # Get the rows with missing features.
        missing_features = fp_data[fp_data.isnull().any(axis=1)]

        # Print the rows with missing features.
        logging.error(missing_features)
        sys.exit(1)

    # Add the labels.
    tp_data['label'] = 1
    fp_data['label'] = 0

    # Combine the true positive and false positive data.
    data = pd.concat([tp_data, fp_data])

    # Get the features and labels.
    features = data[["chrom", "start", "end", "sv_length", "sv_type", "read_depth", "clipped_bases"]]
    labels = data["label"]

    # Check if any features are missing.
    if features.isnull().values.any():
        logging.error('Features are missing.')

        # Get the rows with missing features.
        missing_features = features[features.isnull().any(axis=1)]

        # Print the rows with missing features.
        logging.error(missing_features)
        sys.exit(1)

    # Check if any labels are missing.
    if labels.isnull().values.any():
        logging.error('Labels are missing.')
        sys.exit(1)

    # Train the model.
    model = LogisticRegression()
    model.fit(features, labels)

    # Return the model.
    return model

# Run the program.
def run(true_positives_filepath, false_positives_filepath, output_directory):
    """Run the program."""
    # Train the model.
    model = train(true_positives_filepath, false_positives_filepath)

    # Save the model
    model_path = os.path.join(output_directory, "model.pkl")
    joblib.dump(model, model_path)

    # Print the model.
    print(model)

    # Return the model.
    return model

def score(model, cnv_data):
    """Score the structural variants."""
    # Get the features.
    features = cnv_data[["lrr", "baf"]]

    # Score the structural variants.
    scores = model.predict_proba(features)

    # Return the scores.
    return scores


if __name__ == '__main__':
    # Get the command line arguments.

    # Input VCF of true positive SV calls obtained from a benchmarking dataset.
    tp_filepath = sys.argv[1]

    # Input VCF of false positive SV calls obtained from running the caller on
    # data that is known to be negative for SVs. This data can be obtained by
    # running the caller on a normal sample with known SVs accounted for in the
    # reference genome.
    fp_filepath = sys.argv[2]
    output_dir = sys.argv[3]

    # Run the program.
    logging.info('Training the model.')
    run(tp_filepath, fp_filepath, output_dir)
    logging.info('Done.')
