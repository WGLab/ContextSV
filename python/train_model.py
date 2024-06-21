"""
train_model.py - Train the binary classification model.

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
import logging
import numpy as np
import joblib
import pandas as pd
from sklearn.linear_model import LogisticRegression
import matplotlib.pyplot as plt

from extract_features import extract_features

# Set up the logger.
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def train(true_positives_filepath, false_positives_filepath):
    """Train the binary classification model."""

    # Extract the features from the VCF files.
    logging.info('Extracting features from the true positive VCF file.')
    tp_data = extract_features(true_positives_filepath)

    # Check if any features are missing.
    if tp_data.isnull().values.any():
        logging.error('Features are missing.')

        # Get the rows with missing features.
        missing_features = tp_data[tp_data.isnull().any(axis=1)]

        # Print the rows with missing features.
        logging.error(missing_features)
        sys.exit(1)

    logging.info('Extracting features from the false positive VCF file.')
    fp_data = extract_features(false_positives_filepath)

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

    # Print the number of true positives and false positives.
    logging.info('Number of true labels: %d', tp_data.shape[0])
    logging.info('Number of false labels: %d', fp_data.shape[0])

    # Combine the true positive and false positive data.
    data = pd.concat([tp_data, fp_data])

    # Get the features and labels.
    features = data[["chrom", "start", "sv_length", "sv_type", "read_support", "clipped_bases"]]
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

    # Create the output directory if it does not exist.
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Save the model
    model_path = os.path.join(output_directory, "model.pkl")
    joblib.dump(model, model_path)

    # Print the model.
    print(model)

    # Return the model.
    # return model


if __name__ == '__main__':
    # Get the command line arguments.
    if len(sys.argv) != 4:
        logging.error('Usage: python train_model.py <true_positives_filepath> <false_positives_filepath> <output_directory>\n')
        sys.exit(1)

    # Input VCF of true positive SV calls obtained from a benchmarking dataset.
    tp_filepath = sys.argv[1]

    # Input VCF of false positive SV calls obtained from running the caller on
    # data that is known to be negative for SVs. This data can be obtained by
    # running the caller on a normal sample with known SVs accounted for in the
    # reference genome.
    fp_filepath = sys.argv[2]
    output_dir = sys.argv[3]

    # Run the program.
    logging.info('Training the model...')
    run(tp_filepath, fp_filepath, output_dir)
    logging.info('done.')
