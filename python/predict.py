"""
scoring_model.py: Score the structural variants using the binary classification
model.

Usage:
    scoring_model.py <input> <output> <model>

Arguments:
    <input>     Path to the input VCF file.
    <model>     Path to the model file.
"""

import os
import sys
import logging
import numpy as np
import joblib
import pandas as pd

import matplotlib.pyplot as plt

from extract_features import extract_features


def score(model, input_vcf, output_vcf):
    """Score the structural variants using the binary classification model.

    Args:
        model (str): Path to the model file.
        input_vcf (str): Path to the input VCF file.
        output_vcf (str): Path to the output VCF file.
    """
    # Load the model
    clf = joblib.load(model)

    # Extract the features from the VCF file
    X = extract_features(input_vcf)

    # Predict the labels and get the probabilities
    y_pred = clf.predict_proba(X)

    # logging.info('Predicted labels:\n%s', y_pred)

    # Plot a histogram of the probabilities
    plt.hist(y_pred[:, 1], bins=20)
    plt.xlabel('Probability')
    plt.ylabel('Count')

    # # Save the plot to the input VCF file's directory
    # output_dir = os.path.dirname(output_vcf)
    # output_filepath = os.path.join(output_dir, 'probabilities.png')
    # plt.savefig(output_filepath)
    # logging.info('Saved the plot of the probabilities to %s.', output_filepath)

    # Save the plot to the working directory
    plt.savefig('output/probabilities.png')


if __name__ == '__main__':

    # Model file
    model = sys.argv[1]

    # Input VCF file to score
    input_vcf = sys.argv[2]

    # Output VCF file
    output_vcf = sys.argv[3]

    # Set up the logger
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
    
    # Score the structural variants
    score(model, input_vcf, output_vcf)
    