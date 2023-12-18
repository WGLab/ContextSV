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
    python train_model.py data/NA12878_SV_calls.vcf data/NA12878_SV_calls.vcf
    output/
"""

import os
import sys
import joblib
from sklearn.linear_model import LogisticRegression
import pandas as pd


def train(true_positives_filepath, false_positives_filepath):
    """Train the binary classification model."""
    # Read in the true positive and false positive data.
    tp_data = pd.read_csv(true_positives_filepath, sep="\t")
    fp_data = pd.read_csv(false_positives_filepath, sep="\t")

    # Get the training data.
    tp_data["label"] = 1
    fp_data["label"] = 0
    training_data = pd.concat([tp_data, fp_data])

    # Get the features and labels.
    features = training_data[["lrr", "baf"]]
    labels = training_data["label"]

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
    run(tp_filepath, fp_filepath, output_dir)
