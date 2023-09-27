"""
train_model.py: Train the binary classification model.
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
    tp_filepath = sys.argv[1]
    fp_filepath = sys.argv[2]
    output_dir = sys.argv[3]

    # Run the program.
    run(tp_filepath, fp_filepath, output_dir)
