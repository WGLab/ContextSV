"""
train_model.py: Train the binary classification model.
"""

import os
import sys
import joblib
import sklearn as sk
from sklearn.linear_model import LogisticRegression
import pandas as pd
import os

# Train the model.
def train(tp_filepath, fp_filepath):
    # Read in the true positive and false positive data.
    tp_data = pd.read_csv(tp_filepath, sep="\t")
    fp_data = pd.read_csv(fp_filepath, sep="\t")

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
def run(tp_filepath, fp_filepath, output_dir):
    # Train the model.
    model = train(tp_filepath, fp_filepath)

    # Save the model
    model_path = os.path.join(output_dir, "model.pkl")
    joblib.dump(model, model_path)

    # Print the model.
    print(model)

    # Return the model.
    return model

# Score the structural variants.
def score(model, cnv_data):
    # Get the features.
    features = cnv_data[["lrr", "baf"]]

    # Score the structural variants.
    scores = model.predict_proba(features)

    # Return the scores.
    return scores

# If this is the main program, run it.
if __name__ == '__main__':
    # Get the command line arguments.
    tp_filepath = sys.argv[1]
    fp_filepath = sys.argv[2]
    output_dir = sys.argv[3]

    # Run the program.
    run(tp_filepath, fp_filepath, output_dir)
