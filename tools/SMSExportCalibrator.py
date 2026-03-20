import joblib
import numpy as np
import argparse
import torch

#Note, due to the enviroment and package requirements for this script, it works most reliability in a fresh environment 
#(i.e just after login) which has all the necessary packages and little more. I recommend making a venv with he the reqs
# from requirements.txt and just activating that immediately after login and then running


def plot_sms_calibration_matrix(matrix, outpath, ax=None, cmap='viridis' ):
    """
    Plots the calibration matrix from a fitted probmetrics SMSCalibrator object 
    as a 2D heatmap.
    
    Parameters:
    - calibrator: A fitted SMSCalibrator instance.
    - ax: Optional matplotlib Axes object to plot on.
    - cmap: String representing the matplotlib colormap.
    """
    import matplotlib.pyplot as plt
    
    # Ensure standard format
    matrix = np.asarray(matrix)
    if matrix.ndim != 2:
        raise ValueError(f"Expected a 2D matrix, but got shape {matrix.shape}")
    print(f"matrix.shape: {matrix.shape}")
    # 2. Setup the plot
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    else:
        fig = ax.figure
        
    #cax = ax.imshow(matrix, cmap=cmap, aspect='auto', vmin=-0.05, vmax=0.15)

    #from matplotlib.colors import LogNorm
    #cax = ax.imshow(matrix, cmap=cmap, aspect='auto', norm=LogNorm())

    from matplotlib.colors import SymLogNorm
    cax = ax.imshow(matrix, cmap=cmap, aspect='auto', norm=SymLogNorm(linthresh=0.01, linscale=1.0, vmin=-0.05, vmax=0.15, base=10))

    fig.colorbar(cax, ax=ax, label='Weight / Scaling Value')
    
    # 4. Formatting
    ax.set_title("SMS Calibration Matrix")
    ax.set_xlabel("Input Class Logits / Probs")
    ax.set_ylabel("Calibrated Output Classes")
    
    num_rows, num_cols = matrix.shape

    xticks = np.arange(0, num_cols, 25)
    yticks = np.arange(0, num_rows, 25)
    
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    
    ax.set_xticklabels(xticks, rotation=45, ha='right')
    ax.set_yticklabels(yticks, rotation=45)

    # Optional: flip the y-axis so class 0 starts at the top
    ax.invert_yaxis()
    ax.figure.savefig(outpath, bbox_inches='tight', dpi=300)
    return fig, ax

def softmax(x):
    """Numerically stable softmax."""
    e_x = np.exp(x - np.max(x))
    return e_x / e_x.sum(axis=0)

parser = argparse.ArgumentParser(description="SMSSaveToBinary - This script takes a SMS .pkl file and saves the model parameters to a binary format, interpretable by C++ analysis code")
parser.add_argument("-i", "--input", required=True, type=str, help="Path to .pkl file")
parser.add_argument("-o", "--output", required=True, type=str, help="Output directory")

args = parser.parse_args()

input = args.input
outbase = args.output

# 1. Load your calibrator
with open(input, 'rb') as f:
    calibrator = joblib.load(f)

#testX = [88, 87, 89]
#testX_Probs= [0.934575, 0.0634049, 0.00182167]

#testProbs = [0.906584, 0.0906537, 0.00170366]
testProbs = [0.934575, 0.0634049, 0.00182167]
testSum = np.sum(testProbs)
testProbToDistribute = (1-testSum)/212
print("testProbToDistribute: ", testProbToDistribute)
testArr = np.full(215, testProbToDistribute)
#testArr[199] = 0.906584
#testArr[200] = 0.0906537
#testArr[198] = 0.00170366
testArr[88] = 0.934575
testArr[87] = 0.0634049
testArr[89] = 0.00182167

testArr_np = []
#testArr_np.append(testArr)
testArr_np.append(testArr)
testArr_np2 = np.array(testArr_np)

print("Logits: ", calibrator._get_logits(testArr_np2))
test_calib = calibrator.predict_proba(testArr_np2)
print("test_calib ", test_calib)
print(np.max(test_calib), np.argmax(test_calib))


# 2. Extract weights and biases
# Note: Inspect your 'calibrator' object with dir() if you aren't sure of the attribute names.
# Commonly they are stored as .weights, .W, .coef_, or similar.
# We assume 'W' is the 215x215 matrix and 'b' is the 215-length vector.

W_full = calibrator.W
print(f"W_full.shape: {W_full.shape}")
weights = W_full[:, :-1]
bias = W_full[:, -1]
##plot_sms_calibration_matrix(W_full, f"{outbase}/sms_matrix_full_without_temp.png")
##plot_sms_calibration_matrix(W_full, f"{outbase}/sms_matrix_full_without_temp.pdf")

classes = calibrator.classes_
if hasattr(calibrator.ts.cal_, 'invtemp_'):
    print("Collapsing temp scaling into matrix...")
    alpha = calibrator.ts.cal_.invtemp_ 
    weights = weights * alpha
print("Here0")
#plot_sms_calibration_matrix(weights, f"{outbase}/sms_matrix_without_bias.png")
#plot_sms_calibration_matrix(weights, f"{outbase}/sms_matrix_without_bias.pdf")
W_fused = np.hstack([weights, bias.reshape(-1, 1)])
W_fused = W_fused.astype(np.float64)
print("Here1")
#plot_sms_calibration_matrix(W_fused, f"{outbase}/sms_matrix_with_bias.png")
#plot_sms_calibration_matrix(W_fused, f"{outbase}/sms_matrix_with_bias.pdf")
print("Here2")


#b = calibrator.b  # Replace with actual attribute name
print("Here3")
#W_fused = W_fused.transpose()
# Ensure they are float32 (standard for C++ performance) or float64
W_fused = W_fused.astype(np.float64)
#b = b.astype(np.float64)
print("Here4")


print("--- Python Check ---")
print("Top-left 3x3 (Weights):")
print(W_fused[:3, :3])

print("\nLast 3 Bias terms (Last Column):")
print(W_fused[:3, -1])

#Get results, if they exist
import os
results_path = input.replace("calibrator.pkl","results.txt")
print("results path:", results_path)
ECE_1 = -1
ECE_2 = -1
ECE_3 = -1
if os.path.isfile(results_path):
    print("Found results file")
    import csv
    with open(results_path, newline='') as f:
        reader = csv.reader(f)
        row = next(reader) # Only need the first line in a results file
        ECE_1 = row[0]
        ECE_2 = row[1]
        ECE_3 = row[2]

# 3. Save to raw binary
#W_fused.tofile(f"{outbase}/sms_calibrator.bin")
print(f"Saved to binary file {outbase}/sms_calibrator.bin")
import json
# 3. Create a dictionary for export
export_data = {
    "metadata": {
        "n_classes": W_fused.shape[1],
        "description": "Fused SMS Calibrator (Weights * Alpha | Bias)",
        "fit_hyperparams":{
            "rho": calibrator.rho,
            "tau": calibrator.tau,
            "lambda_intercept": calibrator.lambda_intercept,
            "lambda_diagonal": calibrator.lambda_diagonal,
            "lambda_off_diagonal": calibrator.lambda_off_diagonal,
        },
        "ECE_scores":{
            "1": ECE_1,
            "2": ECE_2,
            "3": ECE_3,
        }
    },
    # .tolist() converts numpy array to standard Python list of lists
    "matrix": W_fused.tolist() 
}

# 4. Save to JSON
#output_path = f"{outbase}/sms_calibrator.json"
#with open(output_path, 'w') as f:
#    json.dump(export_data, f, indent=4) # indent=4 makes it pretty/readable
print(f"Saved to json file {outbase}/sms_calibrator.bin")
print("Export complete.")

#import shutil
#shutil.copyfile(input, f"{outbase}/sms_calibrator.pkl")





#inarray = [-12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -6.37498, -0.0980716, -2.40071, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121, -12.2121]
#log_probs = np.asarray(inarray)

print(f"Input to log: {testArr_np2}")

log_probs = np.log(testArr_np2)

print(f"Log probs: {log_probs}")



# .shape[-1] safely gets the number of classes (e.g., 215) whether 1D or 2D
n = log_probs.shape[-1]
print(n)

weights = W_fused[:, :n]
bias = W_fused[:, -1]

if log_probs.ndim == 1:
    # Standard W * x + b for a single 1D sample
    calibrated_logits = (weights @ log_probs) + bias
else:
    # Batched x * W.T + b for 2D arrays
    calibrated_logits = (log_probs @ weights.T)
    print(f"calibrated_logits pre-biasing: {calibrated_logits}")
    calibrated_logits = calibrated_logits + bias
    print(f"calibrated_logits post-biasing: {calibrated_logits}")

print(f"calibrated_logits: {calibrated_logits}")

print("softmax: ", softmax(calibrated_logits))