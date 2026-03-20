import pandas as pd
import numpy as np
from tqdm import tqdm
import os, joblib, argparse
import ROOT
from array import array
from probmetrics.calibrators import SMSCalibrator
# ==========================================
# CONFIGURATION
# ==========================================
CHUNK_SIZE = 100000  # Process 100k rows at a time

LABEL_COL_IDX = 4 #Col 5 of the DANN raw prediction csv-like text file (index 4) corresponds to the true plane
PROB_COL_START = 6 #Col 7 (index 6) is the start of the actual plane probabilites, with this being plane 0
PROB_COL_END = 220 #Col 221 (index 220) is the end of the actual plane probabilities corresponding to plane 214

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="This script creates a root file for a given set of input prediction files, calibrates them, then stores a reduced subset of the calibrated output in a root file.")
    parser.add_argument("-i", "--input", nargs='+', required=True, type=str, help="Path to text files containing raw DANN output")
    parser.add_argument("-m", "--model", required=True, type=str, help="Path to pre-existing model")
    parser.add_argument("-o", "--output", required=True, type=str, help="Required output file")

    args = parser.parse_args()
    file_paths = args.input
    model = args.model
    output = args.output

    if not os.path.isfile(model):
        print(f"Model file not found at {model}. Quitting...")
        exit()

    print(f"Loading model from {model}...")
    calibrator = joblib.load(model)
    calibratorstr = model.split("/")[-1]
    calibratorstr = calibratorstr.split("_")[0]
    print("==============================================")
    print(f"Using premade {calibratorstr} calibrator")
    print("==============================================")

    file = ROOT.TFile(output, "NEW")
    tree = ROOT.TTree("calibrated", "calibrated")

    probs_out = array('f', [0.0] * 5)  # 'f' corresponds to float (32-bit)
    planes_out   = array('i', [0] * 5)    # 'i' corresponds to signed int

    run_out   = array('i', [0] )    # 'i' corresponds to signed int
    subrun_out   = array('i', [0] )    # 'i' corresponds to signed int
    gate_out   = array('i', [0] )    # 'i' corresponds to signed int
    physEvt_out   = array('i', [0] )    # 'i' corresponds to signed int

    tree.Branch("probs", probs_out, "probs[5]/F")
    tree.Branch("planes", planes_out, "planes[5]/I")
    tree.Branch("run", run_out, "run/I")
    tree.Branch("subrun", subrun_out, "subrun/I")
    tree.Branch("gate", gate_out, "gate/I")
    tree.Branch("physEvt", physEvt_out, "physEvt/I")

    for path in file_paths:
        total_bytes = os.path.getsize(path)
        #for i, chunk in enumerate(pd.read_csv(path, header=None, chunksize=CHUNK_SIZE)):
        with open(path, 'r') as f:
            with tqdm(total=total_bytes, unit='B', unit_scale=True, desc="Processing") as pbar:
                reader = pd.read_csv(f, header=None, chunksize=CHUNK_SIZE)
                
                for chunk in reader:
                    
                    raw_labels = chunk.iloc[:, LABEL_COL_IDX].values
                    cnn_labels = chunk.iloc[:, LABEL_COL_IDX+1].values
                    run_labels_raw = chunk.iloc[:, 0].values
                    subrun_labels_raw = chunk.iloc[:, 1].values
                    gate_labels_raw = chunk.iloc[:, 2].values
                    physEvt_labels_raw = chunk.iloc[:, 3].values
                    valid_mask = (cnn_labels > 0) # A mask to keep only the entries that have a CNN prediction - the CNN prediction is stored the in row after the truth plane
                    
                    # Apply filter masks
                    labels = raw_labels[valid_mask].astype(np.int64)
                    run_labels = run_labels_raw[valid_mask].astype(np.int64) 
                    subrun_labels = subrun_labels_raw[valid_mask].astype(np.int64)
                    gate_labels = gate_labels_raw[valid_mask].astype(np.int64)
                    physEvt_labels = physEvt_labels_raw[valid_mask].astype(np.int64)

                    probs = chunk.iloc[valid_mask, PROB_COL_START : PROB_COL_END + 1].values.astype(np.float32)




                    labels_list = []
                    probs_list = []
                    if len(labels) > 0:
                        labels_list.append(labels)
                        probs_list.append(probs)

                    y = np.concatenate(labels_list) #We don't need y, we're not training
                    X_probs = np.concatenate(probs_list)
                    
                    #do stuff to X_probs
                    calibrated_data = calibrator.predict_proba(X_probs)
                    for c in range(len(calibrated_data[:,1])):
                        tmp = calibrated_data[c,:]
                        for n in range(5):
                            nth_indices = np.argsort(tmp)[-(n+1)]
                            nth_prob = np.sort(tmp)[-(n+1)]
                            planes_out[n] = nth_indices
                            probs_out[n] = nth_prob
                        run_out[0] = run_labels[c]
                        subrun_out[0] = subrun_labels[c]
                        gate_out[0] = gate_labels[c]
                        physEvt_out[0] = physEvt_labels[c]
                        tree.Fill()
                    pbar.update(f.tell() - pbar.n)

    file.Write()
    file.Close()

    print(f"File {output} created.")