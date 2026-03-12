#This script trains a structures matrix scaling (or other specified type) based calibrator for MINERvA DANN predictions
#Based on https://arxiv.org/pdf/2511.03685

import numpy as np
from tqdm import tqdm
import os, joblib, argparse
import matplotlib.pyplot as plt

import sys


from threading import Thread

import ROOT

#Note, due to the enviroment and package requirements for this script, it works most reliability in a fresh environment 
#(i.e just after login) which has all the necessary packages and little more. I recommend making a venv with he the reqs
# from requirements.txt and just activating that immediately after login and then running


ROOT.gInterpreter.Declare("""
#include <cmath>
#include "ROOT/RVec.hxx"
using namespace ROOT::VecOps;

// --- Helper: Theta3D ---
double theta3D(double thetaX, double thetaY) {
    double sec_thetaX = 1.0 / std::cos(thetaX);
    double sec_thetaY = 1.0 / std::cos(thetaY);
    double denom = (sec_thetaX * sec_thetaX) + (sec_thetaY * sec_thetaY) - 1.0;
    if (denom <= 0) return 0.0;
    double inter = std::sqrt(1.0 / denom);
    double angle = std::acos(inter);
    if (thetaX > 1.57079 || thetaY > 1.57079) return (3.14159 - angle);
    return angle;
}

// --- Optimization: Split Geometry Checks ---
// Only needs scalar inputs (fast to read)
bool passes_geometry(RVec<float> vtx, float thetaX, float thetaY) {
    // Apothem
    double hexSlope = -0.57735; 
    double apothem = 850.0;
    if ( !((std::abs(vtx[1]) < hexSlope * std::abs(vtx[0]) + 2.0 * apothem / 1.73205) && (std::abs(vtx[0]) < apothem)) ) return false;
    
    // Angle
    if ( !(theta3D(thetaX, thetaY) < (17.0 * 3.14159 / 180.0)) ) return false;
    return true;
}

// --- Optimization: Energy Checks ---
// Needs Vector inputs (SLOW to read, do this last)
bool passes_energy(RVec<float> corrected_p, float trk_qp, float trk_eqp_qp, int used_curve, bool anuMode) {
    // Energy Range
    if (corrected_p.size() < 4) return false;
    double emu = corrected_p[3] / 1000.0;
    if (emu < 2.0 || emu > 20.0) return false;

    // Neutrino Mode
    if (anuMode) {
        if (trk_qp <= 0) return false;
    } else {
        if (trk_qp >= 0) return false;
    }
    
    // Significance
    if (used_curve == 1 && trk_eqp_qp != 0) {
        double relativeErr = 1.0 / trk_eqp_qp;
        if (anuMode) {
            if (relativeErr < 5.0) return false;
        } else {
            if (relativeErr > -5.0) return false;
        }
    }
    return true;
}
""")

ROOT.gInterpreter.Declare("""
using namespace ROOT::VecOps;
RVec<float> map_probs(RVec<float> segments, RVec<float> probs, int n_classes) {
    float probSum = probs[0]+probs[1]+probs[2];
    float probRemainder = (probSum<1.0) ? (1.0 - probSum)/(n_classes-3) : 0;
    
    RVec<float> result(n_classes, probRemainder);
    for (size_t i = 0; i < 3; ++i) {
        int idx = (int)segments[i];
        if (idx >= 0 && idx < n_classes) {
            result[idx] = probs[i];
        }
    }
    return result;
}
""")


ROOT.gEnv.SetValue("XNet.RequestTimeout", 30)
ROOT.gEnv.SetValue("XNet.ConnectTimeout", 30)
ROOT.gEnv.SetValue("XNet.MaxRedirectCount", 2)
ROOT.gEnv.SetValue("TFile.AsyncPrefetching", 1)

import math
import random

# Col 5 (index 4) = True Label
# Col 6-179 (indices 5-179) = Probabilities
LABEL_COL_IDX = 4 #Col 5 of the DANN raw prediction csv-like text file (index 4) corresponds to the true plane
PROB_COL_START = 6 #Col 7 (index 6) is the start of the actual plane probabilites, with this being plane 0
PROB_COL_END = 220 #Col 221 (index 220) is the end of the actual plane probabilities corresponding to plane 214

#Note, in some earlier iterations of the DANN there were only 174 planes (No DSECAL), if using these for whatever reason make sure to change PROB_COL_END appropriately




class ThreadWithReturn(Thread):
    def __init__(self, group=None, target=None, name=None, args=(), kwargs={}, verbose=None):
        # Initializing the Thread class
        super().__init__(group, target, name, args, kwargs)
        self._return = None

    # Overriding the Thread.run function
    def run(self):
        if self._target is not None:
            self._return = self._target(*self._args, **self._kwargs)

    def join(self):
        super().join()
        return self._return




def get_truth_segment(mdl, plane, tgt_id):

    if tgt_id == 1: return 9
    if tgt_id == 2: return 18
    if tgt_id == 3: return 27
    if tgt_id == 4: return 45
    if tgt_id == 5: return 50

    # --- Module & Plane Lookup Table ---
    # format: { module_id: { plane_id: return_value, ... }, ... }
    segment_map = {
        -5: {1: 1, 2: 2},
        -4: {1: 3, 2: 4},
        -3: {1: 5, 2: 6},
        -2: {1: 7, 2: 8},
        -1: {1: 9, 2: 9},
        0:  {1: 10, 2: 11},
        1:  {1: 12, 2: 13},
        2:  {1: 14, 2: 15},
        3:  {1: 16, 2: 17},
        4:  {1: 18},
        5:  {1: 19, 2: 20},
        6:  {1: 21, 2: 22},
        7:  {1: 23, 2: 24},
        8:  {1: 25, 2: 26},
        9:  {1: 27, 2: 27},
        # mdl 10 is commented out in source
        11: {1: 28, 2: 29},
        12: {1: 30, 2: 31},
        13: {1: 32, 2: 33},
        14: {1: 34, 2: 35},
        -999: {-999: 36},
        15: {1: 37, 2: 38},
        16: {1: 39, 2: 40},
        17: {1: 41, 2: 42},
        18: {1: 43, 2: 44},
        19: {1: 45, 2: 45},
        20: {1: 46, 2: 47},
        21: {1: 48, 2: 49},
        22: {1: 50, 2: 50},
        23: {1: 51, 2: 52},
        24: {1: 53, 2: 54},
        25: {1: 55, 2: 56},
        26: {1: 57, 2: 58},
        27: {1: 59, 2: 60},
        28: {1: 61, 2: 62},
        29: {1: 63, 2: 64},
        30: {1: 65, 2: 66},
        31: {1: 67, 2: 68},
        32: {1: 69, 2: 70},
        33: {1: 71, 2: 72},
        34: {1: 73, 2: 74},
        35: {1: 75, 2: 76},
        36: {1: 77, 2: 78},
        37: {1: 79, 2: 80},
        38: {1: 81, 2: 82},
        39: {1: 83, 2: 84},
        40: {1: 85, 2: 86},
        41: {1: 87, 2: 88},
        42: {1: 89, 2: 90},
        43: {1: 91, 2: 92},
        44: {1: 93, 2: 94},
        45: {1: 95, 2: 96},
        46: {1: 97, 2: 98},
        47: {1: 99, 2: 100},
        48: {1: 101, 2: 102},
        49: {1: 103, 2: 104},
        50: {1: 105, 2: 106},
        51: {1: 107, 2: 108},
        52: {1: 109, 2: 110},
        53: {1: 111, 2: 112},
        54: {1: 113, 2: 114},
        55: {1: 115, 2: 116},
        56: {1: 117, 2: 118},
        57: {1: 119, 2: 120},
        58: {1: 121, 2: 122},
        59: {1: 123, 2: 124},
        60: {1: 125, 2: 126},
        61: {1: 127, 2: 128},
        62: {1: 129, 2: 130},
        63: {1: 131, 2: 132},
        64: {1: 133, 2: 134},
        65: {1: 135, 2: 136},
        66: {1: 137, 2: 138},
        67: {1: 139, 2: 140},
        68: {1: 141, 2: 142},
        69: {1: 143, 2: 144},
        70: {1: 145, 2: 146},
        71: {1: 147, 2: 148},
        72: {1: 149, 2: 150},
        73: {1: 151, 2: 152},
        74: {1: 153, 2: 154},
        75: {1: 155, 2: 156},
        76: {1: 157, 2: 158},
        77: {1: 159, 2: 160},
        78: {1: 161, 2: 162},
        79: {1: 163, 2: 164},
        80: {1: 165, 2: 166},
        81: {1: 167, 2: 168},
        82: {1: 169, 2: 170},
        83: {1: 171, 2: 172},
        84: {1: 173, 2: 174},
        85: {1: 175, 2: 176},
        86: {1: 177, 2: 178},
        87: {1: 179, 2: 180},
        88: {1: 181, 2: 182},
        89: {1: 183, 2: 184},
        90: {1: 185, 2: 186},
        91: {1: 187, 2: 188},
        92: {1: 189, 2: 190},
        93: {1: 191, 2: 192},
        94: {1: 193, 2: 194},
        95: {1: 195, 2: 195},
        96: {1: 196, 2: 196},
        97: {1: 197, 2: 197},
        98: {1: 198, 2: 198},
        99: {1: 199, 2: 199},
        # For mdl >= 100, the C++ code only checked plane == 2
        100: {2: 200},
        101: {2: 201},
        102: {2: 202},
        103: {2: 203},
        104: {2: 204},
        105: {2: 205},
        106: {2: 206},
        107: {2: 207},
        108: {2: 208},
        109: {2: 209},
        110: {2: 210},
        111: {2: 211},
        112: {2: 212},
        113: {2: 213},
        114: {2: 214},
        115: {2: 2},
        116: {2: 0},
        117: {2: 1},
        118: {2: 0},
        119: {2: 2},
    }

    if mdl in segment_map:
        result = segment_map[mdl].get(plane)
        if result is not None:
            return result

    return -999

def getProbDist(plane0, plane1, plane2, prob0, prob1, prob2, n_classes = 215):
    nclasses = n_classes
    probsum = prob0+prob1+prob2
    probspread = (1 - probsum)/(nclasses-3)
    if probsum < 1:
        outarr = np.full(nclasses, probspread)
    else:
        outarr = np.zeros(nclasses, probspread)
    outarr[plane0] = prob0
    outarr[plane1] = prob1
    outarr[plane2] = prob2
    return outarr

def passesCuts(event, anuMode = False): #Unused function since we're using the cling interpreter to perform our cuts in the interest of reading speed
    #Apothem cut
    hexSlope = -1./math.sqrt(3.)
    apothem = 850
    vtx_x = float(event.vtx[0])
    vtx_y = float(event.vtx[1])
    in_apothem = (abs(vtx_y) < hexSlope*abs(vtx_x) + 2.*apothem/math.sqrt(3.)) and (abs(vtx_x) < apothem)
    if not in_apothem:
        #print("Failed Apothem Cut")
        return False
    #Max Muon Angle
    max_angle = 17 #degrees
    def theta3D(thetaX, thetaY):
        sec_thetaX = 1./math.cos(thetaX)
        sec_thetaY = 1./math.cos(thetaY)
        inter = math.sqrt(1./(sec_thetaX*sec_thetaX + sec_thetaY*sec_thetaY -1.))
        angle = math.acos(inter)

        if(thetaX > 3.14159265/2 or thetaY > 3.14159265/2 ):
            return (3.14159265 - angle)
        else:
            return angle
    muon_thetaX = float(event.muon_thetaX)
    muon_thetaY = float(event.muon_thetaY)
    theta_mu = theta3D(muon_thetaX, muon_thetaY)
    if not (theta_mu < (max_angle*3.14159265/180)):
        #print("Failed Muon Angle Cut")
        return False
    #MINOS Match Cut:
    if int(event.has_interaction_vertex) != 1:
        #print("Failed MINOS Match Cut")
        return False
    #Deadtime cut
    if int(event.phys_n_dead_discr_pair_upstream_prim_track_proj) > 1:
        #print("Failed Deadtime Cut")
        return False
    #Is neutrino cut
    if anuMode:
        if not float(event.MasterAnaDev_minos_trk_qp) > 0:
            #print("Failed is anti-neutrino cut")
            return False
    else:
        if not float(event.MasterAnaDev_minos_trk_qp) < 0:
            #print("Failed is neutrino cut")
            return False
    #Muon Curve Significance cut
    significanceMin = 5
    if(int(event.MasterAnaDev_minos_used_curvature) == 1): #Don't make a significance cut if we reconstructed by range
        relativeErr = 1/float(event.MasterAnaDev_minos_trk_eqp_qp)
        if anuMode:
            if not ( relativeErr >= significanceMin ):
                #print("Failed Muon Curve Significance Cut")
                return False
        else:
            if not ( relativeErr <= -significanceMin ):
                #print("Failed Muon Curve Significance Cut")
                return False
    #ANN Muon Energy Min and Max GeV
    GetANNEmuGev = float(event.muon_corrected_p[3])/1000
    if GetANNEmuGev < 2:
        #print("Failed Min Muon Energy Cut")
        return False
    if GetANNEmuGev > 20:
        #print("Failed Max Muon Energy Cut")
        return False
    #Rock muon cut
    if int(event.rock_muons_removed)==1:
        #print("Failed rock muon cut")
        return False
    #Veto Wall cut
    if int(event.VetoWall_event_IsVeto)==1:
        #print("Failed veto wall cut")
        return False
    #print("Event passed all cuts")
    return True


def load_data(file_list_path, subsample_rate=50000, anuMode=False):
    ROOT.EnableImplicitMT()  # Enable threading - NOTE: If files are over slow network (XRootD) it might cause thread contention waiting for data.
    nfilesread = 0
    chain = ROOT.TChain("MasterAnaDev")
    print(f"Reading root file list: {file_list_path}")
    with open(file_list_path, 'r') as f:
        for line in f:
            root_file = line.strip()
            if root_file: 
                chain.Add(root_file)
                nfilesread+=1
    print(f"Finished reading {nfilesread} input root files")

    df = ROOT.RDataFrame(chain)
    #chain.SetCacheSize(30000000) 
    
    #ROOT.RDF.Experimental.AddProgressBar(df)

    #------------------------------------------------------------------------------------------
    #----------------------Begin-Reimplemetation of cuts from my analysis----------------------
    #------------------------------------------------------------------------------------------
    #Alternatively could do these checks in native python, I believe this is a lot faster though


    # Ive grouped the filters by complexity/time to run. I presume RDataFrame does short circuit execution? and filters
    # in the order the filters were added in? We run the quicker checks first so we dont' have to do the more expensive
    # filter checks for those events


    # Add this filter first as it has a very low I/O Cost
    df = df.Filter("hasMLPrediction == 1", "ML Pred Check") \
           .Filter("has_interaction_vertex == 1", "Vertex Check") \
           .Filter("rock_muons_removed == 0", "Rock Muon Check") \
           .Filter("VetoWall_event_IsVeto == 0", "Veto Check") \
           .Filter("phys_n_dead_discr_pair_upstream_prim_track_proj <= 1", "Deadtime Check")

    # Low I/O Cost
    df = df.Filter("passes_geometry(vtx, muon_thetaX, muon_thetaY)", "Geometry Check")

    # Add this filter last as it has a high I/O Cost so we only read the big arrays for events that survived steps 1-3
    anu_str = "true" if anuMode else "false"
    df = df.Filter(f"""
        passes_energy(muon_corrected_p, MasterAnaDev_minos_trk_qp, 
                      MasterAnaDev_minos_trk_eqp_qp, MasterAnaDev_minos_used_curvature, 
                      {anu_str})
    """, "Energy/Significance Check")
    print("Finished energy checks")

    #------------------------------------------------------------------------------------------
    #-----------------------End-Reimplemetation of cuts from my analysis-----------------------
    #------------------------------------------------------------------------------------------
    n_classes = 215 

    #There probably a betwer way to do this
    df = df.Define("ANN_segments_0", "ANN_segments[0]")
    df = df.Define("ANN_segments_1", "ANN_segments[1]")
    df = df.Define("ANN_segments_2", "ANN_segments[2]")
    df = df.Define("ANN_plane_probs_0", "ANN_plane_probs[0]")
    df = df.Define("ANN_plane_probs_1", "ANN_plane_probs[1]")
    df = df.Define("ANN_plane_probs_2", "ANN_plane_probs[2]")
    df = df.Filter("ANN_plane_probs_0 > 0 && ANN_plane_probs_0 < 1", "Valid probabilites") \

    #df = df.Define("X_probs_vec", f"map_probs(ANN_segments, ANN_plane_probs, {n_classes})")
    cols_to_read = ["ANN_segments_0", "ANN_segments_1", "ANN_segments_2", "ANN_plane_probs_0", "ANN_plane_probs_1", "ANN_plane_probs_2"]
    print("Produced X_probs")

    is_mc = hasattr(chain, "truth_vtx_module")
    if is_mc:
        cols_to_read.extend(["truth_vtx_module", "truth_vtx_plane", "truth_targetID"])

    import time
    start = time.time()
    evCount = df.Count().GetValue()
    if evCount > subsample_rate and subsample_rate>0:
        subsample_rate=subsample_rate/evCount
        df = df.Filter(f"gRandom->Rndm() < {subsample_rate}", f"Subsampling to approx {subsample_rate} events")

    print("Processing events in RDataFrame...")
    data_dict = df.AsNumpy(cols_to_read)
    report = df.Report()
    report.Print() 
    end = time.time()
    print(f"Time elapsed: {(end-start)}")
    print("Generating probability distributions")
    X_probs = []

    for seg0, seg1, seg2, prob0, prob1, prob2 in zip(data_dict["ANN_segments_0"],data_dict["ANN_segments_1"], data_dict["ANN_segments_2"], data_dict["ANN_plane_probs_0"],data_dict["ANN_plane_probs_1"], data_dict["ANN_plane_probs_2"] ):
        X_probs.append(getProbDist(seg0, seg1, seg2, prob0, prob1, prob2, n_classes))
    print("Flattening to np array")
    X_probs = np.array(X_probs)
    print(f"Data Loaded. Shape: {X_probs.shape}")
    out_of_bounds_mask = (X_probs < 0) | (X_probs > 1)


    bad_values = X_probs[out_of_bounds_mask]
    if bad_values.size > 0:
        print(f"Found {bad_values.size} values outside the [0, 1] range!")
        row_indices, col_indices = np.where(out_of_bounds_mask)
        #for r, c, val in zip(row_indices, col_indices, bad_values):
        #    print(f"Row {r}, Col {c}: {val}")
    else:
        print("All values are safely within the 0 to 1 range.")
    
    Y = []
    if is_mc:
        for m, p, t in zip(data_dict["truth_vtx_module"],data_dict["truth_vtx_plane"], data_dict["truth_targetID"]):
            Y.append(get_truth_segment(m, p, t))
        pass
    return X_probs, np.array(Y)

def calc_ECE(probs, y_true, n_ranks=3, n_bins=10):
    print(f"Calculating ECE")
    top_preds_indices = np.argsort(probs, axis=1)[:, -n_ranks:][:, ::-1] # Get the indices of the top N predictions, argsort sorts ascending, so we take the last n_ranks columns and reverse
    ECE_out = []
    for rank in range(n_ranks):
        rank_class_indices = top_preds_indices[:, rank]
        
        confidences = probs[np.arange(len(probs)), rank_class_indices]

        accuracies = (rank_class_indices == y_true).astype(float)
        
        # --- Binning ---
        bin_boundaries = np.linspace(0, 1, n_bins + 1)
        bin_lowers = bin_boundaries[:-1]
        bin_uppers = bin_boundaries[1:]
        
        ece = 0.0
        bin_accs = []
        bin_confs = []
        
        for bin_lower, bin_upper in zip(bin_lowers, bin_uppers):
            in_bin = (confidences > bin_lower) & (confidences <= bin_upper)
            prop_in_bin = in_bin.mean()
            if prop_in_bin > 0:
                accuracy = accuracies[in_bin].mean()
                avg_conf = confidences[in_bin].mean()
                #print(f'accuracy: {accuracy}, avg_conf: {avg_conf}')
                ece += np.abs(avg_conf - accuracy) * prop_in_bin
                
                bin_accs.append(accuracy)
                bin_confs.append(avg_conf)
            else:
                bin_accs.append(np.nan)
                bin_confs.append(np.nan)
        ECE_out.append(ece)

    print("Finished ECE calculation")
    return ECE_out

PLANE_CENTRES = [
    4293.04, 4313.68, 4337.25, 4357.9, 4381.47, 4402.11, 4425.68, 4446.33,
    4514.11, 4534.76, 4558.33, 4578.97, 4602.54, 4623.19, 4646.76, 4667.4,
    4735.19, 4755.83, 4779.4, 4800.05, 4823.62, 4844.26, 4867.83, 4888.48,
    5000.48, 5021.12, 5044.69, 5065.34, 5088.91, 5109.55, 5133.12, 5153.77,
    5456.74, 5477.38, 5500.95, 5521.6, 5545.17, 5565.81, 5589.38, 5610.02,
    5677.81, 5698.45, 5722.03, 5742.67, 5810.45, 5831.1, 5855.68, 5876.33,
    5900.91, 5921.56, 5946.14, 5966.79, 5991.37, 6012.01, 6036.6, 6057.24,
    6081.83, 6102.47, 6127.06, 6147.7, 6172.29, 6192.93, 6217.52, 6238.16,
    6262.74, 6283.39, 6307.97, 6328.62, 6353.2, 6373.85, 6398.43, 6419.08,
    6443.66, 6464.3, 6488.89, 6509.53, 6534.12, 6554.76, 6579.35, 6599.99,
    6624.58, 6645.22, 6669.81, 6690.45, 6715.03, 6735.68, 6760.26, 6780.91,
    6805.49, 6826.14, 6850.72, 6871.37, 6895.95, 6916.59, 6941.18, 6961.82,
    6986.41, 7007.05, 7031.64, 7052.28, 7076.87, 7097.51, 7122.1, 7142.74,
    7167.32, 7187.97, 7212.55, 7233.2, 7257.78, 7278.43, 7303.01, 7323.66,
    7348.24, 7368.88, 7393.47, 7414.11, 7438.7, 7459.34, 7483.93, 7504.57,
    7529.16, 7549.8, 7574.39, 7595.03, 7619.61, 7640.26, 7664.84, 7685.49,
    7710.07, 7730.72, 7755.3, 7775.95, 7800.53, 7821.17, 7845.76, 7866.4,
    7890.99, 7911.63, 7936.22, 7956.86, 7981.45, 8002.09, 8026.68, 8047.32,
    8071.9, 8092.55, 8117.13, 8137.78, 8162.36, 8183.01, 8207.59, 8228.24,
    8252.82, 8273.46, 8298.05, 8318.69, 8343.28, 8363.92, 8388.51, 8409.15,
    8433.74, 8454.38, 8478.97, 8499.61, 8524.19, 8544.84, 8569.42, 8590.07,
    8614.65, 8635.3, 8659.46, 8680.1, 8704.26, 8724.9, 8749.06, 8769.71,
    8793.86, 8814.51, 8838.67, 8859.31, 8883.47, 8904.11, 8928.27, 8948.92,
    8973.08, 8993.72, 9017.88, 9038.52, 9088.08, 9135.41, 9182.75, 9230.08,
    9277.41, 9324.74, 9372.08, 9419.41, 9466.74, 9514.07, 9561.41, 9608.74,
    9656.07, 9703.4, 9750.74, 9798.07, 9845.4, 9892.73, 9940.07, 9987.4
]

def get_nearest_segment_to_zpos(zpos):
    util = PlotUtils.TargetUtils()

    if util.InTarget1ZMC(zpos, 26): return 9
    elif util.InTarget2ZMC(zpos, 26): return 18
    elif util.InTarget3ZMC(zpos, 6): return 27
    elif util.InWaterTargetZMC(zpos): return 36
    elif util.InTarget4ZMC(zpos, 26): return 45
    elif util.InTarget5ZMC(zpos, 26): return 50

    if zpos < 4290.11: return 0
    if zpos > 10034.73: return 0

    # bisect_left returns the insertion point index in the sorted list
    idx = bisect.bisect_left(PLANE_CENTRES, zpos)

    if idx == 0:
        return 0
    if idx == len(PLANE_CENTRES):
        return len(PLANE_CENTRES) - 1


    val_above = PLANE_CENTRES[idx]
    val_below = PLANE_CENTRES[idx - 1]
    if abs(zpos - val_below) < abs(zpos - val_above): # The element BEFORE the iterator is closer
        return idx
    else: # The element AT the iterator is closer
        return idx + 1

get_nearest_segment_to_zpos_vec = np.vectorize(get_nearest_segment_to_zpos)

# Defined globally so it's loaded into memory exactly once
_ZPOS_LOOKUP = (
    0.0,        # Index 0 (Padding for invalid/missing segment)
    4293.04,    # 1
    4313.68,    # 2
    4337.25,    # 3
    4357.9,     # 4
    4381.47,    # 5
    4402.11,    # 6
    4425.68,    # 7
    4446.33,    # 8
    4481.21,    # 9: Target 1
    4514.11,    # 10
    4534.76,    # 11
    4558.33,    # 12
    4578.97,    # 13
    4602.54,    # 14
    4623.19,    # 15
    4646.76,    # 16
    4667.4,     # 17
    4702.29,    # 18: Target 2
    4735.19,    # 19
    4755.83,    # 20
    4779.4,     # 21
    4800.05,    # 22
    4823.62,    # 23
    4844.26,    # 24
    4867.83,    # 25
    4888.48,    # 26
    4923.36,    # 27: Target 3
    5000.48,    # 28
    5021.12,    # 29
    5044.69,    # 30
    5065.34,    # 31
    5088.91,    # 32
    5109.55,    # 33
    5133.12,    # 34
    5153.77,    # 35
    5310.0,     # 36: Water Target
    5456.74,    # 37
    5477.38,    # 38
    5500.95,    # 39
    5521.6,     # 40
    5545.17,    # 41
    5565.81,    # 42
    5589.38,    # 43
    5610.02,    # 44
    5644.91,    # 45: Target 4
    5677.81,    # 46
    5698.45,    # 47
    5722.03,    # 48
    5742.67,    # 49
    5777.55,    # 50: Target 5
    # Standard planes continue
    5810.45, 5831.1, 5855.68, 5876.33, 5900.91, 5921.56, 5946.14, 5966.79,
    5991.37, 6012.01, 6036.6, 6057.24, 6081.83, 6102.47, 6127.06, 6147.7,
    6172.29, 6192.93, 6217.52, 6238.16, 6262.74, 6283.39, 6307.97, 6328.62,
    6353.2, 6373.85, 6398.43, 6419.08, 6443.66, 6464.3, 6488.89, 6509.53,
    6534.12, 6554.76, 6579.35, 6599.99, 6624.58, 6645.22, 6669.81, 6690.45,
    6715.03, 6735.68, 6760.26, 6780.91, 6805.49, 6826.14, 6850.72, 6871.37,
    6895.95, 6916.59, 6941.18, 6961.82, 6986.41, 7007.05, 7031.64, 7052.28,
    7076.87, 7097.51, 7122.1, 7142.74, 7167.32, 7187.97, 7212.55, 7233.2,
    7257.78, 7278.43, 7303.01, 7323.66, 7348.24, 7368.88, 7393.47, 7414.11,
    7438.7, 7459.34, 7483.93, 7504.57, 7529.16, 7549.8, 7574.39, 7595.03,
    7619.61, 7640.26, 7664.84, 7685.49, 7710.07, 7730.72, 7755.3, 7775.95,
    7800.53, 7821.17, 7845.76, 7866.4, 7890.99, 7911.63, 7936.22, 7956.86,
    7981.45, 8002.09, 8026.68, 8047.32, 8071.9, 8092.55, 8117.13, 8137.78,
    8162.36, 8183.01, 8207.59, 8228.24, 8252.82, 8273.46, 8298.05, 8318.69,
    8343.28, 8363.92, 8388.51, 8409.15, 8433.74, 8454.38, 8478.97, 8499.61,
    8524.19, 8544.84, 8569.42, 8590.07, 8614.65, 8635.3, 8659.46, 8680.1,
    8704.26, 8724.9, 8749.06, 8769.71, 8793.86, 8814.51, 8838.67, 8859.31,
    8883.47, 8904.11, 8928.27, 8948.92, 8973.08, 8993.72, 9017.88, 9038.52,
    9088.08, 9135.41, 9182.75, 9230.08, 9277.41, 9324.74, 9372.08, 9419.41,
    9466.74, 9514.07, 9561.41, 9608.74, 9656.07, 9703.4, 9750.74, 9798.07,
    9845.4, 9892.73, 9940.07, 9987.4
)
def get_zpos_from_segment(segment: int) -> float:
    if 1 <= segment <= 214:
        return _ZPOS_LOOKUP[segment]
    return 0.0

get_zpos_from_segment_vec = np.vectorize(get_zpos_from_segment)

# ==========================================
# Begin - Plotting functions
# ==========================================

def plot_sms_calibration_matrix_old(calibrator, outpath, ax=None, cmap='viridis' ): #Using viridis colourmap for color-blind friendliness
    import matplotlib
    matplotlib.use('agg')

    if hasattr(calibrator, 'W'):
        matrix = calibrator.W
    else:
        raise AttributeError(
            "Could not find the matrix attribute on the SMSCalibrator object. "
            "Ensure the calibrator has been fitted via .fit() before plotting."
        )
    
    matrix = np.asarray(matrix)
    if matrix.ndim != 2:
        raise ValueError(f"Expected a 2D matrix, but got shape {matrix.shape}")

    if ax is None:
        fig, ax = matplotlib.pyplot.subplots(figsize=(8, 6))
    else:
        fig = ax.figure
        
    cax = ax.imshow(matrix, cmap=cmap, aspect='auto')
    fig.colorbar(cax, ax=ax, label='Weight / Scaling Value')
    
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

    ax.invert_yaxis()
    ax.figure.savefig(outpath, bbox_inches='tight', dpi=300)
    plt.clf()


def plot_sms_calibration_matrix(matrix, outpath, ax=None, cmap='viridis' ):  #Using viridis colourmap for color-blind friendliness
    """
    Plots the calibration matrix from a fitted probmetrics SMSCalibrator object 
    as a 2D heatmap.
    
    Parameters:
    - calibrator: A fitted SMSCalibrator instance.
    - ax: Optional matplotlib Axes object to plot on.
    - cmap: String representing the matplotlib colormap.
    """

    matrix = np.asarray(matrix)
    if matrix.ndim != 2:
        raise ValueError(f"Expected a 2D matrix, but got shape {matrix.shape}")
    print(f"matrix.shape: {matrix.shape}")
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

    ax.invert_yaxis()
    ax.figure.savefig(outpath, bbox_inches='tight', dpi=300)
    return fig, ax

def plot_multirank_reliability(probs, y_true, title, filename, n_ranks=3, n_bins=10): #Plots N reliability diagrams for the top N predicted ranks.

    print(f"Generating reliability diagram: {title}")
    
    plt.figure(figsize=(6 * n_ranks, 5))
    
    # Get the indices of the top N predictions
    # argsort sorts ascending, so we take the last n_ranks columns and reverse
    top_preds_indices = np.argsort(probs, axis=1)[:, -n_ranks:][:, ::-1]
    for rank in range(n_ranks):
        plt.subplot(1, n_ranks, rank + 1)
        
        rank_class_indices = top_preds_indices[:, rank]
        
        # Get the probability assigned to this rank
        # Use numpy advanced indexing to grab the prob for the specific class at this rank
        confidences = probs[np.arange(len(probs)), rank_class_indices]

        accuracies = (rank_class_indices == y_true).astype(float)
        
        # --- Binning ---
        bin_boundaries = np.linspace(0, 1, n_bins + 1)
        bin_lowers = bin_boundaries[:-1]
        bin_uppers = bin_boundaries[1:]
        
        ece = 0.0
        bin_accs = []
        bin_confs = []
        
        for bin_lower, bin_upper in zip(bin_lowers, bin_uppers):
            in_bin = (confidences > bin_lower) & (confidences <= bin_upper)
            prop_in_bin = in_bin.mean()
            
            if prop_in_bin > 0:
                accuracy = accuracies[in_bin].mean()
                avg_conf = confidences[in_bin].mean()
                
                ece += np.abs(avg_conf - accuracy) * prop_in_bin
                
                bin_accs.append(accuracy)
                bin_confs.append(avg_conf)
            else:
                bin_accs.append(np.nan)
                bin_confs.append(np.nan)

        # Plot - reproducing reliability plot style often seen in ML papers
        plt.plot([0, 1], [0, 1], "k--", label="Perfectly Calibrated")
        plt.scatter(bin_confs, bin_accs, color='blue')
        plt.plot(bin_confs, bin_accs, color='blue', alpha=0.5)
        plt.bar(bin_lowers, bin_accs, width=1.0/n_bins, align='edge', color='blue', alpha=0.2, edgecolor='black')
        
        plt.xlabel("Confidence")
        plt.ylabel("Accuracy")
        nth = rank+1
        suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(nth if nth < 20 else nth % 10, 'th')
        plt.title(f"{nth}{suffix} Most Proable Plane (ECE: {ece:.4f})")
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.ylim(0, 1/(rank+1))
        plt.xlim(0, 1/(rank+1))

    plt.suptitle(title, fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(filename)
    plt.close()

    print("Finished Plot Generation")

def plot_diff_distribution(uncalib_probs, calib_probs, avg = False, n=1, num_classes=None, ax=None, label=None, color=None, save_filename=None, title=None):
    if num_classes is None:
        num_classes = calib_probs.shape[1]


    calib_indices1 = np.argsort(calib_probs, axis=1)[:, -1]
    calib_indices2 = np.argsort(calib_probs, axis=1)[:, -2]
    calib_indices3 = np.argsort(calib_probs, axis=1)[:, -3]

    calib_probs1 = np.sort(calib_probs, axis=1)[:, -1]
    calib_probs2 = np.sort(calib_probs, axis=1)[:, -2]
    calib_probs3 = np.sort(calib_probs, axis=1)[:, -3]

    calib_zpos1 = get_zpos_from_segment_vec(calib_indices1)
    calib_zpos2 = get_zpos_from_segment_vec(calib_indices2)
    calib_zpos3 = get_zpos_from_segment_vec(calib_indices3)

    if avg:
        calib_z = (calib_probs1*calib_zpos1+calib_probs2*calib_zpos2+calib_probs3*calib_zpos3)/(calib_probs1+calib_probs2+calib_probs3)
    else:
        calib_z=calib_zpos1

    uncalib_indices1 = np.argsort(uncalib_probs, axis=1)[:, -1]
    uncalib_indices2 = np.argsort(uncalib_probs, axis=1)[:, -2]
    uncalib_indices3 = np.argsort(uncalib_probs, axis=1)[:, -3]

    uncalib_probs1 = np.sort(uncalib_probs, axis=1)[:, -1]
    uncalib_probs2 = np.sort(uncalib_probs, axis=1)[:, -2]
    uncalib_probs3 = np.sort(uncalib_probs, axis=1)[:, -3]

    uncalib_zpos1 = get_zpos_from_segment_vec(calib_indices1)
    uncalib_zpos2 = get_zpos_from_segment_vec(calib_indices2)
    uncalib_zpos3 = get_zpos_from_segment_vec(calib_indices3)

    diffz = calib_z - uncalib_zpos1



    # Create figure if axis not provided otherwise use existing figure, this allows drawing multiple plots overlaid, if ever desired
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    bins = np.arange(num_classes + 1) - 0.5 # Center bins on integers
    
    ax.hist(diffz, bins=bins, histtype='step', linewidth=2, label=label, color=color)
    suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(n if n < 20 else n % 10, 'th')
    ax.set_xlabel('Distance (mm)')
    ax.set_ylabel(f'Number of Events ({n}{suffix} Choice)')
    if title is None:
        ax.set_title(f'Distribution of the {n}{suffix} Most Probable Plane')
    else:
        ax.set_title(title)
    ax.set_xlim(-30, 30)
    
    if label:
        ax.legend()

    if save_filename:
        print(f"Saving plot to {save_filename}...")
        ax.figure.savefig(save_filename, bbox_inches='tight', dpi=300)
        
    return ax

def plot_nth_plane_distribution(probs, n=1, num_classes=None, ax=None, label=None, color=None, save_filename=None, title=None):
    """
    Plots a histogram of the plane numbers that correspond to the n-th highest 
    probability for each event and saves it to a file.

    Parameters:
    - probs: (numpy.ndarray) Matrix of shape [N_events, N_classes] containing probabilities.
    - n: (int) The rank to plot. 1 = Most Probable, 2 = 2nd Most Probable, etc.
    - num_classes: (int) Total number of planes. If None, inferred from probs.
    - ax: (matplotlib.axes.Axes) Optional axis to plot on. If None, creates a new figure.
    - label: (str) Label for the legend.
    - color: (str) Color of the histogram line.
    - save_filename: (str) Full path/filename to save the plot (e.g., 'rank1_calibrated.png').
                     If None, the plot is not saved automatically.
    """
    
    if num_classes is None:
        num_classes = probs.shape[1]

    # Sort and pick the n-th highest index (-1 is highest, -2 is 2nd highest, etc.)
    nth_indices = np.argsort(probs, axis=1)[:, -n]

    # Create figure if axis not provided otherwise use existing figure, this allows drawing multiple plots overlaid, if ever desired
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    bins = np.arange(num_classes + 1) - 0.5 # Center bins on integers
    
    ax.hist(nth_indices, bins=bins, histtype='step', linewidth=2, label=label, color=color)
    suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(n if n < 20 else n % 10, 'th')
    ax.set_xlabel('Plane Number')
    ax.set_ylabel(f'Number of Events ({n}{suffix} Choice)')
    if title is None:
        ax.set_title(f'Distribution of the {n}{suffix} Most Probable Plane')
    else:
        ax.set_title(title)
    ax.set_xlim(0, num_classes)
    
    if label:
        ax.legend()

    if save_filename:
        print(f"Saving plot to {save_filename}...")
        ax.figure.savefig(save_filename, bbox_inches='tight', dpi=300)
        
    return ax

def plot_weighted_avg_plane_distribution(probs, avg = False, num_classes=None, ax=None, label=None, color=None, save_filename=None, title=None):
    if not avg:
        return plot_nth_plane_distribution(probs, color=color, label=label, title=title)
    if num_classes is None:
        num_classes = probs.shape[1]

    # Sort and pick the n-th highest index (-1 is highest, -2 is 2nd highest, etc.)
    indices1 = np.argsort(probs, axis=1)[:, -1]
    indices2 = np.argsort(probs, axis=1)[:, -2]
    indices3 = np.argsort(probs, axis=1)[:, -3]

    probs1 = np.sort(probs, axis=1)[:, -1]
    probs2 = np.sort(probs, axis=1)[:, -2]
    probs3 = np.sort(probs, axis=1)[:, -3]

    zpos1 = get_zpos_from_segment_vec(indices1)
    zpos2 = get_zpos_from_segment_vec(indices2)
    zpos3 = get_zpos_from_segment_vec(indices3)

    weightedZ = (probs1*zpos1+probs2*zpos2+probs3*zpos3)/(probs1+probs2+probs3)
    weightedSegment = get_nearest_segment_to_zpos_vec(weightedZ)

    # Create figure if axis not provided otherwise use existing figure, this allows drawing multiple plots overlaid, if ever desired
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    bins = np.arange(num_classes + 1) - 0.5 # Center bins on integers
    
    ax.hist(weightedSegment, bins=bins, histtype='step', linewidth=2, label=label, color=color)
    suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(n if n < 20 else n % 10, 'th')
    ax.set_xlabel('Plane Number')
    ax.set_ylabel(f'Number of Events ({n}{suffix} Choice)')
    if title is None:
        ax.set_title(f'Distribution of the {n}{suffix} Most Probable Plane')
    else:
        ax.set_title(title)
    ax.set_xlim(0, num_classes)
    
    if label:
        ax.legend()

    # Save file
    if save_filename:
        print(f"Saving plot to {save_filename}...")
        ax.figure.savefig(save_filename, bbox_inches='tight', dpi=300)
        
    return ax

def plot_prob_distributions(probs, n=1, num_bins=10, ax=None, label=None, color=None, save_filename=None, title=None):
    """
    Plots a histogram of the plane numbers that correspond to the n-th highest 
    probability for each event and optionally saves it to a file.

    Parameters:
    - probs: (numpy.ndarray) Matrix of shape [N_events, N_classes] containing probabilities.
    - n: (int) The probability to plot. 1 = Most Probable, 2 = 2nd Most Probable, etc.
    - num_classes: (int) Total number of planes. If None, inferred from probs.
    - ax: (matplotlib.axes.Axes) Optional axis to plot on. If None, creates a new figure.
    - label: (str) Label for the legend.
    - color: (str) Color of the histogram line.
    - save_filename: (str) Full path/filename to save the plot (e.g., 'rank1_calibrated.png').
                     If None, the plot is not saved automatically.
    """
    
    # Sort and pick the n-th highest probability (-1 is highest, -2 is 2nd highest, etc.)
    nth_prob = np.sort(probs, axis=1)[:, -n]

    # Create figure if axis not provided otherwise use existing figure, this allows drawing multiple plots overlaid, if ever desired
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    bins = np.linspace(0, 1, num_bins + 1)
    num_evts = len(nth_prob)
    weights = np.full(num_evts, 1/num_evts)
    ax.hist(nth_prob, bins=bins, histtype='step', linewidth=2, label=label, color=color, weights=weights)
    ax.set_xlabel('CNN Probability Score')
    ax.set_ylabel(f'Number of Events')
    if title is None:
        ax.set_title('Event Probability Distribtuion')
    else:
        ax.set_title(title)
    
    if label:
        ax.legend()

    if save_filename:
        print(f"Saving plot to {save_filename}...")
        ax.figure.savefig(save_filename, bbox_inches='tight', dpi=300)
        
    return ax

def plot_prob_distributions_data_vs_mc(mc_uncalib, mc_calib, data_uncalib, data_calib, n=1, num_bins=10, ax=None, label=None, save_filename=None, title=None):
    
    fig, (ax_main, ax_ratio, ax_ratio2) = plt.subplots(
        nrows=3, 
        sharex=True, 
        gridspec_kw={'height_ratios': [3, 1, 1]}, # Top plot is 3x taller
        figsize=(8, 8)
    )
    plt.xlim(0, 1/n)
    plt.subplots_adjust(hspace=0)

    #The ECAL is not used in my (Akeem's) analysis and the DANN performs very differently in that region
    #Focus on the targets and tracker (first 174 modules)
    mc_uncalib =mc_uncalib[:, :173]
    mc_calib =mc_calib[:, :173]
    data_uncalib =data_uncalib[:, :173]
    data_calib =data_calib[:, :173]

    ax_main = plot_prob_distributions(mc_uncalib, n=n, num_bins=num_bins, ax=ax_main, color='#d2c271', label="MC Uncalibrated")
    ax_main = plot_prob_distributions(mc_calib, n=n, num_bins=num_bins, ax=ax_main, color='#CA5454', label="MC Calibrated")
    ax_main = plot_prob_distributions(data_uncalib, n=n, num_bins=num_bins, ax=ax_main, color='#40B0A6', label="Data Uncalibrated")
    ax_main = plot_prob_distributions(data_calib, n=n, num_bins=num_bins, ax=ax_main, color='#9285d5', label="Data Calibrated")

    if n==1:
        ax_main.axvline(0.4, color='gray', linestyle='--', linewidth=1)
    
    suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(n if n < 20 else n % 10, 'th')
    if title is None:
        ax_main.set_title(f'Probability Distribution of the {n}{suffix} Most Probable Plane')

    ax_ratio.axhline(1, color='gray', linestyle='--')
    ax_ratio.set_ylim(0, 2)
    bins = np.linspace(0, 1, num_bins + 1)

    nth_prob_data_uncalib = np.sort(data_uncalib, axis=1)[:, -n]
    nth_prob_mc_uncalib = np.sort(mc_uncalib, axis=1)[:, -n]
    nth_prob_data_calib = np.sort(data_calib, axis=1)[:, -n]
    nth_prob_mc_calib = np.sort(mc_calib, axis=1)[:, -n]

    data_counts_uncalib, _ = np.histogram(nth_prob_data_uncalib, bins=bins) 
    mc_counts_uncalib, _   = np.histogram(nth_prob_mc_uncalib, bins=bins)
    data_counts_calib, _ = np.histogram(nth_prob_data_calib, bins=bins) 
    mc_counts_calib, _   = np.histogram(nth_prob_mc_calib, bins=bins)

    ratio_uncalib = np.divide(data_counts_uncalib, mc_counts_uncalib, out=np.zeros_like(mc_counts_uncalib, dtype=float), where=data_counts_uncalib!=0)
    ratio_calib = np.divide(data_counts_calib, mc_counts_calib, out=np.zeros_like(mc_counts_calib, dtype=float), where=data_counts_calib!=0)

    ratio_uncalib=ratio_uncalib*np.sum(mc_counts_uncalib)/np.sum(data_counts_uncalib)
    ratio_calib=ratio_calib*np.sum(mc_counts_calib)/np.sum(data_counts_calib)

    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    ax_ratio.errorbar(bin_centers, ratio_uncalib, yerr=0, fmt='-', drawstyle='steps-mid', markersize=4, color='red', label='Uncalibrated Data/MC')
    ax_ratio.errorbar(bin_centers, ratio_calib, yerr=0, fmt='-', drawstyle='steps-mid', markersize=4, color='blue', label='Calibrated Data/MC')
    ax_ratio.legend()
    ax_ratio.axhline(1, color='gray', linestyle='--', linewidth=1)
    if n==1:
        ax_ratio.axvline(0.4, color='gray', linestyle='--', linewidth=1)

    ratio_data = np.divide(data_counts_calib, data_counts_uncalib, out=np.zeros_like(data_counts_uncalib, dtype=float), where=data_counts_uncalib!=0)
    ratio_mc = np.divide(mc_counts_calib, mc_counts_uncalib, out=np.zeros_like(mc_counts_uncalib, dtype=float), where=data_counts_calib!=0)

    ratio_data=ratio_data*np.sum(data_counts_uncalib)/np.sum(data_counts_calib)
    ratio_mc=ratio_mc*np.sum(mc_counts_uncalib)/np.sum(mc_counts_calib)

    ax_ratio2.errorbar(bin_centers, ratio_data, yerr=0, fmt='-', drawstyle='steps-mid', markersize=4, color='orange', label='Calibrated/Uncalibrated Data')
    ax_ratio2.errorbar(bin_centers, ratio_mc, yerr=0, fmt='-', drawstyle='steps-mid', markersize=4, color='purple', label='Calibrated/Uncalibrated MC')
    ax_ratio2.legend()
    ax_ratio2.axhline(1, color='gray', linestyle='--', linewidth=1)
    if n==1:
        ax_ratio2.axvline(0.4, color='gray', linestyle='--', linewidth=1)
    ax_ratio2.set_ylim(0, 2)

    if save_filename:
        print(f"Saving plot to {save_filename}...")
        plt.savefig(save_filename, bbox_inches='tight', dpi=300)

def plot_cdf_data_vs_mc(mc_uncalib, mc_calib, data_uncalib, data_calib, n=1, num_bins=10, save_filename=None, title=None, cdfDirection = 1):
    
    cdf_ax = plot_prob_cdf(mc_uncalib, n=n, num_bins=num_bins, color='#d2c271', label="MC Uncalibrated", cdfDirection=cdfDirection)
    cdf_ax = plot_prob_cdf(mc_calib, n=n, num_bins=num_bins, color='#CA5454', ax = cdf_ax, save_filename=save_filename, label="MC Calibrated", cdfDirection=cdfDirection)
    cdf_ax = plot_prob_cdf(data_uncalib, n=n, num_bins=num_bins, color='#40B0A6', ax = cdf_ax, label="Data Uncalibrated", cdfDirection=cdfDirection)
    cdf_ax = plot_prob_cdf(data_calib, n=n, num_bins=num_bins, color='#9285d5', ax = cdf_ax, save_filename=save_filename, label="Data Calibrated", cdfDirection=cdfDirection)
    cdf_ax.set_xlim(0, 1/n)

    suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(n if n < 20 else n % 10, 'th')
    if title is None:
        cdf_ax.set_title(f'Probability Distribution of the {n}{suffix} Most Probable Plane')

    if save_filename:
        print(f"Saving plot to {save_filename}...")
        plt.savefig(save_filename, bbox_inches='tight', dpi=300)

def plot_prob_cdf(probs, n=1, num_bins=10, ax=None, label=None, color=None, save_filename=None, title=None, cdfDirection = 1):
    
    # Sort and pick the n-th highest probability (-1 is highest, -2 is 2nd highest, etc.)
    nth_prob = np.sort(probs, axis=1)[:, -n]

    # Create figure if axis not provided otherwise use existing figure, this allows drawing multiple plots overlaid, if ever desired
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    bins = np.linspace(0, 1, num_bins + 1)
    num_evts = len(nth_prob)
    weights = np.full(num_evts, 1/num_evts)
    ax.hist(nth_prob, bins=bins, histtype='step', cumulative=cdfDirection, linewidth=2, label=label, color=color, weights=weights)
    ax.set_ylim(0, 1)
    ax.set_xlabel('CNN Probability Score')
    ax.set_ylabel(f'Number of Events')
    if title is None:
        ax.set_title('Event Probability Distribtuion')
    else:
        ax.set_title(title)
    
    if label:
        ax.legend()

    if save_filename:
        print(f"Saving plot to {save_filename}...")
        ax.figure.savefig(save_filename, bbox_inches='tight', dpi=300)
        
    return ax

# ==========================================
# End - Plotting functions
# ==========================================

import probmetrics.utils
def validate_probabilities_monkey_patch(X: np.ndarray) -> None: #The validate_probabilities method built into the probmetrics package is too strict for us here so we're loosening it with this monkey patched method
    """
    Checks if X is a valid 2D array of probability distributions.
    Raises ValueError if the validation fails.
    Monkey patched to loosen probability sum requirement
    """
    if X.ndim != 2:
        raise ValueError(f"Input array must be 2D, but got {X.ndim} dimensions.")
    if np.any(X < 0) or np.any(X > 1):
        raise ValueError("All probability values must be within the [0, 1] range.")
    if not np.all(np.isclose(X.sum(axis=1), 1.0, atol=0.01)):
        raise ValueError("Each row of the input array must sum to 1 to be a valid probability distribution.")

probmetrics.utils.validate_probabilities = validate_probabilities_monkey_patch


def running_on_grid(): # Check for running on fermigrid
    if os.getenv("JOBSUBJOBID") is not None:
        return True
    else:
        return False



def SaveCalibratorInfo(calibrator, filebase, ece_results = None):
    onGrid = running_on_grid()
    savedir = filebase
    if onGrid:
        print("Detected that we're running on Fermigrid, will save calibrated info locally and ifdh cp them to {filebase}_*...")
        savedir = "./"
    else:
        print(f"Saving calibrated info to {filebase}_*...")
        #savedir = filebase

    pkl_path = savedir+"_calibrator.pkl"
    matrix_png_path = savedir+"_calibratorMatrix.png"
    matrix_pdf_path = savedir+"_calibratorMatrix.pdf"

    matrix_json_path = savedir+"_calibratorMatrix.json"
    matrix_bin_path = savedir+"_calibratorMatrix.bin"

    # Save calibrator object as pkl
    joblib.dump(calibrator, pkl_path)
    print("Model saved successfully.")
    #Save calibrator matrix image if using matrix scaling
    if isinstance(calibrator, probmetrics.calibrators.SMSCalibrator):
        W_full = calibrator.W
        print(f"W_full.shape: {W_full.shape}")
        weights = W_full[:, :-1]
        bias = W_full[:, -1]
        if hasattr(calibrator.ts.cal_, 'invtemp_'):
            print("Collapsing temp scaling into matrix...")
            alpha = calibrator.ts.cal_.invtemp_ 
            weights = weights * alpha
        matrix = np.hstack([weights, bias.reshape(-1, 1)])
        matrix = matrix.astype(np.float64)

        plot_sms_calibration_matrix(matrix, outpath = matrix_png_path)
        plot_sms_calibration_matrix(matrix, outpath = matrix_pdf_path)

        ECE_1 = -1
        ECE_2 = -1
        ECE_3 = -1
        if ece_results is not None:
            ECE_1 = ece_results[1]
            ECE_2 = ece_results[2]
            ECE_3 = ece_results[3]
        matrix.tofile(matrix_bin_path)
        print(f"Saved to binary file {matrix_bin_path}")
        import json
        export_data = {
            "metadata": {
                "n_classes": matrix.shape[0],
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
            "matrix": matrix.tolist() 
        }

        # Save to JSON
        output_path = matrix_json_path
        with open(output_path, 'w') as f:
            json.dump(export_data, f, indent=4) # indent=4 makes it pretty/readable
        print(f"Saved to json file {matrix_json_path}")
        print("Export complete.")




    if onGrid:
        pkl_path_ifdh = filebase+"_calibrator.pkl"
        matrix_png_path_ifdh = filebase+"_calibratorMatrix.png"
        matrix_pdf_path_ifdh = filebase+"_calibratorMatrix.pdf"
        matrix_json_path_ifdh  = filebase + "+calibratorMatrix.json"
        matrix_bin_path_ifdh  = filebase + "+calibratorMatrix.bin"
        os.system(f"ifdh cp {pkl_path} {pkl_path_ifdh}")
        #Maybe add a quick check to make sure these latter files actually exist becuase they are saved within a conditional block
        os.system(f"ifdh cp {matrix_png_path} {matrix_png_path_ifdh}")
        os.system(f"ifdh cp {matrix_pdf_path} {matrix_pdf_path_ifdh}")
        os.system(f"ifdh cp {matrix_json_path} {matrix_json_path_ifdh}")
        os.system(f"ifdh cp {matrix_bin_path} {matrix_bin_path_ifdh}")

def trainCalibrator(calibratorstr, train_X, train_Y, test_X=None, test_Y=None, rho_train=1, tau_train=1, lambda_intercept_train=1, lambda_diagonal_train=1, lambda_off_diagonal_train=1, filebase = None):
    print(f"Training with hyper-parameter values: lambda_intercept = {lambda_intercept_train}, lambda_diagonal = {lambda_diagonal_train}, lambda_off_diagonal={lambda_off_diagonal_train}, rho = {rho_train}, tau={tau_train}")
    if calibratorstr == "sms":
        from probmetrics.calibrators import SMSCalibrator
        calibrator = SMSCalibrator (
        penalty = 'ridge' ,
        rho = rho_train ,
        tau = tau_train ,
        lambda_intercept = lambda_intercept_train ,
        lambda_diagonal = lambda_diagonal_train ,
        lambda_off_diagonal = lambda_off_diagonal_train,
        #opt="saga"
        )
    elif calibratorstr == "svs":
        from probmetrics.calibrators import SVSCalibrator
        calibrator = SVSCalibrator (
        penalty = 'ridge' ,
        rho = rho_train ,
        tau = tau_train ,
        lambda_intercept = lambda_intercept_train ,
        lambda_diagonal = lambda_diagonal_train,
        #opt="saga"
        )
    else:
        print("Invalid calibrator: {calibratorstr}")
        exit(1)


    print("==============================================")
    print(f"Using {calibratorstr} calibrator")
    print("==============================================")
    print(f"Fitting Calibrator on {len(train_Y)} samples (this may take time)...")
    print(f"train_X: {train_X.shape}")
    print(f"train_Y: {train_Y.shape}")
    calibrator.fit(train_X, train_Y)
    print(f"Finished training")
    test_calibrated = None
    ece=None
    if test_X is not None and test_Y is not None:
        print(f"Calibrating test MC sample with shape {test_X.shape}")
        test_calibrated = calibrator.predict_proba(test_X)
        print("Finihsed calibrating test MC sample")
        uncalib_ece = calc_ECE(test_X, test_Y, n_ranks=3, n_bins=10)
        ece = calc_ECE(test_calibrated, test_Y, n_ranks=3, n_bins=10)
        if filebase is not None:
            onGrid = running_on_grid()
            savedir = filebase
            if onGrid:
                print(f"Detected that we're running on Fermigrid, will save calibrated info locally and ifdh cp them to {filebase}_*...")
                savedir = "./"
            else:
                print(f"Saving calibrated info to {filebase}_*...")
                #savedir = filebase
            results_path = savedir+"_results.txt"
            with open(results_path, "w") as myfile:
                myfile.write(f"{ece[0]},{ece[1]},{ece[2]},{uncalib_ece[0]},{uncalib_ece[1]},{uncalib_ece[2]},{lambda_intercept_train},{lambda_diagonal_train},{lambda_off_diagonal_train},{rho_train},{tau_train}\n")
            if onGrid:
                results_path_ifdh = filebase+"_results.txt"
                os.system(f"ifdh cp {results_path} {results_path_ifdh}")
        else:
            print(f"Uncalibrated ECE: {uncalib_ece}")
            print(f"Calibrated ECE: {ece}")
    if filebase is not None:
        SaveCalibratorInfo(calibrator, filebase, ece)

    return calibrator, test_calibrated, filebase



def produceCalibPlots(calibrator, calibratorstr, test_x, test_calibrated_x, test_y, filebase, data = None):
    
    onGrid = running_on_grid()
    savedir = filebase
    if onGrid:
        print("Detected that we're running on Fermigrid, will save plots locally and ifdh cp them to {filebase}_*...")
        savedir = "./"
    else:
        print(f"Saving plots to {filebase}_*...")
        #savedir = filebase

    #Reliability plots
    plot_multirank_reliability( test_x, test_y, title="Uncalibrated DANN Output - MC", filename=f"{savedir}_reliability_uncalibrated.png", n_ranks=3)
    plot_multirank_reliability(test_calibrated_x, test_y, title="SMS Calibrated Output - MC", filename=f"{savedir}_reliability_calibrated.png", n_ranks=3)
    if onGrid:
        os.system(f"ifdh cp {savedir}_reliability_uncalibrated.png {filebase}_reliability_uncalibrated.png")
        os.system(f"ifdh cp {savedir}_reliability_uncalibrated.png {filebase}_reliability_uncalibrated.png")

    #Plot uncalibrated plane distributions
    avg_ax = plot_weighted_avg_plane_distribution(test_x, color="red", label="Uncalibrated", title=' ')
    n1_ax = plot_nth_plane_distribution(test_x, color="red", label="Uncalibrated", title=f'Distribution of the most probable Plane Before Calibration')
    n2_ax = plot_nth_plane_distribution(test_x, color="red", label="Uncalibrated", n=2, title=f'Distribution of the second most probable Plane Before Calibration')
    n3_ax = plot_nth_plane_distribution(test_x, color="red", label="Uncalibrated", n=3, title=f'Distribution of the third most probable Plane Before Calibration')
    #Plot calibrated plane distributions
    plot_weighted_avg_plane_distribution(test_calibrated_x, color="blue", ax=avg_ax, label="Uncalibrated", save_filename=f"{savedir}_avg_plane.png", title=f'Distribution of the weighted average position')
    plot_nth_plane_distribution(test_calibrated_x, color="blue", label="Calibrated", ax=n1_ax, save_filename=f"{savedir}_top_plane.png", title=f'Distribution of the most probable Plane')
    plot_nth_plane_distribution(test_calibrated_x, color="blue", label="Calibrated", ax=n2_ax, n=2, save_filename=f"{savedir}_second_plane.png", title=f'Distribution of the second most probable Plane')
    plot_nth_plane_distribution(test_calibrated_x, color="blue", label="Calibrated", ax=n3_ax, n=3, save_filename=f"{savedir}_third_plane.png", title=f'Distribution of the third most probable Plane')
    if onGrid:
        os.system(f"ifdh cp {savedir}_avg_plane.png {filebase}_avg_plane.png")
        os.system(f"ifdh cp {savedir}_top_plane.png {filebase}_top_plane.png")
        os.system(f"ifdh cp {savedir}_second_plane.png {filebase}_second_plane.png")
        os.system(f"ifdh cp {savedir}_third_plane.png {filebase}_third_plane.png")

    #Plot vertex residuals
    zdiff_ax = plot_diff_distribution(test_x, test_calibrated_x, color="red", label="Unweighted Calibrated Z[0] - Uncalibrated Z[0]", title=f'Difference between the most probable z position calib - uncalib')
    plot_diff_distribution(test_x, test_calibrated_x, avg = True, ax = zdiff_ax, color="blue", label="Weighted Calibrated Z[0] - Uncalibrated Z[0]", save_filename=f"{savedir}_calibVsUncalib_zdiff.png", title=f'Difference between the most probable z position calib - uncalib')
    if onGrid:
        os.system(f"ifdh cp {savedir}_calibVsUncalib_zdiff.png {filebase}_calibVsUncalib_zdiff.png")

    if data is not None:
        print(f"Calibrating data sample with shape {data.shape}")
        data_calibrated = calibrator.predict_proba(data)
        print("Finihsed calibrating data MC sample")
        #Plot data/MC probability distributions
        plot_prob_distributions_data_vs_mc(test_x, test_calibrated_x, data, data_calibrated, n=1, num_bins=100, save_filename=f"{savedir}_prob_distribution_n1.png")
        plot_prob_distributions_data_vs_mc(test_x, test_calibrated_x, data, data_calibrated, n=2, num_bins=100, save_filename=f"{savedir}_prob_distribution_n2.png")
        plot_prob_distributions_data_vs_mc(test_x, test_calibrated_x, data, data_calibrated, n=3, num_bins=100, save_filename=f"{savedir}_prob_distribution_n3.png")
        if onGrid:
            os.system(f"ifdh cp {savedir}_prob_distribution_n1.png {filebase}_prob_distribution_n1.png")
            os.system(f"ifdh cp {savedir}_prob_distribution_n2.png {filebase}_prob_distribution_n2.png")
            os.system(f"ifdh cp {savedir}_prob_distribution_n3.png {filebase}_prob_distribution_n3.png")

        #Plotting CDFs
        plot_cdf_data_vs_mc(test_x, test_calibrated_x, data, data_calibrated, n=1, num_bins=100, save_filename=f"{savedir}_cdf_n1.png")
        plot_cdf_data_vs_mc(test_x, test_calibrated_x, data, data_calibrated, n=2, num_bins=100, save_filename=f"{savedir}_cdf_n2.png")
        plot_cdf_data_vs_mc(test_x, test_calibrated_x, data, data_calibrated, n=3, num_bins=100, save_filename=f"{savedir}_cdf_n3.png")
        if onGrid:
            os.system(f"ifdh cp {savedir}_cdf_n1.png {filebase}_cdf_n1.png")
            os.system(f"ifdh cp {savedir}_cdf_n2.png {filebase}_cdf_n2.png")
            os.system(f"ifdh cp {savedir}_cdf_n3.png {filebase}_cdf_n3.png")

        #CDF reversed, i.e  proportion of events passing prob>x axis value 
        plot_cdf_data_vs_mc(test_x, test_calibrated_x, data, data_calibrated, n=1, num_bins=100, title = "Proportion of events passing probaility cut", save_filename=f"{savedir}_cdf_n1_rev.png", cdfDirection = -1)
        plot_cdf_data_vs_mc(test_x, test_calibrated_x, data, data_calibrated, n=2, num_bins=100, title = "Proportion of events passing probaility cut", save_filename=f"{savedir}_cdf_n2_rev.png", cdfDirection = -1)
        plot_cdf_data_vs_mc(test_x, test_calibrated_x, data, data_calibrated, n=3, num_bins=100, title = "Proportion of events passing probaility cut", save_filename=f"{savedir}_cdf_n3_rev.png", cdfDirection = -1)
        if onGrid:
            os.system(f"ifdh cp {savedir}_cdf_n1_rev.png {filebase}_cdf_n1_rev.png")
            os.system(f"ifdh cp {savedir}_cdf_n2_rev.png {filebase}_cdf_n2_rev.png")
            os.system(f"ifdh cp {savedir}_cdf_n3_rev.png {filebase}_cdf_n3_rev.png")
    else:
        ax = plot_prob_distributions(test_x, n=1, num_bins=100, color='#d2c271', label="MC Uncalibrated")
        ax = plot_prob_distributions(test_calibrated_x, n=1, num_bins=100, ax=ax, color='#CA5454', save_filename=f"{savedir}_prob_distribution.png", label="MC Calibrated")
        if onGrid:
            os.system(f"ifdh cp {savedir}_prob_distribution.png {filebase}_prob_distribution.png")


def main():
    parser = argparse.ArgumentParser(description="""SMS Calibrator. 
                                     This scirpt is multi-purpose. When no model is supplied, it calibrates it's own. And it also produces a number of relevant plots. 
                                     When supplied with multiple hyper-paremeter values lambdas, rhos, taus it perfoms a grid sampling across all of them with  mulithreading to speed o- this is used
                                     capabilities to determine the optimal hyperparemeter set. 
                                     """)
    parser.add_argument("-i", "--input", required=True, type=str, help="Path to text file containing the list of MC playlists over which to train and test the calibrator")
    parser.add_argument("-d", "--datafile", type=str, help="Path to text file containing the list of data playlists over which to evaluate the calibrator")
    parser.add_argument("-o", "--outpath", type=str, help="The path to copy the outputs to. Note, if this is detected to be running on the grid ifdh is used")
    parser.add_argument("-s", "--samplerate", default=50000, type=float, help="The approximate number of events to use to fit the calibrator. Note that 20% of this will be separated into a testing subset")
    parser.add_argument("-c", "--calibrator", type=str, default="sms", help="Type of calibrator to use. Options: sms or svs")
    parser.add_argument("-t", "--threads", default=1, type=int, help="The number of threads over which to split the job")
    parser.add_argument("-m", "--model", type=str, help="Path to an existing, already fitted model in .pkl formal")
    parser.add_argument('--lambda_int', nargs='+', help='The values to sample for the lambda_intercept')
    parser.add_argument('--lambda_diag', nargs='+', help='The values to sample for the lambda_diagonal')
    parser.add_argument('--lambda_offd', nargs='+', help='The values to sample for the lambda_off_diagonal (Only used in SMS, not SVS calibration)')
    parser.add_argument('--rho', nargs='+', help='The values to sample for the rho variables')
    parser.add_argument('--tau', nargs='+', help='The values to sample for the tau variables')
    parser.add_argument("--rhc",action='store_true', help="Flag determining if we're looking at a RHC/anti-nu playlist'")
    parser.add_argument("--justTrain",action='store_true', help="Flag determining if we should only train/fit the calibrator and not draw the extra plots")
    parser.add_argument("--profile",action='store_true', help="Profile memory usage and runtime of the calibrator fitting")

    args = parser.parse_args()

    print("Beginning calibration")

    file_path = args.input
    calibratorstr = args.calibrator
    outpath = args.outpath
    samplerate = args.samplerate
    model = args.model
    n_threads = args.threads
    justTrain = args.justTrain
    profiling = args.profile

    print(f"Output path: {outpath}")
    #Creating output directory
    #if not os.path.exists(outpath):
    #    os.makedirs(outpath)

    lambda_int_samples=[]
    lambda_diag_samples=[]
    lambda_offd_samples=[]
    rho_samples=[]
    tau_samples=[]
    
    if args.lambda_int is not None:
        lambda_int_samples = list(map(float, args.lambda_int))
        print(f"Using lambda_intercept values: {lambda_int_samples}")
    if args.lambda_diag is not None:
        lambda_diag_samples = list(map(float, args.lambda_diag))
        print(f"Using lambda_diagonal values: {lambda_diag_samples}")
    if args.lambda_offd is not None:
        lambda_offd_samples = list(map(float, args.lambda_offd))
        print(f"Using lambda_off_diagonal values: {lambda_offd_samples}")
    if args.rho is not None:
        rho_samples = list(map(float, args.rho))
        print(f"Using lambda_intercept values: {rho_samples}")
    if args.tau is not None:
        tau_samples = list(map(float, args.tau))
        print(f"Using lambda_intercept values: {tau_samples}")


    anyHyperParms = False
    if (len(lambda_int_samples)>0 or len(lambda_diag_samples)>0 or len(lambda_offd_samples)>0 or len(rho_samples)>0 or len(tau_samples)>0):
        anyHyperParms = True
    fullSetOfHyperParms = False
    if (len(lambda_int_samples)>0 and len(lambda_diag_samples)>0 and len(lambda_offd_samples)>0 and len(rho_samples)>0 and len(tau_samples)>0):
        fullSetOfHyperParms = True
    doGridSearch = False
    if fullSetOfHyperParms and (len(lambda_int_samples)>1 or len(lambda_diag_samples)>1 or len(lambda_offd_samples)>1 or len(rho_samples)>1 or len(tau_samples)>1):
        doGridSearch = True

    #Input checks
    if justTrain is not None and args.model is not None:
        print("Cannot have --justTrain and --model/-m simultaneously set. Since just train tells the script to only perform calibrator training while --model supplies an already trained calibrator")
        exit(1)
    if args.model is not None and anyHyperParms:
        print("Cannot have --model/-m set simultaneously with hyper-parameter options. The hyperparemeters are used when fitting a calibrator but --model supplies an already fitted calibrator")
        exit(1)
    if args.model is None and not fullSetOfHyperParms:
        print("If you're not providing a prefitted calibrator with --model/-m you need to provide a set of hyper-paremeters for a new model/grid search of models to be fitted with")
        exit(1)
    if n_threads<1:
        print("Cannot have fewer than 1 thread being used")
        exit(1)
    if args.model is not None and n_threads>1:
        print("Cannot have --model/-m set simultaneously with --threads since the --threads option is to enable multithreading during a grid search to find an optimally fitted calibrator but --model supplies an already fitted calibrator")
        exit(1)
    if args.model is not None and profiling:
        print("Cannot have --model/-m set simultaneously with --profile since the --profile option is to profile calibrator fitting but --model supplies an already fitted calibrator")
        exit(1)

    isRHC = args.rhc
    if isRHC:
        print("Treating as an anti-neutrino sample")
    

    testSizeRatio = 0.2 # Proportion of original dataset to use for test subset
    if model is not None:
        testSizeRatio = 1 # No need to split into testing and training sets if we're not training
    #Try to load data
    try:
        X_probs, y = load_data(file_path, subsample_rate=samplerate, anuMode=isRHC)
    except Exception as e:
        print(f"Error loading data: {e}")
        exit()
    if args.datafile is not None:
        try:
            X_data, _ = load_data(args.datafile, subsample_rate=(samplerate*testSizeRatio), anuMode=isRHC) #Downscaling by testSizeRatio so we have approximately matching number of entries between data and MC test samples
        except Exception as e:
            print(f"Error loading data: {e}")
            exit()

    #print(f"Label range: {y.min()} to {y.max()}")
    print(f"Probability matrix shape: {X_probs.shape}")

    if model is None:
        print("Splitting data into Train and Test sets...")
        from sklearn.model_selection import train_test_split
        train_X, test_X, train_Y, test_Y = train_test_split(X_probs, y, test_size=testSizeRatio, random_state=42)
        del X_probs, y # Free up memory (delete the full original arrays which can be quite big and eat a lot of RAM)
    else:
        test_X = X_probs
        test_Y = y
    if (doGridSearch):
        print("Supplied with multiple hyper-parameter values. Performing grid search.")
        intercept_range = lambda_int_samples
        diagonal_range = lambda_diag_samples
        off_diagonal_range = lambda_offd_samples
        rho_range = rho_samples
        tau_range = tau_samples
        threads = []
        total_vals = len(intercept_range)*len(diagonal_range)*len(off_diagonal_range)*len(rho_range)*len(tau_range)
        print(f"Total variations to try {total_vals}")
        #Create header for merged .txt results files
        onGrid = running_on_grid()
        savedir = outpath
        if onGrid:
            savedir = "./"
        #savedir = filebase
        headerpath = f"{savedir}/_resultsHeader.txt"
        with open(headerpath, "w") as myfile:
            myfile.write(f"Top class ECE,2nd Class ECE,3rd Class ECE,Top class ECE Uncalibrated,2nd Class ECE Uncalibrated,3rd Class ECE Uncalibrated,Λ_intercept,Λ_diagonal,Λ_off_diagonal,ρ,τ\n")
        if running_on_grid():
            os.system(f"ifdh cp {headerpath} {outpath}/_resultsHeader.txt")
        for inter in intercept_range:
            for diag in diagonal_range:
                for off_diag in off_diagonal_range:
                    for rho_test in rho_range:
                        for tau_test in tau_range:
                            outstr = f"{outpath}/{calibratorstr}_{str(inter)}_{str(diag)}_{str(off_diag)}_{str(rho_test)}_{str(tau_test)}"
                            t = ThreadWithReturn(target=trainCalibrator, args=(calibratorstr, train_X, train_Y, test_X, test_Y, rho_test, tau_test, inter, diag, off_diag, outstr))
                            threads.append(t)
        count = 1
        len_threads = len(threads)
        from more_itertools import batched
        for thread_batch in batched(threads, n_threads):
            for t in thread_batch:
                print(f"Working on #{count} of {len_threads} hyper-paremeter variations")
                t.start()
            for t in thread_batch:
                cal, test_X_calibrated, path = t.join()
                print(f"Completed #{count} of {total_vals} hyper-paremeter variations")
                if not justTrain and args.datafile is not None:
                    produceCalibPlots(cal, calibratorstr, test_X, test_X_calibrated, test_Y, path, data=X_data)
                count=count+1
    else:
        path = f"{outpath}/{calibratorstr}"
        if model is None:
            if profiling:
                import time
                from memory_profiler import memory_usage
                start = time.time()
                mem_usage, (cal, test_X_calibrated, _)= memory_usage( proc = (trainCalibrator, (calibratorstr, train_X, train_Y, test_X, test_Y, rho_samples[0], tau_samples[0], lambda_int_samples[0], lambda_diag_samples[0], lambda_offd_samples[0], path)), retval=True)
                end = time.time()
                print("=" * 30)
                print(f"Fitting calibrator over {samplerate} training events took {(start-end)} seconds and used a total of {max(mem_usage)} MiB")
                print("=" * 30)
            else:
                cal, test_X_calibrated, _= trainCalibrator(calibratorstr, train_X, train_Y, test_X, test_Y, rho_samples[0], tau_samples[0], lambda_int_samples[0], lambda_diag_samples[0], lambda_offd_samples[0], path)
        else:
            print(f"Loading prefitted model from {model}...")
            cal = joblib.load(model)
            print(f"Succesfully loaded prefitted model from {model}...")
            test_X_calibrated = cal.predict_proba(test_X)
        if not justTrain and args.datafile is not None:
            print("X_data: ", X_data)
            produceCalibPlots(cal, calibratorstr, test_X, test_X_calibrated, test_Y, path, data=X_data)

if __name__ == "__main__":
    main()