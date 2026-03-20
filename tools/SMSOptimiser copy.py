#This script trains a structures matrix scaling (or other specified type) based calibrator for MINERvA DANN predictions
#Based on https://arxiv.org/pdf/2511.03685

import numpy as np
from tqdm import tqdm
import os, joblib, argparse


import sys


import threading

import ROOT


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

OUTPUT_DIR = "calibration_results"

# Col 5 (index 4) = True Label
# Col 6-179 (indices 5-179) = Probabilities
LABEL_COL_IDX = 4 #Col 5 of the DANN raw prediction csv-like text file (index 4) corresponds to the true plane
PROB_COL_START = 6 #Col 7 (index 6) is the start of the actual plane probabilites, with this being plane 0
PROB_COL_END = 220 #Col 221 (index 220) is the end of the actual plane probabilities corresponding to plane 214

#Note, in some earlier iterations of the DANN there were only 174 planes (No DSECAL), if using these for whatever reason make sure to change PROB_COL_END appropriately

#Creating output directory
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)







def get_truth_segment(mdl, plane, tgt_id):

    # --- Target ID Checks (Priority) ---
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
        -999: {-999: 36},  # Fixed C++ bug: assumed comparison, not assignment
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
        # For mdl >= 100, the C++ only checked plane == 2
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

    # Perform Lookup
    if mdl in segment_map:
        # We use .get() to return None if the specific plane isn't found for that module
        result = segment_map[mdl].get(plane)
        if result is not None:
            return result

    return -999

def getProbDist(plane0, plane1, plane2, prob0, prob1, prob2, n_classes = 215):
    nclasses = n_classes+1
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

def passesCuts(event, anuMode = False):
    #print("TestCuts1")
    #Apothem cut
    hexSlope = -1./math.sqrt(3.)
    #print("TestCuts1a")
    apothem = 850
    #print("TestCuts1b")
    vtx_x = float(event.vtx[0])
    #print("TestCuts1c")
    vtx_y = float(event.vtx[1])
    #print("TestCuts1d")
    in_apothem = (abs(vtx_y) < hexSlope*abs(vtx_x) + 2.*apothem/math.sqrt(3.)) and (abs(vtx_x) < apothem)
    #print("TestCuts1e")
    if not in_apothem:
        #print("TestCuts1f")
        #print("Failed Apothem Cut")
        return False
    #print("TestCuts1g")
    #Max Muon Angle
    max_angle = 17 #degrees
    #print("TestCuts2")
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
    #print("TestCuts3")
    if int(event.has_interaction_vertex) != 1:
        #print("Failed MINOS Match Cut")
        return False
    #Deadtime cut
    if int(event.phys_n_dead_discr_pair_upstream_prim_track_proj) > 1:
        #print("Failed Deadtime Cut")
        return False
    #Is neutrino cut
    #print("TestCuts4")
    if anuMode:
        if not float(event.MasterAnaDev_minos_trk_qp) > 0:
            #print("Failed is anti-neutrino cut")
            return False
    else:
        if not float(event.MasterAnaDev_minos_trk_qp) < 0:
            #print("Failed is neutrino cut")
            return False
    #Muon Curve Significance cut
    #print("TestCuts5")
    significanceMin = 5
    #print("TestCuts5a")
    if(int(event.MasterAnaDev_minos_used_curvature) == 1): #Don't make a significance cut if we reconstructed by range
        #print("TestCuts5b")
        relativeErr = 1/float(event.MasterAnaDev_minos_trk_eqp_qp)
        if anuMode:
            #print("TestCuts5c")
            if not ( relativeErr >= significanceMin ):
                #print("Failed Muon Curve Significance Cut")
                return False
        else:
            #print("TestCuts5d")
            if not ( relativeErr <= -significanceMin ):
                #print("Failed Muon Curve Significance Cut")
                return False
    #ANN Muon Energy Min and Max GeV
    #print("TestCuts6")
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
    #print("TestCuts7")
    if int(event.VetoWall_event_IsVeto)==1:
        #print("Failed veto wall cut")
        return False
    #print("Event passed all cuts")
    return True




def load_root_file_rdf(file_list_path, subsample_rate=1.0, anuMode = False):
    # 1. Enable implicit multi-threading (uses all cores)
    #ROOT.EnableImplicitMT()
    # 2. C++ Helper function to handle the specific array mapping logic
    # We compile this once via the interpreter
    ROOT.gInterpreter.Declare("""
    #include <cmath>
    #include <vector>
    #include <iostream>
    #include "ROOT/RVec.hxx"

    using namespace ROOT::VecOps;

    // Helper function for the angle calculation
    double theta3D(double thetaX, double thetaY) {
        double sec_thetaX = 1.0 / std::cos(thetaX);
        double sec_thetaY = 1.0 / std::cos(thetaY);
        // Avoid division by zero if needed, though exact logic is copied below
        double denom = (sec_thetaX * sec_thetaX) + (sec_thetaY * sec_thetaY) - 1.0;
        if (denom <= 0) return 0.0; // Safety check
        
        double inter = std::sqrt(1.0 / denom);
        double angle = std::acos(inter);

        if (thetaX > 3.14159265/2.0 || thetaY > 3.14159265/2.0) {
            return (3.14159265 - angle);
        } else {
            return angle;
        }
    }

    // The main cut function
    bool passesCuts_CPP(
        RVec<double> vtx, 
        double muon_thetaX, 
        double muon_thetaY,
        int has_interaction_vertex,
        int n_dead_discr, 
        double trk_qp,
        int used_curvature,
        double trk_eqp_qp,
        RVec<double> muon_corrected_p,
        int rock_muons_removed,
        int veto_wall_is_veto,
        bool anuMode
    ) {
        // 1. Apothem Cut
        double hexSlope = -1.0 / std::sqrt(3.0);
        double apothem = 850.0;
        double vtx_x = vtx[0];
        double vtx_y = vtx[1];
        
        bool in_apothem = (std::abs(vtx_y) < hexSlope * std::abs(vtx_x) + 2.0 * apothem / std::sqrt(3.0)) 
                        && (std::abs(vtx_x) < apothem);
        
        if (!in_apothem) return false;

        // 2. Max Muon Angle Cut
        double max_angle = 17.0; // degrees
        double theta_mu = theta3D(muon_thetaX, muon_thetaY);
        if (!(theta_mu < (max_angle * 3.14159265 / 180.0))) return false;

        // 3. MINOS Match Cut
        if (has_interaction_vertex != 1) return false;

        // 4. Deadtime Cut
        if (n_dead_discr > 1) return false;

        // 5. Is Neutrino Cut
        if (anuMode) {
            if (!(trk_qp > 0)) return false;
        } else {
            if (!(trk_qp < 0)) return false;
        }

        // 6. Muon Curve Significance Cut
        double significanceMin = 5.0;
        if (used_curvature == 1) {
            // Protect against divide by zero if trk_eqp_qp is 0
            if (trk_eqp_qp == 0) return false; 
            
            double relativeErr = 1.0 / trk_eqp_qp;
            if (anuMode) {
                if (!(relativeErr >= significanceMin)) return false;
            } else {
                if (!(relativeErr <= -significanceMin)) return false;
            }
        }

        // 7. ANN Muon Energy
        // Assuming muon_corrected_p has at least 4 elements
        if (muon_corrected_p.size() < 4) return false; 
        double GetANNEmuGev = muon_corrected_p[3] / 1000.0;
        
        if (GetANNEmuGev < 2.0) return false;
        if (GetANNEmuGev > 20.0) return false;

        // 8. Rock Muon Cut
        if (rock_muons_removed == 1) return false;

        // 9. Veto Wall Cut
        if (veto_wall_is_veto == 1) return false;

        return true;
    }
    """)

    ROOT.gInterpreter.Declare("""
    using namespace ROOT::VecOps;
    RVec<float> map_probs(RVec<double> segments, RVec<float> probs, int n_classes) {
        RVec<float> result(n_classes, 0.0);
        for (size_t i = 0; i < 3; ++i) {
            int idx = (int)segments[i];
            if (idx >= 0 && idx < n_classes) {
                result[idx] = probs[i];
            }
        }
        return result;
    }
    """)

    # 3. Initialize DataFrame
    # Read file list into a std::vector or simple list for RDataFrame
    print("Reading in files")
    chain = ROOT.TChain("MasterAnaDev")
    with open(file_list_path, 'r') as f:
        for line in f:
            if line.strip(): chain.Add(line.strip())
    print("Loading data into a RDataFrame")
    df = ROOT.RDataFrame(chain)
    #ROOT.RDF.Experimental.AddProgressBar(df)
    # 4. Define 'Cuts' Logic
    # Note: 'passesCuts' must be converted to a string expression like "lep_pt > 20 && abs(eta) < 2.5"
    # If passesCuts is complex Python, this is the hardest part to migrate.
    # Assuming simple cuts for this example:
    df = df.Filter("hasMLPrediction == 1 && Sum(ANN_plane_probs) > 0")

    anu_mode_str = "true" if anuMode else "false"

    # Define the filter string calling the C++ function
    # Note: We list the Column Names in the order the C++ function expects them
    filter_expression = f"""
    passesCuts_CPP(
        vtx, 
        muon_thetaX, 
        muon_thetaY, 
        has_interaction_vertex, 
        phys_n_dead_discr_pair_upstream_prim_track_proj,
        MasterAnaDev_minos_trk_qp,
        MasterAnaDev_minos_used_curvature,
        MasterAnaDev_minos_trk_eqp_qp,
        muon_corrected_p,
        rock_muons_removed,
        VetoWall_event_IsVeto,
        {anu_mode_str}
    )
    """

    # Apply the filter
    df_filtered = df.Filter(filter_expression)

    # 5. Handle Subsampling
    # We use gRandom to filter inside C++, avoiding I/O for skipped events
    if subsample_rate < 1.0:
        df = df.Filter(f"gRandom->Rndm() < {subsample_rate}")

    # 6. Define the Target Column (X_probs) using the C++ helper
    # Assuming n_classes is known (e.g., 4)
    n_classes = 215 
    df = df.Define("X_probs", f"map_probs(ANN_segments, ANN_plane_probs, {n_classes})")

    # 7. Extract Data to Numpy directly
    # This executes the event loop in C++ and moves memory directly to Numpy
    columns_to_read = ["X_probs"]
    
    # Check if MC
    is_mc = hasattr(chain, "truth_vtx_module") # simplified check
    if is_mc:
        # Define Y using a similar C++ helper if get_truth_segment is complex
        # Or read the raw columns and process Y in numpy if it's fast
        columns_to_read.extend(["truth_vtx_module", "truth_vtx_plane", "truth_targetID"])
    
    # Run the graph
    print("Processing events in RDataFrame...")
    data_dict = df.AsNumpy(columns_to_read)
    print("Here1")
    # 8. Final Formatting
    X_probs = np.vstack(data_dict["X_probs"])
    
    Y = []
    if is_mc:
        # If your Y logic is very complex Python, run it here on the reduced numpy arrays
        # This is fast because you are looping over numpy arrays, not ROOT TChain
        for m, p, t in zip(data_dict["truth_vtx_module"],data_dict["truth_vtx_plane"], data_dict["truth_targetID"]):
            Y.append(get_truth_segment(m, p, t))
        pass

    return X_probs, np.array(Y)




def load_root_file_optimized_old(file_list_path, targetNSamples=50000, anuMode=False):
    # Enable threading - NOTE: If files are over slow network (XRootD), 
    # disable this as it might cause thread contention waiting for data.
    #ROOT.EnableImplicitMT() 
    import random
    # We ONLY keep the math-heavy geometric cuts here. 
    # The simple flag checks will be done in Python strings to save I/O.

    prescaledownrate = 1
    chain = ROOT.TChain("MasterAnaDev")
    # Read the .txt file and add files to the chain
    print(f"Reading root file list: {file_list_path}")
    nfilesprocessed = 0
    nfilesread = 0
    with open(file_list_path, 'r') as f:
        for line in f:
            root_file = line.strip()
            if root_file:  # distinct from empty lines
                nfilesprocessed+=1
                #if random.randint(0,prescaledownrate)==1:
                chain.Add(root_file)
                nfilesread+=1
    print(f"Finished reading root file list, read {nfilesread} of {nfilesprocessed} files")

    df = ROOT.RDataFrame(chain)
    #chain.SetCacheSize(30000000) 
    
    #ROOT.RDF.Experimental.AddProgressBar(df)

    print("Doing low latency checks")
    # 2. INT/BOOL FLAGS (Very Low I/O Cost)
    # Check these first. If an event fails here, we never read the heavy vectors.
    df = df.Filter("hasMLPrediction == 1", "ML Pred Check") \
           .Filter("has_interaction_vertex == 1", "Vertex Check") \
           .Filter("rock_muons_removed == 0", "Rock Muon Check") \
           .Filter("VetoWall_event_IsVeto == 0", "Veto Check") \
           .Filter("phys_n_dead_discr_pair_upstream_prim_track_proj <= 1", "Deadtime Check")
    print("Finished low latency checks")

    # 3. SCALAR FLOATS (Low I/O Cost)
    # Only now do we read vtx and theta
    print("Doing geometry checks")
    df = df.Filter("passes_geometry(vtx, muon_thetaX, muon_thetaY)", "Geometry Check")
    print("Finished geometry checks")

    # 4. HEAVY VECTORS (High I/O Cost)
    # Only read the big arrays for events that survived steps 1-3
    print("Doing energy checks")
    anu_str = "true" if anuMode else "false"
    df = df.Filter(f"""
        passes_energy(muon_corrected_p, MasterAnaDev_minos_trk_qp, 
                      MasterAnaDev_minos_trk_eqp_qp, MasterAnaDev_minos_used_curvature, 
                      {anu_str})
    """, "Energy/Significance Check")
    print("Finished energy checks")

    n_classes = 215 
    #print("Producing X_probs")
    #df = df.Define("X_probs", f"map_probs(ANN_segments, ANN_plane_probs, {n_classes})")
    #print("Produced X_probs")


    # Explode the vector into scalar columns
    df = df.Define("X_probs_vec", f"map_probs(ANN_segments, ANN_plane_probs, {n_classes})")
    cols_to_read = ["X_probs_vec"]
    #for i in range(n_classes):
    #    col_name = f"prob_{i}"
    #    # Extract the i-th element. The casting (float) ensures it's a simple type.
    #    df = df.Define(col_name, f"(float)X_probs_vec[{i}]")
    #    cols_to_read.append(col_name)
    #print(f"cols_to_read {cols_to_read}")
    print("Produced X_probs")


    # This report will tell you exactly how many events passed each stage!
    import time
    start = time.time()
    report = df.Report()
    report.Print() 
    end = time.time()
    print(f"Time elapsed: {(end-start)}")
    print("Here1")
    # 7. Extract Data to Numpy directly
    # This executes the event loop in C++ and moves memory directly to Numpy
    #columns_to_read = ["X_probs"]
    print("Here2")
    # Check if MCsubsample_rate
    is_mc = hasattr(chain, "truth_vtx_module") # simplified check
    if is_mc:
        # Define Y using a similar C++ helper if get_truth_segment is complex
        # Or read the raw columns and process Y in numpy if it's fast
        cols_to_read.extend(["truth_vtx_module", "truth_vtx_plane", "truth_targetID"])

    evCount = df.Count().GetValue()
    print("Surviving events: ", evCount)
    if evCount > targetNSamples:
        subsample_rate=targetNSamples/evCount
        df = df.Filter(f"gRandom->Rndm() < {subsample_rate}", f"Subsampling to approx {targetNSamples} events")

    print("Processing events in RDataFrame...")
    data_dict = df.AsNumpy(cols_to_read)
    report = df.Report()
    print("Test1")
    report.Print() 
    print("Test2")

    #data_dict.sample(targetNSamples)
    print("data_dict: ", data_dict)
    print("data_dict[\"X_probs_vec\"]: ", data_dict["X_probs_vec"])
    print("Test1.1")
    #X_probs = np.vstack(list(data_dict["X_probs_vec"])).astype(np.float32)
    X_probs = np.array([np.asarray(vec) for vec in data_dict["X_probs_vec"]], dtype=np.float32)
    #X_probs = [np.array(v).tolist() for v in data_dict["X_probs_vec"]]
    print("X_probs: ", X_probs)
    #print(prob_cols)
    #X_probs = np.stack(prob_cols, axis=1)

    X_probs = np.vstack(list(data_dict["X_probs_vec"])).astype(np.float64)
    #print(prob_cols)
    #X_probs = np.stack(prob_cols, axis=1)
    print(f"Data Loaded. Shape: {X_probs.shape}")

    print(f"Data Loaded. Shape: {X_probs.shape}")


    out_of_bounds_mask = (X_probs < 0) | (X_probs > 1)

    # 2. Extract the actual values
    bad_values = X_probs[out_of_bounds_mask]

    # 3. Check if there are any bad values and print them
    if bad_values.size > 0:
        print(f"Found {bad_values.size} values outside the [0, 1] range!")
        
        # Optional but highly recommended: get the exact row and column indices
        row_indices, col_indices = np.where(out_of_bounds_mask)
        
        for r, c, val in zip(row_indices, col_indices, bad_values):
            print(f"Row {r}, Col {c}: {val}")
    else:
        print("All values are safely within the 0 to 1 range.")
    
    Y = []
    if is_mc:
        # If your Y logic is very complex Python, run it here on the reduced numpy arrays
        # This is fast because you are looping over numpy arrays, not ROOT TChain
        for m, p, t in zip(data_dict["truth_vtx_module"],data_dict["truth_vtx_plane"], data_dict["truth_targetID"]):
            Y.append(get_truth_segment(m, p, t))
        pass
    print("type(X_probs): ", type(X_probs))
    return X_probs, np.array(Y)


def load_root_file_optimized(file_list_path, targetNSamples=50000, anuMode=False):
    # Enable threading - NOTE: If files are over slow network (XRootD), 
    # disable this as it might cause thread contention waiting for data.
    ROOT.EnableImplicitMT() 
    import random
    # We ONLY keep the math-heavy geometric cuts here. 
    # The simple flag checks will be done in Python strings to save I/O.

    nfilesread = 0
    chain = ROOT.TChain("MasterAnaDev")
    # Read the .txt file and add files to the chain
    print(f"Reading root file list: {file_list_path}")
    with open(file_list_path, 'r') as f:
        for line in f:
            root_file = line.strip()
            if root_file:  # distinct from empty lines
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

    print("Doing low latency checks")
    # Add this filter first as it has a very low I/O Cost - I presume RDataFrame does short circuit execution?
    # Check these first. If an event fails here, we never read the heavy vectors.
    df = df.Filter("hasMLPrediction == 1", "ML Pred Check") \
           .Filter("has_interaction_vertex == 1", "Vertex Check") \
           .Filter("rock_muons_removed == 0", "Rock Muon Check") \
           .Filter("VetoWall_event_IsVeto == 0", "Veto Check") \
           .Filter("phys_n_dead_discr_pair_upstream_prim_track_proj <= 1", "Deadtime Check")
    print("Finished low latency checks")

    # Low I/O Cost
    print("Doing geometry checks")
    df = df.Filter("passes_geometry(vtx, muon_thetaX, muon_thetaY)", "Geometry Check")
    print("Finished geometry checks")

    # Add this filter last as it has a high I/O Cost - I presume RDataFrame does short circuit execution?
    # Only read the big arrays for events that survived steps 1-3
    print("Doing energy checks")
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

    df = df.Define("ANN_segments_0", "ANN_segments[0]")
    df = df.Define("ANN_segments_1", "ANN_segments[1]")
    df = df.Define("ANN_segments_2", "ANN_segments[2]")
    df = df.Define("ANN_plane_probs_0", "ANN_plane_probs[0]")
    df = df.Define("ANN_plane_probs_1", "ANN_plane_probs[1]")
    df = df.Define("ANN_plane_probs_2", "ANN_plane_probs[2]")
    df = df.Filter("ANN_plane_probs_0 > 0 && ANN_plane_probs_0 < 1", "Valid probabilites") \

    # Explode the vector into scalar columns
    #df = df.Define("X_probs_vec", f"map_probs(ANN_segments, ANN_plane_probs, {n_classes})")
    cols_to_read = ["ANN_segments_0", "ANN_segments_1", "ANN_segments_2", "ANN_plane_probs_0", "ANN_plane_probs_1", "ANN_plane_probs_2"]
    print("Produced X_probs")

    is_mc = hasattr(chain, "truth_vtx_module")
    if is_mc:
        cols_to_read.extend(["truth_vtx_module", "truth_vtx_plane", "truth_targetID"])

    import time
    start = time.time()
    evCount = df.Count().GetValue()
    if evCount > targetNSamples:
        subsample_rate=targetNSamples/evCount
        df = df.Filter(f"gRandom->Rndm() < {subsample_rate}", f"Subsampling to approx {targetNSamples} events")

    print("Processing events in RDataFrame...")
    data_dict = df.AsNumpy(cols_to_read)
    report = df.Report()
    report.Print() 
    end = time.time()
    print(f"Time elapsed: {(end-start)}")
    print("Generating probability distributions")
    X_probs = []

    for seg0, seg1, seg2, prob0, prob1, prob2 in zip(data_dict["ANN_segments_0"],data_dict["ANN_segments_1"], data_dict["ANN_segments_2"], data_dict["ANN_plane_probs_0"],data_dict["ANN_plane_probs_1"], data_dict["ANN_plane_probs_2"] ):
        #print(seg0, seg1, seg2, prob0, prob1, prob2)
        #if seg0 > 0 and seg0 < n_classes+1 and seg1 > 0 and seg1 < n_classes+1 and seg2 > 0 and seg2 < n_classes+1 and prob0 > 0 and prob1 > 0 and prob2 > 0:
        X_probs.append(getProbDist(seg0, seg1, seg2, prob0, prob1, prob2, n_classes))
    print("Flattening to np array")
    X_probs = np.array(X_probs)
    print(f"Data Loaded. Shape: {X_probs.shape}")
    out_of_bounds_mask = (X_probs < 0) | (X_probs > 1)

    # 2. Extract the actual values
    bad_values = X_probs[out_of_bounds_mask]

    #Check if there are any bad values and print them - note no longer bother printing them
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
    print("type(X_probs): ", type(X_probs))
    return X_probs, np.array(Y)



def load_root_file(file_list_path, targetNSamples=50000, subsample_rate=1.0):
    """
    Reads a list of ROOT files and returns processed arrays X_probs and Y.

    Args:
        file_list_path (str): Path to .txt file containing list of .root files.
        fill_mode (int): Determines how to fill the values that don't correspond to an ANN_segment number 0 = set all to 0, 1 = distribute probability evenly, 3 = scale up 3 probabilities we're given to fill probability space and put 0 elsewhere

    Returns:
        tuple: (X_probs, Y)
            X_probs: Numpy array of shape (N_events, N_Classes).
            Y: List of truth segments (empty if is_mc=False).
    """
    ROOT.EnableImplicitMT() 
    # 1. Initialize TChain
    chain = ROOT.TChain("MasterAnaDev")
    
    # Read the .txt file and add files to the chain
    print(f"Reading root file list: {file_list_path}")
    with open(file_list_path, 'r') as f:
        for line in f:
            root_file = line.strip()
            if root_file:  # distinct from empty lines
                chain.Add(root_file)
    #print(f"Finished reading root file list, read {nfilesread} of {nfilesprocessed} files")
    #Speed things up by turning off unused branches
    chain.SetBranchStatus("*", 0)
    required_branches = ["ANN_plane_probs", "ANN_segments", "hasMLPrediction", "truth_vtx_module", "truth_vtx_plane", "truth_targetID", "vtx", "muon_thetaX", "muon_thetaY", "has_interaction_vertex", "phys_n_dead_discr_pair_upstream_prim_track_proj", "MasterAnaDev_minos_trk_qp", "MasterAnaDev_minos_trk_eqp_qp", "MasterAnaDev_minos_used_curvature", "muon_corrected_p", "rock_muons_removed", "VetoWall_event_IsVeto"]
    for b in required_branches:
        chain.SetBranchStatus(b, 1)

    n_entries = chain.GetEntries()
    print(f"Total events found: {n_entries}")
    if subsample_rate != 1.0:
        n_entries_interest = int(n_entries*subsample_rate)
        print(f"Subsampling rate: {subsample_rate * 100}%")
    else:
        n_entries_interest = n_entries
    is_mc = False
    if chain.GetBranch("truth_vtx_module"):
        is_mc = True
    print(f"Processing type: {'MC' if is_mc else 'Data'}")
    # 2. Initialize Output Containers
    # Initialize X_probs with zeros. Shape: (Total Events, N_classes)
    n_classes = PROB_COL_END-PROB_COL_START+1
    X_probs_arr = []
    Y_arr = []
    #print("Test1")
    # 3. Loop over events
    # We enumerate to keep track of the event index 'i' for filling X_probs
    subsample_rate = targetNSamples/n_entries
    for i, event in enumerate(tqdm(chain, total=n_entries, unit="event", desc="Loading Data")):
        if subsample_rate<1.0 and random.random()>(subsample_rate*10):
            continue
        #print("Test1.a")
        prob = [float(event.ANN_plane_probs[0]), float(event.ANN_plane_probs[1]), float(event.ANN_plane_probs[2])]
        probsum = prob[0]+prob[1]+prob[2]

        if int(event.hasMLPrediction) == 1 and probsum > 0 and passesCuts(event):
            #print("Test1.b")
            # --- Fill X_probs ---
            # We assume ANN_segments and ANN_plane_probs are arrays of length 3
            # Note: ANN_segments are doubles, so we cast to int for array indexing.
            tmp_X = [(1-probsum)/(n_classes-3)]*n_classes
            #print("Test1.ba")
            for j in range(3):
                #try:
                # Get the segment index (column) and the probability (value)
                seg_idx = int(event.ANN_segments[j])
                seg_prob = prob[j]
                # Safety check: ensure index is within the range
                if 0 <= seg_idx < n_classes:
                    tmp_X[seg_idx] = seg_prob
            #print("Test1.bb")
            #if any(x < 0 for x in tmp_X) or any(x > 1 for x in tmp_X) or not np.isclose(np.sum(tmp_X), 1.0):
            #    print(np.sum(tmp_X), tmp_X, chain.GetReadEntry(), chain.GetFile().GetName(), i, j, seg_idx, seg_prob, event.mc_run, event.mc_subrun, event.mc_nthEvtInFile)
            #    if(any(x < 0 for x in tmp_X)):
            #        print("Case1")
            #    if(any(x > 1 for x in tmp_X)):
            #        print("Case2")
            #    if(not np.isclose(np.sum(tmp_X), 1.0)):
            #        print("Case3: ", np.sum(tmp_X))
            #    continue #This shouldn't really happen
            #print("Test1.c")
            #if i >= n_entries_interest:
            #    break
            X_probs_arr.append(tmp_X)
            #if not np.isclose(np.sum(tmp_X), 1.0): #Rarely we have cases where due to rounding the probabilites sum to more than 1. This is typically when 1 plane is assigned essentiall 1 as probability
                #print(np.sum(tmp_X), chain.GetReadEntry(), chain.GetFile().GetName(), i, j, seg_idx, seg_prob, event.mc_run, event.mc_subrun, event.mc_nthEvtInFile)
            #print("Test1.d")
            # --- Fill Y (MC Only) ---
            if is_mc:
                # Get the truth label
                mdl = event.truth_vtx_module
                plane = event.truth_vtx_plane
                tgt_id = event.truth_targetID
                truth_val = get_truth_segment(mdl, plane, tgt_id)
                Y_arr.append(truth_val)
                #print("truth: ", truth_val, "\treco: ", np.argmax(tmp_X))
    #print("Test2")
    X_probs = np.array(X_probs_arr)
    Y = np.array(Y_arr)
    print(f"Data Loaded. Shape: {X_probs.shape}")
    print("Processing complete.")
    return X_probs, Y


def load_root_file_awk(file_list_path, targetNSamples=50000):
    """
    Reads a list of ROOT files and returns processed arrays X_probs and Y.

    Args:
        file_list_path (str): Path to .txt file containing list of .root files.
        fill_mode (int): Determines how to fill the values that don't correspond to an ANN_segment number 0 = set all to 0, 1 = distribute probability evenly, 3 = scale up 3 probabilities we're given to fill probability space and put 0 elsewhere

    Returns:
        tuple: (X_probs, Y)
            X_probs: Numpy array of shape (N_events, N_Classes).
            Y: List of truth segments (empty if is_mc=False).
    """
    ROOT.EnableImplicitMT() 
    # 1. Initialize TChain
    chain = ROOT.TChain("MasterAnaDev")
    required_branches = ["ANN_plane_probs", "ANN_segments", "hasMLPrediction", "truth_vtx_module", "truth_vtx_plane", "truth_targetID", "vtx", "muon_thetaX", "muon_thetaY", "has_interaction_vertex", "phys_n_dead_discr_pair_upstream_prim_track_proj", "MasterAnaDev_minos_trk_qp", "MasterAnaDev_minos_trk_eqp_qp", "MasterAnaDev_minos_used_curvature", "muon_corrected_p", "rock_muons_removed", "VetoWall_event_IsVeto"]

    treedict = {}
    with open(file_list_path, 'r') as f:
        for line in f:
            root_file = line.strip()
            if root_file:  # distinct from empty lines
                chain.Add(root_file)
                treedict[root_file] = "MasterAnaDev"
    import uproot
    # Read the .txt file and add files to the chain
    print(f"Reading root file list: {file_list_path}")
    tree = uproot.concatenate(treedict, filter_name=required_branches)
    return
    with open(file_list_path, 'r') as f:
        for line in f:
            root_file = line.strip()
            if root_file:  # distinct from empty lines
                chain.Add(root_file)
    #print(f"Finished reading root file list, read {nfilesread} of {nfilesprocessed} files")
    #Speed things up by turning off unused branches
    chain.SetBranchStatus("*", 0)
    required_branches = ["ANN_plane_probs", "ANN_segments", "hasMLPrediction", "truth_vtx_module", "truth_vtx_plane", "truth_targetID", "vtx", "muon_thetaX", "muon_thetaY", "has_interaction_vertex", "phys_n_dead_discr_pair_upstream_prim_track_proj", "MasterAnaDev_minos_trk_qp", "MasterAnaDev_minos_trk_eqp_qp", "MasterAnaDev_minos_used_curvature", "muon_corrected_p", "rock_muons_removed", "VetoWall_event_IsVeto"]
    for b in required_branches:
        chain.SetBranchStatus(b, 1)

    n_entries = chain.GetEntries()
    print(f"Total events found: {n_entries}")
    if subsample_rate != 1.0:
        n_entries_interest = int(n_entries*subsample_rate)
        print(f"Subsampling rate: {subsample_rate * 100}%")
    else:
        n_entries_interest = n_entries
    is_mc = False
    if chain.GetBranch("truth_vtx_module"):
        is_mc = True
    print(f"Processing type: {'MC' if is_mc else 'Data'}")
    # 2. Initialize Output Containers
    # Initialize X_probs with zeros. Shape: (Total Events, N_classes)
    n_classes = PROB_COL_END-PROB_COL_START+1
    X_probs_arr = []
    Y_arr = []
    #print("Test1")
    # 3. Loop over events
    # We enumerate to keep track of the event index 'i' for filling X_probs
    subsample_rate = targetNSamples/n_entries
    for i, event in enumerate(tqdm(chain, total=n_entries, unit="event", desc="Loading Data")):
        if subsample_rate<1.0 and random.random()>(subsample_rate*10):
            continue
        #print("Test1.a")
        prob = [float(event.ANN_plane_probs[0]), float(event.ANN_plane_probs[1]), float(event.ANN_plane_probs[2])]
        probsum = prob[0]+prob[1]+prob[2]

        if int(event.hasMLPrediction) == 1 and probsum > 0 and passesCuts(event):
            #print("Test1.b")
            # --- Fill X_probs ---
            # We assume ANN_segments and ANN_plane_probs are arrays of length 3
            # Note: ANN_segments are doubles, so we cast to int for array indexing.
            tmp_X = [(1-probsum)/(n_classes-3)]*n_classes
            #print("Test1.ba")
            for j in range(3):
                #try:
                # Get the segment index (column) and the probability (value)
                seg_idx = int(event.ANN_segments[j])
                seg_prob = prob[j]
                # Safety check: ensure index is within the range
                if 0 <= seg_idx < n_classes:
                    tmp_X[seg_idx] = seg_prob
            #print("Test1.bb")
            #if any(x < 0 for x in tmp_X) or any(x > 1 for x in tmp_X) or not np.isclose(np.sum(tmp_X), 1.0):
            #    print(np.sum(tmp_X), tmp_X, chain.GetReadEntry(), chain.GetFile().GetName(), i, j, seg_idx, seg_prob, event.mc_run, event.mc_subrun, event.mc_nthEvtInFile)
            #    if(any(x < 0 for x in tmp_X)):
            #        print("Case1")
            #    if(any(x > 1 for x in tmp_X)):
            #        print("Case2")
            #    if(not np.isclose(np.sum(tmp_X), 1.0)):
            #        print("Case3: ", np.sum(tmp_X))
            #    continue #This shouldn't really happen
            #print("Test1.c")
            #if i >= n_entries_interest:
            #    break
            X_probs_arr.append(tmp_X)
            #if not np.isclose(np.sum(tmp_X), 1.0): #Rarely we have cases where due to rounding the probabilites sum to more than 1. This is typically when 1 plane is assigned essentiall 1 as probability
                #print(np.sum(tmp_X), chain.GetReadEntry(), chain.GetFile().GetName(), i, j, seg_idx, seg_prob, event.mc_run, event.mc_subrun, event.mc_nthEvtInFile)
            #print("Test1.d")
            # --- Fill Y (MC Only) ---
            if is_mc:
                # Get the truth label
                mdl = event.truth_vtx_module
                plane = event.truth_vtx_plane
                tgt_id = event.truth_targetID
                truth_val = get_truth_segment(mdl, plane, tgt_id)
                Y_arr.append(truth_val)
                #print("truth: ", truth_val, "\treco: ", np.argmax(tmp_X))
    #print("Test2")
    X_probs = np.array(X_probs_arr)
    Y = np.array(Y_arr)
    print(f"Data Loaded. Shape: {X_probs.shape}")
    print("Processing complete.")
    return X_probs, Y

def load_text_file(filepath, subsample_rate=1.0, cap_entries=None):
    """
    Because the raw DANN prediction files is HUGE, I need to use a buffered and chunked approach to prevent this from crashing due to RAM overuse on file loading
    Additionally I apply a cut here, I remove events that don't have an uncalibrated ML prediction since there's no point incuding them and contaminating the data we use to train our calibrator
    Also performs random subsampling
    """
    print(f"Reading large file: {filepath}")
    print(f"Subsampling rate: {subsample_rate * 100}%")
    
    labels_list = []
    probs_list = []
    
    import pandas as pd
    
    file_size = os.path.getsize(filepath)
    with open(filepath, 'r') as f:
        with tqdm(total=file_size, unit='B', unit_scale=True, desc="Loading Data") as pbar: #Using tqdm to have some idea of progress as this can take a while
            
            #Using chunksize to break up input file reading
            reader = pd.read_csv(f, header=None, chunksize=100000)  # Process 100k rows at a time
            
            for chunk in reader:
                raw_labels = chunk.iloc[:, LABEL_COL_IDX].values
                cnn_labels = chunk.iloc[:, LABEL_COL_IDX+1].values
                valid_mask = (cnn_labels > 0) # A mask to keep only the entries that have a CNN prediction - the CNN prediction is stored the in row after the truth plane

                
                # Using another mask to subsample to save RAM
                if subsample_rate < 1.0:
                    # Generate a random mask of the same size
                    # Works by creating an array of numbers randomly distrubuted between 0-1 equal in length to the number of text lines/entries in this chunk of data and then evalues each of these array values with < subsample_rate 
                    random_mask = np.random.rand(len(chunk)) < subsample_rate

                    final_mask = valid_mask & random_mask # Combined mask
                else: #If not subsampling can skip the np.random stuff
                    final_mask = valid_mask

                labels = raw_labels[final_mask].astype(np.int64) # Apply filter masks
                probs = chunk.iloc[final_mask, PROB_COL_START : PROB_COL_END + 1].values.astype(np.float32)
                np.set_printoptions(threshold=sys.maxsize)
                if len(labels) > 0:
                    labels_list.append(labels)
                    probs_list.append(probs)

                    #Just top 3
                    #top_k_idx = np.argpartition(probs, -3, axis=1)[:, -3:]
                    #probs_top3 = np.zeros_like(probs)
                    #row_idx = np.arange(probs.shape[0])[:, None]
                    #probs_top3[row_idx, top_k_idx] = probs[row_idx, top_k_idx]
                    #probs_list.append(probs_top3)

                    #Scales top 3
                    #top_k_idx = np.argpartition(probs, -3, axis=1)[:, -3:]
                    #probs_top3 = np.zeros_like(probs)
                    #row_idx = np.arange(probs.shape[0])[:, None]
                    #probs_top3[row_idx, top_k_idx] = probs[row_idx, top_k_idx]
                    #row_sums = probs_top3.sum(axis=1, keepdims=True)
                    #with np.errstate(divide='ignore', invalid='ignore'):
                    #    probs_top3 = probs_top3 / row_sums
                    #probs_list.append(probs_top3)

                    #Applies flat norm to non-top3 values
                    #top_k_idx = np.argpartition(probs, -3, axis=1)[:, -3:]
                    #row_idx = np.arange(probs.shape[0])[:, None]
                    #fill_value = 1e-9
                    #probs_filled = np.full_like(probs, fill_value)
                    #probs_filled[row_idx, top_k_idx] = probs[row_idx, top_k_idx]
                    #row_sums = probs_filled.sum(axis=1, keepdims=True)
                    #with np.errstate(divide='ignore', invalid='ignore'):
                    #    probs_filled = probs_filled / row_sums
                    #probs_list.append(probs_filled)

                    #Applies exponential decay to non-top-3 values centred on top1
                    #decay_rate = 0.5
                    #top3_idx = np.argpartition(probs, -3, axis=1)[:, -3:]
                    #top1_idx = np.argmax(probs, axis=1)[:, None] 
                    #row_idx = np.arange(probs.shape[0])[:, None]
                    #top3_values = probs[row_idx, top3_idx]
                    #current_top3_sum = top3_values.sum(axis=1, keepdims=True)
                    #residual_budget = np.maximum(0, 1.0 - current_top3_sum)
                    #col_indices = np.arange(probs.shape[1])
                    #distances = np.abs(col_indices - top1_idx)
                    #decay_map = np.exp(-decay_rate * distances)
                    #decay_map[row_idx, top3_idx] = 0
                    #decay_sums = decay_map.sum(axis=1, keepdims=True)
                    #with np.errstate(divide='ignore', invalid='ignore'):
                    #    scale_factors = residual_budget / decay_sums
                    #    scale_factors = np.nan_to_num(scale_factors)
                    #background_values = decay_map * scale_factors
                    #final_probs = background_values.copy()
                    #final_probs[row_idx, top3_idx] = top3_values
                    #probs_list.append(final_probs)

                
                pbar.update(f.tell() - pbar.n)
                #break

    print("Concatenating chunks...")
    if not labels_list:
        raise ValueError("No data found! Check your filter, subsample rate or input file.")

    if cap_entries is not None: #Setting a hard, numerical cap on the number of entries to keep, mainly used for Dircal training where the memory usage get's humongous really quickly
        labels_list = labels_list[:cap_entries]
        probs_list = probs_list[:cap_entries]
    y = np.concatenate(labels_list)
    X_probs = np.concatenate(probs_list)
    
    random_mask = np.random.rand(len(chunk)) < subsample_rate

    # Force garbage collection
    del labels_list
    del probs_list
    import gc
    gc.collect()
    
    print(f"Data Loaded. Shape: {X_probs.shape} (Subsampled)")
    return X_probs, y


def load_data_buffered(filepath, subsample_rate=500000, cap_entries=None, anuMode = False):
    #Open first line to check data input type
    first_line = ""
    with open(filepath, "r") as f:
        first_line = f.readline().strip()
    if ".root" in first_line:
        print("Detected a set of root files - Loading in appropriately")
        return load_root_file_optimized(filepath, subsample_rate)
        #return load_root_file_awk(filepath, 50000)
        #return load_root_file_rdf(filepath, subsample_rate)
        #return load_root_file(filepath, 50000)
    else:
        print("Assuming you are passing the raw prediction files - Loading in appropriately")
        return load_text_file(filepath, subsample_rate, cap_entries)



def calc_ECE(probs, y_true, n_ranks=3, n_bins=10):
    """
    Calculate ECE
    """
    print(f"Calculating ECE")

    # Get the indices of the top N predictions
    # argsort sorts ascending, so we take the last n_ranks columns and reverse
    top_preds_indices = np.argsort(probs, axis=1)[:, -n_ranks:][:, ::-1]
    ECE_out = []
    for rank in range(n_ranks):
        # Get the class index predicted at this rank
        rank_class_indices = top_preds_indices[:, rank]
        
        # Get the probability assigned to this rank
        # Use numpy advanced indexing to grab the prob for the specific class at this rank
        confidences = probs[np.arange(len(probs)), rank_class_indices]

        # Determine if this rank was the correct class
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


def plot_sms_calibration_matrix_old(calibrator, outpath, ax=None, cmap='viridis' ):
    """
    Plots the calibration matrix from a fitted probmetrics SMSCalibrator object 
    as a 2D heatmap.
    
    Parameters:
    - calibrator: A fitted SMSCalibrator instance.
    - ax: Optional matplotlib Axes object to plot on.
    - cmap: String representing the matplotlib colormap.
    """
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('agg')
    # 1. Extract the learned matrix. 
    if hasattr(calibrator, 'W'):
        matrix = calibrator.W
    else:
        raise AttributeError(
            "Could not find the matrix attribute on the SMSCalibrator object. "
            "Ensure the calibrator has been fitted via .fit() before plotting."
        )
    
    # Ensure standard format
    matrix = np.asarray(matrix)
    if matrix.ndim != 2:
        raise ValueError(f"Expected a 2D matrix, but got shape {matrix.shape}")

    # 2. Setup the plot
    if ax is None:
        fig, ax = matplotlib.pyplot.subplots(figsize=(8, 6))
    else:
        fig = ax.figure
        
    # 3. Create the heatmap
    cax = ax.imshow(matrix, cmap=cmap, aspect='auto')
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
    plt.clf()

def plot_sms_calibration_matrix(calibrator, outpath, ax=None, cmap='viridis' ):
    """
    Plots the calibration matrix from a fitted probmetrics SMSCalibrator object 
    as a 2D heatmap.
    
    Parameters:
    - calibrator: A fitted SMSCalibrator instance.
    - ax: Optional matplotlib Axes object to plot on.
    - cmap: String representing the matplotlib colormap.
    """
    import matplotlib.pyplot as plt
    W_full = calibrator.W
    weights = W_full[:, :-1]
    bias = W_full[:, -1]
    if hasattr(calibrator.ts.cal_, 'invtemp_'):
        print("Collapsing temp scaling into matrix...")
        alpha = calibrator.ts.cal_.invtemp_ 
        weights = weights * alpha
    W_fused = np.hstack([weights, bias.reshape(-1, 1)])
    matrix = W_fused.astype(np.float64)

    # Ensure standard format
    matrix = np.asarray(matrix)
    if matrix.ndim != 2:
        raise ValueError(f"Expected a 2D matrix, but got shape {matrix.shape}")

    # 2. Setup the plot
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    else:
        fig = ax.figure
        
    #from matplotlib.colors import LogNorm
    # 3. Create the heatmap
    cax = ax.imshow(matrix, cmap=cmap, aspect='auto', vmin=-0.05, vmax=0.15)#, norm=LogNorm())
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


import probmetrics.utils
def validate_probabilities_monkey_patch(X: np.ndarray) -> None:
    """
    Checks if X is a valid 2D array of probability distributions.
    Raises ValueError if the validation fails.
    Monkey patch to loosen probability sum requirement
    """
    if X.ndim != 2:
        raise ValueError(f"Input array must be 2D, but got {X.ndim} dimensions.")
    if np.any(X < 0) or np.any(X > 1):
        raise ValueError("All probability values must be within the [0, 1] range.")
    if not np.all(np.isclose(X.sum(axis=1), 1.0, atol=0.01)):
        raise ValueError("Each row of the input array must sum to 1 to be a valid probability distribution.")


# 2. Overwrite the function in the base module
probmetrics.utils.validate_probabilities = validate_probabilities_monkey_patch

def trainCalibrator(calibratorstr, rho_train, tau_train, lambda_intercept_train, lambda_diagonal_train, lambda_off_diagonal_train):
    print(f"Current hyper-parameter values: lambda_intercept = {lambda_intercept_train}, lambda_diagonal = {lambda_diagonal_train}, lambda_off_diagonal={lambda_off_diagonal_train}, rho = {rho_train}, tau={tau_train}")
    if calibratorstr == "sms":
        lambda_off_diagonal = 0.5
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
    print("================y_train==============================")
    print(f"Fitting Calibrator on {len(y_train)} samples (this may take time)...")
    print(f"X_train_probs: {X_train_probs.shape}")
    print(f"y_train: {y_train.shape}")
    calibrator.fit(X_train_probs, y_train)
    print(f"Finished training")

    test_calibrated = calibrator.predict_proba(X_test_probs)
    print("Model saved successfully.")
    print(f'test_calibrated len {len(test_calibrated)}')
    print(f'y_test len {len(y_test)}')    
    ece = calc_ECE(test_calibrated, y_test, n_ranks=3, n_bins=10)
    print('ECE: ', ece)

    SAVE_PATH = os.path.join(OUTPUT_DIR, f"{calibratorstr}_{lambda_intercept_train}_{lambda_diagonal_train}_{lambda_off_diagonal_train}_{rho_train}_{tau_train}_calibrator.pkl")
    print(f"Saving calibrated model to {SAVE_PATH}...")

    # This saves the entire calibrator object, including the learned matrix W
    joblib.dump(calibrator, SAVE_PATH)
    if calibratorstr == "sms":
        matrix_path = os.path.join(OUTPUT_DIR, f"{calibratorstr}_{lambda_intercept_train}_{lambda_diagonal_train}_{lambda_off_diagonal_train}_{rho_train}_{tau_train}_calibratorMatrix.png")
        plot_sms_calibration_matrix(calibrator, outpath = matrix_path)


    RESULTS_PATH = os.path.join(OUTPUT_DIR, f"{calibratorstr}_{lambda_intercept_train}_{lambda_diagonal_train}_{lambda_off_diagonal_train}_{rho_train}_{tau_train}_results.txt")
    with open(RESULTS_PATH, "a") as myfile:
        myfile.write(f"{ece[0]},{ece[1]},{ece[2]},{lambda_intercept_train},{lambda_diagonal_train},{lambda_off_diagonal_train},{rho_train},{tau_train}\n")
    if LOCALPATH is not None:
        print(f"Saving from grid node to local path {LOCALPATH}")
        os.system(f"ifdh cp {RESULTS_PATH} {LOCALPATH}")
        os.system(f"ifdh cp {matrix_path} {LOCALPATH}")
        os.system(f"ifdh cp {SAVE_PATH} {LOCALPATH}")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="SMS Calibrator. This scirpt is dual purpose. When no model is supplied, it calibrates it's own. And it also produces a number of relevant plots")
    parser.add_argument("-i", "--input", required=True, type=str, help="Path to text file containing raw DANN output")
    parser.add_argument("-r", "--rootfile", type=str, help="Path to root files")
    parser.add_argument("-o", "--output", type=str, help="Output base")
    parser.add_argument("-s", "--samplerate", default=0.1, type=float, help="Proportion of input files to sample, default is to 0.1")
    parser.add_argument("-sr", "--rootsamplerate", default=0.1, type=float, help="Proportion of input root files to sample, default is to 0.1")
    parser.add_argument("-c", "--calibrator", type=str, default="sms", help="Calibrator. Recommended options: sms or dircal")
    parser.add_argument("-l", "--localpath", type=str, help="The path to copy the outputs to (only needed if running on the grid)")
    parser.add_argument("-p", "--processes", default=2, type=int, help="The number of processes over which to split the job")
    parser.add_argument('--lambdas', nargs='+', help='The values to sample for the lambda variables')
    parser.add_argument('--rhos', nargs='+', help='The values to sample for the rho variables')
    parser.add_argument('--taus', nargs='+', help='The values to sample for the tau variables')
    parser.add_argument("--rhc",action='store_true', help="Flag determining if we're looking at a RHC/anti-nu playlist'")



    args = parser.parse_args()

    print("Beginning calibration")

    file_path = args.input
    calibratorstr = args.calibrator
    outbase = args.output
    samplerate = args.samplerate
    rootsamplerate = args.rootsamplerate
    rootfile = args.rootfile
    LOCALPATH = args.localpath
    isRHC = args.rhc
    if isRHC:
        print("Treating as an anti-neutrino sample")
    print(f"Local path: {LOCALPATH}")
    n_processes = args.processes
    if args.lambdas is not None:
        lambda_samples = list(map(float, args.lambdas))
        print(f"Using user specified lambda samples: {lambda_samples}")
    else:
        lambda_samples = [0.01, 0.1, 1, 10, 100]
        print(f"Using default lambda samples: {lambda_samples}")
    if args.rhos is not None:
        rho_samples = list(map(float, args.rhos))
        print(f"Using user specified rho samples: {rho_samples}")
    else:
        rho_samples = [0., 0.5, 1]
        print(f"Using default rho samples: {rho_samples}")
    if args.taus is not None:
        tau_samples = list(map(float, args.taus))
        print(f"Using user specified lambda samples: {tau_samples}")
    else:
        tau_samples = [0.5, 1, 1.5, 2]
        print(f"Using default lambda samples: {tau_samples}")
    
    justTrain = False #flag to determine if we are only training and not testing on the supplied data
    #Try to load data
    try:
        hardcap = None
        if calibratorstr == "dircal":
            hardcap = 10000
        X_probs, y = load_data_buffered(file_path,  subsample_rate=samplerate, cap_entries=hardcap, anuMode=isRHC)
    except Exception as e:
        print(f"Error loading data: {e}")
        exit()

    #print(f"Label range: {y.min()} to {y.max()}")
    print(f"Probability matrix shape: {X_probs.shape}")

    # 5. Standard SMS Calibration
    print("\n--- Starting Calibration ---")
    print("Initializing Probmetrics Calibrator...")
    


    if rootfile is None:
        print("Splitting data into Train and Test sets...")
        from sklearn.model_selection import train_test_split
        X_train_probs, X_test_probs, y_train, y_test = train_test_split(X_probs, y, test_size=0.2, random_state=42)
        del X_probs, y # Free up memory (delete the full original arrays which can be quite big and eat a lot of RAM)
        
        first_line = ""
        with open(file_path, "r") as f:
            first_line = f.readline().strip()
        if not ".root" in first_line:
            #Truncate input planes from text input to just the top 3 most probable to reproduce what we see in our root files
            top_k_idx = np.argpartition(X_test_probs, -3, axis=1)[:, -3:]
            probs_top3 = np.zeros_like(X_test_probs)
            row_idx = np.arange(X_test_probs.shape[0])[:, None]
            probs_top3[row_idx, top_k_idx] = X_test_probs[row_idx, top_k_idx]
            X_test_probs = probs_top3
    else:
        print("Found .root files")
        X_test_probs, y_test = load_root_file(rootfile, fill_mode=None, subsample_rate=rootsamplerate, anuMode=isRHC) #This is wrong
        X_train_probs=X_probs
        y_train = y

    intercept_range = lambda_samples
    diagonal_range = lambda_samples
    off_diagonal_range = lambda_samples
    rho_range = rho_samples
    tau_range = tau_samples
    threads = []
    total_vals = len(intercept_range)*len(diagonal_range)*len(off_diagonal_range)*len(rho_range)*len(tau_range)
    print(f"Total variations to try {total_vals}")
    for inter in intercept_range:
        for diag in diagonal_range:
            for off_diag in off_diagonal_range:
                for rho_test in rho_range:
                    for tau_test in tau_range:
                        #trainCalibrator(calibratorstr, rho_test, tau_test, inter, diag, off_diag,)
                        t = threading.Thread(target=trainCalibrator, args=(calibratorstr, rho_test, tau_test, inter, diag, off_diag,))
                        threads.append(t)
    #RESULTS_PATH = os.path.join(OUTPUT_DIR, f"{calibratorstr}_hyperparam_optimisation_results.txt")
    #if os.path.isfile(RESULTS_PATH):
    #    os.remove(RESULTS_PATH)
    #baselineECE = calc_ECE(X_test_probs, y_test)
    #with open(RESULTS_PATH, "w") as myfile:
    #    myfile.write("ece[0],ece[1],ece[2],lambda_intercept,lambda_diagonal,lambda_off_diagonal,rho,tau\n")
    #    myfile.write(f"{baselineECE[0]},{baselineECE[1]},{baselineECE[2]},BASELINE,BASELINE,BASELINE,BASELINE,BASELINE\n")
    count = 1
    len_threads = len(threads)
    from more_itertools import batched
    for thread_batch in batched(threads, n_processes):
        for t in thread_batch:
            print(f"Working on #{count} of {len_threads} hyper-paremeter variations")
            t.start()
        for t in thread_batch:
            t.join()
            print(f"Completed #{count} of {total_vals} hyper-paremeter variations")
            count=count+1
