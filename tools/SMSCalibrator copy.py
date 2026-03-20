#This script trains a structures matrix scaling (or other specified type) based calibrator for MINERvA DANN predictions
#Based on https://arxiv.org/pdf/2511.03685
import ROOT#, PlotUtils
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import log_loss
from tqdm import tqdm
import os, joblib, argparse

import sys
from probmetrics.calibrators import get_calibrator, SMSCalibrator, SVSCalibrator, Calibrator
import dirichletcal
from typing import Optional, List

import bisect

#import torch, torchmin


#PLOTUTILSROOT = os.path.expandvars("$PLOTUTILSROOT")
#ROOT.gROOT.SetBatch(1)
#plotter = PlotUtils.MnvPlotter()



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
RVec<double> map_probs(RVec<double> segments, RVec<double> probs, int n_classes) {
    double probSum = probs[0]+probs[1]+probs[2];
    double probRemainder = (probSum<1.0) ? (1.0 - probSum)/(n_classes-3) : 0;
    
    RVec<double> result(n_classes, probRemainder);
    for (size_t i = 0; i < 3; ++i) {
        int idx = (int)segments[i];
        if (idx >= 0 && idx < n_classes) {
            result[idx] = probs[i];
        }
    }
    return result;
}
""")

import numpy as np
import os


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

    # Case: zpos is smaller than or equal to the first element
    if idx == 0:
        return 0

    # Case: zpos is strictly larger than the last element
    if idx == len(PLANE_CENTRES):
        return len(PLANE_CENTRES) - 1

    # Case: zpos is between two elements
    val_above = PLANE_CENTRES[idx]
    val_below = PLANE_CENTRES[idx - 1]

    if abs(zpos - val_below) < abs(zpos - val_above):
        # The element BEFORE the iterator is closer
        return idx  # Equivalent to C++ logic: idx - 1 + 1
    else:
        # The element AT the iterator is closer
        return idx + 1 # Equivalent to C++ logic: idx + 1

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
    """
    Returns the Z position for a given segment ID. 
    Expects segment IDs between 1 and 214.
    """
    if 1 <= segment <= 214:
        return _ZPOS_LOOKUP[segment]
    
    # Acts as your fallback 'else' block
    return 0.0
get_zpos_from_segment_vec = np.vectorize(get_zpos_from_segment)

def passesCuts(event, anuMode = False):
    import math
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

def load_root_file_optimized(file_list_path, subsample_rate=1.0, anuMode=False):
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
    
    # Tell ROOT which branches to cache (Auto-learning is usually fine, but this is safer)
    chain.AddBranchToCache("*", True)
    try:
        ROOT.RDF.Experimental.AddProgressBar(df)
    except AttributeError:
        pass
    # 1. SUBSAMPLING (Zero I/O Cost)
    subsample_rate*=prescaledownrate
    if subsample_rate < 1.0:
        df = df.Filter(f"gRandom->Rndm() < {subsample_rate}", "Subsampling")

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

    # 5. FINAL CALCULATION
    # Check the sum last, as it requires reading the ANN_plane_probs vector
    print("Doing prob check")
    df = df.Filter("Sum(ANN_plane_probs) > 0", "Sum Prob Check")
    print("Finished prob check")

    n_classes = 215 
    #print("Producing X_probs")
    #df = df.Define("X_probs", f"map_probs(ANN_segments, ANN_plane_probs, {n_classes})")
    #print("Produced X_probs")

    # Explode the vector into scalar columns
    df = df.Define("X_probs_vec", f"map_probs(ANN_segments, ANN_plane_probs, {n_classes})")
    cols_to_read = ["X_probs_vec"]
    for i in range(n_classes):
        col_name = f"prob_{i}"
        # Extract the i-th element. The casting (float) ensures it's a simple type.
        df = df.Define(col_name, f"(float)X_probs_vec[{i}]")
        cols_to_read.append(col_name)
    #print(f"cols_to_read {cols_to_read}")
    print("Produced X_probs")


    # This report will tell you exactly how many events passed each stage!
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
    
    # Run the graph
    #print(f'{df.Count().GetValue()} events in the RDataFrame')
    print("Processing events in RDataFrame...")
    data_dict = df.AsNumpy(cols_to_read)
    print("Here1")
    report = df.Report()
    report.Print()
    print("Here2")

    #print(data_dict)
    #print(data_dict["X_probs_vec"])

    #X_probs = np.vstack(list(data_dict["X_probs_vec"])).astype(np.float64)
    flat_array = np.concatenate(data_dict["X_probs_vec"])
    print("Here2")

    # 2. Reshape it back into a 2D matrix (N_events, N_classes). 
    # The '-1' tells NumPy to automatically calculate the number of columns (e.g., 215).
    n_events = len(data_dict["X_probs_vec"])
    X_probs = flat_array.reshape(n_events, -1)
    print("Here4")

    # 3. Only cast to float64 if your downstream calibrator absolutely requires it.
    # If it can handle float32, delete this next line to save 50% of your RAM.
    X_probs = X_probs.astype(np.float32)

    #print(prob_cols)
    #X_probs = np.stack(prob_cols, axis=1)

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

    return X_probs, np.array(Y)


def MINERvA-101-Cross-Section/tools/SMSCalibrator.pyload_root_file(file_list_path, subsample_rate=1.0):
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
    
    fill_mode = None
    # 1. Initialize TChain
    chain = ROOT.TChain("MasterAnaDev")
    
    # Read the .txt file and add files to the chain
    print(f"Reading root file list: {file_list_path}")
    with open(file_list_path, 'r') as f:
        for line in f:
            root_file = line.strip()
            if root_file:  # distinct from empty lines
                chain.Add(root_file)

    n_entries = chain.GetEntries()
    print(f"Total events found: {n_entries}")
    if subsample_rate != 1.0:
        print(f"Subsampling rate: {subsample_rate * 100}%")
    is_mc = False
    if chain.GetBranch("truth_vtx_module"):
        is_mc = True
    print(f"Processing type: {'MC' if is_mc else 'Data'}")
    print(f"Using fill_mode {fill_mode}")
    # 2. Initialize Output Containers
    # Initialize X_probs with zeros. Shape: (Total Events, N_classes)
    n_classes = PROB_COL_END-PROB_COL_START+1
    X_probs_arr = []
    Y_arr = []
    import random
    # 3. Loop over events
    # We enumerate to keep track of the event index 'i' for filling X_probs
    for i, event in enumerate(tqdm(chain, total=n_entries, unit="event", desc="Loading Data")):
        if random.random() > subsample_rate:
            continue
        if passesCuts(event):
            continue
        prob = [float(event.ANN_plane_probs[0]), float(event.ANN_plane_probs[1]), float(event.ANN_plane_probs[2])]
        probsum = prob[0]+prob[1]+prob[2]
        if fill_mode == 2 or probsum > 1:
            # With fill_mode 2 we're scaling the probabilites of the 3 planes were given to fill the full probaility space
            # Rarely we have cases where due to rounding the probabilites sum to more than 1. This is typically when 1 plane is assigned essentially 1 as probability, in this case we're scaling down to ensure normalisation
            prob[0] = prob[0] / probsum
            prob[1] = prob[1] / probsum
            prob[2] = prob[2] / probsum
            probsum = prob[0]+prob[1]+prob[2]
        if int(event.hasMLPrediction) == 1 and probsum > 0:
            # --- Fill X_probs ---
            # We assume ANN_segments and ANN_plane_probs are arrays of length 3
            # Note: ANN_segments are doubles, so we cast to int for array indexing.
            tmp_X = [0]*n_classes
            if fill_mode == 1 and probsum<1: # In fill type 1 we distribute remaining probability evenly amongst remaining classes
                tmp_X = [(1-probsum)/(n_classes-3)]*n_classes
            if fill_mode == 3 and probsum<1: # In fill type 3 we put all remaining probability density in an overflow class
                tmp_X[n_classes-1] = (1-probsum)
            for j in range(3):
                #try:
                # Get the segment index (column) and the probability (value)
                seg_idx = int(event.ANN_segments[j])
                seg_prob = prob[j]
                # Safety check: ensure index is within the range
                if 0 <= seg_idx < n_classes:
                    tmp_X[seg_idx] = seg_prob
            #if any(x < 0 for x in tmp_X) or any(x > 1 for x in tmp_X) or not np.isclose(np.sum(tmp_X), 1.0):
            #    print(np.sum(tmp_X), tmp_X, chain.GetReadEntry(), chain.GetFile().GetName(), i, j, seg_idx, seg_prob, event.mc_run, event.mc_subrun, event.mc_nthEvtInFile)
            #    if(any(x < 0 for x in tmp_X)):
            #        print("Case1")
            #    if(any(x > 1 for x in tmp_X)):
            #        print("Case2")
            #    if(not np.isclose(np.sum(tmp_X), 1.0)):
            #        print("Case3: ", np.sum(tmp_X))
            #    continue #This shouldn't really happen
            
            #if i >= n_entries_interest:
            #    break
            X_probs_arr.append(tmp_X)
            #if not np.isclose(np.sum(tmp_X), 1.0): #Rarely we have cases where due to rounding the probabilites sum to more than 1. This is typically when 1 plane is assigned essentiall 1 as probability
                #print(np.sum(tmp_X), chain.GetReadEntry(), chain.GetFile().GetName(), i, j, seg_idx, seg_prob, event.mc_run, event.mc_subrun, event.mc_nthEvtInFile)

            # --- Fill Y (MC Only) ---
            if is_mc:
                # Get the truth label
                mdl = event.truth_vtx_module
                plane = event.truth_vtx_plane
                tgt_id = event.truth_targetID
                truth_val = get_truth_segment(mdl, plane, tgt_id)
                Y_arr.append(truth_val)
                #print("truth: ", truth_val, "\treco: ", np.argmax(tmp_X))
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


def load_data_buffered(filepath, subsample_rate=1.0, cap_entries=None, fill_mode = 2):
    #Open first line to check data input type
    first_line = ""
    with open(filepath, "r") as f:
        first_line = f.readline().strip()
    if ".root" in first_line:
        print("Detected a set of root files - Loading in appropriately")
        #return load_root_file(filepath, subsample_rate)
        return load_root_file_optimized(filepath, subsample_rate)
    else:
        print("Assuming you are passing the raw prediction files - Loading in appropriately")
        return load_text_file(filepath, subsample_rate, cap_entries)

# ==========================================
#This ended up not being needed, the probmetrics package actually does this implicitly
# ==========================================
def probs_to_logits(probs, epsilon=1e-15):
    """
    Inverts Softmax to retrieve 'Pseudo-Logits'.
    x = ln(p)
    """
    print("Converting probabilities to pseudo-logits...")
    probs_clipped = np.clip(probs, epsilon, 1.0) # Clipping probabilities to avoid ln(0)
    logits = np.log(probs_clipped)
    return logits

def plot_multirank_reliability(probs, y_true, title, filename, n_ranks=3, n_bins=10):
    """
    Plots N reliability diagrams for the top N predicted ranks.
    """
    print(f"Generating reliability diagram: {title}")
    
    plt.figure(figsize=(6 * n_ranks, 5))
    
    # Get the indices of the top N predictions
    # argsort sorts ascending, so we take the last n_ranks columns and reverse
    top_preds_indices = np.argsort(probs, axis=1)[:, -n_ranks:][:, ::-1]
    for rank in range(n_ranks):
        plt.subplot(1, n_ranks, rank + 1)
        
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
    plt.savefig(os.path.join(OUTPUT_DIR, filename))
    plt.close()

    print("Finished Plot Generation")

def plot_diff_distribution(uncalib_probs, calib_probs, avg = False, n=1, num_classes=None, ax=None, label=None, color=None, save_filename=None, title=None):
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

    # Save file
    if save_filename:
        print(f"Saving plot to {save_filename}...")
        ax.figure.savefig(save_filename, bbox_inches='tight', dpi=300)
        
    return ax


def plot_weighted_avg_plane_distribution(probs, avg = False, num_classes=None, ax=None, label=None, color=None, save_filename=None, title=None):
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

    mc_uncalib =mc_uncalib[:, :173]
    mc_calib =mc_uncalib[:, :173]
    data_uncalib =mc_uncalib[:, :173]
    data_calib =mc_uncalib[:, :173]

    plot_prob_distributions(mc_uncalib, n=n, num_bins=num_bins, ax=ax_main, color='#d2c271', label="MC Uncalibrated")
    plot_prob_distributions(mc_calib, n=n, num_bins=num_bins, ax=ax_main, color='#CA5454', label="MC Calibrated")
    plot_prob_distributions(data_uncalib, n=n, num_bins=num_bins, ax=ax_main, color='#40B0A6', label="Data Uncalibrated")
    plot_prob_distributions(data_calib, n=n, num_bins=num_bins, ax=ax_main, color='#9285d5', label="Data Calibrated")

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

def plot_cdf_data_vs_mc(mc_uncalib, mc_calib, data_uncalib, data_calib, n=1, num_bins=10, save_filename=None, title=None):
    
    cdf_ax = plot_prob_cdf(mc_uncalib, n=n, num_bins=num_bins, color='#d2c271', label="MC Uncalibrated")
    cdf_ax = plot_prob_cdf(mc_calib, n=n, num_bins=num_bins, color='#CA5454', ax = cdf_ax, save_filename=save_filename, label="MC Calibrated")
    cdf_ax = plot_prob_cdf(data_uncalib, n=n, num_bins=num_bins, color='#40B0A6', ax = cdf_ax, label="Data Uncalibrated")
    cdf_ax = plot_prob_cdf(data_calib, n=n, num_bins=num_bins, color='#9285d5', ax = cdf_ax, save_filename=save_filename, label="Data Calibrated")
    cdf_ax.set_xlim(0, 1/n)

    suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(n if n < 20 else n % 10, 'th')
    if title is None:
        cdf_ax.set_title(f'Probability Distribution of the {n}{suffix} Most Probable Plane')

    if save_filename:
        print(f"Saving plot to {save_filename}...")
        plt.savefig(save_filename, bbox_inches='tight', dpi=300)


def plot_prob_cdf(probs, n=1, num_bins=10, ax=None, label=None, color=None, save_filename=None, title=None):
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
    ax.hist(nth_prob, bins=bins, histtype='step', cumulative=True, linewidth=2, label=label, color=color, weights=weights)
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

def calc_ECE_by_segment_old(probs, y_true, segment, n_bins=10):

    # Get the indices of the top N predictions
    # argsort sorts ascending, so we take the last n_ranks columns and reverse
    confidences = probs[:,segment]

    # Determine if this rank was the correct class
    accuracies = (segment == y_true).astype(float)
    
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
    return ece

def calc_ECE_by_segment(probs, y_true, segment, n_bins=10, threshold=0.01):
    confidences = probs[:, segment]
    is_true_class = (y_true == segment).astype(float)
    
    # 1. Drop the trivial true negatives
    active_mask = confidences > threshold
    if not np.any(active_mask):
        return np.nan
        
    conf_active = confidences[active_mask]
    acc_active = is_true_class[active_mask]
    
    # 2. Bin only the active probabilities
    # Note: Bins now start at the threshold, not 0
    bin_boundaries = np.linspace(threshold, 1, n_bins + 1)
    bin_lowers = bin_boundaries[:-1]
    bin_uppers = bin_boundaries[1:]
    
    ece = 0.0
    for bin_lower, bin_upper in zip(bin_lowers, bin_uppers):
        in_bin = (conf_active > bin_lower) & (conf_active <= bin_upper)
        prop_in_bin = in_bin.mean()
        
        if prop_in_bin > 0:
            accuracy = acc_active[in_bin].mean()
            avg_conf = conf_active[in_bin].mean()
            ece += np.abs(avg_conf - accuracy) * prop_in_bin
            
    return ece

def calc_ECE_by_segment_top1(probs, y_true, segment, n_bins=10):
    # 1. Find what the model actually predicted for all events
    predictions = np.argmax(probs, axis=1)
    
    # 2. Isolate ONLY the events where the model chose your target segment
    class_mask = (predictions == segment)
    
    # If the model never predicted this class, ECE is undefined
    if not np.any(class_mask):
        return np.nan 
        
    # 3. Get the confidences and accuracies ONLY for those events
    confidences = np.max(probs[class_mask], axis=1)
    accuracies = (y_true[class_mask] == segment).astype(float)
    
    # --- Standard Binning ---
    bin_boundaries = np.linspace(0, 1, n_bins + 1)
    bin_lowers = bin_boundaries[:-1]
    bin_uppers = bin_boundaries[1:]
    
    ece = 0.0
    
    for bin_lower, bin_upper in zip(bin_lowers, bin_uppers):
        if bin_lower == 0.0:
            in_bin = (confidences >= bin_lower) & (confidences <= bin_upper)
        else:
            in_bin = (confidences > bin_lower) & (confidences <= bin_upper)
            
        prop_in_bin = in_bin.mean()
        
        if prop_in_bin > 0:
            accuracy = accuracies[in_bin].mean()
            avg_conf = confidences[in_bin].mean()
            ece += np.abs(avg_conf - accuracy) * prop_in_bin

    return ece

""" def plot_sms_calibration_matrix_old(calibrator, outpath, ax=None, cmap='viridis' ):
    
    Plots the calibration matrix from a fitted probmetrics SMSCalibrator object 
    as a 2D heatmap.
    
    Parameters:
    - calibrator: A fitted SMSCalibrator instance.
    - ax: Optional matplotlib Axes object to plot on.
    - cmap: String representing the matplotlib colormap.
    
    import matplotlib.pyplot as plt
    # 1. Extract the learned matrix. 
    if hasattr(calibrator, 'W'):
        matrix = calibrator.W
    else:
        raise AttributeError(
            "Could not find the matrix attribute on the SMSCalibrator object. "
            "Ensure the calibrator has been fitted via .fit() before plotting."
        )
    #import pprint
    #pprint.pprint(calibrator.__dict__)
    classes = calibrator.classes_
    if hasattr(calibrator.ts.calibrator, 'invtemp_'):
        print("Adding const scale to matrix")
        alpha = calibrator.ts.calibrator.invtemp_ 
        I = torch.eye(len(classes))
        I = I * alpha
        matrix = matrix+I

    # Ensure standard format
    matrix = np.asarray(matrix)
    if matrix.ndim != 2:
        raise ValueError(f"Expected a 2D matrix, but got shape {matrix.shape}")

    # 2. Setup the plot
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
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
    print(f"Plotted calibration matrix {outpath}")
    return fig, ax
 """

#BEGIN - WORK IN PROGESS
# Below are early steps taken to implement a diagonal Dirichlet calibration as an alternative to SMS calibration

# A local,  version of DiagonalDirichletCalibrator since versioning issues mean we're stuck with an older version with a known bug
# so we're here implementing the fixed version here ourselves
# In future, if we get access to a newer python version, enabling an upgrade to a newer dirichletcal version, this may be removed
class DiagonalDirichletCalibrator(dirichletcal.calib.fulldirichlet.FullDirichletCalibrator):
    def fit(self, X, y, X_val=None, y_val=None, *args, **kwargs):

        self.weights_ = self.weights_init

        if X_val is None:
            X_val = X.copy()
            y_val = y.copy()

        _X = np.copy(X)
        _X = np.log(dirichletcal.utils.clip_for_log(_X))
        _X_val = np.copy(X_val)
        _X_val = np.log(dirichletcal.utils.clip_for_log(X_val))

        self.calibrator_ = dirichletcal.calib.multinomial.MultinomialRegression(
            method='Diag', reg_lambda=self.reg_lambda, reg_mu=self.reg_mu,
            reg_norm=self.reg_norm, ref_row=self.ref_row,
            optimizer=self.optimizer)
        print(f"Probability matrix shape: {_X.shape}")
        print(f"y shape: {y.shape}")
        self.calibrator_.fit(_X, y, labels = np.arange(0, 215), *args, **kwargs)
        self.final_loss = log_loss(y_val, self.calibrator_.predict_proba(_X_val))

        return self

# Manual extension of the probmetrics or probmetrics style classes for diagonal Dirichlet calibration
class DiagDirichletCalibrator(Calibrator):
    def __init__(
        self,
        n_cv: int = 0,
        use_odir: bool = False,
        reg_lambda: float = 0.0,
        reg_mu: Optional[float] = None,
        reg_lambda_grid: Optional[List[float]] = None,
        reg_mu_grid: Optional[List[float]] = None,
    ):
        super().__init__()

        self.n_cv = n_cv
        self.use_odir = use_odir
        self.reg_lambda = reg_lambda
        self.reg_mu = reg_mu
        self.reg_lambda_grid = reg_lambda_grid
        self.reg_mu_grid = reg_mu_grid

    def _fit_impl(self, X: np.ndarray, y: np.ndarray) -> None:
        self.cal_ = DiagonalDirichletCalibrator()
        self.cal_.l2 = self.reg_lambda
        self.cal_.fit(X, y)

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        return self.cal_.predict_proba(X)



""" class SVSCalibratorMSE(Calibrator):
    
    Multiclass post-hoc calibration with a structured scaling scheme $softmax((aI + diag(v))x + b)$ on the logits x.

    Numba functions are warmed-up the first time the class is initialized with the 'saga' optimizer.

    A penalty (one of [None, "mcp", "lasso", "ridge"]) is applied to the intercept (b) and vector scaling (v) parameters,
    with respective regularization strength lambda_intercept*(k**rho)/(n**tau) and lambda_diagonal*(k**rho)/(n**tau).

    Instead of fitting the global scaling parameter jointly with the other parameters, the logits are pre-processed with
    temperature scaling (with no regularization) and a is then fixed to a=1.

    [!] solver = 'bfgs' only supports 'ridge' penalty.
    [!] For solver = 'bfgs', max_iter and tol are ignored and default torchmin parameters are used.
    

    _warmed_up_saga = False
    _ALLOWED_PENALTIES = [None, "mcp", "lasso", "ridge"]
    _ALLOWED_OPTIMIZERS = ["saga", "bfgs"]

    def __init__(
        self,
        penalty: str = "ridge",
        rho: float = 1.0,
        tau: float = 1.0,
        lambda_intercept: float = 1.0,
        lambda_diagonal: float = 1.0,
        opt: str = "bfgs",
        max_iter: int = 500,
        tol: float = 1e-5,
        print_init_info: bool = True,
    ):
        super().__init__()

        if penalty not in self._ALLOWED_PENALTIES:
            raise ValueError(
                f"Penalty '{penalty}' not recognized. Choose from {self._ALLOWED_PENALTIES}."
            )

        if opt not in self._ALLOWED_OPTIMIZERS:
            raise ValueError(
                f"Optimizer '{opt}' not recognized. Choose from {self._ALLOWED_OPTIMIZERS}."
            )

        if opt == "bfgs" and penalty != "ridge":
            raise ValueError(
                "The 'bfgs' optimizer only supports the 'ridge' penalty, use 'saga' instead."
            )
        self.penalty = penalty
        self.rho = rho
        self.tau = tau
        self.lambda_intercept = lambda_intercept
        self.lambda_diagonal = lambda_diagonal
        self.max_iter = max_iter
        self.tol = tol
        self.opt = opt
        self.print_init_info = print_init_info

        self.v = None
        self.b = None
        self.ts = None

    def _get_logits(self, X: np.ndarray) -> tuple[np.ndarray, int, int]:
        from probmetrics.utils import multiclass_probs_to_logits
        logits = multiclass_probs_to_logits(X)
        return logits
    def _fit_impl(self, X: np.ndarray, y: np.ndarray) -> None:
        #validate_probabilities(X)

        n, k = X.shape

        self.ts = get_calibrator("ts-mix")
        self.ts.fit(X, y)
        scaled_probs = self.ts.predict_proba(X)

        reg_intercept = self.lambda_intercept * (k**self.rho) / (n**self.tau)
        reg_diagonal = self.lambda_diagonal * (k**self.rho) / (n**self.tau)

        logits = self._get_logits(scaled_probs)

        if self.opt == "bfgs":
            self.v, self.b = self._fit_bfgs(
                torch.tensor(logits),
                torch.tensor(y),
                num_classes=k,
                reg_intercept=reg_intercept,
                reg_diagonal=reg_diagonal,
            )
    def _fit_bfgs(
        self,
        logits: torch.Tensor,
        y: torch.Tensor,
        num_classes,
        reg_intercept,
        reg_diagonal,
    ):
        K = num_classes

        initial_params = torch.zeros(2 * K, device=logits.device, dtype=logits.dtype)

        def loss(params):
            v_delta = params[:K]
            b = params[K:]

            # 1. Scale the logits (same as before)
            scaled_logits = logits * (1.0 + v_delta) + b

            # 2. Convert scaled logits to probabilities
            probs = torch.nn.functional.softmax(scaled_logits, dim=-1)

            # 3. Convert target indices 'y' to one-hot encoding
            # Ensure it matches the shape of probs and is a float tensor
            y_onehot = torch.nn.functional.one_hot(y, num_classes=logits.shape[-1]).to(probs.dtype)

            # 4. Calculate Mean Squared Error (Brier Score)
            mse = torch.nn.functional.mse_loss(probs, y_onehot)

            r_b = reg_intercept * b.pow(2).sum()
            r_v = reg_diagonal * v_delta.pow(2).sum()

            return mse + r_b + r_v

        method = "l-bfgs" if initial_params.numel() > 1000 else "bfgs"
        res = torchmin.minimize(loss, initial_params, method=method)

        v_delta_final = res.x[:K]
        b_final = res.x[K:]

        return (
            1.0 + v_delta_final
        ).detach().cpu().numpy(), b_final.detach().cpu().numpy()

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        #validate_probabilities(X)

        if self.v is None:
            raise RuntimeError(
                "Calibrator has not been fitted yet. Call fit() before predict_proba()."
            )

        scaled_probs = self.ts.predict_proba(X)
        logits = self._get_logits(scaled_probs)
        calibrated_logits = self.v.reshape(1, -1) * logits + self.b.reshape(1, -1)
        calibrated_logits -= np.max(calibrated_logits, axis=1, keepdims=True)
        exp_logits = np.exp(calibrated_logits)
        return exp_logits / np.sum(exp_logits, axis=1, keepdims=True)
 """




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="SMS Calibrator. This scirpt is dual purpose. When no model is supplied, it calibrates it's own. And it also produces a number of relevant plots")
    parser.add_argument("-i", "--input", required=True, type=str, help="Path to text file containing raw DANN output")
    parser.add_argument("-m", "--model", required=False, type=str, help="Path to pre-existing model")
    parser.add_argument("-o", "--output", type=str, help="Output directory")
    parser.add_argument("-d", "--data", type=str, help="Input data file")
    parser.add_argument("-s", "--samplerate", default=0.1, type=float, help="Proportion of input files to sample, default is to 0.1")
    parser.add_argument("--split_data",action='store_true', help="Flag determining if we split our data into training and testing")
    parser.add_argument("-c", "--calibrator", type=str, default="sms", help="Calibrator. Recommended options: sms or dircal")
    parser.add_argument("-f", "--fillmode", type=int, help="Fill mode for filling probability space for .root inputs")

    args = parser.parse_args()

    print("Beginning calibration")

    file_path = args.input
    model = args.model
    datafile = args.data
    calibratorstr = args.calibrator
    outbase = args.output
    samplerate = args.samplerate
    fillmode = args.fillmode
    
    justTrain = False #flag to determine if we are only training and not testing on the supplied data
    #Try to load data
    try:
        hardcap = None
        if calibratorstr == "dircal":
            hardcap = 10000
        X_probs, y = load_data_buffered(file_path,  subsample_rate=samplerate, cap_entries=hardcap, fill_mode = fillmode)
    except Exception as e:
        print(f"Error loading data: {e}")
        exit()

    #print(f"Label range: {y.min()} to {y.max()}")
    print(f"Probability matrix shape: {X_probs.shape}")

    # 5. Standard SMS Calibration
    print("\n--- Starting Calibration ---")
    print("Initializing Probmetrics Calibrator...")
    

    #If an existing calibrator is provided, use that and don't waste a lot of time retraining a new one
    #Otherwise, train one
    if model is None:
        print("No pretrained calibration model provided, making one using the supplied data")
        #Split Data. Note to self - maybe we shouldn't do this. Maybe we should train on a whole playlist and then apply it to other playlists
        if args.split_data:
            print("Splitting data into Train and Test sets...")
            X_train_probs, X_test_probs, y_train, y_test = train_test_split(
            X_probs, y, test_size=0.2, random_state=42
            )
            del X_probs, y # Free up memory (delete the full original arrays which can be quite big and eat a lot of RAM)
        else:
            print("You haven't provided a premade calibrator and also not set the --split-data flag to split your data into test/train sets. So I will assume that you don't wish to test on this data, and instead train on all of the supplied data.")
            justTrain = True
            X_train_probs = X_probs
            y_train = y
        #print("X_test_probs: ", X_test_probs)
        #print("y_test: ", y_test)
        #print("y_train: ", y_train)
        #print(f"X test max: {max(X_test_probs[0])}")
        
        if calibratorstr == "sms":
            lambda_off_diagonal = 1
            calibrator = SMSCalibrator (
            penalty = 'ridge' ,

            rho = 0.5 ,
            tau = 0.5 ,
            lambda_intercept = 1 ,
            lambda_diagonal = 1 ,
            lambda_off_diagonal = 0.1)
            print("lambda_off_diagonal: ", lambda_off_diagonal)
        elif calibratorstr == "svs":
            calibrator = SVSCalibrator (
            penalty = 'ridge' ,
            rho = 1.0 ,
            tau = 1.0 ,
            lambda_intercept = 0 ,
            lambda_diagonal = 0)
        elif calibratorstr == "ts":
            calibrator = get_calibrator("temp-scaling") 
        elif calibratorstr == "dircal":
            calibrator = DiagDirichletCalibrator(
                reg_lambda = 1.0
            )
            """calibrator = DiagDirichletCalibrator (
            n_cv = 0,
            use_odir = True,
            reg_lambda = 0.0,
            reg_mu = None,
            reg_lambda_grid = None,
            reg_mu_grid = None) """
            print("Setting ODIR for diriclet calibrator")
        else:
            calibrator = get_calibrator(calibratorstr)


        print("==============================================")
        print(f"Using {calibratorstr} calibrator")
        print("==============================================")
        print(f"Fitting Calibrator on {len(y_train)} samples (this may take time)...")
        print(f"X_train_probs: {X_train_probs.shape}")
        print(f"y_train: {y_train.shape}")
        

        out_of_bounds_mask = (X_train_probs < 0) | (X_train_probs > 1)

        # 2. Extract the actual values
        bad_values = X_train_probs[out_of_bounds_mask]

        # 3. Check if there are any bad values and print them
        if bad_values.size > 0:
            print(f"Found {bad_values.size} values outside the [0, 1] range!")
            
            # Optional but highly recommended: get the exact row and column indices
            row_indices, col_indices = np.where(out_of_bounds_mask)
            
            for r, c, val in zip(row_indices, col_indices, bad_values):
                print(f"Row {r}, Col {c}: {val}")
        else:
            print("All values are safely within the 0 to 1 range.")
    
        row_sums = np.sum(X_train_probs, axis=1)

        # 2. Find where the sum is NOT close to 1.0 (using a standard tolerance)
        # atol=1e-5 means it flags anything outside [0.99999, 1.00001]
        bad_rows_mask = ~np.isclose(row_sums, 1.0, atol=0.01)

        # 3. Get the specific row index numbers where this is True
        bad_indices = np.where(bad_rows_mask)[0]

        # 4. Loop ONLY over the bad rows to print them
        if bad_indices.size > 0:
            print(f"Found {bad_indices.size} rows that do not sum to 1.0!")
            for idx in bad_indices:
                print(f"--- Row {idx} ---")
                print(f"Sum: {row_sums[idx]}")
                print(f"Data: {X_train_probs[idx]}\n")
        else:
            print("All rows successfully sum to 1.0 (within tolerance).")


        calibrator.fit(X_train_probs, y_train)
        print(f"Finished training")
        if args.split_data:
            test_probs = X_test_probs
            y_probs = y_test
    else:
        print(f"Loading model from {model}...")
        calibrator = joblib.load(model)
        test_probs = X_probs
        y_probs = y
        calibratorstr = model.split("/")[-1]
        calibratorstr = calibratorstr.split("_")[0]
        print("==============================================")
        print(f"Using premade {calibratorstr} calibrator")
        print("==============================================")


    #Testing and plotting
    if not justTrain:
        print("Calibrating Test Data...")
        import matplotlib.pyplot as plt #We only need plotting functionality in this code block so lazy-loading it here rather than at the top
        test_calibrated = calibrator.predict_proba(test_probs)
        np.set_printoptions(threshold=sys.maxsize)
        #for i in range(len(test_calibrated)):
        #    print("Uncalibrated: ",np.argmax(test_probs[i]))
        #    print("Calibrated: ",np.argmax(test_calibrated[i]))
        #    print("Truth: ", y_probs[i])
        #Plot uncalibrated reliability



        print("Pre-calibration target ECEs")
        print(f"Target 1 ECE {calc_ECE_by_segment(test_probs, y_probs, 9)}")
        print(f"Target 2 ECE {calc_ECE_by_segment(test_probs, y_probs, 18)}")
        print(f"Target 3 ECE {calc_ECE_by_segment(test_probs, y_probs, 27)}")
        print(f"Water ECE {calc_ECE_by_segment(test_probs, y_probs, 36)}")
        print(f"Target 4 ECE {calc_ECE_by_segment(test_probs, y_probs, 45)}")
        print(f"Target 5 ECE {calc_ECE_by_segment(test_probs, y_probs, 50)}")
        print(f"Segment 22 ECE {calc_ECE_by_segment(test_probs, y_probs, 22)}")
        print(f"Segment 40 ECE {calc_ECE_by_segment(test_probs, y_probs, 40)}")
        print(f"Segment 70 ECE {calc_ECE_by_segment(test_probs, y_probs, 70)}")
        print(f"Segment 90 ECE {calc_ECE_by_segment(test_probs, y_probs, 90)}")

        print("Post-calibration target ECEs")
        print(f"Target 1 ECE {calc_ECE_by_segment(test_calibrated, y_probs, 9)}")
        print(f"Target 2 ECE {calc_ECE_by_segment(test_calibrated, y_probs, 18)}")
        print(f"Target 3 ECE {calc_ECE_by_segment(test_calibrated, y_probs, 27)}")
        print(f"Water ECE {calc_ECE_by_segment(test_calibrated, y_probs, 36)}")
        print(f"Target 4 ECE {calc_ECE_by_segment(test_calibrated, y_probs, 45)}")
        print(f"Target 5 ECE {calc_ECE_by_segment(test_calibrated, y_probs, 50)}")
        print(f"Segment 22 ECE {calc_ECE_by_segment(test_calibrated, y_probs, 22)}")
        print(f"Segment 40 ECE {calc_ECE_by_segment(test_calibrated, y_probs, 40)}")
        print(f"Segment 70 ECE {calc_ECE_by_segment(test_calibrated, y_probs, 70)}")
        print(f"Segment 90 ECE {calc_ECE_by_segment(test_calibrated, y_probs, 90)}")
        print("Top1-------------")
        print("Pre-calibration target ECEs")
        print(f"Target 1 ECE {calc_ECE_by_segment_top1(test_probs, y_probs, 9)}")
        print(f"Target 2 ECE {calc_ECE_by_segment_top1(test_probs, y_probs, 18)}")
        print(f"Target 3 ECE {calc_ECE_by_segment_top1(test_probs, y_probs, 27)}")
        print(f"Water ECE {calc_ECE_by_segment_top1(test_probs, y_probs, 36)}")
        print(f"Target 4 ECE {calc_ECE_by_segment_top1(test_probs, y_probs, 45)}")
        print(f"Target 5 ECE {calc_ECE_by_segment_top1(test_probs, y_probs, 50)}")
        print(f"Segment 22 ECE {calc_ECE_by_segment_top1(test_probs, y_probs, 22)}")
        print(f"Segment 40 ECE {calc_ECE_by_segment_top1(test_probs, y_probs, 40)}")
        print(f"Segment 70 ECE {calc_ECE_by_segment_top1(test_probs, y_probs, 70)}")
        print(f"Segment 90 ECE {calc_ECE_by_segment_top1(test_probs, y_probs, 90)}")
        print("Post-calibration target ECEs")
        print(f"Target 1 ECE {calc_ECE_by_segment_top1(test_calibrated, y_probs, 9)}")
        print(f"Target 2 ECE {calc_ECE_by_segment_top1(test_calibrated, y_probs, 18)}")
        print(f"Target 3 ECE {calc_ECE_by_segment_top1(test_calibrated, y_probs, 27)}")
        print(f"Water ECE {calc_ECE_by_segment_top1(test_calibrated, y_probs, 36)}")
        print(f"Target 4 ECE {calc_ECE_by_segment_top1(test_calibrated, y_probs, 45)}")
        print(f"Target 5 ECE {calc_ECE_by_segment_top1(test_calibrated, y_probs, 50)}")
        print(f"Segment 22 ECE {calc_ECE_by_segment_top1(test_calibrated, y_probs, 22)}")
        print(f"Segment 40 ECE {calc_ECE_by_segment_top1(test_calibrated, y_probs, 40)}")
        print(f"Segment 70 ECE {calc_ECE_by_segment_top1(test_calibrated, y_probs, 70)}")
        print(f"Segment 90 ECE {calc_ECE_by_segment_top1(test_calibrated, y_probs, 90)}")


        plot_multirank_reliability(
            test_probs, 
            y_probs, 
            title="Uncalibrated DANN Output - MC", 
            filename=f"{calibratorstr}_reliability_uncalibrated.png",
            n_ranks=3
        )

        #Plot calibrated reliability
        plot_multirank_reliability(
            test_calibrated, 
            y_probs, 
            title="SMS Calibrated Output - MC", 
            filename=f"{calibratorstr}_reliability_calibrated.png",
            n_ranks=3
        )
        #Plot uncalibrated vs calibrated plane distributions
        avg_ax = plot_weighted_avg_plane_distribution(test_probs, color="red", label="Uncalibrated", title=' ')
        n1_ax = plot_nth_plane_distribution(test_probs, color="red", label="Uncalibrated", title=f'Distribution of the most probable Plane Before Calibration')
        n2_ax = plot_nth_plane_distribution(test_probs, color="red", label="Uncalibrated", n=2, title=f'Distribution of the second most probable Plane Before Calibration')
        n3_ax = plot_nth_plane_distribution(test_probs, color="red", label="Uncalibrated", n=3, title=f'Distribution of the third most probable Plane Before Calibration')

        #Plot calibrated plane distributions
        plot_weighted_avg_plane_distribution(test_probs, color="blue", ax=avg_ax, label="Uncalibrated", save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_avg_plane.png"), title=f'Distribution of the weighted average position')
        plot_nth_plane_distribution(test_calibrated, color="blue", label="Calibrated", ax=n1_ax, save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_top_plane.png"), title=f'Distribution of the most probable Plane')
        plot_nth_plane_distribution(test_calibrated, color="blue", label="Calibrated", ax=n2_ax, n=2, save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_second_plane.png"), title=f'Distribution of the second most probable Plane')
        plot_nth_plane_distribution(test_calibrated, color="blue", label="Calibrated", ax=n3_ax, n=3, save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_third_plane.png"), title=f'Distribution of the third most probable Plane')


        zdiff_ax = plot_diff_distribution(test_probs, test_calibrated, color="red", label="Unweighted Calibrated Z[0] - Uncalibrated Z[0]", title=f'Difference between the most probable z position calib - uncalib')
        plot_diff_distribution(test_probs, test_calibrated, avg = True, ax = zdiff_ax, color="blue", label="Weighted Calibrated Z[0] - Uncalibrated Z[0]", save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_calibVsUncalib_zdiff.png"), title=f'Difference between the most probable z position calib - uncalib')


        #Plot calibrated and uncalibrated probability distributions for most probable plane (not POT normalised)
        cdf_ax = plot_prob_cdf(test_probs, n=1, num_bins=100, color='#d2c271', label="MC Uncalibrated")
        cdf_ax = plot_prob_cdf(test_calibrated, n=1, num_bins=100, color='#CA5454', ax = cdf_ax, save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_cdf_n1.png"), label="MC Calibrated")

        if datafile is None:
            ax = plot_prob_distributions(test_probs, n=1, num_bins=100, color='#d2c271', label="MC Uncalibrated")
            ax = plot_prob_distributions(test_calibrated, n=1, num_bins=100, ax=ax, color='#CA5454', save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_prob_distribution.png"), label="MC Calibrated")
        else:
            X_probs_data, y_none = load_data_buffered(datafile, samplerate, fill_mode = fillmode)
            test_calibrated_data = calibrator.predict_proba(X_probs_data)
            plot_prob_distributions_data_vs_mc(test_probs, test_calibrated, X_probs_data, test_calibrated_data, n=1, num_bins=100, save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_prob_distribution_n1.png"))
            plot_prob_distributions_data_vs_mc(test_probs, test_calibrated, X_probs_data, test_calibrated_data, n=2, num_bins=100, save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_prob_distribution_n2.png"))
            plot_prob_distributions_data_vs_mc(test_probs, test_calibrated, X_probs_data, test_calibrated_data, n=3, num_bins=100, save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_prob_distribution_n3.png"))
            plot_cdf_data_vs_mc(test_probs, test_calibrated, X_probs_data, test_calibrated_data, n=1, num_bins=100, save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_cdf_n1.png"))
            plot_cdf_data_vs_mc(test_probs, test_calibrated, X_probs_data, test_calibrated_data, n=2, num_bins=100, save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_cdf_n2.png"))
            plot_cdf_data_vs_mc(test_probs, test_calibrated, X_probs_data, test_calibrated_data, n=3, num_bins=100, save_filename=os.path.join(OUTPUT_DIR, f"{calibratorstr}_cdf_n3.png"))
    print(f"\nProcessing Complete. Results saved to: {os.path.abspath(OUTPUT_DIR)}")
    
    if model is None:
        SAVE_PATH = os.path.join(OUTPUT_DIR, f"{calibratorstr}_calibrator.pkl")
        print(f"Saving calibrated model to {SAVE_PATH}...")

        # This saves the entire calibrator object, including the learned matrix W
        joblib.dump(calibrator, SAVE_PATH)
        if calibratorstr == "sms":
            matrix_path  = os.path.join(OUTPUT_DIR, f"{calibratorstr}_calibrator_matrix.png")
            plot_sms_calibration_matrix(calibrator, outpath = matrix_path)


        print("Model saved successfully.")