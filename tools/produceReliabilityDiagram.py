import ROOT
import argparse
import sys
import os

def calculate_truth_segment(mdl, plane, tgt_id):
    """
    Converts C++ GetTruthSegment logic to Python.
    """

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

# ==========================================
# MAIN SCRIPT
# ==========================================

def create_reliability_diagram(file_list_path, tree_name, output_name):
    # 1. Setup Data using TChain
    chain = ROOT.TChain(tree_name)
    
    if not os.path.exists(file_list_path):
        print(f"Error: File list '{file_list_path}' not found.")
        return

    print(f"Reading file list: {file_list_path}")
    files_added = 0
    with open(file_list_path, 'r') as f:
        for line in f:
            root_file = line.strip()
            if root_file:
                chain.Add(root_file)
                files_added += 1

    if chain.GetEntries() == 0:
        print(f"Error: Chain is empty or no files loaded.")
        return

    print(f"Chain loaded with {files_added} files and {chain.GetEntries()} total events.")

    # 2. Optimization: Disable all branches, enable only what we need
    chain.SetBranchStatus("*", 0) 
    
    # Enable Truth Branches
    chain.SetBranchStatus("truth_vtx_module", 1)
    chain.SetBranchStatus("truth_vtx_plane", 1)
    chain.SetBranchStatus("truth_target_code", 1)
    
    # Enable Prediction Branches
    chain.SetBranchStatus("ANN_segments", 1)
    chain.SetBranchStatus("ANN_plane_probs", 1) 

    # 3. Create Histograms
    n_bins = 20
    
    # Standard histograms for counts
    h_total = ROOT.TH1F("h_total", "", n_bins, 0.0, 1.0)
    h_correct = ROOT.TH1F("h_correct", "", n_bins, 0.0, 1.0)
    
    # TProfile to calculate the actual average confidence per bin (Crucial for ECE)
    h_avg_conf = ROOT.TProfile("h_avg_conf", "", n_bins, 0.0, 1.0)
    
    h_total.SetStats(0)
    h_correct.SetStats(0)

    # 4. Event Loop
    print("Processing events...")
    
    count = 0
    total_events = chain.GetEntries()
    
    for event in chain:
        count += 1
        if count % 10000 == 0:
            sys.stdout.write(f"\rProcessed {count}/{total_events} events")
            sys.stdout.flush()

        # A. Get Truth Segment
        true_seg = get_truth_segment(event.truth_vtx_module, 
                                       event.truth_vtx_plane, 
                                       event.truth_target_code)
        
        # B. Get Prediction Info
        pred_seg = event.ANN_segments[0]
        confidence = event.ANN_plane_probs[0]

        # C. Fill Histograms
        h_total.Fill(confidence)
        
        # Track the average confidence for this bin
        h_avg_conf.Fill(confidence, confidence)
        
        if pred_seg == true_seg:
            h_correct.Fill(confidence)
        #if count > 20000:
        #    break

    print("\nEvent loop complete.")

    # 5. Calculate Accuracy and ECE
    h_acc = h_correct.Clone("h_acc")
    h_acc.Divide(h_total)

    # --- ECE CALCULATION START ---
    ece = 0.0
    total_samples = h_total.GetEntries()
    
    if total_samples > 0:
        for i in range(1, n_bins + 1):
            n_i = h_total.GetBinContent(i)
            if n_i == 0: continue
            
            # Accuracy in this bin
            acc_i = h_acc.GetBinContent(i)
            
            # Average Confidence in this bin (from TProfile)
            conf_i = h_avg_conf.GetBinContent(i)
            
            # Weighted absolute difference
            ece += (n_i / total_samples) * abs(acc_i - conf_i)
            
    print(f"\n========================================")
    print(f"Calculated ECE: {ece:.4f}")
    print(f"========================================")
    # --- ECE CALCULATION END ---

    # 6. Styling
    c1 = ROOT.TCanvas("c1", "Reliability Diagram", 800, 800)
    c1.SetGrid()
    
    h_acc.SetTitle("Reliability Diagram") 
    h_acc.GetXaxis().SetTitle("Confidence")
    h_acc.GetYaxis().SetTitle("Accuracy")
    h_acc.GetYaxis().SetRangeUser(0.0, 1.0)
    
    h_acc.SetFillColor(ROOT.kAzure-4)
    h_acc.SetLineColor(ROOT.kBlack)
    h_acc.SetLineWidth(1)
    
    h_acc.Draw("HIST BAR")

    # 7. Create Calibration Gap Boxes
    boxes = [] 
    for i in range(1, n_bins + 1):
        count_in_bin = h_total.GetBinContent(i)
        if count_in_bin <= 0:
            continue
            
        accuracy = h_acc.GetBinContent(i)
        
        # Note: For the Visual Diagram, it is cleaner to draw the box to the 
        # Average Confidence (h_avg_conf) rather than the Bin Center.
        confidence_target = h_avg_conf.GetBinContent(i)
        
        bin_width = h_acc.GetBinWidth(i)
        x_min = h_acc.GetBinLowEdge(i)
        x_max = x_min + bin_width
        
        box = ROOT.TBox(x_min, accuracy, x_max, confidence_target)
        box.SetFillColor(ROOT.kRed-7)
        box.SetFillStyle(3354) 
        box.SetLineColor(ROOT.kRed-7)
        box.SetLineWidth(1)
        
        box.Draw("same")
        boxes.append(box)

    h_acc.Draw("HIST BAR same")

    # 8. Diagonal Line
    line = ROOT.TLine(0.0, 0.0, 1.0, 1.0)
    line.SetLineStyle(2)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineWidth(2)
    line.Draw("same")

    # 9. Legend and Labels
    legend = ROOT.TLegend(0.2, 0.75, 0.5, 0.9)
    legend.SetBorderSize(1)
    legend.SetFillStyle(1001)
    legend.SetFillColor(ROOT.kWhite)
    
    legend.AddEntry(h_acc, "Outputs", "f")
    
    gap_entry = ROOT.TH1F("gap_legend", "", 1, 0, 1)
    gap_entry.SetFillColor(ROOT.kRed-7)
    gap_entry.SetFillStyle(3354)
    gap_entry.SetLineColor(ROOT.kRed-7)
    legend.AddEntry(gap_entry, "Calibration gap", "f")
    
    legend.Draw()

    # Draw ECE on the plot
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.035)
    latex.DrawLatex(0.6, 0.2, f"ECE = {ece:.3f}")

    c1.SaveAs(output_name)
    print(f"Reliability diagram saved to {output_name}")

if __name__ == "__main__":
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(True)

    parser = argparse.ArgumentParser(description="Generate reliability diagram with ECE calculation.")
    parser.add_argument("-i", "--input", required=True, type=str, help="Path to text file containing list of ROOT files")
    parser.add_argument("-t", "--tree", default="mytree", type=str, help="Name of the TTree")
    parser.add_argument("-o", "--output", default="reliability_diagram.png", type=str, help="Output image filename")

    args = parser.parse_args()

    create_reliability_diagram(args.input, args.tree, args.output)