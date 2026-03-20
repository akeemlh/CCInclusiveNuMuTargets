import ROOT
import os, argparse


# Set ROOT to batch mode to avoid popping up windows during processing
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)  # Remove stat box

# ----------------------------------------------------------------------
# HELPER FUNCTIONS
# ----------------------------------------------------------------------

def get_chain_from_list(list_file, tree_name):
    """
    Reads a text file of root file paths and adds them to a TChain.
    This effectively 'iterates' through all files in the list.
    """
    chain = ROOT.TChain(tree_name)
    if not os.path.exists(list_file):
        print(f"Error: File {list_file} not found.")
        return None
        
    with open(list_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            # Strip whitespace/newlines and add to chain
            file_path = line.strip()
            if file_path:
                chain.Add(file_path)
    
    print(f"Loaded {chain.GetEntries()} entries from {list_file}")
    return chain

def create_and_save_plot(data_chain, mc_chain, cut_string, plot_name, title_suffix, output_dir):
    """
    Projects data and MC into histograms, normalizes them, and saves the overlay.
    """
    # Define the variable to plot (index 0 of the array)
    var_name = "ANN_plane_probs[0]"
    
    # Create unique names for ROOT internal registry
    h_name_data = f"h_data_{plot_name}"
    h_name_mc = f"h_mc_{plot_name}"
    
    # Define Histograms: 100 bins from 0 to 1
    h_data = ROOT.TH1D(h_name_data, f"ANN_plane_probs[0] ({title_suffix})", 100, 0, 1)
    h_mc = ROOT.TH1D(h_name_mc, f"ANN_plane_probs[0] ({title_suffix})", 100, 0, 1)
    
    # --- Projecting (Iterating) ---
    # TChain.Project loops over all entries in all files passing the 'cut_string'
    data_chain.Project(h_name_data, var_name, cut_string)
    mc_chain.Project(h_name_mc, var_name, cut_string)
    dataint = h_data.Integral()
    mcint = h_mc.Integral()
    # --- Normalization ---
    # Scale so the sum of all bins equals 1 (Probability Density)
    if h_data.Integral() > 0:
        h_data.Scale(1.0 / h_data.Integral())
    if h_mc.Integral() > 0:
        h_mc.Scale(1.0 / h_mc.Integral())

    # --- Styling ---
    # Data: Black markers
    h_data.SetLineColor(ROOT.kBlack)
    h_data.SetMarkerColor(ROOT.kBlack)
    h_data.SetMarkerStyle(20) # Full circle
    h_data.SetMarkerSize(0.8)
    # MC: Red line/fill
    h_mc.SetLineColor(ROOT.kRed)
    h_mc.SetLineWidth(2)
    h_mc.SetFillColorAlpha(ROOT.kRed, 0.3) # Transparent red fill
    
    # Calculate Max Y for plotting range (add 20% headroom for legend)
    max_y = max(h_data.GetMaximum(), h_mc.GetMaximum()) * 1.2
    h_data.SetMaximum(max_y)
    
    # --- Drawing ---
    canvas = ROOT.TCanvas(f"c_{plot_name}", plot_name, 800, 600)
    
    # Draw Data first (to set axes), then MC
    # "E1" = Error bars, "HIST" = Histogram line/fill
    h_data.Draw("E1") 
    h_mc.Draw("HIST SAME")
    h_data.Draw("E1 SAME") # Redraw data on top
    
    # Axis Labels
    h_data.GetXaxis().SetTitle("ANN_plane_probs[0]")
    h_data.GetYaxis().SetTitle("Probability Density")
    
    # Legend
    legend = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.AddEntry(h_data, "Data", "lp")
    legend.AddEntry(h_mc, "Simulation (MC)", "lf")
    legend.Draw()
    
    # Save
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    canvas.SaveAs(f"{output_dir}/{plot_name}.png")

    kscore = h_data.KolmogorovTest(h_mc, "M")
    print(f'{title_suffix} Kolmogorov Score: {kscore}')
    return kscore
    # canvas.SaveAs(f"{output_dir}/{plot_name}.pdf") # Uncomment if you want PDFs

# ----------------------------------------------------------------------
# MAIN EXECUTION
# ----------------------------------------------------------------------

if __name__ == "__main__":
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(True)

    parser = argparse.ArgumentParser(description="Generate reliability diagram with ECE calculation.")
    parser.add_argument("-d", "--datafile", required=True, type=str, help="Data file")
    parser.add_argument("-m", "--mc", required=True, type=str, help="MC file")
    parser.add_argument("-o", "--output", default="./", type=str, help="Output directory")

    args = parser.parse_args()

    datafile = args.datafile
    mcfile = args.mc
    output_dir=args.output

    treename = "MasterAnaDev"


        # 1. Load Chains
    chain_data = get_chain_from_list(datafile, treename)
    chain_mc = get_chain_from_list(mcfile, treename)
    
    if not chain_data or not chain_mc:
        print("Error loading chains. Exiting.")
        exit(0)

    # slice_n_hits
    # ------------------------------------------------
    print("Processing slice_n_hits plots...")
    
    # 0 to 20000 in 10 bins -> Step size 2000
    max_n = 1000
    n_bins_n = 20
    step_n = max_n / n_bins_n
    results = []
    
    for i in range(n_bins_n):
        low_edge = i * step_n
        high_edge = (i + 1) * step_n
        
        # Cut: low <= E < high
        cut = f"slice_n_hits >= {low_edge} && slice_n_hits < {high_edge}"
        name = f"slice_n_hits_bin_{i+1}"
        title = f"slice_n_hits [{int(low_edge)}-{int(high_edge)}]"
        
        results.append((title, create_and_save_plot(chain_data, chain_mc, cut, name, title, output_dir)))
    # ------------------------------------------------
    # 2. Multiplicity Plots
    # ------------------------------------------------
    print("Processing Multiplicity plots...")
    
    # multiplicity 1 to 5
    for m in range(1, 6):
        cut = f"multiplicity == {m}"
        name = f"mult_{m}"
        title = f"Mult={m}"
        results.append((title, create_and_save_plot(chain_data, chain_mc, cut, name, title, output_dir)))

        
    # multiplicity >= 6
    cut = "multiplicity >= 6"
    name = "mult_ge6"
    title = "Mult>=6"
    #create_and_save_plot(chain_data, chain_mc, cut, name, "Mult>=6", output_dir)
    results.append((title, create_and_save_plot(chain_data, chain_mc, cut, name, title, output_dir)))


    # ------------------------------------------------
    # 3. Recoil Energy Plots (MasterAnaDev_ANN_recoil_E)
    # ------------------------------------------------
    print("Processing Recoil Energy plots...")
    
    # 0 to 20000 in 10 bins -> Step size 2000
    max_e = 20000
    n_bins_e = 10
    step_e = max_e / n_bins_e
    
    for i in range(n_bins_e):
        low_edge = i * step_e
        high_edge = (i + 1) * step_e
        
        # Cut: low <= E < high
        cut = f"MasterAnaDev_ANN_recoil_E >= {low_edge} && MasterAnaDev_ANN_recoil_E < {high_edge}"
        name = f"recoil_bin_{i+1}"
        title = f"Recoil E [{int(low_edge)}-{int(high_edge)}]"
        
        results.append((title, create_and_save_plot(chain_data, chain_mc, cut, name, title, output_dir)))

    # ------------------------------------------------
    # 4. Interaction Type Plots (MasterAnaDev_intType)
    # ------------------------------------------------
    print("Processing IntType plots...")
    
    # intType 1 to 10
    for t in range(1, 11):
        cut = f"MasterAnaDev_intType == {t}"
        name = f"intType_{t}"
        create_and_save_plot(chain_data, chain_mc, cut, name, f"IntType={t}", output_dir)

    # ------------------------------------------------
    # MasterAnaDev_E
    # ------------------------------------------------
    print("Processing MasterAnaDev_E plots...")
    
    # 0 to 20000 in 10 bins -> Step size 2000
    max_e = 20000
    n_bins_e = 10
    step_e = max_e / n_bins_e
    
    for i in range(n_bins_e):
        low_edge = i * step_e
        high_edge = (i + 1) * step_e
        
        # Cut: low <= E < high
        cut = f"MasterAnaDev_E >= {low_edge} && MasterAnaDev_E < {high_edge}"
        name = f"E_bin_{i+1}"
        title = f"E [{int(low_edge)}-{int(high_edge)}]"
        
        results.append((title, create_and_save_plot(chain_data, chain_mc, cut, name, title, output_dir)))

    # Bjorken X
    # ------------------------------------------------
    print("Processing MasterAnaDev_x plots...")
    
    # 0 to 20000 in 10 bins -> Step size 2000
    max_x = 1.2
    n_bins_x = 6
    step_x = max_x / n_bins_x
    
    for i in range(n_bins_x):
        low_edge = i * step_x
        high_edge = (i + 1) * step_x
        
        # Cut: low <= E < high
        cut = f"MasterAnaDev_x >= {low_edge} && MasterAnaDev_x < {high_edge}"
        name = f"X_bin_{i+1}"
        title = f"X [{int(low_edge)}-{int(high_edge)}]"
        
        results.append((title, create_and_save_plot(chain_data, chain_mc, cut, name, title, output_dir)))

    # MasterAnaDev_W
    # ------------------------------------------------
    print("Processing MasterAnaDev_W plots...")
    
    # 0 to 20000 in 10 bins -> Step size 2000
    max_w = 20000
    n_bins_w = 10
    step_w = max_w / n_bins_w
    
    for i in range(n_bins_w):
        low_edge = i * step_w
        high_edge = (i + 1) * step_w
        
        # Cut: low <= E < high
        cut = f"MasterAnaDev_W >= {low_edge} && MasterAnaDev_W < {high_edge}"
        name = f"W_bin_{i+1}"
        title = f"W [{int(low_edge)}-{int(high_edge)}]"
        
        results.append((title, create_and_save_plot(chain_data, chain_mc, cut, name, title, output_dir)))

    print(f"Done. Plots saved to directory: {output_dir}")
    results.sort(key=lambda x: x[1])
    
    print(f"\n{'Variable':<30} | {'KS Score':<10}")
    print("-" * 45)
    for name, score in results:
        print(f"{name:<60} | {score:.6f}")