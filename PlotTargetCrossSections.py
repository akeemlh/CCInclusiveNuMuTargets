import argparse
import sys
import ROOT
from array import array

ROOT.gROOT.SetBatch(1)

def main():
    parser = argparse.ArgumentParser(description="Plot Target Cross Sections and Event Rates from MINERvA ROOT files.")
    parser.add_argument("--indir", type=str, required=True, help="Path to the input director containing histograms.")
    parser.add_argument("--prefix", type=str, nargs='*', default=["pTmu", "pZmu", "BjorkenX", "Emu", "Enu"], help="Kinematic variables (e.g., pTmu, pZmu, BjorkenX).")
    parser.add_argument("--target", type=str, nargs='*', default=["Iron", "Lead", "Carbon", "Water"], help="Target strings (e.g., Iron, Lead, Water).")
    parser.add_argument("--nutype", type=str, nargs='*', default=["Neutrino"], help="Neutrino or Antineutrino.")
    parser.add_argument("--pot_data", type=float, default=-1, help="Data POT for scaling. (If left as default -1 will try to read POT from file)")
    parser.add_argument("--pot_mc", type=float, default=-1, help="MC POT for scaling. (If left as default -1 will try to read POT from file)")
    args = parser.parse_args()

    # Prevent PyROOT from destroying histograms when the file closes
    ROOT.TH1.AddDirectory(False)
    
    # Load PlotUtils (Requires MINERvA MAT environment to be sourced)
    try:
        ROOT.gSystem.Load("libPlotUtils")
    except Exception as e:
        print(f"Warning: Could not load libPlotUtils. Ensure your MAT environment is sourced. Error: {e}")

    # Instantiate the MINERvA Plotter
    plotter = ROOT.PlotUtils.MnvPlotter()
    plotter.ApplyStyle(ROOT.PlotUtils.kDefaultStyle)
    
    # Configure common axis logic
    plotter.axis_label_font = 43
    plotter.axis_label_size = 25
    plotter.headroom = 1.2
    plotter.mc_line_width=0
    plotter.ApplyAxisStyle()

    # ---------------------------------------------------------
    # Plotting Cross Sections
    # ---------------------------------------------------------
    for prefix in args.prefix:
        for nutype in args.nutype:
            for target in args.target:
                filename = args.indir  + "/"+nutype+"_"+target+"_"+prefix+"_crossSection.root"
                f = ROOT.TFile.Open(filename, "READ")
                if not f or f.IsZombie():
                    print(f"Error: Could not open {filename}")
                    sys.exit(1)
                print("Plotting Cross Sections...")
                plot_cross_sections(f, plotter, prefix, nutype, target)
                print("Plotting Event Rates...")
                #plot_event_rates(f, plotter, prefix, nutype, target)
                print("Plotting Migration Matrix...")
                plot_migration(f, plotter, prefix, nutype, target)
                print("Plotting Sidebands...")
                #plot_sidebands(f, plotter, prefix, nutype, target)

    # ---------------------------------------------------------
    # Plotting Event Rates (Optional based on your C++ structure)
    # ---------------------------------------------------------
    


    f.Close()

def plot_cross_sections(f, plotter, prefix, nutype, target):
    # 1. Fetch histograms (matching the names saved in your C++ script)
    tempxsec = f.Get("crossSectionTuned")
    tempsimxsec = f.Get("simulatedCrossSection")
    temp2p2h = f.Get("simulatedCrossSection2P2H")
    tempdis = f.Get("simulatedCrossSectionDIS")
    tempres = f.Get("simulatedCrossSectionRES")
    tempqe = f.Get("simulatedCrossSectionQE")
    tempother = f.Get("simulatedCrossSectionOther")
    mcPOT = f.Get("mcPOT")
    dataPOT = f.Get("dataPOT")

    if not all([tempxsec, tempsimxsec, temp2p2h]):
        print("Error: Missing required cross section histograms in the ROOT file.")
        return

    # Clone to avoid modifying the original file objects in memory
    tempxsec = tempxsec.Clone()
    tempsimxsec = tempsimxsec.Clone()
    
    # 2. Scale by 1e39
    scale_factor = 1e39

    # 3. Build TObjArray for stacked MC
    arr = ROOT.TObjArray()
    arr.Add(temp2p2h)
    arr.Add(tempdis)
    arr.Add(tempres)
    arr.Add(tempqe)
    arr.Add(tempother)

    # 4. Set colors (matching your C++ HEX codes)
    channelColArr = array('i', [
        ROOT.TColor.GetColor("#d2c271"), 
        ROOT.TColor.GetColor("#CA5454"), 
        ROOT.TColor.GetColor("#40B0A6"), 
        ROOT.TColor.GetColor("#9285d5"), 
        ROOT.TColor.GetColor("#c2c2c2")
    ])

    # 5. Determine titles and axis ranges
    xaxislabel = ""
    y_title = "Events#times10^{3}"
    xmax = -1111

    if prefix == "pTmu":
        tmp = "Muon p_{T}"
        xaxislabel = "p_{T, #mu} [GeV/c]"
        y_title = "d#sigma/dp_{T} (x10^{-39}) (cm^{2}/(GeV/c)/Nucleon)"
        xmax = 13
    elif prefix == "pZmu":
        tmp = "Muon p_{Z}"
        xaxislabel = "p_{Z, #mu} [GeV/c]"
        y_title = "d#sigma/dp_{Z} (x10^{-39}) (cm^{2}/(GeV/c)/Nucleon)"
        xmax = 13
    elif "BjorkenX" in prefix: #Since we're testing other BjorkenX binnings
        tmp = "Bjorken X"
        xaxislabel = "Bjorken X"
        y_title = "d#sigma/dX (x10^{-39}) (cm^{2}/Nucleon)"
    elif prefix=="Erecoil":
        tmp = "E_{recoil}"
        y_title="d#sigma/dE_{recoil} (x10^{-39}) (cm^{2}/GeV/Nucleon)"
        xaxislabel = "E_{recoil} [GeV]"
    elif prefix=="Emu":
        tmp = "E_{#mu}"
        y_title="d#sigma/dE_{#mu} (x10^{-39}) (cm^{2}/GeV/Nucleon)"
        xaxislabel = "E_{#mu} [GeV]"
        xmax = 11
    elif prefix=="Enu":
        tmp = "E_{#nu}"
        y_title="d#sigma/dE_{#nu} (x10^{-39}) (cm^{2}/GeV/Nucleon)"
        xaxislabel = "E_{#nu} [GeV]"
    else:
        tmp = prefix

    #if xmax != -1111:
    #        tempsimxsec.GetXaxis().SetRange(0, xmax)
    #        tempxsec.GetXaxis().SetRange(0, xmax)

    # 6. Setup Canvas and Pads (Top and Bottom)
    can = ROOT.TCanvas("can", "Cross Section", 800, 800)
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.35, 1, 1.0)
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.35)
    
    pad1.SetBottomMargin(0)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.25)
    
    pad1.Draw()
    pad2.Draw()

    plotter.draw_normalized_to_bin_width = False
    if target == "Water":
        dataPOTEmpty = f.Get("dataPOTWaterEmpty").GetVal()
        dataPOTFilled = f.Get("dataPOTWaterFilled").GetVal()
        potstr_filled = f"Water Filled POT {dataPOTFilled:.2e}"
        potstr_empty = f"Water Empty POT: {dataPOTEmpty:.2e}"

        # Standard console output
        print(f"waterFilledPOTData: {dataPOTFilled}")
        print(f"waterEmptyPOTData: {dataPOTEmpty}")

        # Adding labels to the plotter
        plotter.AddPlotLabel(potstr_filled, 0.18, 0.93, 0.025, 1, 42, 13, 0)
        plotter.AddPlotLabel(potstr_empty, 0.18, 0.90, 0.025, 1, 42, 13, 0)
    else:
        dataPOT = f.Get("dataPOT").GetVal()
        potstr = f"Data POT: {dataPOT:.2e}"
        plotter.AddPlotLabel(potstr , 0.18, 0.93, 0.025, 1, 42, 13, 0)

    if nutype == "Neutrino":
        horncurrent = "FHC"
    elif nutype == "Antieutrino":
        horncurrent = "RHC"
    title = f"{horncurrent} {tmp} Cross-Section on {target}"

    if prefix == "BjorkenX":
        pad1.SetLogx()
        pad2.SetLogx()

    hists_to_format = [tempxsec, tempsimxsec, temp2p2h, tempdis, tempres, tempqe, tempother]
    for hist in hists_to_format:
        hist.Scale(scale_factor)
        hist.GetXaxis().SetRange(1, xmax)

    # 7. Draw Stacked MC
    pad1.cd()
    plotter.axis_minimum = 0.001
    plotter.DrawDataStackedMC(
        tempxsec, 
        arr, 
        channelColArr, 
        1, 
        "TR", 
        "Data", 
        1001, 
        "", 
        "", 
        False, 
        1, 
        xmax
    )



    # 8. Draw Ratio
    pad2.cd()
    plotter.axis_minimum = 0.0 # reset
    #if xmax != -1111:
    #    tempsimxsec.GetXaxis().SetRangeUser(0, xmax)
    #    tempxsec.GetXaxis().SetRangeUser(0, xmax)


    plotter.DrawDataMCRatio(tempxsec, tempsimxsec, 1, True, True, 0.5, 1.5, "")
    

    # Draw axis back on top in pad1 for clarity
    pad1.cd()
    tempxsec.Draw("AXIS SAME")

    can.cd()

    plotter.AddPlotLabel(title, 0.5, 0.97, 25, 1, 43, 22, 0)
    plotter.AddPlotLabel(y_title, 0.07, 0.9, 0.03, 1, 42, 33, 90)
    plotter.AddPlotLabel("Data/MC", 0.07, 0.3, 0.03, 1, 42, 33, 90)
    plotter.AddPlotLabel(xaxislabel, 0.9, 0.08, 0.03, 1, 42, 33, 0)
    can.Update()
    # Save output
    outname = f"{target}_{prefix}_{nutype}_CrossSection"
    plotter.MultiPrint(can, outname, "png,pdf,C")


def plot_event_rates(f, plotter, prefix, nutype, target):
    # 1. Fetch histograms (matching the names saved in your C++ script)
    tempxsec = f.Get("crossSectionTuned")
    tempsimxsec = f.Get("simulatedCrossSection")
    temp2p2h = f.Get("simulatedCrossSection2P2H")
    tempdis = f.Get("simulatedCrossSectionDIS")
    tempres = f.Get("simulatedCrossSectionRES")
    tempqe = f.Get("simulatedCrossSectionQE")
    tempother = f.Get("simulatedCrossSectionOther")
    mcPOT = f.Get("mcPOT")
    dataPOT = f.Get("dataPOT")

    if not all([tempxsec, tempsimxsec, temp2p2h]):
        print("Error: Missing required cross section histograms in the ROOT file.")
        return

    # Clone to avoid modifying the original file objects in memory
    tempxsec = tempxsec.Clone()
    tempsimxsec = tempsimxsec.Clone()
    
    # 2. Scale by 1e39
    scale_factor = 1e39
    for hist in [tempxsec, tempsimxsec, temp2p2h, tempdis, tempres, tempqe, tempother]:
        hist.Scale(scale_factor)

    # 3. Build TObjArray for stacked MC
    arr = ROOT.TObjArray()
    arr.Add(temp2p2h)
    arr.Add(tempdis)
    arr.Add(tempres)
    arr.Add(tempqe)
    arr.Add(tempother)

    # 4. Set colors (matching your C++ HEX codes)
    channelColArr = array('i', [
        ROOT.TColor.GetColor("#d2c271"), 
        ROOT.TColor.GetColor("#CA5454"), 
        ROOT.TColor.GetColor("#40B0A6"), 
        ROOT.TColor.GetColor("#9285d5"), 
        ROOT.TColor.GetColor("#c2c2c2")
    ])

    # 5. Determine titles and axis ranges
    xaxislabel = ""
    y_title = "Events#times10^{3}"
    xmax = -1111

    if prefix == "pTmu":
        tmp = "Muon p_{T}"
        xaxislabel = "p_{T, #mu} [GeV/c]"
        y_title = "d#sigma/dp_{t} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)"
        xmax = 13
    elif prefix == "BjorkenX":
        tmp = "Bjorken X"
        xaxislabel = "Bjorken X"
        y_title = "d#sigma/dX (x10^{-39}) (cm^{2}/Nucleon)"
    else:
        tmp = prefix

    if xmax != -1111:
        tempsimxsec.GetXaxis().SetRangeUser(0, xmax)
        tempxsec.GetXaxis().SetRangeUser(0, xmax)

    if target == "Water":
        dataPOTEmpty = f.Get("dataPOTWaterEmpty")
        dataPOTFilled = f.Get("dataPOTWaterFilled")
        potstr_filled = f"Water Filled POT {dataPOTFilled:.2e}"
        potstr_empty = f"Water Empty POT: {dataPOTEmpty:.2e}"

        # Standard console output
        print(f"waterFilledPOTData: {dataPOTFilled}")
        print(f"waterEmptyPOTData: {dataPOTEmpty}")

        # Adding labels to the plotter
        plotter.AddPlotLabel(potstr_filled, 0.18, 0.93, 0.025, 1, 42, 13, 0)
        plotter.AddPlotLabel(potstr_empty, 0.18, 0.90, 0.025, 1, 42, 13, 0)
    else:
        potstr = f"Data POT: {dataPOT:.2e}"
        plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.93, 0.025, 1, 42, 13, 0)

    title = f"{nutype} {tmp} {target}"

    # 6. Setup Canvas and Pads (Top and Bottom)
    can = ROOT.TCanvas("can", "Cross Section", 800, 800)
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.35, 1, 1.0)
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.35)
    
    pad1.SetBottomMargin(0)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.25)
    
    pad1.Draw()
    pad2.Draw()

    if prefix == "BjorkenX":
        pad1.SetLogx()
        pad2.SetLogx()

    # 7. Draw Stacked MC
    pad1.cd()
    plotter.axis_minimum = 0.001
    plotter.DrawDataStackedMC(
        tempxsec, 
        arr, 
        channelColArr, 
        1, 
        "TR", 
        "Data", 
        1001, 
        "", 
        "", 
        False, 
        1, 
        xmax if xmax != -1111 else -1
    )

    plotter.AddPlotLabel(title, 0.5, 0.97, 25, 1, 43, 22, 0)
    plotter.AddPlotLabel(y_title, 0.08, 0.9, 0.03, 1, 42, 33, 90)

    # 8. Draw Ratio
    pad2.cd()
    plotter.axis_minimum = 0.0 # reset
    plotter.DrawDataMCRatio(tempxsec, tempsimxsec, 1, True, True, 0.5, 1.5, "")
    
    plotter.AddPlotLabel("Data/MC", 0.08, 0.3, 0.03, 1, 42, 33, 90)
    plotter.AddPlotLabel(xaxislabel, 0.9, 0.08, 0.03, 1, 42, 33, 0)

    # Draw axis back on top in pad1 for clarity
    pad1.cd()
    tempxsec.Draw("AXIS SAME")

    can.cd()
    can.Update()
    
    # Save output
    outname = f"{target}_{prefix}_{nutype}_CrossSection"
    plotter.MultiPrint(can, outname, "png,pdf,C")


def plot_event_rates(f, plotter, args):
    """
    Similar logic to plot_cross_sections, but pulling the Event Rate histograms 
    (e.g., MCData2p2h, BackgroundUSPlastic, etc.), applying the dataPOT/mcPOT scaling, 
    and building your background stacks.
    """
    pass



def plot_migration(f, plotter, prefix, nutype, target):
    # 1. Fetch histograms (matching the names saved in your C++ script)
    migration = f.Get("migration")
    mcPOT = f.Get("mcPOT")
    dataPOT = f.Get("dataPOT")

    if not migration:
        print("Error: Missing required migration histograms in the ROOT file.")
        return

    # Clone to avoid modifying the original file objects in memory
    tempxsec = tempxsec.Clone()
    tempsimxsec = tempsimxsec.Clone()
    
    # 2. Scale by 1e39
    scale_factor = 1e39

    # 5. Determine titles and axis ranges
    xaxislabel = ""
    y_title = "Events#times10^{3}"
    xmax = -1111

    if prefix == "pTmu":
        tmp = "Muon p_{T}"
        yaxislabel = "True p_{T, #mu} [GeV/c]"
        xaxislabel = "p_{T, #mu} [GeV/c]"
        xmax = 13
    if prefix == "pZmu":
        tmp = "Muon p_{Z}"
        yaxislabel = "True p_{Z, #mu} [GeV/c]"
        xaxislabel = "p_{Z, #mu} [GeV/c]"
        xmax = 13
    elif "BjorkenX" in prefix: #Since we're testing other BjorkenX binnings
        tmp = "Bjorken X"
        yaxislabel = "True X"
        xaxislabel = "X"
        can.SetLogx()
        can.SetLogy()
    elif prefix=="Erecoil":
        tmp = "E_{recoil}"
        yaxislabel = "Truth E_{recoil} [GeV]"
        xaxislabel = "E_{recoil} [GeV]"
    elif prefix=="Emu":
        tmp = "E_{#mu}"
        y_title="d#sigma/dE_{#mu} (x10^{-39}) (cm^{2}/GeV/Nucleon)"
        yaxislabel = "True E_{#mu} [GeV]"
        xaxislabel = "E_{#mu} [GeV]"
        xmax = 11
    elif prefix=="Enu":
        tmp = "E_{#nu}"
        yaxislabel = "True E_{#mu} [GeV]"
        xaxislabel = "E_{#mu} [GeV]"
    else:
        tmp = prefix

    if xmax != -1111:
            migration.GetXaxis().SetRangeUser(0, xmax)
            migration.GetYaxis().SetRangeUser(0, xmax)

    # 6. Setup Canvas and Pads (Top and Bottom)
    can = ROOT.TCanvas("can", "Cross Section", 800, 800)

    plotter.draw_normalized_to_bin_width = False
    if target == "Water":
        dataPOTEmpty = f.Get("dataPOTWaterEmpty").GetVal()
        dataPOTFilled = f.Get("dataPOTWaterFilled").GetVal()
        potstr_filled = f"Water Filled POT {dataPOTFilled:.2e}"
        potstr_empty = f"Water Empty POT: {dataPOTEmpty:.2e}"

        # Standard console output
        print(f"waterFilledPOTData: {dataPOTFilled}")
        print(f"waterEmptyPOTData: {dataPOTEmpty}")

        # Adding labels to the plotter
        plotter.AddPlotLabel(potstr_filled, 0.18, 0.93, 0.025, 1, 42, 13, 0)
        plotter.AddPlotLabel(potstr_empty, 0.18, 0.90, 0.025, 1, 42, 13, 0)
    else:
        dataPOT = f.Get("dataPOT").GetVal()
        potstr = f"Data POT: {dataPOT:.2e}"
        plotter.AddPlotLabel(potstr , 0.18, 0.93, 0.025, 1, 42, 13, 0)

    title = f"{nutype} {tmp} {target}"


    # 7. Draw Stacked MC

    plotter.AddPlotLabel(title, 0.5, 0.97, 25, 1, 43, 22, 0)

    plotter.axis_title_offset_z = 1.25
    plotter.DrawNormalizedMigrationHistogram(migration, False, False, True, False, 0.02)
    plotter.WritePreliminary( 0.65, 0.2, 0.035, True)
    plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 25, 1, 43, 22, 0)
    plotter.AddPlotLabel(yaxislabel.c_str() , 0.08, 0.9, 0.03, 1, 62, 33, 90)
    plotter.AddPlotLabel(xaxislabel.c_str() , 0.9, 0.08, 0.03, 1, 62, 33, 0)

    can.cd()
    can.Update()
    # Save output
    outname = f"{target}_{prefix}_{nutype}_Migration"
    plotter.MultiPrint(can, outname, "png,pdf,C")
    can.SetLogx(0)
    can.SetLogy(0)


if __name__ == "__main__":
    main()