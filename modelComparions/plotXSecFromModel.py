import ROOT
import array
import sys
import numpy as np
import argparse

ROOT.ROOT.EnableImplicitMT()
ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser( prog='plotXSecFromModel',
                description='Extract and save the cross section from a given input model')
parser.add_argument('input', type=str, help="the input file")
parser.add_argument('-p', '--pdg', default=14, type=int, help="the neutrino pdg to look for (Default: 14)")
parser.add_argument('-d', '--draw', action="store_true", help="Draw plots of model binned cross sections to a file?")
args = parser.parse_args()

input_files = args.input
pdg = args.pdg
draw = args.draw
chain = ROOT.TChain("FlatTree_VARS")
for f in input_files:
    chain.Add(f)
ntrees = chain.GetNtrees()
print(f"Number of trees: {ntrees}")

# Derive prefix from first input file name (before extension)
import os
first_file = input_files
first_file_basename = os.path.basename(first_file)
first_file_prefix = os.path.splitext(first_file_basename)[0]

# Get immediate parent directory of the file
parent_folder = os.path.basename(os.path.dirname(first_file))
material=parent_folder.capitalize()
if material == "Ch":
    material = "Tracker"
model = ""
if "G18_02a_02_11a" in first_file_prefix:
    model = "GENIE_G18_02a_02_11a"
if "G18_02b_02_11a" in first_file_prefix:
    model = "GENIE_G18_02b_02_11a"
if "G18_10a_02_11a" in first_file_prefix:
    model = "GENIE_G18_10a_02_11a"
if "G18_10b_02_11a" in first_file_prefix:
    model = "GENIE_G18_10b_02_11a"
if "SF" in first_file_prefix:
    model = "NEUT_SF_ma103"
if "LFG" in first_file_prefix:
    model = "NEUT_LFG_ma105"
out_prefix = f"{material}_{model}"
print(f"out_prefix {out_prefix}, file: {os.path.dirname(first_file)}, parent_folder {parent_folder}, first_file_prefix {first_file_prefix}")
# Add input files
df = ROOT.RDataFrame("FlatTree_VARS", input_files)

# Bins
ptbins = [0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1, 1.25, 1.5, 2.5, 4.5]
pzbins = [1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20, 40, 60]
#ptbins = [0.0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1.0, 1.25, 1.5, 2.5]
xbins  = [0.001, 0.02, 0.05, 0.1, 0.2, 0.4, 1 , 2.2]

PTBins = [0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1, 1.25, 1.5, 2.5, 4.5]
PzBins = [1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20]
EmuBins = [2,3,4,5,7,9,12,15,20] #[/* 0,1, */2,3,4,5,7,9,12,15,20]
EnuBins = [ 2.0, 3.0, 4.0, 5.0, 6.25, 7.5, 8.75, 10., 12.5, 15., 20.]
bjorkenXbins = [0.001, 0.02, 0.05, 0.1, 0.2, 0.4, 1 , 2.2]
bjorkenYbins = [0.001, 0.02, 0.05, 0.1, 0.2, 0.4, 1 , 2.2]
Erecoilbins = [0, 0.2, 0.75, 1.5, 2.5, 3.5, 5, 8, 11, 14.5, 20]

# Flux correction function in C++
flux_code = r"""
double flux_corr(double enu) {
    if (enu < 7.5) return 1.0;  // no correction below 7.5 GeV

    double ratio = 1.0;

    if (enu < 4.5) {
        ratio = 1.095;
    }
    else if (enu < 12.78) {
        // Polynomial shape from 4.5 to 12.78 GeV
        double enupoly[6] = {
            -4.11608956e+00, 3.66971952e+00, -9.73158862e-01,
             1.21054305e-01, -7.13034149e-03, 1.60919755e-04
        };
        double tempratio = enupoly[0];
        double powenu = enu;
        for (int i = 1; i < 6; ++i) {
            tempratio += enupoly[i] * powenu;
            powenu *= enu;
        }
        ratio = tempratio;
    }
    else if (enu < 26.0) {
        ratio = 1.17;
    }
    else if (enu < 29.0) {
        ratio = 1.17 + (1.09 - 1.17) * (enu - 26.0) / (29.0 - 26.0);
    }
    else if (enu < 48.0) {
        ratio = 1.09;
    }
    else if (enu < 80.0) {
        ratio = 1.09 + (2.5 - 1.09) * (enu - 48.0) / (80.0 - 48.0);
    }
    else if (enu < 100.0) {
        ratio = 2.5 + (1.5 - 2.5) * (enu - 80.0) / (100.0 - 80.0);
    }
    else {
        ratio = 1.5;
    }

    return ratio;
}
"""
ROOT.gInterpreter.Declare(flux_code)

# Define columns for selection
if pdg > 0:
    df = df.Filter("PDGnu == 14 && PDGLep == 13") #Incoming numu and outgoing mu
elif pdg < 0:
    df = df.Filter("PDGnu == -14 && PDGLep == -13")  #Incoming antinumu and outgoing antimu
df = (df.Filter("ELep > 2 && ELep < 20")          # Muon energy cut
        .Filter("acos(CosLep)*180/3.14159 < 17") # Muon angle cut
        .Define("P", "sqrt(ELep*ELep - 0.10566*0.10566)")
        .Define("Pt", "sqrt(1 - CosLep*CosLep)*P")
        .Define("Pz", "CosLep*P")
        #.Define("M_nucleon", "(939.5654133+938.2720813)/2000.0") # Get actual nucleon
        .Define("M_nucleon", 
                "pdg_vert[1] == 2112 ? 939.5654133 / 1000.0 : "
                "(pdg_vert[1] == 2212 ? 938.2720813 / 1000.0 : "
                "(939.5654133 + 938.2720813) / 2000.0)")
        .Define("xcalc", "Q2 / (2*(Enu_true - ELep)*M_nucleon)")
        .Define("ycalc", "((Enu_true - ELep) / Enu_true)")
        .Define("weight_nom", "fScaleFactor")
        .Define("weight_rw", "fScaleFactor * flux_corr(Enu_true)")
)

# Book histograms (weighted & unweighted)
h_pt_nom = df.Histo1D(( "pTmu_nom","pTmu Nominal",len(PTBins)-1,array.array("d",PTBins) ), "Pt", "weight_nom")
h_pt_rw  = df.Histo1D(( "pTmu_rw", "pTmu Reweighted",len(PTBins)-1,array.array("d",PTBins) ), "Pt", "weight_rw")

h_pz_nom = df.Histo1D(( "pZmu_nom","pZmu Nominal",len(PzBins)-1,array.array("d",PzBins) ), "Pz", "weight_nom")
h_pz_rw  = df.Histo1D(( "pZmu_rw", "pZmu Reweighted",len(PzBins)-1,array.array("d",PzBins) ), "Pz", "weight_rw")

h_x_nom  = df.Histo1D(( "BjorkenX_nom","X Nominal",len(bjorkenXbins)-1,array.array("d",bjorkenXbins) ), "xcalc", "weight_nom")
h_x_rw   = df.Histo1D(( "BjorkenX_rw", "X Reweighted",len(bjorkenXbins)-1,array.array("d",bjorkenXbins) ), "xcalc", "weight_rw")

h_y_nom  = df.Histo1D(( "BjorkenY_nom","Y Nominal",len(bjorkenXbins)-1,array.array("d",bjorkenXbins) ), "xcalc", "weight_nom")
h_y_rw   = df.Histo1D(( "BjorkenY_rw", "Y Reweighted",len(bjorkenXbins)-1,array.array("d",bjorkenXbins) ), "xcalc", "weight_rw")

h_Emu_nom  = df.Histo1D(( "Emu_nom","Emu Nominal",len(EmuBins)-1,array.array("d",EmuBins) ), "ELep", "weight_nom")
h_Emu_rw   = df.Histo1D(( "Emu_rw", "Emu Reweighted",len(EmuBins)-1,array.array("d",EmuBins) ), "ELep", "weight_rw")

h_Enu_nom  = df.Histo1D(( "Enu_nom","Enu Nominal",len(EnuBins)-1,array.array("d",EnuBins) ), "Enu_true", "weight_nom")
h_Enu_rw   = df.Histo1D(( "Enu_rw", "Enu Reweighted",len(EnuBins)-1,array.array("d",EnuBins) ), "Enu_true", "weight_rw")

h_Erecoil_nom  = df.Histo1D(( "Erecoil_nom","Erecoil Nominal",len(Erecoilbins)-1,array.array("d",Erecoilbins) ), "q0", "weight_nom")
h_Erecoil_rw   = df.Histo1D(( "Erecoil_rw", "Erecoil Reweighted",len(Erecoilbins)-1,array.array("d",Erecoilbins) ), "q0", "weight_rw")

h_pt_pz_nom = df.Histo2D(( "pTmu_pZmu_nom","pTmu vs pZmu Nominal",len(PTBins)-1,array.array("d",PTBins) ,len(PzBins)-1,array.array("d",PzBins) ), "Pt", "Pz", "weight_nom")
h_pt_pz_rw  = df.Histo2D(( "pTmu_pZmu_rw","pTmu vs pZmu Nominal",len(PTBins)-1,array.array("d",PTBins) ,len(PzBins)-1,array.array("d",PzBins) ), "Pt", "Pz", "weight_rw")

h_emu_recoil_nom = df.Histo2D(("Emu_Erecoil_nom", "Emu vs Erecoil", len(EmuBins) - 1,array.array("d", EmuBins), len(Erecoilbins) - 1, array.array("d", Erecoilbins)), "ELep", "q0", "weight_nom")
h_emu_recoil_rw = df.Histo2D(("Emu_Erecoil_rw", "Emu vs Erecoil", len(EmuBins) - 1, array.array("d", EmuBins), len(Erecoilbins) - 1, array.array("d", Erecoilbins)), "ELep", "q0", "weight_rw")


("emu_erecoil_nom", "Emu vs Erecoil", len(EnuBins) - 1, EnuBins, len(Erecoilbins) - 1, Erecoilbins)
print("Processing...")
df.Count().GetValue()  # Trigger event loop
print("Histograms filled!")


# Save ROOT file
outfile = ROOT.TFile(f"{out_prefix}_cross_sections.root","RECREATE")
h_pt_nom.Write()
h_pt_rw.Write()
h_pz_nom.Write()
h_pz_rw.Write()
h_x_nom.Write()
h_x_rw.Write()
h_y_nom.Write()
h_y_rw.Write()
h_Emu_nom.Write()
h_Emu_rw.Write()
h_Enu_nom.Write()
h_Enu_rw.Write()
h_Erecoil_nom.Write()
h_Erecoil_rw.Write()
h_pt_pz_nom.Write()
h_pt_pz_rw.Write()
h_emu_recoil_nom.Write()
h_emu_recoil_rw.Write()
outfile.Close()
print(f"✅ Saved ROOT file: {out_prefix}_cross_sections.root")
if draw:
    # Normalize by bin width
    all_histos = [h_pt_nom, h_pt_rw, h_x_nom, h_x_rw]

    # Normalize by number of trees
    if ntrees > 0:
        for h in all_histos:
            h.Scale(1.0 / ntrees)
    # Normalize histos for plotting
    for h in all_histos:
        for b in range(1, h.GetNbinsX() + 1):
            bw = h.GetBinWidth(b)
            h.SetBinContent(b, h.GetBinContent(b) / bw)
            h.SetBinError(b, h.GetBinError(b) / bw)

    # Convert to numpy for plotting
    def hist_to_numpy(h):
        nbins = h.GetNbinsX()
        bins = np.array([h.GetBinLowEdge(i) for i in range(1, nbins+2)])
        vals = np.array([h.GetBinContent(i) for i in range(1, nbins+1)])
        return bins, vals

    bins_pt, pt_nom = hist_to_numpy(h_pt_nom)
    _, pt_rw = hist_to_numpy(h_pt_rw)
    bins_x, x_nom = hist_to_numpy(h_x_nom)
    _, x_rw = hist_to_numpy(h_x_rw)

    bins_pz, pz_nom = hist_to_numpy(h_pz_nom)
    _, pz_rw = hist_to_numpy(h_pz_rw)

    ratio_pt = np.divide(pt_rw, pt_nom, out=np.zeros_like(pt_rw), where=pt_nom != 0)
    ratio_x  = np.divide(x_rw, x_nom, out=np.zeros_like(x_rw), where=x_nom != 0)

    # Plot
    import matplotlib.gridspec as gridspec
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(12, 6), constrained_layout=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=[3, 1], hspace=0.05, wspace=0.05, figure=fig)

    ax1 = fig.add_subplot(gs[0, 0])
    ax1.step(bins_pt[:-1], pt_nom, where="mid", label="Nominal", color="blue")
    ax1.step(bins_pt[:-1], pt_rw, where="mid", label="Reweighted", color="red")
    ax1.set_ylabel(r"$d\sigma/dp_T$")
    ax1.set_title("Muon $p_T$ with Low-$nu$ Reweighting")
    ax1.legend()
    ax1.tick_params(labelbottom=False)

    """ ax1 = fig.add_subplot(gs[0, 0])
    ax1.step(bins_pz[:-1], pz_nom, where="mid", label="Nominal", color="blue")
    ax1.step(bins_pz[:-1], pz_rw, where="mid", label="Reweighted", color="red")
    ax1.set_ylabel(r"$d\sigma/dp_Z$")
    ax1.set_title("Muon $p_Z$ with Low-$nu$ Reweighting")
    ax1.legend()
    ax1.tick_params(labelbottom=False) """

    ax2 = fig.add_subplot(gs[0, 1])
    ax2.step(bins_x[:-1], x_nom, where="mid", label="Nominal", color="blue")
    ax2.step(bins_x[:-1], x_rw, where="mid", label="Reweighted", color="red")
    ax2.set_ylabel(r"$d\sigma/dx$")
    ax2.set_title("Bjorken $x$ with Low-$nu$ Reweighting")
    ax2.legend()
    ax2.tick_params(labelbottom=False)

    ax3 = fig.add_subplot(gs[1, 0], sharex=ax1)
    ax3.step(bins_pt[:-1], ratio_pt, where="mid", color="green")
    ax3.axhline(1.0, color="black", linestyle="--")
    ax3.set_xlabel(r"$p_T$ [GeV/c]")
    ax3.set_ylabel("RW/Nom")

    ax4 = fig.add_subplot(gs[1, 1], sharex=ax2)
    ax4.step(bins_x[:-1], ratio_x, where="mid", color="green")
    ax4.axhline(1.0, color="black", linestyle="--")
    ax4.set_xlabel(r"Bjorken $x$")
    ax4.set_ylabel("RW/Nom")

    plt.savefig(f"{out_prefix}_cross_section_comparison.png", dpi=300)
    plt.close()
    print(f"✅ Saved plot: {out_prefix}_cross_section_comparison.png")

    avg_nom = df.Define("w", "fScaleFactor").Mean("w").GetValue()
    avg_rw = df.Define("w", "fScaleFactor * flux_corr(Enu_true)").Mean("w").GetValue()
    print(f"Average nominal weight: {avg_nom}")
    print(f"Average reweighted weight: {avg_rw}")