import ROOT
import PlotUtils
import argparse
import sys
import os

PLOTUTILSROOT = os.path.expandvars("$PLOTUTILSROOT")
plotter = PlotUtils.MnvPlotter()


# 1. Setup Argument Parsing
parser = argparse.ArgumentParser(description="Plot the ratio of a histogram from two ROOT files.")
parser.add_argument("file1", help="Path to the first ROOT file (Numerator)")
parser.add_argument("file2", help="Path to the second ROOT file (Denominator)")
parser.add_argument("-n", "--name", required=True, help="Name of the histogram (must be same in both files)")
parser.add_argument("-o", "--output", default="ratio_plot.pdf", help="Output file name (default: ratio_plot.pdf)")

args = parser.parse_args()
print("Here1")
# 2. Basic Validation
if not os.path.isfile(args.file1) or not os.path.isfile(args.file2):
    print("Error: One or both input files do not exist.")
    sys.exit(1)

# 3. ROOT Setup
ROOT.gROOT.SetBatch(True) # Run in batch mode (no window pops up)
ROOT.gStyle.SetOptStat(0) # Remove statistics box for cleaner plot

# 4. Load Files and Histograms
f1 = ROOT.TFile.Open(args.file1)
f2 = ROOT.TFile.Open(args.file2)

h1 = f1.Get(args.name)
h2 = f2.Get(args.name)

# Check if histograms exist
if not h1 or not h2:
    print(f"Error: Histogram '{args.name}' not found in one of the files.")
    sys.exit(1)

# Detach histograms from the open files so they persist after files might close
h1.SetDirectory(0)
h2.SetDirectory(0)

# Close files to keep memory clean
f1.Close()
f2.Close()

# 5. Styling & Calculation
# Ensure errors are calculated correctly (Sum of weights squared)
#h1.Sumw2()
#h2.Sumw2()

h1.SetTitle(f"Ratio of {args.name}")
h1.SetLineColor(ROOT.kRed)
h2.SetLineColor(ROOT.kBlue)

# 6. Plotting
c1 = ROOT.TCanvas("c1", "Ratio Plot", 800, 800)
plotter.DrawDataMCRatio(h1, h2, 1, True, 0.9, 1.1, "P6/P7 Cross Section Ratio for Muon Transverse Momentum on Iron")
c1.SaveAs(args.output)
# TRatioPlot