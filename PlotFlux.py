import argparse
import sys
import ROOT

def main(ref_file, test_files, hist_name, output_file):
    # Run in batch mode to prevent X11 windows from opening
    ROOT.gROOT.SetBatch(True)
    
    # 1. Extract the reference histogram
    f_ref = ROOT.TFile.Open(ref_file, "READ")
    if not f_ref or f_ref.IsZombie():
        sys.exit(f"Error: Could not open reference file '{ref_file}'")
        
    h_ref = f_ref.Get(hist_name)
    if not h_ref:
        sys.exit(f"Error: Could not find TH1D '{hist_name}' in '{ref_file}'")
        
    # Detach from file directory so it survives after the file closes
    h_ref.SetDirectory(0) 
    f_ref.Close()

    # 2. Setup canvas and legend
    c = ROOT.TCanvas("c", "Ratio Plot", 800, 600)
    c.SetGridy()
    
    # Adjust legend coordinates as needed
    legend = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
    legend.SetBorderSize(0)

    # 3. Process test files and calculate ratios
    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2]
    ratios = []
    open_files = [] # Keep files open until drawing is done

    for i, file_path in enumerate(test_files):
        f_test = ROOT.TFile.Open(file_path, "READ")
        open_files.append(f_test)
        
        h_test = f_test.Get(hist_name)
        if not h_test:
            print(f"Warning: Could not find '{hist_name}' in '{file_path}'. Skipping.")
            continue
            
        # Clone to create the ratio histogram
        h_ratio = h_test.Clone(f"ratio_{i}")
        h_ratio.SetDirectory(0)
        
        # Divide test by reference (automatically handles error propagation)
        h_ratio.Divide(h_ratio, h_ref)

        h_ratio.Smooth(5)
        
        # Styling
        color = colors[i % len(colors)]
        h_ratio.SetLineColor(color)
        h_ratio.SetMarkerColor(color)
        h_ratio.SetMarkerStyle(20 + i)
        
        # Axis labels and title
        h_ratio.SetTitle("Flux Ratio to Reference;X-axis;Ratio")
        h_ratio.SetStats(0) # Turn off the stat box
        
        ratios.append((file_path, h_ratio))

    # 4. Draw histograms
    if not ratios:
        sys.exit("Error: No valid test histograms found to plot.")

    # We determine a sensible Y-axis range based on the data
    max_val = max([h.GetMaximum() for _, h in ratios])
    min_val = min([h.GetMinimum() for _, h in ratios])
    
    for i, (path, h_ratio) in enumerate(ratios):
        if i == 0:
            h_ratio.GetYaxis().SetRangeUser(min_val * 0.9, max_val * 1.1)
            h_ratio.Draw("HIST C") # EP = Draw error bars and markers
        else:
            h_ratio.Draw("HIST C SAME")
            
        legend.AddEntry(h_ratio, path.split('/')[-1], "lep")

    legend.Draw()
    
    # 5. Save and cleanup
    c.SaveAs(output_file)
    print(f"Plot saved successfully to {output_file}")
    
    for f in open_files:
        f.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot TH1D ratios relative to a reference file.")
    parser.add_argument("-r", "--ref", required=True, help="Path to the reference .root file")
    parser.add_argument("-f", "--files", nargs='+', required=True, help="Paths to the other .root files")
    parser.add_argument("-n", "--name", default="flux", help="Name of the histogram (default: 'flux')")
    parser.add_argument("-o", "--output", default="ratio.pdf", help="Output plot filename (default: ratio.pdf)")

    args = parser.parse_args()
    
    if len(args.files) != 3:
        print(f"Note: You provided {len(args.files)} files to compare. The script supports any number, but 3 were expected.")
        
    main(args.ref, args.files, args.name, args.output)