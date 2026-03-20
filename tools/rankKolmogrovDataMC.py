import ROOT
from tqdm import tqdm

def get_ranked_variables(data_txt, mc_txt, tree_name):
    ROOT.gROOT.SetBatch(True)
    ROOT.gErrorIgnoreLevel = ROOT.kWarning

    ch_data = ROOT.TChain(tree_name)
    ch_mc = ROOT.TChain(tree_name)
    
    with open(data_txt) as f:
        for line in f: ch_data.Add(line.strip())
    with open(mc_txt) as f:
        for line in f: ch_mc.Add(line.strip())

    results = []
    data_branches = list(ch_data.GetListOfBranches())
    mc_branches = list(ch_mc.GetListOfBranches())

    mc_branches_names = {obj.GetName() for obj in mc_branches}
    #Only look over branches common to both mc and data
    branches = [obj.GetName() for obj in data_branches if obj.GetName() in mc_branches_names]

    #branches = branches[:100]
    print("Branches to investigate:", branches)
    # Disable global stat box to save a tiny bit of rendering memory
    ROOT.gStyle.SetOptStat(0)

    for branch in tqdm(branches, desc="Scanning Branches", unit="var"):
        b_name = branch#GetName()
        #if "vector" in branch.GetClassName(): continue 
        #print("Branch: ", b_name)

        h_name_data = f"h_data_{b_name}"
        h_name_mc = f"h_mc_{b_name}"
        
        # --- STEP 1: Draw Data ---
        # "goff" = graphics off
        ch_data.Draw(f"{b_name}>>{h_name_data}", "", "goff")
        
        # Retrieve from directory
        h_data = ROOT.gDirectory.Get(h_name_data)
        
        # If it doesn't exist or is empty, skip
        if not h_data or h_data.GetEntries() == 0: 
            # If it exists in directory but is empty, remove it to prevent buildup
            if h_data: 
                ROOT.gDirectory.Delete(f"{h_name_data};*")
            continue

        # *** THE FIX: Detach from ROOT global management ***
        h_data.SetDirectory(0)
        
        # --- STEP 2: Clone for MC ---
        # Clone creates a new object with same bins
        h_mc = h_data.Clone(h_name_mc)

        """ axis = h_data.GetXaxis()
        nbins = axis.GetNbins()

        edges = [axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins() + 1)]
        edges.append(axis.GetBinUpEdge(axis.GetNbins())) # Add the final upper edge
        print("data: ", edges)

        axis = h_mc.GetXaxis()
        nbins = axis.GetNbins()

        edges = [axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins() + 1)]
        edges.append(axis.GetBinUpEdge(axis.GetNbins())) # Add the final upper edge
        print("mc: ", edges) 

        h_mc.Reset() # Empty it
        h_mc.SetDirectory(0) # *** THE FIX: Detach MC too ***

        axis = h_mc.GetXaxis()
        nbins = axis.GetNbins()

        edges = [axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins() + 1)]
        edges.append(axis.GetBinUpEdge(axis.GetNbins())) # Add the final upper edge
        print("mc after: ", edges)  """

        # --- STEP 3: Fill MC ---
        # We project into the EXISTING h_mc object
        ch_mc.Draw(f"{b_name}>>{h_name_mc}", "", "goff")
        #ch_mc.Project(h_name_mc, b_name, "")
        h_mc = ROOT.gDirectory.Get(h_name_mc)
        
        """ axis = h_data.GetXaxis()
        nbins = axis.GetNbins()

        edges = [axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins() + 1)]
        edges.append(axis.GetBinUpEdge(axis.GetNbins())) # Add the final upper edge
        print(edges)

        axis = h_mc.GetXaxis()
        nbins = axis.GetNbins()

        edges = [axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins() + 1)]
        edges.append(axis.GetBinUpEdge(axis.GetNbins())) # Add the final upper edge
        print(edges) """

        cont = [h_data.GetBinContent(i) for i in range(1, h_data.GetXaxis().GetNbins() + 1)]
        #print("data: ", cont)

        cont = [h_mc.GetBinContent(i) for i in range(1, h_mc.GetXaxis().GetNbins() + 1)]
        #print("mc: ", cont)

        
        # --- STEP 4: Test ---
        intedata = h_data.Integral()
        inte = h_mc.Integral()
        if h_mc.Integral() > 0:
            h_mc.Scale(h_data.Integral() / h_mc.Integral())
            
            # Use 'M' option for Kolmogorov (return max deviation)
            ks_score = h_data.KolmogorovTest(h_mc, "M")
            results.append((b_name, ks_score))

        can = ROOT.TCanvas()
        h_data.SetLineColor(3)
        h_data.Draw("HIST")
        h_mc.SetLineColor(2)
        h_mc.Draw("HIST SAME")
        can.SaveAs(f"{b_name}.png")
        #break
        # NO h_data.Delete() HERE! 
        # Python will garbage collect h_data and h_mc when the loop restarts.
        # Since we used SetDirectory(0), ROOT won't try to double-delete them.

    results.sort(key=lambda x: x[1])
    
    print(f"\n{'Variable':<30} | {'KS Score':<10}")
    print("-" * 45)
    for name, score in results:
        print(f"{name:<60} | {score:.6f}")

# Usage
get_ranked_variables("/exp/minerva/app/users/alhart/MAT_AL9/MINERvA-101-Cross-Section/PlaylistFiles/me-playlists/DataP7/Test1A-Data.txt", "/exp/minerva/app/users/alhart/MAT_AL9/MINERvA-101-Cross-Section/PlaylistFiles/me-playlists/MCP7/Test1A-MC.txt", "MasterAnaDev")