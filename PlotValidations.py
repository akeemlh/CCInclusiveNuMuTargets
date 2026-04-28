import ROOT, PlotUtils
import math
import os
import argparse
import numpy as np
from array import array

ROOT.gROOT.SetBatch(True) # Don't open window
ROOT.gStyle.SetPalette(ROOT.kBird)

filledplaylistsused = ""
emptyplaylistsused = ""

materials = {"Iron": 26, "Lead": 82, "Water": -999, "Carbon": 6}
TgtCodeLabels = {
    1026: "1026", 1082: "1082", 2026: "2026", 2082: "2082", 
    3006: "3006", 3026: "3026", 3082: "3082", 4082: "4082", 
    5026: "5026", 5082: "5082", -999: "Water"
}

rebinNum = 10
nbins = 3400.0
xlow = 4200.0
xhigh = 5900.0

def GetIngredient(dir_obj, ingredient):
    obj = dir_obj.Get(ingredient)
    if not obj:
        raise RuntimeError(f"Failed to get {ingredient} in {dir_obj.GetName()}")
    return obj

def find_paths(inpath, file_name, recursive_search=True, skip_test=True, curdepth=0, maxdepth=2):
    # Using a set automatically prevents duplicate directory paths
    filepaths = set()
    if os.path.isfile(inpath) and file_name in inpath: #If fed a single, valid file, return that
        return [inpath]
    try:
        with os.scandir(inpath) as it:
            for entry in it:
                if entry.is_file():
                    # The C++ code searches the full path string
                    if f"Validations" in entry.path:
                        if not (skip_test and "/Test" in inpath):
                            filepaths.add(entry.path)
                            
                # follow_symlinks=True mimics the C++ check for both directory and symlink
                elif entry.is_dir(follow_symlinks=True) and recursive_search and curdepth <= maxdepth:
                    # Recursive call
                    subdirs = find_paths(
                        entry.path, 
                        file_name, 
                        recursive_search, 
                        skip_test, 
                        curdepth + 1, 
                        maxdepth
                    )
                    filepaths.update(subdirs)
                    
    except PermissionError:
        # It's good practice in Python to handle permission errors during directory traversal
        pass

    # Convert the set back to a list and sort it to match the C++ std::sort behavior
    return sorted(list(filepaths))

NU_PLAYLISTS = {
    "minervame1A", "minervame1B", "minervame1C", "minervame1D", 
    "minervame1E", "minervame1F", "minervame1G", "minervame1L", 
    "minervame1M", "minervame1N", "minervame1O", "minervame1P"
}

ANU_PLAYLISTS = {
    "minervame5A", "minervame6A", "minervame6B", "minervame6C", 
    "minervame6D", "minervame6E", "minervame6F", "minervame6G", 
    "minervame6H", "minervame6I", "minervame6J"
}

FILLED_PLAYLISTS = {
    "minervame1L", "minervame1M", "minervame1N", "minervame1O", "minervame1P", 
    "1L", "1M", "1N", "1O", "1P", "minervame5A", "minervame6A", "minervame6B", 
    "minervame6C", "minervame6D", "minervame6E", "minervame6F", "minervame6G", 
    "minervame6H", "5A", "6A", "6B", "6C", "6D", "6E", "6F", "6G", "6H"
}

EMPTY_PLAYLISTS = {
    "minervame1A", "minervame1B", "minervame1C", "minervame1D", "minervame1E", 
    "minervame1F", "minervame1G", "1A", "1B", "1C", "1D", "1E", "1F", "1G", 
    "minervame6J", "6J"
}

def nu_or_anti_nu_mode(playlist: str) -> int:
    if playlist in NU_PLAYLISTS:
        return 1
    elif playlist in ANU_PLAYLISTS:
        return 2
    return 0

def filled_or_empty_me_playlist(playlist: str) -> int:
    print("playlist: ", playlist)
    if playlist in FILLED_PLAYLISTS:
        return 1
    elif playlist in EMPTY_PLAYLISTS:
        return 2
    return 0

def add_hist(hist1, hist2, scale=1.0):
    if hist1 is None or hist1.GetEntries() == 0:
        hist1 = hist2.Clone()
        hist1.Scale(scale)
        hist1.SetDirectory(0) # Detach from the file to prevent PyROOT from garbage collecting it
        return hist1
    # If it already has data, just add the new one
    hist1.Add(hist2, scale)
    return hist1

def add_hist(hist1, file, hist2name, scale=1.0):
    hist2 = file.Get(hist2name)
    if hist1 is None or hist1.GetEntries() == 0:
        hist1 = hist2.Clone()
        hist1.Scale(scale)
        hist1.SetDirectory(0) # Detach from the file to prevent PyROOT from garbage collecting it
        print("hist1 entries: ", hist1.GetEntries())
        return hist1
    # If it already has data, just add the new one
    hist1.Add(hist2, scale)
    return hist1

def normalize_by_row(h_orig):
    if not h_orig: return None
    
    h_norm = h_orig.Clone(f"{h_orig.GetName()}_rowNorm")
    h_norm.SetDirectory(0) # Prevent garbage collection issues
    
    nBinsX = h_norm.GetNbinsX()
    nBinsY = h_norm.GetNbinsY()
    
    for iy in range(1, nBinsY + 1):
        # Integral(x_start, x_end, y_start, y_end)
        row_sum = h_orig.Integral(1, nBinsX, iy, iy)
        
        for ix in range(1, nBinsX + 1):
            if row_sum > 0:
                h_norm.SetBinContent(ix, iy, h_orig.GetBinContent(ix, iy) / row_sum)
                h_norm.SetBinError(ix, iy, h_orig.GetBinError(ix, iy) / row_sum)
            else:
                h_norm.SetBinContent(ix, iy, 0.0)
                h_norm.SetBinError(ix, iy, 0.0)
                
    return h_norm

def normalize_by_column(h_orig):
    if not h_orig: return None
    
    h_norm = h_orig.Clone(f"{h_orig.GetName()}_colNorm")
    h_norm.SetDirectory(0) # Prevent garbage collection issues
    
    nBinsX = h_norm.GetNbinsX()
    nBinsY = h_norm.GetNbinsY()
    
    for ix in range(1, nBinsX + 1):
        # Integral(x_start, x_end, y_start, y_end)
        col_sum = h_orig.Integral(ix, ix, 1, nBinsY)
        
        for iy in range(1, nBinsY + 1):
            if col_sum > 0:
                h_norm.SetBinContent(ix, iy, h_orig.GetBinContent(ix, iy) / col_sum)
                h_norm.SetBinError(ix, iy, h_orig.GetBinError(ix, iy) / col_sum)
            else:
                h_norm.SetBinContent(ix, iy, 0.0)
                h_norm.SetBinError(ix, iy, 0.0)
                
    return h_norm

def PlotRegions(inYmin = None, inYmax = None, mcpot_input = None, datapot_input = None):
    print("Plotting regions test")
    if (inYmin == None):
        ymin = ROOT.gPad.GetFrame().GetY1()
    else:
        ymin = inYmin
    if (inYmax == None):
        ymax = ROOT.gPad.GetFrame().GetY2()
    else:
        ymax = inYmax

    tgt1centre = 4481.0 
    tgt1thickness = 25.78 
    tgt2centre = 4702.0 
    tgt2thickness = 25.81 
    tgt3centre = 4944.48 
    tgt3thickness = 76.2 
    tgt4centre = 5641.68 
    tgt4thickness = 7.95 
    tgt5centre = 5778.0 
    tgt5thickness = 13.17 

    scintillatorthickness = 10.0

    boxes = []

    boxes.append(ROOT.TBox(tgt1centre-tgt1thickness/2, ymin, tgt1centre+tgt1thickness/2, ymax))
    boxes.append(ROOT.TBox(tgt2centre-tgt2thickness/2, ymin, tgt2centre+tgt2thickness/2, ymax))
    boxes.append(ROOT.TBox(tgt3centre-tgt3thickness/2, ymin, tgt3centre+tgt3thickness/2, ymax))
    boxes.append(ROOT.TBox(tgt4centre-tgt4thickness/2, ymin, tgt4centre+tgt4thickness/2, ymax))
    boxes.append(ROOT.TBox(tgt5centre-tgt5thickness/2, ymin, tgt5centre+tgt5thickness/2, ymax))
    boxes.append(ROOT.TBox(5200, ymin, 5420, ymax))

    for i in range(5):
        boxes[i].SetFillColorAlpha(ROOT.kBlack, 0.25)
        boxes[i].SetFillStyle(3144)
        boxes[i].SetLineWidth(0)
        boxes[i].SetLineColor(1)
        boxes[i].Draw("same")
        ROOT.SetOwnership(boxes[i], False)
        
    boxes[5].SetFillColorAlpha(ROOT.TColor.GetColor("#78C3ED"), 0.25)
    boxes[5].SetFillStyle(3144)
    boxes[5].SetLineWidth(0)
    boxes[5].SetLineColor(1)
    boxes[5].Draw("same")
    ROOT.SetOwnership(boxes[5], False)

    # US Boxes
    us_boxes = [
        ROOT.TBox(4446.33-scintillatorthickness, ymin, 4446.33+scintillatorthickness, ymax),
        ROOT.TBox(4667.4-scintillatorthickness, ymin, 4667.4+scintillatorthickness, ymax),
        ROOT.TBox(4888.48-scintillatorthickness, ymin, 4888.48+scintillatorthickness, ymax),
        ROOT.TBox(5610.02-scintillatorthickness, ymin, 5610.02+scintillatorthickness, ymax),
        ROOT.TBox(5742.67-scintillatorthickness, ymin, 5742.67+scintillatorthickness, ymax),
        ROOT.TBox(5153.77-scintillatorthickness, ymin, 5153.77+scintillatorthickness, ymax)
        #ROOT.TBox(5153.77-scintillatorthickness, ymin, 5200, ymax)
    ]

    for b in us_boxes:
        b.SetFillColorAlpha(ROOT.TColor.GetColor("#D694B9"), 0.5)
        b.SetLineWidth(0)
        b.SetLineColor(1)
        b.Draw("same")
        ROOT.SetOwnership(b, False)

    # DS Boxes
    ds_boxes = [
        ROOT.TBox(4514.11-scintillatorthickness, ymin, 4514.11+scintillatorthickness, ymax),
        ROOT.TBox(4735.19-scintillatorthickness, ymin, 4735.19+scintillatorthickness, ymax),
        ROOT.TBox(5000.48-scintillatorthickness, ymin, 5000.48+scintillatorthickness, ymax),
        ROOT.TBox(5677.81-scintillatorthickness, ymin, 5677.81+scintillatorthickness, ymax),
        ROOT.TBox(5810.45-scintillatorthickness, ymin, 5810.45+scintillatorthickness, ymax),
        ROOT.TBox(5456.74-scintillatorthickness, ymin, 5456.74+scintillatorthickness, ymax)
        #ROOT.TBox(5420-scintillatorthickness, ymin, 5456.74+scintillatorthickness, ymax)
    ]

    for b in ds_boxes:
        b.SetFillColorAlpha(ROOT.TColor.GetColor("#33B18F"), 0.5)
        b.SetLineWidth(0)
        b.SetLineColor(1)
        b.Draw("same")
        ROOT.SetOwnership(b, False)

    # To avoid garbage collection issues in PyROOT, we keep references
    ROOT.SetOwnership(boxes[0], False)

    if datapot_input is not None:
        datapotstr = f"Data POT: {datapot_input:.2f}"
        potstr = f"Normalised POT: {datapot_input:.2f}"

        potlabel = ROOT.TLatex(xlow, ymax*1.04, potstr)
        potlabel.SetTextAlign(12)
        potlabel.SetTextFont(43)
        potlabel.SetTextSize(20)

    resolution = ((xhigh - xlow )/ nbins) * rebinNum
    resStr = f"Resolution/Bin Width: {resolution:.1f} mm"
    resLabel = ROOT.TLatex(xhigh, ymax*1.04, resStr)
    resLabel.SetTextAlign(32)
    resLabel.SetTextFont(43)
    resLabel.SetTextSize(0.035) 
    resLabel.Draw("SAME")
    
    titlelabel = ROOT.TLatex(5000, ymax*1.08, "True Vertex Positions of Interactions reconstructed inside the ")
    titlelabel.SetTextAlign(22)
    titlelabel.SetTextFont(43)
    titlelabel.SetTextSize(35) 
    if datapot_input is not None:
        potlabel.Draw("SAME")
    titlelabel.Draw("same")

def DrawDataMCRatio(data, MC, mcScale):
    tmpData = data.Clone()
    tmpData.Divide(data, MC)
    tmpData.Draw("P")
    ROOT.SetOwnership(tmpData, False)

def DrawDataMC(DataIn, MCIn, name, mcpot, datapot, outbase, zoomIntoTargets = False):
    canvas = ROOT.TCanvas("canvas", "canvas", 1500, 1000)
    DataHist = DataIn.Clone()
    MCHist = MCIn.Clone()
    MCHist.SetTitle("MC")
    DataHist.SetTitle("Data")
    
    MCHist.SetUserRange(800, 825)
    DataHist.SetUserRange(800, 825)
    name = name.replace(" ", "").replace("-", "")

    plotter = ROOT.PlotUtils.MnvPlotter()
    plotter.ApplyStyle(ROOT.PlotUtils.kDefaultStyle)
    
    plotter.DrawDataMC(DataHist, MCHist, datapot/mcpot, "TL", True)
    ymin = ROOT.gPad.GetFrame().GetY1()
    ymax = ROOT.gPad.GetFrame().GetY2()
    PlotRegions(ymin, ymax, mcpot, datapot)
    canvas.SaveAs(f"{outbase}{name}.pdf")
    
    canvas.Clear()
    plotter.DrawDataMCRatio(DataHist, MCHist, datapot/mcpot, True, 0, 2, "Data/MC")
    ymin = ROOT.gPad.GetFrame().GetY1()
    ymax = ROOT.gPad.GetFrame().GetY2()
    PlotRegions(ymin, ymax, mcpot, datapot)
    canvas.SaveAs(f"{outbase}{name}Ratio.pdf")

    canvas.Clear()
    canvas.SetLogy()
    DataHist.SetMinimum(10)
    MCHist.SetMinimum(10)
    print("Literally just set min to 1")
    plotter.DrawDataMC(DataHist, MCHist, datapot/mcpot, "TL", True)
    ymin = ROOT.gPad.GetFrame().GetY1()
    ymax = ROOT.gPad.GetFrame().GetY2()
    PlotRegions(10**ymin, 10**ymax, mcpot, datapot)
    print(f"LogY set: {ROOT.gPad.GetLogy()}")
    
    tmpMC = canvas.GetPrimitive("tmpMC") 
    print(f"minimum: {tmpMC.GetMinimum()}") 
    canvas.SaveAs(f"{outbase}{name}Log.pdf")


def GetHistogramQuartiles(h):
    quartiles = [0.0] * 4
    if not h: return quartiles

    total = h.Integral(0, h.GetNbinsX() + 1)  
    if total <= 0: return quartiles

    cumulative = 0.0
    nBins = h.GetNbinsX()

    targets = [0.25, 0.50, 0.75, 1.0]
    t = 0 

    for i in range(1, nBins + 1):
        if t >= 4: break
        cumulative += h.GetBinContent(i)
        frac = cumulative / total

        if frac >= targets[t]:
            xlow = h.GetBinLowEdge(i)
            xhigh = h.GetBinLowEdge(i+1)
            prevFrac = (cumulative - h.GetBinContent(i)) / total
            binFrac = frac - prevFrac

            interp = 0.0
            if binFrac > 0:
                interp = (targets[t] - prevFrac) / binFrac

            quartiles[t] = xlow + interp * (xhigh - xlow)
            t += 1

    return quartiles


""" def PlotRegions(inYmin = None, inYmax = None, mcpot = None, datapot = None):

    if (inYmin == None):
        ymin = ROOT.gPad.GetFrame().GetY1()
    else:
        ymin = inYmin
    if (inYmax == None):
        ymax = ROOT.gPad.GetFrame().GetY2()
    else:
        ymax = inYmax
    print("PlotRegions ymin: ", ymin, "ymax ", ymax)
    #These were all taken from https://cdcvs.fnal.gov/redmine/projects/minerva/wiki/Z_Positions_of_Planes_in_the_Full_MINERvA_Detector
    #For target 3 the centre was taken as the mid point between the US and DS planes, the stated centre gives weird looking results I suspect because it's a 3 material target the reported centre may not be actual centre
    tgt1centre = 4481
    tgt1thickness = 25.78
    tgt2centre = 4702
    tgt2thickness = 25.81
    tgt3centre = 4944.48
    tgt3thickness = 76.2
    tgt4centre = 5641.68
    tgt4thickness = 7.95
    tgt5centre = 5778
    tgt5thickness = 13.17

    tgt1 = ROOT.TBox(tgt1centre-tgt1thickness/2, ymin, tgt1centre+tgt1thickness/2, ymax)
    tgt2 = ROOT.TBox(tgt2centre-tgt2thickness/2, ymin, tgt2centre+tgt2thickness/2, ymax)
    tgt3 = ROOT.TBox(tgt3centre-tgt3thickness/2, ymin, tgt3centre+tgt3thickness/2, ymax)
    tgt4 = ROOT.TBox(tgt4centre-tgt4thickness/2, ymin, tgt4centre+tgt4thickness/2, ymax)
    tgt5 = ROOT.TBox(tgt5centre-tgt5thickness/2, ymin, tgt5centre+tgt5thickness/2, ymax)
    tgt6 = ROOT.TBox(5200, ymin, 5420, ymax)

    tgt1.SetFillColorAlpha(ROOT.kBlack,0.5)
    tgt1.SetFillStyle(3004)
    tgt1.SetLineWidth(0)
    tgt1.SetLineColor(1)
    tgt1.DrawBox(tgt1centre-tgt1thickness/2, ymin, tgt1centre+tgt1thickness/2, ymax)
    tgt2.SetFillColorAlpha(ROOT.kBlack,0.5)
    tgt2.SetFillStyle(3004)
    tgt2.SetLineWidth(0)
    tgt2.SetLineColor(1)
    tgt2.DrawBox(tgt2centre-tgt2thickness/2, ymin, tgt2centre+tgt2thickness/2, ymax)
    tgt3.SetFillColorAlpha(ROOT.kBlack,0.5)
    tgt3.SetFillStyle(3004)
    tgt3.SetLineWidth(0)
    tgt3.SetLineColor(1)
    tgt3.DrawBox(tgt3centre-tgt3thickness/2, ymin, tgt3centre+tgt3thickness/2, ymax)
    tgt4.SetFillColorAlpha(ROOT.kBlack,0.5)
    tgt4.SetFillStyle(3004)
    tgt4.SetLineWidth(0)
    tgt4.SetLineColor(1)
    tgt4.DrawBox(tgt4centre-tgt4thickness/2, ymin, tgt4centre+tgt4thickness/2, ymax)
    tgt5.SetFillColorAlpha(ROOT.kBlack,0.5)
    tgt5.SetFillStyle(3004)
    tgt5.SetLineWidth(0)
    tgt5.SetLineColor(1)
    tgt5.DrawBox(tgt5centre-tgt5thickness/2, ymin, tgt5centre+tgt5thickness/2, ymax)
    tgt6.SetFillColorAlpha(ROOT.kBlue,0.5)
    tgt6.SetFillStyle(3004)
    tgt6.SetLineWidth(0)
    tgt6.SetLineColor(1)
    tgt6.DrawBox(5200, ymin, 5420, ymax)
    print("Finished plotting regions") """


def DrawFullAndEmptyDataMC(canvas, DataInFull, MCInFull, DataInEmpty, MCInEmpty, title, plotType, filledDataPOT = 0, filledMCPOT = 0, emptyDataPOT = 0, emptyMCPOT = 0, yminFull = -1.0, ymaxFull = -1.0, yminEmpty = -1.0, ymaxEmpty = -1.0, titlesize= 0.032, zoomIntoTargets = False):

        plotterfull = PlotUtils.MnvPlotter()
        plotterempty = PlotUtils.MnvPlotter()
        tmpDataFullHist = DataInFull.Clone()
        tmpMCFullHist = MCInFull.Clone()
        tmpDataEmptyHist = DataInEmpty.Clone()
        tmpMCEmptyHist = MCInEmpty.Clone()



        if filledDataPOT > 1e21:
            filleddatapotstr = f"{(filledDataPOT/1e21):.2f}"+" #times 10^{21}"
        if filledDataPOT > 1e20:
            filleddatapotstr = f"{(filledDataPOT/1e20):.2f}"+" #times 10^{20}"
        if filledDataPOT > 1e19:
            filleddatapotstr = f"{(filledDataPOT/1e19):.2f}"+" #times 10^{19}"
        if emptyDataPOT > 1e21:
            emptydatapotstr = f"{(emptyDataPOT/1e21):.2f}"+" #times 10^{21}"
        if emptyDataPOT > 1e20:
            emptydatapotstr = f"{(emptyDataPOT/1e20):.2f}"+" #times 10^{20}"
        if emptyDataPOT > 1e19:
            emptydatapotstr = f"{(emptyDataPOT/1e19):.2f}"+" #times 10^{19}"
        if filledMCPOT > 1e21:
            filledMCpotstr = f"{(filledMCPOT/1e21):.2f}"+" #times 10^{21}"
        if filledMCPOT > 1e20:
            filledMCpotstr = f"{(filledMCPOT/1e20):.2f}"+" #times 10^{20}"
        if filledMCPOT > 1e19:
            filledMCpotstr = f"{(filledMCPOT/1e19):.2f}"+" #times 10^{19}"
        if emptyMCPOT > 1e21:
            emptyMCpotstr = f"{(emptyMCPOT/1e21):.2f}"+" #times 10^{21}"
        if emptyMCPOT > 1e20:
            emptyMCpotstr = f"{(emptyMCPOT/1e20):.2f}"+" #times 10^{20}"
        if emptyMCPOT > 1e19:
            emptyMCpotstr = f"{(emptyMCPOT/1e19):.2f}"+" #times 10^{19}"
        for hist in [tmpDataFullHist, tmpMCFullHist, tmpDataEmptyHist, tmpMCEmptyHist]:
            if zoomIntoTargets:
                hist.GetXaxis().SetRangeUser(-1, 6000)
            else:
                hist.GetXaxis().SetRangeUser(-1, 8300)
            hist.GetXaxis().SetNoExponent(True)

        if (yminFull != -1.0):
            print("here1")
            plotterfull.axis_minimum = yminFull
        if (ymaxFull != -1.0):
            print("here3")
            plotterfull.axis_maximum = ymaxFull
        if (yminEmpty != -1.0):
            print("here3")
            plotterempty.axis_minimum = yminEmpty
        if (ymaxEmpty != -1.0):
            print("here4")
            plotterempty.axis_maximum = ymaxEmpty
        tmpMCFullHist.SetTitle("MC")
        tmpDataFullHist.SetTitle("Data")
        tmpMCEmptyHist.SetTitle("MC")
        tmpDataEmptyHist.SetTitle("Data")
        #Drawing split plots, with both filled and empty in the same plot
        #Preparing canvas
        ROOT.gPad.Update()
        canvas.SetLogy(False); #Unsetting log axis
        p1 = ROOT.TPad("p1","top",0,0.45,1,0.9) # xlow,yloy, xup,yup
        p2 = ROOT.TPad("p2","bottom",0,0,1,0.45)
        tex3 = ROOT.TLatex(.5,.95,"Title2")
        tex3.SetTextAlign(22)
        tex3.SetTextSize(titlesize)
        tex3.DrawLatexNDC(.5,.95, title)
        #Draw linear plots
        canvas.cd()
        p1.SetFillStyle(0)
        p2.SetFillStyle(0)
        p1.Draw()
        p1.cd()
        plotterfull.ratio_marker_size = 0.5
        plotterempty.ratio_marker_size = 0.5
        plotterfull.ratio_line_width = 0
        plotterempty.ratio_line_width = 0
        
        if (plotType == "lin"): #Draw linear plots
            #Draw filled linear plot
            print("Drawing linear plots")
            plotterfull.DrawDataMC(tmpDataFullHist, tmpMCFullHist, filledDataPOT/filledMCPOT, "TL", True)
            ymin = ROOT.gPad.GetFrame().GetY1()
            ymax = ROOT.gPad.GetFrame().GetY2()
            PlotRegions()
            axistext = ROOT.TLatex(.5,.95,"axistext1")
            axistext.SetTextAlign(22)
            axistext.SetTextSize(0.05)
            axistext.DrawLatexNDC(.85,.05, "Vertex Z Position (mm)")
            if (filledDataPOT!=0 and filledMCPOT!=0):
                filledPOTtxt = ROOT.TLatex(.5,.95,"filledPOTtext1")
                filledPOTtxt.SetTextAlign(22)
                filledPOTtxt.SetTextSize(0.05)
                filledPOTtxt.DrawLatexNDC(.5,.95, f"Water Filled - Normalised POT: {filleddatapotstr}")
            canvas.cd()
            p2.Draw()
            p2.cd()
            #Draw empty linear plot
            plotterempty.DrawDataMC(tmpDataEmptyHist, tmpMCEmptyHist, emptyDataPOT/emptyMCPOT, "TL", True)
            ymin = ROOT.gPad.GetFrame().GetY1()
            ymax = ROOT.gPad.GetFrame().GetY2()
            print("ymin: ", ymin, " ymax: ", ymax)
            PlotRegions()
            axistext = ROOT.TLatex(.5,.95,"axistext2")
            axistext.SetTextAlign(22)
            axistext.SetTextSize(0.05)
            axistext.DrawLatexNDC(.85,.05, "Vertex Z Position (mm)")
            if (emptyDataPOT!=0 and emptyMCPOT!=0):
                emptyPOTtxt = ROOT.TLatex(.5,.95,"emptyPOTtext1")
                emptyPOTtxt.SetTextAlign(22)
                emptyPOTtxt.SetTextSize(0.05)
                emptyPOTtxt.DrawLatexNDC(.5,.95, f"Water Empty - Normalised POT: {emptydatapotstr}")
        elif (plotType == "rat"): #Draw ratio plots
            #Draw filled ratio plot
            print("Drawing ratio plots")
            #DrawDataMCRatio(tmpDataFullHist, tmpMCFullHist, filledDataPOT/filledMCPOT)
            plotterfull.DrawDataMCRatio(tmpDataFullHist, tmpMCFullHist, filledDataPOT/filledMCPOT, True, 0.5, 1.5, "")
            ymin = ROOT.gPad.GetFrame().GetY1()
            ymax = ROOT.gPad.GetFrame().GetY2()
            PlotRegions()
            axistext1 = ROOT.TLatex(.5,.95,"axistext1")
            axistext1.SetTextAlign(22)
            axistext1.SetTextSize(0.05)
            axistext1.DrawLatexNDC(.85,.05, "Vertex Z Position (mm)")
            axistext2 = ROOT.TLatex(.5,.95,"axistext2")
            axistext2.SetTextAlign(22)
            axistext2.SetTextSize(0.05)
            axistext2.SetTextAngle(90)
            axistext2.DrawLatexNDC(.1,.5, "Data/MC")
            if (filledDataPOT!=0 and filledMCPOT!=0):
                filledPOTtxt = ROOT.TLatex(.5,.95,"filledPOTtext1")
                filledPOTtxt.SetTextAlign(22)
                filledPOTtxt.SetTextSize(0.05)
                filledPOTtxt.DrawLatexNDC(.5,.95, f"Water Filled - Normalised POT: {filleddatapotstr}")
            #Draw empty ratio plot
            canvas.cd()
            p2.Draw()
            p2.cd()
            #DrawDataMCRatio(tmpDataEmptyHist, tmpMCEmptyHist, emptyDataPOT/emptyMCPOT)
            plotterempty.DrawDataMCRatio(tmpDataEmptyHist, tmpMCEmptyHist, emptyDataPOT/emptyMCPOT, True, 0.5, 1.5, "")
            ymin = ROOT.gPad.GetFrame().GetY1()
            ymax = ROOT.gPad.GetFrame().GetY2()
            PlotRegions(ymin, ymax)
            axistext1 = ROOT.TLatex(.5,.95,"axistext1")
            axistext1.SetTextAlign(22)
            axistext1.SetTextSize(0.05)
            axistext1.DrawLatexNDC(.85,.05, "Vertex Z Position (mm)")
            axistext2 = ROOT.TLatex(.5,.95,"axistext2")
            axistext2.SetTextAlign(22)
            axistext2.SetTextSize(0.05)
            axistext2.SetTextAngle(90)
            axistext2.DrawLatexNDC(.1,.5, "Data/MC")
            if (emptyDataPOT!=0 and emptyMCPOT!=0):
                emptyPOTtxt = ROOT.TLatex(.5,.95,"emptyPOTtext1")
                emptyPOTtxt.SetTextAlign(22)
                emptyPOTtxt.SetTextSize(0.05)
                emptyPOTtxt.DrawLatexNDC(.5,.95, f"Water Empty - Normalised POT: {emptydatapotstr}")
        elif (plotType == "log"): #Draw ratio plots
            #Draw filled log plot
            print("Drawing log plots")
            p1.SetLogy()

            if (yminFull == -1.0):
                tmpDataFullHist.SetMinimum(1)
                tmpMCFullHist.SetMinimum(1)

            plotterfull.DrawDataMC(tmpDataFullHist, tmpMCFullHist, filledDataPOT/filledMCPOT, "TL", True)
            ymin = ROOT.gPad.GetFrame().GetY1()
            ymax = ROOT.gPad.GetFrame().GetY2()
            #std::cout << "Filled ymin: " << ymin << " ymax: " << ymax << std::endl
            PlotRegions(np.power(10,ymin), np.power(10,ymax), 0, 0)
            axistext = ROOT.TLatex(.5,.95,"axistext1")
            axistext.SetTextAlign(22)
            axistext.SetTextSize(0.05)
            axistext.DrawLatexNDC(.85,.05, "Vertex Z Position (mm)")
            if (filledDataPOT!=0 and filledMCPOT!=0):
                filledPOTtxt = ROOT.TLatex(.5,.95,"filledPOTtext1")
                filledPOTtxt.SetTextAlign(22)
                filledPOTtxt.SetTextSize(0.05)
                filledPOTtxt.DrawLatexNDC(.5,.95, f"Water Filled - Normalised POT: {filleddatapotstr}")
            #Draw empty log plot
            canvas.cd()
            p2.Draw()
            p2.cd()
            p2.SetLogy()
            #canvas.SetLogy()
            if (yminEmpty == -1.0):
                tmpDataEmptyHist.SetMinimum(1)
                tmpMCEmptyHist.SetMinimum(1)

            plotterempty.DrawDataMC(tmpDataEmptyHist, tmpMCEmptyHist, emptyDataPOT/emptyMCPOT, "TL", True)
            ymin = ROOT.gPad.GetFrame().GetY1()
            ymax = ROOT.gPad.GetFrame().GetY2()
            #std::cout << "Empty ymin: " << ymin << " ymax: " << ymax << std::endl
            PlotRegions(np.power(10,ymin), np.power(10,ymax), 0, 0)
            axistext = ROOT.TLatex(.5,.95,"axistext2")
            axistext.SetTextAlign(22)
            axistext.SetTextSize(0.05)
            axistext.DrawLatexNDC(.85,.05, "Vertex Z Position (mm)")
            if (emptyDataPOT!=0 and emptyMCPOT!=0):
                emptyPOTtxt = ROOT.TLatex(.5,.95,"emptyPOTtext1")
                emptyPOTtxt.SetTextAlign(22)
                emptyPOTtxt.SetTextSize(0.05)
                emptyPOTtxt.DrawLatexNDC(.5,.95, f"Water Empty - Normalised POT: {emptydatapotstr}")
        canvas.cd()
        return


def DrawSinglePlotsDataMC(canvas, DataIn, MCIn, title, plotType, DataPOT = 0, MCPOT = 0, ymin = -1.0, ymax = -1.0, titlesize= 0.032):

        plotter = PlotUtils.MnvPlotter()
        tmpDataHist = DataIn.Clone()
        tmpMCHist = MCIn.Clone()
        if (ymin != -1.0):
            plotter.axis_minimum = ymin
        if (ymax != -1.0):
            plotter.axis_maximum = ymax
        tmpMCHist.SetTitle("MC")
        tmpDataHist.SetTitle("Data")
        #Drawing split plots, with both filled and empty in the same plot
        #Preparing canvas
        ROOT.gPad.Update()
        canvas.SetLogy(False); #Unsetting log axis
        p1 = ROOT.TPad("p1","main",0,0,1,0.9) # xlow,yloy, xup,yup
        tex3 = ROOT.TLatex(.5,.95,"Title2")
        tex3.SetTextAlign(22)
        tex3.SetTextSize(titlesize)
        tex3.DrawLatexNDC(.5,.95, title)
        #Draw linear plots
        canvas.cd()
        p1.SetFillStyle(0)
        p1.Draw()
        p1.cd()
        if (plotType == "lin"): #Draw linear plots
            #Draw filled linear plot
            print("Drawing linear plots")
            plotter.DrawDataMC(tmpDataHist, tmpMCHist, 1.0, "TL", True)
            ymin = ROOT.gPad.GetFrame().GetY1()
            ymax = ROOT.gPad.GetFrame().GetY2()
            PlotRegions()
            axistext = ROOT.TLatex(.5,.95,"axistext1")
            axistext.SetTextAlign(22)
            axistext.SetTextSize(0.05)
            axistext.DrawLatexNDC(.85,.05, "Vertex Z Position (mm)")
            if (DataPOT!=0 and MCPOT!=0):
                filledPOTtxt = ROOT.TLatex(.5,.95,"filledPOTtext1")
                filledPOTtxt.SetTextAlign(22)
                filledPOTtxt.SetTextSize(0.05)
                filledPOTtxt.DrawLatexNDC(.5,.95, "Normalised POT: " + "{0:e}".format(DataPOT))
            canvas.cd()
        elif (plotType == "rat"): #Draw ratio plots
            #Draw filled ratio plot
            print("Drawing ratio plots")
            plotter.DrawDataMCRatio(tmpDataHist, tmpMCHist, 1.0, True,ymin, ymax, "")
            ymin = ROOT.gPad.GetFrame().GetY1()
            ymax = ROOT.gPad.GetFrame().GetY2()
            PlotRegions()
            axistext = ROOT.TLatex(.5,.95,"axistext1")
            axistext.SetTextAlign(22)
            axistext.SetTextSize(0.05)
            axistext.DrawLatexNDC(.85,.05, "Vertex Z Position (mm)")
            if (DataPOT!=0 and MCPOT!=0):
                filledPOTtxt = ROOT.TLatex(.5,.95,"filledPOTtext1")
                filledPOTtxt.SetTextAlign(22)
                filledPOTtxt.SetTextSize(0.05)
                filledPOTtxt.DrawLatexNDC(.5,.95, "Normalised POT: " + "{0:e}".format(DataPOT))
            #Draw empty ratio plot
            canvas.cd()
        elif (plotType == "log"): #Draw ratio plots
            #Draw filled log plot
            print("Drawing log plots")
            p1.SetLogy()

            if (ymin == -1.0):
                tmpDataHist.SetMinimum(1)
                tmpMCHist.SetMinimum(1)

            plotter.DrawDataMC(tmpDataHist, tmpMCHist, 1.0, "TL", True)
            ymin = ROOT.gPad.GetFrame().GetY1()
            ymax = ROOT.gPad.GetFrame().GetY2()
            #std::cout << "Filled ymin: " << ymin << " ymax: " << ymax << std::endl
            PlotRegions(np.power(10,ymin), np.power(10,ymax), 0, 0)
            axistext = ROOT.TLatex(.5,.95,"axistext1")
            axistext.SetTextAlign(22)
            axistext.SetTextSize(0.05)
            axistext.DrawLatexNDC(.85,.05, "Vertex Z Position (mm)")
            if (DataPOT!=0 and MCPOT!=0):
                filledPOTtxt = ROOT.TLatex(.5,.95,"filledPOTtext1")
                filledPOTtxt.SetTextAlign(22)
                filledPOTtxt.SetTextSize(0.05)
                filledPOTtxt.DrawLatexNDC(.5,.95, "Normalised POT: " + "{0:e}".format(DataPOT))
            #Draw empty log plot
            canvas.cd()
        canvas.cd()
        return

def get_quantile_bins(h3, num_bins, axis = 2):
    """
    Calculates the Z-axis (var3) bin boundaries corresponding to equal-probability 
    quantiles (chunks) of the total integral.
    
    Args:
        h3 (ROOT.TH3): The 3D histogram to process.
        num_bins (int): The number of quantile bins to divide the distribution into 
                        (e.g., 10 for deciles, 4 for quartiles, 100 for percentiles).
                        
    Returns:
        list of int: Bin indices representing the boundaries 
        (including underflow and overflow).
    """
    if num_bins < 2:
        raise ValueError("num_bins must be at least 2 to create boundaries.")
        
    # 1. Project to Z using a unique name to avoid ROOT directory namespace collisions
    proj_name = f"h_var3_proj_temp"
    if axis == 0:
        h_var3 = h3.ProjectionX(proj_name)
    elif axis == 1:
        h_var3 = h3.ProjectionY(proj_name)
    elif axis == 2:
        h_var3 = h3.ProjectionZ(proj_name)
    
    # 2. Calculate the quantile thresholds dynamically
    num_boundaries = num_bins - 1
    
    # Generate probabilities: e.g. for num_bins=4, generates [0.25, 0.50, 0.75]
    probs = [float(i) / num_bins for i in range(1, num_bins)]
    
    probabilities = array('d', probs)
    quantiles = array('d', [0.0] * num_boundaries)
    
    h_var3.GetQuantiles(num_boundaries, quantiles, probabilities)
    
    # 3. Map exact quantile values to discrete Z-axis bin indices
    z_axis = h3.GetZaxis()
    bin_boundaries = [0] 
    
    for q in quantiles:
        bin_boundaries.append(z_axis.FindBin(q))
        
    bin_boundaries.append(z_axis.GetNbins() + 1)
    
    # 4. Clean up the temporary 1D projection to prevent memory leaks in ROOT
    h_var3.Delete()
    
    return bin_boundaries

def get_var3_decile_bins(h3):
    """
    Calculates the Z-axis (var3) bin boundaries corresponding to deciles 
    (10% chunks) of the total integral.
    
    Returns:
        list of int: 11 bin indices representing the 10 decile boundaries 
        (including underflow and overflow).
    """
    # 1. Project to Z to get the var3 weight distribution
    h_var3 = h3.ProjectionZ("h_var3_proj_temp")
    
    # 2. Calculate the 9 quantile thresholds
    probabilities = array('d', [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    quantiles = array('d', [0.0] * 9)
    h_var3.GetQuantiles(9, quantiles, probabilities)
    
    # 3. Map exact quantile values to discrete Z-axis bin indices
    z_axis = h3.GetZaxis()
    bin_boundaries = [0] 
    
    for q in quantiles:
        bin_boundaries.append(z_axis.FindBin(q))
        
    bin_boundaries.append(z_axis.GetNbins() + 1)
    
    # 4. Clean up the temporary 1D projection to prevent memory leaks in ROOT
    h_var3.Delete()
    
    return bin_boundaries

def extract_var1_decile_hists(h3, bin_boundaries):
    """
    Projects a TH3D onto the Y-axis (var2) across 10 discrete Z-axis (var3) 
    bin ranges, integrating over the X-axis (var1).
    
    Args:
        h3 (ROOT.TH3D): The source histogram.
        bin_boundaries (list): 11 bin boundaries generated by get_var3_decile_bins.
        
    Returns:
        list of ROOT.TH1D: 10 histograms of var2.
    """
    #if len(bin_boundaries) != 11:
    #    raise ValueError(f"Expected 11 bin boundaries, got {len(bin_boundaries)}.")
        
    decile_hists = []
    
    for i in range(len(bin_boundaries)-1):
        start_bin = bin_boundaries[i]
        
        # Prevent bin overlap, except for the upper bound of the final decile
        end_bin = bin_boundaries[i+1] - 1 if i < 9 else bin_boundaries[i+1]
        
        # Failsafe for highly peaked data falling into a single bin
        if start_bin > end_bin:
            end_bin = start_bin
            
        h_name = f"h_var2_decile_{i+1}"
        
        # X bins: 0 to -1 (flattens var1)
        # Z bins: start_bin to end_bin (slices var3)
        h_y = h3.ProjectionX(h_name, 0, -1, start_bin, end_bin)
        h_y.SetTitle(f"Var 2 (Var 3 Decile {i+1});Var 2;Events")
        
        decile_hists.append(h_y)
        
    return decile_hists

def extract_var2_decile_hists(h3, bin_boundaries):
    """
    Projects a TH3D onto the Y-axis (var2) across 10 discrete Z-axis (var3) 
    bin ranges, integrating over the X-axis (var1).
    
    Args:
        h3 (ROOT.TH3D): The source histogram.
        bin_boundaries (list): 11 bin boundaries generated by get_var3_decile_bins.
        
    Returns:
        list of ROOT.TH1D: 10 histograms of var2.
    """
    if len(bin_boundaries) != 11:
        raise ValueError(f"Expected 11 bin boundaries, got {len(bin_boundaries)}.")
        
    decile_hists = []
    
    for i in range(10):
        start_bin = bin_boundaries[i]
        
        # Prevent bin overlap, except for the upper bound of the final decile
        end_bin = bin_boundaries[i+1] - 1 if i < 9 else bin_boundaries[i+1]
        
        # Failsafe for highly peaked data falling into a single bin
        if start_bin > end_bin:
            end_bin = start_bin
            
        h_name = f"h_var2_decile_{i+1}"
        
        # X bins: 0 to -1 (flattens var1)
        # Z bins: start_bin to end_bin (slices var3)
        h_y = h3.ProjectionY(h_name, 0, -1, start_bin, end_bin)
        h_y.SetTitle(f"Var 2 (Var 3 Decile {i+1});Var 2;Events")
        
        decile_hists.append(h_y)
        
    return decile_hists

def PlotValidation(inpath, outbase, nu=True):
    print("outbase: ", outbase)
    global filledplaylistsused, emptyplaylistsused
    ROOT.gStyle.SetOptTitle(False)

    foundPaths = find_paths(inpath, "VertexValidations.root")
    
    #Initialising histograms to be plotted later

    #Water Filled Histograms
    ANNZ_Data_Filled = None
    ANNZ_MC_Filled = None
    ANNZUncalibrated_Data_Filled = None
    ANNZUncalibrated_MC_Filled = None
    ANNZ_vs_Conf_vs_Erecoil_MC_Filled = None
    ANNZ_vs_Conf_vs_Erecoil_Data_Filled = None

    #Water Empty Histograms
    ANNZ_Data_Empty = None
    ANNZ_MC_Empty = None
    ANNZUncalibrated_Data_Empty = None
    ANNZUncalibrated_MC_Empty = None
    ANNZ_vs_Conf_vs_Erecoil_MC_Empty = None
    ANNZ_vs_Conf_vs_Erecoil_Data_Empty = None

    ANNUncalibratedZResidual_vs_Conf_vs_Erecoil = None
    ANNZResidual_vs_Conf_vs_Erecoil = None
    ANNZResidual_vs_Conf_vs_Erecoil_Tracker = None
    ANNZResidual_vs_Conf_vs_Erecoil_Tgt2Iron = None
    ANNZResidual_vs_Conf_vs_Erecoil_Tgt2Lead = None
    ANNZResidual_vs_Conf_vs_Erecoil_Tgt3Iron = None
    ANNZResidual_vs_Conf_vs_Erecoil_Tgt3Lead = None
    ANNZResidual_vs_Conf_vs_Erecoil_Tgt3Carbon = None
    ANNZResidual_vs_Conf_vs_Erecoil_Tgt4Lead = None
    ANNZResidual_vs_Conf_vs_Erecoil_Tg5Iron = None
    ANNZResidual_vs_Conf_vs_Erecoil_Tg5Lead = None
    ANNZResidual_vs_Conf_vs_Erecoil_Water = None

    #Histograms that contain both filled + empty
    ANNConfidence_Data = None
    ANNConfidence_MC = None
    ANNUnalibratedConfidence_Data = None
    ANNUnalibratedConfidence_MC = None

    ANNCalibratedConfidence_Selected1 = None
    ANNCalibratedConfidence_SelectedSignal1 = None
    ANNCalibratedConfidence_Signal1 = None
    ANNCalibratedConfidence_Selected2 = None
    ANNCalibratedConfidence_SelectedSignal2 = None
    ANNCalibratedConfidence_Signal2 = None
    ANNCalibratedConfidence_Selected3 = None
    ANNCalibratedConfidence_SelectedSignal3 = None
    ANNCalibratedConfidence_Signal3 = None

    ANNConfidence_vs_NHitsMC = None
    ANNConfidence_vs_NHitsQE_MC = None
    ANNConfidence_vs_NHitsRES_MC = None
    ANNConfidence_vs_NHitsDIS_MC = None
    ANNConfidence_vs_NHitsCOHPI_MC = None
    ANNConfidence_vs_NHitsAMNUGAMMA_MC = None
    ANNConfidence_vs_NHitsIMD_MC = None
    ANNConfidence_vs_NHitsQENUEEL_MC = None
    ANNConfidence_vs_NHitsQEMEC_MC = None
    ANNConfidence_vs_NHitsQEUnknownInt_MC = None
    ANNConfidence_vs_NHitsData = None
    ANNConfidence_vs_NTracks_vs_NHitsMC = None
    ANNConfidence_vs_NTracks_vs_NHitsData = None

    ANNZ_Data_Combined = None
    ANNZ_MC_Combined = None
    ANNZUncalibrated_Data_Combined = None
    ANNZUncalibrated_MC_Combined = None

    ANNZ_vs_Conf_vs_Erecoil_MC = None
    ANNZ_vs_Conf_vs_Erecoil_Data = None

    filledMCPOT = 0.0
    emptyMCPOT = 0.0
    filledDataPOT = 0.0
    emptyDataPOT = 0.0
    totalMCPOT = 0.0
    totalDataPOT = 0.0

    for path in foundPaths:
        print("path: ", path)
        validationFile = ROOT.TFile.Open(path)
        playlistUsed = validationFile.Get("PlaylistUsed").GetTitle()

        print("filled_or_empty_me_playlist(playlistUsed): ", filled_or_empty_me_playlist(playlistUsed))
        print("nu_or_anti_nu_mode(playlistUsed): ", nu_or_anti_nu_mode(playlistUsed))
        if nu_or_anti_nu_mode(playlistUsed) == 1 and not nu: #neutrino mode
            print("Skipping1")
            continue
        elif nu_or_anti_nu_mode(playlistUsed) == 2 and nu: #antineutrino mode
            print("Skipping2")
            continue
        print("Here1")
        if filled_or_empty_me_playlist(playlistUsed) == 1: #filled mode
            print("In filled file")
            ANNZ_Data_Filled = add_hist(ANNZ_Data_Filled, validationFile, "ANNZ_Data")
            ANNZ_MC_Filled = add_hist(ANNZ_MC_Filled, validationFile, "ANNZ_MC")
            ANNZUncalibrated_Data_Filled = add_hist(ANNZUncalibrated_Data_Filled, validationFile, "ANNZUncalibrated_Data")
            ANNZUncalibrated_MC_Filled = add_hist(ANNZUncalibrated_MC_Filled, validationFile, "ANNZUncalibrated_MC")
            ANNZ_vs_Conf_vs_Erecoil_MC_Filled=add_hist(ANNZ_vs_Conf_vs_Erecoil_MC_Filled, validationFile, "ANNZ_vs_Conf_vs_Erecoil_MC")
            ANNZ_vs_Conf_vs_Erecoil_Data_Filled=add_hist(ANNZ_vs_Conf_vs_Erecoil_Data_Filled, validationFile, "ANNZ_vs_Conf_vs_Erecoil_Data")
            filledMCPOT += validationFile.Get("MCPOT").GetVal()
            print("testVal: ", validationFile.Get("MCPOT").GetVal())
            filledDataPOT += validationFile.Get("DataPOT").GetVal()


        elif filled_or_empty_me_playlist(playlistUsed) == 2: #empty mode
            print("In empty file")
            ANNZ_Data_Empty = add_hist(ANNZ_Data_Empty, validationFile, "ANNZ_Data")
            ANNZ_MC_Empty = add_hist(ANNZ_MC_Empty, validationFile, "ANNZ_MC")
            ANNZUncalibrated_Data_Empty = add_hist(ANNZUncalibrated_Data_Empty, validationFile, "ANNZUncalibrated_Data")
            ANNZUncalibrated_MC_Empty = add_hist(ANNZUncalibrated_MC_Empty, validationFile, "ANNZUncalibrated_MC")
            ANNZ_vs_Conf_vs_Erecoil_MC_Empty=add_hist(ANNZ_vs_Conf_vs_Erecoil_MC_Empty, validationFile, "ANNZ_vs_Conf_vs_Erecoil_MC")
            ANNZ_vs_Conf_vs_Erecoil_Data_Empty=add_hist(ANNZ_vs_Conf_vs_Erecoil_Data_Empty, validationFile, "ANNZ_vs_Conf_vs_Erecoil_Data")
            emptyMCPOT += validationFile.Get("MCPOT").GetVal()
            emptyDataPOT += validationFile.Get("DataPOT").GetVal()

        ANNZ_Data_Combined=add_hist(ANNZ_Data_Combined, validationFile, "ANNZ_Data")
        ANNZ_MC_Combined=add_hist(ANNZ_MC_Combined, validationFile, "ANNZ_MC")
        ANNZUncalibrated_Data_Combined=add_hist(ANNZUncalibrated_Data_Combined, validationFile, "ANNZUncalibrated_Data")
        ANNZUncalibrated_MC_Combined=add_hist(ANNZUncalibrated_MC_Combined, validationFile, "ANNZUncalibrated_MC")

        ANNZ_vs_Conf_vs_Erecoil_MC=add_hist(ANNZ_vs_Conf_vs_Erecoil_MC, validationFile, "ANNZ_vs_Conf_vs_Erecoil_MC")
        ANNZ_vs_Conf_vs_Erecoil_Data=add_hist(ANNZ_vs_Conf_vs_Erecoil_Data, validationFile, "ANNZ_vs_Conf_vs_Erecoil_Data")

        ANNConfidence_Data = add_hist(ANNConfidence_Data, validationFile, "ANNConfidence_Data")
        ANNConfidence_MC = add_hist(ANNConfidence_MC, validationFile, "ANNConfidence_MC")
        ANNUnalibratedConfidence_Data = add_hist(ANNUnalibratedConfidence_Data, validationFile, "ANNUnalibratedConfidence_Data")
        ANNUnalibratedConfidence_MC = add_hist(ANNUnalibratedConfidence_MC, validationFile, "ANNUnalibratedConfidence_MC")
        
        ANNCalibratedConfidence_Selected1 = add_hist(ANNCalibratedConfidence_Selected1, validationFile, "ANNCalibratedConfidence_Selected1")
        ANNCalibratedConfidence_SelectedSignal1 = add_hist(ANNCalibratedConfidence_SelectedSignal1, validationFile, "ANNCalibratedConfidence_SelectedSignal1")
        ANNCalibratedConfidence_Signal1 = add_hist(ANNCalibratedConfidence_Signal1, validationFile, "ANNCalibratedConfidence_Signal1")
        ANNCalibratedConfidence_Selected2 = add_hist(ANNCalibratedConfidence_Selected2, validationFile, "ANNCalibratedConfidence_Selected2")
        ANNCalibratedConfidence_SelectedSignal2 = add_hist(ANNCalibratedConfidence_SelectedSignal2, validationFile, "ANNCalibratedConfidence_SelectedSignal2")
        ANNCalibratedConfidence_Signal2 = add_hist(ANNCalibratedConfidence_Signal2, validationFile, "ANNCalibratedConfidence_Signal2")
        ANNCalibratedConfidence_Selected3 = add_hist(ANNCalibratedConfidence_Selected3, validationFile, "ANNCalibratedConfidence_Selected3")
        ANNCalibratedConfidence_SelectedSignal3 = add_hist(ANNCalibratedConfidence_SelectedSignal3, validationFile, "ANNCalibratedConfidence_SelectedSignal3")
        ANNCalibratedConfidence_Signal3 = add_hist(ANNCalibratedConfidence_Signal3, validationFile, "ANNCalibratedConfidence_Signal3")


        ANNConfidence_vs_NHitsMC = add_hist(ANNConfidence_vs_NHitsMC, validationFile, "ANNConfidence_vs_NHitsMC")
        ANNConfidence_vs_NHitsQE_MC = add_hist(ANNConfidence_vs_NHitsQE_MC, validationFile, "ANNConfidence_vs_NHitsQE_MC")
        ANNConfidence_vs_NHitsRES_MC = add_hist(ANNConfidence_vs_NHitsRES_MC, validationFile, "ANNConfidence_vs_NHitsRES_MC")
        ANNConfidence_vs_NHitsDIS_MC = add_hist(ANNConfidence_vs_NHitsDIS_MC, validationFile, "ANNConfidence_vs_NHitsDIS_MC")
        ANNConfidence_vs_NHitsCOHPI_MC = add_hist(ANNConfidence_vs_NHitsCOHPI_MC, validationFile, "ANNConfidence_vs_NHitsCOHPI_MC")
        ANNConfidence_vs_NHitsAMNUGAMMA_MC = add_hist(ANNConfidence_vs_NHitsAMNUGAMMA_MC, validationFile, "ANNConfidence_vs_NHitsAMNUGAMMA_MC")
        ANNConfidence_vs_NHitsIMD_MC = add_hist(ANNConfidence_vs_NHitsIMD_MC, validationFile, "ANNConfidence_vs_NHitsIMD_MC")
        ANNConfidence_vs_NHitsQENUEEL_MC = add_hist(ANNConfidence_vs_NHitsQENUEEL_MC, validationFile, "ANNConfidence_vs_NHitsQENUEEL_MC")
        ANNConfidence_vs_NHitsQEMEC_MC = add_hist(ANNConfidence_vs_NHitsQEMEC_MC, validationFile, "ANNConfidence_vs_NHitsQEMEC_MC")
        ANNConfidence_vs_NHitsQEUnknownInt_MC = add_hist(ANNConfidence_vs_NHitsQEUnknownInt_MC, validationFile, "ANNConfidence_vs_NHitsQEUnknownInt_MC")
        ANNConfidence_vs_NHitsData = add_hist(ANNConfidence_vs_NHitsData, validationFile, "ANNConfidence_vs_NHitsData")
        ANNConfidence_vs_NTracks_vs_NHitsMC = add_hist(ANNConfidence_vs_NTracks_vs_NHitsMC, validationFile, "ANNConfidence_vs_NTracks_vs_NHitsMC")
        ANNConfidence_vs_NTracks_vs_NHitsData = add_hist(ANNConfidence_vs_NTracks_vs_NHitsData, validationFile, "ANNConfidence_vs_NTracks_vs_NHitsData")

        ANNUncalibratedZResidual_vs_Conf_vs_Erecoil = add_hist(ANNUncalibratedZResidual_vs_Conf_vs_Erecoil, validationFile, "ANNUncalibratedZResidual_vs_Conf_vs_Erecoil")
        ANNZResidual_vs_Conf_vs_Erecoil = add_hist(ANNZResidual_vs_Conf_vs_Erecoil, validationFile, "ANNZResidual_vs_Conf_vs_Erecoil")

        totalMCPOT += validationFile.Get("MCPOT").GetVal()
        totalDataPOT += validationFile.Get("DataPOT").GetVal()
        validationFile.Close()
    # To avoid histograms being destroyed when files are closed,
    # we'll store them in a dictionary and SetDirectory(0) on each.

    print(f"empty pot: {emptyDataPOT}/{emptyMCPOT}")
    print(f"filled pot: {filledDataPOT}/{filledMCPOT}")


    canvas = ROOT.TCanvas()
    Proj_ANNUncalibratedZResidual_vs_Conf_vs_Erecoil = ANNUncalibratedZResidual_vs_Conf_vs_Erecoil.ProjectionX("Uncalibrated", 0, -1)
    Proj_ANNZResidual_vs_Conf_vs_Erecoil = ANNZResidual_vs_Conf_vs_Erecoil.ProjectionX("Calibrated", 0, -1)
    Proj_ANNUncalibratedZResidual_vs_Conf_vs_Erecoil.SetLineColor(2)
    Proj_ANNUncalibratedZResidual_vs_Conf_vs_Erecoil.SetStats(False)
    Proj_ANNUncalibratedZResidual_vs_Conf_vs_Erecoil.Draw("HIST")
    Proj_ANNZResidual_vs_Conf_vs_Erecoil.SetLineColor(4)
    Proj_ANNZResidual_vs_Conf_vs_Erecoil.Draw("HIST SAME")
    legend = ROOT.TLegend(0.55,0.7,0.9,0.95)
    legend.SetBorderSize(0)
    #legend.SetNColumns(4)
    legend.AddEntry(Proj_ANNZResidual_vs_Conf_vs_Erecoil,"Calibrated","l")
    legend.AddEntry(Proj_ANNUncalibratedZResidual_vs_Conf_vs_Erecoil,"Uncalibrated","l")
    legend.Draw()
    filename = f"{outbase}/ANN_Residuals.pdf"
    canvas.SaveAs(filename)

    canvas.Clear()
    ANNConfidence_Data.SetLineColor(2)
    ANNConfidence_MC.SetLineColor(3)
    ANNUnalibratedConfidence_Data.SetLineColor(4)
    ANNUnalibratedConfidence_MC.SetLineColor(5)

    if ANNConfidence_Data.Integral() > 0:
        ANNConfidence_Data.Scale(1.0 / ANNConfidence_Data.Integral())
    if ANNConfidence_MC.Integral() > 0:
        ANNConfidence_MC.Scale(1.0 / ANNConfidence_MC.Integral())
    if ANNUnalibratedConfidence_Data.Integral() > 0:
        ANNUnalibratedConfidence_Data.Scale(1.0 / ANNUnalibratedConfidence_Data.Integral())
    if ANNUnalibratedConfidence_MC.Integral() > 0:
        ANNUnalibratedConfidence_MC.Scale(1.0 / ANNUnalibratedConfidence_MC.Integral())
    ANNConfidence_Data.SetStats(False)
    ANNConfidence_Data.Draw("HIST")
    ANNConfidence_MC.Draw("HIST SAME")
    ANNUnalibratedConfidence_Data.Draw("HIST SAME")
    ANNUnalibratedConfidence_MC.Draw("HIST SAME")
    legend = ROOT.TLegend(0.15,0.4,0.5,0.95)
    legend.SetBorderSize(0)
    #legend.SetNColumns(4)
    legend.AddEntry(ANNConfidence_Data,"Calibrated Data","l")
    legend.AddEntry(ANNConfidence_MC,"Calibrated MC","l")
    legend.AddEntry(ANNUnalibratedConfidence_Data,"Uncalibrated Data","l")
    legend.AddEntry(ANNUnalibratedConfidence_MC,"Uncalibrated MC","l")
    legend.Draw()
    filename = f"{outbase}/ANN_Confidences.pdf"
    canvas.SaveAs(filename)

    canvas.Clear()
    CumANNConfidence_Data= ANNConfidence_Data.Clone()
    CumANNConfidence_MC= ANNConfidence_MC.Clone()
    CumANNUnalibratedConfidence_Data= ANNUnalibratedConfidence_Data.Clone()
    CumANNUnalibratedConfidence_MC= ANNUnalibratedConfidence_MC.Clone()
    CumANNConfidence_Data = ANNConfidence_Data.GetCumulative(True)
    CumANNConfidence_MC = ANNConfidence_MC.GetCumulative(True)
    CumANNUnalibratedConfidence_Data = ANNUnalibratedConfidence_Data.GetCumulative(True)
    CumANNUnalibratedConfidence_MC = ANNUnalibratedConfidence_MC.GetCumulative(True)
    if CumANNConfidence_Data.GetMaximum() > 0:
        CumANNConfidence_Data.Scale(1.0 / CumANNConfidence_Data.GetMaximum())
    if CumANNConfidence_MC.GetMaximum() > 0:
        CumANNConfidence_MC.Scale(1.0 / CumANNConfidence_MC.GetMaximum())
    if CumANNUnalibratedConfidence_Data.GetMaximum() > 0:
        CumANNUnalibratedConfidence_Data.Scale(1.0 / CumANNUnalibratedConfidence_Data.GetMaximum())
    if CumANNUnalibratedConfidence_MC.GetMaximum() > 0:
        CumANNUnalibratedConfidence_MC.Scale(1.0 / CumANNUnalibratedConfidence_MC.GetMaximum())
    CumANNConfidence_Data.SetLineColor(2)
    CumANNConfidence_MC.SetLineColor(3)
    CumANNUnalibratedConfidence_Data.SetLineColor(4)
    CumANNUnalibratedConfidence_MC.SetLineColor(5)
    CumANNConfidence_Data.SetStats(False)
    CumANNConfidence_Data.Draw("HIST")
    CumANNConfidence_MC.Draw("HIST SAME")
    CumANNUnalibratedConfidence_Data.Draw("HIST SAME")
    CumANNUnalibratedConfidence_MC.Draw("HIST SAME")
    legend = ROOT.TLegend(0.55,0.7,0.9,0.95)
    legend.SetBorderSize(0)
    #legend.SetNColumns(4)
    legend.AddEntry(CumANNConfidence_Data,"Uncalibrated Data","l")
    legend.AddEntry(CumANNConfidence_MC,"Uncalibrated MC","l")
    legend.AddEntry(CumANNUnalibratedConfidence_Data,"Calibrated Data","l")
    legend.AddEntry(CumANNUnalibratedConfidence_MC,"Calibrated MC","l")
    legend.Draw()
    filename = f"{outbase}/ANN_CDFConfidences.pdf"
    canvas.SaveAs(filename)


    canvas.Clear()
    CumANNConfidence_Data= ANNConfidence_Data.Clone()
    CumANNConfidence_MC= ANNConfidence_MC.Clone()
    CumANNUnalibratedConfidence_Data= ANNUnalibratedConfidence_Data.Clone()
    CumANNUnalibratedConfidence_MC= ANNUnalibratedConfidence_MC.Clone()
    CumANNConfidence_Data = ANNConfidence_Data.GetCumulative(False)
    CumANNConfidence_MC = ANNConfidence_MC.GetCumulative(False)
    CumANNUnalibratedConfidence_Data = ANNUnalibratedConfidence_Data.GetCumulative(False)
    CumANNUnalibratedConfidence_MC = ANNUnalibratedConfidence_MC.GetCumulative(False)
    if CumANNConfidence_Data.GetMaximum() > 0:
        CumANNConfidence_Data.Scale(1.0 / CumANNConfidence_Data.GetMaximum())
    if CumANNConfidence_MC.GetMaximum() > 0:
        CumANNConfidence_MC.Scale(1.0 / CumANNConfidence_MC.GetMaximum())
    if CumANNUnalibratedConfidence_Data.GetMaximum() > 0:
        CumANNUnalibratedConfidence_Data.Scale(1.0 / CumANNUnalibratedConfidence_Data.GetMaximum())
    if CumANNUnalibratedConfidence_MC.GetMaximum() > 0:
        CumANNUnalibratedConfidence_MC.Scale(1.0 / CumANNUnalibratedConfidence_MC.GetMaximum())
    CumANNConfidence_Data.SetLineColor(2)
    CumANNConfidence_MC.SetLineColor(3)
    CumANNUnalibratedConfidence_Data.SetLineColor(4)
    CumANNUnalibratedConfidence_MC.SetLineColor(5)
    CumANNConfidence_Data.SetStats(False)
    CumANNConfidence_Data.Draw("HIST")
    CumANNConfidence_MC.Draw("HIST SAME")
    CumANNUnalibratedConfidence_Data.Draw("HIST SAME")
    CumANNUnalibratedConfidence_MC.Draw("HIST SAME")
    legend = ROOT.TLegend(0.55,0.7,0.9,0.95)
    legend.SetBorderSize(0)
    #legend.SetNColumns(4)
    legend.AddEntry(CumANNConfidence_Data,"Uncalibrated Data","l")
    legend.AddEntry(CumANNConfidence_MC,"Uncalibrated MC","l")
    legend.AddEntry(CumANNUnalibratedConfidence_Data,"Calibrated Data","l")
    legend.AddEntry(CumANNUnalibratedConfidence_MC,"Calibrated MC","l")
    legend.Draw()
    filename = f"{outbase}/ANN_CDFRevConfidences.pdf"
    canvas.SaveAs(filename)


    selected1 = ANNCalibratedConfidence_Selected1.GetCumulative(False)
    selected2 = ANNCalibratedConfidence_Selected2.GetCumulative(False)
    selected3 = ANNCalibratedConfidence_Selected3.GetCumulative(False)

    selectedSignal1 = ANNCalibratedConfidence_SelectedSignal1.GetCumulative(False)
    selectedSignal2 = ANNCalibratedConfidence_SelectedSignal2.GetCumulative(False)
    selectedSignal3 = ANNCalibratedConfidence_SelectedSignal3.GetCumulative(False)

    signal = ANNCalibratedConfidence_Signal1.GetSumOfWeights()

    purity1 = selectedSignal1.Clone("Purity")
    purity2 = selectedSignal2.Clone("Purity")
    purity3 = selectedSignal3.Clone("Purity")

    recall1 = selectedSignal1.Clone("Recall")
    recall2 = selectedSignal2.Clone("Recall")
    recall3 = selectedSignal3.Clone("Recall")

    F1Score1 = selectedSignal1.Clone("F1 Score")
    F1Score2 = selectedSignal2.Clone("F1 Score")
    F1Score3 = selectedSignal3.Clone("F1 Score")

    purity1.Divide(selected1)
    purity2.Divide(selected2)
    purity3.Divide(selected3)

    recall1.Scale(1/signal)
    recall2.Scale(1/signal)
    recall3.Scale(1/signal)

    PrecisionRecallSum1 = recall1.Clone("Precision + Recall")
    PrecisionRecallSum2 = recall2.Clone("Precision + Recall")
    PrecisionRecallSum3 = recall3.Clone("Precision + Recall")

    PrecisionRecallSum1.Add(purity1)
    PrecisionRecallSum2.Add(purity2)
    PrecisionRecallSum3.Add(purity3)

    F1Score1.Multiply(purity1, recall1, 2)
    F1Score2.Multiply(purity2, recall2, 2)
    F1Score3.Multiply(purity3, recall3, 2)

    F1Score1.Divide(PrecisionRecallSum1)
    F1Score2.Divide(PrecisionRecallSum2)
    F1Score3.Divide(PrecisionRecallSum3)


    #Plotting test
    canvas.Clear()
    ANNCalibratedConfidence_SelectedSignal1.SetTitle("Test")
    ANNCalibratedConfidence_SelectedSignal1.SetYTitle("Number of events")
    ANNCalibratedConfidence_SelectedSignal1.SetXTitle("Calibrated ANN Confidence")
    ANNCalibratedConfidence_SelectedSignal1.SetLineWidth(3)
    ANNCalibratedConfidence_SelectedSignal1.SetLineColor(ROOT.kBlue)
    ANNCalibratedConfidence_SelectedSignal1.SetStats(0)
    ANNCalibratedConfidence_SelectedSignal1.Draw("HIST")

    ANNCalibratedConfidence_SelectedSignal2.SetLineWidth(3)
    ANNCalibratedConfidence_SelectedSignal2.SetLineColor(ROOT.kRed)
    ANNCalibratedConfidence_SelectedSignal2.Draw("HIST SAME")

    ANNCalibratedConfidence_SelectedSignal3.SetLineWidth(3)
    ANNCalibratedConfidence_SelectedSignal3.SetLineColor(ROOT.kGreen)
    ANNCalibratedConfidence_SelectedSignal3.Draw("HIST SAME")


    legend = ROOT.TLegend(0.8,0.7,0.9,0.8)
    legend.SetBorderSize(0)
    legend.AddEntry(ANNCalibratedConfidence_SelectedSignal1,"Top Plane","l")
    legend.AddEntry(ANNCalibratedConfidence_SelectedSignal2,"2nd Plane","l")
    legend.AddEntry(ANNCalibratedConfidence_SelectedSignal3,"3rd Plane","l")
    legend.Draw()
    canvas.SaveAs("test.pdf")

    #Plotting selected signal
    canvas.Clear()
    selectedSignal1.SetTitle("Selected Signal")
    selectedSignal1.SetYTitle("Number of events")
    selectedSignal1.SetXTitle("Calibrated ANN Confidence")
    selectedSignal1.SetLineWidth(3)
    selectedSignal1.SetLineColor(ROOT.kBlue)
    selectedSignal1.SetStats(0)
    selectedSignal1.Draw("HIST")

    selectedSignal2.SetLineWidth(3)
    selectedSignal2.SetLineColor(ROOT.kRed)
    selectedSignal2.Draw("HIST SAME")

    selectedSignal3.SetLineWidth(3)
    selectedSignal3.SetLineColor(ROOT.kGreen)
    selectedSignal3.Draw("HIST SAME")


    legend = ROOT.TLegend(0.8,0.7,0.9,0.8)
    legend.SetBorderSize(0)
    legend.AddEntry(selectedSignal1,"Top Plane","l")
    legend.AddEntry(selectedSignal2,"2nd Plane","l")
    legend.AddEntry(selectedSignal3,"3rd Plane","l")
    legend.Draw()
    canvas.SaveAs("SelectedSignal.pdf")

    #Plotting signal
    canvas.Clear()
    selected1.SetTitle("Selected")
    selected1.SetYTitle("Number of events")
    selected1.SetXTitle("Calibrated ANN Confidence")
    selected1.SetLineWidth(3)
    selected1.SetLineColor(ROOT.kBlue)
    selected1.SetStats(0)
    selected1.Draw("HIST")

    selected2.SetLineWidth(3)
    selected2.SetLineColor(ROOT.kRed)
    selected2.Draw("HIST SAME")

    selected3.SetLineWidth(3)
    selected3.SetLineColor(ROOT.kGreen)
    selected3.Draw("HIST SAME")


    legend = ROOT.TLegend(0.8,0.7,0.9,0.8)
    legend.SetBorderSize(0)
    legend.AddEntry(selected1,"Top Plane","l")
    legend.AddEntry(selected2,"2nd Plane","l")
    legend.AddEntry(selected3,"3rd Plane","l")
    legend.Draw()
    canvas.SaveAs("Selected.pdf")

    #Plotting purity
    canvas.Clear()
    purity1.SetTitle("Purity")
    purity1.SetYTitle("Number of events")
    purity1.SetXTitle("Calibrated ANN Confidence")
    purity1.SetLineWidth(3)
    purity1.SetLineColor(ROOT.kBlue)
    purity1.SetStats(0)
    purity1.Draw("HIST")
    purity1.SetMinimum(0.95)
    purity1.SetMaximum(1)

    purity2.SetLineWidth(3)
    purity2.SetLineColor(ROOT.kRed)
    purity2.Draw("HIST SAME")

    purity3.SetLineWidth(3)
    purity3.SetLineColor(ROOT.kGreen)
    purity3.Draw("HIST SAME")


    legend = ROOT.TLegend(0.8,0.7,0.9,0.8)
    legend.SetBorderSize(0)
    legend.AddEntry(purity1,"Top Plane","l")
    legend.AddEntry(purity2,"2nd Plane","l")
    legend.AddEntry(purity3,"3rd Plane","l")
    legend.Draw()
    canvas.SaveAs("Purity.pdf")

    #Plotting recall
    canvas.Clear()
    recall1.SetTitle("Recall")
    recall1.SetYTitle("Number of events")
    recall1.SetXTitle("Calibrated ANN Confidence")
    recall1.SetLineWidth(3)
    recall1.SetLineColor(ROOT.kBlue)
    recall1.SetStats(0)
    recall1.Draw("HIST")

    recall2.SetLineWidth(3)
    recall2.SetLineColor(ROOT.kRed)
    recall2.Draw("HIST SAME")

    recall3.SetLineWidth(3)
    recall3.SetLineColor(ROOT.kGreen)
    recall3.Draw("HIST SAME")

    legend = ROOT.TLegend(0.8,0.7,0.9,0.8)
    legend.SetBorderSize(0)
    legend.AddEntry(recall1,"Top Plane","l")
    legend.AddEntry(recall2,"2nd Plane","l")
    legend.AddEntry(recall3,"3rd Plane","l")
    legend.Draw()
    canvas.SaveAs("Recall.pdf")

    #Plotting F1 Score
    canvas.Clear()
    F1Score1.SetTitle("F1Score")
    F1Score1.SetYTitle("Number of events")
    F1Score1.SetXTitle("Calibrated ANN Confidence")
    F1Score1.SetLineWidth(3)
    F1Score1.SetLineColor(ROOT.kBlue)
    F1Score1.SetStats(0)
    F1Score1.Draw("HIST")

    F1Score2.SetLineWidth(3)
    F1Score2.SetLineColor(ROOT.kRed)
    F1Score2.Draw("HIST SAME")

    F1Score3.SetLineWidth(3)
    F1Score3.SetLineColor(ROOT.kGreen)
    F1Score3.Draw("HIST SAME")

    legend = ROOT.TLegend(0.8,0.7,0.9,0.8)
    legend.SetBorderSize(0)
    legend.AddEntry(F1Score1,"Top Plane","l")
    legend.AddEntry(F1Score2,"2nd Plane","l")
    legend.AddEntry(F1Score3,"3rd Plane","l")
    legend.Draw()
    canvas.SaveAs("F1Score.pdf")


    print(f"Num entries {ANNConfidence_vs_NTracks_vs_NHitsMC.GetEntries()}")

    ROOT.gStyle.SetOptTitle(True)
    NTracks_vs_NHitsMC = ANNConfidence_vs_NTracks_vs_NHitsMC.Project3D("yz")
    NTracks_vs_NHitsMC.SetStats(False)
    NTracks_vs_NHitsMC = normalize_by_row(NTracks_vs_NHitsMC)
    NTracks_vs_NHitsMC.SetTitle("Row normalised multiplicity vs number of hits (MC)")
    NTracks_vs_NHitsMC.GetXaxis().SetTitle("Number of hits")
    NTracks_vs_NHitsMC.GetYaxis().SetTitle("Multiplicity")
    NTracks_vs_NHitsMC.SetMinimum(0.0)
    NTracks_vs_NHitsMC.SetMaximum(1.0)
    NTracks_vs_NHitsMC.Draw("COLZ")
    canvas.SaveAs(f"{outbase}/NTracks_vs_NHitsMC.pdf")
    NTracks_vs_NHitsData = ANNConfidence_vs_NTracks_vs_NHitsData.Project3D("yz")
    NTracks_vs_NHitsData.SetStats(False)
    NTracks_vs_NHitsData = normalize_by_row(NTracks_vs_NHitsData)
    NTracks_vs_NHitsData.SetTitle("Row normalised multiplicity vs number of hits (Data)")
    NTracks_vs_NHitsData.GetXaxis().SetTitle("Number of hits")
    NTracks_vs_NHitsData.GetYaxis().SetTitle("Multiplicity")
    NTracks_vs_NHitsData.SetMinimum(0.0)
    NTracks_vs_NHitsData.SetMaximum(1.0)
    NTracks_vs_NHitsData.Draw("COLZ")
    canvas.SaveAs(f"{outbase}/NTracks_vs_NHitsData.pdf")

    n_bins_y = ANNConfidence_vs_NHitsMC.GetNbinsY()
    y_axis = ANNConfidence_vs_NHitsMC.GetYaxis()
    ROOT.gErrorIgnoreLevel = ROOT.kWarning #Otherwise you'll get a bunch of "INFO" messages on each file save

    #================================================================================================================
    #===============================  Producing plots broken down by number of hits  #===============================
    #================================================================================================================
    for iy in range(1, n_bins_y + 1):
        # Extract the bin edges
        y_lower = y_axis.GetBinLowEdge(iy)
        y_higher = y_axis.GetBinUpEdge(iy)
        
        plotter = PlotUtils.MnvPlotter()

        Proj_ANNConfidence_vs_NHitsQE_MC = ANNConfidence_vs_NHitsQE_MC.ProjectionX("QE", iy, iy)
        Proj_ANNConfidence_vs_NHitsRES_MC = ANNConfidence_vs_NHitsRES_MC.ProjectionX("RES", iy, iy)
        Proj_ANNConfidence_vs_NHitsDIS_MC = ANNConfidence_vs_NHitsDIS_MC.ProjectionX("DIS", iy, iy)
        Proj_ANNConfidence_vs_NHitsCOHPI_MC = ANNConfidence_vs_NHitsCOHPI_MC.ProjectionX("COHPI", iy, iy)
        Proj_ANNConfidence_vs_NHitsAMNUGAMMA_MC = ANNConfidence_vs_NHitsAMNUGAMMA_MC.ProjectionX("AMNUGAMMA", iy, iy)
        Proj_ANNConfidence_vs_NHitsIMD_MC = ANNConfidence_vs_NHitsIMD_MC.ProjectionX("IMD", iy, iy)
        Proj_ANNConfidence_vs_NHitsQENUEEL_MC = ANNConfidence_vs_NHitsQENUEEL_MC.ProjectionX("QENUEEL", iy, iy)
        Proj_ANNConfidence_vs_NHitsQEMEC_MC = ANNConfidence_vs_NHitsQEMEC_MC.ProjectionX("QEMEC", iy, iy)

        Proj_ANNConfidence_vs_NHitsMC = ANNConfidence_vs_NHitsMC.ProjectionX("MC", iy, iy)
        Proj_ANNConfidence_vs_NHitsData = ANNConfidence_vs_NHitsData.ProjectionX("Data", iy, iy)

        # Set the custom title. 
        # The :g format specifier removes unnecessary trailing zeros (e.g., 5.000 -> 5)
        #title = f"{y_lower:g} < n_hits < {y_higher:g}"
        #h1_proj.SetTitle(title)
        
        # Optional: Set axis titles if needed
        # h1_proj.GetXaxis().SetTitle("X Axis Label")
        # h1_proj.GetYaxis().SetTitle("Counts")

        Cumulative_Proj_ANNConfidence_vs_NHitsQE_MC = Proj_ANNConfidence_vs_NHitsQE_MC.GetCumulative(False)
        Cumulative_Proj_ANNConfidence_vs_NHitsRES_MC = Proj_ANNConfidence_vs_NHitsRES_MC.GetCumulative(False)
        Cumulative_Proj_ANNConfidence_vs_NHitsDIS_MC = Proj_ANNConfidence_vs_NHitsDIS_MC.GetCumulative(False)
        Cumulative_Proj_ANNConfidence_vs_NHitsCOHPI_MC = Proj_ANNConfidence_vs_NHitsCOHPI_MC.GetCumulative(False)
        Cumulative_Proj_ANNConfidence_vs_NHitsAMNUGAMMA_MC = Proj_ANNConfidence_vs_NHitsAMNUGAMMA_MC.GetCumulative(False)
        Cumulative_Proj_ANNConfidence_vs_NHitsIMD_MC = Proj_ANNConfidence_vs_NHitsIMD_MC.GetCumulative(False)
        Cumulative_Proj_ANNConfidence_vs_NHitsQENUEEL_MC = Proj_ANNConfidence_vs_NHitsQENUEEL_MC.GetCumulative(False)
        Cumulative_Proj_ANNConfidence_vs_NHitsQEMEC_MC = Proj_ANNConfidence_vs_NHitsQEMEC_MC.GetCumulative(False)

        Cumulative_Proj_ANNConfidence_vs_NHitsMC = Proj_ANNConfidence_vs_NHitsMC.GetCumulative(False)
        Cumulative_Proj_ANNConfidence_vs_NHitsData = Proj_ANNConfidence_vs_NHitsData.GetCumulative(False)
        
        # Normalise
        max_MC_val = Cumulative_Proj_ANNConfidence_vs_NHitsMC.GetMaximum()
        max_Data_val = Cumulative_Proj_ANNConfidence_vs_NHitsData.GetMaximum()
        if max_Data_val > 0:
            Cumulative_Proj_ANNConfidence_vs_NHitsData.Scale(1.0/max_Data_val)

        Cumulative_Proj_ANNConfidence_vs_NHitsQE_MC = PlotUtils.MnvH1D(Cumulative_Proj_ANNConfidence_vs_NHitsQE_MC)
        Cumulative_Proj_ANNConfidence_vs_NHitsRES_MC = PlotUtils.MnvH1D(Cumulative_Proj_ANNConfidence_vs_NHitsRES_MC)
        Cumulative_Proj_ANNConfidence_vs_NHitsDIS_MC = PlotUtils.MnvH1D(Cumulative_Proj_ANNConfidence_vs_NHitsDIS_MC)
        Cumulative_Proj_ANNConfidence_vs_NHitsCOHPI_MC = PlotUtils.MnvH1D(Cumulative_Proj_ANNConfidence_vs_NHitsCOHPI_MC)
        Cumulative_Proj_ANNConfidence_vs_NHitsAMNUGAMMA_MC = PlotUtils.MnvH1D(Cumulative_Proj_ANNConfidence_vs_NHitsAMNUGAMMA_MC)
        Cumulative_Proj_ANNConfidence_vs_NHitsIMD_MC = PlotUtils.MnvH1D(Cumulative_Proj_ANNConfidence_vs_NHitsIMD_MC)
        Cumulative_Proj_ANNConfidence_vs_NHitsQENUEEL_MC = PlotUtils.MnvH1D(Cumulative_Proj_ANNConfidence_vs_NHitsQENUEEL_MC)
        Cumulative_Proj_ANNConfidence_vs_NHitsQEMEC_MC = PlotUtils.MnvH1D(Cumulative_Proj_ANNConfidence_vs_NHitsQEMEC_MC)
        Cumulative_Proj_ANNConfidence_vs_NHitsData = PlotUtils.MnvH1D(Cumulative_Proj_ANNConfidence_vs_NHitsData)

        Proj_ANNConfidence_vs_NHitsQE_MC = PlotUtils.MnvH1D(Proj_ANNConfidence_vs_NHitsQE_MC)
        Proj_ANNConfidence_vs_NHitsRES_MC = PlotUtils.MnvH1D(Proj_ANNConfidence_vs_NHitsRES_MC)
        Proj_ANNConfidence_vs_NHitsDIS_MC = PlotUtils.MnvH1D(Proj_ANNConfidence_vs_NHitsDIS_MC)
        Proj_ANNConfidence_vs_NHitsCOHPI_MC = PlotUtils.MnvH1D(Proj_ANNConfidence_vs_NHitsCOHPI_MC)
        Proj_ANNConfidence_vs_NHitsAMNUGAMMA_MC = PlotUtils.MnvH1D(Proj_ANNConfidence_vs_NHitsAMNUGAMMA_MC)
        Proj_ANNConfidence_vs_NHitsIMD_MC = PlotUtils.MnvH1D(Proj_ANNConfidence_vs_NHitsIMD_MC)
        Proj_ANNConfidence_vs_NHitsQENUEEL_MC = PlotUtils.MnvH1D(Proj_ANNConfidence_vs_NHitsQENUEEL_MC)
        Proj_ANNConfidence_vs_NHitsQEMEC_MC = PlotUtils.MnvH1D(Proj_ANNConfidence_vs_NHitsQEMEC_MC)
        Proj_ANNConfidence_vs_NHitsData = PlotUtils.MnvH1D(Proj_ANNConfidence_vs_NHitsData)

        Proj_ANNConfidence_vs_NHitsQE_MC.SetTitle("QE")
        Proj_ANNConfidence_vs_NHitsRES_MC.SetTitle("RES")
        Proj_ANNConfidence_vs_NHitsDIS_MC.SetTitle("DIS")
        Proj_ANNConfidence_vs_NHitsCOHPI_MC.SetTitle("COHPI")
        Proj_ANNConfidence_vs_NHitsAMNUGAMMA_MC.SetTitle("AMNUGAMMA")
        Proj_ANNConfidence_vs_NHitsIMD_MC.SetTitle("IMD")
        Proj_ANNConfidence_vs_NHitsQENUEEL_MC.SetTitle("QENUEEL")
        Proj_ANNConfidence_vs_NHitsQEMEC_MC.SetTitle("QEMEC")

        Cumulative_Proj_ANNConfidence_vs_NHitsQE_MC.SetTitle("QE")
        Cumulative_Proj_ANNConfidence_vs_NHitsRES_MC.SetTitle("RES")
        Cumulative_Proj_ANNConfidence_vs_NHitsDIS_MC.SetTitle("DIS")
        Cumulative_Proj_ANNConfidence_vs_NHitsCOHPI_MC.SetTitle("COHPI")
        Cumulative_Proj_ANNConfidence_vs_NHitsAMNUGAMMA_MC.SetTitle("AMNUGAMMA")
        Cumulative_Proj_ANNConfidence_vs_NHitsIMD_MC.SetTitle("IMD")
        Cumulative_Proj_ANNConfidence_vs_NHitsQENUEEL_MC.SetTitle("QENUEEL")
        Cumulative_Proj_ANNConfidence_vs_NHitsQEMEC_MC.SetTitle("QEMEC")

        Cumulative_Proj_ANNConfidence_vs_NHitsQE_MC.SetMinimum(0.0)
        Cumulative_Proj_ANNConfidence_vs_NHitsQE_MC.SetMaximum(1.0)

        cumulative_arr = ROOT.TObjArray()
        cumulative_arr.Add(Cumulative_Proj_ANNConfidence_vs_NHitsQE_MC)
        cumulative_arr.Add(Cumulative_Proj_ANNConfidence_vs_NHitsRES_MC)
        cumulative_arr.Add(Cumulative_Proj_ANNConfidence_vs_NHitsDIS_MC)
        cumulative_arr.Add(Cumulative_Proj_ANNConfidence_vs_NHitsCOHPI_MC)
        #Unused entirely
        #cumulative_arr.Add(Cumulative_Proj_ANNConfidence_vs_NHitsAMNUGAMMA_MC)
        #cumulative_arr.Add(Cumulative_Proj_ANNConfidence_vs_NHitsIMD_MC)
        #cumulative_arr.Add(Cumulative_Proj_ANNConfidence_vs_NHitsQENUEEL_MC)
        cumulative_arr.Add(Cumulative_Proj_ANNConfidence_vs_NHitsQEMEC_MC)
        
        arr = ROOT.TObjArray()
        arr.Add(Proj_ANNConfidence_vs_NHitsQE_MC)
        arr.Add(Proj_ANNConfidence_vs_NHitsRES_MC)
        arr.Add(Proj_ANNConfidence_vs_NHitsDIS_MC)
        arr.Add(Proj_ANNConfidence_vs_NHitsCOHPI_MC)
        #Unused entirely
        #arr.Add(Proj_ANNConfidence_vs_NHitsAMNUGAMMA_MC)
        #arr.Add(Proj_ANNConfidence_vs_NHitsIMD_MC)
        #arr.Add(Proj_ANNConfidence_vs_NHitsQENUEEL_MC)
        arr.Add(Proj_ANNConfidence_vs_NHitsQEMEC_MC)


        plotter.axis_title_font_x = 43
        plotter.axis_title_font_y = 43
        plotter.axis_title_size_x = 15
        plotter.axis_title_size_y = 15
        plotter.axis_title_offset_x = 1.5
        plotter.axis_title_offset_y = 2

        plotter.legend_text_size = 0.02
        plotter.DrawDataStackedMC(Proj_ANNConfidence_vs_NHitsData, arr, (totalDataPOT/totalMCPOT), "TR", "Data", 2, 1, 3001, "ANN Calibrated Confidence", "ANN Calibrated Confidence")
        plotter.AddPlotLabel(f"ANN confidence (calibrated) for {y_lower:g} < n_hits < {y_higher:g}" , 0.5, 0.97, 20, 1, 43, 22, 0)
        filename = f"{outbase}/ANN_ProbabilityDistribution_Confidence_vs_NHitsMC_{iy}.pdf"
        canvas.SaveAs(filename)
        print("*"*30)
        tmpMC = Proj_ANNConfidence_vs_NHitsMC.Clone()
        tmpMC.Scale(totalDataPOT/totalMCPOT)
        chisq = Proj_ANNConfidence_vs_NHitsData.Chi2Test(tmpMC, "UW")
        print(f"ChiSq for {y_lower:g} < n_hits < {y_higher:g} = {chisq}")

        plotter.axis_minimum=0
        plotter.axis_maximum=1

        plotter.legend_text_size = 0.02
        print(f"max_MC_val: {max_MC_val}")
        plotter.DrawDataStackedMC(Cumulative_Proj_ANNConfidence_vs_NHitsData, cumulative_arr, (1/max_MC_val), "TR", "Data", 2, 1, 3001, "ANN Calibrated Confidence", "Proportion of events passing")
        plotter.AddPlotLabel(f"Proportion of events passing ANN confidence cut for {y_lower:g} < n_hits < {y_higher:g}" , 0.5, 0.97, 20, 1, 43, 22, 0)
        filename = f"{outbase}/ANN_PassingEvents_Confidence_vs_NHitsMC_{iy}.pdf"
        canvas.SaveAs(filename)

        #plotter.DrawDataStackedMC(Cumulative_Proj_ANNConfidence_vs_NHitsData, cumulative_arr, (1/max_MC_val), "TR", "Data", 2, 1, 3001, "ANN Calibrated Confidence", "Proportion of events passing")
        #plotter.AddPlotLabel(f"Proportion of events passing ANN confidence cut for {y_lower:g} < n_hits < {y_higher:g}" , 0.5, 0.97, 20, 1, 43, 22, 0)
        #filename = f"{outbase}/ANN_PassingEvents_Confidence_vs_NHitsMC_{iy}.pdf"
        #canvas.SaveAs(filename)

        summedMC = Proj_ANNConfidence_vs_NHitsQE_MC.Clone()
        summedMC.Add(Proj_ANNConfidence_vs_NHitsRES_MC)
        summedMC.Add(Proj_ANNConfidence_vs_NHitsDIS_MC)
        summedMC.Add(Proj_ANNConfidence_vs_NHitsCOHPI_MC)
        summedMC.Add(Proj_ANNConfidence_vs_NHitsQEMEC_MC)

        plotter.DrawDataMCRatio(Proj_ANNConfidence_vs_NHitsData, summedMC, (totalDataPOT/totalMCPOT), True, True, 0, 2)
        plotter.AddPlotLabel(f"Data MC Ratio of ANN confidences for {y_lower:g} < n_hits < {y_higher:g}" , 0.5, 0.97, 20, 1, 43, 22, 0)
        filename = f"{outbase}/ANN_DataMCRatio_Confidence_vs_NHitsMC_{iy}.pdf"
        canvas.SaveAs(filename)

    #================================================================================================================
    #End - =========================  Producing plots broken down by number of hits  #===============================
    #================================================================================================================
    #=============================  Producing confidence distributons by ERecoil  #=================================
    #================================================================================================================

    erecoil_decile_limits = get_quantile_bins(ANNZ_vs_Conf_vs_Erecoil_Data, 10, 2) # These will also be useful later on
    # You can print erecoil_decile_limits here to verify the bin distributions 
    # before generating the plots.

    data_by_Erecoil_deciles = extract_var2_decile_hists(ANNZ_vs_Conf_vs_Erecoil_Data, erecoil_decile_limits)
    mc_by_Erecoil_deciles = extract_var2_decile_hists(ANNZ_vs_Conf_vs_Erecoil_MC, erecoil_decile_limits)
    for i in range(len(data_by_Erecoil_deciles)):
        dataHist = data_by_Erecoil_deciles[i]
        mcHist = mc_by_Erecoil_deciles[i]
        if dataHist.Integral()>0:
            dataHist.Scale(1.0/dataHist.Integral())
        if mcHist.Integral()>0:
            mcHist.Scale(1.0/mcHist.Integral())
        dataHist.SetTitle(f"Distribution of event confidences at Erecoil bin {i}")
        dataHist.SetLineColor(ROOT.kBlue)
        mcHist.SetLineColor(ROOT.kRed)
        #dataHist.Draw("HIST")
        mcHist.Draw("HIST")
        filename = f"{outbase}/EHad_bin{i}_ProbDist.pdf"
        canvas.SaveAs(filename)
        passingData = dataHist.GetCumulative(False)
        passingMC = mcHist.GetCumulative(False)
        passingData.SetTitle(f"Distribution of event confidences at Erecoil bin {i}")
        passingData.SetLineColor(ROOT.kBlue)
        passingMC.SetLineColor(ROOT.kRed)
        passingData.Draw("HIST")
        passingMC.Draw("HIST SAME")
        filename = f"{outbase}/EHad_bin{i}_PassRateDist.pdf"
        canvas.SaveAs(filename)
    #================================================================================================================
    #End - =======================  Producing confidencen distributons by ERecoil  #=================================
    #================================================================================================================
    #===================================  Producing vertex distribution plots  #=====================================
    #================================================================================================================

    empty_vtxZ_data = ANNZ_vs_Conf_vs_Erecoil_Data_Empty.ProjectionX("VtxZDataEmpty", -1, -1, -1, -1)
    empty_vtxZ_mc = ANNZ_vs_Conf_vs_Erecoil_MC_Empty.ProjectionX("VtxZMCEmpty", -1, -1, -1, -1)
    filled_vtxZ_data = ANNZ_vs_Conf_vs_Erecoil_Data_Filled.ProjectionX("VtxZDataFilled", -1, -1, -1, -1)
    filled_vtxZ_mc = ANNZ_vs_Conf_vs_Erecoil_MC_Filled.ProjectionX("VtxZMCFilled", -1, -1, -1, -1)
    for hist in [empty_vtxZ_data, empty_vtxZ_mc, filled_vtxZ_data, filled_vtxZ_mc]:
        if hist.Integral()>0:
            hist.Scale(1.0/dataHist.Integral())
    canvas.Clear()
    title = f"Z Vertex Distribution Data/MC Ratios"
    DrawFullAndEmptyDataMC(canvas, filled_vtxZ_data, filled_vtxZ_mc, empty_vtxZ_data, empty_vtxZ_mc, title, "lin", filledDataPOT, filledMCPOT, emptyDataPOT, emptyMCPOT)
    filename = f"{outbase}/Overall_Z_vertex_distributions.pdf"
    canvas.SaveAs(filename)
    canvas.Clear()
    DrawFullAndEmptyDataMC(canvas, filled_vtxZ_data, filled_vtxZ_mc, empty_vtxZ_data, empty_vtxZ_mc, title, "rat", filledDataPOT, filledMCPOT, emptyDataPOT, emptyMCPOT)
    filename = f"{outbase}/Overall_Z_vertex_distributions_ratio.pdf"
    canvas.SaveAs(filename)
    canvas.Clear()
    DrawFullAndEmptyDataMC(canvas, filled_vtxZ_data, filled_vtxZ_mc, empty_vtxZ_data, empty_vtxZ_mc, title, "lin", filledDataPOT, filledMCPOT, emptyDataPOT, emptyMCPOT, zoomIntoTargets=True)
    filename = f"{outbase}/Overall_Z_vertex_distributions_zoomed.pdf"
    canvas.SaveAs(filename)
    canvas.Clear()
    DrawFullAndEmptyDataMC(canvas, filled_vtxZ_data, filled_vtxZ_mc, empty_vtxZ_data, empty_vtxZ_mc, title, "rat", filledDataPOT, filledMCPOT, emptyDataPOT, emptyMCPOT, zoomIntoTargets=True)
    filename = f"{outbase}/Overall_Z_vertex_distributions_ratio_zoomed.pdf"
    canvas.SaveAs(filename)

    erecoil_quantile_limits = get_quantile_bins(ANNZ_vs_Conf_vs_Erecoil_Data, 4, 2) # These will also be useful later on

    empty_vtxZ_data_by_Erecoil_deciles = extract_var1_decile_hists(ANNZ_vs_Conf_vs_Erecoil_Data_Empty, erecoil_quantile_limits)
    empty_vtxZ_mc_by_Erecoil_deciles = extract_var1_decile_hists(ANNZ_vs_Conf_vs_Erecoil_MC_Empty, erecoil_quantile_limits)
    filled_vtxZ_data_by_Erecoil_deciles = extract_var1_decile_hists(ANNZ_vs_Conf_vs_Erecoil_Data_Filled, erecoil_quantile_limits)
    filled_vtxZ_mc_by_Erecoil_deciles = extract_var1_decile_hists(ANNZ_vs_Conf_vs_Erecoil_MC_Filled, erecoil_quantile_limits)
    for i in range(len(empty_vtxZ_mc_by_Erecoil_deciles)):
        canvas.Clear()
        emptydataHist = empty_vtxZ_data_by_Erecoil_deciles[i]
        emptymcHist = empty_vtxZ_mc_by_Erecoil_deciles[i]
        filleddataHist = filled_vtxZ_data_by_Erecoil_deciles[i]
        filledmcHist = filled_vtxZ_mc_by_Erecoil_deciles[i]
        for hist in [emptydataHist, emptymcHist, filleddataHist, filledmcHist]:
            if hist.Integral()>0:
                hist.Scale(1.0/dataHist.Integral())
        binLowerLim = -1
        binUpperLim = -1
        title = f"Z vertex distribution for "
        if i != 0:
            title += f"{erecoil_quantile_limits[i]:.2f} < "
        title += "E_{recoil}"
        if i != len(erecoil_quantile_limits)-1:
            title += f" < {erecoil_quantile_limits[i+1]:.2f}"
        dataHist.SetTitle(title)
        canvas.Clear()
        DrawFullAndEmptyDataMC(canvas, filleddataHist, filledmcHist, emptydataHist, emptymcHist, title, "lin", filledDataPOT, filledMCPOT, emptyDataPOT, emptyMCPOT)
        filename = f"{outbase}/EHad_bin{i}_Z_vertex_distributions.pdf"
        canvas.SaveAs(filename)
        canvas.Clear()
        DrawFullAndEmptyDataMC(canvas, filleddataHist, filledmcHist, emptydataHist, emptymcHist, title, "rat", filledDataPOT, filledMCPOT, emptyDataPOT, emptyMCPOT, yminFull=0.5, yminEmpty=0.5, ymaxFull=1.5, ymaxEmpty=1.5)
        filename = f"{outbase}/EHad_bin{i}_Z_vertex_distributions_ratio.pdf"
        canvas.SaveAs(filename)
        canvas.Clear()
        DrawFullAndEmptyDataMC(canvas, filleddataHist, filledmcHist, emptydataHist, emptymcHist, title, "log", filledDataPOT, filledMCPOT, emptyDataPOT, emptyMCPOT)
        filename = f"{outbase}/EHad_bin{i}_Z_vertex_distributions_log.pdf"
        canvas.SaveAs(filename)

if __name__ == "__main__":
    parser = argparse.ArgumentParser( prog='PlotValidations',
                    description='Extracts and produces a number of useful plot from the VertexValidation.root files')
    parser.add_argument('input', help="the path of a directory containing a number of input files or the path of a single input file")
    parser.add_argument('out', help="the output directory")
    parser.add_argument('nu', help="Neutrino or antineutrino mode")
    args = parser.parse_args()
    PlotValidation(args.input, args.out, bool(int(args.nu)))