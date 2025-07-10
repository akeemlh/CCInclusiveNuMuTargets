## import common python modules
import ROOT,PlotUtils,os,sys, argparse
import csv
import numpy as np
import ctypes

import ROOT
from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, THStack, TObjArray, TRatioPlot, TLegend,TLine, TParameter
from ROOT import gROOT
import copy
from array import array

## modify $PYTHONPATH to know about custom python modules
PLOTUTILSROOT = os.path.expandvars("$PLOTUTILSROOT")
ROOT.gROOT.SetBatch(1)

rebinNum = 10
nbins = 3400
xlow = 4200
#xhigh = 5900 # Zoom into target region
xhigh = 8410 # Full detector

def PlotRegions(inYmin = None, inYmax = None, mcpot = None, datapot = None):

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
    print("Finished plotting regions")

def getChi2(val):
    ndof = ctypes.c_int(0)
    """ chisq = ctypes.c_double(0)
    ndof = ctypes.c_int(0)
    igood = ctypes.c_int(0)
    mc = MChist.Clone("temp_plastic_float")
    print("val[0] ", val[0])
    print("mc.GetBinContent(10) ", mc.GetBinContent(15))
    mc = mc*val[0]
    print("after mc.GetBinContent(10) ", mc.GetBinContent(15))
    data = DataHist
    data.Chi2TestX(mc, chisq, ndof ,igood , "UW")
    #return chisq.value
    return data.Chi2Test(mc, "UW") """
    data = DataHist.Clone("data")
    mc = MChist.Clone("mc")
    """ mc.Scale(val[0])
    nbinsdata= data.GetNbinsX()
    chisq = 0
    for i in range(nbinsdata):
        dataval = data.GetBinContent(i)
        mcval = mc.GetBinContent(i)
        if mcval == 0:
            continue
        chisq = chisq + (dataval-mcval)**2/mcval
    return chisq """
    plotter = PlotUtils.MnvPlotter(6)
    chisq = plotter.Chi2DataMC(data, mc, ndof, val[0])
    #print("chisq ", chisq)

    #print("ndof ", ndof.value)
    return chisq

global MChist
MChist = 0
global DataHist
DataHist = 0

def PlotSidebands(dataSideband, MCSidebandSignal, MCSidebandUS, MCSidebandDS, MCSidebandOther, scaling, title, filepath):

    print("Plotting sidebands")
    canvas = ROOT.TCanvas("canvas", "canvas", 1500, 1500)
    """ 
    p1 = ROOT.TPad("p1","top",0,0.3,1,0.95) # xlow,yloy, xup,yup
    p2 = ROOT.TPad("p2","bottom",0,0,1,0.3)
    tex3 = ROOT.TLatex(.5,.95,"Title2")
    tex3.SetTextAlign(22)
    tex3.SetTextSize(200)
    tex3.DrawLatexNDC(.5,.95, title)
    #ROOT.gStyle.SetOptTile(False)

    ROOT.gPad.Update()
    canvas.cd()
    p1.SetFillStyle(0)
    p2.SetFillStyle(0)
    p1.Draw()
    p1.cd()
    stack = THStack("hs","")
    #DataUSSideband.Draw("P SAME")
    plotter = PlotUtils.MnvPlotter(6) #kCCQENuStyle
    arr = TObjArray()
    arr.Add(MCSidebandSignal)
    arr.Add(MCSidebandUS)
    arr.Add(MCSidebandDS)
    arr.Add(MCSidebandOther)
    summedSideband = (USSidebandSignal+USSidebandUS+USSidebandDS+USSidebandOther)*dataMCScale
    colours = (ctypes.c_int*4)(1,2,3,4)
    plotter.DrawDataStackedMC(dataSideband, arr, colours,scaling, 'R', "Data", 3001, "x test", "y test")

    canvas.cd()
    p2.Draw()
    p2.cd()
    plotter.DrawDataMCRatio(dataSideband, summedSideband) """

    stack = THStack("hs",title)
    MCSidebandSignalTemp = MCSidebandSignal.Clone("MCSidebandSignalTemp")
    MCSidebandUSTemp = MCSidebandUS.Clone("MCSidebandUSTemp")
    MCSidebandDSTemp = MCSidebandDS.Clone("MCSidebandDSTemp")
    MCSidebandOtherTemp = MCSidebandOther.Clone("MCSidebandOtherTemp")

    MCSidebandSignalTemp.Scale(scaling)
    MCSidebandUSTemp.Scale(scaling)
    MCSidebandDSTemp.Scale(scaling)
    MCSidebandOtherTemp.Scale(scaling)

    MCSidebandSignalTemp.SetFillColor(1)
    MCSidebandUSTemp.SetFillColor(2)
    MCSidebandDSTemp.SetFillColor(3)
    MCSidebandOtherTemp.SetFillColor(4)

    stack.Add(MCSidebandSignalTemp)
    stack.Add(MCSidebandUSTemp)
    stack.Add(MCSidebandDSTemp)
    stack.Add(MCSidebandOtherTemp)

    dataSideband.SetMarkerStyle(8)

    ratioPlot = TRatioPlot(stack, dataSideband)

    ratioPlot.SetSeparationMargin(0.0)
    ratioPlot.Draw("hideup")
    ratioPlot.GetUpperRefXaxis().SetRangeUser(0, 120)
    ratioPlot.GetLowerRefYaxis().SetRangeUser(0.5, 1.5)

    ratioPlot.GetLowerRefGraph().SetMinimum(0.5)
    ratioPlot.GetLowerRefGraph().SetMaximum(1.5)
    ratioPlot.GetLowerRefYaxis().SetTitle("ratio")
    ratioPlot.GetUpperRefYaxis().SetTitle("entries")
    ratioPlot.GetLowerRefGraph().SetMarkerStyle(8)
    ratioPlot.GetLowerRefGraph().SetMarkerSize(2)
    ratioPlot.GetLowerRefGraph().GetYaxis().SetNdivisions(5)
    ratioPlot.GetLowerRefGraph().GetXaxis().SetLabelSize(0.03)
    #ratioPlot.GetLowerPad().cd()
    #l = TLine(0,1,120,1)
    #l.Draw()
    ratioPlot.GetUpperPad().cd()
    legend = TLegend(0.6, 0.7, 0.8, 0.85)
    
    legend.AddEntry("MCSidebandSignalTemp", "Signal", "f")
    legend.AddEntry("MCSidebandUSTemp", "Upstream Sideband", "f")
    legend.AddEntry("MCSidebandDSTemp", "Downstream Sideband", "f")
    legend.AddEntry("MCSidebandOtherTemp", "Other", "f")
    legend.Draw()

    canvas.Update()

    canvas.SaveAs(filepath)

    canvas.Clear()


def minimiseChiSq():
    minimizer = ROOT.Math.Factory.CreateMinimizer("Minuit2", "")
    if (not minimizer):
        raise RuntimeError(
            "Cannot create minimizer Minuit2. Maybe the required library was not built?")
    
    # Set tolerance and other minimizer parameters, one can also use default
    # values

    minimizer.SetMaxFunctionCalls(1000000)  # working for Minuit/Minuit2
    # for GSL minimizers - no effect in Minuit/Minuit2
    minimizer.SetMaxIterations(10000)
    minimizer.SetTolerance(0.001)
    minimizer.SetPrintLevel(1)

    # Create function wrapper for minimizer

    f = ROOT.Math.Functor(getChi2, 1)

    minimizer.SetFunction(f)
    # Set the free variables to be minimized !
    minimizer.SetVariable(0, "scale", 1, 0.001)
    
    # Do the minimization
    ret = minimizer.Minimize()

    xs = minimizer.X()
    #print("Minimum: f({}) = {}".format(xs[0],minimizer.MinValue()))
    return xs[0]

if __name__ == '__main__':
    parser = argparse.ArgumentParser( prog='SidebandTuning',
                    description='automates the sideband tuning process')
    parser.add_argument('data', nargs=1, help="the data input file")
    parser.add_argument('mc', nargs=1, help="the data input directory")
    #parser.add_argument('out', nargs=1, help="Where to put the scale factors and assosciated chisq")
    args = parser.parse_args()


    MCFile = args.mc[0] #"/pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/7July2025P4/1A/lead_mergedMC.root"
    DataFile = args.data[0] #"/pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/7July2025P4/1A/lead_mergedData.root"
    mcfile = TFile(MCFile)
    datafile = TFile(DataFile)
    #outbase = args.out[0]
    mcextension = MCFile.find(".root")
    dataextension = DataFile.find(".root")
    mcbase = MCFile[:mcextension]
    database = DataFile[:dataextension]

    USSidebandSignal = mcfile.Get("nuke_segment_US_sideband_Signal")
    USSidebandDS = mcfile.Get("nuke_segment_US_sideband_DS")
    USSidebandUS = mcfile.Get("nuke_segment_US_sideband_US")
    USSidebandOther = mcfile.Get("nuke_segment_US_sideband_Other")
    DSSidebandSignal = mcfile.Get("nuke_segment_DS_sideband_Signal")
    DSSidebandDS = mcfile.Get("nuke_segment_DS_sideband_DS")
    DSSidebandUS = mcfile.Get("nuke_segment_DS_sideband_US")
    DSSidebandOther = mcfile.Get("nuke_segment_DS_sideband_Other")
    DataUSSideband = datafile.Get("nuke_segment_US_Sideband")
    DataDSSideband = datafile.Get("nuke_segment_DS_Sideband")
    DataSignal = datafile.Get("nuke_segment_data")
    dataPOT = datafile.Get("POTUsed").GetVal()
    mcPOT = mcfile.Get("POTUsed").GetVal()
    dataMCScale = dataPOT/mcPOT

    PlotSidebands(DataUSSideband, USSidebandSignal, USSidebandUS, USSidebandDS, USSidebandOther, dataMCScale, "Upstream Sideband - Untuned", mcbase+"USUntuned.png")
    PlotSidebands(DataDSSideband, DSSidebandSignal, DSSidebandUS, DSSidebandDS, DSSidebandOther, dataMCScale, "Downstream Sideband - Untuned", mcbase+"DSUntuned.png")
    ##Scaling stuff
    ##Find US Scaling
    combined = USSidebandSignal.Clone("Combined")
    combined.Add(USSidebandUS)
    combined.Add(USSidebandDS)
    combined.Add(USSidebandOther)
    combined.Scale(dataMCScale)
    DataHist = DataUSSideband
    MChist = combined
    data = DataHist.Clone("data")
    mc = MChist.Clone("mc")
    usscaling = minimiseChiSq()

    ##Find US Scaling
    combined = DSSidebandSignal.Clone("Combined")
    combined.Add(DSSidebandUS)
    combined.Add(DSSidebandDS)
    combined.Add(DSSidebandOther)
    combined.Scale(dataMCScale)
    DataHist = DataDSSideband
    MChist = combined
    data = DataHist.Clone("data")
    mc = MChist.Clone("mc")
    dsscaling = minimiseChiSq()
    PlotSidebands(DataUSSideband, USSidebandSignal, USSidebandUS, USSidebandDS, USSidebandOther, dataMCScale*usscaling, "Upstream Sideband - Tuned", mcbase+"USTuned.png")
    PlotSidebands(DataDSSideband, DSSidebandSignal, DSSidebandUS, DSSidebandDS, DSSidebandOther, dataMCScale*dsscaling, "Downstream Sideband - Tuned", mcbase+"DSTuned.png")

    mcoutpath = mcbase+"Scaled.root"
    print("mcoutpath: ",mcoutpath)
    mcoutfile = TFile(mcoutpath, "RECREATE")
    mcoutfile.cd()
    
    for key in mcfile.GetListOfKeys():
        #print("key ", key)
        _name = key.GetName()
        obj = mcfile.Get(_name)
        if "background_US" in _name:
            obj.Scale(usscaling)
        if "background_DS" in _name:
            obj.Scale(dsscaling)
        mcoutfile.WriteTObject(obj,_name)

    mcoutfile.Close()

    dataoutpath = database+"Scaled.root"
    print("dataoutpath: ",dataoutpath)
    dataoutfile = TFile(dataoutpath, "RECREATE")
    dataoutfile.cd()

    for key in datafile.GetListOfKeys():
        #print("key ", key)
        _name = key.GetName()
        obj = datafile.Get(_name)
        if "background_US" in _name:
            obj.Scale(usscaling)
        if "background_DS" in _name:
            obj.Scale(dsscaling)
        dataoutfile.WriteTObject(obj,_name)

    dataoutfile.Close()

    mcfile.Close()
    datafile.Close()

