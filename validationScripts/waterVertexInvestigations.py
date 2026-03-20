import os, argparse, math
import ROOT
from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, THStack, TObjArray, TRatioPlot, TLegend,TLine, TParameter
from ROOT import gROOT
from ROOT import PlotUtils
import ctypes

import numpy as np
#import PlotUtils does this work too?
from array import array   


def PlotRegions(bySegment=False, inYmin = None, inYmax = None, mcpot = None, datapot = None):

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

    if not bySegment:
        tgt1 = ROOT.TBox(tgt1centre-tgt1thickness/2, ymin, tgt1centre+tgt1thickness/2, ymax)
        tgt2 = ROOT.TBox(tgt2centre-tgt2thickness/2, ymin, tgt2centre+tgt2thickness/2, ymax)
        tgt3 = ROOT.TBox(tgt3centre-tgt3thickness/2, ymin, tgt3centre+tgt3thickness/2, ymax)
        tgt4 = ROOT.TBox(tgt4centre-tgt4thickness/2, ymin, tgt4centre+tgt4thickness/2, ymax)
        tgt5 = ROOT.TBox(tgt5centre-tgt5thickness/2, ymin, tgt5centre+tgt5thickness/2, ymax)
        tgt6 = ROOT.TBox(5200, ymin, 5420, ymax)
    else:
        tgt1 = ROOT.TBox(8, ymin, 9, ymax)
        tgt2 = ROOT.TBox(17, ymin, 18, ymax)
        tgt3 = ROOT.TBox(26, ymin, 27, ymax)
        tgt4 = ROOT.TBox(44, ymin, 45, ymax)
        tgt5 = ROOT.TBox(49, ymin, 50, ymax)
        tgt6 = ROOT.TBox(35, ymin, 36, ymax)

    tgt1.SetFillColorAlpha(ROOT.kBlack,0.5)
    tgt1.SetFillStyle(3004)
    tgt1.SetLineWidth(0)
    tgt1.SetLineColor(1)
    if not bySegment:
        tgt1.DrawBox(tgt1centre-tgt1thickness/2, ymin, tgt1centre+tgt1thickness/2, ymax)
    else:
        tgt1.DrawBox(8.5, ymin, 9.5, ymax)
    tgt2.SetFillColorAlpha(ROOT.kBlack,0.5)
    tgt2.SetFillStyle(3004)
    tgt2.SetLineWidth(0)
    tgt2.SetLineColor(1)
    if not bySegment:
        tgt2.DrawBox(tgt2centre-tgt2thickness/2, ymin, tgt2centre+tgt2thickness/2, ymax)
    else:
        tgt1.DrawBox(17.5, ymin, 18.5, ymax)
    tgt3.SetFillColorAlpha(ROOT.kBlack,0.5)
    tgt3.SetFillStyle(3004)
    tgt3.SetLineWidth(0)
    tgt3.SetLineColor(1)
    if not bySegment:
        tgt3.DrawBox(tgt3centre-tgt3thickness/2, ymin, tgt3centre+tgt3thickness/2, ymax)
    else:
        tgt1.DrawBox(26.5, ymin, 27.5, ymax)
    tgt4.SetFillColorAlpha(ROOT.kBlack,0.5)
    tgt4.SetFillStyle(3004)
    tgt4.SetLineWidth(0)
    tgt4.SetLineColor(1)
    if not bySegment:
        tgt4.DrawBox(tgt4centre-tgt4thickness/2, ymin, tgt4centre+tgt4thickness/2, ymax)
    else:
        tgt4.DrawBox(44.5, ymin, 45.5, ymax)
    tgt5.SetFillColorAlpha(ROOT.kBlack,0.5)
    tgt5.SetFillStyle(3004)
    tgt5.SetLineWidth(0)
    tgt5.SetLineColor(1)
    if not bySegment:
        tgt5.DrawBox(tgt5centre-tgt5thickness/2, ymin, tgt5centre+tgt5thickness/2, ymax)
    else:
        tgt5.DrawBox(49.5, ymin, 50.5, ymax)
    tgt6.SetFillColorAlpha(ROOT.kBlue,0.5)
    tgt6.SetFillStyle(3004)
    tgt6.SetLineWidth(0)
    tgt6.SetLineColor(1)        
    if not bySegment:
        tgt6.DrawBox(5200, ymin, 5420, ymax)
    else:
        tgt6.DrawBox(35.5, ymin, 36.5, ymax)
    print("Finished plotting regions")
        

def zPosFromModAndPlane(module, plane):
    if (module == -5 and plane == 1):
        return 4293.04
    elif (module == -5 and plane == 2):
        return 4313.68
    elif (module == -4 and plane == 1):
        return 4337.25
    elif (module == -4 and plane == 2):
        return 4357.9
    elif (module == -3 and plane == 1):
        return 4381.47
    elif (module == -3 and plane == 2):
        return 4402.11
    elif (module == -2 and plane == 1):
        return 4425.68
    elif (module == -2 and plane == 2):
        return 4446.33
    elif (module == 0 and plane == 1):
        return 4514.11
    elif (module == 0 and plane == 2):
        return 4534.76
    elif (module == 1 and plane == 1):
        return 4558.33
    elif (module == 1 and plane == 2):
        return 4578.97
    elif (module == 2 and plane == 1):
        return 4602.54
    elif (module == 2 and plane == 2):
        return 4623.19
    elif (module == 3 and plane == 1):
        return 4646.76
    elif (module == 3 and plane == 2):
        return 4667.4
    elif (module == 5 and plane == 1):
        return 4735.19
    elif (module == 5 and plane == 2):
        return 4755.83
    elif (module == 6 and plane == 1):
        return 4779.4
    elif (module == 6 and plane == 2):
        return 4800.05
    elif (module == 7 and plane == 1):
        return 4823.62
    elif (module == 7 and plane == 2):
        return 4844.26
    elif (module == 8 and plane == 1):
        return 4867.83
    elif (module == 8 and plane == 2):
        return 4888.48
    elif (module == 11 and plane == 1):
        return 5000.48
    elif (module == 11 and plane == 2):
        return 5021.12
    elif (module == 12 and plane == 1):
        return 5044.69
    elif (module == 12 and plane == 2):
        return 5065.34
    elif (module == 13 and plane == 1):
        return 5088.91
    elif (module == 13 and plane == 2):
        return 5109.55
    elif (module == 14 and plane == 1):
        return 5133.12
    elif (module == 14 and plane == 2):
        return 5153.77
    elif (module == 15 and plane == 1):
        return 5456.74
    elif (module == 15 and plane == 2):
        return 5477.38
    elif (module == 16 and plane == 1):
        return 5500.95
    elif (module == 16 and plane == 2):
        return 5521.6
    elif (module == 17 and plane == 1):
        return 5545.17
    elif (module == 17 and plane == 2):
        return 5565.81
    elif (module == 18 and plane == 1):
        return 5589.38
    elif (module == 18 and plane == 2):
        return 5610.02
    elif (module == 20 and plane == 1):
        return 5677.81
    elif (module == 20 and plane == 2):
        return 5698.45
    elif (module == 21 and plane == 1):
        return 5722.03
    elif (module == 21 and plane == 2):
        return 5742.67
    elif (module == 23 and plane == 1):
        return 5810.45
    elif (module == 23 and plane == 2):
        return 5831.1
    elif (module == 24 and plane == 1):
        return 5855.68
    elif (module == 24 and plane == 2):
        return 5876.33
    elif (module == 25 and plane == 1):
        return 5900.91
    elif (module == 25 and plane == 2):
        return 5921.56
    elif (module == 26 and plane == 1):
        return 5946.14
    elif (module == 26 and plane == 2):
        return 5966.79
    elif (module == 27 and plane == 1):
        return 5991.37
    elif (module == 27 and plane == 2):
        return 6012.01
    elif (module == 28 and plane == 1):
        return 6036.6
    elif (module == 28 and plane == 2):
        return 6057.24
    elif (module == 29 and plane == 1):
        return 6081.83
    elif (module == 29 and plane == 2):
        return 6102.47
    elif (module == 30 and plane == 1):
        return 6127.06
    elif (module == 30 and plane == 2):
        return 6147.7
    elif (module == 31 and plane == 1):
        return 6172.29
    elif (module == 31 and plane == 2):
        return 6192.93
    elif (module == 32 and plane == 1):
        return 6217.52
    elif (module == 32 and plane == 2):
        return 6238.16
    elif (module == 33 and plane == 1):
        return 6262.74
    elif (module == 33 and plane == 2):
        return 6283.39
    elif (module == 34 and plane == 1):
        return 6307.97
    elif (module == 34 and plane == 2):
        return 6328.62
    elif (module == 35 and plane == 1):
        return 6353.2
    elif (module == 35 and plane == 2):
        return 6373.85
    elif (module == 36 and plane == 1):
        return 6398.43
    elif (module == 36 and plane == 2):
        return 6419.08
    elif (module == 37 and plane == 1):
        return 6443.66
    elif (module == 37 and plane == 2):
        return 6464.3
    elif (module == 38 and plane == 1):
        return 6488.89
    elif (module == 38 and plane == 2):
        return 6509.53
    elif (module == 39 and plane == 1):
        return 6534.12
    elif (module == 39 and plane == 2):
        return 6554.76
    elif (module == 40 and plane == 1):
        return 6579.35
    elif (module == 40 and plane == 2):
        return 6599.99
    elif (module == 41 and plane == 1):
        return 6624.58
    elif (module == 41 and plane == 2):
        return 6645.22
    elif (module == 42 and plane == 1):
        return 6669.81
    elif (module == 42 and plane == 2):
        return 6690.45
    elif (module == 43 and plane == 1):
        return 6715.03
    elif (module == 43 and plane == 2):
        return 6735.68
    elif (module == 44 and plane == 1):
        return 6760.26
    elif (module == 44 and plane == 2):
        return 6780.91
    elif (module == 45 and plane == 1):
        return 6805.49
    elif (module == 45 and plane == 2):
        return 6826.14
    elif (module == 46 and plane == 1):
        return 6850.72
    elif (module == 46 and plane == 2):
        return 6871.37
    elif (module == 47 and plane == 1):
        return 6895.95
    elif (module == 47 and plane == 2):
        return 6916.59
    elif (module == 48 and plane == 1):
        return 6941.18
    elif (module == 48 and plane == 2):
        return 6961.82
    elif (module == 49 and plane == 1):
        return 6986.41
    elif (module == 49 and plane == 2):
        return 7007.05
    elif (module == 50 and plane == 1):
        return 7031.64
    elif (module == 50 and plane == 2):
        return 7052.28
    elif (module == 51 and plane == 1):
        return 7076.87
    elif (module == 51 and plane == 2):
        return 7097.51
    elif (module == 52 and plane == 1):
        return 7122.1
    elif (module == 52 and plane == 2):
        return 7142.74
    elif (module == 53 and plane == 1):
        return 7167.32
    elif (module == 53 and plane == 2):
        return 7187.97
    elif (module == 54 and plane == 1):
        return 7212.55
    elif (module == 54 and plane == 2):
        return 7233.2
    elif (module == 55 and plane == 1):
        return 7257.78
    elif (module == 55 and plane == 2):
        return 7278.43
    elif (module == 56 and plane == 1):
        return 7303.01
    elif (module == 56 and plane == 2):
        return 7323.66
    elif (module == 57 and plane == 1):
        return 7348.24
    elif (module == 57 and plane == 2):
        return 7368.88
    elif (module == 58 and plane == 1):
        return 7393.47
    elif (module == 58 and plane == 2):
        return 7414.11
    elif (module == 59 and plane == 1):
        return 7438.7
    elif (module == 59 and plane == 2):
        return 7459.34
    elif (module == 60 and plane == 1):
        return 7483.93
    elif (module == 60 and plane == 2):
        return 7504.57
    elif (module == 61 and plane == 1):
        return 7529.16
    elif (module == 61 and plane == 2):
        return 7549.8
    elif (module == 62 and plane == 1):
        return 7574.39
    elif (module == 62 and plane == 2):
        return 7595.03
    elif (module == 63 and plane == 1):
        return 7619.61
    elif (module == 63 and plane == 2):
        return 7640.26
    elif (module == 64 and plane == 1):
        return 7664.84
    elif (module == 64 and plane == 2):
        return 7685.49
    elif (module == 65 and plane == 1):
        return 7710.07
    elif (module == 65 and plane == 2):
        return 7730.72
    elif (module == 66 and plane == 1):
        return 7755.3
    elif (module == 66 and plane == 2):
        return 7775.95
    elif (module == 67 and plane == 1):
        return 7800.53
    elif (module == 67 and plane == 2):
        return 7821.17
    elif (module == 68 and plane == 1):
        return 7845.76
    elif (module == 68 and plane == 2):
        return 7866.4
    elif (module == 69 and plane == 1):
        return 7890.99
    elif (module == 69 and plane == 2):
        return 7911.63
    elif (module == 70 and plane == 1):
        return 7936.22
    elif (module == 70 and plane == 2):
        return 7956.86
    elif (module == 71 and plane == 1):
        return 7981.45
    elif (module == 71 and plane == 2):
        return 8002.09
    elif (module == 72 and plane == 1):
        return 8026.68
    elif (module == 72 and plane == 2):
        return 8047.32
    elif (module == 73 and plane == 1):
        return 8071.9
    elif (module == 73 and plane == 2):
        return 8092.55
    elif (module == 74 and plane == 1):
        return 8117.13
    elif (module == 74 and plane == 2):
        return 8137.78
    elif (module == 75 and plane == 1):
        return 8162.36
    elif (module == 75 and plane == 2):
        return 8183.01
    elif (module == 76 and plane == 1):
        return 8207.59
    elif (module == 76 and plane == 2):
        return 8228.24
    elif (module == 77 and plane == 1):
        return 8252.82
    elif (module == 77 and plane == 2):
        return 8273.46
    elif (module == 78 and plane == 1):
        return 8298.05
    elif (module == 78 and plane == 2):
        return 8318.69
    elif (module == 79 and plane == 1):
        return 8343.28
    elif (module == 79 and plane == 2):
        return 8363.92
    elif (module == 80 and plane == 1):
        return 8388.51
    elif (module == 80 and plane == 2):
        return 8409.15
    elif (module == 81 and plane == 1):
        return 8433.74
    elif (module == 81 and plane == 2):
        return 8454.38
    elif (module == 82 and plane == 1):
        return 8478.97
    elif (module == 82 and plane == 2):
        return 8499.61
    elif (module == 83 and plane == 1):
        return 8524.19
    elif (module == 83 and plane == 2):
        return 8544.84
    elif (module == 84 and plane == 1):
        return 8569.42
    elif (module == 84 and plane == 2):
        return 8590.07
    elif (module == 85 and plane == 1):
        return 8614.65
    elif (module == 85 and plane == 2):
        return 8635.3
    elif (module == 86 and plane == 1):
        return 8659.46
    elif (module == 86 and plane == 2):
        return 8680.1
    elif (module == 87 and plane == 1):
        return 8704.26
    elif (module == 87 and plane == 2):
        return 8724.9
    elif (module == 88 and plane == 1):
        return 8749.06
    elif (module == 88 and plane == 2):
        return 8769.71
    elif (module == 89 and plane == 1):
        return 8793.86
    elif (module == 89 and plane == 2):
        return 8814.51
    elif (module == 90 and plane == 1):
        return 8838.67
    elif (module == 90 and plane == 2):
        return 8859.31
    elif (module == 91 and plane == 1):
        return 8883.47
    elif (module == 91 and plane == 2):
        return 8904.11
    elif (module == 92 and plane == 1):
        return 8928.27
    elif (module == 92 and plane == 2):
        return 8948.92
    elif (module == 93 and plane == 1):
        return 8973.08
    elif (module == 93 and plane == 2):
        return 8993.72
    elif (module == 94 and plane == 1):
        return 9017.88
    elif (module == 94 and plane == 2):
        return 9038.52
    elif (module == 95 and plane == 2):
        return 9088.08
    elif (module == 96 and plane == 2):
        return 9135.41
    elif (module == 97 and plane == 2):
        return 9182.75
    elif (module == 98 and plane == 2):
        return 9230.08
    elif (module == 99 and plane == 2):
        return 9277.41
    elif (module == 1 and plane == 2):
        return 9324.74
    elif (module == 1 and plane == 2):
        return 9372.08
    elif (module == 1 and plane == 2):
        return 9419.41
    elif (module == 1 and plane == 2):
        return 9466.74
    elif (module == 1 and plane == 2):
        return 9514.07
    elif (module == 1 and plane == 2):
        return 9561.41
    elif (module == 1 and plane == 2):
        return 9608.74
    elif (module == 1 and plane == 2):
        return 9656.07
    elif (module == 1 and plane == 2):
        return 9703.4
    elif (module == 1 and plane == 2):
        return 9750.74
    elif (module == 1 and plane == 2):
        return 9798.07
    elif (module == 1 and plane == 2):
        return 9845.4
    elif (module == 1 and plane == 2):
        return 9892.73
    elif (module == 1 and plane == 2):
        return 9940.07
    elif (module == 1 and plane == 2):
        return 9987.4


def planecode(mdl, plane):
    if (mdl == -5):
            
        if (plane == 1):
            return 1
        if (plane == 2):
            return 2
    if (mdl == -4):
            
        if (plane == 1):
            return 3
        if (plane == 2):
            return 4
    if (mdl == -3):
            
        if (plane == 1):
            return 5
        if (plane == 2):
            return 6
    if (mdl == -2):
            
        if (plane == 1):
            return 7
        if (plane == 2):
            return 8
    if (mdl == -1):
            
        if (plane == 1):
            return 9 #target 1
        if (plane == 2):
            return 9
    if (mdl == 0):
            
        if (plane == 1):
            return 10
        if (plane == 2):
            return 11
    if (mdl == 1):
            
        if (plane == 1):
            return 12
        if (plane == 2):
            return 13
    if (mdl == 2):
            
        if (plane == 1):
            return 14
        if (plane == 2):
            return 15
    if (mdl == 3):
            
        if (plane == 1):
            return 16
        if (plane == 2):
            return 17
    if (mdl == 4):
            
        if (plane == 1):
            return 18 #Target 2!
    if (mdl == 5):
            
        if (plane == 1):
            return 19
        if (plane == 2):
            return 20
    if (mdl == 6):
            
        if (plane == 1):
            return 21
        if (plane == 2):
            return 22
    
    if (mdl == 7):
            
        if (plane == 1):
            return 23
        if (plane == 2):
            return 24
    if (mdl == 8):
            
        if (plane == 1):
            return 25
        if (plane == 2):
            return 26
    if (mdl == 9):
            
        if (plane == 1):
            return 27 #target 3!
        if (plane == 2):
            return 27
    #water?
    #if (mdl == 10):
            
    #    if (plane == 1):
            
    #        return 1 #weird!!
    #    
    #    if (plane == 2):
            
    #        return 0 #weird!!
    #    
    #

    if (mdl == 11):
            
        if (plane == 1):
            return 28
        if (plane == 2):
            return 29
    if (mdl == 12):
            
        if (plane == 1):
            return 30
        if (plane == 2):
            return 31
    if (mdl == 13):
            
        if (plane == 1):
            return 32
        if (plane == 2):
            return 33
    if (mdl == 14):
            
        if (plane == 1):
            return 34
        if (plane == 2):
            return 35
    #Here goes Water!!!
    if (mdl == -999):
            
        if(plane == -999):
            return 36
    if (mdl == 15):
            
        if (plane == 1):
            return 37
        if (plane == 2):
            return 38
    if (mdl == 16):
            
        if (plane == 1):
            return 39
        if (plane == 2):
            return 40
    if (mdl == 17):
            
        if (plane == 1):
            return 41
        if (plane == 2):
            return 42
    if (mdl == 18):
            
        if (plane == 1):
            return 43
        if (plane == 2):
            return 44
    if (mdl == 19):
            
        if (plane == 1):
            return 45 #Target 4!!
        if (plane == 2):
            return 45 #Target 4!!
    if (mdl == 20):
            
        if (plane == 1):
            return 46
        if (plane == 2):
            return 47
    if (mdl == 21):
            
        if (plane == 1):
            return 48
        if (plane == 2):
            return 49  
    if (mdl == 22):
            
        if (plane == 1):
            return 50 #target 5!!
        if (plane == 2):
            return 50 #target 5!!  
    if (mdl == 23):
            
        if (plane == 1):
            return 51
        if (plane == 2):
            return 52
    if (mdl == 24):
            
        if (plane == 1):
            return 53
        if (plane == 2):
            return 54  
    if (mdl == 25):
            
        if (plane == 1):
            return 55
        if (plane == 2):
            return 56  
    if (mdl == 26):
            
        if (plane == 1):
            return 57
        if (plane == 2):
            return 58  
    if (mdl == 27):
            
        if (plane == 1):
            return 59
        if (plane == 2):
            return 60  
    if (mdl == 28):
            
        if (plane == 1):
            return 61
        if (plane == 2):
            return 62  
    if (mdl == 29):
            
        if (plane == 1):
            return 63
        if (plane == 2):
            return 64  
    if (mdl == 30):
            
        if (plane == 1):
            return 65
        if (plane == 2):
            return 66  
    if (mdl == 31):
            
        if (plane == 1):
            return 67
        if (plane == 2):
            return 68  
    if (mdl == 32):
            
        if (plane == 1):
            return 69
        if (plane == 2):
            return 70  
    if (mdl == 33):
            
        if (plane == 1):
            return 71
        if (plane == 2):
            return 72  
    if (mdl == 34):
            
        if (plane == 1):
            return 73
        if (plane == 2):
            return 74  
    if (mdl == 35):
            
        if (plane == 1):
            return 75
        if (plane == 2):
            return 76  
    if (mdl == 36):
            
        if (plane == 1):
            return 77
        if (plane == 2):
            return 78  
    if (mdl == 37):
            
        if (plane == 1):
            return 79
        if (plane == 2):
            return 80  
    if (mdl == 38):
            
        if (plane == 1):
            return 81
        if (plane == 2):
            return 82  
    if (mdl == 39):
            
        if (plane == 1):
            return 83
        if (plane == 2):
            return 84  
    if (mdl == 40):
            
        if (plane == 1):
            return 85
        if (plane == 2):
            return 86  
    if (mdl == 41):
            
        if (plane == 1):
            return 87
        if (plane == 2):
            return 88  
    if (mdl == 42):
            
        if (plane == 1):
            return 89
        if (plane == 2):
            return 90  
    if (mdl == 43):
            
        if (plane == 1):
            return 91
        if (plane == 2):
            return 92  
    if (mdl == 44):
            
        if (plane == 1):
            return 93
        if (plane == 2):
            return 94  
    if (mdl == 45):
            
        if (plane == 1):
            return 95
        if (plane == 2):
            return 96  
    if (mdl == 46):
            
        if (plane == 1):
            return 97
        if (plane == 2):
            return 98  
    if (mdl == 47):
            
        if (plane == 1):
            return 99
        if (plane == 2):
            return 100  
    if (mdl == 48):
            
        if (plane == 1):
            return 101
        if (plane == 2):
            return 102  
    if (mdl == 49):
            
        if (plane == 1):
            return 103
        if (plane == 2):
            return 104  
    if (mdl == 50):
            
        if (plane == 1):
            return 105
        if (plane == 2):
            return 106  
    if (mdl == 51):
            
        if (plane == 1):
            return 107
        if (plane == 2):
            return 108  
    if (mdl == 52):
            
        if (plane == 1):
            return 109
        if (plane == 2):
            return 110  
    if (mdl == 53):
            
        if (plane == 1):
            return 111
        if (plane == 2):
            return 112  
    if (mdl == 54):
            
        if (plane == 1):
            return 113
        if (plane == 2):
            return 114  
    if (mdl == 55):
            
        if (plane == 1):
            return 115
        if (plane == 2):
            return 116  
    if (mdl == 56):
            
        if (plane == 1):
            return 117
        if (plane == 2):
            return 118  
    if (mdl == 57):
            
        if (plane == 1):
            return 119
        if (plane == 2):
            return 120  
    if (mdl == 58):
            
        if (plane == 1):
            return 121
        if (plane == 2):
            return 122  
    if (mdl == 59):
            
        if (plane == 1):
            return 123
        if (plane == 2):
            return 124  
    if (mdl == 60):
            
        if (plane == 1):
            return 125
        if (plane == 2):
            return 126  
    if (mdl == 61):
            
        if (plane == 1):
            return 127
        if (plane == 2):
            return 128  
    if (mdl == 62):
            
        if (plane == 1):
            return 129
        if (plane == 2):
            return 130  
    if (mdl == 63):
            
        if (plane == 1):
            return 131
        if (plane == 2):
            return 132  
    if (mdl == 64):
            
        if (plane == 1):
            return 133
        if (plane == 2):
            return 134  
    if (mdl == 65):
            
        if (plane == 1):
            return 135
        if (plane == 2):
            return 136  
    if (mdl == 66):
            
        if (plane == 1):
            return 1
        if (plane == 2):
            return 138  
    if (mdl == 67):
            
        if (plane == 1):
            return 139
        if (plane == 2):
            return 140  
    if (mdl == 68):
            
        if (plane == 1):
            return 141
        if (plane == 2):
            return 142  
    if (mdl == 69):
            
        if (plane == 1):
            return 143
        if (plane == 2):
            return 144  
    if (mdl == 70):
            
        if (plane == 1):
            return 145
        if (plane == 2):
            return 146  
    if (mdl == 71):
            
        if (plane == 1):
            return 147
        if (plane == 2):
            return 148  
    if (mdl == 72):
            
        if (plane == 1):
            return 149
        if (plane == 2):
            return 150  
    if (mdl == 73):
            
        if (plane == 1):
            return 151
        if (plane == 2):
            return 152  
    if (mdl == 74):
            
        if (plane == 1):
            return 153
        if (plane == 2):
            return 154  
    if (mdl == 75):
            
        if (plane == 1):
            return 155
        if (plane == 2):
            return 156
    
    if (mdl == 76):
            
        if (plane == 1):
            return 157
        if (plane == 2):
            return 158  
    if (mdl == 77):
            
        if (plane == 1):
            return 159
        if (plane == 2):
            return 160  
    if (mdl == 78):
            
        if (plane == 1):
            return 161
        if (plane == 2):
            return 162  
    if (mdl == 79):
            
        if (plane == 1):
            return 163
        if (plane == 2):
            return 164  
    if (mdl == 80):
            
        if (plane == 1):
            return 165
        if (plane == 2):
            return 166  
    if (mdl == 81):
            
        if (plane == 1):
            return 167
        if (plane == 2):
            return 168  
    if (mdl == 82):
            
        if (plane == 1):
            return 169
        if (plane == 2):
            return 170  
    if (mdl == 83):
            
        if (plane == 1):
            return 171
        if (plane == 2):
            return 172  
    if (mdl == 84):
            
        if (plane == 1):
            return 173
        if (plane == 2):
            return 174  
    if (mdl == 85):
            
        if (plane == 1):
            return 175
        if (plane == 2):
            return 176  
    if (mdl == 86):
            
        if (plane == 1):
            return 177
        if (plane == 2):
            return 178  
    if (mdl == 87):
            
        if (plane == 1):
            return 179
        if (plane == 2):
            return 180  
    if (mdl == 88):
            
        if (plane == 1):
            return 181
        if (plane == 2):
            return 182  
    if (mdl == 89):
            
        if (plane == 1):
            return 183
        if (plane == 2):
            return 184  
    if (mdl == 90):
            
        if (plane == 1):
            return 185
        if (plane == 2):
            return 186  
    if (mdl == 91):
            
        if (plane == 1):
            return 187
        if (plane == 2):
            return 188  
    if (mdl == 92):
            
        if (plane == 1):
            return 189
        if (plane == 2):
            return 190  
    if (mdl == 93):
            
        if (plane == 1):
            return 191
        if (plane == 2):
            return 192  
    if (mdl == 94):
            
        if (plane == 1):
            return 193
        if (plane == 2):
            return 194  
    if (mdl == 95):
            
        if (plane == 1):
            return 195
        if (plane == 2):
            return 195  
    if (mdl == 96):
            
        if (plane == 1):
            return 196
        if (plane == 2):
            return 196  
    if (mdl == 97):
            
        if (plane == 1):
            return 197
        if (plane == 2):
            return 197  
    if (mdl == 98):
            
        if (plane == 1):
            return 198 #not needed
        if (plane == 2):
            return 198  
    if (mdl == 99):
            
        if (plane == 1):
            return 199 #not needed
        if (plane == 2):
            return 199  
    if (mdl == 100):
            
        if (plane == 2):
            return 200  
    if (mdl == 101):
            
        if (plane == 2):
            return 201  
    if (mdl == 102):
            
        if (plane == 2):
            return 202  
    if (mdl == 103):
            
        if (plane == 2):
            return 203  
    if (mdl == 104):
            
        if (plane == 2):
            return 204  
    if (mdl == 105):
            
        if (plane == 2):
            return 205  
    if (mdl == 106):
            
        if (plane == 2):
            return 206  
    if (mdl == 107):
            
        if (plane == 2):
            return 207  
    if (mdl == 108):
            
        if (plane == 2):
            return 208  
    if (mdl == 109):
            
        if (plane == 2):
            return 209  
    if (mdl == 110):
            
        if (plane == 2):
            return 210  
    if (mdl == 111):
            
        if (plane == 2):
            return 211  
    if (mdl == 112):
            
        if (plane == 2):
            return 212  
    if (mdl == 113):
            
        if (plane == 2):
            return 213  
    if (mdl == 114):
            
        if (plane == 2):
            return 214  
    if (mdl == 115):
            
        if (plane == 2):
            return 2  
    if (mdl == 116):
            
        if (plane == 2):
            return 0  
    if (mdl == 117):
            
        if (plane == 2):
            return 1  
    if (mdl == 118):
            
        if (plane == 2):
            return 0  
    if (mdl == 119):
            
        if (plane == 2):
            return 2  
    


def CountPOT(fileName):
  meta = PlotUtils.ChainWrapper("Meta")
  meta.Add(fileName)

  nEntries = meta.GetEntries()
  pot_used = 0
  for entry in range(nEntries):
      pot_used += meta.GetValue("POT_Used", entry)
  return pot_used



if __name__ == '__main__':

    parser = argparse.ArgumentParser( prog='waterVertexInvestigations',
                    description='Produces several plots relating to the water target ANN vertexing')
    parser.add_argument('datapath', nargs=1, help="the data input file")
    parser.add_argument('mcpath', nargs=1, help="the data input file")
    parser.add_argument('outpath', nargs=1, help="the output file basepath")
    args = parser.parse_args()
    
        
    datapath=args.datapath[0]
    mcpath=args.mcpath[0]
    outpath=args.outpath[0]


    emptyplaylists = ["1A","1B", "1C", "1D", "1E", "1F", "1G"]
    filledplaylists = ["1L","1M", "1N", "1O", "1P"]

    recochw = PlotUtils.ChainWrapper("MasterAnaDev")
    mcchw = PlotUtils.ChainWrapper("MasterAnaDev")     
    recochwfilled = PlotUtils.ChainWrapper("MasterAnaDev")
    mcchwfilled = PlotUtils.ChainWrapper("MasterAnaDev")     
    for empty in emptyplaylists:
        datapathstr = datapath+"/"+empty+"-Data.txt"
        mcpathstr = mcpath+"/"+empty+"-MC.txt"
        recochw.Add(datapathstr)
        mcchw.Add(mcpathstr)
    for filled in filledplaylists:
        datapathstr = datapath+"/"+filled+"-Data.txt"
        mcpathstr = mcpath+"/"+filled+"-MC.txt"
        recochwfilled.Add(datapathstr)
        mcchwfilled.Add(mcpathstr)


    nDataEntries = recochw.GetEntries()
    nMCEntries = mcchw.GetEntries()
    nDataEntriesFilled = recochwfilled.GetEntries()
    nMCEntriesFilled = mcchwfilled.GetEntries()
    print("nDataEntries", nDataEntries)
    print("nMCEntries", nMCEntries)
    print("nDataEntriesFilled", nDataEntriesFilled)
    print("nMCEntriesFilled", nMCEntriesFilled)
    nbins = 10
    tgtUtil=PlotUtils.TargetUtils()

    """ trueZ = ROOT.MnvH1D("trueWaterZ", "trueWaterZ", 80, 0, 80)
    recoZ = ROOT.MnvH1D("recoZ", "Reconstructed Water Z", 80, 0, 80)
    recoInThisZTrueInThisZ = ROOT.MnvH1D("recoInThisZTrueInThisZ", "Reconstructed Water Z", 80, 0, 80)
    trueInThisZRecoInThisZ = ROOT.MnvH1D("trueInThisZRecoInThisZ", "trueInThisZRecoInThisZ", 80, 0, 80) """


    vertexBins = [4293.04, 4337.25, 4381.47, 4425.68, 4514.11, 4558.33, 4602.54, 4646.76, 4735.19, 4779.4, 4823.62, 4867.83, 5000.48, 5044.69, 5088.91, 5133.12, 5200, 5420, 5456.74, 5500.95, 5545.17, 5589.38, 5677.81, 5722.03, 5810.45, 5855.68, 5900.91, 5946.14, 5991.37, 6036.6, 6081.83, 6127.06, 6172.29, 6217.52, 6262.74, 6307.97, 6353.2, 6398.43, 6443.66, 6488.89, 6534.12, 6579.35, 6624.58, 6669.81, 6715.03, 6760.26, 6805.49, 6850.72, 6895.95, 6941.18, 6986.41, 7031.64, 7076.87, 7122.1, 7167.32, 7212.55, 7257.78, 7303.01, 7348.24, 7393.47, 7438.7, 7483.93, 7529.16, 7574.39, 7619.61, 7664.84, 7710.07, 7755.3, 7800.53, 7845.76, 7890.99, 7936.22, 7981.45, 8026.68, 8071.9, 8117.13, 8162.36, 8207.59, 8252.82, 8298.05, 8343.28, 8388.51, 8433.74, 8478.97, 8524.19, 8569.42, 8614.65]
    print("len(vertexBins): ", len(vertexBins))
    binsC = (ctypes.c_double*len(vertexBins))(*vertexBins)
    binsPointer = ctypes.cast(binsC, ctypes.POINTER(ctypes.c_double))

    """ trueZ = ROOT.MnvH1D("trueWaterZ", "trueWaterZ", len(vertexBins)-1, binsPointer)
    recoZ = ROOT.MnvH1D("recoZ", "Reconstructed Water Z", len(vertexBins)-1, binsPointer)
    recoInThisZTrueInThisZ = ROOT.MnvH1D("recoInThisZTrueInThisZ", "Reconstructed Water Z", len(vertexBins)-1, binsPointer)
    trueInThisZRecoInThisZ = ROOT.MnvH1D("trueInThisZRecoInThisZ", "trueInThisZRecoInThisZ", len(vertexBins)-1, binsPointer)
    trueInWaterRecoWhere = ROOT.MnvH1D("trueInWaterRecoWhere", "trueInWaterRecoWhere", len(vertexBins)-1, binsPointer)
    recoInWaterTrueWhere = ROOT.MnvH1D("recoInWaterTrueWhere", "recoInWaterTrueWhere", len(vertexBins)-1, binsPointer)

    vertexingConfusionMatrix = ROOT.MnvH2D("vertexingConfusionMatrix", "vertexingConfusionMatrix", len(vertexBins)-1, binsPointer, len(vertexBins)-1, binsPointer)
    """

    nModulesToPlot = 90
    trueZ = ROOT.MnvH1D("trueWaterZ", "trueWaterZ", nModulesToPlot, -.5, nModulesToPlot-.5)
    recoZ = ROOT.MnvH1D("recoZ", "Reconstructed Water Z", nModulesToPlot, -.5, nModulesToPlot-.5)
    recoInThisZTrueInThisZ = ROOT.MnvH1D("recoInThisZTrueInThisZ", "Reconstructed Water Z", nModulesToPlot, -.5, nModulesToPlot-.5)
    trueInThisZRecoInThisZ = ROOT.MnvH1D("trueInThisZRecoInThisZ", "trueInThisZRecoInThisZ", nModulesToPlot, -.5, nModulesToPlot-.5)
    trueInWaterRecoWhere = ROOT.MnvH1D("trueInWaterRecoWhere", "trueInWaterRecoWhere", nModulesToPlot, -.5, nModulesToPlot-.5)
    recoInWaterTrueWhere = ROOT.MnvH1D("recoInWaterTrueWhere", "recoInWaterTrueWhere", nModulesToPlot, -.5, nModulesToPlot-.5)

    vertexingConfusionMatrix = ROOT.MnvH2D("vertexingConfusionMatrix", "vertexingConfusionMatrix", nModulesToPlot, -.5, nModulesToPlot-.5, nModulesToPlot, -.5, nModulesToPlot-.5)


    trueWaterRecoInWater = ROOT.MnvH1D("trueWaterRecoInWater", "Reconstructed In Water", nbins, PlotUtils.TargetProp.WaterTarget.Face, PlotUtils.TargetProp.WaterTarget.Back)
    trueWaterRecoDS = ROOT.MnvH1D("trueWaterRecoDS", "Reconstructed Downstream", nbins, PlotUtils.TargetProp.WaterTarget.Face, PlotUtils.TargetProp.WaterTarget.Back)
    trueWaterRecoUS = ROOT.MnvH1D("trueWaterRecoUS", "Reconstructed Upstream", nbins, PlotUtils.TargetProp.WaterTarget.Face, PlotUtils.TargetProp.WaterTarget.Back)
    trueWater = ROOT.MnvH1D("trueWater", "trueWater", nbins, PlotUtils.TargetProp.WaterTarget.Face, PlotUtils.TargetProp.WaterTarget.Back)

    recoWaterTrueInWater = ROOT.MnvH1D("recoWaterTrueInWater", "Truth In Water", nbins, PlotUtils.TargetProp.WaterTarget.Face, PlotUtils.TargetProp.WaterTarget.Back)
    recoWaterTrueDS = ROOT.MnvH1D("recoWaterTrueDS", "Truth Downstream", nbins, PlotUtils.TargetProp.WaterTarget.Face, PlotUtils.TargetProp.WaterTarget.Back)
    recoWaterTrueUS = ROOT.MnvH1D("recoWaterTrueUS", "Truth Upstream", nbins, PlotUtils.TargetProp.WaterTarget.Face, PlotUtils.TargetProp.WaterTarget.Back)
    recoWater = ROOT.MnvH1D("recoWater", "recoWater", nbins, PlotUtils.TargetProp.WaterTarget.Face, PlotUtils.TargetProp.WaterTarget.Back)


    trueWaterRecoInWaterApothem = ROOT.MnvH1D("trueWaterRecoInWaterApothem", "Reconstructed In Water Apothem", nbins, 0, 900)
    trueWaterRecoFalseByApothem = ROOT.MnvH1D("######", "#####", nbins, 0, 900)
    trueWaterApothem = ROOT.MnvH1D("trueWaterApothem", "trueWater", nbins, 0, 900)

    recoWaterTrueInWaterApothem = ROOT.MnvH1D("recoWaterTrueInWaterApothem", "Truth In Water Apothem", nbins, 0, 900)
    recoWaterTrueFalseByApothem = ROOT.MnvH1D("######", "#######", nbins, 0, 900)
    recoWaterApothem = ROOT.MnvH1D("recoWaterApothem", "recoWater", nbins, 0, 900)



    rebinNum = 10
    nbins = 3400
    xlow = 4200
    #xhigh = 5900 # Zoom into target region
    xhigh = 8410 # Full detector


    TruthVerticesMCHighResEmpty = ROOT.TH1D("TruthVerticesMCHighResEmpty", "TruthVerticesMCHighResEmpty", nbins, xlow, xhigh)
    TBVerticesMCHighResEmpty = ROOT.TH1D("TBVerticesMCHighResEmpty", "TBVerticesMCHighResEmpty", nbins, xlow, xhigh)
    ANNVerticesMCHighResEmpty = ROOT.TH1D("ANNVerticesMCHighResEmpty", "ANNVerticesMCHighResEmpty", nbins, xlow, xhigh)
    TBVerticesDataHighResEmpty = ROOT.TH1D("TBVerticesDataHighResEmpty", "TBVerticesDataHighResEmpty", nbins, xlow, xhigh)
    ANNVerticesDataHighResEmpty = ROOT.TH1D("ANNVerticesDataHighResEmpty", "ANNVerticesDataHighResEmpty", nbins, xlow, xhigh)

    TruthVerticesMCHighResFilled = ROOT.TH1D("TruthVerticesMCHighResFilled", "TruthVerticesMCHighResFilled", nbins, xlow, xhigh)
    TBVerticesMCHighResFilled = ROOT.TH1D("TBVerticesMCHighResFilled", "TBVerticesMCHighResFilled", nbins, xlow, xhigh)
    ANNVerticesMCHighResFilled = ROOT.TH1D("ANNVerticesMCHighResFilled", "ANNVerticesMCHighResFilled", nbins, xlow, xhigh)
    TBVerticesDataHighResFilled = ROOT.TH1D("TBVerticesDataHighResFilled", "TBVerticesDataHighResFilled", nbins, xlow, xhigh)
    ANNVerticesDataHighResFilled = ROOT.TH1D("ANNVerticesDataHighResFilled", "ANNVerticesDataHighResFilled", nbins, xlow, xhigh)


    for i in range(nMCEntries):
        if (i%100000==0):
            print(i,"/",nMCEntries, flush=True)
        annvtxsize=mcchw.GetValue("ANN_vtx_sz",i, 0)
        mcvtxx = mcchw.GetValue("vtx",i, 0)
        mcvtxy = mcchw.GetValue("vtx",i, 1)
        mcvtxz = mcchw.GetValue("vtx",i, 2)
        mc_plane = mcchw.GetValue("truth_vtx_plane",i, 0)
        mc_module = mcchw.GetValue("truth_vtx_module",i, 0)
        tbvtxz = mcchw.GetValue("MasterAnaDev_vtx", i, 2)
        annsegment = int(mcchw.GetValue("ANN_segments",i, 0))
        mcsegment = planecode(mc_module, mc_plane)
        annProbs = mcchw.GetValue("ANN_plane_probs",i, 0)
        mc_tgtid = mcchw.GetValue("truth_targetID",i, 0)


        TruthVerticesMCHighResEmpty.Fill(mcvtxz)
        TBVerticesMCHighResEmpty.Fill(tbvtxz)
        if annvtxsize>=3:
            annvtxz = mcchw.GetValue("ANN_vtx",i, 2)
            ANNVerticesMCHighResEmpty.Fill(annvtxz)


    for i in range(nDataEntries):
        if (i%100000==0):
            print(i,"/",nDataEntries, flush=True)
        annvtxsize = recochw.GetValue("ANN_vtx_sz",i, 0)
        tbvtxz = recochw.GetValue("MasterAnaDev_vtx", i, 2)
        annProbs = recochw.GetValue("ANN_plane_probs",i, 0)

        TBVerticesMCHighResEmpty.Fill(tbvtxz)
        if annvtxsize>=3:
            annvtxz = recochw.GetValue("ANN_vtx",i, 2)
            ANNVerticesMCHighResEmpty.Fill(annvtxz)

    #MC loop
    mcsum = 0
    for i in range(nMCEntries):        ##Add a progress update?
        if (i%100000==0):
            print(i,"/",nMCEntries, flush=True)
        annvtxsize=mcchw.GetValue("ANN_vtx_sz",i, 0)
        mcvtxx = mcchw.GetValue("vtx",i, 0)
        mcvtxy = mcchw.GetValue("vtx",i, 1)
        mcvtxz = mcchw.GetValue("vtx",i, 2)
        mc_plane = mcchw.GetValue("truth_vtx_plane",i, 0)
        mc_module = mcchw.GetValue("truth_vtx_module",i, 0)
        tbvtxz = mcchw.GetValue("MasterAnaDev_vtx", i, 2)
        annsegment = int(mcchw.GetValue("ANN_segments",i, 0))
        mcsegment = planecode(mc_module, mc_plane)
        annProbs = mcchw.GetValue("ANN_plane_probs",i, 0)
        mc_tgtid = mcchw.GetValue("truth_targetID",i, 0)

        if mc_tgtid == 1:
            mcsegment = 9
        elif mc_tgtid == 2:
            mcsegment = 18
        elif mc_tgtid == 3:
            mcsegment = 27
        elif mc_tgtid == 4:
            mcsegment = 45
        elif mc_tgtid == 5:
            mcsegment = 50
        #print("truebin: ",truebin)
        if mcsegment:
            truebin = trueZ.Fill(mcsegment)

        if annvtxsize>=3:
            annvtxx = mcchw.GetValue("ANN_vtx",i, 0)
            annvtxy = mcchw.GetValue("ANN_vtx",i, 1)
            annvtxz = mcchw.GetValue("ANN_vtx",i, 2)

            waterMidZ=(PlotUtils.TargetProp.WaterTarget.Face+PlotUtils.TargetProp.WaterTarget.Back)/2
            waterBin = trueZ.FindBin(waterMidZ)
            #if (trueZ.FindBin(annvtxz) == waterBin): #Debugging, if reconstructed in water
                #print("Water bin: ", waterBin, ", true bin: ", truebin)


            #trueZ.Fill(mcvtxz, 1)
            recobin = recoZ.Fill(annsegment)
            if (truebin==recobin):
                recoInThisZTrueInThisZ.Fill(annsegment, 1)
                trueInThisZRecoInThisZ.Fill(mcsegment, 1)


            vertexingConfusionMatrix.Fill(mcsegment, annsegment)

            if (truebin==17):
                trueInWaterRecoWhere.Fill(annvtxz)
            if (recobin==17):
                recoInWaterTrueWhere.Fill(mcvtxz)
            if (mcvtxz>7500 and recobin==17):
                print("What the helly:", i)
                print("mc_run:", mcchw.GetValue("mc_run",i, 0))
                print("mc_subrun:", mcchw.GetValue("mc_subrun",i, 0))

            annapothem = math.sqrt(annvtxx**2 + annvtxy**2)
            mcapothem = math.sqrt(mcvtxx**2 + mcvtxy**2)
            ##print("annvtxz: ", annvtxz)
            if tgtUtil.InWaterTargetVolMC(mcvtxx,mcvtxy,mcvtxz): #For each event truly in the water, where was it reconstructed
                trueWater.Fill(mcvtxz,1)
                trueWaterApothem.Fill(mcapothem,1)
                if tgtUtil.InWaterTargetVolMC(annvtxx,annvtxy,annvtxz): #If it was truly in water and reconstructed in water z
                    trueWaterRecoInWater.Fill(mcvtxz, 1)
                elif annvtxz> PlotUtils.TargetProp.WaterTarget.Back: #If it was truly in water and reconstructed in water z
                    trueWaterRecoDS.Fill(mcvtxz, 1)
                elif annvtxz< PlotUtils.TargetProp.WaterTarget.Face: #If it was truly in water and reconstructed in water z
                    trueWaterRecoUS.Fill(mcvtxz, 1)
                ##Else something else idk

                if tgtUtil.InWaterTargetVolMC(annvtxx,annvtxy,annvtxz): #If it was truly in water and reconstructed in water z
                    trueWaterRecoInWaterApothem.Fill(mcapothem, 1)
                else: #If it was truly in water and reconstructed in water z
                    trueWaterRecoFalseByApothem.Fill(mcapothem, 1)

            if tgtUtil.InWaterTargetZMC(annvtxz): #For each event truly in the water, where was it reconstructed
                recoWater.Fill(annvtxz,1)
                recoWaterApothem.Fill(annapothem,1)
                if tgtUtil.InWaterTargetZMC(mcvtxz): #If it was truly in water and reconstructed in water
                    recoWaterTrueInWater.Fill(annvtxz, 1)
                elif mcvtxz> PlotUtils.TargetProp.WaterTarget.Back: #If it was truly in water and reconstructed in water
                    recoWaterTrueDS.Fill(annvtxz, 1)
                elif mcvtxz< PlotUtils.TargetProp.WaterTarget.Face: #If it was truly in water and reconstructed in water
                    recoWaterTrueUS.Fill(annvtxz, 1)

                if tgtUtil.InWaterTargetVolMC(mcvtxx,mcvtxy,mcvtxz): #If it was truly in water and reconstructed in water
                    recoWaterTrueInWaterApothem.Fill(annapothem, 1)
                else: #If it was truly in water and reconstructed in water
                    recoWaterTrueFalseByApothem.Fill(annapothem, 1)


    recall =  trueWaterRecoInWater.GetSumOfWeights()/trueWater.GetSumOfWeights()* 100
    accuracy =  recoWaterTrueInWater.GetSumOfWeights()/recoWater.GetSumOfWeights()* 100
    recallapothem =  trueWaterRecoInWaterApothem.GetSumOfWeights()/trueWater.GetSumOfWeights()* 100
    accuracyapothem =  recoWaterTrueInWaterApothem.GetSumOfWeights()/recoWater.GetSumOfWeights()* 100

    trueWaterRecoInWater.Divide(trueWaterRecoInWater,trueWater)
    trueWaterRecoDS.Divide(trueWaterRecoDS, trueWater)
    trueWaterRecoUS.Divide(trueWaterRecoUS,trueWater)
    #trueWater.Divide(trueWater, trueWater)

    recoWaterTrueInWater.Divide(recoWaterTrueInWater,recoWater)
    recoWaterTrueDS.Divide(recoWaterTrueDS, recoWater)
    recoWaterTrueUS.Divide(recoWaterTrueUS,recoWater)
    #recoWater.Divide(recoWater, recoWater)

    canvas1 = ROOT.TCanvas() # have to declare canvas before calling mnvplotter :))
    mnv = PlotUtils.MnvPlotter()
    mnv.legend_text_size=0.02

    objarray = ROOT.TObjArray()

    recallstr = "Overall Reco Recall: " + '%.4g' % recall + "%"
    objarray.Add(trueWaterRecoInWater)
    objarray.Add(trueWaterRecoDS)
    objarray.Add(trueWaterRecoUS)
    mnv.DrawStackedMC(objarray, 1, "TR")
    mnv.AddPlotLabel(recallstr, 0.85, 0.95, 0.03)
    canvas1.Print(outpath+"trueWaterRecoWhere.png")

    objarray2 = ROOT.TObjArray()
    accuracystr = "Overall Reco Accuracy: " + '%.4g' % accuracy + "%"
    objarray2.Add(recoWaterTrueInWater)
    objarray2.Add(recoWaterTrueDS)
    objarray2.Add(recoWaterTrueUS)
    mnv.DrawStackedMC(objarray2, 1, "TR")
    mnv.AddPlotLabel(accuracystr, 0.85, 0.95, 0.03)
    canvas1.Print(outpath+"recoWaterTrueWhere.png")



    trueWaterRecoInWaterApothem.Divide(trueWaterRecoInWaterApothem,trueWaterApothem)
    trueWaterRecoFalseByApothem.Divide(trueWaterRecoFalseByApothem, trueWaterApothem)
    #trueWater.Divide(trueWater, trueWater)

    recoWaterTrueInWaterApothem.Divide(recoWaterTrueInWaterApothem,recoWaterApothem)
    recoWaterTrueFalseByApothem.Divide(recoWaterTrueFalseByApothem, recoWaterApothem)
    #recoWater.Divide(recoWater, recoWater)

    objarray3 = ROOT.TObjArray()
    recallstr = "Overall Reco Recall: " + '%.4g' % recallapothem + "%"
    objarray3.Add(trueWaterRecoInWaterApothem)
    objarray3.Add(trueWaterRecoFalseByApothem)
    mnv.DrawStackedMC(objarray3, 1, "TR")
    mnv.AddPlotLabel(recallstr, 0.85, 0.97, 0.03)
    canvas1.Print(outpath+"trueWaterRecoWhereApothem.png")

    objarray4 = ROOT.TObjArray()
    accuracystr = "Overall Reco Accuracy: " + '%.4g' % accuracyapothem + "%"
    objarray4.Add(recoWaterTrueInWaterApothem)
    objarray4.Add(recoWaterTrueFalseByApothem)
    mnv.DrawStackedMC(objarray4, 1, "TR")
    mnv.AddPlotLabel(accuracystr, 0.85, 0.97, 0.03)
    canvas1.Print(outpath+"recoWaterTrueWhereApothem.png")

    canvas1.Clear()
    
    #print(trueZ.GetSumOfWeights())

    #Calculate recall and accuracy
    #Overall
    recall =  trueInThisZRecoInThisZ.GetSumOfWeights()/trueZ.GetSumOfWeights()* 100
    accuracy =  recoInThisZTrueInThisZ.GetSumOfWeights()/recoZ.GetSumOfWeights()* 100
    accuracystr = "Overall Reco Accuracy: " + '%.4g' % accuracy + "%"
    recallstr = "Overall Reco Recall: " + '%.4g' % recall + "%"
    #Specifically in nuclear targets
    #Get nuclear target bins
    tgt1bin = trueZ.FindBin(PlotUtils.TargetProp.CenterZMC.Tgt1.center_nooffset)
    tgt2bin = trueZ.FindBin(PlotUtils.TargetProp.CenterZMC.Tgt2.center_nooffset)
    tgt3bin = trueZ.FindBin(PlotUtils.TargetProp.CenterZMC.Tgt3.center_nooffset)
    tgt4bin = trueZ.FindBin(PlotUtils.TargetProp.CenterZMC.Tgt4.center_nooffset)
    tgt5bin = trueZ.FindBin(PlotUtils.TargetProp.CenterZMC.Tgt5.center_nooffset)
    waterMidZ=(PlotUtils.TargetProp.WaterTarget.Face+PlotUtils.TargetProp.WaterTarget.Back)/2
    waterBin = trueZ.FindBin(waterMidZ)
    #Get appropriate values
    #Get the total number of events reconstructed in each target
    recoInTgt1 = recoZ.GetBinContent(tgt1bin)
    recoInTgt2 = recoZ.GetBinContent(tgt2bin)
    recoInTgt3 = recoZ.GetBinContent(tgt3bin)
    recoInTgt4 = recoZ.GetBinContent(tgt4bin)
    recoInTgt5 = recoZ.GetBinContent(tgt5bin)
    recoInWater = recoZ.GetBinContent(waterBin)
    recoInTargetsSum = recoInTgt1 + recoInTgt2 +recoInTgt3 +recoInTgt4 + recoInTgt5 + recoInWater
    #Get the total number of events truly in each target
    trueInTgt1 = trueZ.GetBinContent(tgt1bin)
    trueInTgt2 = trueZ.GetBinContent(tgt2bin)
    trueInTgt3 = trueZ.GetBinContent(tgt3bin)
    trueInTgt4 = trueZ.GetBinContent(tgt4bin)
    trueInTgt5 = trueZ.GetBinContent(tgt5bin)
    trueInWater = trueZ.GetBinContent(waterBin)
    trueInTargetsSum = trueInTgt1 + trueInTgt2 + trueInTgt3 + trueInTgt4 + trueInTgt4 + trueInWater
    #Get the number of events reconstructed in each target and truly in that target
    recoInAndTrueInTgt1 = recoInThisZTrueInThisZ.GetBinContent(tgt1bin)
    recoInAndTrueInTgt2 = recoInThisZTrueInThisZ.GetBinContent(tgt2bin)
    recoInAndTrueInTgt3 = recoInThisZTrueInThisZ.GetBinContent(tgt3bin)
    recoInAndTrueInTgt4 = recoInThisZTrueInThisZ.GetBinContent(tgt4bin)
    recoInAndTrueInTgt5 = recoInThisZTrueInThisZ.GetBinContent(tgt5bin)
    recoInAndTrueInWater = recoInThisZTrueInThisZ.GetBinContent(waterBin)
    recoInAndTrueInTargetsSum = recoInAndTrueInTgt1 + recoInAndTrueInTgt2 +recoInAndTrueInTgt3 +recoInAndTrueInTgt4 + recoInAndTrueInTgt5 + recoInAndTrueInWater

    #Get the number of events truly in each target and also reconstructed there
    trueInAndRecoInTgt1 = trueInThisZRecoInThisZ.GetBinContent(tgt1bin)
    trueInAndRecoInTgt2 = trueInThisZRecoInThisZ.GetBinContent(tgt2bin)
    trueInAndRecoInTgt3 = trueInThisZRecoInThisZ.GetBinContent(tgt3bin)
    trueInAndRecoInTgt4 = trueInThisZRecoInThisZ.GetBinContent(tgt4bin)
    trueInAndRecoInTgt5 = trueInThisZRecoInThisZ.GetBinContent(tgt5bin)
    trueInAndRecoInWater = trueInThisZRecoInThisZ.GetBinContent(waterBin)
    trueInAndRecoInTargetsSum = trueInAndRecoInTgt1 + trueInAndRecoInTgt2 +trueInAndRecoInTgt3 +trueInAndRecoInTgt4 + trueInAndRecoInTgt5 + trueInAndRecoInWater
    
    #Calculate total targets recall and accuracy
    recalltargets =  trueInAndRecoInTargetsSum/trueInTargetsSum * 100
    accuracytargets =  recoInAndTrueInTargetsSum/recoInTargetsSum * 100
    accuracytargetsstr = "Targets Reco Accuracy: " + '%.4g' % accuracytargets + "%"
    recalltargetsstr = "Targets Reco Recall: " + '%.4g' % recalltargets + "%"



    recall =  trueInThisZRecoInThisZ.GetSumOfWeights()/trueZ.GetSumOfWeights()* 100
    accuracy =  recoInThisZTrueInThisZ.GetSumOfWeights()/recoZ.GetSumOfWeights()* 100
    print("================================")
    print("========OVERALL ACCURACY========")
    print(recoInThisZTrueInThisZ.GetSumOfWeights(), "/",recoZ.GetSumOfWeights())
    print("=========OVERALL RECALL=========")
    print(trueInThisZRecoInThisZ.GetSumOfWeights(), "/",trueZ.GetSumOfWeights())
    print("================================")
    print("=========WATER ACCURACY=========")
    print(recoInAndTrueInWater, "/",recoInWater)
    print("=========WATER RECALL===========")
    print(trueInAndRecoInWater, "/",trueInWater)
    print("================================")

    mnv.DrawDataMCRatio(recoInThisZTrueInThisZ, recoZ)
    mnv.AddHistoTitle ("Machine learning vertexing accuracy", 0.05)
    PlotRegions(True)
    mnv.AddPlotLabel(accuracystr, 0.9, 0.95, 0.02)
    mnv.AddPlotLabel(accuracytargetsstr, 0.9, 0.92, 0.02)
    canvas1.Print(outpath+"test123.png")
    mnv.DrawDataMCRatio(trueInThisZRecoInThisZ, trueZ)
    mnv.AddHistoTitle ("Machine learning vertexing recall", 0.05)
    PlotRegions(True)
    mnv.AddPlotLabel(recallstr, 0.9, 0.95, 0.02)
    mnv.AddPlotLabel(recalltargetsstr, 0.9, 0.92, 0.02)
    canvas1.Print(outpath+"test1232.png")

    mnv.DrawDataMC(recoZ, trueZ)
    PlotRegions(True)
    canvas1.Print(outpath+"test1233.png")

    vertexingConfusionMatrix.Draw("COLZ")
    vertexingConfusionMatrix.GetXaxis().SetTitle("True vertex module/cm")
    vertexingConfusionMatrix.GetYaxis().SetTitle("Reconstruction vertex module/cm")
    PlotRegions(True, inYmax=nModulesToPlot)
    canvas1.Print(outpath+"test1234.png")

    trueInWaterRecoWhere.Draw("HIST")
    PlotRegions(True)
    canvas1.Print(outpath+"test1235.png")
    recoInWaterTrueWhere.Draw("HIST")
    PlotRegions(True)
    canvas1.Print(outpath+"test1236.png")

