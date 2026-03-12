import ROOT
import array
import os,sys


def getMuonMomentum(mytree):
    muon_mom = ROOT.TVector3()

    nfsp = mytree.nfsp
    Efsp = mytree.E
    pdg  = mytree.pdg
    px = mytree.px
    py = mytree.py
    pz = mytree.pz
    
    for p in range(0,nfsp):
        if(pdg[p]==13):
            muon_mom.SetX(px[p])
            muon_mom.SetY(py[p])
            muon_mom.SetZ(pz[p])
            break
    return muon_mom

def isSignal(mytree):

    isCC = mytree.cc
    isnumu = mytree.PDGnu==14
    mymuon = getMuonMomentum(mytree)
    isTheta = mymuon.Theta()*180/3.1415 < 20
    #print isCC,isnumu,isTheta
    if(isCC and isnumu and isTheta): return True
    return False


inputrw = ROOT.TFile("nuwro_flux_rw.root")
rwhist = inputrw.Get("flux2")
inputfiles = os.popen("ls %s/*%s*.root"%(sys.argv[1],sys.argv[2])).readlines()
mytree = ROOT.TChain("FlatTree_VARS")
for f in inputfiles:
    fn = f.rstrip("\n")
    print fn
    mytree.Add(fn)


#histogram

ptbins = [0,0.075,0.15,0.25,0.325,0.4,0.475,0.55,0.7,0.85,1,1.25,1.5,2.5,4.5]
pzbins = [1.5,2.,2.5,3.,3.5,4.,4.5,5.,6.,7.,8.,9.,10.,15.,20.,40.,60.]

myptpz = ROOT.TH2D("ptpz","ptpz",len(pzbins)-1, array.array("d",pzbins),len(ptbins)-1,array.array("d",ptbins))
count = 0
for e in mytree:
    count +=1
    if(count%50000==0): print count
    if not isSignal(e): continue
    
    fScaleFactor = e.fScaleFactor/len(inputfiles)
    muonmom = getMuonMomentum(e)
    pl = muonmom.Z()
    pt = muonmom.Perp()
    extra = rwhist.GetBinContent(rwhist.FindBin(e.Enu_true))
    myptpz.Fill(pl,pt,fScaleFactor*extra)






myptpz.Draw("COLZ")

myout = ROOT.TFile("NuWro_%s.root"%(sys.argv[2]),"RECREATE")
myptpz.Write("nuwro_%s"%(sys.argv[2]))
myout.Close()

raw_input("done")
    
    