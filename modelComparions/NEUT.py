import ROOT
import array
import sys

def isCCQELike(mytree):

    #need 1 muon and no mesons and photons > 10 MeV in final state

    nfsp = mytree.nfsp
    Efsp = mytree.E
    pdg  = mytree.pdg

    n_muon   = 0
    n_meson  = 0
    n_baryon = 0
    n_gamma  = 0

    badMesons = [211,321,323,111,130,310,311,313,411,421]
    badBaryons = [3112,3122,3212,3222,4112,4122,4212,4222]

    for p in range(0,nfsp):
        if(pdg[p]==13): n_muon+=1
        if(abs(pdg[p]) in badMesons): n_meson+=1
        if(abs(pdg[p]) in badBaryons): n_baryon+=1
        if(pdg[p]==22 and Efsp[p] > 0.010): n_gamma+=1

    if n_muon == 1 and n_meson==0 and n_gamma==0 and n_baryon==0 : return True
    else: return False

def getProtonMomentum(mytree):
    nfsp = mytree.nfsp
    Efsp = mytree.E
    pdg  = mytree.pdg
    px = mytree.px
    py = mytree.py
    pz = mytree.pz
    
    prot_mom = ROOT.TVector3()
    prot_e = 0

    for p in range(0,nfsp):
        if(pdg[p]==2212 and Efsp[p]>prot_e):#highest energy proton
            prot_mom.SetX(px[p])
            prot_mom.SetY(py[p])
            prot_mom.SetZ(pz[p])
    return prot_mom

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

def isInclusiveFHC(mytree):
    nu_pdg = mytree.PDGnu
    lep_pdg = mytree.PDGLep

    if(nu_pdg == 14 and lep_pdg == 13):
        return True
    
    return False

def isTKI(mytree):
    #<70 degree proton
    #<20 degree muon
    #1.5 to 10 muon mom
    #0.45 to 1.2 proton mom

    muonmom = getMuonMomentum(mytree)
    protonmom = getProtonMomentum(mytree)

    goodProtonMom = protonmom.Mag()>0.45 and protonmom.Mag()<1.2
    goodMuonMom   = muonmom.Mag()>1.5 and muonmom.Mag()<10

    goodProtonAngle = protonmom.Theta()*180/3.1415
    goodMuonAngle = muonmom.Theta()*180/3.1415

    if goodProtonMom and goodMuonMom and goodProtonAngle and goodMuonAngle: return True
    
    return False

def getPn(e):

    promom = getProtonMomentum(e)
    mumom = getMuonMomentum(e)

    protonPtVect = getProtonMomentum(e)
    muonPtVect = getMuonMomentum(e)

    protonPtVect.SetZ(0)
    muonPtVect.SetZ(0)


    pn = -99;
    Mn = 939.5654133/1000;
    Mp = 938.2720813/1000;
    Mmu = 105.6583745/1000;
    MA = 6*Mn + 6*Mp - 92.162/1000;
    Bin=27.13/1000;
    MAstar = MA - Mn + Bin;
    
    Epprim = ROOT.TMath.Sqrt(Mp*Mp + promom.Mag2());
    kprimL = mumom.Z();
    pprimL = promom.Z();
    Eprim = ROOT.TMath.Sqrt(Mmu*Mmu+mumom.Mag2());
    
    factor = MA - Eprim - Epprim + kprimL + pprimL;
    pT = (protonPtVect+muonPtVect).Mag();
    pL = -(MAstar*MAstar + pT*pT-factor*factor)/2.0/factor;
    

    pn = ROOT.TMath.Sqrt(pL*pL + pT*pT);
    return pn

inputFile = ROOT.TFile(sys.argv[1])

mytree = inputFile.Get("FlatTree_VARS")



ptbins = [0,0.075,0.15,0.25,0.325,0.4,0.475,0.55,0.7,0.85,1,1.25,1.5,2.5,4.5]
pzbins = [1.5,2.,2.5,3.,3.5,4.,4.5,5.,6.,7.0,8.,9.0,10.,15.,20.,40.,60.]

q3bins = [0,0.2,0.3,0.4,0.6,0.9,1.2]
eabins = [0,0.04,0.08,0.12,0.16,0.24,0.32,0.4,0.6,0.8,1.0]

myptpz = ROOT.TH2D("ptpz","ptpz",len(pzbins)-1, array.array("d",pzbins),len(ptbins)-1,array.array("d",ptbins))
mypt = ROOT.TH1D("pt","pt",len(ptbins)-1,array.array("d",ptbins))
mypz = ROOT.TH1D("pz","pz",len(pzbins)-1,array.array("d",pzbins))
myptpz_rate = ROOT.TH2D("ptpz_rate","ptpz_rate",len(pzbins)-1, array.array("d",pzbins),len(ptbins)-1,array.array("d",ptbins))
mytheta = ROOT.TH1D("theta","theta",360,0,180)

myptpz_bymode = []
for i in range(0,60):
    tmphist = ROOT.TH2D("ptpz-%d"%i,"ptpz-%d"%i,len(pzbins)-1, array.array("d",pzbins),len(ptbins)-1,array.array("d",ptbins))
    myptpz_bymode.append(tmphist)

myptpz_qelike = ROOT.TH2D("ptpz_qelike","ptpz_qelike",len(pzbins)-1, array.array("d",pzbins),len(ptbins)-1,array.array("d",ptbins))
counter = 0
for e in mytree:
    if not isInclusiveFHC(e): continue
    coslep = e.CosLep
    #20 degree cut
    if coslep<0.93969262078: continue
    elep= e.ELep
    fScaleFactor = e.fScaleFactor
  
    P = ROOT.TMath.Sqrt(elep*elep-0.105*0.105)
    Pl = coslep*P
    Pt = ROOT.TMath.Sqrt(1-coslep*coslep)*P

#    print P,Pl,Pt,ROOT.TMath.Sqrt(Pl*Pl+Pt*Pt)/P,e.Mode

    myptpz.Fill(Pl,Pt,fScaleFactor)
    mypt.Fill(Pt,fScaleFactor)
    mypz.Fill(Pl,fScaleFactor)
    myptpz_rate.Fill(Pl,Pt)
    myptpz_bymode[e.Mode].Fill(Pl,Pt,fScaleFactor)
    if(Pl>1.5 and Pl<60 and Pt>0 and Pt<4.5): mytheta.Fill(ROOT.TMath.ACos(coslep)*180.0/3.14159,fScaleFactor)
    counter+=1
#    if(counter==100):
#        break
    if isCCQELike(e): myptpz_qelike.Fill(Pl,Pt,fScaleFactor)

#myptpz.Scale(1,"width")
#mypt.Scale(1,"width")
#mypz.Scale(1,"width")

myptpz.GetXaxis().SetTitle("Muon p_{||} (GeV)")
myptpz.GetYaxis().SetTitle("Muon p_{t} (GeV)")

mypz.GetXaxis().SetTitle("Muon p_{||} (GeV)")
mypz.GetYaxis().SetTitle("d#sigma/dp_{||} cm^2/GeV/nucleon")

mypt.GetXaxis().SetTitle("Muon p_{t} (GeV)")
mypt.GetYaxis().SetTitle("d#sigma/dp_{t} cm^2/GeV/nucleon")


myoutput = ROOT.TFile(sys.argv[2],"RECREATE")
myptpz.Write();
mypz.Write()
mypt.Write()
myptpz_rate.Write()
mytheta.Write()
for i in range(60):
    myptpz_bymode[i].Write()

#can = ROOT.TCanvas()
#myptpz.Draw("COLZ")

#can2 = ROOT.TCanvas()
#mypt.Draw()

#can3 = ROOT.TCanvas()
#mypz.Draw()


#raw_input("done")

    