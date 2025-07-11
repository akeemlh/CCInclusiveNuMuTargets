#include "GENIEXSecExtract/XSecLooper.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvH2D.h>
#include "PlotUtils/TargetUtils.h"
#include "event/CVUniverse.h"
#include "util/NukeUtils.h"

#include <cstdlib>
typedef unsigned int uint;


class CCInclTargets : public XSec
{
public:
  CCInclTargets(const char* name, std::string TargetCode) : 
  XSec(name), 
  fTargetCode(TargetCode){};
  
  std::string fTargetCode;

  double getVariableValue(ChainWrapper& chw, int entry)
  {
      return getValue(chw, entry, m_variable_x);
  }

  double getValue(ChainWrapper& chw, int entry, XSec::EVariable var)
  {
      switch(var){
          case XSec::kEHad: //GeV
            return ( chw.GetValue("mc_incomingE", entry) - chw.GetValue("mc_primFSLepton", entry, 3) );
      }
      return XSec::getValue(chw, entry, var);
  }


  bool isCCInclusiveSignal( PlotUtils::ChainWrapper& chw, int entry )
  {
    double theta              = 0.;
    double true_muon_px   = (double)chw.GetValue("mc_primFSLepton",entry,0)/1000;
    double true_muon_py   = (double)chw.GetValue("mc_primFSLepton",entry,1)/1000;
    double true_muon_pz   = (double)chw.GetValue("mc_primFSLepton",entry,2)/1000;
    double true_muon_vtx_x   = (double)chw.GetValue("mc_vtx",entry,0);
    double true_muon_vtx_y   = (double)chw.GetValue("mc_vtx",entry,1);
    double true_muon_vtx_z   = (double)chw.GetValue("mc_vtx",entry,2);
    double numi_beam_angle_rad = -0.05887;
    double muon_E     = (double)chw.GetValue("mc_primFSLepton",entry,3);
    double pyprime = -1.0*sin(numi_beam_angle_rad)*true_muon_pz + cos(numi_beam_angle_rad)*true_muon_py;
    double pzprime =  1.0*cos(numi_beam_angle_rad)*true_muon_pz + sin(numi_beam_angle_rad)*true_muon_py;
    double pSquare = pow(true_muon_px,2) + pow(pyprime,2) + pow(pzprime,2);
    theta = acos( pzprime / sqrt(pSquare) );
    theta *= 180./3.14159265;
    double truth_muon_E = (double)chw.GetValue("truth_muon_E",entry);

    TVector3 p3lep((double)chw.GetValue("mc_primFSLepton",entry,0),
                  (double)chw.GetValue("mc_primFSLepton",entry,1),
                  (double)chw.GetValue("mc_primFSLepton",entry,2));
    p3lep.RotateX(numi_beam_angle_rad);
    theta = p3lep.Theta();
    //std::cout<<"entry: " << entry <<std::endl;
    //std::cout<<"theta: " << theta << std::endl;
    //theta = p3lep.Theta();
    bool inAngle = theta <= (17.0*M_PI/180); //Used only when using this TVector3 p3lep approach
    //bool inAngle = theta <= 17.0;



    //bool inApothem = (fabs(true_muon_vtx_y) < (-1./sqrt(3.))*fabs(true_muon_vtx_x) + 2.*850/sqrt(3.)) && (fabs(true_muon_vtx_x) < 850);

    //std::cout<<"runxsecloop apothem check: x " << true_muon_vtx_x << " y " << true_muon_vtx_y <<std::endl;

    //bool inZRange = true_muon_vtx_z > PlotUtils::TargetProp::NukeRegion::Face && true_muon_vtx_z < PlotUtils::TargetProp::NukeRegion::Back;

    bool inERange = (truth_muon_E>2000) && (truth_muon_E<50000);
    if(inAngle  && inERange ) return true;
    return false;

  }
  // Override this method from the base class to decide what events to
  // include in this selection
  virtual bool passesCuts(PlotUtils::ChainWrapper& chw, int entry)
  {

    double true_muon_vtx_x   = (double)chw.GetValue("mc_vtx",entry,0);
    double true_muon_vtx_y   = (double)chw.GetValue("mc_vtx",entry,1);
    double true_muon_vtx_z   = (double)chw.GetValue("mc_vtx",entry,2);


    if((int)chw.GetValue("mc_incoming", entry)!=14) return false;
    if((int)chw.GetValue("mc_current", entry)!=1) return false;
    if(!isCCInclusiveSignal  ( chw, entry ) ) 
    { 
      //std::cout<<"Failed signal check\n";
      return false;
    }
    if (!inTarget( chw, entry ))
    {
      return false;
    }
    /* if ((int)chw.GetValue("truth_target_code", entry) != fTargetCode)
    {
      if (fTargetCode == 6008 ||fTargetCode == 6001)
      {
        if (!PlotUtils::TargetUtils::Get().InWaterTargetVolMC(true_muon_vtx_x, true_muon_vtx_y, true_muon_vtx_z)) return false;
      }
      else return false;
    }    */

    //if(!PassTrueDistToDivisionCut( chw, entry, target, nucleus, 25.0)) return false;
    return true;
  }


    bool PassTrueDistToDivisionCut( ChainWrapper& chw, int entry, int targetID, int num_targetZ, double xySep /* = 25. */ )
    {
        double true_target_dist_to_division = (double)chw.GetValue("truth_target_dist_to_division",entry);
	
        int true_targetID = (int)chw.GetValue("truth_targetID",entry);
        int true_targetZ  = (int)chw.GetValue("truth_targetZ",entry);
        int true_module   = (int)chw.GetValue("truth_vtx_module",entry);
        if( targetID < 10 && true_targetID != targetID) return false;
        if( targetID < 10 && true_targetZ != num_targetZ) return false;
	if(targetID >10 && false) {
	  if( true_module < 27 || true_module > 80 ) return false;
	}//NukeCC really does cut off at 80 for all tracker, but 81 for 94 it looks like. But I compared the two and they were the same?? I think there might be a number of targets issue here. (maybe result of the plane code number change faiza did?). No I'm pretty sure that target 94 is just only supposed to go up to 80. Its 1 module bigger than every other one. 
	else{
	  if( targetID == 14 && ( true_module < 27 || true_module > 32)) return false;
	  if( targetID == 24 && ( true_module < 33 || true_module > 38)) return false;
	  if( targetID == 34 && ( true_module < 39 || true_module > 44)) return false;
	  if( targetID == 44 && ( true_module < 45 || true_module > 50)) return false;
	  if( targetID == 54 && ( true_module < 51 || true_module > 56)) return false;
	  if( targetID == 64 && ( true_module < 57 || true_module > 62)) return false;
	  if( targetID == 74 && ( true_module < 63 || true_module > 68)) return false;
	  if( targetID == 84 && ( true_module < 69 || true_module > 74)) return false;
	  if( targetID == 94 && ( true_module < 75 || true_module > 80)) return false;//true_module > 81)) return false;
	}
        //            only relevant for passive targets 1235
	if( targetID<10 && targetID!=4)
        if( 0 < true_targetID && true_targetID < 10 && 4 != true_targetID)
            return ( xySep < true_target_dist_to_division );
        //        cout<<"fiducial event "<<targetID<<endl;
        
        return true;
    }

    bool inTarget(PlotUtils::ChainWrapper& chw, int entry)
    {
      double vtx_x   = (double)chw.GetValue("mc_vtx",entry,0);
      double vtx_y   = (double)chw.GetValue("mc_vtx",entry,1);
      double vtx_z   = (double)chw.GetValue("mc_vtx",entry,2);
      int mod   = (int)chw.GetValue("truth_vtx_module",entry);
      int plane   = (int)chw.GetValue("truth_vtx_plane",entry);
      int tgtCode = util::getTargetCodeFromVtxInfo(vtx_x, vtx_y, vtx_z, mod, plane);
      bool inTgt = false;
      if (fTargetCode=="Iron")
      {
        inTgt = (tgtCode == 2026) || (tgtCode == 3026) || (tgtCode == 5026);
      }
      else if (fTargetCode=="Lead")
      {
        inTgt = (tgtCode == 2082) || (tgtCode == 3082) || (tgtCode == 4082) || (tgtCode == 5082);
      }
      else if (fTargetCode=="Carbon")
      {
        inTgt = (tgtCode == 3006);
      }
      else
      {
        inTgt = std::stoi(fTargetCode) == tgtCode;
      }
      return inTgt;
    }
};

double GetNormValue( std::string tgt)
{
    std::cout<<"tgtCode : "<<tgt<<std::endl;
    double passiveNucleons = 0;
    PlotUtils::TargetUtils tgtUtil;
    if (tgt=="Iron")
    {
      passiveNucleons +=tgtUtil.GetPassiveTargetNNucleons(2, 26, true, 850); 
      passiveNucleons +=tgtUtil.GetPassiveTargetNNucleons(3, 26, true, 850); 
      passiveNucleons +=tgtUtil.GetPassiveTargetNNucleons(5, 26, true, 850); 
    }
    else if (tgt=="Lead")
    {
      passiveNucleons +=tgtUtil.GetPassiveTargetNNucleons(2, 82, true, 850); 
      passiveNucleons +=tgtUtil.GetPassiveTargetNNucleons(3, 82, true, 850); 
      passiveNucleons +=tgtUtil.GetPassiveTargetNNucleons(4, 82, true, 850); 
      passiveNucleons +=tgtUtil.GetPassiveTargetNNucleons(5, 82, true, 850); 
    }
    else if (tgt=="Carbon")
    {
      passiveNucleons +=tgtUtil.GetPassiveTargetNNucleons(3, 6, true, 850); 
    }
    else
    {
      int tgtnum =  std::stoi(tgt);
      double targetZ = tgtnum % 1000;
      double targetID = (tgtnum - targetZ) / 1000;

      if (targetID==6) passiveNucleons = PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(6, 1, true, 850); 
      else if (targetID<7) passiveNucleons = PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(targetID, targetZ, true, 850); 
      else
      {
        if (targetID==7) passiveNucleons = PlotUtils::TargetUtils::Get().GetTrackerNNucleons( 7, true); 
        if (targetID==8) passiveNucleons = PlotUtils::TargetUtils::Get().GetTrackerNNucleons( 6, true); 
        if (targetID==9) passiveNucleons = PlotUtils::TargetUtils::Get().GetTrackerNNucleons( 6, true); 
        if (targetID==10) passiveNucleons = PlotUtils::TargetUtils::Get().GetTrackerNNucleons( 6, true); 
        if (targetID==11) passiveNucleons = PlotUtils::TargetUtils::Get().GetTrackerNNucleons( 6, true); 
        if (targetID==12) passiveNucleons = PlotUtils::TargetUtils::Get().GetTrackerNNucleons( 2, true); 

      }
      if( passiveNucleons < 0 )
          assert( false && "Normalizations only known for Z = 1,6,26,82 and 0(scint)" );

    }
    std::cout<<"passiveNucleons: " <<passiveNucleons<<std::endl;
    
    double trackerAtomsC = TargetUtils::Get().GetTrackerElementNAtoms( 6,  5970, 8450, true, 850.0);
   std::cout<<"trackerAtomsC: " <<trackerAtomsC<<std::endl;
    
    std::cout<<"The normalization factor for this analysis is "<<trackerAtomsC / passiveNucleons<<std::endl;
    
    return trackerAtomsC / passiveNucleons;
    //return 1.0;
}

int main(const int argc, const char** argv)
{
  //Read a playlist file from the command line
  if(argc != 2)
  {
    std::cerr << "Expected exactly 1 command line argument, but got " << argc - 1 << ".\n\n"
              << "USAGE: runXSecLooper <MCPlaylist.txt>\n\n"
              << "MCPlaylist.txt shall contain one .root file per line that has a Truth tree in it.\n"
              << "This program returns 0 when it suceeds.  It produces a .root file with GENIEXSECEXTRACT in its name.\n";
    return 1;
  }

  const std::string playlistFile = argv[1]; //argv[0] is the name of the executable

  //std::vector<int> targetCodes = /* { 1026, 1082, */ {/* 2026, 2082, 3006, 3026, 3082, 4082, 5026, 5082, 6001, 6008,  */7, 8, 9, 10, 11, 12};
  //std::vector<int> targetCodes = { 1026, 1082, 2026, 2082, 3006, 3026, 3082, 4082, 5026, 5082, 6000};
  //std::vector<int> targetCodes = { 2026, 2082, 3006, 3026, 3082, 4082, 5026, 5082, 6000};
  //std::vector<string> targetCodes = { "2026", "2082", "3006", "3026", "3082", "4082", "5026", "5082", "6000 ,7", "8", "9", "10", "11", "12","Iron", "Lead", "Carbon"};
  std::vector<std::string> targetCodes = { "Iron", "Lead", "Carbon"};
  for (auto const& tgtCode : targetCodes)
  {
    std::cout<< "Target: " << tgtCode << std::endl;
    // Create the XSecLooper and tell it the input files
    // Inputs should be the merged ntuples:
    XSecLooper loop(playlistFile.c_str());

    double norm = GetNormValue( tgtCode);
  /* 
    std::string playlistname = "minervame1A";
    size_t mc_label = playlistFile.find("/MC/");
    if (mc_label == string::npos) std::cout << "FAILED TO DETERMINE WHICH PLAYLIST THIS IS. DEFAULTING TO ME1A FOR FLUX\n";
    else if (playlistFile.substr(mc_label+4, 4) == "Test" ) std::cout << "IDENTIFIED TEST PLAYLIST. DEFAULTING TO ME1A FOR FLUX\n";
    else
    {
      std::string playlistkey = playlistFile.substr(mc_label+4, 2);
      std::cout<<"Identified Playlist " << playlistkey << std::endl;
      playlistname = "minervame"+playlistkey;
    }

    loop.setPlaylist(FluxReweighter::GetPlaylistEnum(playlistname));
    */

    // Tell the XSecLooper which neutrino type we're considering (mandatory)
    loop.setNuPDG(14);

    // Setting the number of Universes in the GENIE error band (default 100, put 0 if you do not want to include the universes)
    loop.setNumUniv(0); 

    loop.setFiducial(5970, 8450);
    //loop.setFiducial(PlotUtils::TargetProp::Tracker::Face, PlotUtils::TargetProp::Tracker::Back, 850.0);
    //loop.setFiducial(PlotUtils::TargetProp::NukeRegion::Face, PlotUtils::TargetProp::NukeRegion::Back);
    loop.setPlaylist(PlotUtils::FluxReweighter::minervame1A);

    // Add the differential cross section dsigma/ds_dpT
    double pt_edges[] = { 0.0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1.0, 1.25, 1.5, 2.5, 4.5 };
    int pt_nbins = 14;

    double pz_edges[] = {1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20, 40, 60};
    int pz_nbins = 16; 

    double eMu_edges[] = {0,1,2,3,4,5,7,9,12,15,18,22,36,50,75,100,120};
    int eMu_nbins = 16; 
    double bjorken_edges[] = {0.0, 0.1, 0.3, 0.5, 0.7, 0.9 , 1.1, 1.5};
    int bjorken_nbins = 7; 


    double Erecoil_edges[101];
    int Erecoil_nbins = 100;

    const double robsRecoilBinWidth = 50;
    for(int whichBin = 0; whichBin < Erecoil_nbins+1 ; whichBin++)
    {
      Erecoil_edges[whichBin] = robsRecoilBinWidth * whichBin;
    }
  
    // Flux-integrated over the range 0.0 to 100.0 GeV
    CCInclTargets* ds_dpT = new CCInclTargets("pTmu", tgtCode);
    ds_dpT->setBinEdges(pt_nbins, pt_edges);
    ds_dpT->setVariable(XSec::kPTLep);
    ds_dpT->setIsFluxIntegrated(true);
    ds_dpT->setDimension(1);
    ds_dpT->setFluxIntLimits(0, 100.0);
    ds_dpT->setNormalizationValue( norm );
    ds_dpT->setNormalizationType(XSec::kSelfNorm);
    ds_dpT->setUniverses(0); //default value, put 0 if you do not want universes to be included.
    loop.addXSec(ds_dpT);


    CCInclTargets* ds_dpZ = new CCInclTargets("pZmu", tgtCode);
    ds_dpZ->setBinEdges(pz_nbins, pz_edges);
    ds_dpZ->setVariable(XSec::kPZLep);
    ds_dpZ->setIsFluxIntegrated(true);
    ds_dpZ->setDimension(1);
    ds_dpZ->setFluxIntLimits(0.0, 100.0);
    ds_dpZ->setNormalizationValue( norm );
    ds_dpZ->setNormalizationType(XSec::kSelfNorm);
    ds_dpZ->setUniverses(0); //default value, put 0 if you do not want universes to be included.
    loop.addXSec(ds_dpZ);

    CCInclTargets* ds_dpEMu = new CCInclTargets("Emu", tgtCode);
    ds_dpEMu->setBinEdges(eMu_nbins, eMu_edges);
    ds_dpEMu->setVariable(XSec::kELep);
    ds_dpEMu->setIsFluxIntegrated(true);
    ds_dpEMu->setDimension(1);
    ds_dpEMu->setFluxIntLimits(0.0, 100.0);
    ds_dpEMu->setNormalizationValue( norm );
    ds_dpEMu->setNormalizationType(XSec::kSelfNorm);
    ds_dpEMu->setUniverses(0); //default value, put 0 if you do not want universes to be included.
    loop.addXSec(ds_dpEMu);

    CCInclTargets* ds_dpERecoil = new CCInclTargets("Erecoil", tgtCode);
    ds_dpERecoil->setBinEdges(Erecoil_nbins, Erecoil_edges);
    ds_dpERecoil->setVariable(XSec::kEHad);
    ds_dpERecoil->setIsFluxIntegrated(true);
    ds_dpERecoil->setDimension(1);
    ds_dpERecoil->setFluxIntLimits(0.0, 100.0);
    ds_dpERecoil->setNormalizationValue( norm );
    ds_dpERecoil->setNormalizationType(XSec::kSelfNorm);
    ds_dpERecoil->setUniverses(0); //default value, put 0 if you do not want universes to be included.
    loop.addXSec(ds_dpERecoil);

    CCInclTargets* ds_dpMeasX = new CCInclTargets("measX", tgtCode);
    ds_dpMeasX->setBinEdges(bjorken_nbins, bjorken_edges);
    ds_dpMeasX->setVariable(XSec::kxExp);
    ds_dpMeasX->setIsFluxIntegrated(true);
    ds_dpMeasX->setDimension(1);
    ds_dpMeasX->setFluxIntLimits(0.0, 100.0);
    ds_dpMeasX->setNormalizationValue( norm );
    ds_dpMeasX->setNormalizationType(XSec::kSelfNorm);
    ds_dpMeasX->setUniverses(0); //default value, put 0 if you do not want universes to be included.
    loop.addXSec(ds_dpMeasX);

    CCInclTargets* ds_dpX = new CCInclTargets("BjorkenX", tgtCode);
    ds_dpX->setBinEdges(bjorken_nbins, bjorken_edges);
    ds_dpX->setVariable(XSec::kx);
    ds_dpX->setIsFluxIntegrated(true);
    ds_dpX->setDimension(1);
    ds_dpX->setFluxIntLimits(0.0, 100.0);
    ds_dpX->setNormalizationValue( norm );
    ds_dpX->setNormalizationType(XSec::kSelfNorm);
    ds_dpX->setUniverses(0); //default value, put 0 if you do not want universes to be included.
    loop.addXSec(ds_dpX);

    CCInclTargets* ds_dpTdpZ = new CCInclTargets("pTpZ", tgtCode);
    ds_dpTdpZ->setBinEdges(pz_nbins, pz_edges, pt_nbins, pt_edges);
    ds_dpTdpZ->setVariable(XSec::kPZLep, XSec::kPTLep);
    ds_dpTdpZ->setIsFluxIntegrated(true);
    ds_dpTdpZ->setDimension(2);
    ds_dpTdpZ->setFluxIntLimits(0.0, 100.0);
    ds_dpTdpZ->setNormalizationValue( norm );
    ds_dpTdpZ->setNormalizationType(XSec::kSelfNorm);
    ds_dpTdpZ->setUniverses(0); //default value, put 0 if you do not want universes to be included.
    loop.addXSec(ds_dpTdpZ);

    loop.runLoop();

    // Get the output histograms and save them to file
    std::string geniefilename =  "GENIEXSECEXTRACT_Tgt"+tgtCode+".root";
    TFile fout(geniefilename.c_str(), "RECREATE");
    for(uint i=0; i<loop.getXSecs().size(); ++i)
    {
      if (loop.getXSecs()[i]->getDimension()==1)
      {
        loop.getXSecs()[i]->getXSecHist()->Write();
        loop.getXSecs()[i]->getEvRateHist()->Write();
        loop.getXSecs()[i]->integratedFlux->Write();
        loop.getXSecs()[i]->flux->Write();
        loop.getFluxHist()->Write();
      }
      else if (loop.getXSecs()[i]->getDimension()==2)
      {
        loop.getXSecs()[i]->get2DXSecHist()->Write();
      }
    }
  }
  std::cout<<"SUCCESS\n";
  return 0;
}
