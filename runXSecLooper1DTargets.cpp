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
  CCInclTargets(const char* name, int targetZ) : 
  XSec(name), 
  fTargetZ(targetZ){};
  
  int fTargetZ;

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
    //theta = p3lep.Theta();
    //bool inAngle = theta <= (17.0*M_PI/180); //Used only when using this TVector3 p3lep approach

    bool inAngle = theta <= 17.0;



    bool inApothem = (fabs(true_muon_vtx_y) < (-1./sqrt(3.))*fabs(true_muon_vtx_x) + 2.*850/sqrt(3.))
          && (fabs(true_muon_vtx_x) < 850);
    bool inZRange = true_muon_vtx_z > PlotUtils::TargetProp::NukeRegion::Face && true_muon_vtx_z < PlotUtils::TargetProp::NukeRegion::Back;

    bool inERange = (truth_muon_E>2000) && (truth_muon_E<50000);
    if(inAngle  && inERange /* && inZRange && inApothem */ ) return true;
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
    if(!isCCInclusiveSignal  ( chw, entry ) ) return false;

    if ((int)chw.GetValue("truth_targetID", entry) == 1) return false; //Discarding all target 1 interactions
    if (fTargetZ == 6)
    {
      //if ((int)chw.GetValue("truth_targetZ", entry) != 6) return false;
      if ((int)chw.GetValue("truth_target_code", entry) != 6) return false;
    }
    if (fTargetZ == 26)
    {
      //if ((int)chw.GetValue("truth_targetZ", entry) != 26) return false;
      if ((int)chw.GetValue("truth_target_code", entry) != 26) return false;
    }
    if (fTargetZ == 82)
    {
      //if ((int)chw.GetValue("truth_targetZ", entry) != 82) return false;
      if ((int)chw.GetValue("truth_target_code", entry) != 82) return false;
    }
    if (fTargetZ == 18)
    {
      if (!PlotUtils::TargetUtils::Get().InWaterTargetVolMC(true_muon_vtx_x, true_muon_vtx_y, true_muon_vtx_z)) return false;
    }
    return true;
  }
};

double GetNormValue( int targetZ )
{
    
    double trackerAtomsC = TargetUtils::Get().GetTrackerElementNAtoms( 6, PlotUtils::TargetProp::NukeRegion::Face, PlotUtils::TargetProp::NukeRegion::Back, true, 850.0);
    
   std::cout<<"trackerAtomsC: " <<trackerAtomsC<<std::endl;
    
    double passiveNucleons = 0;
    if (targetZ == 6)
    {
      passiveNucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(3, 6, true);
    }
    if (targetZ == 26)
    {
      //passiveNucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(1, 26, true);
      passiveNucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(2, 26, true);
      passiveNucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(3, 26, true);
      passiveNucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(5, 26, true);
    }
    if (targetZ == 82)
    {
      //passiveNucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(1, 82, true);
      passiveNucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(2, 82, true);
      passiveNucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(3, 82, true);
      passiveNucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(4, 82, true);
      passiveNucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(5, 82, true);
    }   
    if (targetZ == 18)
    {
      //passiveNucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(1, 82, true);
      passiveNucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(6, 1, true);
    }   
    if( passiveNucleons < 0 )
        assert( false && "Normalizations only known for Z = 1,6,26,82 and 0(scint)" );
    
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


  std::vector<int> targetZs = { 6, 18, 26, 82};
  for (auto const& tgtZ : targetZs)
  {
    std::cout<< "Target: " << tgtZ << std::endl;
    // Create the XSecLooper and tell it the input files
    // Inputs should be the merged ntuples:
    XSecLooper loop(playlistFile.c_str());

    double norm = GetNormValue( tgtZ );
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

    loop.setPlaylist(FluxReweighter::minervame1A);


    // Tell the XSecLooper which neutrino type we're considering (mandatory)
    loop.setNuPDG(14);

    // Setting the number of Universes in the GENIE error band (default 100, put 0 if you do not want to include the universes)
    loop.setNumUniv(0); 
    loop.setFiducial(PlotUtils::TargetProp::NukeRegion::Face, PlotUtils::TargetProp::NukeRegion::Back, 850);

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

    const double robsRecoilBinWidth = 50; //MeV
    for(int whichBin = 0; whichBin < Erecoil_nbins+1 ; ++whichBin)
    {
      Erecoil_edges[whichBin-1] = robsRecoilBinWidth * whichBin;
      std::cout << "whichBin: " << whichBin << std::endl;
    }
  
    // Flux-integrated over the range 0.0 to 100.0 GeV
    CCInclTargets* ds_dpT = new CCInclTargets("pTmu", tgtZ);
    ds_dpT->setBinEdges(pt_nbins, pt_edges);
    ds_dpT->setVariable(XSec::kPTLep);
    ds_dpT->setIsFluxIntegrated(true);
    ds_dpT->setDimension(1);
    ds_dpT->setFluxIntLimits(0, 100.0);
    ds_dpT->setNormalizationValue( norm );
    ds_dpT->setNormalizationType(XSec::kSelfNorm);
    ds_dpT->setUniverses(0); //default value, put 0 if you do not want universes to be included.

    loop.addXSec(ds_dpT);


    CCInclTargets* ds_dpZ = new CCInclTargets("pZmu", tgtZ);
    ds_dpZ->setBinEdges(pz_nbins, pz_edges);
    ds_dpZ->setVariable(XSec::kPZLep);
    ds_dpZ->setIsFluxIntegrated(true);
    ds_dpZ->setDimension(1);
    ds_dpZ->setFluxIntLimits(0.0, 100.0);
    ds_dpZ->setNormalizationValue( norm );
    ds_dpZ->setNormalizationType(XSec::kSelfNorm);
    ds_dpZ->setUniverses(0); //default value, put 0 if you do not want universes to be included.
    loop.addXSec(ds_dpZ);

    CCInclTargets* ds_dpEMu = new CCInclTargets("Emu", tgtZ);
    ds_dpEMu->setBinEdges(eMu_nbins, eMu_edges);
    ds_dpEMu->setVariable(XSec::kELep);
    ds_dpEMu->setIsFluxIntegrated(true);
    ds_dpEMu->setDimension(1);
    ds_dpEMu->setFluxIntLimits(0.0, 100.0);
    ds_dpEMu->setNormalizationValue( norm );
    ds_dpEMu->setNormalizationType(XSec::kSelfNorm);
    ds_dpEMu->setUniverses(0); //default value, put 0 if you do not want universes to be included.
    loop.addXSec(ds_dpEMu);

    CCInclTargets* ds_dpERecoil = new CCInclTargets("Erecoil", tgtZ);
    ds_dpERecoil->setBinEdges(Erecoil_nbins, Erecoil_edges);
    ds_dpERecoil->setVariable(XSec::kEHad);
    ds_dpERecoil->setIsFluxIntegrated(true);
    ds_dpERecoil->setDimension(1);
    ds_dpERecoil->setFluxIntLimits(0.0, 100.0);
    ds_dpERecoil->setNormalizationValue( norm );
    ds_dpERecoil->setNormalizationType(XSec::kSelfNorm);
    ds_dpERecoil->setUniverses(0); //default value, put 0 if you do not want universes to be included.
    loop.addXSec(ds_dpERecoil);

    CCInclTargets* ds_dpMeasX = new CCInclTargets("measX", tgtZ);
    ds_dpMeasX->setBinEdges(bjorken_nbins, bjorken_edges);
    ds_dpMeasX->setVariable(XSec::kxExp);
    ds_dpMeasX->setIsFluxIntegrated(true);
    ds_dpMeasX->setDimension(1);
    ds_dpMeasX->setFluxIntLimits(0.0, 100.0);
    ds_dpMeasX->setNormalizationValue( norm );
    ds_dpMeasX->setNormalizationType(XSec::kSelfNorm);
    ds_dpMeasX->setUniverses(0); //default value, put 0 if you do not want universes to be included.
    loop.addXSec(ds_dpMeasX);

    CCInclTargets* ds_dpX = new CCInclTargets("BjorkenX", tgtZ);
    ds_dpX->setBinEdges(bjorken_nbins, bjorken_edges);
    ds_dpX->setVariable(XSec::kx);
    ds_dpX->setIsFluxIntegrated(true);
    ds_dpX->setDimension(1);
    ds_dpX->setFluxIntLimits(0.0, 100.0);
    ds_dpX->setNormalizationValue( norm );
    ds_dpX->setNormalizationType(XSec::kSelfNorm);
    ds_dpX->setUniverses(0); //default value, put 0 if you do not want universes to be included.
    loop.addXSec(ds_dpX);

    CCInclTargets* ds_dpTdpZ = new CCInclTargets("pTpZ", tgtZ);
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
    std::string geniefilename =  "GENIEXSECEXTRACT_" + playlistFile.substr(playlistFile.rfind("/")+1, playlistFile.find(".")) + "Tgt"+std::to_string(tgtZ)+".root";
    TFile fout(geniefilename.c_str(), "RECREATE");
    for(uint i=0; i<loop.getXSecs().size(); ++i)
    {
      if (loop.getXSecs()[i]->getDimension()==1)
      {
        loop.getXSecs()[i]->getXSecHist()->Write();
        loop.getXSecs()[i]->getEvRateHist()->Write();
        loop.getXSecs()[i]->integratedFlux->Write();
        loop.getXSecs()[i]->flux->Write();
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
