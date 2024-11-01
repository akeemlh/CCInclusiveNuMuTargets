#include "GENIEXSecExtract/XSecLooper.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvH2D.h>
#include "PlotUtils/TargetUtils.h"


#include <cstdlib>
typedef unsigned int uint;

class MinModDepCCQEXSec : public XSec
{
public:
  MinModDepCCQEXSec(const char* name)
    :XSec(name)
  {
  };

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
  bool inZRange = true_muon_vtx_z > PlotUtils::TargetProp::Tracker::Face && true_muon_vtx_z < PlotUtils::TargetProp::Tracker::Back;

  bool inERange = (truth_muon_E>2000) && (truth_muon_E<50000);
  std::cout<<"truth_muon_E: " <<truth_muon_E <<"inERange: "<< inERange << std::endl;
  std::cout<<"inApothem: "<< inApothem << std::endl;
  std::cout<<"true_muon_vtx_z: "<< true_muon_vtx_z << std::endl;

  //if(!chw.GetValue("truth_is_fiducial",entry)) return false; //Doesn't work for MasterAnaDev tuples.  What does this even mean in the targets anyway? :(
  if(inAngle  && inERange /* && inZRange && inApothem */ ) return true;
  return false;

}
  // Override this method from the base class to decide what events to
  // include in this selection
  virtual bool passesCuts(PlotUtils::ChainWrapper& chw, int entry)
  {
    if((int)chw.GetValue("mc_incoming", entry)!=14) return false;
    if((int)chw.GetValue("mc_current", entry)!=1) return false;
    if(!isCCInclusiveSignal  ( chw, entry ) ) return false;
    
    return true;
  }
};

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

  // Create the XSecLooper and tell it the input files
  // Inputs should be the merged ntuples:
  XSecLooper loop(playlistFile.c_str());

  // Tell the XSecLooper which neutrino type we're considering (mandatory)
  loop.setNuPDG(14);

  // Setting the number of Universes in the GENIE error band (default 100, put 0 if you do not want to include the universes)
  loop.setNumUniv(0); 
  loop.setFiducial(PlotUtils::TargetProp::Tracker::Face, PlotUtils::TargetProp::Tracker::Back, 850);

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
  for(int whichBin = 0; whichBin < 100 + 1; ++whichBin)
  {
    Erecoil_edges[whichBin-1] = robsRecoilBinWidth * whichBin;
  }
 
  // Flux-integrated over the range 0.0 to 100.0 GeV
  MinModDepCCQEXSec* ds_dpT = new MinModDepCCQEXSec("pT");
  ds_dpT->setBinEdges(pt_nbins, pt_edges);
  ds_dpT->setVariable(XSec::kPTLep);
  ds_dpT->setIsFluxIntegrated(true);
  ds_dpT->setDimension(1);
  ds_dpT->setFluxIntLimits(0, 100.0);
  ds_dpT->setNormalizationType(XSec::kPerNucleon);  
  ds_dpT->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  //loop.addXSec(ds_dpT);

  MinModDepCCQEXSec* ds_dpTdpZ = new MinModDepCCQEXSec("pT");
  ds_dpTdpZ->setBinEdges(pt_nbins, pt_edges, pz_nbins, pz_edges);
  ds_dpTdpZ->setVariable(XSec::kPTLep, XSec::kPZLep);
  ds_dpTdpZ->setIsFluxIntegrated(true);
  ds_dpTdpZ->setDimension(2);
  ds_dpTdpZ->setFluxIntLimits(0.0, 100.0);
  ds_dpTdpZ->setNormalizationType(XSec::kPerNucleon);  
  ds_dpTdpZ->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpTdpZ);

  // Get the output histograms and save them to file
  std::string geniefilename =  "GENIEXSECEXTRACT_" + playlistFile.substr(playlistFile.rfind("/")+1, playlistFile.find(".")) + ".root";
  TFile fout(geniefilename.c_str(), "RECREATE");
  for(uint i=0; i<loop.getXSecs().size(); ++i)
  {
    loop.getXSecs()[i]->getXSecHist()->Write();
    loop.getXSecs()[i]->getEvRateHist()->Write();
  }
  //loop.getFluxHist()->Write();
  std::cout<<"SUCCESS\n";
  return 0;
}
