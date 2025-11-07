#define HELP \
"\n*** Help: ***\n"\
" File: Extract1DCrossSectionTrackerDaisy.cpp\n"\
" Brief: Given data and MC files produced by runEventLoopTracker.cpp, extract a 1D differential cross section.\n"\
"        Automatically merges input root files. Directories must have structure of /playlist/runEventLoopTracker<MC-or-Data>.root\n"\
"        Example directory structure:\n\n"\
"        /path/to/dirs/\n"\
"        ├──1A\n"\
"        │   ├──runEventLoopTrackerData.root\n"\
"        │   ├──runEventLoopTrackerMC.root\n"\
"        │   ├──runEventLoopTrackerData_petal_0.root\n"\
"        │   ├──runEventLoopTrackerMC_petal_0.root\n"\
"        │   ├──...\n"\
"        │   ├──runEventLoopTrackerData_petal_11.root\n"\
"        │   ├──runEventLoopTrackerMC_petal_11.root\n"\
"        ├──1B\n"\
"        │   ├──runEventLoopTrackerData.root\n"\
"        │   ├──runEventLoopTrackerMC.root\n"\
"        │   ├──runEventLoopTrackerData_petal_0.root\n"\
"        │   ├──runEventLoopTrackerMC_petal_0.root\n"\
"        │   ├──...\n"\
"        │   ├──runEventLoopTrackerData_petal_11.root\n"\
"        │   ├──runEventLoopTrackerMC_petal_11.root\n"\
"        └──.....\n\n"\
"        Automatically identifies and searches all playlist directories within /path/to/dirs/ and\n"\
"        merges files for appropriate neutrino beam mode\n"\
"        Subtracts backgrounds, performs unfolding, applies efficiency x acceptance correction,\n"\
"        divides by flux and number of nucleons and if option is selected and performs daisy reweight.\n"\
"        Writes a .root file with the cross section histograms\n"\
"        To run for a single playlist (for example 1A) simply pass /path/to/dirs/1A as the directory path\n\n"\
" Usage: ExtractCrossSectionTrackerDaisy <unfolding iterations> <directory> <pdg>\n"\
"        e.g:   ExtractCrossSectionTrackerDaisy 5 /path/to/dirs 14 -- to extract xsecs for neutrino-mode playlists with 5 iterations\n"\
"        e.g:   ExtractCrossSectionTrackerDaisy 10 /path/to/dirs -14 -- to extract xsecs for antineutrino-mode playlists with 10 iterations\n\n"\

// Author: Akeem Hart a.l.hart@qmul.ac.uk based on original ExtractCrossSection.cpp by Andrew Olivier aolivier@ur.rochester.edu

// util includes
#include "util/GetIngredient.h"

// UnfoldUtils includes
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#include "MinervaUnfold/MnvUnfold.h"

// PlotUtils includes
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvPlotter.h"
#include "util/NukeUtils.h"
#include "PlotUtils/TargetUtils.h"
#include "PlotUtils/FluxReweighter.h"
#pragma GCC diagnostic pop

// ROOT includes
#include "TH1D.h"
#include "TFile.h"
#include "TKey.h"
#include "TParameter.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLatex.h"
#include "string.h"

// Cintex is only needed for older ROOT versions like the GPVMs.
////Let CMake decide whether it's needed.
#ifndef NCINTEX
#include "Cintex/Cintex.h"
#endif

// c++ includes
#include <iostream>
#include <exception>
#include <algorithm>
#include <numeric>

// Convince the STL to talk to TIter so I can use std::find_if()
namespace std
{
  template <>
  struct iterator_traits<TIter>
  {
    using value_type = TObject;
    using pointer = TObject *;
    using reference = TObject &;
    using iterator_category = forward_iterator_tag;
  };
}

double GetTotalScatteringCenters(int targetZ, bool isMC)
{
  // TARGET INFO
  PlotUtils::TargetUtils targetInfo;
  double Nucleons = 0;

  // Target 1 is generally excluded due to rock muon contamination (in the inclusive analysis)
  if (targetZ == 6)
  {
    Nucleons = targetInfo.GetPassiveTargetNNucleons(3, targetZ, isMC); // Target 3
  }

  if (targetZ == 26)
  {
    Nucleons = targetInfo.GetPassiveTargetNNucleons(2, targetZ, isMC)    // Target 2
               + targetInfo.GetPassiveTargetNNucleons(3, targetZ, isMC)  // Target 3
               + targetInfo.GetPassiveTargetNNucleons(5, targetZ, isMC); // Target 5
  }

  if (targetZ == 82)
  {
    Nucleons = targetInfo.GetPassiveTargetNNucleons(2, targetZ, isMC)    // Target 2
               + targetInfo.GetPassiveTargetNNucleons(3, targetZ, isMC)  // Target 3
               + targetInfo.GetPassiveTargetNNucleons(4, targetZ, isMC)  // Target 4
               + targetInfo.GetPassiveTargetNNucleons(5, targetZ, isMC); // Target 5
  }
  if (targetZ > 90)
  {
    Nucleons = targetInfo.GetTrackerNNucleons(5980, 8422, isMC, 850);
    // double TargetUtils::GetTrackerNNucleons( double minZ, double maxZ, bool isMC, double apothem /* = 850. */ ) const
  }

  return Nucleons;
}

// Plot a step in cross section extraction.
//Doesn't yet work for 2D
void Plot(PlotUtils::MnvH2D &hist, const std::string &stepName, const std::string &prefix, const std::string &target)
{
  /* bool plotpngs = false;
  TCanvas can(stepName.c_str());
  hist.GetCVHistoWithError().Clone()->Draw();
  if (plotpngs) can.Print((target + "_" + prefix + "_" + stepName + ".png").c_str());

  // Uncertainty summary
  PlotUtils::MnvPlotter plotter;
  plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
  plotter.axis_maximum = 0.4;

  plotter.DrawErrorSummary(&hist);
  if (plotpngs) can.Print((target + "_" + prefix + "_" + stepName + "_uncertaintySummary.png").c_str());

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "Other");
  if (plotpngs) can.Print((target + "_" + prefix + "_" + stepName + "_otherUncertainties.png").c_str()); */
}

// Plot a step in cross section extraction.
void Plot(PlotUtils::MnvH1D &hist, const std::string &stepName, const std::string &prefix)
{
  bool plotpngs = false;
  TCanvas can(stepName.c_str());
  hist.GetCVHistoWithError().Clone()->Draw();
  if (plotpngs) can.Print((prefix + "_" + stepName + ".png").c_str());

  // Uncertainty summary
  PlotUtils::MnvPlotter plotter;
  plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
  plotter.axis_maximum = 0.4;

  plotter.DrawErrorSummary(&hist);
  if (plotpngs) can.Print((prefix + "_" + stepName + "_uncertaintySummary.png").c_str());

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "Other");
  if (plotpngs) can.Print((prefix + "_" + stepName + "_otherUncertainties.png").c_str());
}

// Unfolding function from Aaron Bercelle
// TODO: Trim it down a little?  Remove that static?
PlotUtils::MnvH1D *UnfoldHist(PlotUtils::MnvH1D *h_folded, PlotUtils::MnvH2D *h_migration, int num_iter)
{
  static MinervaUnfold::MnvUnfold unfold;
  PlotUtils::MnvH1D *h_unfolded = nullptr;

  // bool bUnfolded = false;

  TMatrixD dummyCovMatrix;
  if (!unfold.UnfoldHisto(h_unfolded, dummyCovMatrix, h_migration, h_folded, RooUnfold::kBayes, num_iter, true, false))
    return nullptr;

  /////////////////////////////////////////////////////////////////////////////////////////
  // No idea if this is still needed
  // Probably.  This gets your stat unfolding covariance matrix
  TMatrixD unfoldingCovMatrixOrig;
  int correctNbins;
  int matrixRows;
  TH1D *hUnfoldedDummy = new TH1D(h_unfolded->GetCVHistoWithStatError());
  TH1D *hRecoDummy = new TH1D(h_migration->ProjectionX()->GetCVHistoWithStatError());
  TH1D *hTruthDummy = new TH1D(h_migration->ProjectionY()->GetCVHistoWithStatError());
  TH1D *hBGSubDataDummy = new TH1D(h_folded->GetCVHistoWithStatError());
  TH2D *hMigrationDummy = new TH2D(h_migration->GetCVHistoWithStatError());
  unfold.UnfoldHisto(hUnfoldedDummy, unfoldingCovMatrixOrig, hMigrationDummy, hRecoDummy, hTruthDummy, hBGSubDataDummy, RooUnfold::kBayes, num_iter); // Stupid RooUnfold.  This is dummy, we don't need iterations

  correctNbins = hUnfoldedDummy->fN;
  matrixRows = unfoldingCovMatrixOrig.GetNrows();
  if (correctNbins != matrixRows)
  {
    std::cout << "****************************************************************************" << std::endl;
    std::cout << "*  Fixing unfolding matrix size because of RooUnfold bug. From " << matrixRows << " to " << correctNbins << std::endl;
    std::cout << "****************************************************************************" << std::endl;
    // It looks like this, since the extra last two bins don't have any content
    unfoldingCovMatrixOrig.ResizeTo(correctNbins, correctNbins);
  }

  for (int i = 0; i < unfoldingCovMatrixOrig.GetNrows(); ++i)
    unfoldingCovMatrixOrig(i, i) = 0;
  delete hUnfoldedDummy;
  delete hMigrationDummy;
  delete hRecoDummy;
  delete hTruthDummy;
  delete hBGSubDataDummy;
  h_unfolded->PushCovMatrix("unfoldingCov", unfoldingCovMatrixOrig);

  /////////////////////////////////////////////////////////////////////////////////////////
  return h_unfolded;
}

// The final step of cross section extraction: normalize by flux, bin width, POT, and number of targets
PlotUtils::MnvH1D *normalize(PlotUtils::MnvH1D *efficiencyCorrected, PlotUtils::MnvH1D *fluxIntegral, const double nNucleons, const double POT)
{
  efficiencyCorrected->Divide(efficiencyCorrected, fluxIntegral);

  efficiencyCorrected->Scale(1. / nNucleons / POT);
  efficiencyCorrected->Scale(1.e4); // Flux histogram is in m^-2, but convention is to report cm^2
  efficiencyCorrected->Scale(1., "width");

  return efficiencyCorrected;
}


//Unfolding function from Aaron Bercelle
//TODO: Trim it down a little?  Remove that static?
PlotUtils::MnvH2D* UnfoldHist( PlotUtils::MnvH2D* h_folded, PlotUtils::MnvH2D* h_migration, PlotUtils::MnvH2D* h_reco, PlotUtils::MnvH2D* h_truth,  int num_iter )
{
  static MinervaUnfold::MnvUnfold unfold;
  PlotUtils::MnvH2D* h_unfolded = nullptr;

  //bool bUnfolded = false;



  if(!unfold.UnfoldHisto2D( h_unfolded, h_migration, h_reco, h_truth, h_folded, num_iter, true, false ))
    return nullptr;

  /////////////////////////////////////////////////////////////////////////////////////////  
  //No idea if this is still needed
  //Probably.  This gets your stat unfolding covariance matrix
  TMatrixD unfoldingCovMatrixOrig; 
  int correctNbins;
  int matrixRows;  
  TH2D* hUnfoldedDummy  = new TH2D(h_unfolded->GetCVHistoWithStatError());
  TH2D* hRecoDummy      = new TH2D(h_reco->GetCVHistoWithStatError());
  TH2D* hTruthDummy     = new TH2D(h_truth->GetCVHistoWithStatError());
  TH2D* hBGSubDataDummy = new TH2D(h_folded->GetCVHistoWithStatError());
  TH2D* hMigrationDummy = new TH2D(h_migration->GetCVHistoWithStatError());
  unfold.UnfoldHisto2D(hUnfoldedDummy ,unfoldingCovMatrixOrig, hMigrationDummy, hRecoDummy, hTruthDummy, hBGSubDataDummy, num_iter);

  correctNbins=hUnfoldedDummy->fN;
  matrixRows=unfoldingCovMatrixOrig.GetNrows();
  if(correctNbins!=matrixRows){
    std::cout << "****************************************************************************" << std::endl;
    std::cout << "*  Fixing unfolding matrix size because of RooUnfold bug. From " << matrixRows << " to " << correctNbins << std::endl;
    std::cout << "****************************************************************************" << std::endl;
    // It looks like this, since the extra last two bins don't have any content
    unfoldingCovMatrixOrig.ResizeTo(correctNbins, correctNbins);
  }

  for(int i=0; i<unfoldingCovMatrixOrig.GetNrows(); ++i) unfoldingCovMatrixOrig(i,i)=0;
  delete hUnfoldedDummy;
  delete hMigrationDummy;
  delete hRecoDummy;
  delete hTruthDummy;
  delete hBGSubDataDummy;
  h_unfolded->PushCovMatrix("unfoldingCov",unfoldingCovMatrixOrig);

  /////////////////////////////////////////////////////////////////////////////////////////  
  return h_unfolded;
}

// The final step of cross section extraction: normalize by flux, bin width, POT, and number of targets
PlotUtils::MnvH2D *normalize(PlotUtils::MnvH2D *efficiencyCorrected, PlotUtils::MnvH2D *fluxIntegral, const double nNucleons, const double POT)
{
  efficiencyCorrected->Divide(efficiencyCorrected, fluxIntegral);

  efficiencyCorrected->Scale(1. / nNucleons / POT);
  efficiencyCorrected->Scale(1.e4); // Flux histogram is in m^-2, but convention is to report cm^2
  efficiencyCorrected->Scale(1., "width");

  return efficiencyCorrected;
}

int main(const int argc, const char **argv)
{
#ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable(); // Needed to look up dictionaries for PlotUtils classes like MnvH1D
#endif

  TH1::AddDirectory(kFALSE); // Needed so that MnvH1D gets to clean up its own MnvLatErrorBands (which are TH1Ds).

  if (argc != 4)
  {
    std::cerr << "Expected 4 arguments, but I got " << argc - 1 << ".\n"  << HELP <<  std::endl;  
    return 1;
  }

  const int nIterations = std::stoi(argv[1]);
  std::string indir = std::string(argv[2]);
  int pdg = std::stoi(argv[3]);

  bool doDaisy = true;

  std::vector<std::string> dirs = util::findContainingDirectories(indir, "Tracker", true);

  if (dirs.size()==0) dirs = {"./"};
  
  std::vector<std::string> crossSectionPrefixes = {"pTmu_pZmu"};
  //std::vector<std::string> crossSectionPrefixes = {"pTmu_pZmu", "Emu_Erecoil"};

  double mcPOT = 0;
  double dataPOT = 0;

  std::vector<std::pair<std::string, double>> playlistPOTpair; //Used for daisy reweight POT normalising fluxes

  for (const auto &prefix : crossSectionPrefixes)
  {
    //if (!(prefix == "Erecoil" || prefix == "pTmu")) continue; //Used for testing with only subset of prefixes
    std::cout << "Currently working on variable: " << prefix << std::endl;
    try
    {
      //Flux parameters
      int n_flux_universes = 100; // Is this right
      const bool use_nue_constraint = true;
      const std::string project_dir = "targets_2345_jointNueIMD";
      double min_energy = 0;
      double max_energy = 100;
      
      PlotUtils::MnvH2D *fluxIntReweighted = new PlotUtils::MnvH2D();
      std::shared_ptr<PlotUtils::MnvH2D> DaisyEffNum[12], DaisyFolded[12];
      std::shared_ptr<PlotUtils::MnvH2D> DaisyMigration[12], DaisyMigrationReco[12], DaisyMigrationTruth[12];
      std::vector<PlotUtils::MnvH2D*> DaisyBackgrounds[12];
      std::shared_ptr<PlotUtils::MnvH2D> DaisyEffDenom[12], DaisyEffDenom2P2H[12], DaisyEffDenomDIS[12], DaisyEffDenomRES[12], DaisyEffDenomQE[12], DaisyEffDenomOther[12];

      for (int c = 0; c<12; c++)
      {
        DaisyEffNum[c] = std::shared_ptr<PlotUtils::MnvH2D>(new PlotUtils::MnvH2D());
        DaisyFolded[c] = std::shared_ptr<PlotUtils::MnvH2D>(new PlotUtils::MnvH2D());
        DaisyMigration[c] = std::shared_ptr<PlotUtils::MnvH2D>(new PlotUtils::MnvH2D());
        DaisyMigrationReco[c] = std::shared_ptr<PlotUtils::MnvH2D>(new PlotUtils::MnvH2D());
        DaisyMigrationTruth[c] = std::shared_ptr<PlotUtils::MnvH2D>(new PlotUtils::MnvH2D());
        DaisyEffDenom[c] = std::shared_ptr<PlotUtils::MnvH2D>(new PlotUtils::MnvH2D());
        DaisyEffDenom2P2H[c] = std::shared_ptr<PlotUtils::MnvH2D>(new PlotUtils::MnvH2D());
        DaisyEffDenomDIS[c] = std::shared_ptr<PlotUtils::MnvH2D>(new PlotUtils::MnvH2D());
        DaisyEffDenomRES[c] = std::shared_ptr<PlotUtils::MnvH2D>(new PlotUtils::MnvH2D());
        DaisyEffDenomQE[c] = std::shared_ptr<PlotUtils::MnvH2D>(new PlotUtils::MnvH2D());
        DaisyEffDenomOther[c] = std::shared_ptr<PlotUtils::MnvH2D>(new PlotUtils::MnvH2D());
      }

      PlotUtils::MnvH2D* flux = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* folded = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* migration = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* migration_reco = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* migration_truth = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* effNum = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* effDenom = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* effDenom2P2H = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* effDenomDIS = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* effDenomRES = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* effDenomQE = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* effDenomOther = new PlotUtils::MnvH2D();

      PlotUtils::MnvH2D* BackgroundWrongSign = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* BackgroundNC = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* BackgroundOther = new PlotUtils::MnvH2D();

      PlotUtils::MnvH2D* SelectedSignalReco = new PlotUtils::MnvH2D();

      PlotUtils::MnvH2D* MCData = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* MCData2p2h = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* MCDataDIS = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* MCDataRES = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* MCDataQE = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* MCDataOther = new PlotUtils::MnvH2D();

      std::vector<PlotUtils::MnvH2D *> backgrounds;

      for (int c = 0; c<dirs.size(); c++)
      {
        std::cout<<"Investigating directory " << dirs[c] << " which is " << c+1 <<"/"<<dirs.size() <<" playlists identified" <<std::endl;

        //***********************************************
        //Daisy reweight
        //Getting ingredients
        //***********************************************
        doDaisy = true;
        for (int petal=0; petal<12; petal++){
          std::string datapathdaisy = dirs[c] + "/runEventLoopTrackerData_petal_"+std::to_string(petal)+".root";
          std::string mcpathdaisy = dirs[c] + "/runEventLoopTrackerMC_petal_"+std::to_string(petal)+".root";
          std::string migrationpathdaisy = dirs[c] + "/runEventLoopTracker2DMigration_petal_"+std::to_string(petal)+".root";
          bool petalFilesExist = std::filesystem::exists(datapathdaisy) && std::filesystem::exists(mcpathdaisy) && std::filesystem::exists(migrationpathdaisy);
          std::cout<< "petalFilesExist: "<< petalFilesExist<<std::endl;
          doDaisy = doDaisy && petalFilesExist;
        }
        //doDaisy=false;#
        if (doDaisy) std::cout<<"Will do daisy petal reweight in this analysis\n";
        else std::cout<<"Will not do daisy petal reweight in this analysis\n";
        for (int petal=-1; petal<12 && doDaisy; petal++){
          std::string datapathdaisy = dirs[c] + "/runEventLoopTrackerData_petal_"+std::to_string(petal)+".root";
          std::string mcpathdaisy = dirs[c] + "/runEventLoopTrackerMC_petal_"+std::to_string(petal)+".root";
          std::string migrationpathdaisy = dirs[c] + "/runEventLoopTracker2DMigration_petal_"+std::to_string(petal)+".root";
          std::cout<< "datapathdaisy: "<<datapathdaisy <<std::endl;
          auto dataDaisyFile = TFile::Open(datapathdaisy.c_str(), "READ");
          if (!dataDaisyFile)
          {
            std::cerr << "Failed to open data file " << datapathdaisy.c_str() << ".\n";
            return 2;
          }

          auto mcDaisyFile = TFile::Open(mcpathdaisy.c_str(), "READ");
          if (!mcDaisyFile)
          {
            std::cerr << "Failed to open MC file " << mcpathdaisy.c_str() << ".\n";
            return 3;
          }

          auto migFile = TFile::Open(migrationpathdaisy.c_str(), "READ");
          if (!migFile)
          {
            std::cerr << "Failed to open migration file " << migrationpathdaisy.c_str() << ".\n";
            return 3;
          }

          //std::cout<< "Here1"<<std::endl;

          std::string playlistUsed = util::GetIngredient<TNamed>(*mcDaisyFile, "PlaylistUsed")->GetTitle();
          int filledorempty = util::filledOrEmptyMEPlaylist(playlistUsed);

          int nuoranu = util::nuOrAntiNuMode(playlistUsed);
          int nupdg;
          if (nuoranu==1) nupdg = 14;
          else if (nuoranu==2) nupdg = -14;

          if (pdg!=nupdg)
          {
            std::cout<<"Skipping this set of files because this playlist pdg is " << nupdg << " but you specified " << pdg << std::endl;
            continue;
          }

          //std::cout<< "Here2"<<std::endl;

          auto tmpEffNum = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, (std::string("efficiency_numerator")), prefix);
          auto tmpEffDenom = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, (std::string("efficiency_denominator")), prefix);
          auto tmpMigration = util::GetIngredient<PlotUtils::MnvH2D>(*migFile, (std::string("migration")), prefix);
          auto tmpMigration_reco = util::GetIngredient<PlotUtils::MnvH2D>(*migFile, (std::string("reco")), prefix);
          auto tmpMigration_truth = util::GetIngredient<PlotUtils::MnvH2D>(*migFile, (std::string("truth")), prefix);
          auto tmpFolded = util::GetIngredient<PlotUtils::MnvH2D>(*dataDaisyFile, (std::string("_data_")+prefix));
          auto tmpEffDenom2P2H = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, (std::string("efficiency_denominator_intChannels_2p2h")), prefix);
          auto tmpEffDenomDIS = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, (std::string("efficiency_denominator_intChannels_DIS")), prefix);
          auto tmpEffDenomRES = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, (std::string("efficiency_denominator_intChannels_RES")), prefix);
          auto tmpEffDenomQE = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, (std::string("efficiency_denominator_intChannels_QE")), prefix);
          auto tmpEffDenomOther = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, (std::string("efficiency_denominator_intChannels_Other")), prefix);

          auto tmpBkgWS = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, (std::string("by_BKG_Label_Bkg_Wrong_Sign")), prefix);
          auto tmpBkgNC = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, (std::string("by_BKG_Label_NC_Bkg")), prefix);
          auto tmpBkgOther = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, (std::string("by_BKG_Label_Other")), prefix);

          auto tmpSelectedSignalReco = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, (std::string("selected_signal_reco")), prefix);
          
          auto tmpMCData = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, (std::string("data")), prefix);
          auto tmpMCData2p2h = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, (std::string("intChannels_2p2h")), prefix);
          auto tmpMCDataDIS = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, (std::string("intChannels_DIS")), prefix);
          auto tmpMCDataRES = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, (std::string("intChannels_RES")), prefix);
          auto tmpMCDataQE = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, (std::string("intChannels_QE")), prefix);
          auto tmpMCDataOther = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, (std::string("intChannels_Other")), prefix);

          if (petal>=0)
          {
            //std::cout<< "Here2.1"<<std::endl;
            std::cout<< "petal:"<<petal<<std::endl;
            util::AddHist(*DaisyEffNum[petal], tmpEffNum);
            //std::cout<< "Here2.11"<<std::endl;

            util::AddHist(*DaisyEffDenom[petal], tmpEffDenom);
            util::AddHist(*DaisyMigration[petal], tmpMigration);
            util::AddHist(*DaisyMigrationReco[petal], tmpMigration_reco);
            util::AddHist(*DaisyMigrationTruth[petal], tmpMigration_truth);
            util::AddHist(*DaisyFolded[petal], tmpFolded);
            util::AddHist(*DaisyEffDenom2P2H[petal], tmpEffDenom2P2H);
            util::AddHist(*DaisyEffDenomDIS[petal], tmpEffDenomDIS);
            util::AddHist(*DaisyEffDenomRES[petal], tmpEffDenomRES);
            util::AddHist(*DaisyEffDenomQE[petal], tmpEffDenomQE);
            util::AddHist(*DaisyEffDenomOther[petal], tmpEffDenomOther);

            //std::cout<< "Here2.2"<<std::endl;

            DaisyBackgrounds[petal].push_back(tmpBkgWS->Clone());
            DaisyBackgrounds[petal].push_back(tmpBkgNC->Clone());
            DaisyBackgrounds[petal].push_back(tmpBkgOther->Clone());
          }

          backgrounds.push_back(tmpBkgWS->Clone());
          backgrounds.push_back(tmpBkgNC->Clone());
          backgrounds.push_back(tmpBkgOther->Clone());

          //std::cout<< "Here2.3"<<std::endl;

          if (petal >=-1){
            util::AddHist(*folded, tmpFolded);
            util::AddHist(*migration,tmpMigration);
            util::AddHist(*migration_reco,tmpMigration_reco);
            util::AddHist(*migration_truth, tmpMigration_truth);
            util::AddHist(*effNum, tmpEffNum);
            util::AddHist(*effDenom, tmpEffDenom);
            util::AddHist(*effDenom2P2H, tmpEffDenom2P2H);
            util::AddHist(*effDenomDIS, tmpEffDenomDIS);
            util::AddHist(*effDenomRES, tmpEffDenomRES);
            util::AddHist(*effDenomQE, tmpEffDenomQE);
            util::AddHist(*effDenomOther, tmpEffDenomOther);

            util::AddHist(*BackgroundWrongSign, tmpBkgWS);
            util::AddHist(*BackgroundNC, tmpBkgNC);
            util::AddHist(*BackgroundOther,tmpBkgOther);

            util::AddHist(*SelectedSignalReco,tmpSelectedSignalReco);
                      
            util::AddHist(*MCData, tmpMCData);
            util::AddHist(*MCData2p2h, tmpMCData2p2h);
            util::AddHist(*MCDataDIS,tmpMCDataDIS);
            util::AddHist(*MCDataRES, tmpMCDataRES);
            util::AddHist(*MCDataQE, tmpMCDataQE);
            util::AddHist(*MCDataOther,tmpMCDataOther);
          }
          //std::cout<< "Here2.4"<<std::endl;


          if (petal ==0)
          {
            mcPOT += util::GetIngredient<TParameter<double>>(*mcDaisyFile, "POTUsed")->GetVal();
            double tempDataPOT = util::GetIngredient<TParameter<double>>(*dataDaisyFile, "POTUsed")->GetVal();
            dataPOT+=tempDataPOT;
            playlistPOTpair.push_back(std::make_pair(playlistUsed, tempDataPOT));
            PlotUtils::FluxReweighter frw_temp = PlotUtils::FluxReweighter( nupdg, use_nue_constraint, playlistUsed, PlotUtils::FluxReweighter::gen2thin, PlotUtils::FluxReweighter::g4numiv6, n_flux_universes );
            auto tempIntFlux = frw_temp.GetIntegratedFluxReweighted(nupdg, effDenom, min_energy, max_energy, true)->Clone();
            tempIntFlux->Scale(tempDataPOT);
            util::AddHist(*fluxIntReweighted,tempIntFlux);
          }

          delete tmpEffNum;
          delete tmpEffDenom;
          delete tmpMigration;
          delete tmpFolded;
          delete tmpEffDenom2P2H;
          delete tmpEffDenomDIS;
          delete tmpEffDenomRES;
          delete tmpEffDenomQE;
          delete tmpEffDenomOther;

          delete tmpBkgWS;
          delete tmpBkgNC;
          delete tmpBkgOther;

          delete tmpSelectedSignalReco;
          
          delete tmpMCData;
          delete tmpMCData2p2h;
          delete tmpMCDataDIS;
          delete tmpMCDataRES;
          delete tmpMCDataQE;
          delete tmpMCDataOther;

          dataDaisyFile->Close();
          mcDaisyFile->Close();

          std::cout<< "Here3"<<std::endl;

        }
      }
      doDaisy = true; //Hacky
      //Normalising the integrated flux across different playlists by dataPOT
      fluxIntReweighted->Scale(1/dataPOT);
      auto simEventRate = effDenom->Clone(); // Make a copy for later
      auto simEventRate2P2H = effDenom2P2H->Clone(); // Make a copy for later
      auto simEventRateDIS = effDenomDIS->Clone(); // Make a copy for later
      auto simEventRateRES = effDenomRES->Clone(); // Make a copy for later
      auto simEventRateQE = effDenomQE->Clone(); // Make a copy for later
      auto simEventRateOther = effDenomOther->Clone(); // Make a copy for later
      // There are no error bands in the data, but I need somewhere to put error bands on the results I derive from it.
      folded->AddMissingErrorBandsAndFillWithCV(*migration);
      std::cout<< "Here5"<<std::endl;
      for (int petal=0; petal<12 && doDaisy; petal++) DaisyFolded[petal]->AddMissingErrorBandsAndFillWithCV(*(DaisyMigration[petal]));
      std::cout<< "Here5.1"<<std::endl;

      // Basing my unfolding procedure for a differential cross section on Alex's MINERvA 101 talk at https://minerva-docdb.fnal.gov/cgi-bin/private/RetrieveFile?docid=27438&filename=whatsACrossSection.pdf&version=1

      // TODO: Remove these debugging plots when done
      auto toSubtract = std::accumulate(std::next(backgrounds.begin()), backgrounds.end(), (*backgrounds.begin())->Clone(),
                                        [](auto sum, const auto hist)
                                        {
                                          sum->Add(hist);
                                          return sum;
                                        });
      //Plot(*toSubtract, "BackgroundSum", prefix);
      std::cout<< "Here5.2"<<std::endl;

      auto bkgSubtracted = std::accumulate(backgrounds.begin(), backgrounds.end(), folded->Clone(),
                                            [mcPOT, dataPOT](auto sum, const auto hist)
                                            {
                                              std::cout << "Subtracting " << hist->GetName() << " scaled by " << -dataPOT / mcPOT << " from " << sum->GetName() << "\n";
                                              sum->Add(hist, -dataPOT / mcPOT);
                                              return sum;
                                            });
      //Plot(*bkgSubtracted, "backgroundSubtracted", prefix);
      //std::cout<< "Here6"<<std::endl;
      auto outFile = TFile::Open((prefix + "_crossSection.root").c_str(), "RECREATE");
      if (!outFile)
      {
        std::cerr << "Could not create a file called " << prefix + "_crossSection.root" << ".  Does it already exist?\n";
        return 5;
      }

      bkgSubtracted->Write("backgroundSubtracted");
      // d'Aogstini unfolding
      migration->Write("migration");
      //return 0;
      std::cout<<migration_reco->GetName()<<", "<<migration_truth->GetName()<<", "<<migration->GetName()<< std::endl;
      std::cout<<(migration_reco->GetNbinsX()+2)<<" "<<(migration_reco->GetNbinsY()+2)<<" "<<migration->GetNbinsX()<<" "<<(migration_truth->GetNbinsX()+2)<<" "<<(migration_truth->GetNbinsY()+2) <<" "<< migration->GetNbinsY()<< std::endl;
      std::cout<< "Here5.3"<<std::endl;
      auto unfolded = UnfoldHist(bkgSubtracted, migration, migration_reco, migration_truth, nIterations);
      if (!unfolded)
        throw std::runtime_error(std::string("Failed to unfold ") + folded->GetName() + " using " + migration->GetName());
      //Plot(*unfolded, "unfolded", prefix);
      unfolded->Clone()->Write("unfolded"); // TODO: Seg fault first appears when I uncomment this line
      std::cout << "Survived writing the unfolded histogram.\n"
                << std::endl; // This is evidence that the problem is on the final file Write() and not unfolded->Clone()->Write().
      effNum->Divide(effNum, effDenom); // Only the 2 parameter version of MnvH1D::Divide()
                                        // handles systematics correctly.
      //Plot(*effNum, "efficiency", prefix);
      effNum->Write(("efficiency"));
      unfolded->Divide(unfolded, effNum);
      //Plot(*unfolded, "efficiencyCorrected", prefix);


      PlotUtils::TargetUtils targetInfo;
      double nnucleons = targetInfo.GetTrackerNNucleons(PlotUtils::TargetProp::Tracker::Face, PlotUtils::TargetProp::Tracker::Back, true);
      double nnucleonsData = targetInfo.GetTrackerNNucleons(PlotUtils::TargetProp::Tracker::Face, PlotUtils::TargetProp::Tracker::Back, false);

      PlotUtils::MnvH2D *flux2;
      PlotUtils::MnvH2D *fluxIntegral;
      PlotUtils::MnvH2D *fluxRebinned;

      outFile->cd();
      auto crossSection = normalize(unfolded, fluxIntReweighted, nnucleonsData, dataPOT);
      //Plot(*crossSection, "crossSection", prefix);
      crossSection->Clone()->Write("crossSection");
      simEventRate->Write("simulatedEventRate");
      simEventRate2P2H->Write("simulatedEventRate2P2H");
      simEventRateDIS->Write("simulatedEventRateDIS");
      simEventRateRES->Write("simulatedEventRateRES");
      simEventRateQE->Write("simulatedEventRateQE");
      simEventRateOther->Write("simulatedEventRateOther");
      MCData->Write("mcData");
      MCData2p2h->Write("MCData2p2h");
      MCDataDIS->Write("MCDataDIS");
      MCDataRES->Write("MCDataRES");
      MCDataQE->Write("MCDataQE");
      MCDataOther->Write("MCDataOther");
      folded->Write("folded");
      fluxIntReweighted->Write("fluxIntReweighted");

      // Write a "simulated cross section" to compare to the data I just extracted.
      // If this analysis passed its closure test, this should be the same cross section as
      // what GENIEXSecExtract would produce.
      //fluxIntegral about 1pc off
      //fluxIntReweighted passes
      auto simulatedCrossSection = normalize(simEventRate, fluxIntReweighted, nnucleons, mcPOT);
      //Plot(*simulatedCrossSection, "simulatedCrossSection", prefix);
      simulatedCrossSection->Write("simulatedCrossSection");
      auto simulatedCrossSection2P2H = normalize(simEventRate2P2H, fluxIntReweighted, nnucleons, mcPOT);
      simulatedCrossSection2P2H->Write("simulatedCrossSection2P2H");
      auto simulatedCrossSectionDIS = normalize(simEventRateDIS, fluxIntReweighted, nnucleons, mcPOT);
      simulatedCrossSectionDIS->Write("simulatedCrossSectionDIS");
      auto simulatedCrossSectionRES = normalize(simEventRateRES, fluxIntReweighted, nnucleons, mcPOT);
      simulatedCrossSectionRES->Write("simulatedCrossSectionRES");
      auto simulatedCrossSectionQE = normalize(simEventRateQE, fluxIntReweighted, nnucleons, mcPOT);
      simulatedCrossSectionQE->Write("simulatedCrossSectionQE");
      auto simulatedCrossSectionOther = normalize(simEventRateOther, fluxIntReweighted, nnucleons, mcPOT);
      simulatedCrossSectionOther->Write("simulatedCrossSectionOther");

      outFile->Close();

      //Daisy reweight-----------------------------------------
      //map of petal distributions

      // ---------------------------------------------------------------------
      // Flux reweighter information, get reweighted daisy sum according to a material
      // ---------------------------------------------------------------------
      if (doDaisy)
      {
        auto outFileDaisy = TFile::Open((prefix + "_Daisy_crossSection.root").c_str(), "RECREATE");
        if (!outFileDaisy)
        {
          std::cerr << "Could not create a file called " << prefix + "_crossSection.root" << ".  Does it already exist?\n";
          return 5;
        }

        std::map<int, PlotUtils::MnvH2D*> daisy_petal_hists;
        std::map<int, PlotUtils::MnvH2D*> daisy_petal_hists_eff_denom;
        std::map<int, PlotUtils::MnvH2D*> daisy_petal_hists_eff_denom_2p2h;
        std::map<int, PlotUtils::MnvH2D*> daisy_petal_hists_eff_denom_dis;
        std::map<int, PlotUtils::MnvH2D*> daisy_petal_hists_eff_denom_res;
        std::map<int, PlotUtils::MnvH2D*> daisy_petal_hists_eff_denom_qe;
        std::map<int, PlotUtils::MnvH2D*> daisy_petal_hists_eff_denom_other;



        for (int petal=0; petal<12; petal++){
          std::vector<PlotUtils::MnvH2D*> DaisyBackgroundsTemp = DaisyBackgrounds[petal];
          std::cout << "Petal " <<petal <<std::endl;
          auto toSubtractDaisy = std::accumulate(std::next(DaisyBackgroundsTemp.begin()), DaisyBackgroundsTemp.end(), (*DaisyBackgroundsTemp.begin())->Clone(),
                                          [](auto sum, const auto hist)
                                          {
                                            sum->Add(hist);
                                            return sum;
                                          });
          //Plot(*toSubtract, "BackgroundSum", prefix);
          outFileDaisy->cd();
          toSubtractDaisy->Write((prefix+"_toSubtractDaisy_"+petal));
          auto bkgSubtractedDaisy = std::accumulate(DaisyBackgrounds[petal].begin(), DaisyBackgrounds[petal].end(), DaisyFolded[petal]->Clone(),
                                              [mcPOT, dataPOT](auto sum, const auto hist)
                                              {
                                                std::cout << "Subtracting " << hist->GetName() << " scaled by " << -dataPOT/mcPOT << " from " << sum->GetName() << "\n";
                                                sum->Add(hist, -dataPOT/mcPOT);
                                                return sum;
                                              });
          outFileDaisy->cd();
          DaisyEffNum[petal]->Write((prefix+"_DaisyEffNum_"+petal));
          outFileDaisy->cd();
          DaisyEffDenom[petal]->Write((prefix+"_DaisyEffDenom_"+petal));
          outFileDaisy->cd();
          DaisyFolded[petal]->Write((prefix+"_DaisyFolded_"+petal));
          outFileDaisy->cd();
          bkgSubtractedDaisy->Write((prefix+"_bkgSubtractedDaisy_"+petal));
          auto unfoldedDaisy = UnfoldHist(bkgSubtractedDaisy, DaisyMigration[petal].get(), DaisyMigrationReco[petal].get(), DaisyMigrationTruth[petal].get(), nIterations);
          outFileDaisy->cd();
          unfoldedDaisy->Write((prefix+"_unfoldedDaisy_"+petal));
          if(!unfoldedDaisy) throw std::runtime_error(std::string("Failed to unfold ") + DaisyFolded[petal]->GetName() + " using " + DaisyMigration[petal]->GetName());

          DaisyEffNum[petal]->Divide(DaisyEffNum[petal].get(),DaisyEffDenom[petal].get());
          DaisyEffNum[petal]->Write((prefix+"_DaisyEfficiency_"+petal));
          unfoldedDaisy->Divide(unfoldedDaisy, DaisyEffNum[petal].get());
          unfoldedDaisy->Write((prefix+"_unfoldedDaisyEffCorrected_"+petal));
          daisy_petal_hists[petal]=unfoldedDaisy->Clone();

          daisy_petal_hists_eff_denom[petal]= DaisyEffDenom[petal].get();
          daisy_petal_hists_eff_denom_2p2h[petal]= DaisyEffDenom2P2H[petal].get();
          daisy_petal_hists_eff_denom_dis[petal]= DaisyEffDenomDIS[petal].get();
          daisy_petal_hists_eff_denom_res[petal]= DaisyEffDenomRES[petal].get();
          daisy_petal_hists_eff_denom_qe[petal]= DaisyEffDenomQE[petal].get();
          daisy_petal_hists_eff_denom_other[petal]= DaisyEffDenomOther[petal].get();

        }
        //std::unique_ptr<PlotUtils::MnvH1D> fluxIntegralCtmp = std::unique_ptr<PlotUtils::MnvH1D>(new PlotUtils::MnvH1D());
        PlotUtils::MnvH2D *fluxIntegralC = new PlotUtils::MnvH2D();
        PlotUtils::MnvH2D *fluxIntegralFe = new PlotUtils::MnvH2D();
        PlotUtils::MnvH2D *fluxIntegralPb = new PlotUtils::MnvH2D(); 
        PlotUtils::MnvH2D *DaisyCorrectedC = new PlotUtils::MnvH2D();
        PlotUtils::MnvH2D *DaisyCorrectedFe = new PlotUtils::MnvH2D();
        PlotUtils::MnvH2D *DaisyCorrectedPb = new PlotUtils::MnvH2D(); 
        PlotUtils::MnvH2D *DaisyCorrectedSimC = new PlotUtils::MnvH2D();
        PlotUtils::MnvH2D *DaisyCorrectedSimFe = new PlotUtils::MnvH2D();
        PlotUtils::MnvH2D *DaisyCorrectedSimPb = new PlotUtils::MnvH2D(); 
        PlotUtils::MnvH2D *DaisyCorrectedSim2p2hC = new PlotUtils::MnvH2D();
        PlotUtils::MnvH2D *DaisyCorrectedSim2p2hFe = new PlotUtils::MnvH2D();
        PlotUtils::MnvH2D *DaisyCorrectedSim2p2hPb = new PlotUtils::MnvH2D(); 
        PlotUtils::MnvH2D *DaisyCorrectedSimDISC = new PlotUtils::MnvH2D();
        PlotUtils::MnvH2D *DaisyCorrectedSimDISFe = new PlotUtils::MnvH2D();
        PlotUtils::MnvH2D *DaisyCorrectedSimDISPb = new PlotUtils::MnvH2D(); 
        PlotUtils::MnvH2D *DaisyCorrectedSimRESC = new PlotUtils::MnvH2D();
        PlotUtils::MnvH2D *DaisyCorrectedSimRESFe = new PlotUtils::MnvH2D();
        PlotUtils::MnvH2D *DaisyCorrectedSimRESPb = new PlotUtils::MnvH2D(); 
        PlotUtils::MnvH2D *DaisyCorrectedSimQEC = new PlotUtils::MnvH2D();
        PlotUtils::MnvH2D *DaisyCorrectedSimQEFe = new PlotUtils::MnvH2D();
        PlotUtils::MnvH2D *DaisyCorrectedSimQEPb = new PlotUtils::MnvH2D();
        PlotUtils::MnvH2D *DaisyCorrectedSimOtherC = new PlotUtils::MnvH2D();
        PlotUtils::MnvH2D *DaisyCorrectedSimOtherFe = new PlotUtils::MnvH2D();
        PlotUtils::MnvH2D *DaisyCorrectedSimOtherPb = new PlotUtils::MnvH2D(); 
        std::cout<<"ABC1\n";
        for (auto pp : playlistPOTpair)
        {
          std::cout<<"ABC1.1\n";
          PlotUtils::FluxReweighter *frw = new PlotUtils::FluxReweighter( pdg, use_nue_constraint, pp.first, PlotUtils::FluxReweighter::gen2thin, PlotUtils::FluxReweighter::g4numiv6, n_flux_universes );

          util::AddHist(*DaisyCorrectedC, frw->GetReweightedDaisySum(14, "carbon", daisy_petal_hists, project_dir ), pp.second);
          util::AddHist(*fluxIntegralC,frw->GetIntegratedTargetFlux(14, "carbon", DaisyCorrectedC, min_energy, max_energy, project_dir), pp.second);
          util::AddHist(*DaisyCorrectedFe, frw->GetReweightedDaisySum(14, "iron", daisy_petal_hists, project_dir ), pp.second);
          util::AddHist(*fluxIntegralFe,frw->GetIntegratedTargetFlux(14, "iron", DaisyCorrectedFe, min_energy, max_energy, project_dir), pp.second);
          util::AddHist(*DaisyCorrectedPb, frw->GetReweightedDaisySum(14, "lead", daisy_petal_hists, project_dir ), pp.second);
          util::AddHist(*fluxIntegralPb,frw->GetIntegratedTargetFlux(14, "lead", DaisyCorrectedPb, min_energy, max_energy, project_dir), pp.second);

          util::AddHist(*DaisyCorrectedSimC, frw->GetReweightedDaisySum(14, "carbon", daisy_petal_hists_eff_denom, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimFe, frw->GetReweightedDaisySum(14, "iron", daisy_petal_hists_eff_denom, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimPb, frw->GetReweightedDaisySum(14, "lead", daisy_petal_hists_eff_denom, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSim2p2hC, frw->GetReweightedDaisySum(14, "carbon", daisy_petal_hists_eff_denom_2p2h, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSim2p2hFe, frw->GetReweightedDaisySum(14, "iron", daisy_petal_hists_eff_denom_2p2h, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSim2p2hPb, frw->GetReweightedDaisySum(14, "lead", daisy_petal_hists_eff_denom_2p2h, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimDISC, frw->GetReweightedDaisySum(14, "carbon", daisy_petal_hists_eff_denom_dis, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimDISFe, frw->GetReweightedDaisySum(14, "iron", daisy_petal_hists_eff_denom_dis, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimDISPb, frw->GetReweightedDaisySum(14, "lead", daisy_petal_hists_eff_denom_dis, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimRESC, frw->GetReweightedDaisySum(14, "carbon", daisy_petal_hists_eff_denom_res, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimRESFe, frw->GetReweightedDaisySum(14, "iron", daisy_petal_hists_eff_denom_res, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimRESPb, frw->GetReweightedDaisySum(14, "lead", daisy_petal_hists_eff_denom_res, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimQEC, frw->GetReweightedDaisySum(14, "carbon", daisy_petal_hists_eff_denom_qe, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimQEFe, frw->GetReweightedDaisySum(14, "iron", daisy_petal_hists_eff_denom_qe, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimQEPb, frw->GetReweightedDaisySum(14, "lead", daisy_petal_hists_eff_denom_qe, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimOtherC, frw->GetReweightedDaisySum(14, "carbon", daisy_petal_hists_eff_denom_other, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimOtherFe, frw->GetReweightedDaisySum(14, "iron", daisy_petal_hists_eff_denom_other, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimOtherPb, frw->GetReweightedDaisySum(14, "lead", daisy_petal_hists_eff_denom_other, project_dir ), pp.second);
        }
        DaisyCorrectedC->Scale(1/dataPOT);
        fluxIntegralC->Scale(1/dataPOT);
        DaisyCorrectedFe->Scale(1/dataPOT);
        fluxIntegralFe->Scale(1/dataPOT);
        DaisyCorrectedPb->Scale(1/dataPOT);
        fluxIntegralPb->Scale(1/dataPOT);

        DaisyCorrectedSimC->Scale(1/dataPOT);
        DaisyCorrectedSimFe->Scale(1/dataPOT);
        DaisyCorrectedSimPb->Scale(1/dataPOT);
        DaisyCorrectedSim2p2hC->Scale(1/dataPOT);
        DaisyCorrectedSim2p2hFe->Scale(1/dataPOT);
        DaisyCorrectedSim2p2hPb->Scale(1/dataPOT);
        DaisyCorrectedSimDISC->Scale(1/dataPOT);
        DaisyCorrectedSimDISFe->Scale(1/dataPOT);
        DaisyCorrectedSimDISPb->Scale(1/dataPOT);
        DaisyCorrectedSimRESC->Scale(1/dataPOT);
        DaisyCorrectedSimRESFe->Scale(1/dataPOT);
        DaisyCorrectedSimRESPb->Scale(1/dataPOT);
        DaisyCorrectedSimQEC->Scale(1/dataPOT);
        DaisyCorrectedSimQEFe->Scale(1/dataPOT);
        DaisyCorrectedSimQEPb->Scale(1/dataPOT);
        DaisyCorrectedSimOtherC->Scale(1/dataPOT);
        DaisyCorrectedSimOtherFe->Scale(1/dataPOT);
        DaisyCorrectedSimOtherPb->Scale(1/dataPOT);

        std::cout<<"ABC2\n";

        /* fluxIntegralC->Scale(1/dataPOT);
        fluxIntegralFe->Scale(1/dataPOT);
        fluxIntegralPb->Scale(1/dataPOT); */
        /* DaisyCorrectedC->Scale(1/dataPOT);
        DaisyCorrectedFe->Scale(1/dataPOT);
        DaisyCorrectedPb->Scale(1/dataPOT); */
        
        /* DaisyCorrected2p2hC->Scale(1/dataPOT);
        DaisyCorrected2p2hFe->Scale(1/dataPOT);
        DaisyCorrected2p2hPb->Scale(1/dataPOT);

        DaisyCorrectedDISC->Scale(1/dataPOT);
        DaisyCorrectedDISFe->Scale(1/dataPOT);
        DaisyCorrectedDISPb->Scale(1/dataPOT);

        DaisyCorrectedRESC->Scale(1/dataPOT);
        DaisyCorrectedRESFe->Scale(1/dataPOT);
        DaisyCorrectedRESPb->Scale(1/dataPOT);

        DaisyCorrectedQEC->Scale(1/dataPOT);
        DaisyCorrectedQEFe->Scale(1/dataPOT);
        DaisyCorrectedQEPb->Scale(1/dataPOT);

        DaisyCorrectedOtherC->Scale(1/dataPOT);
        DaisyCorrectedOtherFe->Scale(1/dataPOT);
        DaisyCorrectedOtherPb->Scale(1/dataPOT); */

        PlotUtils::MnvH2D *crossSectionC, *crossSectionFe, *crossSectionPb;
        
        PlotUtils::MnvH2D *crossSectionSimC, *crossSectionSimFe, *crossSectionSimPb;
        PlotUtils::MnvH2D *crossSectionSimC2p2h, *crossSectionSimFe2p2h, *crossSectionSimPb2p2h;
        PlotUtils::MnvH2D *crossSectionSimCDIS, *crossSectionSimFeDIS, *crossSectionSimPbDIS;
        PlotUtils::MnvH2D *crossSectionSimCRes, *crossSectionSimFeRes, *crossSectionSimPbRes;
        PlotUtils::MnvH2D *crossSectionSimCQE, *crossSectionSimFeQE, *crossSectionSimPbQE;
        PlotUtils::MnvH2D *crossSectionSimCOther, *crossSectionSimFeOther, *crossSectionSimPbOther;
        std::cout<<"ABC3\n";







        //Carbon
        {
          outFileDaisy->cd();
          DaisyCorrectedC->Write((prefix+"_DaisyCorrectedC").c_str());
          fluxIntegralC->Write((prefix+"_fluxIntegralC").c_str());
          DaisyCorrectedSimC->Write((prefix+"_DaisyCorrectedSimC").c_str());
          DaisyCorrectedSim2p2hC->Write((prefix+"_DaisyCorrectedSim2p2hC").c_str());
          crossSectionC = normalize(DaisyCorrectedC, fluxIntegralC, nnucleonsData, dataPOT);
          crossSectionC->Write((prefix+"_C_CrossSection").c_str());
          crossSectionSimC = normalize(DaisyCorrectedSimC, fluxIntegralC, nnucleons, mcPOT);
          crossSectionSimC->Write((prefix+"_C_CrossSectionSimulated").c_str());
          crossSectionSimC2p2h = normalize(DaisyCorrectedSim2p2hC, fluxIntegralC, nnucleons, mcPOT);
          crossSectionSimC2p2h->Write((prefix+"_C_CrossSectionSimulated2p2h").c_str());
          crossSectionSimCDIS = normalize(DaisyCorrectedSimDISC, fluxIntegralC, nnucleons, mcPOT);
          crossSectionSimCDIS->Write((prefix+"_C_CrossSectionSimulatedDIS").c_str());
          crossSectionSimCRes = normalize(DaisyCorrectedSimRESC, fluxIntegralC, nnucleons, mcPOT);
          crossSectionSimCRes->Write((prefix+"_C_CrossSectionSimulatedRES").c_str());
          crossSectionSimCQE = normalize(DaisyCorrectedSimQEC, fluxIntegralC, nnucleons, mcPOT);
          crossSectionSimCQE->Write((prefix+"_C_CrossSectionSimulatedQE").c_str());
          crossSectionSimCOther= normalize(DaisyCorrectedSimOtherC, fluxIntegralC, nnucleons, mcPOT);
          crossSectionSimCOther->Write((prefix+"_C_CrossSectionSimulatedOther").c_str());
        }
        std::cout<<"ABC4\n";

        //Iron
        {
          outFileDaisy->cd();
          DaisyCorrectedFe->Write((prefix+"DaisyCorrectedFe").c_str());
          crossSectionFe = normalize(DaisyCorrectedFe, fluxIntegralFe, nnucleonsData, dataPOT);
          crossSectionFe->Write((prefix+"_Fe_CrossSection").c_str());
          crossSectionSimFe = normalize(DaisyCorrectedSimFe, fluxIntegralFe, nnucleons, mcPOT);
          crossSectionSimFe->Write((prefix+"_Fe_CrossSectionSimulated").c_str());
          crossSectionSimFe2p2h = normalize(DaisyCorrectedSim2p2hFe, fluxIntegralFe, nnucleons, mcPOT);
          crossSectionSimFe2p2h->Write((prefix+"_Fe_CrossSectionSimulated2p2h").c_str());
          crossSectionSimFeDIS = normalize(DaisyCorrectedSimDISFe, fluxIntegralFe, nnucleons, mcPOT);
          crossSectionSimFeDIS->Write((prefix+"_Fe_CrossSectionSimulatedDIS").c_str());
          crossSectionSimFeRes = normalize(DaisyCorrectedSimRESFe, fluxIntegralFe, nnucleons, mcPOT);
          crossSectionSimFeRes->Write((prefix+"_Fe_CrossSectionSimulatedRES").c_str());
          crossSectionSimFeQE = normalize(DaisyCorrectedSimQEFe, fluxIntegralFe, nnucleons, mcPOT);
          crossSectionSimFeQE->Write((prefix+"_Fe_CrossSectionSimulatedQE").c_str());
          crossSectionSimFeOther= normalize(DaisyCorrectedSimOtherFe, fluxIntegralFe, nnucleons, mcPOT);
          crossSectionSimFeOther->Write((prefix+"_Fe_CrossSectionSimulatedOther").c_str());
        }

        //Lead
        {
          outFileDaisy->cd();
          DaisyCorrectedPb->Write((prefix+"DaisyCorrectedPb").c_str());
          crossSectionPb = normalize(DaisyCorrectedPb, fluxIntegralPb, nnucleonsData, dataPOT);
          crossSectionPb->Write((prefix+"_Pb_CrossSection").c_str());
          crossSectionSimPb = normalize(DaisyCorrectedSimPb, fluxIntegralPb, nnucleons, mcPOT);
          crossSectionSimPb->Write((prefix+"_Pb_CrossSectionSimulated").c_str());
          crossSectionSimPb2p2h = normalize(DaisyCorrectedSim2p2hPb, fluxIntegralPb, nnucleons, mcPOT);
          crossSectionSimPb2p2h->Write((prefix+"_Pb_CrossSectionSimulated2p2h").c_str());
          crossSectionSimPbDIS = normalize(DaisyCorrectedSimDISPb, fluxIntegralPb, nnucleons, mcPOT);
          crossSectionSimPbDIS->Write((prefix+"_Pb_CrossSectionSimulatedDIS").c_str());
          crossSectionSimPbRes = normalize(DaisyCorrectedSimRESPb, fluxIntegralPb, nnucleons, mcPOT);
          crossSectionSimPbRes->Write((prefix+"_Pb_CrossSectionSimulatedRES").c_str());
          crossSectionSimPbQE = normalize(DaisyCorrectedSimQEPb, fluxIntegralPb, nnucleons, mcPOT);
          crossSectionSimPbQE->Write((prefix+"_Pb_CrossSectionSimulatedQE").c_str());
          crossSectionSimPbOther= normalize(DaisyCorrectedSimOtherPb, fluxIntegralPb, nnucleons, mcPOT);
          crossSectionSimPbOther->Write((prefix+"_Pb_CrossSectionSimulatedOther").c_str());
        }

        outFileDaisy->Close();


        //--------------------------------------------------------------------------------------
        //Plotting
        //--------------------------------------------------------------------------------------

        { //Cross-Section - Transverse
          TCanvas *canvas = new TCanvas("canvas", "canvas", 1500, 1000);
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
          plotter.axis_minimum=0.001;
          plotter.draw_normalized_to_bin_width=false;
          //Plotting
          TPad *PanelPad = new TPad("pad1","",0.02,0.05,1,0.95);
          PanelPad->Draw();
          PanelPad->Divide(4, 4, 0, 0);
          PlotUtils::MnvH2D* temp2p2h = simulatedCrossSection2P2H->Clone();
          PlotUtils::MnvH2D* tempdis = simulatedCrossSectionDIS->Clone();
          PlotUtils::MnvH2D* tempres = simulatedCrossSectionRES->Clone();
          PlotUtils::MnvH2D* tempqe = simulatedCrossSectionQE->Clone();
          PlotUtils::MnvH2D* tempother = simulatedCrossSectionOther->Clone();
          PlotUtils::MnvH2D* tempsimxsec = simulatedCrossSection->Clone();
          PlotUtils::MnvH2D* tempxsec = crossSection->Clone();
          //TH2D* tempxsec = tempsimxsec->Clone();

          temp2p2h->Scale(1e39);
          tempdis->Scale(1e39);
          tempres->Scale(1e39);
          tempqe->Scale(1e39);
          tempother->Scale(1e39);
          tempxsec->Scale(1e39);
          tempsimxsec->Scale(1e39);
          int nbinsY = tempxsec->GetNbinsY();

          int nbinsX = tempxsec->GetNbinsX();



          PlotUtils::MnvH1D* hist;
          PlotUtils::MnvH1D* histmcproj;
          PlotUtils::MnvH1D* histmcbkgs;
          PlotUtils::MnvH1D* DISHistProj;
          PlotUtils::MnvH1D* RESHistProj;
          PlotUtils::MnvH1D* HistProj2p2h;
          PlotUtils::MnvH1D* QEHistProj;
          PlotUtils::MnvH1D* OtherHistProj;
          PlotUtils::MnvH1D* NCProj;

          //Autos-scaling
          double maxVal = 0 ;
          for (int i = 0; i < nbinsY; i++){
              TH1D* tempdataproj = tempxsec->ProjectionX("a", i+1, i+1);
              tempdataproj->Scale(1);
              TH1D* tempmcproj = tempsimxsec->ProjectionX("b", i+1, i+1);
              tempmcproj->Scale(dataPOT/mcPOT);
              Int_t datamaxbin = tempdataproj->GetMaximumBin();
              double datamax = tempdataproj->GetBinContent(datamaxbin);
              Int_t mcmaxbin = tempmcproj->GetMaximumBin();
              double mcmax = tempmcproj->GetBinContent(mcmaxbin);
              double overallmax = (datamax > mcmax) ? datamax : mcmax;
              if (overallmax > maxVal) maxVal = overallmax;
          }
          std::cout<<"maxVal: " << maxVal <<std::endl;
          for (int i = 0; i < nbinsX; i++){
              std::cout<<"On pad " << i << "\n";
              PanelPad->cd(i+1);
              hist = tempxsec->ProjectionX("c", i+1, i+1);
              histmcproj = tempsimxsec->ProjectionX("d", i+1, i+1);
              DISHistProj = tempdis->ProjectionX("pxdis", i+1, i+1);
              RESHistProj = tempres->ProjectionX("pxres", i+1, i+1);
              QEHistProj = tempqe->ProjectionX("pxqe", i+1, i+1);
              HistProj2p2h = temp2p2h->ProjectionX("px2p2h", i+1, i+1);
              OtherHistProj = tempother->ProjectionX("pxother", i+1, i+1);
              TH1D* inverseProjection = tempxsec->ProjectionY("py", 1, 1); //used for auto-detecting subplot limits

              std::cout<<"Data/MC pot << " << dataPOT/mcPOT <<std::endl;
              Int_t datamaxbin = hist->GetMaximumBin();
              double datamax = hist->GetBinContent(datamaxbin);
              Int_t mcmaxbin = histmcproj->GetMaximumBin();
              double mcmax = histmcproj->GetBinContent(mcmaxbin);
              std::cout<<"data max " << datamax << std::endl;
              std::cout<<"mc max " << mcmax << std::endl;
              std::cout<<"mc max bin" << mcmaxbin << std::endl;
              double overallmax = (datamax > mcmax) ? datamax : mcmax;
              double scalefactor = (std::round((maxVal*5) / overallmax))/5; //round to nearest 0.2
              std::cout<<"scale factor " << scalefactor << std::endl;
              std::cout<<"max " << overallmax*scalefactor << std::endl;

              hist->Scale(scalefactor);
              std::cout<<"Here\n";
              histmcproj->Scale(scalefactor);

              DISHistProj->Scale(scalefactor);
              RESHistProj->Scale(scalefactor);
              QEHistProj->Scale(scalefactor);
              HistProj2p2h->Scale(scalefactor);
              OtherHistProj->Scale(scalefactor);
              gPad->Update();

              hist->SetLineColor(kBlack);
              hist->SetMarkerStyle(20);

              if (i != 12) histmcproj->GetXaxis()->ChangeLabel(1, -1, 0);
              else histmcproj->GetXaxis()->ChangeLabel(1, -1, -1, 31);
              histmcproj->GetXaxis()->ChangeLabel(6, -1, 0);
              histmcproj->GetXaxis()->SetRangeUser(0., 2.5);
              histmcproj->GetXaxis()->SetLabelFont(54);
              histmcproj->GetYaxis()->SetLabelFont(54);
              histmcproj->GetXaxis()->SetLabelSize(25);
              histmcproj->GetYaxis()->SetLabelSize(25);
              histmcproj->GetXaxis()->SetNdivisions(505);
              //histmcproj->GetYaxis()->LabelsOption("v");
              histmcproj->GetYaxis()->SetNdivisions(505);

              histmcproj->SetLineColor(kRed);
              histmcproj->SetLineWidth(3);
              DISHistProj->SetLineColor(6);
              DISHistProj->SetLineWidth(3);
              RESHistProj->SetLineColor(kOrange);
              RESHistProj->SetLineWidth(3);
              QEHistProj->SetLineColor(kBlue);
              QEHistProj->SetLineWidth(3);
              HistProj2p2h->SetLineColor(kGreen);
              HistProj2p2h->SetLineWidth(3);
              OtherHistProj->SetLineColor(kMagenta);
              OtherHistProj->SetLineWidth(3);
              TObjArray* mcArr = new TObjArray();
              mcArr->Add(histmcproj);
              mcArr->Add(DISHistProj);
              mcArr->Add(RESHistProj);
              mcArr->Add(QEHistProj);
              mcArr->Add(HistProj2p2h);
              mcArr->Add(OtherHistProj);
                
              plotter.mc_error_color=kRed;
              plotter.mc_error_style=3002;
              plotter.DrawDataMCVariations(hist, mcArr, 1.0, "N", true, true, true, false, false, true);
              if (i == (0)) //And check for number of panels/y bins
              {
                  canvas->cd();
                  TLegend *legend = new TLegend(0.55,0.1,0.9,0.25);
                  legend->SetBorderSize(0);
                  legend->SetNColumns(4);
                  legend->AddEntry(hist,"Data","lep");
                  legend->AddEntry(histmcproj,"MINERvA Tune v4.3.0","l");
                  legend->AddEntry(DISHistProj,"DIS","l");
                  legend->AddEntry(RESHistProj,"RES","l");
                  legend->AddEntry(QEHistProj,"QE","l");
                  legend->AddEntry(HistProj2p2h,"2p2h","l");
                  legend->AddEntry(OtherHistProj,"Other","l");
                  legend->Draw();
                  PanelPad->cd(i+1);
              }
                      
              //gROOT->SetSelectedPad(PanelPad->cd(i+1));

              std::stringstream lowerlim, upperlim, multip;
              lowerlim << std::fixed << std::setprecision(3) << inverseProjection->GetBinLowEdge(i+1);
              upperlim << std::fixed << std::setprecision(3) << inverseProjection->GetBinLowEdge(i+1)+inverseProjection->GetBinWidth(i+1);
              multip << std::fixed << std::setprecision(1) << scalefactor;
              std::string lowerlimit = lowerlim.str();
              std::string upperlimit = upperlim.str();
              std::string multiplier = std::string("#bf{#times ") +multip.str()+std::string("}");
              std::string textStr = std::string("#bf{")+lowerlimit+std::string(" < p_{Z} < ") + upperlimit + std::string("}");

              std::cout<<"X2: " << hist->GetXaxis()->GetXmax() << "\n";
              std::cout<<"Y2: " << hist->GetZaxis()->GetXmax() << "\n";
              TLatex limLabel;
              limLabel.SetTextSize(20);
              limLabel.SetTextAlign(33);
              limLabel.SetTextFont(43);
              limLabel.DrawLatexNDC(.95,.9,textStr.c_str());
              if (scalefactor!=1)
              {
                  TLatex multiplierlabel;
                  multiplierlabel.SetTextSize(20);
                  multiplierlabel.SetTextAlign(33);  //align at top
                  multiplierlabel.SetTextFont(43);
                  multiplierlabel.DrawLatexNDC(.95,.8,multiplier.c_str());
              }
              //delete t1;
          }

          canvas->cd();
          
          std::string title = "Tracker";

          std::stringstream POTstr;
          POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOT;
          std::string potstr = std::string("POT Used: ")+POTstr.str();
          plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.97, 0.025, 1, 42, 13, 0);

          plotter.WritePreliminary( 0.75, 0.075, 0.035, true);
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 35, 1, 43, 22, 0);
          //plotter.AddPlotLabel("Data/MC" , 0.08, 0.3, 0.03, 1, 42, 33, 90);

          plotter.AddPlotLabel("d^{2}#sigma/dp_{t}dp_{z} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)", 0.020 ,0.51, 30, 1, 43, 22, 90);
          plotter.AddPlotLabel("Muon Transverse Momentum (Gev/c)", 0.51, 0.025, 30, 1, 43);
          canvas->Print((prefix+"_CrossSectionTransverse.png").c_str());
        }
        { //Cross-Section - Transverse
          TCanvas *canvas = new TCanvas("canvas", "canvas", 1500, 1000);
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
          plotter.axis_minimum=0.001;
          plotter.draw_normalized_to_bin_width=false;
          //Plotting
          TPad *PanelPad = new TPad("pad1","",0.02,0.05,1,0.95);
          PanelPad->Draw();
          PanelPad->Divide(4, 4, 0, 0);
          PlotUtils::MnvH2D* temp2p2h = simulatedCrossSection2P2H->Clone();
          PlotUtils::MnvH2D* tempdis = simulatedCrossSectionDIS->Clone();
          PlotUtils::MnvH2D* tempres = simulatedCrossSectionRES->Clone();
          PlotUtils::MnvH2D* tempqe = simulatedCrossSectionQE->Clone();
          PlotUtils::MnvH2D* tempother = simulatedCrossSectionOther->Clone();
          PlotUtils::MnvH2D* tempsimxsec = simulatedCrossSection->Clone();
          PlotUtils::MnvH2D* tempxsec = crossSection->Clone();
          //TH2D* tempxsec = tempsimxsec->Clone();

          temp2p2h->Scale(1e39);
          tempdis->Scale(1e39);
          tempres->Scale(1e39);
          tempqe->Scale(1e39);
          tempother->Scale(1e39);
          tempxsec->Scale(1e39);
          tempsimxsec->Scale(1e39);
          int nbinsY = tempxsec->GetNbinsY();

          int nbinsX = tempxsec->GetNbinsX();



          PlotUtils::MnvH1D* hist;
          PlotUtils::MnvH1D* histmcproj;
          PlotUtils::MnvH1D* histmcbkgs;
          PlotUtils::MnvH1D* DISHistProj;
          PlotUtils::MnvH1D* RESHistProj;
          PlotUtils::MnvH1D* HistProj2p2h;
          PlotUtils::MnvH1D* QEHistProj;
          PlotUtils::MnvH1D* OtherHistProj;
          PlotUtils::MnvH1D* NCProj;

          //Autos-scaling
          double maxVal = 0 ;
          for (int i = 0; i < nbinsY; i++){
              TH1D* tempdataproj = tempxsec->ProjectionY("a", i+1, i+1);
              tempdataproj->Scale(1);
              TH1D* tempmcproj = tempsimxsec->ProjectionY("b", i+1, i+1);
              tempmcproj->Scale(dataPOT/mcPOT);
              Int_t datamaxbin = tempdataproj->GetMaximumBin();
              double datamax = tempdataproj->GetBinContent(datamaxbin);
              Int_t mcmaxbin = tempmcproj->GetMaximumBin();
              double mcmax = tempmcproj->GetBinContent(mcmaxbin);
              double overallmax = (datamax > mcmax) ? datamax : mcmax;
              if (overallmax > maxVal) maxVal = overallmax;
          }
          std::cout<<"maxVal: " << maxVal <<std::endl;
          for (int i = 0; i < nbinsX; i++){
              std::cout<<"On pad " << i << "\n";
              PanelPad->cd(i+1);
              hist = tempxsec->ProjectionY("c", i+1, i+1);
              histmcproj = tempsimxsec->ProjectionY("d", i+1, i+1);
              DISHistProj = tempdis->ProjectionY("pxdis", i+1, i+1);
              RESHistProj = tempres->ProjectionY("pxres", i+1, i+1);
              QEHistProj = tempqe->ProjectionY("pxqe", i+1, i+1);
              HistProj2p2h = temp2p2h->ProjectionY("px2p2h", i+1, i+1);
              OtherHistProj = tempother->ProjectionY("pxother", i+1, i+1);
              TH1D* inverseProjection = tempxsec->ProjectionX("py", 1, 1); //used for auto-detecting subplot limits

              std::cout<<"Data/MC pot << " << dataPOT/mcPOT <<std::endl;
              Int_t datamaxbin = hist->GetMaximumBin();
              double datamax = hist->GetBinContent(datamaxbin);
              Int_t mcmaxbin = histmcproj->GetMaximumBin();
              double mcmax = histmcproj->GetBinContent(mcmaxbin);
              std::cout<<"data max " << datamax << std::endl;
              std::cout<<"mc max " << mcmax << std::endl;
              std::cout<<"mc max bin" << mcmaxbin << std::endl;
              double overallmax = (datamax > mcmax) ? datamax : mcmax;
              double scalefactor = (std::round((maxVal*5) / overallmax))/5; //round to nearest 0.2
              std::cout<<"scale factor " << scalefactor << std::endl;
              std::cout<<"max " << overallmax*scalefactor << std::endl;

              hist->Scale(scalefactor);
              std::cout<<"Here\n";
              histmcproj->Scale(scalefactor);

              DISHistProj->Scale(scalefactor);
              RESHistProj->Scale(scalefactor);
              QEHistProj->Scale(scalefactor);
              HistProj2p2h->Scale(scalefactor);
              OtherHistProj->Scale(scalefactor);
              gPad->Update();

              hist->SetLineColor(kBlack);
              hist->SetMarkerStyle(20);

              if (i != 12) histmcproj->GetXaxis()->ChangeLabel(1, -1, 0);
              else histmcproj->GetXaxis()->ChangeLabel(1, -1, -1, 31);
              histmcproj->GetXaxis()->ChangeLabel(6, -1, 0);
              histmcproj->GetXaxis()->SetRangeUser(1.5, 15);
              histmcproj->GetXaxis()->SetLabelFont(54);
              histmcproj->GetYaxis()->SetLabelFont(54);
              histmcproj->GetXaxis()->SetLabelSize(25);
              histmcproj->GetYaxis()->SetLabelSize(25);
              histmcproj->GetXaxis()->SetNdivisions(502);
              //histmcproj->GetYaxis()->LabelsOption("v");
              histmcproj->GetYaxis()->SetNdivisions(505);

              histmcproj->SetLineColor(kRed);
              histmcproj->SetLineWidth(3);
              DISHistProj->SetLineColor(6);
              DISHistProj->SetLineWidth(3);
              RESHistProj->SetLineColor(kOrange);
              RESHistProj->SetLineWidth(3);
              QEHistProj->SetLineColor(kBlue);
              QEHistProj->SetLineWidth(3);
              HistProj2p2h->SetLineColor(kGreen);
              HistProj2p2h->SetLineWidth(3);
              OtherHistProj->SetLineColor(kMagenta);
              OtherHistProj->SetLineWidth(3);
              TObjArray* mcArr = new TObjArray();
              mcArr->Add(histmcproj);
              mcArr->Add(DISHistProj);
              mcArr->Add(RESHistProj);
              mcArr->Add(QEHistProj);
              mcArr->Add(HistProj2p2h);
              mcArr->Add(OtherHistProj);
                
              plotter.mc_error_color=kRed;
              plotter.mc_error_style=3002;
              plotter.DrawDataMCVariations(hist, mcArr, 1.0, "N", true, true, true, false, false, true);
              if (i == (0)) //And check for number of panels/y bins
              {
                  canvas->cd();
                  TLegend *legend = new TLegend(0.55,0.1,0.9,0.25);
                  legend->SetBorderSize(0);
                  legend->SetNColumns(4);
                  legend->AddEntry(hist,"Data","lep");
                  legend->AddEntry(histmcproj,"MINERvA Tune v4.3.0","l");
                  legend->AddEntry(DISHistProj,"DIS","l");
                  legend->AddEntry(RESHistProj,"RES","l");
                  legend->AddEntry(QEHistProj,"QE","l");
                  legend->AddEntry(HistProj2p2h,"2p2h","l");
                  legend->AddEntry(OtherHistProj,"Other","l");
                  legend->Draw();
                  PanelPad->cd(i+1);
              }
                      
              //gROOT->SetSelectedPad(PanelPad->cd(i+1));

              std::stringstream lowerlim, upperlim, multip;
              lowerlim << std::fixed << std::setprecision(3) << inverseProjection->GetBinLowEdge(i+1);
              upperlim << std::fixed << std::setprecision(3) << inverseProjection->GetBinLowEdge(i+1)+inverseProjection->GetBinWidth(i+1);
              multip << std::fixed << std::setprecision(1) << scalefactor;
              std::string lowerlimit = lowerlim.str();
              std::string upperlimit = upperlim.str();
              std::string multiplier = std::string("#bf{#times ") +multip.str()+std::string("}");
              std::string textStr = std::string("#bf{")+lowerlimit+std::string(" < p_{t} < ") + upperlimit + std::string("}");

              std::cout<<"X2: " << hist->GetXaxis()->GetXmax() << "\n";
              std::cout<<"Y2: " << hist->GetZaxis()->GetXmax() << "\n";
              TLatex limLabel;
              limLabel.SetTextSize(20);
              limLabel.SetTextAlign(33);
              limLabel.SetTextFont(43);
              limLabel.DrawLatexNDC(.95,.9,textStr.c_str());
              if (scalefactor!=1)
              {
                  TLatex multiplierlabel;
                  multiplierlabel.SetTextSize(20);
                  multiplierlabel.SetTextAlign(33);  //align at top
                  multiplierlabel.SetTextFont(43);
                  multiplierlabel.DrawLatexNDC(.95,.8,multiplier.c_str());
              }
              //delete t1;
          }

          canvas->cd();
          
          std::string title = "Tracker";

          std::stringstream POTstr;
          POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOT;
          std::string potstr = std::string("POT Used: ")+POTstr.str();
          plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.97, 0.025, 1, 42, 13, 0);

          plotter.WritePreliminary( 0.75, 0.075, 0.035, true);
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 35, 1, 43, 22, 0);
          //plotter.AddPlotLabel("Data/MC" , 0.08, 0.3, 0.03, 1, 42, 33, 90);

          plotter.AddPlotLabel("d^{2}#sigma/dp_{t}dp_{z} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)", 0.020 ,0.51, 30, 1, 43, 22, 90);
          plotter.AddPlotLabel("Muon Longitudinal Momentum (Gev/c)", 0.51, 0.025, 30, 1, 43);
          canvas->Print((prefix+"_CrossSectionLongitudinal.png").c_str());
        }





























        { //Data MC ratio- Transverse
          TCanvas *canvas = new TCanvas("canvas", "canvas", 1500, 1000);
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
          plotter.axis_minimum=0.001;
          plotter.draw_normalized_to_bin_width=false;
          //Plotting
          TPad *PanelPad = new TPad("pad1","",0.02,0.05,1,0.95);
          PanelPad->Draw();
          PanelPad->Divide(4, 4, 0, 0);
          PlotUtils::MnvH2D* temp2p2h = simulatedCrossSection2P2H->Clone();
          PlotUtils::MnvH2D* tempdis = simulatedCrossSectionDIS->Clone();
          PlotUtils::MnvH2D* tempres = simulatedCrossSectionRES->Clone();
          PlotUtils::MnvH2D* tempqe = simulatedCrossSectionQE->Clone();
          PlotUtils::MnvH2D* tempother = simulatedCrossSectionOther->Clone();
          PlotUtils::MnvH2D* tempsimxsec = simulatedCrossSection->Clone();
          PlotUtils::MnvH2D* tempxsec = crossSection->Clone();
          //TH2D* tempxsec = tempsimxsec->Clone();

          temp2p2h->Scale(1e39);
          tempdis->Scale(1e39);
          tempres->Scale(1e39);
          tempqe->Scale(1e39);
          tempother->Scale(1e39);
          tempxsec->Scale(1e39);
          tempsimxsec->Scale(1e39);
          int nbinsY = tempxsec->GetNbinsY();

          int nbinsX = tempxsec->GetNbinsX();



          PlotUtils::MnvH1D* hist;
          PlotUtils::MnvH1D* histmcproj;
          PlotUtils::MnvH1D* histmcbkgs;
          PlotUtils::MnvH1D* DISHistProj;
          PlotUtils::MnvH1D* RESHistProj;
          PlotUtils::MnvH1D* HistProj2p2h;
          PlotUtils::MnvH1D* QEHistProj;
          PlotUtils::MnvH1D* OtherHistProj;
          PlotUtils::MnvH1D* NCProj;

          //Autos-scaling
          double maxVal = 0 ;
          for (int i = 0; i < nbinsY; i++){
              TH1D* tempdataproj = tempxsec->ProjectionX("a", i+1, i+1);
              tempdataproj->Scale(1);
              TH1D* tempmcproj = tempsimxsec->ProjectionX("b", i+1, i+1);
              tempmcproj->Scale(dataPOT/mcPOT);
              Int_t datamaxbin = tempdataproj->GetMaximumBin();
              double datamax = tempdataproj->GetBinContent(datamaxbin);
              Int_t mcmaxbin = tempmcproj->GetMaximumBin();
              double mcmax = tempmcproj->GetBinContent(mcmaxbin);
              double overallmax = (datamax > mcmax) ? datamax : mcmax;
              if (overallmax > maxVal) maxVal = overallmax;
          }
          std::cout<<"maxVal: " << maxVal <<std::endl;
          for (int i = 0; i < nbinsX; i++){
              std::cout<<"On pad " << i << "\n";
              PanelPad->cd(i+1);
              hist = tempxsec->ProjectionX("c", i+1, i+1);
              histmcproj = tempsimxsec->ProjectionX("d", i+1, i+1);
              DISHistProj = tempdis->ProjectionX("pxdis", i+1, i+1);
              RESHistProj = tempres->ProjectionX("pxres", i+1, i+1);
              QEHistProj = tempqe->ProjectionX("pxqe", i+1, i+1);
              HistProj2p2h = temp2p2h->ProjectionX("px2p2h", i+1, i+1);
              OtherHistProj = tempother->ProjectionX("pxother", i+1, i+1);
              TH1D* inverseProjection = tempxsec->ProjectionY("py", 1, 1); //used for auto-detecting subplot limits

              std::cout<<"Data/MC pot << " << dataPOT/mcPOT <<std::endl;
              Int_t datamaxbin = hist->GetMaximumBin();
              double datamax = hist->GetBinContent(datamaxbin);
              Int_t mcmaxbin = histmcproj->GetMaximumBin();
              double mcmax = histmcproj->GetBinContent(mcmaxbin);
              std::cout<<"data max " << datamax << std::endl;
              std::cout<<"mc max " << mcmax << std::endl;
              std::cout<<"mc max bin" << mcmaxbin << std::endl;
              double overallmax = (datamax > mcmax) ? datamax : mcmax;
              double scalefactor = (std::round((maxVal*5) / overallmax))/5; //round to nearest 0.2
              std::cout<<"scale factor " << scalefactor << std::endl;
              std::cout<<"max " << overallmax*scalefactor << std::endl;

              hist->Scale(scalefactor);
              std::cout<<"Here\n";
              histmcproj->Scale(scalefactor);

              DISHistProj->Scale(scalefactor);
              RESHistProj->Scale(scalefactor);
              QEHistProj->Scale(scalefactor);
              HistProj2p2h->Scale(scalefactor);
              OtherHistProj->Scale(scalefactor);
              gPad->Update();

              hist->SetLineColor(kBlack);
              hist->SetMarkerStyle(20);

              if (i != 12) histmcproj->GetXaxis()->ChangeLabel(1, -1, 0);
              else histmcproj->GetXaxis()->ChangeLabel(1, -1, -1, 31);
              histmcproj->GetXaxis()->ChangeLabel(6, -1, 0);
              histmcproj->GetXaxis()->SetRangeUser(0., 2.5);
              histmcproj->GetXaxis()->SetLabelFont(54);
              histmcproj->GetYaxis()->SetLabelFont(54);
              histmcproj->GetXaxis()->SetLabelSize(25);
              histmcproj->GetYaxis()->SetLabelSize(25);
              histmcproj->GetXaxis()->SetNdivisions(505);
              //histmcproj->GetYaxis()->LabelsOption("v");
              histmcproj->GetYaxis()->SetNdivisions(505);

              histmcproj->SetLineColor(kRed);
              histmcproj->SetLineWidth(3);
              DISHistProj->SetLineColor(6);
              DISHistProj->SetLineWidth(3);
              RESHistProj->SetLineColor(kOrange);
              RESHistProj->SetLineWidth(3);
              QEHistProj->SetLineColor(kBlue);
              QEHistProj->SetLineWidth(3);
              HistProj2p2h->SetLineColor(kGreen);
              HistProj2p2h->SetLineWidth(3);
              OtherHistProj->SetLineColor(kMagenta);
              OtherHistProj->SetLineWidth(3);
              TObjArray* mcArr = new TObjArray();
              mcArr->Add(histmcproj);
              mcArr->Add(DISHistProj);
              mcArr->Add(RESHistProj);
              mcArr->Add(QEHistProj);
              mcArr->Add(HistProj2p2h);
              mcArr->Add(OtherHistProj);
                
              plotter.mc_error_color=kRed;
              plotter.mc_error_style=3002;
              plotter.DrawDataMCRatio(hist, histmcproj, 1.0, true, true, 0.5, 1.5);
              if (i == (0)) //And check for number of panels/y bins
              {
                  canvas->cd();
                  TLegend *legend = new TLegend(0.55,0.1,0.9,0.25);
                  legend->SetBorderSize(0);
                  legend->SetNColumns(4);
                  legend->AddEntry(hist,"Data","lep");
                  legend->AddEntry(histmcproj,"MINERvA Tune v4.3.0","l");
                  legend->AddEntry(DISHistProj,"DIS","l");
                  legend->AddEntry(RESHistProj,"RES","l");
                  legend->AddEntry(QEHistProj,"QE","l");
                  legend->AddEntry(HistProj2p2h,"2p2h","l");
                  legend->AddEntry(OtherHistProj,"Other","l");
                  legend->Draw();
                  PanelPad->cd(i+1);
              }
                      
              //gROOT->SetSelectedPad(PanelPad->cd(i+1));

              std::stringstream lowerlim, upperlim, multip;
              lowerlim << std::fixed << std::setprecision(3) << inverseProjection->GetBinLowEdge(i+1);
              upperlim << std::fixed << std::setprecision(3) << inverseProjection->GetBinLowEdge(i+1)+inverseProjection->GetBinWidth(i+1);
              multip << std::fixed << std::setprecision(1) << scalefactor;
              std::string lowerlimit = lowerlim.str();
              std::string upperlimit = upperlim.str();
              std::string multiplier = std::string("#bf{#times ") +multip.str()+std::string("}");
              std::string textStr = std::string("#bf{")+lowerlimit+std::string(" < p_{Z} < ") + upperlimit + std::string("}");

              std::cout<<"X2: " << hist->GetXaxis()->GetXmax() << "\n";
              std::cout<<"Y2: " << hist->GetZaxis()->GetXmax() << "\n";
              TLatex limLabel;
              limLabel.SetTextSize(20);
              limLabel.SetTextAlign(33);
              limLabel.SetTextFont(43);
              limLabel.DrawLatexNDC(.95,.9,textStr.c_str());
              if (scalefactor!=1)
              {
                  TLatex multiplierlabel;
                  multiplierlabel.SetTextSize(20);
                  multiplierlabel.SetTextAlign(33);  //align at top
                  multiplierlabel.SetTextFont(43);
                  multiplierlabel.DrawLatexNDC(.95,.8,multiplier.c_str());
              }
              //delete t1;
          }

          canvas->cd();
          
          std::string title = "";

          std::stringstream POTstr;
          POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOT;
          std::string potstr = std::string("POT Used: ")+POTstr.str();
          plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.97, 0.025, 1, 42, 13, 0);

          plotter.WritePreliminary( 0.75, 0.075, 0.035, true);
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 35, 1, 43, 22, 0);
          //plotter.AddPlotLabel("Data/MC" , 0.08, 0.3, 0.03, 1, 42, 33, 90);

          plotter.AddPlotLabel("d^{2}#sigma/dp_{t}dp_{z} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)", 0.020 ,0.51, 30, 1, 43, 22, 90);
          plotter.AddPlotLabel("Muon Transverse Momentum (Gev/c)", 0.51, 0.025, 30, 1, 43);
          canvas->Print((prefix+"_DataMCRatioTransverse.png").c_str());
        }
        { //Data MC ratio - Transverse
          TCanvas *canvas = new TCanvas("canvas", "canvas", 1500, 1000);
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
          plotter.axis_minimum=0.001;
          plotter.draw_normalized_to_bin_width=false;
          //Plotting
          TPad *PanelPad = new TPad("pad1","",0.02,0.05,1,0.95);
          PanelPad->Draw();
          PanelPad->Divide(4, 4, 0, 0);
          PlotUtils::MnvH2D* temp2p2h = simulatedCrossSection2P2H->Clone();
          PlotUtils::MnvH2D* tempdis = simulatedCrossSectionDIS->Clone();
          PlotUtils::MnvH2D* tempres = simulatedCrossSectionRES->Clone();
          PlotUtils::MnvH2D* tempqe = simulatedCrossSectionQE->Clone();
          PlotUtils::MnvH2D* tempother = simulatedCrossSectionOther->Clone();
          PlotUtils::MnvH2D* tempsimxsec = simulatedCrossSection->Clone();
          PlotUtils::MnvH2D* tempxsec = crossSection->Clone();
          //TH2D* tempxsec = tempsimxsec->Clone();

          temp2p2h->Scale(1e39);
          tempdis->Scale(1e39);
          tempres->Scale(1e39);
          tempqe->Scale(1e39);
          tempother->Scale(1e39);
          tempxsec->Scale(1e39);
          tempsimxsec->Scale(1e39);
          int nbinsY = tempxsec->GetNbinsY();

          int nbinsX = tempxsec->GetNbinsX();



          PlotUtils::MnvH1D* hist;
          PlotUtils::MnvH1D* histmcproj;
          PlotUtils::MnvH1D* histmcbkgs;
          PlotUtils::MnvH1D* DISHistProj;
          PlotUtils::MnvH1D* RESHistProj;
          PlotUtils::MnvH1D* HistProj2p2h;
          PlotUtils::MnvH1D* QEHistProj;
          PlotUtils::MnvH1D* OtherHistProj;
          PlotUtils::MnvH1D* NCProj;

          //Autos-scaling
          double maxVal = 0 ;
          for (int i = 0; i < nbinsY; i++){
              TH1D* tempdataproj = tempxsec->ProjectionY("a", i+1, i+1);
              tempdataproj->Scale(1);
              TH1D* tempmcproj = tempsimxsec->ProjectionY("b", i+1, i+1);
              tempmcproj->Scale(dataPOT/mcPOT);
              Int_t datamaxbin = tempdataproj->GetMaximumBin();
              double datamax = tempdataproj->GetBinContent(datamaxbin);
              Int_t mcmaxbin = tempmcproj->GetMaximumBin();
              double mcmax = tempmcproj->GetBinContent(mcmaxbin);
              double overallmax = (datamax > mcmax) ? datamax : mcmax;
              if (overallmax > maxVal) maxVal = overallmax;
          }
          std::cout<<"maxVal: " << maxVal <<std::endl;
          for (int i = 0; i < nbinsX; i++){
              std::cout<<"On pad " << i << "\n";
              PanelPad->cd(i+1);
              hist = tempxsec->ProjectionY("c", i+1, i+1);
              histmcproj = tempsimxsec->ProjectionY("d", i+1, i+1);
              DISHistProj = tempdis->ProjectionY("pxdis", i+1, i+1);
              RESHistProj = tempres->ProjectionY("pxres", i+1, i+1);
              QEHistProj = tempqe->ProjectionY("pxqe", i+1, i+1);
              HistProj2p2h = temp2p2h->ProjectionY("px2p2h", i+1, i+1);
              OtherHistProj = tempother->ProjectionY("pxother", i+1, i+1);
              TH1D* inverseProjection = tempxsec->ProjectionX("py", 1, 1); //used for auto-detecting subplot limits

              std::cout<<"Data/MC pot << " << dataPOT/mcPOT <<std::endl;
              Int_t datamaxbin = hist->GetMaximumBin();
              double datamax = hist->GetBinContent(datamaxbin);
              Int_t mcmaxbin = histmcproj->GetMaximumBin();
              double mcmax = histmcproj->GetBinContent(mcmaxbin);
              std::cout<<"data max " << datamax << std::endl;
              std::cout<<"mc max " << mcmax << std::endl;
              std::cout<<"mc max bin" << mcmaxbin << std::endl;
              double overallmax = (datamax > mcmax) ? datamax : mcmax;
              double scalefactor = (std::round((maxVal*5) / overallmax))/5; //round to nearest 0.2
              std::cout<<"scale factor " << scalefactor << std::endl;
              std::cout<<"max " << overallmax*scalefactor << std::endl;

              hist->Scale(scalefactor);
              std::cout<<"Here\n";
              histmcproj->Scale(scalefactor);

              DISHistProj->Scale(scalefactor);
              RESHistProj->Scale(scalefactor);
              QEHistProj->Scale(scalefactor);
              HistProj2p2h->Scale(scalefactor);
              OtherHistProj->Scale(scalefactor);
              gPad->Update();

              hist->SetLineColor(kBlack);
              hist->SetMarkerStyle(20);

              if (i != 12) histmcproj->GetXaxis()->ChangeLabel(1, -1, 0);
              else histmcproj->GetXaxis()->ChangeLabel(1, -1, -1, 31);
              histmcproj->GetXaxis()->ChangeLabel(6, -1, 0);
              histmcproj->GetXaxis()->SetRangeUser(1.5, 15);
              histmcproj->GetXaxis()->SetLabelFont(54);
              histmcproj->GetYaxis()->SetLabelFont(54);
              histmcproj->GetXaxis()->SetLabelSize(25);
              histmcproj->GetYaxis()->SetLabelSize(25);
              histmcproj->GetXaxis()->SetNdivisions(502);
              //histmcproj->GetYaxis()->LabelsOption("v");
              histmcproj->GetYaxis()->SetNdivisions(505);

              histmcproj->SetLineColor(kRed);
              histmcproj->SetLineWidth(3);
              DISHistProj->SetLineColor(6);
              DISHistProj->SetLineWidth(3);
              RESHistProj->SetLineColor(kOrange);
              RESHistProj->SetLineWidth(3);
              QEHistProj->SetLineColor(kBlue);
              QEHistProj->SetLineWidth(3);
              HistProj2p2h->SetLineColor(kGreen);
              HistProj2p2h->SetLineWidth(3);
              OtherHistProj->SetLineColor(kMagenta);
              OtherHistProj->SetLineWidth(3);
              TObjArray* mcArr = new TObjArray();
              mcArr->Add(histmcproj);
              mcArr->Add(DISHistProj);
              mcArr->Add(RESHistProj);
              mcArr->Add(QEHistProj);
              mcArr->Add(HistProj2p2h);
              mcArr->Add(OtherHistProj);
                
              plotter.mc_error_color=kRed;
              plotter.mc_error_style=3002;

              plotter.DrawDataMCRatio(hist, histmcproj, 1.0, true, true, 0.5, 1.5);
              if (i == (0)) //And check for number of panels/y bins
              {
                  canvas->cd();
                  TLegend *legend = new TLegend(0.55,0.1,0.9,0.25);
                  legend->SetBorderSize(0);
                  legend->SetNColumns(4);
                  legend->AddEntry(hist,"Data","lep");
                  legend->AddEntry(histmcproj,"MINERvA Tune v4.3.0","l");
                  legend->AddEntry(DISHistProj,"DIS","l");
                  legend->AddEntry(RESHistProj,"RES","l");
                  legend->AddEntry(QEHistProj,"QE","l");
                  legend->AddEntry(HistProj2p2h,"2p2h","l");
                  legend->AddEntry(OtherHistProj,"Other","l");
                  legend->Draw();
                  PanelPad->cd(i+1);
              }
                      
              //gROOT->SetSelectedPad(PanelPad->cd(i+1));

              std::stringstream lowerlim, upperlim, multip;
              lowerlim << std::fixed << std::setprecision(3) << inverseProjection->GetBinLowEdge(i+1);
              upperlim << std::fixed << std::setprecision(3) << inverseProjection->GetBinLowEdge(i+1)+inverseProjection->GetBinWidth(i+1);
              multip << std::fixed << std::setprecision(1) << scalefactor;
              std::string lowerlimit = lowerlim.str();
              std::string upperlimit = upperlim.str();
              std::string multiplier = std::string("#bf{#times ") +multip.str()+std::string("}");
              std::string textStr = std::string("#bf{")+lowerlimit+std::string(" < p_{t} < ") + upperlimit + std::string("}");

              std::cout<<"X2: " << hist->GetXaxis()->GetXmax() << "\n";
              std::cout<<"Y2: " << hist->GetZaxis()->GetXmax() << "\n";
              TLatex limLabel;
              limLabel.SetTextSize(20);
              limLabel.SetTextAlign(33);
              limLabel.SetTextFont(43);
              limLabel.DrawLatexNDC(.95,.9,textStr.c_str());
              if (scalefactor!=1)
              {
                  TLatex multiplierlabel;
                  multiplierlabel.SetTextSize(20);
                  multiplierlabel.SetTextAlign(33);  //align at top
                  multiplierlabel.SetTextFont(43);
                  multiplierlabel.DrawLatexNDC(.95,.8,multiplier.c_str());
              }
              //delete t1;
          }

          canvas->cd();
          
          std::string title = "";

          std::stringstream POTstr;
          POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOT;
          std::string potstr = std::string("POT Used: ")+POTstr.str();
          plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.97, 0.025, 1, 42, 13, 0);

          plotter.WritePreliminary( 0.75, 0.075, 0.035, true);
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 35, 1, 43, 22, 0);
          //plotter.AddPlotLabel("Data/MC" , 0.08, 0.3, 0.03, 1, 42, 33, 90);

          plotter.AddPlotLabel("d^{2}#sigma/dp_{t}dp_{z} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)", 0.020 ,0.51, 30, 1, 43, 22, 90);
          plotter.AddPlotLabel("Muon Longitudinal Momentum (Gev/c)", 0.51, 0.025, 30, 1, 43);
          canvas->Print((prefix+"_DataMCRatioLongitudinal.png").c_str());
        }


































        { //Event Selection - Transverse
          TCanvas *canvas = new TCanvas("canvas", "canvas", 1500, 1000);
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
          plotter.axis_minimum=0.001;
          plotter.draw_normalized_to_bin_width=true;
          //Plotting
          TPad *PanelPad = new TPad("pad1","",0.02,0.05,1,0.95);
          PanelPad->Draw();
          PanelPad->Divide(4, 4, 0, 0);
          PlotUtils::MnvH2D* temp2p2h = simEventRate2P2H->Clone();
          PlotUtils::MnvH2D* tempdis = simEventRateDIS->Clone();
          PlotUtils::MnvH2D* tempres = simEventRateRES->Clone();
          PlotUtils::MnvH2D* tempqe = simEventRateQE->Clone();
          PlotUtils::MnvH2D* tempother = simEventRateOther->Clone();
          PlotUtils::MnvH2D* tempsimevrate = simEventRate->Clone();
          PlotUtils::MnvH2D* tempevrate = unfolded->Clone();
          //TH2D* tempevrate = tempsimevrate->Clone();

          temp2p2h->Scale(1e-3);
          tempdis->Scale(1e-3);
          tempres->Scale(1e-3);
          tempqe->Scale(1e-3);
          tempother->Scale(1e-3);
          tempsimevrate->Scale(1e-3);
          tempevrate->Scale(1e-3);

          int nbinsY = tempevrate->GetNbinsY();

          int nbinsX = tempevrate->GetNbinsX();



          PlotUtils::MnvH1D* hist;
          PlotUtils::MnvH1D* histmcproj;
          PlotUtils::MnvH1D* histmcbkgs;
          PlotUtils::MnvH1D* DISHistProj;
          PlotUtils::MnvH1D* RESHistProj;
          PlotUtils::MnvH1D* HistProj2p2h;
          PlotUtils::MnvH1D* QEHistProj;
          PlotUtils::MnvH1D* OtherHistProj;
          PlotUtils::MnvH1D* NCProj;

          //Autos-scaling
          double maxVal = 0 ;
          for (int i = 0; i < nbinsY; i++){
              TH1D* tempdataproj = tempevrate->ProjectionX("a", i+1, i+1);
              tempdataproj->Scale(1);
              TH1D* tempmcproj = tempsimevrate->ProjectionX("b", i+1, i+1);
              tempmcproj->Scale(dataPOT/mcPOT);
              Int_t datamaxbin = tempdataproj->GetMaximumBin();
              double datamax = tempdataproj->GetBinContent(datamaxbin);
              Int_t mcmaxbin = tempmcproj->GetMaximumBin();
              double mcmax = tempmcproj->GetBinContent(mcmaxbin);
              double overallmax = (datamax > mcmax) ? datamax : mcmax;
              if (overallmax > maxVal) maxVal = overallmax;
          }
          std::cout<<"maxVal: " << maxVal <<std::endl;
          for (int i = 0; i < nbinsX; i++){
              std::cout<<"On pad " << i << "\n";
              PanelPad->cd(i+1);
              hist = tempevrate->ProjectionX("c", i+1, i+1);
              histmcproj = tempsimevrate->ProjectionX("d", i+1, i+1);
              DISHistProj = tempdis->ProjectionX("pxdis", i+1, i+1);
              RESHistProj = tempres->ProjectionX("pxres", i+1, i+1);
              QEHistProj = tempqe->ProjectionX("pxqe", i+1, i+1);
              HistProj2p2h = temp2p2h->ProjectionX("px2p2h", i+1, i+1);
              OtherHistProj = tempother->ProjectionX("pxother", i+1, i+1);
              TH1D* inverseProjection = tempevrate->ProjectionY("py", 1, 1); //used for auto-detecting subplot limits

              std::cout<<"Data/MC pot << " << dataPOT/mcPOT <<std::endl;
              Int_t datamaxbin = hist->GetMaximumBin();
              double datamax = hist->GetBinContent(datamaxbin);
              Int_t mcmaxbin = histmcproj->GetMaximumBin();
              double mcmax = histmcproj->GetBinContent(mcmaxbin);
              std::cout<<"data max " << datamax << std::endl;
              std::cout<<"mc max " << mcmax << std::endl;
              std::cout<<"mc max bin" << mcmaxbin << std::endl;
              double overallmax = (datamax > mcmax) ? datamax : mcmax;
              double scalefactor = (std::round((maxVal*5) / overallmax))/5; //round to nearest 0.2
              std::cout<<"scale factor " << scalefactor << std::endl;
              std::cout<<"max " << overallmax*scalefactor << std::endl;

              hist->Scale(scalefactor);
              std::cout<<"Here\n";
              histmcproj->Scale(scalefactor);

              DISHistProj->Scale(scalefactor);
              RESHistProj->Scale(scalefactor);
              QEHistProj->Scale(scalefactor);
              HistProj2p2h->Scale(scalefactor);
              OtherHistProj->Scale(scalefactor);
              gPad->Update();

              hist->SetLineColor(kBlack);
              hist->SetMarkerStyle(20);

              if (i != 12) histmcproj->GetXaxis()->ChangeLabel(1, -1, 0);
              else histmcproj->GetXaxis()->ChangeLabel(1, -1, -1, 31);
              histmcproj->GetXaxis()->ChangeLabel(6, -1, 0);
              histmcproj->GetXaxis()->SetRangeUser(0., 2.5);
              histmcproj->GetXaxis()->SetLabelFont(54);
              histmcproj->GetYaxis()->SetLabelFont(54);
              histmcproj->GetXaxis()->SetLabelSize(25);
              histmcproj->GetYaxis()->SetLabelSize(25);
              histmcproj->GetXaxis()->SetNdivisions(505);
              //histmcproj->GetYaxis()->LabelsOption("v");
              histmcproj->GetYaxis()->SetNdivisions(505);

              histmcproj->SetLineColor(kRed);
              histmcproj->SetLineWidth(3);
              DISHistProj->SetLineColor(6);
              DISHistProj->SetLineWidth(3);
              RESHistProj->SetLineColor(kOrange);
              RESHistProj->SetLineWidth(3);
              QEHistProj->SetLineColor(kBlue);
              QEHistProj->SetLineWidth(3);
              HistProj2p2h->SetLineColor(kGreen);
              HistProj2p2h->SetLineWidth(3);
              OtherHistProj->SetLineColor(kMagenta);
              OtherHistProj->SetLineWidth(3);
              TObjArray* mcArr = new TObjArray();
              mcArr->Add(histmcproj);
              mcArr->Add(DISHistProj);
              mcArr->Add(RESHistProj);
              mcArr->Add(QEHistProj);
              mcArr->Add(HistProj2p2h);
              mcArr->Add(OtherHistProj);
                
              plotter.mc_error_color=kRed;
              plotter.mc_error_style=3002;
              plotter.DrawDataMCVariations(hist, mcArr, 1.0, "N", true, true, true, false, false, true);
              if (i == (0)) //And check for number of panels/y bins
              {
                  canvas->cd();
                  TLegend *legend = new TLegend(0.55,0.1,0.9,0.25);
                  legend->SetBorderSize(0);
                  legend->SetNColumns(4);
                  legend->AddEntry(hist,"Data","lep");
                  legend->AddEntry(histmcproj,"MINERvA Tune v4.3.0","l");
                  legend->AddEntry(DISHistProj,"DIS","l");
                  legend->AddEntry(RESHistProj,"RES","l");
                  legend->AddEntry(QEHistProj,"QE","l");
                  legend->AddEntry(HistProj2p2h,"2p2h","l");
                  legend->AddEntry(OtherHistProj,"Other","l");
                  legend->Draw();
                  PanelPad->cd(i+1);
              }
                      
              //gROOT->SetSelectedPad(PanelPad->cd(i+1));

              std::stringstream lowerlim, upperlim, multip;
              lowerlim << std::fixed << std::setprecision(3) << inverseProjection->GetBinLowEdge(i+1);
              upperlim << std::fixed << std::setprecision(3) << inverseProjection->GetBinLowEdge(i+1)+inverseProjection->GetBinWidth(i+1);
              multip << std::fixed << std::setprecision(1) << scalefactor;
              std::string lowerlimit = lowerlim.str();
              std::string upperlimit = upperlim.str();
              std::string multiplier = std::string("#bf{#times ") +multip.str()+std::string("}");
              std::string textStr = std::string("#bf{")+lowerlimit+std::string(" < p_{Z} < ") + upperlimit + std::string("}");

              std::cout<<"X2: " << hist->GetXaxis()->GetXmax() << "\n";
              std::cout<<"Y2: " << hist->GetZaxis()->GetXmax() << "\n";
              TLatex limLabel;
              limLabel.SetTextSize(20);
              limLabel.SetTextAlign(33);
              limLabel.SetTextFont(43);
              limLabel.DrawLatexNDC(.95,.9,textStr.c_str());
              if (scalefactor!=1)
              {
                  TLatex multiplierlabel;
                  multiplierlabel.SetTextSize(20);
                  multiplierlabel.SetTextAlign(33);  //align at top
                  multiplierlabel.SetTextFont(43);
                  multiplierlabel.DrawLatexNDC(.95,.8,multiplier.c_str());
              }
              //delete t1;
          }

          canvas->cd();
          
          std::string title = "";

          std::stringstream POTstr;
          POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOT;
          std::string potstr = std::string("POT Used: ")+POTstr.str();
          plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.97, 0.025, 1, 42, 13, 0);

          plotter.WritePreliminary( 0.75, 0.075, 0.035, true);
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 35, 1, 43, 22, 0);
          //plotter.AddPlotLabel("Data/MC" , 0.08, 0.3, 0.03, 1, 42, 33, 90);

          plotter.AddPlotLabel("Events#times10^{3} /(Gev/c)^2", 0.020 ,0.51, 30, 1, 43, 22, 90);
          plotter.AddPlotLabel("Muon Transverse Momentum (Gev/c)", 0.51, 0.025, 30, 1, 43);
          canvas->Print((prefix+"_EvRateTransverse.png").c_str());
        }
        { //Cross-Section - Transverse
          TCanvas *canvas = new TCanvas("canvas", "canvas", 1500, 1000);
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
          plotter.axis_minimum=0.001;
          plotter.draw_normalized_to_bin_width=true;
          //Plotting
          TPad *PanelPad = new TPad("pad1","",0.02,0.05,1,0.95);
          PanelPad->Draw();
          PanelPad->Divide(4, 4, 0, 0);
          PlotUtils::MnvH2D* temp2p2h = simEventRate2P2H->Clone();
          PlotUtils::MnvH2D* tempdis = simEventRateDIS->Clone();
          PlotUtils::MnvH2D* tempres = simEventRateRES->Clone();
          PlotUtils::MnvH2D* tempqe = simEventRateQE->Clone();
          PlotUtils::MnvH2D* tempother = simEventRateOther->Clone();
          PlotUtils::MnvH2D* tempsimevrate = simEventRate->Clone();
          PlotUtils::MnvH2D* tempevrate = unfolded->Clone();
          //TH2D* tempevrate = tempsimevrate->Clone();

          temp2p2h->Scale(1e-3);
          tempdis->Scale(1e-3);
          tempres->Scale(1e-3);
          tempqe->Scale(1e-3);
          tempother->Scale(1e-3);
          tempsimevrate->Scale(1e-3);
          tempevrate->Scale(1e-3);

          int nbinsY = tempevrate->GetNbinsY();

          int nbinsX = tempevrate->GetNbinsX();



          PlotUtils::MnvH1D* hist;
          PlotUtils::MnvH1D* histmcproj;
          PlotUtils::MnvH1D* histmcbkgs;
          PlotUtils::MnvH1D* DISHistProj;
          PlotUtils::MnvH1D* RESHistProj;
          PlotUtils::MnvH1D* HistProj2p2h;
          PlotUtils::MnvH1D* QEHistProj;
          PlotUtils::MnvH1D* OtherHistProj;
          PlotUtils::MnvH1D* NCProj;

          //Autos-scaling
          double maxVal = 0 ;
          for (int i = 0; i < nbinsY; i++){
              TH1D* tempdataproj = tempevrate->ProjectionY("a", i+1, i+1);
              tempdataproj->Scale(1);
              TH1D* tempmcproj = tempsimevrate->ProjectionY("b", i+1, i+1);
              tempmcproj->Scale(dataPOT/mcPOT);
              Int_t datamaxbin = tempdataproj->GetMaximumBin();
              double datamax = tempdataproj->GetBinContent(datamaxbin);
              Int_t mcmaxbin = tempmcproj->GetMaximumBin();
              double mcmax = tempmcproj->GetBinContent(mcmaxbin);
              double overallmax = (datamax > mcmax) ? datamax : mcmax;
              if (overallmax > maxVal) maxVal = overallmax;
          }
          std::cout<<"maxVal: " << maxVal <<std::endl;
          for (int i = 0; i < nbinsX; i++){
              std::cout<<"On pad " << i << "\n";
              PanelPad->cd(i+1);
              hist = tempevrate->ProjectionY("c", i+1, i+1);
              histmcproj = tempsimevrate->ProjectionY("d", i+1, i+1);
              DISHistProj = tempdis->ProjectionY("pxdis", i+1, i+1);
              RESHistProj = tempres->ProjectionY("pxres", i+1, i+1);
              QEHistProj = tempqe->ProjectionY("pxqe", i+1, i+1);
              HistProj2p2h = temp2p2h->ProjectionY("px2p2h", i+1, i+1);
              OtherHistProj = tempother->ProjectionY("pxother", i+1, i+1);
              TH1D* inverseProjection = tempevrate->ProjectionX("py", 1, 1); //used for auto-detecting subplot limits

              std::cout<<"Data/MC pot << " << dataPOT/mcPOT <<std::endl;
              Int_t datamaxbin = hist->GetMaximumBin();
              double datamax = hist->GetBinContent(datamaxbin);
              Int_t mcmaxbin = histmcproj->GetMaximumBin();
              double mcmax = histmcproj->GetBinContent(mcmaxbin);
              std::cout<<"data max " << datamax << std::endl;
              std::cout<<"mc max " << mcmax << std::endl;
              std::cout<<"mc max bin" << mcmaxbin << std::endl;
              double overallmax = (datamax > mcmax) ? datamax : mcmax;
              double scalefactor = (std::round((maxVal*5) / overallmax))/5; //round to nearest 0.2
              std::cout<<"scale factor " << scalefactor << std::endl;
              std::cout<<"max " << overallmax*scalefactor << std::endl;

              hist->Scale(scalefactor);
              std::cout<<"Here\n";
              histmcproj->Scale(scalefactor);

              DISHistProj->Scale(scalefactor);
              RESHistProj->Scale(scalefactor);
              QEHistProj->Scale(scalefactor);
              HistProj2p2h->Scale(scalefactor);
              OtherHistProj->Scale(scalefactor);
              gPad->Update();

              hist->SetLineColor(kBlack);
              hist->SetMarkerStyle(20);

              if (i != 12) histmcproj->GetXaxis()->ChangeLabel(1, -1, 0);
              else histmcproj->GetXaxis()->ChangeLabel(1, -1, -1, 31);
              histmcproj->GetXaxis()->ChangeLabel(6, -1, 0);
              histmcproj->GetXaxis()->SetRangeUser(1.5, 15);
              histmcproj->GetXaxis()->SetLabelFont(54);
              histmcproj->GetYaxis()->SetLabelFont(54);
              histmcproj->GetXaxis()->SetLabelSize(25);
              histmcproj->GetYaxis()->SetLabelSize(25);
              histmcproj->GetXaxis()->SetNdivisions(502);
              //histmcproj->GetYaxis()->LabelsOption("v");
              histmcproj->GetYaxis()->SetNdivisions(505);

              histmcproj->SetLineColor(kRed);
              histmcproj->SetLineWidth(3);
              DISHistProj->SetLineColor(6);
              DISHistProj->SetLineWidth(3);
              RESHistProj->SetLineColor(kOrange);
              RESHistProj->SetLineWidth(3);
              QEHistProj->SetLineColor(kBlue);
              QEHistProj->SetLineWidth(3);
              HistProj2p2h->SetLineColor(kGreen);
              HistProj2p2h->SetLineWidth(3);
              OtherHistProj->SetLineColor(kMagenta);
              OtherHistProj->SetLineWidth(3);
              TObjArray* mcArr = new TObjArray();
              mcArr->Add(histmcproj);
              mcArr->Add(DISHistProj);
              mcArr->Add(RESHistProj);
              mcArr->Add(QEHistProj);
              mcArr->Add(HistProj2p2h);
              mcArr->Add(OtherHistProj);
                
              plotter.mc_error_color=kRed;
              plotter.mc_error_style=3002;
              plotter.DrawDataMCVariations(hist, mcArr, 1.0, "N", true, true, true, false, false, true);
              if (i == (0)) //And check for number of panels/y bins
              {
                  canvas->cd();
                  TLegend *legend = new TLegend(0.55,0.1,0.9,0.25);
                  legend->SetBorderSize(0);
                  legend->SetNColumns(4);
                  legend->AddEntry(hist,"Data","lep");
                  legend->AddEntry(histmcproj,"MINERvA Tune v4.3.0","l");
                  legend->AddEntry(DISHistProj,"DIS","l");
                  legend->AddEntry(RESHistProj,"RES","l");
                  legend->AddEntry(QEHistProj,"QE","l");
                  legend->AddEntry(HistProj2p2h,"2p2h","l");
                  legend->AddEntry(OtherHistProj,"Other","l");
                  legend->Draw();
                  PanelPad->cd(i+1);
              }
                      
              //gROOT->SetSelectedPad(PanelPad->cd(i+1));

              std::stringstream lowerlim, upperlim, multip;
              lowerlim << std::fixed << std::setprecision(3) << inverseProjection->GetBinLowEdge(i+1);
              upperlim << std::fixed << std::setprecision(3) << inverseProjection->GetBinLowEdge(i+1)+inverseProjection->GetBinWidth(i+1);
              multip << std::fixed << std::setprecision(1) << scalefactor;
              std::string lowerlimit = lowerlim.str();
              std::string upperlimit = upperlim.str();
              std::string multiplier = std::string("#bf{#times ") +multip.str()+std::string("}");
              std::string textStr = std::string("#bf{")+lowerlimit+std::string(" < p_{t} < ") + upperlimit + std::string("}");

              std::cout<<"X2: " << hist->GetXaxis()->GetXmax() << "\n";
              std::cout<<"Y2: " << hist->GetZaxis()->GetXmax() << "\n";
              TLatex limLabel;
              limLabel.SetTextSize(20);
              limLabel.SetTextAlign(33);
              limLabel.SetTextFont(43);
              limLabel.DrawLatexNDC(.95,.9,textStr.c_str());
              if (scalefactor!=1)
              {
                  TLatex multiplierlabel;
                  multiplierlabel.SetTextSize(20);
                  multiplierlabel.SetTextAlign(33);  //align at top
                  multiplierlabel.SetTextFont(43);
                  multiplierlabel.DrawLatexNDC(.95,.8,multiplier.c_str());
              }
              //delete t1;
          }

          canvas->cd();
          
          std::string title = "";

          std::stringstream POTstr;
          POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOT;
          std::string potstr = std::string("POT Used: ")+POTstr.str();
          plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.97, 0.025, 1, 42, 13, 0);

          plotter.WritePreliminary( 0.75, 0.075, 0.035, true);
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 35, 1, 43, 22, 0);
          //plotter.AddPlotLabel("Data/MC" , 0.08, 0.3, 0.03, 1, 42, 33, 90);

          plotter.AddPlotLabel("Events#times10^{3} /(Gev/c)^2", 0.020 ,0.51, 30, 1, 43, 22, 90);
          plotter.AddPlotLabel("Muon Longitudinal Momentum (Gev/c)", 0.51, 0.025, 30, 1, 43);
          canvas->Print((prefix+"_EvRateLongitudinal.png").c_str());
        }

        
      }
      // return 0;
    }
    catch (const std::runtime_error &e)
    {
      std::cerr << "Failed to extract a cross section for prefix " << prefix <<  " : " << e.what() << "\n";
      return 4;
      // break;
    }
  }
  return 0;
}
