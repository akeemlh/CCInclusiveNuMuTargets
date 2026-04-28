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

bool verbose = false;

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

int main(const int argc, const char **argv)
{
#ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable(); // Needed to look up dictionaries for PlotUtils classes like MnvH1D
#endif

  TH1::AddDirectory(kFALSE); // Needed so that MnvH1D gets to clean up its own MnvLatErrorBands (which are TH1Ds).

  if (argc != 4 && argc != 5)
  {
    std::cerr << "Expected 4 or 5 arguments, but I got " << argc - 1 << ".\n"  << HELP <<  std::endl;  
    return 1;
  }

  const int nIterations = std::stoi(argv[1]);
  std::string indir = std::string(argv[2]);
  int pdg = std::stoi(argv[3]);
  std::string variable = std::string(argv[4]);

  bool doDaisy = true;

  std::vector<std::string> dirs = util::findContainingDirectories(indir, "Tracker", true, true);

  if (dirs.size()==0) dirs = {"./"};
  
  //std::vector<std::string> crossSectionPrefixes = {"pTmu", "pZmu", "BjorkenX", "BjorkenY", "Erecoil", "Emu" , "beamAngle", "segment"};
  std::vector<std::string> crossSectionPrefixes = {variable};
  //std::vector<std::string> crossSectionPrefixes = {"pTmu", "pZmu", "BjorkenX", "Erecoil", "Emu"};
  
  PlotUtils::MinervaUniverse::SetNuEConstraint(true);
  PlotUtils::MinervaUniverse::SetNFluxUniverses(100);

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
      double max_energy = 120;
      
      PlotUtils::MnvH1D *fluxIntReweighted = new PlotUtils::MnvH1D();
      std::shared_ptr<PlotUtils::MnvH1D> DaisyEffNum[12], DaisyFolded[12];
      std::shared_ptr<PlotUtils::MnvH2D> DaisyMigration[12];
      std::vector<PlotUtils::MnvH1D*> DaisyBackgrounds[12];
      std::shared_ptr<PlotUtils::MnvH1D> DaisyEffDenom[12], DaisyEffDenom2P2H[12], DaisyEffDenomDIS[12], DaisyEffDenomRES[12], DaisyEffDenomQE[12], DaisyEffDenomOther[12];

      for (int c = 0; c<12; c++)
      {
        DaisyEffNum[c] = std::shared_ptr<PlotUtils::MnvH1D>(new PlotUtils::MnvH1D());
        DaisyFolded[c] = std::shared_ptr<PlotUtils::MnvH1D>(new PlotUtils::MnvH1D());
        DaisyMigration[c] = std::shared_ptr<PlotUtils::MnvH2D>(new PlotUtils::MnvH2D());
        DaisyEffDenom[c] = std::shared_ptr<PlotUtils::MnvH1D>(new PlotUtils::MnvH1D());
        DaisyEffDenom2P2H[c] = std::shared_ptr<PlotUtils::MnvH1D>(new PlotUtils::MnvH1D());
        DaisyEffDenomDIS[c] = std::shared_ptr<PlotUtils::MnvH1D>(new PlotUtils::MnvH1D());
        DaisyEffDenomRES[c] = std::shared_ptr<PlotUtils::MnvH1D>(new PlotUtils::MnvH1D());
        DaisyEffDenomQE[c] = std::shared_ptr<PlotUtils::MnvH1D>(new PlotUtils::MnvH1D());
        DaisyEffDenomOther[c] = std::shared_ptr<PlotUtils::MnvH1D>(new PlotUtils::MnvH1D());
      }

      PlotUtils::MnvH1D* flux = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* folded = new PlotUtils::MnvH1D();
      PlotUtils::MnvH2D* migration = new PlotUtils::MnvH2D();
      PlotUtils::MnvH1D* effNum = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* effDenom = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* effDenom2P2H = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* effDenomDIS = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* effDenomRES = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* effDenomQE = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* effDenomOther = new PlotUtils::MnvH1D();

      PlotUtils::MnvH1D* BackgroundWrongSign = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* BackgroundNC = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* BackgroundOther = new PlotUtils::MnvH1D();

      PlotUtils::MnvH1D* SelectedSignalReco = new PlotUtils::MnvH1D();

      PlotUtils::MnvH1D* MCData = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* MCData2p2h = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* MCDataDIS = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* MCDataRES = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* MCDataQE = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* MCDataOther = new PlotUtils::MnvH1D();

      std::vector<PlotUtils::MnvH1D *> backgrounds;

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
          bool petalFilesExist = std::filesystem::exists(datapathdaisy) && std::filesystem::exists(mcpathdaisy);
          doDaisy = doDaisy && petalFilesExist;
        }
        //doDaisy=false;
        if (doDaisy) std::cout<<"Will do daisy petal reweight in this analysis\n";
        else std::cout<<"Will not do daisy petal reweight in this analysis\n";
        for (int petal=0; petal<12 && doDaisy; petal++){
          std::string datapathdaisy = dirs[c] + "/runEventLoopTrackerData_petal_"+std::to_string(petal)+".root";
          std::string mcpathdaisy = dirs[c] + "/runEventLoopTrackerMC_petal_"+std::to_string(petal)+".root";
          if (verbose) std::cout<< "datapathdaisy: "<<datapathdaisy <<std::endl;
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
            break;
          }

          //std::cout<< "Here2"<<std::endl;

          auto tmpEffNum = util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("efficiency_numerator")), prefix);
          auto tmpEffDenom = util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("efficiency_denominator")), prefix);
          auto tmpMigration = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, (std::string("migration")), prefix);
          auto tmpFolded = util::GetIngredient<PlotUtils::MnvH1D>(*dataDaisyFile, (std::string("data")), prefix);
          auto tmpEffDenom2P2H = util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("efficiency_denominator_intChannels_2p2h")), prefix);
          auto tmpEffDenomDIS = util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("efficiency_denominator_intChannels_DIS")), prefix);
          auto tmpEffDenomRES = util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("efficiency_denominator_intChannels_RES")), prefix);
          auto tmpEffDenomQE = util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("efficiency_denominator_intChannels_QE")), prefix);
          auto tmpEffDenomOther = util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("efficiency_denominator_intChannels_Other")), prefix);

          auto tmpBkgWS = util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("background_Wrong_Sign")), prefix);
          auto tmpBkgNC = util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("background_NC")), prefix);
          auto tmpBkgOther = util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("background_Other")), prefix);

          auto tmpSelectedSignalReco = util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("selected_signal_reco")), prefix);
          
          auto tmpMCData = util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("data")), prefix);
          auto tmpMCData2p2h = util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("intChannels_2p2h")), prefix);
          auto tmpMCDataDIS = util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("intChannels_DIS")), prefix);
          auto tmpMCDataRES = util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("intChannels_RES")), prefix);
          auto tmpMCDataQE = util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("intChannels_QE")), prefix);
          auto tmpMCDataOther = util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("intChannels_Other")), prefix);

          if (petal>=0)
          {
            //std::cout<< "Here2.1"<<std::endl;
            if (verbose) std::cout<< "petal:"<<petal<<std::endl;
            util::AddHist(*DaisyEffNum[petal], tmpEffNum);
            //std::cout<< "Here2.11"<<std::endl;

            util::AddHist(*DaisyEffDenom[petal], tmpEffDenom);
            util::AddHist(*DaisyMigration[petal], tmpMigration);
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

          if (petal >=0){
            util::AddHist(*folded, tmpFolded);
            util::AddHist(*migration,tmpMigration);
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
            
            //bool useMuCorr = true;
            //if (nupdg<14) useMuCorr == false;
            //PlotUtils::FluxReweighter frw_temp = PlotUtils::flux_reweighter(playlistUsed, nupdg, use_nue_constraint, n_flux_universes);
            //auto tempIntFlux = frw_temp.GetIntegratedFluxReweighted(nupdg, effDenom, min_energy, max_energy, useMuCorr)->Clone();

            PlotUtils::FluxReweighter frw_temp = PlotUtils::FluxReweighter( nupdg, use_nue_constraint, playlistUsed, PlotUtils::FluxReweighter::gen2thin, PlotUtils::FluxReweighter::g4numiv6, n_flux_universes );
            auto tempIntFlux = frw_temp.GetIntegratedFluxReweighted(nupdg, effDenom, min_energy, max_energy, true)->Clone();
            tempIntFlux->Scale(tempDataPOT);
            util::AddHist(*fluxIntReweighted,tempIntFlux);
            delete tempIntFlux;
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
      Plot(*toSubtract, "BackgroundSum", prefix);
      std::cout<< "Here5.2"<<std::endl;

      auto bkgSubtracted = std::accumulate(backgrounds.begin(), backgrounds.end(), folded->Clone(),
                                            [mcPOT, dataPOT](auto sum, const auto hist)
                                            {
                                              std::cout << "Subtracting " << hist->GetName() << " scaled by " << -dataPOT / mcPOT << " from " << sum->GetName() << "\n";
                                              sum->Add(hist, -dataPOT / mcPOT);
                                              return sum;
                                            });
      Plot(*bkgSubtracted, "backgroundSubtracted", prefix);
      //std::cout<< "Here6"<<std::endl;

      std::string nutypestr;
      if (pdg>0) nutypestr = "Neutrino";
      else nutypestr = "Antineutrino";

      auto outFile = TFile::Open((nutypestr+"_"+prefix + "_crossSection.root").c_str(), "RECREATE");
      if (!outFile)
      {
        std::cerr << "Could not create a file called " << prefix + "_crossSection.root" << ".  Does it already exist?\n";
        return 5;
      }

      bkgSubtracted->Write("backgroundSubtracted");
      // d'Aogstini unfolding
      migration->Write("migration");
      //return 0;
      auto unfolded = UnfoldHist(bkgSubtracted, migration, nIterations);
      if (!unfolded)
        throw std::runtime_error(std::string("Failed to unfold ") + folded->GetName() + " using " + migration->GetName());
      //Plot(*unfolded, "unfolded", prefix);
      unfolded->Clone()->Write("unfolded"); // TODO: Seg fault first appears when I uncomment this line
      std::cout << "Survived writing the unfolded histogram.\n"
                << std::endl; // This is evidence that the problem is on the final file Write() and not unfolded->Clone()->Write().
      effNum->Divide(effNum, effDenom); // Only the 2 parameter version of MnvH1D::Divide()
                                        // handles systematics correctly.
      Plot(*effNum, "efficiency", prefix);
      effNum->Write(("efficiency"));
      unfolded->Divide(unfolded, effNum);
      Plot(*unfolded, "efficiencyCorrected", prefix);


      PlotUtils::TargetUtils targetInfo;
      double nnucleons = targetInfo.GetTrackerNNucleons(PlotUtils::TargetProp::Tracker::Face, PlotUtils::TargetProp::Tracker::Back, true);
      double nnucleonsData = targetInfo.GetTrackerNNucleons(PlotUtils::TargetProp::Tracker::Face, PlotUtils::TargetProp::Tracker::Back, false);

      PlotUtils::MnvH1D *flux2;
      PlotUtils::MnvH1D *fluxIntegral;
      PlotUtils::MnvH1D *fluxRebinned;

      outFile->cd();
      auto crossSection = normalize(unfolded, fluxIntReweighted, nnucleonsData, dataPOT);
      Plot(*crossSection, "crossSection", prefix);
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
      auto crossSection2 = normalize(simEventRate, fluxIntReweighted, nnucleons, mcPOT);
      Plot(*crossSection2, "simulatedCrossSection", prefix);
      crossSection2->Write("simulatedCrossSection");
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
      delete fluxIntReweighted;
      outFile->Close();

      //Daisy reweight-----------------------------------------
      //map of petal distributions

      // ---------------------------------------------------------------------
      // Flux reweighter information, get reweighted daisy sum according to a material
      // ---------------------------------------------------------------------
      if (doDaisy)
      {
        auto outFileDaisy = TFile::Open((nutypestr+"_"+prefix + "_Daisy_crossSection.root").c_str(), "RECREATE");
        if (!outFileDaisy)
        {
          std::cerr << "Could not create a file called " << prefix + "_crossSection.root" << ".  Does it already exist?\n";
          return 5;
        }

        std::map<int, PlotUtils::MnvH1D*> daisy_petal_hists;
        std::map<int, PlotUtils::MnvH1D*> daisy_petal_hists_eff_denom;
        std::map<int, PlotUtils::MnvH1D*> daisy_petal_hists_eff_denom_2p2h;
        std::map<int, PlotUtils::MnvH1D*> daisy_petal_hists_eff_denom_dis;
        std::map<int, PlotUtils::MnvH1D*> daisy_petal_hists_eff_denom_res;
        std::map<int, PlotUtils::MnvH1D*> daisy_petal_hists_eff_denom_qe;
        std::map<int, PlotUtils::MnvH1D*> daisy_petal_hists_eff_denom_other;



        for (int petal=0; petal<12; petal++){
          std::vector<PlotUtils::MnvH1D*> DaisyBackgroundsTemp = DaisyBackgrounds[petal];
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
          auto unfoldedDaisy = UnfoldHist(bkgSubtractedDaisy, DaisyMigration[petal].get(), nIterations);
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
          delete toSubtractDaisy;
          delete bkgSubtractedDaisy;

        }
        for (int p = 0; p < 12; p++) {
          for (auto hist : DaisyBackgrounds[p]) delete hist;
          DaisyBackgrounds[p].clear();
        }
        delete toSubtract;
        delete bkgSubtracted; 
        //std::unique_ptr<PlotUtils::MnvH1D> fluxIntegralCtmp = std::unique_ptr<PlotUtils::MnvH1D>(new PlotUtils::MnvH1D());
        PlotUtils::MnvH1D *fluxIntegralC = new PlotUtils::MnvH1D();
        PlotUtils::MnvH1D *fluxIntegralFe = new PlotUtils::MnvH1D();
        PlotUtils::MnvH1D *fluxIntegralPb = new PlotUtils::MnvH1D(); 
        PlotUtils::MnvH1D *DaisyCorrectedC = new PlotUtils::MnvH1D();
        PlotUtils::MnvH1D *DaisyCorrectedFe = new PlotUtils::MnvH1D();
        PlotUtils::MnvH1D *DaisyCorrectedPb = new PlotUtils::MnvH1D(); 
        PlotUtils::MnvH1D *DaisyCorrectedSimC = new PlotUtils::MnvH1D();
        PlotUtils::MnvH1D *DaisyCorrectedSimFe = new PlotUtils::MnvH1D();
        PlotUtils::MnvH1D *DaisyCorrectedSimPb = new PlotUtils::MnvH1D(); 
        PlotUtils::MnvH1D *DaisyCorrectedSim2p2hC = new PlotUtils::MnvH1D();
        PlotUtils::MnvH1D *DaisyCorrectedSim2p2hFe = new PlotUtils::MnvH1D();
        PlotUtils::MnvH1D *DaisyCorrectedSim2p2hPb = new PlotUtils::MnvH1D(); 
        PlotUtils::MnvH1D *DaisyCorrectedSimDISC = new PlotUtils::MnvH1D();
        PlotUtils::MnvH1D *DaisyCorrectedSimDISFe = new PlotUtils::MnvH1D();
        PlotUtils::MnvH1D *DaisyCorrectedSimDISPb = new PlotUtils::MnvH1D(); 
        PlotUtils::MnvH1D *DaisyCorrectedSimRESC = new PlotUtils::MnvH1D();
        PlotUtils::MnvH1D *DaisyCorrectedSimRESFe = new PlotUtils::MnvH1D();
        PlotUtils::MnvH1D *DaisyCorrectedSimRESPb = new PlotUtils::MnvH1D(); 
        PlotUtils::MnvH1D *DaisyCorrectedSimQEC = new PlotUtils::MnvH1D();
        PlotUtils::MnvH1D *DaisyCorrectedSimQEFe = new PlotUtils::MnvH1D();
        PlotUtils::MnvH1D *DaisyCorrectedSimQEPb = new PlotUtils::MnvH1D();
        PlotUtils::MnvH1D *DaisyCorrectedSimOtherC = new PlotUtils::MnvH1D();
        PlotUtils::MnvH1D *DaisyCorrectedSimOtherFe = new PlotUtils::MnvH1D();
        PlotUtils::MnvH1D *DaisyCorrectedSimOtherPb = new PlotUtils::MnvH1D(); 
        for (auto pp : playlistPOTpair)
        {
          //PlotUtils::FluxReweighter frw = PlotUtils::flux_reweighter(pp.first, pdg, use_nue_constraint, n_flux_universes);
          PlotUtils::FluxReweighter frw( pdg, use_nue_constraint, pp.first, PlotUtils::FluxReweighter::gen2thin, PlotUtils::FluxReweighter::g4numiv6, n_flux_universes );

          util::AddHist(*DaisyCorrectedC, frw.GetReweightedDaisySum(pdg, "carbon", daisy_petal_hists, project_dir ), pp.second);
          util::AddHist(*fluxIntegralC,frw.GetIntegratedTargetFlux(pdg, "carbon", DaisyCorrectedC, min_energy, max_energy, project_dir), pp.second);
          util::AddHist(*DaisyCorrectedFe, frw.GetReweightedDaisySum(pdg, "iron", daisy_petal_hists, project_dir ), pp.second);
          util::AddHist(*fluxIntegralFe,frw.GetIntegratedTargetFlux(pdg, "iron", DaisyCorrectedFe, min_energy, max_energy, project_dir), pp.second);
          util::AddHist(*DaisyCorrectedPb, frw.GetReweightedDaisySum(pdg, "lead", daisy_petal_hists, project_dir ), pp.second);
          util::AddHist(*fluxIntegralPb,frw.GetIntegratedTargetFlux(pdg, "lead", DaisyCorrectedPb, min_energy, max_energy, project_dir), pp.second);

          util::AddHist(*DaisyCorrectedSimC, frw.GetReweightedDaisySum(pdg, "carbon", daisy_petal_hists_eff_denom, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimFe, frw.GetReweightedDaisySum(pdg, "iron", daisy_petal_hists_eff_denom, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimPb, frw.GetReweightedDaisySum(pdg, "lead", daisy_petal_hists_eff_denom, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSim2p2hC, frw.GetReweightedDaisySum(pdg, "carbon", daisy_petal_hists_eff_denom_2p2h, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSim2p2hFe, frw.GetReweightedDaisySum(pdg, "iron", daisy_petal_hists_eff_denom_2p2h, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSim2p2hPb, frw.GetReweightedDaisySum(pdg, "lead", daisy_petal_hists_eff_denom_2p2h, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimDISC, frw.GetReweightedDaisySum(pdg, "carbon", daisy_petal_hists_eff_denom_dis, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimDISFe, frw.GetReweightedDaisySum(pdg, "iron", daisy_petal_hists_eff_denom_dis, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimDISPb, frw.GetReweightedDaisySum(pdg, "lead", daisy_petal_hists_eff_denom_dis, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimRESC, frw.GetReweightedDaisySum(pdg, "carbon", daisy_petal_hists_eff_denom_res, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimRESFe, frw.GetReweightedDaisySum(pdg, "iron", daisy_petal_hists_eff_denom_res, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimRESPb, frw.GetReweightedDaisySum(pdg, "lead", daisy_petal_hists_eff_denom_res, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimQEC, frw.GetReweightedDaisySum(pdg, "carbon", daisy_petal_hists_eff_denom_qe, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimQEFe, frw.GetReweightedDaisySum(pdg, "iron", daisy_petal_hists_eff_denom_qe, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimQEPb, frw.GetReweightedDaisySum(pdg, "lead", daisy_petal_hists_eff_denom_qe, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimOtherC, frw.GetReweightedDaisySum(pdg, "carbon", daisy_petal_hists_eff_denom_other, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimOtherFe, frw.GetReweightedDaisySum(pdg, "iron", daisy_petal_hists_eff_denom_other, project_dir ), pp.second);
          util::AddHist(*DaisyCorrectedSimOtherPb, frw.GetReweightedDaisySum(pdg, "lead", daisy_petal_hists_eff_denom_other, project_dir ), pp.second);
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

        PlotUtils::MnvH1D *crossSectionC, *crossSectionFe, *crossSectionPb;
        
        PlotUtils::MnvH1D *crossSectionSimC, *crossSectionSimFe, *crossSectionSimPb;
        PlotUtils::MnvH1D *crossSectionSimC2p2h, *crossSectionSimFe2p2h, *crossSectionSimPb2p2h;
        PlotUtils::MnvH1D *crossSectionSimCDIS, *crossSectionSimFeDIS, *crossSectionSimPbDIS;
        PlotUtils::MnvH1D *crossSectionSimCRes, *crossSectionSimFeRes, *crossSectionSimPbRes;
        PlotUtils::MnvH1D *crossSectionSimCQE, *crossSectionSimFeQE, *crossSectionSimPbQE;
        PlotUtils::MnvH1D *crossSectionSimCOther, *crossSectionSimFeOther, *crossSectionSimPbOther;



        int mcColors[12] = {   TColor::GetColor("#332288"), 
                          TColor::GetColor("#117733"),
                          TColor::GetColor("#44AA99"), 
                          TColor::GetColor("#88CCEE"), 
                          TColor::GetColor("#DDCC77"), 
                          TColor::GetColor("#E2AC34"), 
                          TColor::GetColor("#CC6677"), 
                          TColor::GetColor("#AA4499"), 
                          TColor::GetColor("#882255"), 
                          TColor::GetColor("#9E92D8"), 
                          17,
                          12};

        double normPOT = 1.06e+21;
        { 
          TCanvas can("Temp");
          PlotUtils::MnvPlotter plotter(PlotUtils::kDefaultStyle);

          TPad pad("pad1","pad1",0,0,1.0,1.0);
          pad.Draw();
          TLegend leg(0.6,0.5,0.9,0.9);
          leg.SetBorderSize(0);   // no border
          leg.SetFillColorAlpha(0, 0);    // transparent background
          std::vector<PlotUtils::MnvH1D*> dataHists, mcHists; //To keep track of histogram pointers from Clone() below for deletion
          for (int p = 0; p<12; p++)
          {
            auto tmpData = daisy_petal_hists[p]->Clone((std::string("Petal ") + std::to_string(p)).c_str());
            auto tmpMC = daisy_petal_hists_eff_denom[p]->Clone((std::string("Petal ") + std::to_string(p)).c_str());
            dataHists.push_back(tmpData);
            mcHists.push_back(tmpMC);
            tmpData->Scale(1e-3);
            tmpMC->Scale(1e-3);
            tmpMC->Scale(dataPOT/mcPOT);

            tmpData->SetMarkerStyle(20);
            tmpData->SetMarkerSize(1);
            tmpData->SetMarkerColor(mcColors[p]);
            tmpData->SetLineWidth(1);
            tmpData->SetLineStyle(1);
            tmpData->SetLineColor(mcColors[p]);

            tmpMC->SetLineColor(mcColors[p]);
            tmpMC->SetLineWidth(3);
            tmpMC->SetLineStyle(1);

            std::string tmp, xaxislabel, yaxislabel;
            if (prefix=="pTmu")
            {
              tmp = "Tracker Muon p_{T}";
              xaxislabel = "p_{T, #mu} [GeV/c]";
              yaxislabel = "Events #times 10^{3} /(GeV/c)";

              tmpMC->GetXaxis()->SetRange(1, 13);
            }
            else if (prefix=="pZmu")
            {
              tmp = "Tracker  Muon p_{Z}";
              xaxislabel = "p_{Z, #mu} [GeV/c]";
              tmpMC->GetXaxis()->SetRange(1, 13);
              tmpData->GetXaxis()->SetRange(1, 13);
            }
            else if (prefix=="BjorkenX")
            {
              tmp = "Tracker  Bjorken X";
              xaxislabel = "X";
              yaxislabel = "Events #times 10^{3} per unit X";
              tmpMC->GetXaxis()->SetRange(1, 9);
              can.SetLogx();
            }
            else if (prefix=="Erecoil")
            {
              tmp = "Tracker E_{recoil}";
              xaxislabel = "E_{recoil} [GeV]";
              yaxislabel = "Events #times 10^{3} /(GeV/c)";

            } 
            else if (prefix=="Emu")
            {
              tmp = "Tracker  E_{#mu}";
              xaxislabel = "E_{#mu} [GeV]";
              yaxislabel = "Events #times 10^{3} /(GeV/c)";
              tmpMC->GetXaxis()->SetRange(1, 11);
            }
            tmpMC->GetYaxis()->SetTitle(yaxislabel.c_str());
            tmpMC->GetXaxis()->SetTitle(xaxislabel.c_str());
            tmpMC->GetXaxis()->SetTitleSize(0.04);
            tmpMC->SetTitle(tmp.c_str());
            tmpMC->Scale(1, "Width");
            tmpData->Scale(1, "Width");
            tmpMC->SetMaximum(tmpMC->GetMaximum()*1.4);
            tmpMC->Draw("HIST SAME");
            tmpData->Draw("P SAME");

            std::stringstream POTstr;
            POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOT;
            std::string potstr = std::string("POT Normalised: ")+POTstr.str();

            TLatex latex;
            latex.SetNDC();               // use normalized coordinates (0–1 across pad)
            latex.SetTextSize(0.02);      // text size
            latex.SetTextAlign(31);       // align top-left
            //latex.SetTextColor(12);
            latex.SetTextFont(62);
            latex.DrawLatex(0.93, 0.94, potstr.c_str());
            TLatex latex2;
            latex2.SetNDC();               // use normalized coordinates (0–1 across pad)
            latex2.SetTextSize(0.04);      // text size
            latex2.SetTextAlign(22);       // align top-left
            latex2.SetTextFont(62);       // align top-left
            //latex2.SetTextColor(12);
            latex2.DrawLatex(0.5, 0.95, tmp.c_str());

            leg.AddEntry(tmpMC, (std::string("Petal ") + std::to_string(p)).c_str(), "l");
          }
          leg.Draw();
          plotter.WritePreliminary( 0.4, 0.87, 0.035, false);
          can.Print((std::string("Tracker_")+prefix+"_EvRate.png").c_str());
          for (auto h : dataHists) delete h;
          for (auto h : mcHists) delete h;
        }

        for (auto p : daisy_petal_hists) delete p.second;

        { //Plotting event rates
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter(PlotUtils::kDefaultStyle);
          plotter.axis_label_font=43;
          plotter.axis_label_size = 25;
          plotter.headroom = 1.2;
          /* plotter.axis_title_font_y=43;
          plotter.axis_title_size_y = 35;
          plotter.axis_title_font_x=43;
          plotter.axis_title_size_x = 35;
          */
          plotter.ApplyAxisStyle();
          //plotter.axis_maximum = 0.4;
          auto temp2p2h = MCData2p2h->Clone();
          auto tempdis = MCDataDIS->Clone();
          auto tempres = MCDataRES->Clone();
          auto tempqe = MCDataQE->Clone();
          auto tempChannelOther = MCDataOther->Clone();

          auto tempWS = BackgroundWrongSign->Clone();
          auto tempNC = BackgroundNC->Clone();
          auto tempBackgroundOther = BackgroundOther->Clone();
          auto tempsignal = SelectedSignalReco->Clone();

          auto tempevrate = folded->Clone();
          auto tempsim = MCData->Clone();

          tempWS->SetTitle("Wrong Sign");
          tempNC->SetTitle("Neutral Current");
          tempBackgroundOther->SetTitle("Other");
          tempsignal->SetTitle("Signal");


          temp2p2h->Scale(1e-3);
          tempdis->Scale(1e-3);
          tempres->Scale(1e-3);
          tempqe->Scale(1e-3);
          tempChannelOther->Scale(1e-3);

          tempWS->Scale(1e-3);
          tempNC->Scale(1e-3);
          tempBackgroundOther->Scale(1e-3);
          tempsim->Scale(1e-3);
          tempsignal->Scale(1e-3);

          tempevrate->Scale(1e-3);

          //double normPOT = 1.06e+21;

          temp2p2h->Scale(dataPOT/mcPOT);
          tempdis->Scale(dataPOT/mcPOT);
          tempres->Scale(dataPOT/mcPOT);
          tempqe->Scale(dataPOT/mcPOT);
          tempChannelOther->Scale(dataPOT/mcPOT);

          tempWS->Scale(dataPOT/mcPOT);
          tempNC->Scale(dataPOT/mcPOT);
          tempBackgroundOther->Scale(dataPOT/mcPOT);
          tempsignal->Scale(dataPOT/mcPOT);

          tempsim->Scale(dataPOT/mcPOT);
          //std::cout<<"dataScale = " << normPOT << "/" << dataPOT << " = " << normPOT/dataPOT << std::endl;
          //std::cout<<"mcScale = " << normPOT << "/" << mcPOT << " = " << normPOT/mcPOT << std::endl;

          //tempevrate->Scale(normPOT/dataPOT);

          TObjArray channelArr;
          channelArr.Add(temp2p2h);
          channelArr.Add(tempdis);
          channelArr.Add(tempres);
          channelArr.Add(tempqe);
          channelArr.Add(tempChannelOther);

          TObjArray bkgArr;
          bkgArr.Add(tempWS);
          bkgArr.Add(tempNC);
          bkgArr.Add(tempBackgroundOther);
          bkgArr.Add(tempsignal);

          std::string tmp; 
          std::string y_title="Events#times10^{3}";
          std::cout<<"Here2.n";
          double xmax = -1111; //Default value
          std::string xaxislabel = "";
          if (prefix=="pTmu")
          {
            tmp = "Muon p_{T}";
            xaxislabel = "p_{T, #mu} [GeV/c]";
            y_title += " /(GeV/c)";
            xmax = 13;
            tempsim->GetXaxis()->SetRange(1);
          }
          else if (prefix=="pZmu")
          {
            tmp = "Muon p_{Z}";
            xaxislabel = "p_{Z, #mu} [GeV/c]";
            y_title += " /(GeV/c)";
            xmax = 13;
            tempsim->GetXaxis()->SetRange(1);
          }
          else if (prefix=="BjorkenX")
          {
            tmp = "Bjorken X";
            xaxislabel = "Bjorken X";
            y_title += " per unit X";
            xmax = 7;
            tempsim->GetXaxis()->SetRange(1);
            can.SetLogx();
          }
          else if (prefix=="Erecoil")
          {
            tmp = "E_{recoil}";
            y_title += " /GeV";
            xaxislabel = "E_{recoil} [GeV]";
          } 
          else if (prefix=="Emu")
          {
            tmp = "E_{#mu}";
            y_title += " /GeV";
            xaxislabel = "E_{#mu} [GeV]";
            xmax = 11;
            tempsim->GetXaxis()->SetRange(1);
          } 

          std::string title = tmp + " Plastic Tracker";
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;

          plotter.DrawErrorSummary(tempsim);
          can.Print((std::string("Tracker_")+prefix+"_EvRateErrors.png").c_str());
          can.Clear();


          // Adjust pad heights manually (top bigger, bottom smaller)
          TPad pad1("pad1","pad1",0,0.35,1,1.0); // top pad
          TPad pad2("pad2","pad2",0,0.05,1,0.35); // bottom pad
          pad1.SetBottomMargin(0); // remove x-axis label space for top
          pad2.SetTopMargin(0);
          pad2.SetBottomMargin(0.25);
          pad1.Draw();
          pad2.Draw();

          if (prefix=="BjorkenX")
          {
            pad1.SetLogx();
            pad2.SetLogx();
          }

          double defaultmin = plotter.axis_minimum;
          plotter.legend_offset_x = 0.1325;
          plotter.legend_border_size = 0;
          plotter.legend_text_size = 0.03;

          plotter.mc_line_width=0;
          std::stringstream POTstr;
          POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOT;
          std::string potstr = std::string("POT Normalised: ")+POTstr.str();
          plotter.WritePreliminary( 0.5, 0.15, 0.035, true);
          plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.93, 0.025, 1, 42, 13, 0);
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 25, 1, 53, 22, 0);
          plotter.AddPlotLabel(y_title.c_str() , 0.08, 0.9, 0.03, 1, 42, 33, 90);
          plotter.AddPlotLabel("Data/MC" , 0.08, 0.3, 0.03, 1, 42, 33, 90);
          plotter.AddPlotLabel(xaxislabel.c_str() , 0.9, 0.08, 0.03, 1, 42, 33, 0);

          //Interaction Type
          pad1.cd(); 
          plotter.axis_minimum=0.001;

          int channelColArr[5] = {TColor::GetColor("#d2c271"), TColor::GetColor("#CA5454"), TColor::GetColor("#40B0A6"), TColor::GetColor("#9285d5"), TColor::GetColor("#c2c2c2")};
          plotter.DrawDataStackedMC(tempevrate, &channelArr, channelColArr, 1, "TR", "Data (Stat. + Sys.)", 1001, "", "", false, 1);
          pad2.cd();
          plotter.axis_minimum=defaultmin;
          plotter.DrawDataMCRatio(tempevrate, tempsim, 1, true, true, 0.5, 1.5, "");
          can.SetLogx();
          pad1.cd();
          tempevrate->Draw("AXIS SAME");
          can.cd();

          can.Print((std::string("Tracker_")+prefix+"_EvRateByIntType.png").c_str());

          //Backgrounds
          pad1.cd(); 
          plotter.axis_minimum=0.001;
          plotter.DrawDataStackedMC(tempevrate, &bkgArr, channelColArr, 1, "TR", "Data (Stat. Only)", 1001, "", "", false, 1);
          pad2.cd();
          plotter.axis_minimum=defaultmin;
          plotter.DrawDataMCRatio(tempevrate, tempsim, 1, true, true, 0.5, 1.5, "");
          can.SetLogx();
          pad1.cd();
          tempevrate->Draw("AXIS SAME");
          can.cd();
          can.Print((std::string("Tracker_")+prefix+"_EvRateByBackgrounds.png").c_str());
          can.SetLogx(0);
          can.SetLogy(0);

          delete temp2p2h;
          delete tempdis;
          delete tempres;
          delete tempqe;
          delete tempChannelOther;
          delete tempWS;
          delete tempNC;
          delete tempBackgroundOther;
          delete tempsignal;
          delete tempevrate;
          delete tempsim;
        } //End plotting event rates

        //Writing out daisy weighting histograms
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

        //Plotting stuff - Should move this elsewhere long-term
        { //Plotting cross section ratios
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
          //plotter.axis_maximum = 0.4;
          std::string tmp; 
          if (prefix=="pTmu")
          {
            tmp = "Muon p_{T}";
            //tempxsec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix=="pZmu")
          {
            tmp = "Muon p_{Z}";
            //tempxsec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix=="BjorkenX") tmp = "{X}";
          else if (prefix=="Erecoil") tmp = "E_{recoil}";
          else if (prefix=="Emu")
          {
            tmp = "E_{#mu}";
            //tempxsec->GetXaxis()->SetRange(1, 11);
          } 
          plotter.axis_title_size_y=0.025;
          plotter.axis_title_offset_y=1.8;

          std::string title;

          title = tmp + " " + "Tracker Cross Section to Daisy Reweighted Cross Section ratio - Lead";
          plotter.DrawDataMCRatio(crossSection, crossSectionPb, 1, true, true, -0.95, 1.05);
          plotter.AddHistoTitle(title.c_str(), 0.028);
          can.Update();
          can.Print((std::string("Tracker_")+"LeadDaisy"+prefix+"_CrossSectionRatio.pdf").c_str());
          can.Print((std::string("Tracker_")+"LeadDaisy"+prefix+"_CrossSectionRatio.png").c_str());

          title = tmp + " " + "Tracker Cross Section to Daisy Reweighted Cross Section ratio - Iron";
          plotter.DrawDataMCRatio(crossSection, crossSectionFe, 1, true, true, -0.95, 1.05);
          plotter.AddHistoTitle(title.c_str(), 0.028);
          can.Update();
          can.Print((std::string("Tracker_")+"IronDaisy"+prefix+"_CrossSectionRatio.pdf").c_str());
          can.Print((std::string("Tracker_")+"IronDaisy"+prefix+"_CrossSectionRatio.png").c_str());

          title = tmp + " " + "Tracker Cross Section to Daisy Reweighted Cross Section ratio - Carbon";
          plotter.DrawDataMCRatio(crossSection, crossSectionC, 1, true, true, -0.95, 1.05);
          plotter.AddHistoTitle(title.c_str(), 0.028);
          can.Update();
          can.Print((std::string("Tracker_")+"CarbonDaisy"+prefix+"_CrossSectionRatio.pdf").c_str());
          can.Print((std::string("Tracker_")+"CarbonDaisy"+prefix+"_CrossSectionRatio.png").c_str());
        }

        { //Plotting daisy cross-sections - Carbon
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter(PlotUtils::kDefaultStyle);
          plotter.draw_normalized_to_bin_width=false; //Already width normlized by normalise()
          plotter.axis_label_font=43;
          plotter.axis_label_size = 25;
          plotter.headroom = 1.2;
          /* plotter.axis_title_font_y=43;
          plotter.axis_title_size_y = 35;
          plotter.axis_title_font_x=43;
          plotter.axis_title_size_x = 35;
 */
          plotter.ApplyAxisStyle();
          //plotter.axis_maximum = 0.4;



          TObjArray arr;;
          auto temp2p2h = crossSectionSimC2p2h->Clone();
          auto tempdis = crossSectionSimCDIS->Clone();
          auto tempres = crossSectionSimCRes->Clone();
          auto tempqe = crossSectionSimCQE->Clone();
          auto tempother = crossSectionSimCOther->Clone();
          auto tempsimxsec = crossSectionSimC->Clone();
          auto tempxsec = crossSectionC->Clone();
          temp2p2h->Scale(1e39);
          tempdis->Scale(1e39);
          tempres->Scale(1e39);
          tempqe->Scale(1e39);
          tempother->Scale(1e39);
          tempxsec->Scale(1e39);
          tempsimxsec->Scale(1e39);
          double dataPOTforScale = dataPOT;
          double mcPOTforScale = mcPOT;

          arr.Add(temp2p2h);
          arr.Add(tempdis);
          arr.Add(tempres);
          arr.Add(tempqe);
          arr.Add(tempother);

          double defaultmin = plotter.axis_minimum;
          plotter.legend_offset_x = 0.0825;
          plotter.legend_border_size = 0;
          plotter.legend_text_size = 0.03;

          plotter.mc_line_width=0;

          std::string tmp; 
          std::string y_title="Events#times10^{3}";
          std::cout<<"Here2.n";
          double xmax = -1111; //Default value
          std::string xaxislabel = "";
          if (prefix=="pTmu")
          {
            tmp = "Muon p_{T}";
            xaxislabel = "p_{T, #mu} [GeV/c]";
            y_title="d#sigma/dp_{t} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            xmax = 13;
            tempsimxsec->GetXaxis()->SetRange(1, xmax);
          }
          else if (prefix=="pZmu")
          {
            tmp = "Muon p_{Z}";
            xaxislabel = "p_{Z, #mu} [GeV/c]";
            y_title="d#sigma/dp_{Z} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            xmax = 13;
            tempsimxsec->GetXaxis()->SetRange(1, xmax);
          }
          else if (prefix=="BjorkenX")
          {
            tmp = "Bjorken X";
            xaxislabel = "Bjorken X";
            y_title="d#sigma/d{X} (x10^{-39}) (cm^{2}/Nucleon)";
            xmax = 7;
            tempsimxsec->GetXaxis()->SetRange(1, xmax);
            can.SetLogx();
          }
          else if (prefix=="Erecoil")
          {
            tmp = "E_{recoil}";
            y_title="d#sigma/dX (x10^{-39}) (cm^{2}/X/Nucleon)";
            xaxislabel = "E_{recoil} [GeV]";
          } 
          else if (prefix=="Emu")
          {
            tmp = "E_{#mu}";
            y_title="d#sigma/dE_{#mu} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            xaxislabel = "E_{#mu} [GeV]";
            xmax = 11;
            tempsimxsec->GetXaxis()->SetRange(1, xmax);
          } 

          std::string title = "Tracker " + tmp + " Daisy Reweighted to Carbon Flux";
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;

          plotter.DrawErrorSummary(tempsimxsec);
          can.Print((std::string("Tracker_")+prefix+"_CarbonDaisyXSecErrors.png").c_str());
          can.Clear();


          // Adjust pad heights manually (top bigger, bottom smaller)
          TPad pad1("pad1","pad1",0,0.35,1,1.0); // top pad
          TPad pad2("pad2","pad2",0,0.05,1,0.35); // bottom pad
          pad1.SetBottomMargin(0); // remove x-axis label space for top
          pad2.SetTopMargin(0);
          pad2.SetBottomMargin(0.25);
          pad1.Draw();
          pad2.Draw();

          if (prefix=="BjorkenX")
          {
            pad1.SetLogx();
            pad2.SetLogx();
          }


          std::stringstream POTstr;
          POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOTforScale;
          std::string potstr = std::string("POT Used: ")+POTstr.str();
          plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.93, 0.025, 1, 42, 13, 0);

          plotter.WritePreliminary( 0.5, 0.15, 0.035, true);
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 25, 1, 53, 22, 0);
          plotter.AddPlotLabel(y_title.c_str() , 0.08, 0.9, 0.03, 1, 42, 33, 90);
          plotter.AddPlotLabel("Data/MC" , 0.08, 0.3, 0.03, 1, 42, 33, 90);
          plotter.AddPlotLabel(xaxislabel.c_str() , 0.9, 0.08, 0.03, 1, 42, 33, 0);

          //Interaction Type
          pad1.cd(); 
          plotter.axis_minimum=0.001;

          int channelColArr[5] = {TColor::GetColor("#d2c271"), TColor::GetColor("#CA5454"), TColor::GetColor("#40B0A6"), TColor::GetColor("#9285d5"), TColor::GetColor("#c2c2c2")};
          plotter.DrawDataStackedMC(tempxsec, &arr, channelColArr, 1, "TR", "Data", 1001, "", "", false, 1, xmax);
          pad2.cd();
          plotter.axis_minimum=defaultmin;
          plotter.DrawDataMCRatio(tempxsec, tempsimxsec, 1, true, true, 0.5, 1.5, "");
          can.SetLogx();
          pad1.cd();
          tempxsec->Draw("AXIS SAME");
          can.cd();
          can.Print((std::string("Tracker_")+prefix+"_CarbonDaisyCrossSection.png").c_str());
          delete temp2p2h;
          delete tempdis;
          delete tempres;
          delete tempqe;
          delete tempother;
          delete tempsimxsec;
          delete tempxsec;
        }

        { //Plotting daisy cross-sections - Iron
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter(PlotUtils::kDefaultStyle);
          plotter.draw_normalized_to_bin_width=false; //Already width normlized by normalise()
          plotter.axis_label_font=43;
          plotter.axis_label_size = 25;
          plotter.headroom = 1.2;
          /* plotter.axis_title_font_y=43;
          plotter.axis_title_size_y = 35;
          plotter.axis_title_font_x=43;
          plotter.axis_title_size_x = 35;
 */
          plotter.ApplyAxisStyle();
          //plotter.axis_maximum = 0.4;



          TObjArray arr;;
          auto temp2p2h = crossSectionSimFe2p2h->Clone();
          auto tempdis = crossSectionSimFeDIS->Clone();
          auto tempres = crossSectionSimFeRes->Clone();
          auto tempqe = crossSectionSimFeQE->Clone();
          auto tempother = crossSectionSimFeOther->Clone();
          auto tempsimxsec = crossSectionSimFe->Clone();
          auto tempxsec = crossSectionFe->Clone();
          temp2p2h->Scale(1e39);
          tempdis->Scale(1e39);
          tempres->Scale(1e39);
          tempqe->Scale(1e39);
          tempother->Scale(1e39);
          tempxsec->Scale(1e39);
          tempsimxsec->Scale(1e39);
          double dataPOTforScale = dataPOT;
          double mcPOTforScale = mcPOT;

          arr.Add(temp2p2h);
          arr.Add(tempdis);
          arr.Add(tempres);
          arr.Add(tempqe);
          arr.Add(tempother);

          double defaultmin = plotter.axis_minimum;
          plotter.legend_offset_x = 0.0825;
          plotter.legend_border_size = 0;
          plotter.legend_text_size = 0.03;

          plotter.mc_line_width=0;

          std::string tmp; 
          std::string y_title="Events#times10^{3}";
          std::cout<<"Here2.n";
          double xmax = -1111; //Default value
          std::string xaxislabel = "";
          if (prefix=="pTmu")
          {
            tmp = "Muon p_{T}";
            xaxislabel = "p_{T, #mu} [GeV/c]";
            y_title="d#sigma/dp_{t} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            xmax = 13;
            tempsimxsec->GetXaxis()->SetRange(1, xmax);
          }
          else if (prefix=="pZmu")
          {
            tmp = "Muon p_{Z}";
            xaxislabel = "p_{Z, #mu} [GeV/c]";
            y_title="d#sigma/dp_{Z} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            xmax = 13;
            tempsimxsec->GetXaxis()->SetRange(1, xmax);
          }
          else if (prefix=="BjorkenX")
          {
            tmp = "Bjorken X";
            xaxislabel = "Bjorken X";
            y_title="d#sigma/d{X} (x10^{-39}) (cm^{2}/Nucleon)";
            xmax = 7;
            tempsimxsec->GetXaxis()->SetRange(1, xmax);
            can.SetLogx();
          }
          else if (prefix=="Erecoil")
          {
            tmp = "E_{recoil}";
            y_title="d#sigma/dX (x10^{-39}) (cm^{2}/X/Nucleon)";
            xaxislabel = "E_{recoil} [GeV]";
          } 
          else if (prefix=="Emu")
          {
            tmp = "E_{#mu}";
            y_title="d#sigma/dE_{#mu} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            xaxislabel = "E_{#mu} [GeV]";
            xmax = 11;
            tempsimxsec->GetXaxis()->SetRange(1, xmax);
          } 

          std::string title = "Tracker " + tmp + " Daisy Reweighted to Iron Flux";
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;

          plotter.DrawErrorSummary(tempsimxsec);
          can.Print((prefix+"_IronDaisyXSecErrors.png").c_str());
          can.Clear();


          // Adjust pad heights manually (top bigger, bottom smaller)
          TPad pad1("pad1","pad1",0,0.35,1,1.0); // top pad
          TPad pad2("pad2","pad2",0,0.05,1,0.35); // bottom pad
          pad1.SetBottomMargin(0); // remove x-axis label space for top
          pad2.SetTopMargin(0);
          pad2.SetBottomMargin(0.25);
          pad1.Draw();
          pad2.Draw();

          if (prefix=="BjorkenX")
          {
            pad1.SetLogx();
            pad2.SetLogx();
          }


          std::stringstream POTstr;
          POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOTforScale;
          std::string potstr = std::string("POT Used: ")+POTstr.str();
          plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.93, 0.025, 1, 42, 13, 0);

          plotter.WritePreliminary( 0.5, 0.15, 0.035, true);
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 25, 1, 53, 22, 0);
          plotter.AddPlotLabel(y_title.c_str() , 0.08, 0.9, 0.03, 1, 42, 33, 90);
          plotter.AddPlotLabel("Data/MC" , 0.08, 0.3, 0.03, 1, 42, 33, 90);
          plotter.AddPlotLabel(xaxislabel.c_str() , 0.9, 0.08, 0.03, 1, 42, 33, 0);

          //Interaction Type
          pad1.cd(); 
          plotter.axis_minimum=0.001;

          int channelColArr[5] = {TColor::GetColor("#d2c271"), TColor::GetColor("#CA5454"), TColor::GetColor("#40B0A6"), TColor::GetColor("#9285d5"), TColor::GetColor("#c2c2c2")};
          plotter.DrawDataStackedMC(tempxsec, &arr, channelColArr, 1, "TR", "Data", 1001, "", "", false, 1, xmax);
          pad2.cd();
          plotter.axis_minimum=defaultmin;
          plotter.DrawDataMCRatio(tempxsec, tempsimxsec, 1, true, true, 0.5, 1.5, "");
          can.SetLogx();
          pad1.cd();
          tempxsec->Draw("AXIS SAME");
          can.cd();
          can.Print((std::string("Tracker_")+prefix+"_IronDaisyCrossSection.png").c_str());
          delete temp2p2h;
          delete tempdis;
          delete tempres;
          delete tempqe;
          delete tempother;
          delete tempsimxsec;
          delete tempxsec;
        }

        { //Plotting daisy cross-sections - Lead
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter(PlotUtils::kDefaultStyle);
          plotter.draw_normalized_to_bin_width=false; //Already width normlized by normalise()
          plotter.axis_label_font=43;
          plotter.axis_label_size = 25;
          plotter.headroom = 1.2;
          /* plotter.axis_title_font_y=43;
          plotter.axis_title_size_y = 35;
          plotter.axis_title_font_x=43;
          plotter.axis_title_size_x = 35;
 */
          plotter.ApplyAxisStyle();
          //plotter.axis_maximum = 0.4;



          TObjArray arr;
          auto temp2p2h = crossSectionSimPb2p2h->Clone();
          auto tempdis = crossSectionSimPbDIS->Clone();
          auto tempres = crossSectionSimPbRes->Clone();
          auto tempqe = crossSectionSimPbQE->Clone();
          auto tempother = crossSectionSimPbOther->Clone();
          auto tempsimxsec = crossSectionSimPb->Clone();
          auto tempxsec = crossSectionPb->Clone();
          temp2p2h->Scale(1e39);
          tempdis->Scale(1e39);
          tempres->Scale(1e39);
          tempqe->Scale(1e39);
          tempother->Scale(1e39);
          tempxsec->Scale(1e39);
          tempsimxsec->Scale(1e39);
          double dataPOTforScale = dataPOT;
          double mcPOTforScale = mcPOT;

          arr.Add(temp2p2h);
          arr.Add(tempdis);
          arr.Add(tempres);
          arr.Add(tempqe);
          arr.Add(tempother);

          double defaultmin = plotter.axis_minimum;
          plotter.legend_offset_x = 0.0825;
          plotter.legend_border_size = 0;
          plotter.legend_text_size = 0.03;

          plotter.mc_line_width=0;

          std::string tmp; 
          std::string y_title="Events#times10^{3}";
          std::cout<<"Here2.n";
          double xmax = -1111; //Default value
          std::string xaxislabel = "";
          if (prefix=="pTmu")
          {
            tmp = "Muon p_{T}";
            xaxislabel = "p_{T, #mu} [GeV/c]";
            y_title="d#sigma/dp_{t} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            xmax = 13;
            tempsimxsec->GetXaxis()->SetRange(1, xmax);
          }
          else if (prefix=="pZmu")
          {
            tmp = "Muon p_{Z}";
            xaxislabel = "p_{Z, #mu} [GeV/c]";
            y_title="d#sigma/dp_{Z} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            xmax = 13;
            tempsimxsec->GetXaxis()->SetRange(1, xmax);
          }
          else if (prefix=="BjorkenX")
          {
            tmp = "Bjorken X";
            xaxislabel = "Bjorken X";
            y_title="d#sigma/d{X} (x10^{-39}) (cm^{2}/Nucleon)";
            xmax = 7;
            tempsimxsec->GetXaxis()->SetRange(1, xmax);
            can.SetLogx();
          }
          else if (prefix=="Erecoil")
          {
            tmp = "E_{recoil}";
            y_title="d#sigma/dX (x10^{-39}) (cm^{2}/X/Nucleon)";
            xaxislabel = "E_{recoil} [GeV]";
          } 
          else if (prefix=="Emu")
          {
            tmp = "E_{#mu}";
            y_title="d#sigma/dE_{#mu} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            xaxislabel = "E_{#mu} [GeV]";
            xmax = 11;
            tempsimxsec->GetXaxis()->SetRange(1, xmax);
          } 

          std::string title = "Tracker " + tmp + " Daisy Reweighted to Lead Flux";
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;

          plotter.DrawErrorSummary(tempsimxsec);
          can.Print((std::string("Tracker_")+prefix+"_LeadDaisyXSecErrors.png").c_str());
          can.Clear();


          // Adjust pad heights manually (top bigger, bottom smaller)
          TPad pad1("pad1","pad1",0,0.35,1,1.0); // top pad
          TPad pad2("pad2","pad2",0,0.05,1,0.35); // bottom pad
          pad1.SetBottomMargin(0); // remove x-axis label space for top
          pad2.SetTopMargin(0);
          pad2.SetBottomMargin(0.25);
          pad1.Draw();
          pad2.Draw();

          if (prefix=="BjorkenX")
          {
            pad1.SetLogx();
            pad2.SetLogx();
          }


          std::stringstream POTstr;
          POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOTforScale;
          std::string potstr = std::string("POT Used: ")+POTstr.str();
          plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.93, 0.025, 1, 42, 13, 0);

          plotter.WritePreliminary( 0.5, 0.15, 0.035, true);
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 25, 1, 53, 22, 0);
          plotter.AddPlotLabel(y_title.c_str() , 0.08, 0.9, 0.03, 1, 42, 33, 90);
          plotter.AddPlotLabel("Data/MC" , 0.08, 0.3, 0.03, 1, 42, 33, 90);
          plotter.AddPlotLabel(xaxislabel.c_str() , 0.9, 0.08, 0.03, 1, 42, 33, 0);

          //Interaction Type
          pad1.cd(); 
          plotter.axis_minimum=0.001;

          int channelColArr[5] = {TColor::GetColor("#d2c271"), TColor::GetColor("#CA5454"), TColor::GetColor("#40B0A6"), TColor::GetColor("#9285d5"), TColor::GetColor("#c2c2c2")};
          plotter.DrawDataStackedMC(tempxsec, &arr, channelColArr, 1, "TR", "Data", 1001, "", "", false, 1, xmax);
          pad2.cd();
          plotter.axis_minimum=defaultmin;
          plotter.DrawDataMCRatio(tempxsec, tempsimxsec, 1, true, true, 0.5, 1.5, "");
          can.SetLogx();
          pad1.cd();
          tempxsec->Draw("AXIS SAME");
          can.cd();
          can.Print((prefix+"_LeadDaisyCrossSection.png").c_str());
          delete temp2p2h;
          delete tempdis;
          delete tempres;
          delete tempqe;
          delete tempother;
          delete tempsimxsec;
          delete tempxsec;
        }


        { //Plotting cross sections
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
          plotter.draw_normalized_to_bin_width=false; //Already width normlized by normalise()
          plotter.mc_line_width=0;
          //plotter.axis_maximum = 0.4;
          TObjArray arr;
          auto temp2p2h = crossSectionSimPb2p2h->Clone();
          auto tempdis = crossSectionSimPbDIS->Clone();
          auto tempres = crossSectionSimPbRes->Clone();
          auto tempqe = crossSectionSimPbQE->Clone();
          auto tempother = crossSectionSimPbOther->Clone();
          auto tempsimxsec = crossSectionSimPb->Clone();
          auto tempxsec = crossSectionPb->Clone();

          temp2p2h->Scale(1e39);
          tempdis->Scale(1e39);
          tempres->Scale(1e39);
          tempqe->Scale(1e39);
          tempother->Scale(1e39);
          tempxsec->Scale(1e39);
          arr.Add(temp2p2h);
          arr.Add(tempdis);
          arr.Add(tempres);
          arr.Add(tempqe);
          arr.Add(tempother);

          std::string tmp; 
          std::string y_title;
          if (prefix.find("pTmu")!=std::string::npos)
          {
            tmp = "Muon p_{T}";
            y_title="d#sigma/dp_{t} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            tempxsec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix.find("pZmu")!=std::string::npos)
          {
            tmp = "Muon p_{Z}";
            y_title="d#sigma/dp_{Z} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            tempxsec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix.find("BjorkenX")!=std::string::npos)
          {
            tmp = "Bjorken X";
            y_title="d#sigma/dX (x10^{-39}) (cm^{2}/X/Nucleon)";
            temp2p2h->GetXaxis()->SetRange(1, 7);
            can.SetLogx();
          }
          else if (prefix.find("Erecoil")!=std::string::npos)
          {
            tmp = "E_{recoil}";
            y_title="d#sigma/dX (x10^{-39}) (cm^{2}/X/Nucleon)";
          } 
          else if (prefix.find("Emu")!=std::string::npos)
          {
            tmp = "E_{#mu}";
            y_title="d#sigma/dE_{#mu} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            tempxsec->GetXaxis()->SetRange(1, 11);
          } 

          std::string title = tmp + " Daisy Reweighted to Lead Flux";
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;
          int channelColArr[5] = {TColor::GetColor("#d2c271"), TColor::GetColor("#CA5454"), TColor::GetColor("#40B0A6"), TColor::GetColor("#9285d5"), TColor::GetColor("#c2c2c2")};
          plotter.DrawDataStackedMC(tempxsec, &arr, channelColArr, 1, "TR", "Data", 1001, "", "", false, 1);
          plotter.AddPlotLabel(y_title.c_str() , 0.05, 0.9, 0.04, 1, 42, 33, 90);
          //plotter.DrawDataStackedMC(tempxsec, arr, 1, "TR", "Data", 2, 1, 3001, "", y_title.c_str());
          plotter.AddHistoTitle(title.c_str(), 0.039);
          plotter.AddPlotLabel("p_{Z, #mu}" , 0.87, 0.07, 0.04, 1, 42, 33);
          can.Update();
          can.Print((std::string("Tracker_")+prefix+"_PbDaisyCrossSection.pdf").c_str());
          can.Print((std::string("Tracker_")+prefix+"_PbDaisyCrossSection.png").c_str());
          delete temp2p2h;
          delete tempdis;
          delete tempres;
          delete tempqe;
          delete tempother;
          delete tempsimxsec;
          delete tempxsec;
        }
        { //Plotting cross sections
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
          plotter.draw_normalized_to_bin_width=false; //Already width normlized by normalise()
          plotter.mc_line_width=0;
          //plotter.axis_maximum = 0.4;
          TObjArray arr;
          auto temp2p2h = crossSectionSimFe2p2h->Clone();
          auto tempdis = crossSectionSimFeDIS->Clone();
          auto tempres = crossSectionSimFeRes->Clone();
          auto tempqe = crossSectionSimFeQE->Clone();
          auto tempother = crossSectionSimFeOther->Clone();
          auto tempxsec = crossSectionFe->Clone();
          temp2p2h->Scale(1e39);
          tempdis->Scale(1e39);
          tempres->Scale(1e39);
          tempqe->Scale(1e39);
          tempother->Scale(1e39);
          tempxsec->Scale(1e39);


          std::string tmp; 
          std::string y_title;
          if (prefix.find("pTmu")!=std::string::npos)
          {
            tmp = "Muon p_{T}";
            y_title="d#sigma/dp_{t} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            tempxsec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix.find("pZmu")!=std::string::npos)
          {
            tmp = "Muon p_{Z}";
            y_title="d#sigma/dp_{Z} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            tempxsec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix.find("BjorkenX")!=std::string::npos)
          {
            tmp = "Bjorken X";
            y_title="d#sigma/dX (x10^{-39}) (cm^{2}/X/Nucleon)";
            temp2p2h->GetXaxis()->SetRange(1, 7);
            can.SetLogx();
          }
          else if (prefix.find("Erecoil")!=std::string::npos)
          {
            tmp = "E_{recoil}";
            y_title="d#sigma/dX (x10^{-39}) (cm^{2}/X/Nucleon)";
          } 
          else if (prefix.find("Emu")!=std::string::npos)
          {
            tmp = "E_{#mu}";
            y_title="d#sigma/dE_{#mu} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            tempxsec->GetXaxis()->SetRange(1, 11);
          } 
          arr.Add(temp2p2h);
          arr.Add(tempdis);
          arr.Add(tempres);
          arr.Add(tempqe);
          arr.Add(tempother);
          std::string title = tmp + " Daisy Reweighted to Iron Flux";
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;
          int channelColArr[5] = {TColor::GetColor("#d2c271"), TColor::GetColor("#CA5454"), TColor::GetColor("#40B0A6"), TColor::GetColor("#9285d5"), TColor::GetColor("#c2c2c2")};
          plotter.DrawDataStackedMC(tempxsec, &arr, channelColArr, 1, "TR", "Data (Stat. Only)", 1001, "", "", false, 1);
          plotter.AddHistoTitle(title.c_str(), 0.039);
          plotter.AddPlotLabel(y_title.c_str() , 0.05, 0.9, 0.04, 1, 42, 33, 90);
          plotter.AddPlotLabel("p_{Z, #mu}" , 0.87, 0.07, 0.04, 1, 42, 33);
          can.Update();
          can.Print((std::string("Tracker_")+prefix+"_FeDaisyCrossSection.pdf").c_str());
          can.Print((std::string("Tracker_")+prefix+"_FeDaisyCrossSection.png").c_str());
          delete temp2p2h;
          delete tempdis;
          delete tempres;
          delete tempqe;
          delete tempother;
          delete tempxsec;
        }
        { //Plotting cross sections
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
          plotter.draw_normalized_to_bin_width=false; //Already width normlized by normalise()
          plotter.mc_line_width=0;
          //plotter.axis_maximum = 0.4;
          TObjArray arr;
          auto temp2p2h = crossSectionSimC2p2h->Clone();
          auto tempdis = crossSectionSimCDIS->Clone();
          auto tempres = crossSectionSimCRes->Clone();
          auto tempqe = crossSectionSimCQE->Clone();
          auto tempother = crossSectionSimCOther->Clone();
          auto tempxsec = crossSectionC->Clone();
          temp2p2h->Scale(1e39);
          tempdis->Scale(1e39);
          tempres->Scale(1e39);
          tempqe->Scale(1e39);
          tempother->Scale(1e39);
          tempxsec->Scale(1e39);

          std::string tmp; 
          std::string y_title;
          if (prefix.find("pTmu")!=std::string::npos)
          {
            tmp = "Muon p_{T}";
            y_title="d#sigma/dp_{t} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            tempxsec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix.find("pZmu")!=std::string::npos)
          {
            tmp = "Muon p_{Z}";
            y_title="d#sigma/dp_{Z} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            tempxsec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix.find("BjorkenX")!=std::string::npos)
          {
            tmp = "Bjorken X";
            y_title="d#sigma/dX (x10^{-39}) (cm^{2}/X/Nucleon)";
            temp2p2h->GetXaxis()->SetRange(1, 7);
            can.SetLogx();
          }
          else if (prefix.find("Erecoil")!=std::string::npos)
          {
            tmp = "E_{recoil}";
            y_title="d#sigma/dX (x10^{-39}) (cm^{2}/X/Nucleon)";
          } 
          else if (prefix.find("Emu")!=std::string::npos)
          {
            tmp = "E_{#mu}";
            y_title="d#sigma/dE_{#mu} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            tempxsec->GetXaxis()->SetRange(1, 11);
          } 
          arr.Add(temp2p2h);
          arr.Add(tempdis);
          arr.Add(tempres);
          arr.Add(tempqe);
          arr.Add(tempother);
          std::string title = tmp + " Daisy Reweighted to Carbon Flux";
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;
          int channelColArr[5] = {TColor::GetColor("#d2c271"), TColor::GetColor("#CA5454"), TColor::GetColor("#40B0A6"), TColor::GetColor("#9285d5"), TColor::GetColor("#c2c2c2")};
          plotter.DrawDataStackedMC(tempxsec, &arr, channelColArr, 1, "TR", "Data (Stat. Only)", 1001, "", "", false, 1);
          plotter.AddHistoTitle(title.c_str(), 0.039);
          plotter.AddPlotLabel(y_title.c_str() , 0.05, 0.9, 0.04, 1, 42, 33, 90);
          plotter.AddPlotLabel("p_{Z, #mu}" , 0.87, 0.07, 0.04, 1, 42, 33);
          can.Update();
          can.Print((std::string("Tracker_")+prefix+"_CDaisyCrossSection.pdf").c_str());
          can.Print((std::string("Tracker_")+prefix+"_CDaisyCrossSection.png").c_str());
          delete temp2p2h;
          delete tempdis;
          delete tempres;
          delete tempqe;
          delete tempother;
          delete tempxsec;
        }
        { //Plotting cross sections
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
          plotter.draw_normalized_to_bin_width=false; //Already width normlized by normalise()
          plotter.mc_line_width=0;
          //plotter.axis_maximum = 0.4;
          TObjArray arr;

          auto temp2p2h = simulatedCrossSection2P2H->Clone();
          auto tempdis = simulatedCrossSectionDIS->Clone();
          auto tempres = simulatedCrossSectionRES->Clone();
          auto tempqe = simulatedCrossSectionQE->Clone();
          auto tempother = simulatedCrossSectionOther->Clone();
          auto tempxsec = crossSection->Clone();
          temp2p2h->Scale(1e39);
          tempdis->Scale(1e39);
          tempres->Scale(1e39);
          tempqe->Scale(1e39);
          tempother->Scale(1e39);
          tempxsec->Scale(1e39);

          std::string tmp; 
          std::string y_title;
          if (prefix.find("pTmu")!=std::string::npos)
          {
            tmp = "Muon p_{T}";
            y_title="d#sigma/dp_{t} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            tempxsec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix.find("pZmu")!=std::string::npos)
          {
            tmp = "Muon p_{Z}";
            y_title="d#sigma/dp_{Z} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            tempxsec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix.find("BjorkenX")!=std::string::npos)
          {
            tmp = "Bjorken X";
            y_title="d#sigma/dX (x10^{-39}) (cm^{2}/X/Nucleon)";
            tempxsec->GetXaxis()->SetRange(1, 7);
            can.SetLogx();
          }
          else if (prefix.find("Erecoil")!=std::string::npos)
          {
            tmp = "E_{recoil}";
            y_title="d#sigma/dX (x10^{-39}) (cm^{2}/X/Nucleon)";
          } 
          else if (prefix.find("Emu")!=std::string::npos)
          {
            tmp = "E_{#mu}";
            y_title="d#sigma/dE_{#mu} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            tempxsec->GetXaxis()->SetRange(1, 11);
          } 
          
          arr.Add(temp2p2h);
          arr.Add(tempdis);
          arr.Add(tempres);
          arr.Add(tempqe);
          arr.Add(tempother);

          std::string title = tmp + " - No Daisy Reweighting Applied";
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;
          int channelColArr[5] = {TColor::GetColor("#d2c271"), TColor::GetColor("#CA5454"), TColor::GetColor("#40B0A6"), TColor::GetColor("#9285d5"), TColor::GetColor("#c2c2c2")};
          plotter.DrawDataStackedMC(tempxsec, &arr, channelColArr, 1, "TR", "Data (Stat. Only)", 1001, "", "", false, 1);
          plotter.AddHistoTitle(title.c_str(), 0.039);
          plotter.AddPlotLabel(y_title.c_str() , 0.05, 0.9, 0.04, 1, 42, 33, 90);
          plotter.AddPlotLabel("p_{Z, #mu}" , 0.87, 0.07, 0.04, 1, 42, 33);
          can.Update();
          can.Print((std::string("Tracker_")+prefix+"_NoDaisyCrossSection.pdf").c_str());
          can.Print((std::string("Tracker_")+prefix+"_NoDaisyCrossSection.png").c_str());
          delete temp2p2h;
          delete tempdis;
          delete tempres;
          delete tempqe;
          delete tempother;
          delete tempxsec;
        }

        for (int i = 0; i<12; i++)
        { //Plotting migration
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kDefaultStyle);
          //plotter.axis_maximum = 0.4;

          plotter.axis_title_size_y=0.03;
          plotter.axis_title_size_x=0.03;
          plotter.axis_title_offset_x=1.8;
          plotter.axis_title_offset_y=1.8;
          plotter.SetRedHeatPalette();

          std::string tmp = "";
         /*  if (prefix=="pTmu")
          {
            tmp = "Muon p_{T}";
            xaxislabel = "p_{T, #mu} [GeV/c]";
            migration->GetXaxis()->SetRange(1, 13);
            migration->GetYaxis()->SetRange(1, 13);
          }
          else if (prefix=="pZmu")
          {
            tmp = "Muon p_{Z}";
            xaxislabel = "p_{Z, #mu} [GeV/c]";
            migration->GetXaxis()->SetRange(1, 13);
            migration->GetYaxis()->SetRange(1, 13);
          }
          else if (prefix=="BjorkenX")
          {
            tmp = "Bjorken X";
            xaxislabel = "X";
            migration->GetXaxis()->SetRange(1, 9);
            migration->GetYaxis()->SetRange(1, 9);
            can.SetLogx();
            can.SetLogy();
          }
          else if (prefix=="Erecoil")
          {
            tmp = "E_{recoil}";
            xaxislabel = "E_{recoil} [GeV]";
          } 
          else if (prefix=="Emu")
          {
            tmp = "E_{#mu}";
            xaxislabel = "E_{#mu} [GeV]";
            xaxislabel = "X";
            migration->GetXaxis()->SetRange(1, 11);
            migration->GetYaxis()->SetRange(1, 11);
          } */

          std::string xaxislabel = "";
          std::string yaxislabel = "";
          if (prefix=="pTmu")
          {
            tmp = "Muon p_{T}";
            yaxislabel = "True p_{T, #mu} [GeV/c]";
            xaxislabel = "p_{T, #mu} [GeV/c]";
            //migration->GetXaxis()->SetRange(1, 13);
            //migration->GetYaxis()->SetRange(1, 13);
          }
          else if (prefix=="pZmu")
          {
            tmp = "Muon p_{Z}";
            yaxislabel = "True p_{Z, #mu} [GeV/c]";
            xaxislabel = "p_{Z, #mu} [GeV/c]";
            //migration->GetXaxis()->SetRange(1, 13);
            //migration->GetYaxis()->SetRange(1, 13);
          }
          else if (prefix=="BjorkenX")
          {
            tmp = "Bjorken X";
            yaxislabel = "True X";
            xaxislabel = "X";
            can.SetLogx();
            can.SetLogy();
            //migration->GetXaxis()->SetRange(1, 6);
            //migration->GetYaxis()->SetRange(1, 6);

            //tempevrate->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix=="Erecoil")
          {
            tmp = "E_{recoil}";
            yaxislabel = "Truth E_{recoil} [GeV]";
            xaxislabel = "E_{recoil} [GeV]";
          } 
          else if (prefix=="Emu")
          {
            tmp = "E_{#mu}";
            yaxislabel = "True E_{#mu} [GeV]";
            xaxislabel = "E_{#mu} [GeV]";
            //migration->GetXaxis()->SetRange(1, 11);
            //migration->GetYaxis()->SetRange(1, 11);
          } 

          //Int_t colors[] = {0, 1};
          //gStyle->SetPalette(2, colors);
          plotter.axis_title_offset_z = 1.25;
          plotter.DrawNormalizedMigrationHistogram(DaisyMigration[i].get(), false, false, true, false, 0.02);
          plotter.WritePreliminary( 0.65, 0.2, 0.035, true);
          std::string title = "Tracker Petal " +std::to_string(i)+ " " + tmp;
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 25, 1, 53, 22, 0);
          plotter.AddPlotLabel(yaxislabel.c_str() , 0.08, 0.9, 0.03, 1, 62, 33, 90);
          plotter.AddPlotLabel(xaxislabel.c_str() , 0.9, 0.08, 0.03, 1, 62, 33, 0);

          //plotter.AddHistoTitle(title.c_str(), 0.039);
          can.Update();
          //can.Print((tgtname+"_"+prefix+"_DSSidebandScaled.pdf").c_str());
          can.Print(std::string(std::string("Tracker_")+prefix+"_petal_"+i+"_migration.png").c_str());
          can.Print(std::string(std::string("Tracker_")+prefix+"_petal_"+i+"_migration.pdf").c_str());
          can.SetLogx(0);
          can.SetLogy(0);
        }
        { //Plotting migration
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kDefaultStyle);
          //plotter.axis_maximum = 0.4;

          plotter.axis_title_size_y=0.03;
          plotter.axis_title_size_x=0.03;
          plotter.axis_title_offset_x=1.8;
          plotter.axis_title_offset_y=1.8;
          plotter.SetRedHeatPalette();

          std::string tmp = "";
         /*  if (prefix=="pTmu")
          {
            tmp = "Muon p_{T}";
            xaxislabel = "p_{T, #mu} [GeV/c]";
            migration->GetXaxis()->SetRange(1, 13);
            migration->GetYaxis()->SetRange(1, 13);
          }
          else if (prefix=="pZmu")
          {
            tmp = "Muon p_{Z}";
            xaxislabel = "p_{Z, #mu} [GeV/c]";
            migration->GetXaxis()->SetRange(1, 13);
            migration->GetYaxis()->SetRange(1, 13);
          }
          else if (prefix=="BjorkenX")
          {
            tmp = "Bjorken X";
            xaxislabel = "X";
            migration->GetXaxis()->SetRange(1, 9);
            migration->GetYaxis()->SetRange(1, 9);
            can.SetLogx();
            can.SetLogy();
          }
          else if (prefix=="Erecoil")
          {
            tmp = "E_{recoil}";
            xaxislabel = "E_{recoil} [GeV]";
          } 
          else if (prefix=="Emu")
          {
            tmp = "E_{#mu}";
            xaxislabel = "E_{#mu} [GeV]";
            xaxislabel = "X";
            migration->GetXaxis()->SetRange(1, 11);
            migration->GetYaxis()->SetRange(1, 11);
          } */

          std::string xaxislabel = "";
          std::string yaxislabel = "";
          if (prefix=="pTmu")
          {
            tmp = "Muon p_{T}";
            yaxislabel = "True p_{T, #mu} [GeV/c]";
            xaxislabel = "p_{T, #mu} [GeV/c]";
            //migration->GetXaxis()->SetRange(1, 13);
            //migration->GetYaxis()->SetRange(1, 13);
          }
          else if (prefix=="pZmu")
          {
            tmp = "Muon p_{Z}";
            yaxislabel = "True p_{Z, #mu} [GeV/c]";
            xaxislabel = "p_{Z, #mu} [GeV/c]";
            //migration->GetXaxis()->SetRange(1, 13);
            //migration->GetYaxis()->SetRange(1, 13);
          }
          else if (prefix=="BjorkenX")
          {
            tmp = "Bjorken X";
            yaxislabel = "True X";
            xaxislabel = "X";
            can.SetLogx();
            can.SetLogy();
            //migration->GetXaxis()->SetRange(1, 6);
            //migration->GetYaxis()->SetRange(1, 6);

            //tempevrate->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix=="Erecoil")
          {
            tmp = "E_{recoil}";
            yaxislabel = "Truth E_{recoil} [GeV]";
            xaxislabel = "E_{recoil} [GeV]";
          } 
          else if (prefix=="Emu")
          {
            tmp = "E_{#mu}";
            yaxislabel = "True E_{#mu} [GeV]";
            xaxislabel = "E_{#mu} [GeV]";
            //migration->GetXaxis()->SetRange(1, 11);
            //migration->GetYaxis()->SetRange(1, 11);
          } 

          //Int_t colors[] = {0, 1};
          //gStyle->SetPalette(2, colors);
          plotter.axis_title_offset_z = 1.25;
          plotter.DrawNormalizedMigrationHistogram(migration, false, false, true, false, 0.02);
          plotter.WritePreliminary( 0.65, 0.2, 0.035, true);
          std::string title = "Tracker " + tmp;
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 25, 1, 53, 22, 0);
          plotter.AddPlotLabel(yaxislabel.c_str() , 0.08, 0.9, 0.03, 1, 62, 33, 90);
          plotter.AddPlotLabel(xaxislabel.c_str() , 0.9, 0.08, 0.03, 1, 62, 33, 0);

          //plotter.AddHistoTitle(title.c_str(), 0.039);
          can.Update();
          //can.Print((tgtname+"_"+prefix+"_DSSidebandScaled.pdf").c_str());
          can.Print(std::string(std::string("Tracker_")+prefix+"_migration.png").c_str());
          can.Print(std::string(std::string("Tracker_")+prefix+"_migration.pdf").c_str());
          can.SetLogx(0);
          can.SetLogy(0);
        }









        if (prefix!="segment")
        { //Plotting cross-sections
          //GENIE3 G18_02a
          std::string G18_02a_02_11a_string = std::string("/exp/minerva/app/users/alhart/MAT_AL9/MINERvA-101-Cross-Section/modelComparions/")+nutypestr+"/Tracker_GENIE_G18_02a_02_11a_cross_sections.root";
          auto G18_02a_02_11a_file = TFile::Open(G18_02a_02_11a_string.c_str());
          if (!G18_02a_02_11a_file)
          {
            std::cerr << "Failed to open model comparison file " << G18_02a_02_11a_string << "\n";
            return 2;
          }
          //GENIE3 G18_02a
          std::string G18_02b_02_11a_string = std::string("/exp/minerva/app/users/alhart/MAT_AL9/MINERvA-101-Cross-Section/modelComparions/")+nutypestr+"/Tracker_GENIE_G18_02b_02_11a_cross_sections.root";
          auto G18_02b_02_11a_file = TFile::Open(G18_02b_02_11a_string.c_str());
          if (!G18_02b_02_11a_file)
          {
            std::cerr << "Failed to open model comparison file " << G18_02b_02_11a_string << "\n";
            return 2;
          }
          //GENIE3 G18_10a
          std::string G18_10a_02_11a_string = std::string("/exp/minerva/app/users/alhart/MAT_AL9/MINERvA-101-Cross-Section/modelComparions/")+nutypestr+"/Tracker_GENIE_G18_10a_02_11a_cross_sections.root";
          auto G18_10a_02_11a_file = TFile::Open(G18_10a_02_11a_string.c_str());
          if (!G18_10a_02_11a_file)
          {
            std::cerr << "Failed to open model comparison file " << G18_10a_02_11a_string << "\n";
            return 2;
          }
          //GENIE3 G18_10b
          std::string G18_10b_02_11a_string = std::string("/exp/minerva/app/users/alhart/MAT_AL9/MINERvA-101-Cross-Section/modelComparions/")+nutypestr+"/Tracker_GENIE_G18_10b_02_11a_cross_sections.root";
          auto G18_10b_02_11a_file = TFile::Open(G18_10b_02_11a_string.c_str());
          if (!G18_10b_02_11a_file)
          {
            std::cerr << "Failed to open model comparison file " << G18_10b_02_11a_string << "\n";
            return 2;
          }
          //NEUT LFG ma105
          std::string NEUT_tune_LFG_maqe1_05_string = std::string("/exp/minerva/app/users/alhart/MAT_AL9/MINERvA-101-Cross-Section/modelComparions/")+nutypestr+"/Tracker_NEUT_SF_ma103_cross_sections.root";
          auto NEUT_tune_LFG_maqe1_05_file = TFile::Open(NEUT_tune_LFG_maqe1_05_string.c_str());
          if (!NEUT_tune_LFG_maqe1_05_file)
          {
            std::cerr << "Failed to open model comparison file " << NEUT_tune_LFG_maqe1_05_string << "\n";
            return 2;
          }
          //NEUT SFG ma103
          std::string NEUT_tune_SFG_maqe1_03_string = std::string("/exp/minerva/app/users/alhart/MAT_AL9/MINERvA-101-Cross-Section/modelComparions/")+nutypestr+"/Tracker_NEUT_LFG_ma105_cross_sections.root";
          auto NEUT_tune_SFG_maqe1_03_file = TFile::Open(NEUT_tune_SFG_maqe1_03_string.c_str());
          if (!NEUT_tune_SFG_maqe1_03_file)
          {
            std::cerr << "Failed to open model comparison file " << NEUT_tune_SFG_maqe1_03_string << "\n";
            return 2;
          }
          std::cout<<"Model comparison cross sections found - plotting ratios\n";
          std::string titlebase;
          TH1D* G18_02a_02_11a_xsec=(TH1D*)util::GetIngredientSmart<TH1D>(*G18_02a_02_11a_file, prefix+"_rw").get()->Clone("GENIE3 G18_02a");
          TH1D* G18_02b_02_11a_xsec=(TH1D*)util::GetIngredientSmart<TH1D>(*G18_02b_02_11a_file, prefix+"_rw").get()->Clone("GENIE3 G18_02b");
          TH1D* G18_10a_02_11a_xsec=(TH1D*)util::GetIngredientSmart<TH1D>(*G18_10a_02_11a_file, prefix+"_rw").get()->Clone("GENIE3 G18_10a");
          TH1D* G18_10b_02_11a_xsec=(TH1D*)util::GetIngredientSmart<TH1D>(*G18_10b_02_11a_file, prefix+"_rw").get()->Clone("GENIE3 G18_10b");
          TH1D* NEUT_tune_LFG_maqe1_05_xsec=(TH1D*)util::GetIngredientSmart<TH1D>(*NEUT_tune_LFG_maqe1_05_file, prefix+"_rw").get()->Clone("NEUT LFG ma=105");


          G18_02a_02_11a_xsec->Scale(1., "width");
          G18_02b_02_11a_xsec->Scale(1., "width");
          G18_10a_02_11a_xsec->Scale(1., "width");
          G18_10b_02_11a_xsec->Scale(1., "width");
          NEUT_tune_LFG_maqe1_05_xsec->Scale(1., "width");

          std::cout<<"Model comparison cross sections found - plotting ratios1\n";
          PlotUtils::MnvH1D hist_G18_02a(*G18_02a_02_11a_xsec);
          PlotUtils::MnvH1D hist_G18_02b(*G18_02b_02_11a_xsec);
          PlotUtils::MnvH1D hist_G18_10a(*G18_10a_02_11a_xsec);
          PlotUtils::MnvH1D hist_G18_10b(*G18_10b_02_11a_xsec);
          PlotUtils::MnvH1D hist_NEUT_LFG(*NEUT_tune_LFG_maqe1_05_xsec);

          std::cout<<"Model comparison cross sections found - plotting ratios2\n";
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kDefaultStyle);
          plotter.draw_normalized_to_bin_width=false; //Already width normlized by normalise()
          plotter.axis_label_font=43;
          plotter.axis_label_size = 25;
          plotter.headroom = 1.2;
          /* plotter.axis_title_font_y=43;
          plotter.axis_title_size_y = 35;
          plotter.axis_title_font_x=43;
          plotter.axis_title_size_x = 35;
          */
          plotter.ApplyAxisStyle();
          //plotter.axis_maximum = 0.4;

          std::cout<<"1.111\n";
          std::cout<<"1.11.11.\n";
          TObjArray* arr = new TObjArray();
          std::cout<<"1.11.1.\n";
          auto temp2p2h = simulatedCrossSection2P2H->Clone();
          std::cout<<"1.111.\n";
          auto tempdis = simulatedCrossSectionDIS->Clone();
          std::cout<<"1.111..\n";
          auto tempres = simulatedCrossSectionRES->Clone();
          std::cout<<"1.111...\n";
          auto tempqe = simulatedCrossSectionQE->Clone();
          std::cout<<"1.111....\n";
          auto tempother = simulatedCrossSectionOther->Clone();
          std::cout<<"1.112\n";
          auto tempsimxsec = temp2p2h->Clone();
          tempsimxsec->Add(tempdis);
          tempsimxsec->Add(tempres);
          tempsimxsec->Add(tempqe);
          tempsimxsec->Add(tempother);
          auto tempxsec = crossSection->Clone();
          std::cout<<"1.113\n";
          double defaultmin = plotter.axis_minimum;
          plotter.legend_offset_x = 0.0825;
          plotter.legend_border_size = 0;
          plotter.legend_text_size = 0.03;

          std::string tmp; 
          std::string y_title="Events#times10^{3}";
          double xmax = -1111; //Default value
          std::string xaxislabel = "";
          if (prefix=="pTmu")
          {
            tmp = "Muon p_{T}";
            xaxislabel = "p_{T, #mu} [GeV/c]";
            y_title="d#sigma/dp_{t} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            xmax = 13;
            tempsimxsec->GetXaxis()->SetRange(1, xmax);
            tempxsec->GetXaxis()->SetRange(1, xmax);
          }
          else if (prefix=="pZmu")
          {
            tmp = "Muon p_{Z}";
            xaxislabel = "p_{Z, #mu} [GeV/c]";
            y_title="d#sigma/dp_{Z} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            //xmax = 13;
            //tempsimxsec->GetXaxis()->SetRange(1, xmax);
            //tempxsec->GetXaxis()->SetRange(1, xmax);
          }
          else if (prefix=="BjorkenX")
          {
            tmp = "Bjorken X";
            xaxislabel = "Bjorken X";
            y_title="d#sigma/dX (x10^{-39}) (cm^{2}/Nucleon)";
            xmax = 6;
            tempsimxsec->GetXaxis()->SetRange(1, xmax);
            can.SetLogx();
          }
          else if (prefix=="Erecoil")
          {
            tmp = "E_{recoil}";
            y_title="d#sigma/dE_{recoil} (x10^{-39}) (cm^{2}/GeV/Nucleon)";
            xaxislabel = "E_{recoil} [GeV]";
          } 
          else if (prefix=="Emu")
          {
            tmp = "E_{#mu}";
            y_title="d#sigma/dE_{#mu} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            xaxislabel = "E_{#mu} [GeV]";
            xmax = 11;
            tempsimxsec->GetXaxis()->SetRange(1, xmax);
            tempxsec->GetXaxis()->SetRange(1, xmax);
          } 
          std::cout<<"1.114\n";
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;
          tempxsec->SetTitle("Data");
          tempxsec->SetLineWidth(3);
          tempxsec->SetMarkerStyle(20);
          tempsimxsec->SetLineColor(46);
          tempsimxsec->SetLineWidth(2);
          tempsimxsec->SetTitle("MnvTune 4.3.0 (CV)");
          //double max_y = tempsimxsec->GetMaximum();
          //tempsimxsec->SetMaximum(max_y * 1.2);
          //max_y = tempxsec->GetMaximum();
          //tempxsec->SetMaximum(max_y * 1.2);
          //plotter.axis_maximum = max_y * 1.2;
          std::cout<<"1.115\n";

          hist_G18_02a. SetLineColor(TColor::GetColor("#332288"));
          hist_G18_02b. SetLineColor(TColor::GetColor("#44AA99"));
          hist_G18_10a. SetLineColor(TColor::GetColor("#88CCEE"));
          hist_G18_10b. SetLineColor(TColor::GetColor("#DDCC77"));
          hist_NEUT_LFG.SetLineColor(TColor::GetColor("#AA4499"));

          hist_G18_02a.SetLineWidth(3);
          hist_G18_02b.SetLineWidth(3);
          hist_G18_10a.SetLineWidth(3);
          hist_G18_10b.SetLineWidth(3);
          hist_NEUT_LFG.SetLineWidth(3);
          std::cout<<"1.116\n";

          tempxsec->Scale(1e39);
          tempsimxsec->Scale(1e39);
          hist_G18_02a.Scale(1e39);
          hist_G18_02b.Scale(1e39);
          hist_G18_10a.Scale(1e39);
          hist_G18_10b.Scale(1e39);
          hist_NEUT_LFG.Scale(1e39);


          arr->Add(tempsimxsec);

          hist_G18_02a.SetTitle("GENIE3 G18_02a");
          hist_G18_02b.SetTitle("GENIE3 G18_02b");
          hist_G18_10a.SetTitle("GENIE3 G18_10a");
          hist_G18_10b.SetTitle("GENIE3 G18_10b");
          hist_NEUT_LFG.SetTitle("NEUT LFG ma105");

          std::cout<<"1.117\n";
          //hist_G18_02a.SetMaximum(max_y);
          //hist_G18_02b.SetMaximum(max_y);
          //hist_G18_10a.SetMaximum(max_y);
          //hist_G18_10b.SetMaximum(max_y);
          //hist_NEUT_LFG.SetMaximum(max_y);
          //hist_NEUT_SFG.SetMaximum(max_y);

          arr->Add(&hist_G18_02a);
          arr->Add(&hist_G18_02b);
          arr->Add(&hist_G18_10a);
          arr->Add(&hist_G18_10b);
          arr->Add(&hist_NEUT_LFG);

          std::cout<<"1.118\n";
          //Used in ratio drawing later
          PlotUtils::MnvH1D * cv_hist_ratio= (PlotUtils::MnvH1D*)tempxsec->Clone("MnvTune 4.3.0 (CV)");
          PlotUtils::MnvH1D * hist_G18_02a_ratio= (PlotUtils::MnvH1D*)tempxsec->Clone("GENIE3 G18_02a");
          PlotUtils::MnvH1D * hist_G18_02b_ratio= (PlotUtils::MnvH1D*)tempxsec->Clone("GENIE3 G18_02b");
          PlotUtils::MnvH1D * hist_G18_10a_ratio= (PlotUtils::MnvH1D*)tempxsec->Clone("GENIE3 G18_10a");
          PlotUtils::MnvH1D * hist_G18_10b_ratio= (PlotUtils::MnvH1D*)tempxsec->Clone("GENIE3 G18_10b");
          PlotUtils::MnvH1D * hist_NEUT_LFG_ratio= (PlotUtils::MnvH1D*)tempxsec->Clone("NEUT LFG ma105");
          std::cout<<"1.119\n";
          cv_hist_ratio->Divide(cv_hist_ratio, tempsimxsec);
          hist_G18_02a_ratio->Divide(hist_G18_02a_ratio, &hist_G18_02a);
          hist_G18_02b_ratio->Divide(hist_G18_02b_ratio, &hist_G18_02b);
          hist_G18_10a_ratio->Divide(hist_G18_10a_ratio, &hist_G18_10a);
          hist_G18_10b_ratio->Divide(hist_G18_10b_ratio, &hist_G18_10b);
          hist_NEUT_LFG_ratio->Divide(hist_NEUT_LFG_ratio, &hist_NEUT_LFG);
          std::cout<<"1.1110\n";
          cv_hist_ratio->SetLineColor(46);
          hist_G18_02a_ratio->SetLineColor(TColor::GetColor("#332288")); 
          hist_G18_02b_ratio->SetLineColor(TColor::GetColor("#44AA99")); 
          hist_G18_10a_ratio->SetLineColor(TColor::GetColor("#88CCEE")); 
          hist_G18_10b_ratio->SetLineColor(TColor::GetColor("#DDCC77")); 
          hist_NEUT_LFG_ratio->SetLineColor(TColor::GetColor("#AA4499")); 
          std::cout<<"1.1111\n";
          cv_hist_ratio->SetMarkerColor(46);
          hist_G18_02a_ratio->SetMarkerColor(TColor::GetColor("#332288")); 
          hist_G18_02b_ratio->SetMarkerColor(TColor::GetColor("#44AA99")); 
          hist_G18_10a_ratio->SetMarkerColor(TColor::GetColor("#88CCEE")); 
          hist_G18_10b_ratio->SetMarkerColor(TColor::GetColor("#DDCC77")); 
          hist_NEUT_LFG_ratio->SetMarkerColor(TColor::GetColor("#AA4499")); 
          std::cout<<"1.1112\n";
          cv_hist_ratio->SetMinimum(0.0001);
          cv_hist_ratio->SetMaximum(1.9999);
          std::cout<<"1.1113\n";
          plotter.ApplyNextLineStyle(cv_hist_ratio, true, false);
          plotter.ApplyNextLineStyle(cv_hist_ratio, false, false);
          plotter.ApplyNextLineStyle(hist_G18_02a_ratio, false, false);
          plotter.ApplyNextLineStyle(hist_G18_02b_ratio, false, false);
          plotter.ApplyNextLineStyle(hist_G18_10a_ratio, false, false);
          plotter.ApplyNextLineStyle(hist_G18_10b_ratio, false, false);
          plotter.ApplyNextLineStyle(hist_NEUT_LFG_ratio, false, false);
          std::cout<<"1.1114\n";

          TH1D* NEUT_tune_SFG_maqe1_03_xsec=(TH1D*)util::GetIngredientSmart<TH1D>(*NEUT_tune_SFG_maqe1_03_file, prefix+"_rw").get()->Clone("NEUT SF ma=103");
          NEUT_tune_SFG_maqe1_03_xsec->Scale(1., "width");
          PlotUtils::MnvH1D hist_NEUT_SFG = PlotUtils::MnvH1D(*NEUT_tune_SFG_maqe1_03_xsec);
          hist_NEUT_SFG.SetLineColor(15);
          hist_NEUT_SFG.SetLineWidth(3);
          hist_NEUT_SFG.Scale(1e39);
          hist_NEUT_SFG.SetTitle("NEUT SFG ma103");
          PlotUtils::MnvH1D * hist_NEUT_SFG_ratio= (PlotUtils::MnvH1D*)tempxsec->Clone("NEUT SFG ma103");
          hist_NEUT_SFG_ratio->Divide(hist_NEUT_SFG_ratio, &hist_NEUT_SFG);
          hist_NEUT_SFG_ratio->SetLineColor(15);
          hist_NEUT_SFG_ratio->SetMarkerColor(15);
          plotter.ApplyNextLineStyle(hist_NEUT_SFG_ratio, false, false);
          arr->Add(&hist_NEUT_SFG);
          std::cout<<"1.1115\n";

          

          std::string title = tmp + " Tracker Model Comparisons";
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;

          // Adjust pad heights manually (top bigger, bottom smaller)
          TPad *pad1 = new TPad("pad1","pad1",0,0.35,1,1.0); // top pad
          TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.35); // bottom pad
          pad1->SetBottomMargin(0); // remove x-axis label space for top
          pad2->SetTopMargin(0);
          pad2->SetBottomMargin(0.25);
          pad1->Draw();
          pad2->Draw();
          std::cout<<"1.1116\n";
          if (prefix=="BjorkenX")
          {
            pad1->SetLogx();
            pad2->SetLogx();
          }

          std::cout<<"1.1117\n";
        
          std::stringstream POTstr;
          POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOT;
          std::string potstr = std::string("POT Used: ")+POTstr.str();
          plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.93, 0.025, 1, 42, 13, 0);
          plotter.WritePreliminary( 0.5, 0.15, 0.035, true);
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 25, 1, 53, 22, 0);
          plotter.AddPlotLabel(y_title.c_str() , 0.08, 0.9, 0.03, 1, 42, 33, 90);
          plotter.AddPlotLabel("Data/MC" , 0.08, 0.3, 0.03, 1, 42, 33, 90);
          plotter.AddPlotLabel(xaxislabel.c_str() , 0.9, 0.08, 0.03, 1, 42, 33, 0);

          //Interaction Type
          pad1->cd(); 
          plotter.axis_minimum=0.001;
          std::cout<<"1.1118\n";
          //tempxsec->Draw("HIST");
          //tempsimxsec->Draw("HIST SAME");
          //plotter.DrawDataMC(tempxsec, tempsimxsec);
          tempxsec->GetXaxis()->SetLabelFont(43);
          tempxsec->GetXaxis()->SetLabelSize(25);
          std::string legPos = (prefix == "BjorkenX") ? "TL" : "TR";
          plotter.DrawDataMCVariations(tempxsec, arr, 1, legPos, true, true, false, false, true, false);
          pad2->cd();
          plotter.axis_minimum=defaultmin;
          std::cout<<"1.1119\n";
          cv_hist_ratio->GetXaxis()->SetLabelFont(43);
          cv_hist_ratio->GetXaxis()->SetLabelSize(25);
          cv_hist_ratio->Draw("HIST");
          hist_G18_02a_ratio->Draw("HIST SAME");
          hist_G18_02b_ratio->Draw("HIST SAME");
          hist_G18_10a_ratio->Draw("HIST SAME");
          hist_G18_10b_ratio->Draw("HIST SAME");
          hist_NEUT_LFG_ratio->Draw("HIST SAME");
          if (hist_NEUT_SFG_ratio) hist_NEUT_SFG_ratio->Draw("HIST SAME");



          const TAxis *axis = cv_hist_ratio->GetXaxis();

          TLine *line = new TLine(axis->GetBinLowEdge( axis->GetFirst() ), 1.0, axis->GetBinUpEdge(  axis->GetLast() ), 1.0);
          line->SetLineColor(kBlack);
          line->SetLineStyle(2); // 2 = Dashed
          line->SetLineWidth(3);
          line->SetLineColor(36);
          line->Draw();
          std::cout<<"1.1120\n";
          //plotter.DrawDataMCRatio(tempxsec, tempsimxsec, 1, false, true, 0.5, 1.5, "");
          //plotter.DrawDataMCRatio(tempxsec, &hist_G18_02a, 1, false, true, 0.5, 1.5, "");
          //can.SetLogx();
          pad1->cd();
          can.cd();
          can.Print(("Tracker_"+prefix+"_ModelComparisonsCrossSection.png").c_str());
          std::cout<<"1.1121\n";
        }


        if (pdg == -14 && prefix=="pTmu")
        { //Plotting cross-sections
          TCanvas can("Temp", "Temp", 1200, 1000);
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kDefaultStyle);
          plotter.draw_normalized_to_bin_width=false; //Already width normlized by normalise()
          plotter.axis_label_font=43;
          plotter.axis_label_size = 25;
          plotter.headroom = 1.2;
          /* plotter.axis_title_font_y=43;
          plotter.axis_title_size_y = 35;
          plotter.axis_title_font_x=43;
          plotter.axis_title_size_x = 35;
          */
          plotter.ApplyAxisStyle();
          //plotter.axis_maximum = 0.4;


          TObjArray* arr = new TObjArray();
          auto temp2p2h = simulatedCrossSection2P2H->Clone();
          auto tempdis = simulatedCrossSectionDIS->Clone();
          auto tempres = simulatedCrossSectionRES->Clone();
          auto tempqe = simulatedCrossSectionQE->Clone();
          auto tempother = simulatedCrossSectionOther->Clone();

          auto tempsimxsec = temp2p2h->Clone();
          tempsimxsec->Add(tempdis);
          tempsimxsec->Add(tempres);
          tempsimxsec->Add(tempqe);
          tempsimxsec->Add(tempother);

          auto tempxsec = crossSection->Clone();
          temp2p2h->Scale(1e39);
          tempdis->Scale(1e39);
          tempres->Scale(1e39);
          tempqe->Scale(1e39);
          tempother->Scale(1e39);
          tempxsec->Scale(1e39);
          tempsimxsec->Scale(1e39);

          arr->Add(temp2p2h);
          arr->Add(tempdis);
          arr->Add(tempres);
          arr->Add(tempqe);
          arr->Add(tempother);

          double defaultmin = plotter.axis_minimum;
          plotter.legend_offset_x = -0.0825;
          plotter.legend_border_size = 0;
          plotter.legend_text_size = 0.03;

          plotter.mc_line_width=0;  //Why??

          std::string tmp; 
          std::string y_title="Events#times10^{3}";
          double xmax = -1111; //Default value
          std::string xaxislabel = "";
          if (prefix=="pTmu")
          {
            tmp = "Muon p_{T}";
            xaxislabel = "p_{T, #mu} [GeV/c]";
            y_title="d#sigma/dp_{t} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            xmax = 13;
            tempsimxsec->GetXaxis()->SetRange(1, xmax);
            tempxsec->GetXaxis()->SetRange(1, xmax);
          }
          else if (prefix=="pZmu")
          {
            tmp = "Muon p_{Z}";
            xaxislabel = "p_{Z, #mu} [GeV/c]";
            y_title="d#sigma/dp_{Z} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            //xmax = 13;
            //tempsimxsec->GetXaxis()->SetRange(1, xmax);
            //tempxsec->GetXaxis()->SetRange(1, xmax);
          }
          else if (prefix=="BjorkenX")
          {
            tmp = "Bjorken X";
            xaxislabel = "Bjorken X";
            y_title="d#sigma/dX (x10^{-39}) (cm^{2}/Nucleon)";
            //xmax = 6;
            //tempsimxsec->GetXaxis()->SetRange(1, xmax);
            can.SetLogx();
          }
          else if (prefix=="Erecoil")
          {
            tmp = "E_{recoil}";
            y_title="d#sigma/dE_{recoil} (x10^{-39}) (cm^{2}/GeV/Nucleon)";
            xaxislabel = "E_{recoil} [GeV]";
          } 
          else if (prefix=="Emu")
          {
            tmp = "E_{#mu}";
            y_title="d#sigma/dE_{#mu} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            xaxislabel = "E_{#mu} [GeV]";
            xmax = 11;
            tempsimxsec->GetXaxis()->SetRange(1, xmax);
            tempxsec->GetXaxis()->SetRange(1, xmax);
          } 

          std::string title = nutypestr+ " " +tmp;
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;

          plotter.DrawErrorSummary(tempsimxsec, "TR", true, true,  0.00001, false, "", true, "", true);

          can.Print((prefix+"_"+nutypestr+"_XSecErrors.png").c_str());
          can.Clear();


          // Adjust pad heights manually (top bigger, bottom smaller)
          TPad *pad1 = new TPad("pad1","pad1",0,0.45,1,1.0); // top pad
          TPad *pad2 = new TPad("pad2","pad2",0,0.25,1,0.45); // middle pad
          TPad *pad3 = new TPad("pad3","pad3",0,0.05,1,0.25); // bottom pad
          pad1->SetBottomMargin(0); // remove x-axis label space for top
          pad2->SetTopMargin(0);
          pad2->SetBottomMargin(0.0);
          pad3->SetTopMargin(0);
          pad3->SetBottomMargin(0.25);
          pad1->Draw();
          pad2->Draw();
          pad3->Draw();

          if (prefix=="BjorkenX")
          {
            pad1->SetLogx();
            pad2->SetLogx();
            pad3->SetLogx();
          }


          std::stringstream POTstr;
          POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOT;
          std::string potstr = std::string("POT Used: ")+POTstr.str();
          plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.93, 0.025, 1, 42, 13, 0);
          
          plotter.WritePreliminary( 0.5, 0.30, 0.035, true);
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 25, 1, 53, 22, 0);
          plotter.AddPlotLabel(y_title.c_str() , 0.08, 0.9, 0.03, 1, 42, 33, 90);
          plotter.AddPlotLabel("Data/MC" , 0.08, 0.4, 0.03, 1, 42, 33, 90);
          plotter.AddPlotLabel("Mine/Anezkas" , 0.08, 0.25, 0.02, 1, 42, 33, 90);
          plotter.AddPlotLabel(xaxislabel.c_str() , 0.9, 0.08, 0.03, 1, 42, 33, 0);

          //Interaction Type
          pad1->cd(); 
          plotter.axis_minimum=0.001;

          int channelColArr[5] = {TColor::GetColor("#d2c271"), TColor::GetColor("#CA5454"), TColor::GetColor("#40B0A6"), TColor::GetColor("#9285d5"), TColor::GetColor("#c2c2c2")};
          plotter.DrawDataStackedMC(tempxsec, arr, channelColArr, 1, "TR", "Data", 1001, "", "", false, 1, xmax);
          
          Double_t AnezkaBins[15] = {0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1, 1.25, 1.5, 2.5, 4.5}; 
          TH1D *AnezkasResults = new TH1D("AnezkasResults","pTmu1D",14, AnezkaBins);
          TH1D *AnezkasResultsWithErr = new TH1D("AnezkasResultsWithErr","pTmu1D",14, AnezkaBins);
          TH1D *AnezkasMCOther = new TH1D("AnezkasMCOther","pTmu1D",14, AnezkaBins);
          TH1D *AnezkasMC2P2H = new TH1D("AnezkasMC2P2H","pTmu1D",14, AnezkaBins);
          TH1D *AnezkasMCQE = new TH1D("AnezkasMCQE","pTmu1D",14, AnezkaBins);
          TH1D *AnezkasMCDIS = new TH1D("AnezkasMCDIS","pTmu1D",14, AnezkaBins);
          TH1D *AnezkasMCSoftDIS = new TH1D("AnezkasMCSoftDIS","pTmu1D",14, AnezkaBins);
          TH1D *AnezkasMCRES = new TH1D("AnezkasMCRES","pTmu1D",14, AnezkaBins);
          AnezkasMCOther->SetBinContent(15,0);
          AnezkasMC2P2H->SetBinContent(15,0);
          AnezkasMCQE->SetBinContent(15,0);
          AnezkasMCDIS->SetBinContent(15,0);
          AnezkasMCSoftDIS->SetBinContent(15,0);
          AnezkasMCRES->SetBinContent(15,0);
          AnezkasMCOther->SetBinError(15,0);
          AnezkasMC2P2H->SetBinError(15,0);
          AnezkasMCQE->SetBinError(15,0);
          AnezkasMCDIS->SetBinError(15,0);
          AnezkasMCSoftDIS->SetBinError(15,0);
          AnezkasMCRES->SetBinError(15,0);
          AnezkasResults->SetBinContent(15,0);
          AnezkasResults->SetBinError(15,0);
          AnezkasResultsWithErr->SetBinContent(15,0);
          AnezkasResultsWithErr->SetBinError(15,0);
          AnezkasResults->SetBinContent(1,1.157488);
          AnezkasResults->SetBinContent(2,4.114783);
          AnezkasResults->SetBinContent(3,8.795731);
          AnezkasResults->SetBinContent(4,13.88215);
          AnezkasResults->SetBinContent(5,18.00916);
          AnezkasResults->SetBinContent(6,21.14226);
          AnezkasResults->SetBinContent(7,22.71399);
          AnezkasResults->SetBinContent(8,22.08141);
          AnezkasResults->SetBinContent(9,17.32603);
          AnezkasResults->SetBinContent(10,11.46697);
          AnezkasResults->SetBinContent(11,5.350887);
          AnezkasResults->SetBinContent(12,1.606471);
          AnezkasResults->SetBinContent(13,0.2060639);
          AnezkasResults->SetBinContent(14,0.007842739);
          AnezkasResults->SetBinError(1,0.01399283);
          AnezkasResults->SetBinError(2,0.02535665);
          AnezkasResults->SetBinError(3,0.03258475);
          AnezkasResults->SetBinError(4,0.04578264);
          AnezkasResults->SetBinError(5,0.0513595);
          AnezkasResults->SetBinError(6,0.05548039);
          AnezkasResults->SetBinError(7,0.05829278);
          AnezkasResults->SetBinError(8,0.04659434);
          AnezkasResults->SetBinError(9,0.04257159);
          AnezkasResults->SetBinError(10,0.0348274);
          AnezkasResults->SetBinError(11,0.02025414);
          AnezkasResults->SetBinError(12,0.01099306);
          AnezkasResults->SetBinError(13,0.002241809);
          AnezkasResults->SetBinError(14,0.0004377194);
          AnezkasResults->SetEntries(1160151);
          AnezkasResults->SetStats(0);
          AnezkasResults->SetMarkerStyle(20);
          AnezkasResults->SetMarkerStyle(20);
          AnezkasResults->SetLineColor(kGreen);
          AnezkasResults->SetMarkerColor(kGreen);
          AnezkasResults->Draw("SAME HIST P E1 X0");

          AnezkasResultsWithErr->SetBinContent(1,1.157488);
          AnezkasResultsWithErr->SetBinContent(2,4.114783);
          AnezkasResultsWithErr->SetBinContent(3,8.795731);
          AnezkasResultsWithErr->SetBinContent(4,13.88215);
          AnezkasResultsWithErr->SetBinContent(5,18.00916);
          AnezkasResultsWithErr->SetBinContent(6,21.14226);
          AnezkasResultsWithErr->SetBinContent(7,22.71399);
          AnezkasResultsWithErr->SetBinContent(8,22.08141);
          AnezkasResultsWithErr->SetBinContent(9,17.32603);
          AnezkasResultsWithErr->SetBinContent(10,11.46697);
          AnezkasResultsWithErr->SetBinContent(11,5.350887);
          AnezkasResultsWithErr->SetBinContent(12,1.606471);
          AnezkasResultsWithErr->SetBinContent(13,0.2060639);
          AnezkasResultsWithErr->SetBinContent(14,0.007842739);
          AnezkasResultsWithErr->SetBinError(1,0.07010461);
          AnezkasResultsWithErr->SetBinError(2,0.2454011);
          AnezkasResultsWithErr->SetBinError(3,0.5161058);
          AnezkasResultsWithErr->SetBinError(4,0.808047);
          AnezkasResultsWithErr->SetBinError(5,1.026125);
          AnezkasResultsWithErr->SetBinError(6,1.183164);
          AnezkasResultsWithErr->SetBinError(7,1.25153);
          AnezkasResultsWithErr->SetBinError(8,1.206148);
          AnezkasResultsWithErr->SetBinError(9,0.9793414);
          AnezkasResultsWithErr->SetBinError(10,0.7114732);
          AnezkasResultsWithErr->SetBinError(11,0.3832299);
          AnezkasResultsWithErr->SetBinError(12,0.1312486);
          AnezkasResultsWithErr->SetBinError(13,0.01706393);
          AnezkasResultsWithErr->SetBinError(14,0.0004377194);
          AnezkasResultsWithErr->SetEntries(1160151);
          AnezkasResultsWithErr->SetStats(0);
          AnezkasResultsWithErr->SetMarkerStyle(20);
          AnezkasResultsWithErr->SetMarkerStyle(20);
          AnezkasResultsWithErr->SetLineColor(kGreen);
          AnezkasResultsWithErr->SetMarkerColor(kGreen);
          AnezkasResultsWithErr->Draw("HIST E1 P SAME X0");

          AnezkasMCOther->SetBinContent(1,0.3070969);
          AnezkasMCOther->SetBinContent(2,0.8359186);
          AnezkasMCOther->SetBinContent(3,1.120626);
          AnezkasMCOther->SetBinContent(4,1.059088);
          AnezkasMCOther->SetBinContent(5,0.851739);
          AnezkasMCOther->SetBinContent(6,0.6256);
          AnezkasMCOther->SetBinContent(7,0.423567);
          AnezkasMCOther->SetBinContent(8,0.2045737);
          AnezkasMCOther->SetBinContent(9,0.05631914);
          AnezkasMCOther->SetBinContent(10,0.01300191);
          AnezkasMCOther->SetBinContent(11,0.00193429);
          AnezkasMCOther->SetBinContent(12,0.0002543287);
          AnezkasMCOther->SetBinContent(13,4.504711e-06);
          AnezkasMCOther->SetBinError(1,0.003015121);
          AnezkasMCOther->SetBinError(2,0.004899222);
          AnezkasMCOther->SetBinError(3,0.004767944);
          AnezkasMCOther->SetBinError(4,0.005240861);
          AnezkasMCOther->SetBinError(5,0.004662188);
          AnezkasMCOther->SetBinError(6,0.004000697);
          AnezkasMCOther->SetBinError(7,0.003301374);
          AnezkasMCOther->SetBinError(8,0.001644984);
          AnezkasMCOther->SetBinError(9,0.0008844199);
          AnezkasMCOther->SetBinError(10,0.0004419895);
          AnezkasMCOther->SetBinError(11,0.0001402958);
          AnezkasMCOther->SetBinError(12,5.155948e-05);
          AnezkasMCOther->SetBinError(13,3.193204e-06);
          AnezkasMCOther->SetEntries(222994.7);

          AnezkasMC2P2H->SetBinContent(1,0.1808324);
          AnezkasMC2P2H->SetBinContent(2,0.6043056);
          AnezkasMC2P2H->SetBinContent(3,1.254858);
          AnezkasMC2P2H->SetBinContent(4,2.131201);
          AnezkasMC2P2H->SetBinContent(5,2.893184);
          AnezkasMC2P2H->SetBinContent(6,3.133849);
          AnezkasMC2P2H->SetBinContent(7,2.578227);
          AnezkasMC2P2H->SetBinContent(8,1.475026);
          AnezkasMC2P2H->SetBinContent(9,0.6202112);
          AnezkasMC2P2H->SetBinContent(10,0.1718075);
          AnezkasMC2P2H->SetBinContent(11,0.01195846);
          AnezkasMC2P2H->SetBinError(1,0.001747763);
          AnezkasMC2P2H->SetBinError(2,0.003291973);
          AnezkasMC2P2H->SetBinError(3,0.004855);
          AnezkasMC2P2H->SetBinError(4,0.009870802);
          AnezkasMC2P2H->SetBinError(5,0.01467615);
          AnezkasMC2P2H->SetBinError(6,0.01623659);
          AnezkasMC2P2H->SetBinError(7,0.01289437);
          AnezkasMC2P2H->SetBinError(8,0.004667596);
          AnezkasMC2P2H->SetBinError(9,0.002290987);
          AnezkasMC2P2H->SetBinError(10,0.001188522);
          AnezkasMC2P2H->SetBinError(11,0.0002460717);
          AnezkasMC2P2H->SetEntries(280290.8);

          AnezkasMCQE->SetBinContent(1,0.2263757);
          AnezkasMCQE->SetBinContent(2,0.8451885);
          AnezkasMCQE->SetBinContent(3,1.811677);
          AnezkasMCQE->SetBinContent(4,2.820535);
          AnezkasMCQE->SetBinContent(5,3.600091);
          AnezkasMCQE->SetBinContent(6,4.035558);
          AnezkasMCQE->SetBinContent(7,4.19025);
          AnezkasMCQE->SetBinContent(8,3.891123);
          AnezkasMCQE->SetBinContent(9,2.912163);
          AnezkasMCQE->SetBinContent(10,1.80704);
          AnezkasMCQE->SetBinContent(11,0.783444);
          AnezkasMCQE->SetBinContent(12,0.1942666);
          AnezkasMCQE->SetBinContent(13,0.009802562);
          AnezkasMCQE->SetBinContent(14,2.941885e-05);
          AnezkasMCQE->SetBinError(1,0.001663243);
          AnezkasMCQE->SetBinError(2,0.003127107);
          AnezkasMCQE->SetBinError(3,0.004010935);
          AnezkasMCQE->SetBinError(4,0.005939925);
          AnezkasMCQE->SetBinError(5,0.006800773);
          AnezkasMCQE->SetBinError(6,0.007316766);
          AnezkasMCQE->SetBinError(7,0.007685678);
          AnezkasMCQE->SetBinError(8,0.005476894);
          AnezkasMCQE->SetBinError(9,0.00495443);
          AnezkasMCQE->SetBinError(10,0.003972278);
          AnezkasMCQE->SetBinError(11,0.002019631);
          AnezkasMCQE->SetBinError(12,0.0009964447);
          AnezkasMCQE->SetBinError(13,0.0001133221);
          AnezkasMCQE->SetBinError(14,6.612182e-06);
          AnezkasMCQE->SetEntries(2468003);

          AnezkasMCDIS->SetBinContent(3,3.121936e-05);
          AnezkasMCDIS->SetBinContent(4,0.0009299204);
          AnezkasMCDIS->SetBinContent(5,0.005945889);
          AnezkasMCDIS->SetBinContent(6,0.02056209);
          AnezkasMCDIS->SetBinContent(7,0.08148147);
          AnezkasMCDIS->SetBinContent(8,0.545608);
          AnezkasMCDIS->SetBinContent(9,2.024473);
          AnezkasMCDIS->SetBinContent(10,2.326776);
          AnezkasMCDIS->SetBinContent(11,1.23369);
          AnezkasMCDIS->SetBinContent(12,0.517761);
          AnezkasMCDIS->SetBinContent(13,0.1186402);
          AnezkasMCDIS->SetBinContent(14,0.005384884);
          AnezkasMCDIS->SetBinError(3,1.832147e-05);
          AnezkasMCDIS->SetBinError(4,0.0001345316);
          AnezkasMCDIS->SetBinError(5,0.0003335521);
          AnezkasMCDIS->SetBinError(6,0.0006244211);
          AnezkasMCDIS->SetBinError(7,0.001211941);
          AnezkasMCDIS->SetBinError(8,0.002129607);
          AnezkasMCDIS->SetBinError(9,0.004054516);
          AnezkasMCDIS->SetBinError(10,0.004396799);
          AnezkasMCDIS->SetBinError(11,0.002532571);
          AnezkasMCDIS->SetBinError(12,0.0016907);
          AnezkasMCDIS->SetBinError(13,0.0004086432);
          AnezkasMCDIS->SetBinError(14,8.473396e-05);
          AnezkasMCDIS->SetEntries(913868.2);

          AnezkasMCSoftDIS->SetBinContent(1,0.0329431);
          AnezkasMCSoftDIS->SetBinContent(2,0.2457491);
          AnezkasMCSoftDIS->SetBinContent(3,1.050188);
          AnezkasMCSoftDIS->SetBinContent(4,2.494229);
          AnezkasMCSoftDIS->SetBinContent(5,4.112991);
          AnezkasMCSoftDIS->SetBinContent(6,5.738131);
          AnezkasMCSoftDIS->SetBinContent(7,7.043573);
          AnezkasMCSoftDIS->SetBinContent(8,7.433257);
          AnezkasMCSoftDIS->SetBinContent(9,5.000544);
          AnezkasMCSoftDIS->SetBinContent(10,2.63645);
          AnezkasMCSoftDIS->SetBinContent(11,1.257635);
          AnezkasMCSoftDIS->SetBinContent(12,0.3467069);
          AnezkasMCSoftDIS->SetBinContent(13,0.0282341);
          AnezkasMCSoftDIS->SetBinContent(14,0.0005485851);
          AnezkasMCSoftDIS->SetBinError(1,0.0007054644);
          AnezkasMCSoftDIS->SetBinError(2,0.001913068);
          AnezkasMCSoftDIS->SetBinError(3,0.003411352);
          AnezkasMCSoftDIS->SetBinError(4,0.006068694);
          AnezkasMCSoftDIS->SetBinError(5,0.007805679);
          AnezkasMCSoftDIS->SetBinError(6,0.009237756);
          AnezkasMCSoftDIS->SetBinError(7,0.01025576);
          AnezkasMCSoftDIS->SetBinError(8,0.007457671);
          AnezkasMCSoftDIS->SetBinError(9,0.006080943);
          AnezkasMCSoftDIS->SetBinError(10,0.004359117);
          AnezkasMCSoftDIS->SetBinError(11,0.002341241);
          AnezkasMCSoftDIS->SetBinError(12,0.001253425);
          AnezkasMCSoftDIS->SetBinError(13,0.0001866556);
          AnezkasMCSoftDIS->SetBinError(14,2.71622e-05);
          AnezkasMCSoftDIS->SetEntries(3312344);

          AnezkasMCRES->SetBinContent(1,0.6281972);
          AnezkasMCRES->SetBinContent(2,1.909452);
          AnezkasMCRES->SetBinContent(3,3.758868);
          AnezkasMCRES->SetBinContent(4,5.486567);
          AnezkasMCRES->SetBinContent(5,6.568613);
          AnezkasMCRES->SetBinContent(6,7.132694);
          AnezkasMCRES->SetBinContent(7,7.1922);
          AnezkasMCRES->SetBinContent(8,6.332155);
          AnezkasMCRES->SetBinContent(9,4.463852);
          AnezkasMCRES->SetBinContent(10,2.704325);
          AnezkasMCRES->SetBinContent(11,1.179173);
          AnezkasMCRES->SetBinContent(12,0.2981872);
          AnezkasMCRES->SetBinContent(13,0.01888634);
          AnezkasMCRES->SetBinContent(14,0.0001943831);
          AnezkasMCRES->SetBinError(1,0.00317438);
          AnezkasMCRES->SetBinError(2,0.005514471);
          AnezkasMCRES->SetBinError(3,0.006819899);
          AnezkasMCRES->SetBinError(4,0.009705431);
          AnezkasMCRES->SetBinError(5,0.01073782);
          AnezkasMCRES->SetBinError(6,0.0112432);
          AnezkasMCRES->SetBinError(7,0.01131868);
          AnezkasMCRES->SetBinError(8,0.007524324);
          AnezkasMCRES->SetBinError(9,0.006332856);
          AnezkasMCRES->SetBinError(10,0.004951286);
          AnezkasMCRES->SetBinError(11,0.00256773);
          AnezkasMCRES->SetBinError(12,0.001332316);
          AnezkasMCRES->SetBinError(13,0.0001757611);
          AnezkasMCRES->SetBinError(14,2.018734e-05);
          AnezkasMCRES->SetEntries(3339058);

          TLegend* leg = nullptr;
          TList* primitives = gPad->GetListOfPrimitives();
          for (int i = 0; i < primitives->GetSize(); ++i) {
              TObject* obj = primitives->At(i);
              if (obj->InheritsFrom(TLegend::Class())) {
                  leg = static_cast<TLegend*>(obj);
                  break; 
              }
          }
          if (leg) {
              leg->AddEntry(AnezkasResults, "Anezka's Paper Results", "ple"); 
              leg->Draw(); 
              gPad->Update();
          }
        
          pad2->cd();
          plotter.axis_minimum=defaultmin;
          plotter.DrawDataMCRatio(tempxsec, tempsimxsec, 1, true, true, 0.5, 1.5, "");
          TH1D simTH1D =tempsimxsec->GetCVHistoWithStatError();
          TH1D* CloneAnezka = (TH1D*)AnezkasResultsWithErr->Clone();
          CloneAnezka->Divide(CloneAnezka, &simTH1D);
          CloneAnezka->Draw("SAME X0");
          pad3->cd();
          AnezkasResultsWithErr->GetXaxis()->SetRange(1, 13);
          tempxsec->GetXaxis()->SetRange(1, 13);
          plotter.DrawDataMCRatio(tempxsec, AnezkasResultsWithErr, 1, true, 0.75, 1.25, "");
          TH1D *AnezkasMCSum = (TH1D*)AnezkasMCOther->Clone("MCSum");
          AnezkasMCSum->Add(AnezkasMC2P2H);
          AnezkasMCSum->Add(AnezkasMCQE);
          AnezkasMCSum->Add(AnezkasMCDIS);
          AnezkasMCSum->Add(AnezkasMCSoftDIS);
          AnezkasMCSum->Add(AnezkasMCRES);
          TH1D *tempth1simxsec = (TH1D*)tempsimxsec;
          tempth1simxsec->Divide(tempth1simxsec, AnezkasMCSum);
          tempth1simxsec->SetLineColor(kBlack);
          tempth1simxsec->SetLineWidth(2);
          tempth1simxsec->Draw("HIST SAME");
          TH1D *AnezkasMCDISSum = (TH1D*)AnezkasMCDIS->Clone("MCDIS");
          AnezkasMCDISSum->Add(AnezkasMCSoftDIS);
          TH1D *tempTH12p2h = (TH1D*)temp2p2h;
          TH1D *tempTH1dis = (TH1D*)tempdis;
          TH1D *tempTH1res = (TH1D*)tempres;
          TH1D *tempTH1qe = (TH1D*)tempqe;
          TH1D *tempTH1other = (TH1D*)tempother;
          tempTH12p2h->Divide(tempTH12p2h, AnezkasMC2P2H);
          tempTH1dis->Divide(tempTH1dis, AnezkasMCDISSum);
          tempTH1res->Divide(tempTH1res, AnezkasMCRES);
          tempTH1qe->Divide(tempTH1qe, AnezkasMCQE);
          tempTH1other->Divide(tempTH1other, AnezkasMCOther);
          tempTH12p2h->SetLineWidth(2);
          tempTH1dis->SetLineWidth(2);
          tempTH1res->SetLineWidth(2);
          tempTH1qe->SetLineWidth(2);
          tempTH1other->SetLineWidth(2);
          tempTH12p2h->SetLineColor(kBlue);
          tempTH1dis->SetLineColor(kGreen);
          tempTH1res->SetLineColor(kYellow);
          tempTH1qe->SetLineColor(kRed);
          tempTH1other->SetLineColor(kViolet);
          //tempTH12p2h->Draw("HIST SAME");
          //tempTH1dis->Draw("HIST SAME");
          //tempTH1res->Draw("HIST SAME");
          //tempTH1qe->Draw("HIST SAME");
          //tempTH1other->Draw("HIST SAME");
          pad1->cd();
          tempxsec->Draw("SAME HIST P E1 X0"); //Draw it back on top
          can.cd();
          can.Print((prefix+"_"+nutypestr+"_CrossSectionWithAnezkas.png").c_str());
        }


        delete fluxIntegralC;
        delete fluxIntegralFe;
        delete fluxIntegralPb; 
        delete DaisyCorrectedC;
        delete DaisyCorrectedFe;
        delete DaisyCorrectedPb; 
        delete DaisyCorrectedSimC;
        delete DaisyCorrectedSimFe;
        delete DaisyCorrectedSimPb; 
        delete DaisyCorrectedSim2p2hC;
        delete DaisyCorrectedSim2p2hFe;
        delete DaisyCorrectedSim2p2hPb; 
        delete DaisyCorrectedSimDISC;
        delete DaisyCorrectedSimDISFe;
        delete DaisyCorrectedSimDISPb; 
        delete DaisyCorrectedSimRESC;
        delete DaisyCorrectedSimRESFe;
        delete DaisyCorrectedSimRESPb; 
        delete DaisyCorrectedSimQEC;
        delete DaisyCorrectedSimQEFe;
        delete DaisyCorrectedSimQEPb;
        delete DaisyCorrectedSimOtherC;
        delete DaisyCorrectedSimOtherFe;
        delete DaisyCorrectedSimOtherPb; 
        
        delete simEventRate;
        delete simEventRate2P2H;
        delete simEventRateDIS;
        delete simEventRateRES;
        delete simEventRateQE;
        delete simEventRateOther;

        delete unfolded;
        //delete crossSection;
      }


      // return 0;
      for (auto hist : backgrounds) delete hist;


      delete flux;
      delete folded;
      delete migration;
      delete effNum;
      delete effDenom;
      delete effDenom2P2H;
      delete effDenomDIS;
      delete effDenomRES;
      delete effDenomQE;
      delete effDenomOther;
      delete BackgroundWrongSign;
      delete BackgroundNC;
      delete BackgroundOther;
      delete SelectedSignalReco;
      delete MCData;
      delete MCData2p2h;
      delete MCDataDIS;
      delete MCDataRES;
      delete MCDataQE;
      delete MCDataOther;
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
