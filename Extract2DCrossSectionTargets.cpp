
#define HELP \
"\n*** Help: ***\n"\
" File: Extract1DCrossSectionTargets_ByTargetNew.cpp\n"\
" Brief: Given data and MC files produced by runEventLoopTargets.cpp, extract a 1D differential cross section.\n"\
"        Automatically merges input root files. Directories must have structure of /playlist/runEventLoopTargets<MC-or-Data>.root\n"\
"        Example directory structure:\n\n"\
"        /path/to/dirs/\n"\
"        ├──1A\n"\
"        │   ├──runEventLoopTargetsData2026.root\n"\
"        │   ├──runEventLoopTargetsMC2026.root\n"\
"        │   ├──runEventLoopTargetsData2082.root\n"\
"        │   ├──runEventLoopTargetsMC2082.root\n"\
"        │   ├──...\n"\
"        ├──1B\n"\
"        │   ├──runEventLoopTargetsData2026.root\n"\
"        │   ├──runEventLoopTargetsMC2026.root\n"\
"        │   ├──runEventLoopTargetsData2082.root\n"\
"        │   ├──runEventLoopTargetsMC2082.root\n"\
"        │   ├──...\n"\
"        ├──1C\n"\
"        │   ├──runEventLoopTargetsData2026.root\n"\
"        │   ├──runEventLoopTargetsMC2026.root\n"\
"        │   ├──runEventLoopTargetsData2082.root\n"\
"        │   ├──runEventLoopTargetsMC2082.root\n"\
"        │   ├──...\n"\
"        └──.....\n\n"\
"        Automatically identifies and searches all playlist directories within /path/to/dirs/ and\n"\
"        merges files for appropriate neutrino beam mode\n"\
"        Performs plastic sideband scaling, subtracts backgrounds, performs unfolding, applies efficiency x \n"\
"        acceptance correction, divides by flux and number of nucleons and if option is selected .\n"\
"        Writes a .root file with the cross section histograms\n"\
"        To run for a single playlist (for example 1A) simply pass /path/to/dirs/1A as the directory path\n\n"\
" Usage: Extract1DCrossSectionTargets_ByTargetNew <unfolding iterations> <directory> <target> <pdg>\n"\
"        e.g:   Extract1DCrossSectionTargets_ByTargetNew 5 /path/to/dirs 2026 14 -- to extract xsecs for target 2 Iron over neutrino-mode playlists with 5 iterations\n"\
"        e.g:   Extract1DCrossSectionTargets_ByTargetNew 10 /path/to/dirs 6000 -14 -- to extract xsecs for the water target over antineutrino-mode playlists with 10 iterations\n"\
"        e.g:   Extract1DCrossSectionTargets_ByTargetNew 10 /path/to/dirs iron -14 -- to extract xsecs for the combined iron targets over the antineutrino-mode playlists with 10 iterations\n"\
"        e.g:   Extract1DCrossSectionTargets_ByTargetNew 5 /path/to/dirs all 14 -- to extract xsecs for all targets over the neutrino-mode playlists with 5 iterations\n"\
"        e.g:   Extract1DCrossSectionTargets_ByTargetNew 10 /path/to/dirs WaterFull -14 -- to extract xsecs for the water target over the antineutrino-mode playlists with the water target filled with 10 iterations\n"\
"        e.g:   Extract1DCrossSectionTargets_ByTargetNew 10 /path/to/dirs WaterEmpty -14 -- to extract xsecs for the water target over the antineutrino-mode playlists with the water target empty with 10 iterations\n\n"\

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
#include "TLatex.h"
#include "TColor.h"
#include "string.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
//#include "TStyle.h"

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
#include <filesystem>

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
  double Nucleons;

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


PlotUtils::MnvH1D *g_combinedMCHist, *g_dataHist;


double getChi2( const double * val )
{
    double scale  = val[0];
    int ndof;
    PlotUtils::MnvPlotter plotter;
    plotter.ApplyStyle(PlotUtils::kCCQENuInclusiveStyle);
    auto tmpHist = g_dataHist->Clone();
    double chiSq = plotter.Chi2DataMC(tmpHist, g_combinedMCHist, ndof, scale, true, true);
    return chiSq;
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
  std::vector<std::string> filepathBases; //Excluding the ending eg "MC2006.root"
  if (argc != 5)
  {
    std::cerr << "Expected 4 arguments, but I got " << argc - 1 << ".\n" << HELP << std::endl;
    return 1;
  }
  const int nIterations = std::stoi(argv[1]);
  std::string indir = std::string(argv[2]);
  std::string intgt = std::string(argv[3]);
  int pdg = std::stoi(argv[4]);
  
  std::vector<std::string> dirs = util::findContainingDirectories(indir, "Targets", true);

  if (dirs.size()==0) dirs = {"./"};

  std::vector<std::string> targets;
  if (intgt == "ALL" || intgt == "all")
  {
    targets = {"1026", "1082", "2026", "2082", "3006", "3026", "3082", "4082", "5026", "5082", "6000", "Iron", "Carbon", "Lead", "WaterFull", "WaterEmpty"};
  }
  else
  {
    targets = {intgt};
  }
  for (std::string &tgt : targets)
  {
    std::cout<<"Working on target " << tgt << std::endl;
    std::vector<std::string> crossSectionPrefixes = {"pTmu_pZmu"};
    //std::vector<std::string> crossSectionPrefixes = {"pTmu_pZmu", "Emu_Erecoil"};

    for (const auto &prefix : crossSectionPrefixes)
    {
      std::cout<<"Working on prefix " << prefix << std::endl;

      PlotUtils::MnvH2D *fluxIntReweighted = new PlotUtils::MnvH2D();

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

      PlotUtils::MnvH1D* USSidebandSignal = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* USSidebandDS = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* USSidebandUS = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* USSidebandOther = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* DSSidebandSignal = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* DSSidebandDS = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* DSSidebandUS = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* DSSidebandOther = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* DataUSSideband = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* DataDSSideband = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* DataSignal = new PlotUtils::MnvH1D();

      PlotUtils::MnvH2D* SelectedSignalReco = new PlotUtils::MnvH2D();

      PlotUtils::MnvH2D* BackgroundWrongSign = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* BackgroundNC = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* BackgroundOther = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* BackgroundWaterTank = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* BackgroundUSPlastic = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* BackgroundDSPlastic = new PlotUtils::MnvH2D();

      PlotUtils::MnvH2D* MCData = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* MCData2p2h = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* MCDataDIS = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* MCDataRES = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* MCDataQE = new PlotUtils::MnvH2D();
      PlotUtils::MnvH2D* MCDataOther = new PlotUtils::MnvH2D();

      //if (!(prefix == "Erecoil" || prefix == "pTmu")) continue; //Used for testing with only subset of prefixes
      std::cout << "Currently working on prefix: " << prefix << std::endl;
      try
      {
        int waterFilledEmpty = -1; //Only used with water target

        double waterFilledPOTData = 0;
        double waterEmptyPOTData = 0;
        double waterFilledPOTMC = 0;
        double waterEmptyPOTMC = 0;
        std::cout<<"Investigating directory\n";
        std::vector <std::string> targetsInTgt;
        if (tgt=="Iron") targetsInTgt = {"2026", "3026", "5026"};
        else if (tgt=="Carbon") targetsInTgt = {"3006"};
        else if (tgt=="Lead") targetsInTgt = {"2082", "3082", "4082", "5082"};
        else if (tgt=="WaterEmpty")
        {
          targetsInTgt = {"6000"};
          waterFilledEmpty = 2;
        }
        else if (tgt=="WaterFull")
        {
          targetsInTgt = {"6000"};
          waterFilledEmpty = 1;
        }
        else if (tgt=="Water")
        {
          targetsInTgt = {"6000"};
          waterFilledEmpty = 0; //Do subtractive analysis
          for (int c = 0; c<dirs.size(); c++)
          {
            std::cout<<"Investigating directory " << dirs[c] << " which is " << c+1 <<"/"<<dirs.size() <<" playlists identified" <<std::endl;
            std::string datapath = dirs[c] + "/runEventLoopTargetsData" + targetsInTgt[0] + ".root";
            std::string mcpath = dirs[c] + "/runEventLoopTargetsMC" + targetsInTgt[0] + ".root";

            auto dataFile = TFile::Open(datapath.c_str(), "READ");
            if (!dataFile)
            {
              std::cerr << "Failed to open data file " << datapath.c_str() << ".\n";
              continue;
              return 2;
            }
            double data_pot = util::GetIngredient<TParameter<double>>(*dataFile, "POTUsed")->GetVal();
            auto mcFile = TFile::Open(mcpath.c_str(), "READ");
            if (!mcFile)
            {
              std::cerr << "Failed to open MC file " << mcpath.c_str() << ".\n";
              continue;
              return 3;
            }
            double mc_pot = util::GetIngredient<TParameter<double>>(*mcFile, "POTUsed")->GetVal();
            std::string playlistUsed = util::GetIngredient<TNamed>(*mcFile, "PlaylistUsed")->GetTitle();
            int filledorempty = util::filledOrEmptyMEPlaylist(playlistUsed);
            if (filledorempty == 2) //Empty
            {
              waterEmptyPOTData += data_pot;
              waterEmptyPOTMC += mc_pot;
            }
            if (filledorempty == 1) //Filled
            {
              waterFilledPOTData += data_pot;
              waterFilledPOTMC += mc_pot;
            }
          }
          std::cout<<"Total MC POT Filled: "<< waterFilledPOTMC << std::endl;
          std::cout<<"Total MC POT Empty: "<< waterEmptyPOTMC << std::endl;
          std::cout<<"Total Data POT Filled: "<< waterFilledPOTData << std::endl;
          std::cout<<"Total Data POT Empty: "<< waterEmptyPOTData << std::endl;
        }
        else targetsInTgt = {tgt};
        std::cout<<"Investigating directory1\n";
        //Flux parameters
        int n_flux_universes = 100; // Is this right
        const bool use_nue_constraint = true;
        const std::string project_dir = "targets_2345_jointNueIMD";
        double min_energy = 0;
        double max_energy = 100;

        double mcPOT = 0;
        double dataPOT = 0;
        std::cout<<"Investigating directory2\n";
        //Only needed for water subtractive analysis
        double scaleFactor = 1e18;

        double totalMCPOTFull = 0; //Only needed for water subtractive analysis
        double totalDataPOTFull = 0;
        double totalMCPOTEmpty = 0; //Only needed for water subtractive analysis
        double totalDataPOTEMpt = 0;
        double culmulativeMCPOT = 0; //Only needed for water subtractive analysis
        double culmulativeDataPOT = 0;

        std::cout<<"dirs " << dirs.size()<<"\n";
        for (int c = 0; c<dirs.size(); c++)
        {
          std::cout<<"Investigating directory " << dirs[c] << " which is " << c+1 <<"/"<<dirs.size() <<" playlists identified" <<std::endl;
          for (int d = 0; d<targetsInTgt.size(); d++)
          {
            std::cout<<"Investigating target " << targetsInTgt[d] << " which is " << d+1 <<"/"<<targetsInTgt.size() <<" targets selected" <<std::endl;
            std::string datapath = dirs[c] + "/runEventLoopTargetsData" + targetsInTgt[d] + ".root";
            std::string mcpath = dirs[c] + "/runEventLoopTargetsMC" + targetsInTgt[d] + ".root";
            std::string migpath = dirs[c] + "/runEventLoopTargets2DMigration" + targetsInTgt[d] + ".root";

            auto dataFile = TFile::Open(datapath.c_str(), "READ");
            if (!dataFile)
            {
              std::cerr << "Failed to open data file " << datapath.c_str() << ".\n";
              continue;
              return 2;
            }

            auto mcFile = TFile::Open(mcpath.c_str(), "READ");
            if (!mcFile)
            {
              std::cerr << "Failed to open MC file " << mcpath.c_str() << ".\n";
              continue;
              return 3;
            }

            auto migFile = TFile::Open(migpath.c_str(), "READ");
            if (!migFile)
            {
              std::cerr << "Failed to open MC file " << migpath.c_str() << ".\n";
              continue;
              return 3;
            }
            
            std::string playlistUsed = util::GetIngredient<TNamed>(*mcFile, "PlaylistUsed")->GetTitle();
            int filledorempty = util::filledOrEmptyMEPlaylist(playlistUsed);
            std::cout<<"playlistUsed: " << playlistUsed << " filledorempty " << filledorempty <<std::endl;

            int nuoranu = util::nuOrAntiNuMode(playlistUsed);
            int nupdg;
            if (nuoranu==1) nupdg = 14;
            else if (nuoranu==2) nupdg = -14;

            if (waterFilledEmpty>0 && (waterFilledEmpty != filledorempty)) continue;
            if (pdg!=nupdg)
            {
              std::cout<<"Skipping this set of files because this playlist pdg is " << nupdg << " but you specified " << pdg << std::endl;
              continue;
            }

            double mcscale = 1; //Only needed for water subtractive analysis
            double datascale = 1;

              //We scale to a standard POT, perform the necessary addition of subtraction for each playlist, keeping track of the culmulative POT (Standard x number of playlists), then scale back down to the actual total POT
            if (waterFilledEmpty == 0 && filledorempty==2)
            {
              mcscale = -(waterFilledPOTMC/waterEmptyPOTMC);
              datascale = -(waterFilledPOTData/waterEmptyPOTData);
              std::cout<<"datascale: " << datascale << std::endl;
              std::cout<<"mcscale: " << mcscale << std::endl;
            }


            PlotUtils::MnvH2D* tmpfolded = util::GetIngredient<PlotUtils::MnvH2D>(*dataFile, (std::string("data")), prefix);
            PlotUtils::MnvH2D* tmpmigration = util::GetIngredient<PlotUtils::MnvH2D>(*migFile, (std::string("migration")), prefix);
            PlotUtils::MnvH2D* tmpmigration_reco = util::GetIngredient<PlotUtils::MnvH2D>(*migFile, (std::string("reco")), prefix);
            PlotUtils::MnvH2D* tmpmigration_truth = util::GetIngredient<PlotUtils::MnvH2D>(*migFile, (std::string("truth")), prefix);
            PlotUtils::MnvH2D* tmpeffNum = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_numerator")), prefix);
            PlotUtils::MnvH2D* tmpeffDenom = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator")), prefix);
            PlotUtils::MnvH2D* tmpeffDenom2P2H = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_2p2h")), prefix);
            PlotUtils::MnvH2D* tmpeffDenomDIS = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_DIS")), prefix);
            PlotUtils::MnvH2D* tmpeffDenomRES = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_RES")), prefix);
            PlotUtils::MnvH2D* tmpeffDenomQE = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_QE")), prefix);
            PlotUtils::MnvH2D* tmpeffDenomOther = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_Other")), prefix);
            PlotUtils::MnvH1D* tmpUSSidebandSignal = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_US_sideband_Signal")));
            PlotUtils::MnvH1D* tmpUSSidebandDS = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_US_sideband_DS")));
            PlotUtils::MnvH1D* tmpUSSidebandUS = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_US_sideband_US")));
            PlotUtils::MnvH1D* tmpUSSidebandOther = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_US_sideband_Other")));
            PlotUtils::MnvH1D* tmpDSSidebandSignal = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_DS_sideband_Signal")));
            PlotUtils::MnvH1D* tmpDSSidebandDS = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_DS_sideband_DS")));
            PlotUtils::MnvH1D* tmpDSSidebandUS = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_DS_sideband_US")));
            PlotUtils::MnvH1D* tmpDSSidebandOther = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_DS_sideband_Other")));
            PlotUtils::MnvH1D* tmpDataUSSideband = util::GetIngredient<PlotUtils::MnvH1D>(*dataFile, (std::string("segment_US_Sideband")));
            PlotUtils::MnvH1D* tmpDataDSSideband = util::GetIngredient<PlotUtils::MnvH1D>(*dataFile, (std::string("segment_DS_Sideband")));
            PlotUtils::MnvH1D* tmpDataSignal = util::GetIngredient<PlotUtils::MnvH1D>(*dataFile, (std::string("segment_data")));
            PlotUtils::MnvH2D* tmpSelectedSignalReco = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("selected_signal_reco")), prefix);
            PlotUtils::MnvH2D* tmpBackgroundWrongSign = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("by_BKG_Label_Wrong_Sign_Bkg")), prefix);
            PlotUtils::MnvH2D* tmpBackgroundNC = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("by_BKG_Label_NC_Bkg")), prefix);
            PlotUtils::MnvH2D* tmpBackgroundWaterTank = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("by_BKG_Label_Water_Tank_Bkg")), prefix);
            PlotUtils::MnvH2D* tmpBackgroundDSPlastic = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("by_BKG_Label_Downstream_Plastic_Bkg")), prefix);
            PlotUtils::MnvH2D* tmpBackgroundUSPlastic = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("by_BKG_Label_Upstream_Plastic_Bkg")), prefix);
            PlotUtils::MnvH2D* tmpBackgroundOther = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("by_BKG_Label_Other")), prefix);
            PlotUtils::MnvH2D* tmpMCData = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("data")), prefix);
            PlotUtils::MnvH2D* tmpMCData2p2h = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("intType_2p2h")), prefix);
            PlotUtils::MnvH2D* tmpMCDataDIS = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("intType_DIS")), prefix);
            PlotUtils::MnvH2D* tmpMCDataRES = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("intType_RES")), prefix);
            PlotUtils::MnvH2D* tmpMCDataQE = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("intType_QE")), prefix);
            PlotUtils::MnvH2D* tmpMCDataOther = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("intType_Other")), prefix);

            util::AddHist(*folded,tmpfolded, datascale);
            util::AddHist(*migration,tmpmigration, mcscale);
            util::AddHist(*migration_reco,tmpmigration_reco, mcscale);
            util::AddHist(*migration_truth,tmpmigration_truth, mcscale);
            util::AddHist(*effNum,tmpeffNum, mcscale);
            util::AddHist(*effDenom,tmpeffDenom, mcscale);
            util::AddHist(*effDenom2P2H,tmpeffDenom2P2H, mcscale);
            util::AddHist(*effDenomDIS,tmpeffDenomDIS, mcscale);
            util::AddHist(*effDenomRES,tmpeffDenomRES, mcscale);
            util::AddHist(*effDenomQE,tmpeffDenomQE, mcscale);
            util::AddHist(*effDenomOther,tmpeffDenomOther, mcscale);
            util::AddHist(*USSidebandSignal,tmpUSSidebandSignal, mcscale);
            util::AddHist(*USSidebandDS,tmpUSSidebandDS, mcscale);
            util::AddHist(*USSidebandUS,tmpUSSidebandUS, mcscale);
            util::AddHist(*USSidebandOther,tmpUSSidebandOther, mcscale);
            util::AddHist(*DSSidebandSignal,tmpDSSidebandSignal, mcscale);
            util::AddHist(*DSSidebandDS,tmpDSSidebandDS, mcscale);
            util::AddHist(*DSSidebandUS,tmpDSSidebandUS, mcscale);
            util::AddHist(*DSSidebandOther,tmpDSSidebandOther, mcscale);
            util::AddHist(*DataUSSideband,tmpDataUSSideband, datascale);
            util::AddHist(*DataDSSideband,tmpDataDSSideband, datascale);
            util::AddHist(*DataSignal,tmpDataSignal, mcscale);
            util::AddHist(*SelectedSignalReco,tmpSelectedSignalReco, mcscale);
            util::AddHist(*BackgroundWrongSign,tmpBackgroundWrongSign, mcscale);
            util::AddHist(*BackgroundNC,tmpBackgroundNC, mcscale);
            util::AddHist(*BackgroundWaterTank,tmpBackgroundWaterTank, mcscale);
            util::AddHist(*BackgroundDSPlastic,tmpBackgroundDSPlastic, mcscale);
            util::AddHist(*BackgroundUSPlastic,tmpBackgroundUSPlastic, mcscale);
            util::AddHist(*BackgroundOther,tmpBackgroundOther, mcscale);
            util::AddHist(*MCData,tmpMCData, mcscale);
            util::AddHist(*MCData2p2h,tmpMCData2p2h, mcscale);
            util::AddHist(*MCDataDIS,tmpMCDataDIS, mcscale);
            util::AddHist(*MCDataRES,tmpMCDataRES, mcscale);
            util::AddHist(*MCDataQE,tmpMCDataQE, mcscale);
            util::AddHist(*MCDataOther,tmpMCDataOther, mcscale);

            delete tmpfolded;
            delete tmpmigration;
            delete tmpmigration_reco;
            delete tmpmigration_truth;
            delete tmpeffNum;
            delete tmpeffDenom;
            delete tmpeffDenom2P2H;
            delete tmpeffDenomDIS;
            delete tmpeffDenomRES;
            delete tmpeffDenomQE;
            delete tmpeffDenomOther;
            delete tmpUSSidebandSignal;
            delete tmpUSSidebandDS;
            delete tmpUSSidebandUS;
            delete tmpUSSidebandOther;
            delete tmpDSSidebandSignal;
            delete tmpDSSidebandDS;
            delete tmpDSSidebandUS;
            delete tmpDSSidebandOther;
            delete tmpDataUSSideband;
            delete tmpDataDSSideband;
            delete tmpDataSignal;
            delete tmpSelectedSignalReco;
            delete tmpBackgroundWrongSign;
            delete tmpBackgroundNC;
            delete tmpBackgroundWaterTank;
            delete tmpBackgroundDSPlastic;
            delete tmpBackgroundUSPlastic;
            delete tmpBackgroundOther;
            delete tmpMCData;
            delete tmpMCData2p2h;
            delete tmpMCDataDIS;
            delete tmpMCDataRES;
            delete tmpMCDataQE;
            delete tmpMCDataOther;

            //util::AddHist(*flux,util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("reweightedflux_integrated")), prefix), mcscale);
            /* util::AddHist(*folded,util::GetIngredient<PlotUtils::MnvH2D>(*dataFile, (std::string("data")), prefix), datascale);
            util::AddHist(*migration,util::GetIngredient<PlotUtils::MnvH2D>(*migFile, (std::string("migration")), prefix), mcscale);
            util::AddHist(*migration_reco,util::GetIngredient<PlotUtils::MnvH2D>(*migFile, (std::string("reco")), prefix), mcscale);
            util::AddHist(*migration_truth,util::GetIngredient<PlotUtils::MnvH2D>(*migFile, (std::string("truth")), prefix), mcscale);
            util::AddHist(*effNum,util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_numerator")), prefix), mcscale);
            util::AddHist(*effDenom,util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator")), prefix), mcscale);
            util::AddHist(*effDenom2P2H,util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_2p2h")), prefix), mcscale);
            util::AddHist(*effDenomDIS,util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_DIS")), prefix), mcscale);
            util::AddHist(*effDenomRES,util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_RES")), prefix), mcscale);
            util::AddHist(*effDenomQE,util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_QE")), prefix), mcscale);
            util::AddHist(*effDenomOther,util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_Other")), prefix), mcscale);
            util::AddHist(*USSidebandSignal,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_US_sideband_Signal"))), mcscale);
            util::AddHist(*USSidebandDS,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_US_sideband_DS"))), mcscale);
            util::AddHist(*USSidebandUS,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_US_sideband_US"))), mcscale);
            util::AddHist(*USSidebandOther,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_US_sideband_Other"))), mcscale);
            util::AddHist(*DSSidebandSignal,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_DS_sideband_Signal"))), mcscale);
            util::AddHist(*DSSidebandDS,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_DS_sideband_DS"))), mcscale);
            util::AddHist(*DSSidebandUS,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_DS_sideband_US"))), mcscale);
            util::AddHist(*DSSidebandOther,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_DS_sideband_Other"))), mcscale);
            util::AddHist(*DataUSSideband,util::GetIngredient<PlotUtils::MnvH1D>(*dataFile, (std::string("segment_US_Sideband"))), datascale);
            util::AddHist(*DataDSSideband,util::GetIngredient<PlotUtils::MnvH1D>(*dataFile, (std::string("segment_DS_Sideband"))), datascale);
            util::AddHist(*DataSignal,util::GetIngredient<PlotUtils::MnvH1D>(*dataFile, (std::string("segment_data"))), datascale);
            util::AddHist(*SelectedSignalReco,util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("selected_signal_reco")), prefix), mcscale);
            util::AddHist(*BackgroundWrongSign,util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("by_BKG_Label_Wrong_Sign_Bkg")), prefix), mcscale);
            util::AddHist(*BackgroundNC,util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("by_BKG_Label_NC_Bkg")), prefix), mcscale);
            util::AddHist(*BackgroundWaterTank,util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("by_BKG_Label_Water_Tank_Bkg")), prefix), mcscale);
            util::AddHist(*BackgroundDSPlastic,util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("by_BKG_Label_Downstream_Plastic_Bkg")), prefix), mcscale);
            util::AddHist(*BackgroundUSPlastic,util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("by_BKG_Label_Upstream_Plastic_Bkg")), prefix), mcscale);
            util::AddHist(*BackgroundOther,util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("by_BKG_Label_Other")), prefix), mcscale);
            util::AddHist(*MCData,util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("data")), prefix), mcscale);
            util::AddHist(*MCData2p2h,util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("intType_2p2h")), prefix), mcscale);
            util::AddHist(*MCDataDIS,util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("intType_DIS")), prefix), mcscale);
            util::AddHist(*MCDataRES,util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("intType_RES")), prefix), mcscale);
            util::AddHist(*MCDataQE,util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("intType_QE")), prefix), mcscale);
            util::AddHist(*MCDataOther,util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("intType_Other")), prefix), mcscale); */

            if (d==0) //Only get POT once per playlist. E.g we don't wanna double count the POT for the same runs just because we're looking at target 2026 and 3026
            {
              mcPOT += util::GetIngredient<TParameter<double>>(*mcFile, "POTUsed")->GetVal();
              std::cout<<"mcPOT: " << mcPOT << std::endl;
              double tempDataPOT = util::GetIngredient<TParameter<double>>(*dataFile, "POTUsed")->GetVal();
              dataPOT+=tempDataPOT;

              //Scale the integrated flux for this playlist by playlist data POT for appropriate scaling when applying across different playlists
              PlotUtils::FluxReweighter frw = PlotUtils::FluxReweighter( nupdg, use_nue_constraint, playlistUsed, PlotUtils::FluxReweighter::gen2thin, PlotUtils::FluxReweighter::g4numiv6, n_flux_universes );                    
              auto tempIntFlux = frw.GetIntegratedFluxReweighted(nupdg, effDenom, min_energy, max_energy, true)->Clone();
              //tempIntFlux->Scale(tempDataPOT);
              util::AddHist(*fluxIntReweighted,tempIntFlux,tempDataPOT);
            }

            dataFile->Close();
            mcFile->Close();
          }
        }

        std::cout<<"mcPOT: " << mcPOT << std::endl;
        std::cout<<"dataPOT: " << dataPOT << std::endl;

        //Normalising the integrated flux across different playlists by dataPOT
        fluxIntReweighted->Scale(1/dataPOT);
        Plot(*folded, "folded", prefix, tgt);
        //************************************************
        //Sideband chisq minimisation - single scale factor across a given sideband region
        //************************************************

        bool doMinimiser = true;
        double USScaleFactor = 1;
        double DSScaleFactor = 1;
        if (doMinimiser)
        {
          PlotUtils::MnvH1D* combinedUS = USSidebandSignal->Clone("CombinedUSMC");
          combinedUS->Add(USSidebandUS);
          combinedUS->Add(USSidebandDS);
          combinedUS->Add(USSidebandOther);
          combinedUS->Scale(dataPOT/mcPOT);

          PlotUtils::MnvH1D* combinedDS = DSSidebandSignal->Clone("CombinedUSMC");
          combinedDS->Add(DSSidebandUS);
          combinedDS->Add(DSSidebandDS);
          combinedDS->Add(DSSidebandOther);
          combinedDS->Scale(dataPOT/mcPOT);

          g_combinedMCHist = combinedUS;
          g_dataHist = DataUSSideband;

          //Plot(*combinedUS, "combinedUS", prefix, tgt);
          //Plot(*DataUSSideband, "DataUSSideband", prefix, tgt);
          double minValUS, minValDS;
          {          //US Sideband
            ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2");
            minimizer->SetMaxFunctionCalls(1000000);
            minimizer->SetMaxIterations(100000);
            minimizer->SetTolerance(0.01);
            minimizer->SetPrintLevel(1);
            ROOT::Math::Functor f( & getChi2, 1 );
            minimizer->SetFunction(f);
            minimizer->SetVariable(0, "scale", 1, 0.001);
            minimizer->Minimize();
            const double *bestfit = minimizer->X();
            double bestChiSq = minimizer->MinValue();
            USScaleFactor = bestfit[0];
            std::cout<<"Minised chisq for USSideband. Scale factor: " << USScaleFactor << " best chisq: " << bestChiSq << std::endl;
          }

          g_combinedMCHist = combinedDS;
          g_dataHist = DataDSSideband;
          {          //DS Sideband
            ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2");
            minimizer->SetMaxFunctionCalls(1000000);
            minimizer->SetMaxIterations(100000);
            minimizer->SetTolerance(0.01);
            minimizer->SetPrintLevel(1);
            ROOT::Math::Functor f( & getChi2, 1 );
            minimizer->SetFunction(f);
            minimizer->SetVariable(0, "scale", 1, 0.001);
            minimizer->Minimize();
            const double *bestfit = minimizer->X();
            double bestChiSq = minimizer->MinValue();
            DSScaleFactor = bestfit[0];
            std::cout<<"Minised chisq for DSSideband. Scale factor: " << DSScaleFactor << " best chisq: " << bestChiSq << std::endl;
          }

        }
        //************************************************
        //End - Sideband chisq minimisation - single scale factor across a given sideband region
        //************************************************

        auto simEventRate = effDenom->Clone(); // Make a copy for later
        auto simEventRate2P2H = effDenom2P2H->Clone(); // Make a copy for later
        auto simEventRateDIS = effDenomDIS->Clone(); // Make a copy for later
        auto simEventRateRES = effDenomRES->Clone(); // Make a copy for later
        auto simEventRateQE = effDenomQE->Clone(); // Make a copy for later
        auto simEventRateOther = effDenomOther->Clone(); // Make a copy for later
        // There are no error bands in the data, but I need somewhere to put error bands on the results I derive from it.
        folded->AddMissingErrorBandsAndFillWithCV(*migration);

        // Basing my unfolding procedure for a differential cross section on Alex's MINERvA 101 talk at https://minerva-docdb.fnal.gov/cgi-bin/private/RetrieveFile?docid=27438&filename=whatsACrossSection.pdf&version=1

        // TODO: Remove these debugging plots when done


        USScaleFactor = 1.05;
        DSScaleFactor = 1.05;

        BackgroundUSPlastic->Scale(USScaleFactor);
        BackgroundDSPlastic->Scale(DSScaleFactor);
        std::vector<PlotUtils::MnvH2D*> sidebandScaled = {BackgroundUSPlastic, BackgroundDSPlastic};
        std::vector<PlotUtils::MnvH2D*> backgrounds = { BackgroundWrongSign, BackgroundNC, BackgroundOther, BackgroundWaterTank};
        auto toSubtract = std::accumulate(std::next(backgrounds.begin()), backgrounds.end(), (*backgrounds.begin())->Clone(),
                                          [](auto sum, const auto hist)
                                          {
                                            sum->Add(hist);
                                            return sum;
                                          });
        Plot(*toSubtract, "BackgroundSum", prefix, tgt);
        auto toSubtractSideband = std::accumulate(std::next(sidebandScaled.begin()), sidebandScaled.end(), (*sidebandScaled.begin())->Clone(),
                                          [](auto sum, const auto hist)
                                          {
                                            sum->Add(hist);
                                            return sum;
                                          });
        Plot(*toSubtractSideband, "sidebandScaledSum", prefix, tgt);

        //sum up backgrounds without sideband scaling
        auto bkgSubtracted = std::accumulate(backgrounds.begin(), backgrounds.end(), folded->Clone(),
                                              [mcPOT, dataPOT](auto sum, const auto hist)
                                              {
                                                std::cout << "Subtracting " << hist->GetName() << " scaled by " << -dataPOT / mcPOT << " from " << sum->GetName() << "\n";
                                                sum->Add(hist, -dataPOT / mcPOT);
                                                return sum;
                                              });

        //remove scaled sideband
        auto sidebandScaledSubtracted = std::accumulate(sidebandScaled.begin(), sidebandScaled.end(), folded->Clone(),
                                              [mcPOT, dataPOT](auto sum, const auto hist)
                                              {
                                                std::cout << "Subtracting " << hist->GetName() << " scaled by " << -dataPOT / mcPOT << " from " << sum->GetName() << "\n";
                                                sum->Add(hist, -dataPOT / mcPOT);
                                                return sum;
                                              });
        //Both background subtracted and sideband scaled
        auto bkgScaledSubtracted = std::accumulate(sidebandScaled.begin(), sidebandScaled.end(), bkgSubtracted->Clone(),
                                              [mcPOT, dataPOT](auto sum, const auto hist)
                                              {
                                                std::cout << "Subtracting " << hist->GetName() << " scaled by " << -dataPOT / mcPOT << " from " << sum->GetName() << "\n";
                                                sum->Add(hist, -dataPOT / mcPOT);
                                                return sum;
                                              });

        Plot(*bkgSubtracted, "backgroundSubtracted", prefix, tgt);
        Plot(*bkgScaledSubtracted, "backgroundSubtractedScaled", prefix, tgt);
        Plot(*sidebandScaledSubtracted, "sidebandScaledSubtracted", prefix, tgt);

        auto outFile = TFile::Open((tgt + prefix + "_crossSection.root").c_str(), "RECREATE");
        if (!outFile)
        {
          std::cerr << "Could not create a file called " << prefix + "_crossSection.root" << ".  Does it already exist?\n";
          return 5;
        }

        toSubtractSideband->Clone()->Write("toSubtractSideband"); // TODO: Seg fault first appears when I uncomment this line
        toSubtract->Clone()->Write("toSubtract"); // TODO: Seg fault first appears when I uncomment this line

        bkgSubtracted->Clone()->Write("bkgSubtracted"); // TODO: Seg fault first appears when I uncomment this line
        bkgScaledSubtracted->Clone()->Write("bkgScaledSubtracted"); // TODO: Seg fault first appears when I uncomment this line
        sidebandScaledSubtracted->Clone()->Write("sidebandScaledSubtracted"); // TODO: Seg fault first appears when I uncomment this line

        folded->Clone()->Write("folded"); // TODO: Seg fault first appears when I uncomment this line

        DataUSSideband->Write("DataUSSideband");
        DataDSSideband->Write("DataDSSideband");

        USSidebandSignal->Write("USSidebandSignal");
        USSidebandUS->Write("USSidebandUS");
        USSidebandDS->Write("USSidebandDS");
        USSidebandOther->Write("USSidebandOther");
        
        DSSidebandSignal->Write("DSSidebandSignal");
        DSSidebandUS->Write("DSSidebandUS");
        DSSidebandDS->Write("DSSidebandDS");
        DSSidebandOther->Write("DSSidebandOther");

        bkgSubtracted->Write("backgroundSubtracted");
        bkgScaledSubtracted->Write("backgroundSubtractedSidebandScaled");


        // d'Aogstini unfolding
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Temporary to speed up plotting debugging
        auto unfolded = UnfoldHist(bkgSubtracted, migration, migration_reco, migration_truth, nIterations);
        if (!unfolded)
          throw std::runtime_error(std::string("Failed to unfold ") + folded->GetName() + " using " + migration->GetName());
        Plot(*unfolded, "unfolded", prefix, tgt);
        unfolded->Clone()->Write("unfolded"); // TODO: Seg fault first appears when I uncomment this line
        auto unfolded_tuned = UnfoldHist(bkgScaledSubtracted, migration, migration_reco, migration_truth, nIterations);
        if (!unfolded_tuned)
          throw std::runtime_error(std::string("Failed to unfold ") + folded->GetName() + " using " + migration->GetName());
        //////////////////////////////////////////Delete these 2 lines
        //auto unfolded = bkgSubtracted;
        //auto unfolded_tuned = bkgScaledSubtracted;
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! End Temporary to speed up plotting debugging

        Plot(*unfolded_tuned, "unfolded_sidebandTuned", prefix, tgt);
        unfolded_tuned->Clone()->Write("unfolded_sidebandTuned"); // TODO: Seg fault first appears when I uncomment this line
        std::cout << "Survived writing the unfolded histogram.\n"
                  << std::flush; // This is evidence that the problem is on the final file Write() and not unfolded->Clone()->Write().

        effNum->Divide(effNum, effDenom); // Only the 2 parameter version of MnvH1D::Divide()
                                          // handles systematics correctly.
        Plot(*effNum, "efficiency", prefix, tgt);

        unfolded->Divide(unfolded, effNum);
        Plot(*unfolded, "efficiencyCorrected", prefix, tgt);

        unfolded_tuned->Divide(unfolded_tuned, effNum);
        Plot(*unfolded_tuned, "efficiencyCorrected_sidebandTuned", prefix, tgt);

        // double nnucleons = nNucleons->GetVal()/numMergedPlaylists;
        double nnucleons = 0;

        std::string tgtname = tgt;
        const size_t base = tgt.find("Scaled");
        if (base != std::string::npos)
        {
          tgtname= tgt.substr(0, base);
        }

        std::string material = "tracker";
        if(tgtname == "3006") material = "carbon";
        else if(tgtname == "1026" || tgtname == "2026" ||tgtname == "3026" || tgtname == "5026" ) material = "iron";
        else if(tgtname == "1082" || tgtname == "2082" ||tgtname == "3082" || tgtname == "4082" || tgtname == "5082" ) material = "lead";
        else if (tgtname == "Lead") material = "lead";
        else if (tgtname == "Iron") material = "iron";
        else if (tgtname == "Carbon") material = "carbon";
        //Is it ok to just let water use the tracker flux?
        
        double nnucleonsData;
        std::cout<<"material: " << material << std::endl;
        PlotUtils::TargetUtils targetInfo;

        if (tgtname == "Lead")
        {
          nnucleons = targetInfo.GetPassiveTargetNNucleons(2, 82, true);
          nnucleons += targetInfo.GetPassiveTargetNNucleons(3, 82, true);
          nnucleons += targetInfo.GetPassiveTargetNNucleons(4, 82, true);
          nnucleons += targetInfo.GetPassiveTargetNNucleons(5, 82, true);
          nnucleonsData = targetInfo.GetPassiveTargetNNucleons(2, 82, false);
          nnucleonsData += targetInfo.GetPassiveTargetNNucleons(3, 82, false);
          nnucleonsData += targetInfo.GetPassiveTargetNNucleons(4, 82, false);
          nnucleonsData += targetInfo.GetPassiveTargetNNucleons(5, 82, false);
        }
        else if (tgtname == "Iron")
        {
          nnucleons = targetInfo.GetPassiveTargetNNucleons(2, 26, true);
          nnucleons += targetInfo.GetPassiveTargetNNucleons(3, 26, true);
          nnucleons += targetInfo.GetPassiveTargetNNucleons(5, 26, true);
          nnucleonsData = targetInfo.GetPassiveTargetNNucleons(2, 26, false);
          nnucleonsData += targetInfo.GetPassiveTargetNNucleons(3, 26, false);
          nnucleonsData += targetInfo.GetPassiveTargetNNucleons(5, 26, false);
        }
        else if (tgtname == "Carbon")
        {
          nnucleons = targetInfo.GetPassiveTargetNNucleons(3, 6, true);
          nnucleonsData = targetInfo.GetPassiveTargetNNucleons(3, 6, false);
        }
        else if (tgtname.find("Water") != std::string::npos)
        {
          nnucleons = targetInfo.GetPassiveTargetNNucleons(6, 1, true);
          nnucleonsData = targetInfo.GetPassiveTargetNNucleons(6, 1, false);
        }
        else
        {
          std::cout<<"Here1.n";
          int tgtCode = std::stoi(tgtname);
          int tgtMat = tgtCode%1000;
          int tgtNum = (tgtCode-tgtMat)/1000;
          if (tgtCode == 6000)
          {
            tgtMat = 1;
            tgtNum = 6;
          }
          if (tgtNum<7 && tgtCode>1000) nnucleons = targetInfo.GetPassiveTargetNNucleons(tgtNum, tgtMat, true);
          else
          {
            if (tgtNum==7) nnucleons = targetInfo.GetTrackerNNucleons( 7, true); 
            if (tgtNum==8) nnucleons = targetInfo.GetTrackerNNucleons( 6, true); 
            if (tgtNum==9) nnucleons = targetInfo.GetTrackerNNucleons( 6, true); 
            if (tgtNum==10) nnucleons = targetInfo.GetTrackerNNucleons( 6, true); 
            if (tgtNum==11) nnucleons = targetInfo.GetTrackerNNucleons( 6, true); 
            if (tgtNum==12) nnucleons = targetInfo.GetTrackerNNucleons( 2, true); 
          }
          if (tgtNum<7 && tgtCode>1000) nnucleonsData = targetInfo.GetPassiveTargetNNucleons(tgtNum, tgtMat, false);
          else
          {
            if (tgtNum==7) nnucleonsData = targetInfo.GetTrackerNNucleons( 7, false); 
            if (tgtNum==8) nnucleonsData = targetInfo.GetTrackerNNucleons( 6, false); 
            if (tgtNum==9) nnucleonsData = targetInfo.GetTrackerNNucleons( 6, false); 
            if (tgtNum==10) nnucleonsData = targetInfo.GetTrackerNNucleons( 6, false); 
            if (tgtNum==11) nnucleonsData = targetInfo.GetTrackerNNucleons( 6, false); 
            if (tgtNum==12) nnucleonsData = targetInfo.GetTrackerNNucleons( 2, false); 
          }
        }

        std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
        std::cout<<"2082 Nnucleons: " << targetInfo.GetPassiveTargetNNucleons(2, 82, true)<< "\n";
        std::cout<<"3082 Nnucleons: " << targetInfo.GetPassiveTargetNNucleons(3, 82, true)<< "\n";
        std::cout<<"4082 Nnucleons: " << targetInfo.GetPassiveTargetNNucleons(4, 82, true)<< "\n";
        std::cout<<"5082 Nnucleons: " << targetInfo.GetPassiveTargetNNucleons(5, 82, true)<< "\n";
        std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
        std::cout<<"2026 Nnucleons: " << targetInfo.GetPassiveTargetNNucleons(2, 26, true)<< "\n";
        std::cout<<"3026 Nnucleons: " << targetInfo.GetPassiveTargetNNucleons(3, 26, true)<< "\n";
        std::cout<<"5026 Nnucleons: " << targetInfo.GetPassiveTargetNNucleons(5, 26, true)<< "\n";
        std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";




        std::string tgtString = "";
        if (tgt == "WaterEmpty") tgtString = "Water (Empty Playlists - Approx 63% of ME POT)";
        else if (tgt == "WaterFull") tgtString = "Water (Filled Playlists - Approx 37% of ME POT)";
        else if (tgt == "Water") tgtString = "Water (Filled - Empty)";
        else if (tgt == "Iron") tgtString = "Iron";
        else if (tgt == "Carbon") tgtString = "Carbon";
        else if (tgt == "Lead") tgtString = "Lead";
        else
        {
          int tgtInt = std::stoi(tgt);
          int tgtZ = tgtInt%1000;
          int tgtNum = (tgtInt-tgtZ)/1000;
          if (tgtNum==6) tgtString = "Water Target";
          else if (tgtNum==1) tgtString = "Target 1 ";
          else if (tgtNum==2) tgtString = "Target 2 ";
          else if (tgtNum==3) tgtString = "Target 3 ";
          else if (tgtNum==4) tgtString = "Target 4 ";
          else  if (tgtNum==5) tgtString = "Target 5 ";
          if (tgtZ==26) tgtString += "Iron";
          if (tgtZ==6) tgtString += "Carbon";
          if (tgtZ==82) tgtString += "Lead";
        }

        PlotUtils::MnvH1D *flux2;
        PlotUtils::MnvH1D *fluxIntegral;
        PlotUtils::MnvH1D *fluxRebinned;

        //std::cout<<"nupdg: " << nupdg << std::endl;
        std::cout<<"use_nue_constraint: " << use_nue_constraint << std::endl;
        std::cout<<"n_flux_universes: " << n_flux_universes << std::endl;
        std::cout<<"min_energy: " << min_energy << std::endl;
        std::cout<<"max_energy: " << max_energy << std::endl;
        std::cout<<"nnucleons: " << nnucleons << std::endl;
        std::cout<<"mcPOT: " << mcPOT << std::endl;
        //PlotUtils::FluxReweighter frw = PlotUtils::flux_reweighter(platlistname, pdg, use_nue_constraint, n_flux_universes);
        //PlotUtils::FluxReweighter *frw = new PlotUtils::FluxReweighter( pdg, use_nue_constraint, platlistname, PlotUtils::FluxReweighter::gen2thin, PlotUtils::FluxReweighter::g4numiv6, n_flux_universes );
        //std::cout << "ABC123-1 " << std::endl;
        //PlotUtils::MnvH1D *fluxIntReweighted = frw->GetIntegratedFluxReweighted(pdg, simEventRate, min_energy, max_energy, true);
        //fluxIntegral = frw->GetIntegratedTargetFlux(pdg, material, simEventRate, min_energy, max_energy, project_dir);
        //auto &frw2 = PlotUtils::flux_reweighter("minervame1A", pdg, use_nue_constraint, n_flux_universes);
        //flux2 = frw2.GetTargetFluxMnvH1D(pdg, material, project_dir);

        //fluxRebinned = frw->GetRebinnedFluxReweighted_FromInputFlux(flux2, simEventRate); // issue here
        //PlotUtils::MnvH1D *Integrated_fluxGenie = frw->GetIntegratedFluxReweighted_FromInputFlux(flux2, simEventRate, min_energy, max_energy);
        //PlotUtils::MnvH1D *Integrated_fluxGenie2 = frw->GetIntegratedFluxReweighted_FromInputFlux(flux, simEventRate, min_energy, max_energy);
        outFile->cd();
        auto crossSection = normalize(unfolded, fluxIntReweighted, nnucleonsData, dataPOT);
        Plot(*crossSection, "crossSection", prefix, tgt);
        crossSection->Write("crossSection");
        auto crossSectionTuned = normalize(unfolded_tuned, fluxIntReweighted, nnucleonsData, dataPOT);
        Plot(*crossSectionTuned, "crossSectionTuned", prefix, tgt);
        crossSectionTuned->Write("crossSectionTuned");
        simEventRate->Write("simulatedEventRate");
        fluxIntReweighted->Write("fluxIntReweighted");

        //These lines are/were mainly used to debug a scaling issue
        /* flux->Write("flux(FRW/from event loop)");
          //flux2->Write("flux2");
          //flux2->Write("flux2");
        fluxIntegral->Write("fluxIntegral");
        //Integrated_fluxGenie->Write("Integrated_fluxGenie");
        Integrated_fluxGenie2->Write("Integrated_fluxGenie2");
        //fluxRebinned->Write("fluxRebinned");
        std::cout << "GetTotalScatteringCenters(26,true) " << GetTotalScatteringCenters(26, true) << std::endl;
        std::cout << "GetTotalScatteringCenters(26,false) " << GetTotalScatteringCenters(26, false) << std::endl;
        std::cout << "targetInfo.GetPassiveTargetNNucleons( 2, 26, true ): " << targetInfo.GetPassiveTargetNNucleons(2, 26, true) << std::endl;
        std::cout << "targetInfo.GetPassiveTargetNNucleons( 2, 26, false ): " << targetInfo.GetPassiveTargetNNucleons(2, 26, false) << std::endl;
        std::cout << "mcPOT " << mcPOT << std::endl;
        std::cout << "NNucleons " << nnucleons << std::endl; */
        // Write a "simulated cross section" to compare to the data I just extracted.
        // If this analysis passed its closure test, this should be the same cross section as
        // what GENIEXSecExtract would produce.
        //fluxIntegral about 1pc off
        //fluxIntReweighted passes

        auto simulatedCrossSection = normalize(simEventRate, fluxIntReweighted, nnucleons, mcPOT);
        //Plot(*simulatedCrossSection, "simulatedCrossSection", prefix, tgt);
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
        TParameter<double> potMC = TParameter<double>("mcPOT", mcPOT);
        TParameter<double> potData = TParameter<double>("dataPOT", dataPOT);
        potMC.Write();
        potData.Write();

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
          PlotUtils::MnvH2D* tempxsec = crossSectionTuned->Clone();
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
          
          std::string title = tgtString;

          if (tgtname == "Water")
          {
            std::stringstream POTstrFilled, POTstrEmpty;
            POTstrFilled << std::fixed << std::scientific << std::setprecision(2) << waterFilledPOTData;
            POTstrEmpty << std::fixed << std::scientific << std::setprecision(2) << waterEmptyPOTData;
            std::cout << "waterFilledPOTData: "<< waterFilledPOTData <<std::endl;
            std::cout << "waterEmptyPOTData: "<< waterEmptyPOTData <<std::endl;
            std::string potstrFilled = std::string("Water Filled POT ")+POTstrFilled.str();
            std::string potstrEmpty = std::string("Water Empty POT ")+POTstrEmpty.str();// + " (Scaled to match filled)";
            plotter.AddPlotLabel(potstrFilled.c_str() , 0.18, 0.99, 0.025, 1, 42, 13, 0);
            plotter.AddPlotLabel(potstrEmpty.c_str() , 0.18, 0.96, 0.025, 1, 42, 13, 0);
          }
          else
          {
            std::stringstream POTstr;
            POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOT;
            std::string potstr = std::string("POT Used: ")+POTstr.str();
            plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.97, 0.025, 1, 42, 13, 0);
          }

          plotter.WritePreliminary( 0.75, 0.075, 0.035, true);
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 35, 1, 43, 22, 0);
          //plotter.AddPlotLabel("Data/MC" , 0.08, 0.3, 0.03, 1, 42, 33, 90);

          plotter.AddPlotLabel("d^{2}#sigma/dp_{t}dp_{z} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)", 0.020 ,0.51, 30, 1, 43, 22, 90);
          plotter.AddPlotLabel("Muon Transverse Momentum (Gev/c)", 0.51, 0.025, 30, 1, 43);
          canvas->Print((tgtname+"_"+prefix+"_CrossSectionTransverse.png").c_str());
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
          PlotUtils::MnvH2D* tempxsec = crossSectionTuned->Clone();
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
          
          std::string title = tgtString;

          if (tgtname == "Water")
          {
            std::stringstream POTstrFilled, POTstrEmpty;
            POTstrFilled << std::fixed << std::scientific << std::setprecision(2) << waterFilledPOTData;
            POTstrEmpty << std::fixed << std::scientific << std::setprecision(2) << waterEmptyPOTData;
            std::cout << "waterFilledPOTData: "<< waterFilledPOTData <<std::endl;
            std::cout << "waterEmptyPOTData: "<< waterEmptyPOTData <<std::endl;
            std::string potstrFilled = std::string("Water Filled POT: ")+POTstrFilled.str();
            std::string potstrEmpty = std::string("Water Empty POT ")+POTstrEmpty.str();// + " (Scaled to match filled)";
            plotter.AddPlotLabel(potstrFilled.c_str() , 0.18, 0.99, 0.025, 1, 42, 13, 0);
            plotter.AddPlotLabel(potstrEmpty.c_str() , 0.18, 0.96, 0.025, 1, 42, 13, 0);
          }
          else
          {
            std::stringstream POTstr;
            POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOT;
            std::string potstr = std::string("POT Used: ")+POTstr.str();
            plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.97, 0.025, 1, 42, 13, 0);
          }

          plotter.WritePreliminary( 0.75, 0.075, 0.035, true);
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 35, 1, 43, 22, 0);
          //plotter.AddPlotLabel("Data/MC" , 0.08, 0.3, 0.03, 1, 42, 33, 90);

          plotter.AddPlotLabel("d^{2}#sigma/dp_{t}dp_{z} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)", 0.020 ,0.51, 30, 1, 43, 22, 90);
          plotter.AddPlotLabel("Muon Longitudinal Momentum (Gev/c)", 0.51, 0.025, 30, 1, 43);
          canvas->Print((tgtname+"_"+prefix+"_CrossSectionLongitudinal.png").c_str());
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
          PlotUtils::MnvH2D* tempxsec = crossSectionTuned->Clone();
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
          
          std::string title = tgtString;

          if (tgtname == "Water")
          {
            std::stringstream POTstrFilled, POTstrEmpty;
            POTstrFilled << std::fixed << std::scientific << std::setprecision(2) << waterFilledPOTData;
            POTstrEmpty << std::fixed << std::scientific << std::setprecision(2) << waterEmptyPOTData;
            std::cout << "waterFilledPOTData: "<< waterFilledPOTData <<std::endl;
            std::cout << "waterEmptyPOTData: "<< waterEmptyPOTData <<std::endl;
            std::string potstrFilled = std::string("Water Filled POT ")+POTstrFilled.str();
            std::string potstrEmpty = std::string("Water Empty POT ")+POTstrEmpty.str();// + " (Scaled to match filled)";
            plotter.AddPlotLabel(potstrFilled.c_str() , 0.18, 0.99, 0.025, 1, 42, 13, 0);
            plotter.AddPlotLabel(potstrEmpty.c_str() , 0.18, 0.96, 0.025, 1, 42, 13, 0);
          }
          else
          {
            std::stringstream POTstr;
            POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOT;
            std::string potstr = std::string("POT Used: ")+POTstr.str();
            plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.97, 0.025, 1, 42, 13, 0);
          }

          plotter.WritePreliminary( 0.75, 0.075, 0.035, true);
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 35, 1, 43, 22, 0);
          //plotter.AddPlotLabel("Data/MC" , 0.08, 0.3, 0.03, 1, 42, 33, 90);

          plotter.AddPlotLabel("d^{2}#sigma/dp_{t}dp_{z} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)", 0.020 ,0.51, 30, 1, 43, 22, 90);
          plotter.AddPlotLabel("Muon Transverse Momentum (Gev/c)", 0.51, 0.025, 30, 1, 43);
          canvas->Print((tgtname+"_"+prefix+"_DataMCRatioTransverse.png").c_str());
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
          PlotUtils::MnvH2D* tempxsec = crossSectionTuned->Clone();
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
          
          std::string title = tgtString;

          if (tgtname == "Water")
          {
            std::stringstream POTstrFilled, POTstrEmpty;
            POTstrFilled << std::fixed << std::scientific << std::setprecision(2) << waterFilledPOTData;
            POTstrEmpty << std::fixed << std::scientific << std::setprecision(2) << waterEmptyPOTData;
            std::cout << "waterFilledPOTData: "<< waterFilledPOTData <<std::endl;
            std::cout << "waterEmptyPOTData: "<< waterEmptyPOTData <<std::endl;
            std::string potstrFilled = std::string("Water Filled POT: ")+POTstrFilled.str();
            std::string potstrEmpty = std::string("Water Empty POT ")+POTstrEmpty.str();// + " (Scaled to match filled)";
            plotter.AddPlotLabel(potstrFilled.c_str() , 0.18, 0.99, 0.025, 1, 42, 13, 0);
            plotter.AddPlotLabel(potstrEmpty.c_str() , 0.18, 0.96, 0.025, 1, 42, 13, 0);
          }
          else
          {
            std::stringstream POTstr;
            POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOT;
            std::string potstr = std::string("POT Used: ")+POTstr.str();
            plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.97, 0.025, 1, 42, 13, 0);
          }

          plotter.WritePreliminary( 0.75, 0.075, 0.035, true);
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 35, 1, 43, 22, 0);
          //plotter.AddPlotLabel("Data/MC" , 0.08, 0.3, 0.03, 1, 42, 33, 90);

          plotter.AddPlotLabel("d^{2}#sigma/dp_{t}dp_{z} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)", 0.020 ,0.51, 30, 1, 43, 22, 90);
          plotter.AddPlotLabel("Muon Longitudinal Momentum (Gev/c)", 0.51, 0.025, 30, 1, 43);
          canvas->Print((tgtname+"_"+prefix+"_DataMCRatioLongitudinal.png").c_str());
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
          PlotUtils::MnvH2D* tempevrate = unfolded_tuned->Clone();
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
          
          std::string title = tgtString;

          if (tgtname == "Water")
          {
            std::stringstream POTstrFilled, POTstrEmpty;
            POTstrFilled << std::fixed << std::scientific << std::setprecision(2) << waterFilledPOTData;
            POTstrEmpty << std::fixed << std::scientific << std::setprecision(2) << waterEmptyPOTData;
            std::cout << "waterFilledPOTData: "<< waterFilledPOTData <<std::endl;
            std::cout << "waterEmptyPOTData: "<< waterEmptyPOTData <<std::endl;
            std::string potstrFilled = std::string("Water Filled POT ")+POTstrFilled.str();
            std::string potstrEmpty = std::string("Water Empty POT ")+POTstrEmpty.str();// + " (Scaled to match filled)";
            plotter.AddPlotLabel(potstrFilled.c_str() , 0.18, 0.99, 0.025, 1, 42, 13, 0);
            plotter.AddPlotLabel(potstrEmpty.c_str() , 0.18, 0.96, 0.025, 1, 42, 13, 0);
          }
          else
          {
            std::stringstream POTstr;
            POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOT;
            std::string potstr = std::string("POT Used: ")+POTstr.str();
            plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.97, 0.025, 1, 42, 13, 0);
          }

          plotter.WritePreliminary( 0.75, 0.075, 0.035, true);
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 35, 1, 43, 22, 0);
          //plotter.AddPlotLabel("Data/MC" , 0.08, 0.3, 0.03, 1, 42, 33, 90);

          plotter.AddPlotLabel("Events#times10^{3} /(Gev/c)^2", 0.020 ,0.51, 30, 1, 43, 22, 90);
          plotter.AddPlotLabel("Muon Transverse Momentum (Gev/c)", 0.51, 0.025, 30, 1, 43);
          canvas->Print((tgtname+"_"+prefix+"_EvRateTransverse.png").c_str());
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
          PlotUtils::MnvH2D* tempevrate = unfolded_tuned->Clone();
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
          
          std::string title = tgtString;

          if (tgtname == "Water")
          {
            std::stringstream POTstrFilled, POTstrEmpty;
            POTstrFilled << std::fixed << std::scientific << std::setprecision(2) << waterFilledPOTData;
            POTstrEmpty << std::fixed << std::scientific << std::setprecision(2) << waterEmptyPOTData;
            std::cout << "waterFilledPOTData: "<< waterFilledPOTData <<std::endl;
            std::cout << "waterEmptyPOTData: "<< waterEmptyPOTData <<std::endl;
            std::string potstrFilled = std::string("Water Filled POT: ")+POTstrFilled.str();
            std::string potstrEmpty = std::string("Water Empty POT ")+POTstrEmpty.str();// + " (Scaled to match filled)";
            plotter.AddPlotLabel(potstrFilled.c_str() , 0.18, 0.99, 0.025, 1, 42, 13, 0);
            plotter.AddPlotLabel(potstrEmpty.c_str() , 0.18, 0.96, 0.025, 1, 42, 13, 0);
          }
          else
          {
            std::stringstream POTstr;
            POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOT;
            std::string potstr = std::string("POT Used: ")+POTstr.str();
            plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.97, 0.025, 1, 42, 13, 0);
          }

          plotter.WritePreliminary( 0.75, 0.075, 0.035, true);
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 35, 1, 43, 22, 0);
          //plotter.AddPlotLabel("Data/MC" , 0.08, 0.3, 0.03, 1, 42, 33, 90);

          plotter.AddPlotLabel("Events#times10^{3} /(Gev/c)^2", 0.020 ,0.51, 30, 1, 43, 22, 90);
          plotter.AddPlotLabel("Muon Longitudinal Momentum (Gev/c)", 0.51, 0.025, 30, 1, 43);
          canvas->Print((tgtname+"_"+prefix+"_EvRateLongitudinal.png").c_str());
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
            xaxislabel = "Reconstructed p_{T, #mu} [GeV/c]";
            migration->GetXaxis()->SetRange(1, 13);
            migration->GetYaxis()->SetRange(1, 13);
          }
          else if (prefix=="pZmu")
          {
            tmp = "Muon p_{Z}";
            yaxislabel = "True p_{Z, #mu} [GeV/c]";
            xaxislabel = "Reconstructed p_{Z, #mu} [GeV/c]";
            migration->GetXaxis()->SetRange(1, 13);
            migration->GetYaxis()->SetRange(1, 13);
          }
          else if (prefix=="BjorkenX")
          {
            tmp = "Bjorken X";
            yaxislabel = "True X";
            xaxislabel = "Reconstructed X";
            can.SetLogx();
            can.SetLogy();
            migration->GetXaxis()->SetRange(1, 6);
            migration->GetYaxis()->SetRange(1, 6);
            //tempevrate->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix=="Erecoil")
          {
            tmp = "E_{recoil}";
            yaxislabel = "Truth E_{recoil} [GeV]";
            xaxislabel = "Reconstructed E_{recoil} [GeV]";
          } 
          else if (prefix=="Emu")
          {
            tmp = "E_{#mu}";
            yaxislabel = "True E_{#mu} [GeV]";
            xaxislabel = "Reconstructed E_{#mu} [GeV]";
            migration->GetXaxis()->SetRange(1, 11);
            migration->GetYaxis()->SetRange(1, 11);
          } 

          //Int_t colors[] = {0, 1};
          //gStyle->SetPalette(2, colors);
          plotter.axis_title_offset_z = 1.25;
          plotter.DrawNormalizedMigrationHistogram(migration, false, false, true, false, 0.02);
          plotter.WritePreliminary( 0.65, 0.2, 0.035, true);
          std::string title = tgtString + " " + tmp;
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 25, 1, 53, 22, 0);
          plotter.AddPlotLabel(yaxislabel.c_str() , 0.08, 0.9, 0.03, 1, 62, 33, 90);
          plotter.AddPlotLabel(xaxislabel.c_str() , 0.9, 0.08, 0.03, 1, 62, 33, 0);

          /*   void DrawNormalizedMigrationHistogram(
            const TH2D* h_migration,
            const bool drawAsMatrix = false,
            const bool coarseContours = false,
            const bool includeFlows = true,
            const bool no_text = false
            ); */
          //plotter.AddHistoTitle(title.c_str(), 0.039);
          can.Update();
          //can.Print((tgtname+"_"+prefix+"_DSSidebandScaled.pdf").c_str());
          can.Print((tgtname+"_"+prefix+"_migration.png").c_str());
          can.SetLogx(0);
          can.SetLogy(0);
        }





















        //Daisy Comparisons


        std::string daisyFilepath = "./"+prefix+"_Daisy_crossSection.root";
        std::string trackerNoDaisyFilepath = "./"+prefix+"_crossSection.root";
        if (std::filesystem::exists(daisyFilepath))
        {
          auto trackerFile = TFile::Open(daisyFilepath.c_str());
          if (!trackerFile)
          {
            std::cerr << "Failed to open tracker file\n";
            return 2;
          }
          auto trackerFileNoDaisy = TFile::Open(trackerNoDaisyFilepath.c_str());
          if (!trackerFileNoDaisy)
          {
            std::cerr << "Failed to open tracker file\n";
            return 2;
          }
          std::cout<<"Tracker extracted cross sections found - plotting material/CH ratios\n";
          std::string titlebase;
          PlotUtils::MnvH2D* plasticxsec;
          if (tgtname == "Iron" || tgtname == "2026" || tgtname == "3026")
          {
            titlebase = "Fe";
            plasticxsec=util::GetIngredient<PlotUtils::MnvH2D>(*trackerFile, prefix+"_Fe_CrossSection");
          }
          if (tgtname == "Lead")
          {
            titlebase="Pb";
            plasticxsec=util::GetIngredient<PlotUtils::MnvH2D>(*trackerFile, prefix+"_Pb_CrossSection");
          }
          if (tgtname == "Carbon")
          {
            titlebase="C";
            plasticxsec=util::GetIngredient<PlotUtils::MnvH2D>(*trackerFile, prefix+"_C_CrossSection");
          }
           if (tgtname == "Water")
          {
            titlebase="H_{2}O";
            plasticxsec=util::GetIngredient<PlotUtils::MnvH2D>(*trackerFileNoDaisy, "crossSection");
          } 
          titlebase+="/CH Ratio in "+prefix;
          






          











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
            PlotUtils::MnvH2D* tempplastic = plasticxsec->Clone();
            PlotUtils::MnvH2D* tempxsec = crossSectionTuned->Clone();
            //TH2D* tempxsec = tempsimxsec->Clone();

            int nbinsY = tempxsec->GetNbinsY();

            int nbinsX = tempxsec->GetNbinsX();



            PlotUtils::MnvH1D* hist;
            PlotUtils::MnvH1D* histplastic;

            //Autos-scaling
            double maxVal = 0 ;
            for (int i = 0; i < nbinsY; i++){
                TH1D* tempdataproj = tempxsec->ProjectionX("a", i+1, i+1);
                tempdataproj->Scale(1);
                TH1D* tempmcproj = tempplastic->ProjectionX("b", i+1, i+1);
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
                histplastic = tempplastic->ProjectionX("d", i+1, i+1);
                TH1D* inverseProjection = tempxsec->ProjectionY("py", 1, 1); //used for auto-detecting subplot limits

                std::cout<<"Data/MC pot << " << dataPOT/mcPOT <<std::endl;
                Int_t datamaxbin = hist->GetMaximumBin();
                double datamax = hist->GetBinContent(datamaxbin);
                Int_t mcmaxbin = histplastic->GetMaximumBin();
                double mcmax = histplastic->GetBinContent(mcmaxbin);
                std::cout<<"data max " << datamax << std::endl;
                std::cout<<"mc max " << mcmax << std::endl;
                std::cout<<"mc max bin" << mcmaxbin << std::endl;
                double overallmax = (datamax > mcmax) ? datamax : mcmax;
                double scalefactor = (std::round((maxVal*5) / overallmax))/5; //round to nearest 0.2
                std::cout<<"scale factor " << scalefactor << std::endl;
                std::cout<<"max " << overallmax*scalefactor << std::endl;

                hist->Scale(scalefactor);
                std::cout<<"Here\n";
                histplastic->Scale(scalefactor);
                gPad->Update();

                plotter.mc_error_color=kRed;
                plotter.mc_error_style=3002;
                hist->GetXaxis()->SetRangeUser(0, 2.5);
                histplastic->GetXaxis()->SetRangeUser(0, 2.5);
                plotter.DrawDataMCRatio(hist, histplastic, 1.0, true, true, -1, -1, "");
                        
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
            
            std::string title = tgtString;

            if (tgtname == "Water")
            {
              std::stringstream POTstrFilled, POTstrEmpty;
              POTstrFilled << std::fixed << std::scientific << std::setprecision(2) << waterFilledPOTData;
              POTstrEmpty << std::fixed << std::scientific << std::setprecision(2) << waterEmptyPOTData;
              std::cout << "waterFilledPOTData: "<< waterFilledPOTData <<std::endl;
              std::cout << "waterEmptyPOTData: "<< waterEmptyPOTData <<std::endl;
              std::string potstrFilled = std::string("Water Filled POT ")+POTstrFilled.str();
              std::string potstrEmpty = std::string("Water Empty POT ")+POTstrEmpty.str();// + " (Scaled to match filled)";
              plotter.AddPlotLabel(potstrFilled.c_str() , 0.18, 0.99, 0.025, 1, 42, 13, 0);
              plotter.AddPlotLabel(potstrEmpty.c_str() , 0.18, 0.96, 0.025, 1, 42, 13, 0);
            }
            else
            {
              std::stringstream POTstr;
              POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOT;
              std::string potstr = std::string("POT Used: ")+POTstr.str();
              plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.97, 0.025, 1, 42, 13, 0);
            }
            title+="/CH Ratio";
            plotter.WritePreliminary( 0.75, 0.075, 0.035, true);
            plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 35, 1, 43, 22, 0);
            //plotter.AddPlotLabel("Data/MC" , 0.08, 0.3, 0.03, 1, 42, 33, 90);

            plotter.AddPlotLabel("data/MC", 0.020 ,0.51, 30, 1, 43, 22, 90);
            plotter.AddPlotLabel("Muon Transverse Momentum (Gev/c)", 0.51, 0.025, 30, 1, 43);
            canvas->Print((tgtname+"_"+prefix+"_CHRatio_Transverse.png").c_str());
          }
          /* { //Data MC ratio - Transverse
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
            PlotUtils::MnvH2D* tempxsec = crossSectionTuned->Clone();
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
            
            std::string title = tgtString;

            if (tgtname == "Water")
            {
              std::stringstream POTstrFilled, POTstrEmpty;
              POTstrFilled << std::fixed << std::scientific << std::setprecision(2) << waterFilledPOTData;
              POTstrEmpty << std::fixed << std::scientific << std::setprecision(2) << waterEmptyPOTData;
              std::cout << "waterFilledPOTData: "<< waterFilledPOTData <<std::endl;
              std::cout << "waterEmptyPOTData: "<< waterEmptyPOTData <<std::endl;
              std::string potstrFilled = std::string("Water Filled POT: ")+POTstrFilled.str();
              std::string potstrEmpty = std::string("Water Empty POT ")+POTstrEmpty.str();// + " (Scaled to match filled)";
              plotter.AddPlotLabel(potstrFilled.c_str() , 0.18, 0.99, 0.025, 1, 42, 13, 0);
              plotter.AddPlotLabel(potstrEmpty.c_str() , 0.18, 0.96, 0.025, 1, 42, 13, 0);
            }
            else
            {
              std::stringstream POTstr;
              POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOT;
              std::string potstr = std::string("POT Used: ")+POTstr.str();
              plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.97, 0.025, 1, 42, 13, 0);
            }

            plotter.WritePreliminary( 0.75, 0.075, 0.035, true);
            plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 35, 1, 43, 22, 0);
            //plotter.AddPlotLabel("Data/MC" , 0.08, 0.3, 0.03, 1, 42, 33, 90);

            plotter.AddPlotLabel("d^{2}#sigma/dp_{t}dp_{z} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)", 0.020 ,0.51, 30, 1, 43, 22, 90);
            plotter.AddPlotLabel("Muon Longitudinal Momentum (Gev/c)", 0.51, 0.025, 30, 1, 43);
            canvas->Print((tgtname+"_"+prefix+"_DataMCRatioLongitudinal.png").c_str());
          } */

        }
        else std::cout<<"No tracker extracted cross sections found - Not plotting material/CH ratios\n";




        //--------------------------------------------------------------------------------------
        //End - Plotting
        //--------------------------------------------------------------------------------------
        std::cout<<"Closing\n";
        outFile->Close();
        std::cout<<"Closed\n";
      }
      catch (const std::runtime_error &e)
      {
        std::cerr << "Failed to extract a cross section for prefix " << prefix << " and target " << tgt << " : " << e.what() << "\n";
        return 4;
        // break;
      }
      std::cout<<"Deleting\n";
      std::cout<<"Deleted\n";
    }
  }
  return 0;
}
