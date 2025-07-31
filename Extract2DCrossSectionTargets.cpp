// File: ExtractCrossSection.cpp
// Brief: Given data and MC files from analyses/studies/CrossSection.h, extract a 1D differential cross section.
//        Subtracts backgrounds, performs unfolding, applies efficiency x acceptance correction, and
//        divides by flux and number of nucleons.  Writes a .root file with the cross section histogram.
//
// Usage: ExtractCrossSection <unfolding iterations> <data.root> <mc.root>
//
// Author: Andrew Olivier aolivier@ur.rochester.edu

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
#include "string.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

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
    plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
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



std::vector<std::string> findContainingDirectories(std::string dir, bool recursiveSearch = true, int curdepth = 0, int maxdepth = 2)
{
  std::vector<std::string> dirpaths;
  for (const auto &entry : std::filesystem::directory_iterator(dir))
  {
    std::string path = entry.path();
    //std::cout<<"Investigating path -wtf: " << path <<std::endl;
    std::filesystem::file_type ft = std::filesystem::status(path).type();
    if (ft == std::filesystem::file_type::regular)
    {
      //std::cout<<"Found regular file at: " << path <<std::endl;
      const size_t base = path.find("runEventLoopTargets");
      if (base != std::string::npos)
      {
        //std::cout<<"Found matching file: " << path <<std::endl;
        if (std::find(dirpaths.begin(), dirpaths.end(), dir) == dirpaths.end()) dirpaths.push_back(dir);
      }
    }
    else if ((ft == std::filesystem::file_type::directory || ft == std::filesystem::file_type::symlink) && recursiveSearch && curdepth<=maxdepth)
    {
      //std::cout<<"Found directory file at: " << path <<std::endl;
      std::vector<std::string> subdirs = findContainingDirectories(path, recursiveSearch, curdepth+1);
      for (auto subdirpath : subdirs)
      {
        if (std::find(dirpaths.begin(), dirpaths.end(), subdirpath) == dirpaths.end()) dirpaths.push_back(subdirpath);
      }
    }
  }
  std::sort(dirpaths.begin(), dirpaths.end());
  return dirpaths;
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

  if (argc != 5)
  {
    std::cerr << "Expected 4 arguments, but I got " << argc - 1 << ".\n"
              << "USAGE: ExtractCrossSection <unfolding iterations> <directory> <target>\n"
              << "Where <target> is the name of a target or \"all\" to loop over all targets found\n"
              << "or the name of the material(+'_tuned', if applicable), if looking at combined playlists\n"
              << "e.g: ExtractCrossSection 5 ./ 2026 minervame1A  -- to just extract xsecs for tgt2 iron\n"
              << "e.g: ExtractCrossSection 5 ./ all minervame1A  -- to just extract xsecs for all targets\n"
              << "e.g: ExtractCrossSection 5 ./ iron minervame1A  -- to just extract xsecs for combined iron\n"
              << "e.g: ExtractCrossSection 5 ./ iron_tuned minervame1A  -- to just extract xsecs for combined iron after sideband tune\n";
    return 1;
  }

  const int nIterations = std::stoi(argv[1]);
  std::string indir = std::string(argv[2]);
  std::string intgt = std::string(argv[3]);

  bool recursiveSearch = false;
  for (int c = 0; c<argc; c++) if (std::string(argv[c])=="-r") recursiveSearch = true;


  std::vector<std::string> dirs = findContainingDirectories(indir, recursiveSearch);

  std::vector<std::string> targets;
  if (intgt == "ALL" || intgt == "all")
  {
    targets = {"1026", "1082", "2026", "2082", "3006", "3026", "3082", "4082", "5026", "5082", "6000", "Iron", "Carbon", "Lead"};
  }
  else
  {
    targets = {intgt};
  }
  for (std::string &tgt : targets)
  {
    //if (tgt != "3006") continue; //Used for testing with only one target

    std::vector<std::string> crossSectionPrefixes = {"nuke_pTmu", "nuke_pZmu", "nuke_BjorkenX", "nuke_Erecoil", "nuke_Emu"};

    /* for (auto key : *dataFile->GetListOfKeys())
    {
      const std::string keyName = key->GetName();
      if (keyName == "POTUsed")
        continue;
      std::cout << "keyName " << keyName << std::endl;
      const size_t endOfPrefix = keyName.find("_data");
      std::string prefix = keyName.substr(0, endOfPrefix);
      std::cout << "prefix " << prefix << std::endl;
      bool twoDimension = (keyName == "nuke_pTmu_pZmu_data");

      bool alreadyInVector = std::find(crossSectionPrefixes.begin(), crossSectionPrefixes.end(), prefix) != crossSectionPrefixes.end();
      std::cout << "twoDimension " << twoDimension << std::endl;
      std::cout << "alreadyInVector " << alreadyInVector << std::endl;
      if (endOfPrefix != std::string::npos && !alreadyInVector && twoDimension)
        crossSectionPrefixes.push_back(prefix);
    } */
    for (const auto &prefix : crossSectionPrefixes)
    {
      //if (!(prefix == "nuke_Erecoil" || prefix == "nuke_pTmu")) continue; //Used for testing with only subset of prefixes
      std::cout << "Current working on prefix: " << prefix << std::endl;
      try
      {
        PlotUtils::MnvH2D* flux;
        PlotUtils::MnvH2D* folded;
        PlotUtils::MnvH2D* migration;
        PlotUtils::MnvH2D* migration_reco;
        PlotUtils::MnvH2D* migration_truth;
        PlotUtils::MnvH2D* effNum;
        PlotUtils::MnvH2D* effDenom;
        PlotUtils::MnvH2D* effDenom2P2H;
        PlotUtils::MnvH2D* effDenomDIS;
        PlotUtils::MnvH2D* effDenomRES;
        PlotUtils::MnvH2D* effDenomQE;
        PlotUtils::MnvH2D* effDenomOther;

        PlotUtils::MnvH1D* USSidebandSignal;
        PlotUtils::MnvH1D* USSidebandDS;
        PlotUtils::MnvH1D* USSidebandUS;
        PlotUtils::MnvH1D* USSidebandOther;
        PlotUtils::MnvH1D* DSSidebandSignal;
        PlotUtils::MnvH1D* DSSidebandDS;
        PlotUtils::MnvH1D* DSSidebandUS;
        PlotUtils::MnvH1D* DSSidebandOther;
        PlotUtils::MnvH1D* DataUSSideband;
        PlotUtils::MnvH1D* DataDSSideband;
        PlotUtils::MnvH1D* DataSignal;

        PlotUtils::MnvH2D* BackgroundUSPlastic;
        PlotUtils::MnvH2D* BackgroundDSPlastic;
        PlotUtils::MnvH2D* BackgroundWrongSign;
        PlotUtils::MnvH2D* BackgroundNC;
        PlotUtils::MnvH2D* BackgroundOther;

        PlotUtils::MnvH2D *fluxIntReweighted;

        //Flux parameters
        int n_flux_universes = 100; // Is this right
        const bool use_nue_constraint = true;
        const std::string project_dir = "targets_2345_jointNueIMD";
        double min_energy = 0;
        double max_energy = 100;

        double mcPOT = 0;
        double dataPOT = 0;

        std::vector <std::string> targetsInTgt;
        if (tgt=="Iron") targetsInTgt = {"2026", "3026", "5026"};
        else if (tgt=="Carbon") targetsInTgt = {"3006"};
        else if (tgt=="Lead") targetsInTgt = {"2082", "3082", "4082", "5082"};
        else targetsInTgt = {tgt};
    
        {

          std::cout<<"Adding target: " << targetsInTgt[0] << " from directory: " << dirs[0] << std::endl;
          std::string datapath = dirs[0] + "/runEventLoopTargetsData" + targetsInTgt[0] + ".root";
          std::string mcpath = dirs[0] + "/runEventLoopTargetsMC" + targetsInTgt[0] + ".root";
          std::string migrationpath = dirs[0] + "/runEventLoopTargets2DMigration" + targetsInTgt[0] + ".root";

          auto dataFile = TFile::Open(datapath.c_str(), "READ");
          if (!dataFile)
          {
            std::cerr << "Failed to open data file " << datapath.c_str() << ".\n";
            return 2;
          }

          auto mcFile = TFile::Open(mcpath.c_str(), "READ");
          if (!mcFile)
          {
            std::cerr << "Failed to open MC file " << mcpath.c_str() << ".\n";
            return 3;
          }

          auto migrationFile = TFile::Open(migrationpath.c_str(), "READ");
          if (!migrationFile)
          {
            std::cerr << "Failed to open MC file " << migrationpath.c_str() << ".\n";
            return 3;
          }

          folded = util::GetIngredient<PlotUtils::MnvH2D>(*dataFile, (std::string("data")), prefix);
          migration = util::GetIngredient<PlotUtils::MnvH2D>(*migrationFile, (std::string("migration")), prefix);
          migration_reco = util::GetIngredient<PlotUtils::MnvH2D>(*migrationFile, (std::string("reco")), prefix);
          migration_truth = util::GetIngredient<PlotUtils::MnvH2D>(*migrationFile, (std::string("truth")), prefix);

          effNum = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_numerator")), prefix);
          effDenom = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator")), prefix);
          effDenom2P2H = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_2p2h")), prefix);
          effDenomDIS = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_DIS")), prefix);
          effDenomRES = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_RES")), prefix);
          effDenomQE = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_QE")), prefix);
          effDenomOther = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_Other")), prefix);

          USSidebandSignal = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_US_sideband_Signal")));
          USSidebandDS = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_US_sideband_DS")));
          USSidebandUS = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_US_sideband_US")));
          USSidebandOther = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_US_sideband_Other")));
          DSSidebandSignal = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_DS_sideband_Signal")));
          DSSidebandDS = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_DS_sideband_DS")));
          DSSidebandUS = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_DS_sideband_US")));
          DSSidebandOther = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_DS_sideband_Other")));
          DataUSSideband = util::GetIngredient<PlotUtils::MnvH1D>(*dataFile, (std::string("nuke_segment_US_Sideband")));
          DataDSSideband = util::GetIngredient<PlotUtils::MnvH1D>(*dataFile, (std::string("nuke_segment_DS_Sideband")));
          DataSignal = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_data")));
          DataSignal = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_data")));

          BackgroundUSPlastic = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("by_BKG_Label_US_Plastic")), prefix);
          BackgroundDSPlastic = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("by_BKG_Label_DS_Plastic")), prefix);
          BackgroundWrongSign = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("by_BKG_Label_Wrong_Sign_Bkg")), prefix);
          BackgroundNC = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("by_BKG_Label_NC_Bkg")), prefix);
          BackgroundOther = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("by_BKG_Label_Other")), prefix);
          
          
          mcPOT += util::GetIngredient<TParameter<double>>(*mcFile, "POTUsed")->GetVal();
          std::cout<<"mcPOT: " << mcPOT << std::endl;
          double tempDataPOT = util::GetIngredient<TParameter<double>>(*dataFile, "POTUsed")->GetVal();
          dataPOT+=tempDataPOT;

          std::string playlistUsed = util::GetIngredient<TNamed>(*mcFile, "PlaylistUsed")->GetTitle();
          int nuoranu = util::nuOrAntiNuMode(playlistUsed);
          int nupdg;
          if (nuoranu==1) nupdg = 14;
          else if (nuoranu==2) nupdg = -14;

          PlotUtils::FluxReweighter frw = PlotUtils::FluxReweighter( nupdg, use_nue_constraint, playlistUsed, PlotUtils::FluxReweighter::gen2thin, PlotUtils::FluxReweighter::g4numiv6, n_flux_universes );
          fluxIntReweighted = frw.GetIntegratedFluxReweighted(nupdg, effDenom, min_energy, max_energy, true);
          fluxIntReweighted->Scale(tempDataPOT);


          dataFile->Close();
          mcFile->Close();
        }
        for (int c = 0; c<dirs.size(); c++)
        {
          for (int d = 1; d<targetsInTgt.size(); d++)
          {
            std::cout<<"Adding target: " << targetsInTgt[d] << " from directory: " << dirs[c] << std::endl;
            if (c==0 && d==0) continue; //Already covered by the above
            std::string datapath = dirs[0] + "/runEventLoopTargetsData" + targetsInTgt[d] + ".root";
            std::string mcpath = dirs[0] + "/runEventLoopTargetsMC" + targetsInTgt[d] + ".root";
            std::string migrationpath = dirs[0] + "/runEventLoopTargets2DMigration" + targetsInTgt[d] + ".root";

            auto dataFile = TFile::Open(datapath.c_str(), "READ");
            if (!dataFile)
            {
              std::cerr << "Failed to open data file " << datapath.c_str() << ".\n";
              return 2;
            }

            auto mcFile = TFile::Open(mcpath.c_str(), "READ");
            if (!mcFile)
            {
              std::cerr << "Failed to open MC file " << mcpath.c_str() << ".\n";
              return 3;
            }
            
            auto migrationFile = TFile::Open(migrationpath.c_str(), "READ");
            if (!migrationFile)
            {
              std::cerr << "Failed to open MC file " << migrationpath.c_str() << ".\n";
              return 3;
            }
            
            folded->Add(util::GetIngredient<PlotUtils::MnvH2D>(*dataFile, (std::string("data")), prefix));
            migration->Add(util::GetIngredient<PlotUtils::MnvH2D>(*migrationFile, (std::string("migration")), prefix));
            migration_reco->Add(util::GetIngredient<PlotUtils::MnvH2D>(*migrationFile, (std::string("reco")), prefix));
            migration_truth->Add(util::GetIngredient<PlotUtils::MnvH2D>(*migrationFile, (std::string("truth")), prefix));

            effNum->Add(util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_numerator")), prefix));
            effDenom->Add(util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator")), prefix));
            effDenom2P2H->Add(util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_2p2h")), prefix));
            effDenomDIS->Add(util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_DIS")), prefix));
            effDenomRES->Add(util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_RES")), prefix));
            effDenomQE->Add(util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_QE")), prefix));
            effDenomOther->Add(util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_Other")), prefix));

            USSidebandSignal->Add(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_US_sideband_Signal"))));
            USSidebandDS->Add(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_US_sideband_DS"))));
            USSidebandUS->Add(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_US_sideband_US"))));
            USSidebandOther->Add(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_US_sideband_Other"))));
            DSSidebandSignal->Add(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_DS_sideband_Signal"))));
            DSSidebandDS->Add(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_DS_sideband_DS"))));
            DSSidebandUS->Add(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_DS_sideband_US"))));
            DSSidebandOther->Add(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_DS_sideband_Other"))));
            DataUSSideband->Add(util::GetIngredient<PlotUtils::MnvH1D>(*dataFile, (std::string("nuke_segment_US_Sideband"))));
            DataDSSideband->Add(util::GetIngredient<PlotUtils::MnvH1D>(*dataFile, (std::string("nuke_segment_DS_Sideband"))));
            DataSignal->Add(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_data"))));
            DataSignal->Add(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_data"))));

            BackgroundUSPlastic->Add(util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("by_BKG_Label_US_Plastic")), prefix));
            BackgroundDSPlastic->Add(util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("by_BKG_Label_DS_Plastic")), prefix));
            BackgroundWrongSign->Add(util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("by_BKG_Label_Wrong_Sign_Bkg")), prefix));
            BackgroundNC->Add(util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("by_BKG_Label_NC_Bkg")), prefix));
            BackgroundOther->Add(util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("by_BKG_Label_Other")), prefix));

            mcPOT += util::GetIngredient<TParameter<double>>(*mcFile, "POTUsed")->GetVal();
            std::cout<<"mcPOT: " << mcPOT << std::endl;
            double tempDataPOT = util::GetIngredient<TParameter<double>>(*dataFile, "POTUsed")->GetVal();
            dataPOT+=tempDataPOT;

            std::string playlistUsed = util::GetIngredient<TNamed>(*mcFile, "PlaylistUsed")->GetTitle();
            int nuoranu = util::nuOrAntiNuMode(playlistUsed);
            int nupdg;
            if (nuoranu==1) nupdg = 14;
            else if (nuoranu==2) nupdg = -14;

            PlotUtils::FluxReweighter frw = PlotUtils::FluxReweighter( nupdg, use_nue_constraint, playlistUsed, PlotUtils::FluxReweighter::gen2thin, PlotUtils::FluxReweighter::g4numiv6, n_flux_universes );
            fluxIntReweighted = frw.GetIntegratedFluxReweighted(nupdg, effDenom, min_energy, max_energy, true);
            fluxIntReweighted->Scale(tempDataPOT);


            if (d==0) //Only get POT once per playlist. E.g we don't wanna double count the POT for the same runs just because we're looking at target 2026 and 2082
            {
              mcPOT += util::GetIngredient<TParameter<double>>(*mcFile, "POTUsed")->GetVal();
              std::cout<<"mcPOT: " << mcPOT << std::endl;
              double tempDataPOT = util::GetIngredient<TParameter<double>>(*dataFile, "POTUsed")->GetVal();
              dataPOT+=tempDataPOT;


              std::string playlistUsed = util::GetIngredient<TNamed>(*mcFile, "PlaylistUsed")->GetTitle();
              int nuoranu = util::nuOrAntiNuMode(playlistUsed);
              int nupdg;
              if (nuoranu==1) nupdg = 14;
              else if (nuoranu==2) nupdg = -14;

              PlotUtils::FluxReweighter frw = PlotUtils::FluxReweighter( nupdg, use_nue_constraint, playlistUsed, PlotUtils::FluxReweighter::gen2thin, PlotUtils::FluxReweighter::g4numiv6, n_flux_universes );
              fluxIntReweighted->Add(frw.GetIntegratedFluxReweighted(nupdg, effDenom, min_energy, max_energy, true),tempDataPOT);
            }

            dataFile->Close();
            mcFile->Close();
          }
        }

        fluxIntReweighted->Scale(1/dataPOT);


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
            std::cout<<"Made a Minimizer" <<std::endl;
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
            std::cout<<"Made a Minimizer" <<std::endl;
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

        auto simEventRate = effDenom->Clone(); // Make a copy for later
        auto simEventRate2P2H = effDenom2P2H->Clone(); // Make a copy for later
        auto simEventRateDIS = effDenomDIS->Clone(); // Make a copy for later
        auto simEventRateRES = effDenomRES->Clone(); // Make a copy for later
        auto simEventRateQE = effDenomQE->Clone(); // Make a copy for later
        auto simEventRateOther = effDenomOther->Clone(); // Make a copy for later
        // There are no error bands in the data, but I need somewhere to put error bands on the results I derive from it.
        folded->AddMissingErrorBandsAndFillWithCV(*migration);

        // Basing my unfolding procedure for a differential cross section on Alex's MINERvA 101 talk at https://minerva-docdb.fnal.gov/cgi-bin/private/RetrieveFile?docid=27438&filename=whatsACrossSection.pdf&version=1

        std::vector<PlotUtils::MnvH2D*> backgrounds = {BackgroundUSPlastic, BackgroundDSPlastic, BackgroundWrongSign, BackgroundNC, BackgroundOther};
        PlotUtils::MnvH2D* BackgroundUSPlasticScaled = BackgroundUSPlastic->Clone();
        BackgroundUSPlasticScaled->Scale(USScaleFactor);
        PlotUtils::MnvH2D* BackgroundDSPlasticScaled = BackgroundDSPlastic->Clone();
        BackgroundDSPlasticScaled->Scale(DSScaleFactor);
        std::vector<PlotUtils::MnvH2D*> backgroundsScaled = {BackgroundUSPlasticScaled, BackgroundDSPlasticScaled, BackgroundWrongSign, BackgroundNC, BackgroundOther};
        auto toSubtract = std::accumulate(std::next(backgrounds.begin()), backgrounds.end(), (*backgrounds.begin())->Clone(),
                                          [](auto sum, const auto hist)
                                          {
                                            sum->Add(hist);
                                            return sum;
                                          });
        Plot(*toSubtract, "BackgroundSum", prefix, tgt);

        auto bkgSubtracted = std::accumulate(backgrounds.begin(), backgrounds.end(), folded->Clone(),
                                              [mcPOT, dataPOT](auto sum, const auto hist)
                                              {
                                                std::cout << "Subtracting " << hist->GetName() << " scaled by " << -dataPOT / mcPOT << " from " << sum->GetName() << "\n";
                                                sum->Add(hist, -dataPOT / mcPOT);
                                                return sum;
                                              });

        auto bkgScaledSubtracted = std::accumulate(backgroundsScaled.begin(), backgroundsScaled.end(), folded->Clone(),
                                              [mcPOT, dataPOT](auto sum, const auto hist)
                                              {
                                                std::cout << "Subtracting " << hist->GetName() << " scaled by " << -dataPOT / mcPOT << " from " << sum->GetName() << "\n";
                                                sum->Add(hist, -dataPOT / mcPOT);
                                                return sum;
                                              });

        Plot(*bkgSubtracted, "backgroundSubtracted", prefix, tgt);
        Plot(*bkgScaledSubtracted, "backgroundSubtractedScaled", prefix, tgt);

        auto outFile = TFile::Open((tgt + prefix + "_crossSection.root").c_str(), "RECREATE");
        if (!outFile)
        {
          std::cerr << "Could not create a file called " << prefix + "_crossSection.root" << ".  Does it already exist?\n";
          return 5;
        }

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

        USSidebandSignal->Scale(USScaleFactor);
        USSidebandUS->Scale(USScaleFactor);
        USSidebandDS->Scale(USScaleFactor);
        USSidebandOther->Scale(USScaleFactor);

        DSSidebandSignal->Scale(DSScaleFactor);
        DSSidebandUS->Scale(DSScaleFactor);
        DSSidebandDS->Scale(DSScaleFactor);
        DSSidebandOther->Scale(DSScaleFactor);
        
        USSidebandSignal->Write("USSidebandSignalScaled");
        USSidebandUS->Write("USSidebandUSScaled");
        USSidebandDS->Write("USSidebandDSScaled");
        USSidebandOther->Write("USSidebandOtherScaled");
        
        DSSidebandSignal->Write("DSSidebandSignalScaled");
        DSSidebandUS->Write("DSSidebandUSScaled");
        DSSidebandDS->Write("DSSidebandDSScaled");
        DSSidebandOther->Write("DSSidebandOtherScaled");

        bkgSubtracted->Write("backgroundSubtracted");
        bkgScaledSubtracted->Write("backgroundSubtractedSidebandScaled");

        // d'Aogstini unfolding
        auto unfolded = UnfoldHist(bkgSubtracted, migration, migration_reco, migration_truth, nIterations);
        if (!unfolded)
          throw std::runtime_error(std::string("Failed to unfold ") + folded->GetName() + " using " + migration->GetName());
        Plot(*unfolded, "unfolded", prefix, tgt);
        unfolded->Clone()->Write("unfolded"); // TODO: Seg fault first appears when I uncomment this line
        std::cout << "Survived writing the unfolded histogram.\n"
                  << std::flush; // This is evidence that the problem is on the final file Write() and not unfolded->Clone()->Write().
        auto unfolded_tuned = UnfoldHist(bkgScaledSubtracted, migration, migration_reco, migration_truth, nIterations);
        if (!unfolded_tuned)
          throw std::runtime_error(std::string("Failed to unfold ") + folded->GetName() + " using " + migration->GetName());
        Plot(*unfolded_tuned, "unfolded_sidebandTuned", prefix, tgt);
        unfolded_tuned->Clone()->Write("unfolded"); // TODO: Seg fault first appears when I uncomment this line



        effNum->Divide(effNum, effDenom); // Only the 2 parameter version of MnvH1D::Divide()
                                          // handles systematics correctly.
        Plot(*effNum, "efficiency", prefix, tgt);

        unfolded->Divide(unfolded, effNum);
        Plot(*unfolded, "efficiencyCorrected", prefix, tgt);

        unfolded_tuned->Divide(unfolded_tuned, effNum);
        Plot(*unfolded_tuned, "efficiencyCorrected_sidebandTuned", prefix, tgt);


        double nnucleons = 0;
        double nnucleonsData = 0;


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
        else
        {
          int tgtCode = std::stoi(tgtname);
          int tgtMat = tgtCode%1000;
          int tgtNum = (tgtCode-tgtMat)/1000;
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

        PlotUtils::MnvH1D *flux2;
        PlotUtils::MnvH1D *fluxIntegral;
        PlotUtils::MnvH1D *fluxRebinned;

        //std::cout<<"nu_pdg: " << nu_pdg << std::endl;
        std::cout<<"use_nue_constraint: " << use_nue_constraint << std::endl;
        std::cout<<"n_flux_universes: " << n_flux_universes << std::endl;
        std::cout<<"min_energy: " << min_energy << std::endl;
        std::cout<<"max_energy: " << max_energy << std::endl;
        std::cout<<"nnucleons: " << nnucleons << std::endl;
        std::cout<<"nnucleonsData: " << nnucleonsData << std::endl;
        std::cout<<"mcPOT: " << mcPOT << std::endl;
        std::cout<<"dataPOT: " << dataPOT << std::endl;
        //PlotUtils::FluxReweighter frw = PlotUtils::flux_reweighter(platlistname, nu_pdg, use_nue_constraint, n_flux_universes);
        //PlotUtils::FluxReweighter *frw = new PlotUtils::FluxReweighter( nu_pdg, use_nue_constraint, platlistname, PlotUtils::FluxReweighter::gen2thin, PlotUtils::FluxReweighter::g4numiv6, n_flux_universes );
        //std::cout << "ABC123-1 " << std::endl;
        //PlotUtils::MnvH2D *fluxIntReweighted = frw->GetIntegratedFluxReweighted(nu_pdg, simEventRate, min_energy, max_energy, true);
        //fluxIntegral = frw->GetIntegratedTargetFlux(nu_pdg, material, simEventRate, min_energy, max_energy, project_dir);
        //auto &frw2 = PlotUtils::flux_reweighter("minervame1A", nu_pdg, use_nue_constraint, n_flux_universes);
        //flux2 = frw2.GetTargetFluxMnvH1D(nu_pdg, material, project_dir);
        std::cout<<"1234.1\n"; 
        //fluxRebinned = frw->GetRebinnedFluxReweighted_FromInputFlux(flux2, simEventRate); // issue here
        //PlotUtils::MnvH1D *Integrated_fluxGenie = frw->GetIntegratedFluxReweighted_FromInputFlux(flux2, simEventRate, min_energy, max_energy);
        //PlotUtils::MnvH1D *Integrated_fluxGenie2 = frw->GetIntegratedFluxReweighted_FromInputFlux(flux, simEventRate, min_energy, max_energy);
        outFile->cd();
        std::cout<<"1234.2\n"; 
        unfolded->Write("unfolded2");
        std::cout<<"1234.3\n"; 
        fluxIntReweighted->Write("fluxIntReweighted2");
        std::cout<<"1234.4\n"; 
        auto crossSection = normalize(unfolded, fluxIntReweighted, nnucleonsData, dataPOT);
        std::cout<<"1234.5\n"; 
        Plot(*crossSection, "crossSection", prefix, tgt);
        std::cout<<"1234.6\n"; 
        crossSection->Write("crossSection");
        std::cout<<"1234.7\n"; 
        auto crossSectionTuned = normalize(unfolded_tuned, fluxIntReweighted, nnucleonsData, dataPOT);
        //Plot(*crossSectionTuned, "crossSectionTuned", prefix, tgt);
        std::cout<<"1234.8\n"; 
        crossSectionTuned->Write("crossSectionTuned");
        std::cout<<"1234.9\n"; 
        simEventRate->Write("simulatedEventRate");
        std::cout<<"1234.10\n"; 
        fluxIntReweighted->Write("fluxIntReweighted");
        std::cout<<"1234.11\n"; 
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
        std::cout<<"1234.12\n"; 
        //Plot(*simulatedCrossSection, "simulatedCrossSection", prefix, tgt);
        simulatedCrossSection->Write("simulatedCrossSection");
        std::cout<<"1234.13\n";
        auto simulatedCrossSection2P2H = normalize(simEventRate2P2H, fluxIntReweighted, nnucleons, mcPOT);
        std::cout<<"1234.14\n";
        simulatedCrossSection2P2H->Write("simulatedCrossSection2P2H");
        std::cout<<"1234.15\n";
        auto simulatedCrossSectionDIS = normalize(simEventRateDIS, fluxIntReweighted, nnucleons, mcPOT);
        std::cout<<"1234.16\n";
        simulatedCrossSectionDIS->Write("simulatedCrossSectionDIS");
        std::cout<<"1234.17\n";
        auto simulatedCrossSectionRES = normalize(simEventRateRES, fluxIntReweighted, nnucleons, mcPOT);
        simulatedCrossSectionRES->Write("simulatedCrossSectionRES");
        auto simulatedCrossSectionQE = normalize(simEventRateQE, fluxIntReweighted, nnucleons, mcPOT);
        simulatedCrossSectionQE->Write("simulatedCrossSectionQE");
        auto simulatedCrossSectionOther = normalize(simEventRateOther, fluxIntReweighted, nnucleons, mcPOT);
        simulatedCrossSectionOther->Write("simulatedCrossSectionOther");
        std::cout<<"1234.11\n"; 
        outFile->Close();
        // return 0;
      }
      catch (const std::runtime_error &e)
      {
        std::cerr << "Failed to extract a cross section for prefix " << prefix << " and target " << tgt << " : " << e.what() << "\n";
        return 4;
        // break;
      }
    }
  }
  return 0;
}
