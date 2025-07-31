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
void Plot(PlotUtils::MnvH1D &hist, const std::string &stepName, const std::string &prefix, const std::string &target)
{
  return;
  bool plotpngs = true;
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
  if (plotpngs) can.Print((target + "_" + prefix + "_" + stepName + "_otherUncertainties.png").c_str());
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
    //std::vector<std::string> crossSectionPrefixes = {"nuke_pTmu", "nuke_pZmu", "nuke_BjorkenX", "nuke_Erecoil", "nuke_Emu"};
    std::vector<std::string> crossSectionPrefixes = {"nuke_pTmu", "nuke_pZmu", "nuke_BjorkenX", "nuke_Erecoil", "nuke_Emu"};

    for (const auto &prefix : crossSectionPrefixes)
    {
      //if (!(prefix == "nuke_Erecoil" || prefix == "nuke_pTmu")) continue; //Used for testing with only subset of prefixes
      std::cout << "Currently working on prefix: " << prefix << std::endl;
      try
      {
        int waterFilledEmpty = 0; //Only used with water target
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
        else targetsInTgt = {tgt};

        //Flux parameters
        int n_flux_universes = 100; // Is this right
        const bool use_nue_constraint = true;
        const std::string project_dir = "targets_2345_jointNueIMD";
        double min_energy = 0;
        double max_energy = 100;
        
        PlotUtils::MnvH1D *fluxIntReweighted = new PlotUtils::MnvH1D();

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

        PlotUtils::MnvH1D* BackgroundUSPlastic = new PlotUtils::MnvH1D();
        PlotUtils::MnvH1D* BackgroundDSPlastic = new PlotUtils::MnvH1D();
        PlotUtils::MnvH1D* BackgroundWrongSign = new PlotUtils::MnvH1D();
        PlotUtils::MnvH1D* BackgroundNC = new PlotUtils::MnvH1D();
        PlotUtils::MnvH1D* BackgroundOther = new PlotUtils::MnvH1D();

        double mcPOT = 0;
        double dataPOT = 0;

        for (int c = 0; c<dirs.size(); c++)
        {
          std::cout<<"Investigating directory " << dirs[c] << " which is " << c+1 <<"/"<<dirs.size() <<" playlists identified" <<std::endl;
          for (int d = 0; d<targetsInTgt.size(); d++)
          {
            std::cout<<"Investigating target " << targetsInTgt[d] << " which is " << d+1 <<"/"<<targetsInTgt.size() <<" targets selected" <<std::endl;
            std::string datapath = dirs[c] + "/runEventLoopTargetsData" + targetsInTgt[d] + ".root";
            std::string mcpath = dirs[c] + "/runEventLoopTargetsMC" + targetsInTgt[d] + ".root";

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
            
            std::string playlistUsed = util::GetIngredient<TNamed>(*mcFile, "PlaylistUsed")->GetTitle();
            int filledorempty = util::filledOrEmptyMEPlaylist(playlistUsed);

            int nuoranu = util::nuOrAntiNuMode(playlistUsed);
            int nupdg;
            if (nuoranu==1) nupdg = 14;
            else if (nuoranu==2) nupdg = -14;

            if (waterFilledEmpty && (waterFilledEmpty != filledorempty)) continue;
            if (pdg!=nupdg)
            {
              std::cout<<"Skipping this set of files because this playlist pdg is " << nupdg << " but you specified " << pdg << std::endl;
              continue;
            }
            util::AddHist(*flux,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("reweightedflux_integrated")), prefix));
            util::AddHist(*folded,util::GetIngredient<PlotUtils::MnvH1D>(*dataFile, (std::string("data")), prefix));
            util::AddHist(*migration,util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("migration")), prefix));
            util::AddHist(*effNum,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_numerator")), prefix));
            util::AddHist(*effDenom,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator")), prefix));
            util::AddHist(*effDenom2P2H,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator_intChannels_2p2h")), prefix));
            util::AddHist(*effDenomDIS,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator_intChannels_DIS")), prefix));
            util::AddHist(*effDenomRES,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator_intChannels_RES")), prefix));
            util::AddHist(*effDenomQE,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator_intChannels_QE")), prefix));
            util::AddHist(*effDenomOther,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator_intChannels_Other")), prefix));

            util::AddHist(*USSidebandSignal,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_US_sideband_Signal"))));
            util::AddHist(*USSidebandDS,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_US_sideband_DS"))));
            util::AddHist(*USSidebandUS,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_US_sideband_US"))));
            util::AddHist(*USSidebandOther,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_US_sideband_Other"))));
            util::AddHist(*DSSidebandSignal,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_DS_sideband_Signal"))));
            util::AddHist(*DSSidebandDS,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_DS_sideband_DS"))));
            util::AddHist(*DSSidebandUS,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_DS_sideband_US"))));
            util::AddHist(*DSSidebandOther,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_DS_sideband_Other"))));
            util::AddHist(*DataUSSideband,util::GetIngredient<PlotUtils::MnvH1D>(*dataFile, (std::string("nuke_segment_US_Sideband"))));
            util::AddHist(*DataDSSideband,util::GetIngredient<PlotUtils::MnvH1D>(*dataFile, (std::string("nuke_segment_DS_Sideband"))));
            util::AddHist(*DataSignal,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("nuke_segment_data"))));

            util::AddHist(*BackgroundUSPlastic,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_US_Plastic")), prefix));
            util::AddHist(*BackgroundDSPlastic,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_DS_Plastic")), prefix));
            util::AddHist(*BackgroundWrongSign,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_Wrong_Sign_Bkg")), prefix));
            util::AddHist(*BackgroundNC,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_NC_Bkg")), prefix));
            util::AddHist(*BackgroundOther,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_Other")), prefix));

            if (d==0) //Only get POT once per playlist. E.g we don't wanna double count the POT for the same runs just because we're looking at target 2026 and 3026
            {
              mcPOT += util::GetIngredient<TParameter<double>>(*mcFile, "POTUsed")->GetVal();
              std::cout<<"mcPOT: " << mcPOT << std::endl;
              double tempDataPOT = util::GetIngredient<TParameter<double>>(*dataFile, "POTUsed")->GetVal();
              dataPOT+=tempDataPOT;

              //Scale the integrated flux for this playlist by playlist data POT for appropriate scaling when applying across different playlists
              PlotUtils::FluxReweighter frw = PlotUtils::FluxReweighter( nupdg, use_nue_constraint, playlistUsed, PlotUtils::FluxReweighter::gen2thin, PlotUtils::FluxReweighter::g4numiv6, n_flux_universes );                    
              auto tempIntFlux = frw.GetIntegratedFluxReweighted(nupdg, effDenom, min_energy, max_energy, true)->Clone();
              tempIntFlux->Scale(tempDataPOT);
              util::AddHist(*fluxIntReweighted,tempIntFlux);
            }

            dataFile->Close();
            mcFile->Close();
          }
        }

        //Normalising the integrated flux across different playlists by dataPOT
        fluxIntReweighted->Scale(1/dataPOT);

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

        std::vector<PlotUtils::MnvH1D*> backgrounds = {BackgroundUSPlastic, BackgroundDSPlastic, BackgroundWrongSign, BackgroundNC, BackgroundOther};
        PlotUtils::MnvH1D* BackgroundUSPlasticScaled = BackgroundUSPlastic->Clone();
        BackgroundUSPlasticScaled->Scale(USScaleFactor);
        PlotUtils::MnvH1D* BackgroundDSPlasticScaled = BackgroundDSPlastic->Clone();
        BackgroundDSPlasticScaled->Scale(DSScaleFactor);
        std::vector<PlotUtils::MnvH1D*> backgroundsScaled = {BackgroundUSPlasticScaled, BackgroundDSPlasticScaled, BackgroundWrongSign, BackgroundNC, BackgroundOther};
        auto toSubtract = std::accumulate(std::next(backgrounds.begin()), backgrounds.end(), (*backgrounds.begin())->Clone(),
                                          [](auto sum, const auto hist)
                                          {
                                            sum->Add(hist);
                                            return sum;
                                          });
        Plot(*toSubtract, "BackgroundSum", prefix, tgt);

        //sum up backgrounds without sideband scaling
        auto bkgSubtracted = std::accumulate(backgrounds.begin(), backgrounds.end(), folded->Clone(),
                                              [mcPOT, dataPOT](auto sum, const auto hist)
                                              {
                                                std::cout << "Subtracting " << hist->GetName() << " scaled by " << -dataPOT / mcPOT << " from " << sum->GetName() << "\n";
                                                sum->Add(hist, -dataPOT / mcPOT);
                                                return sum;
                                              });

        //sum up backgrounds with sideband scaling
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

        bkgSubtracted->Write("backgroundSubtracted");
        bkgScaledSubtracted->Write("backgroundSubtractedSidebandScaled");

        // d'Aogstini unfolding
        auto unfolded = UnfoldHist(bkgSubtracted, migration, nIterations);
        if (!unfolded)
          throw std::runtime_error(std::string("Failed to unfold ") + folded->GetName() + " using " + migration->GetName());
        Plot(*unfolded, "unfolded", prefix, tgt);
        unfolded->Clone()->Write("unfolded"); // TODO: Seg fault first appears when I uncomment this line
        auto unfolded_tuned = UnfoldHist(bkgScaledSubtracted, migration, nIterations);
        if (!unfolded_tuned)
          throw std::runtime_error(std::string("Failed to unfold ") + folded->GetName() + " using " + migration->GetName());
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
        //Plot(*crossSectionTuned, "crossSectionTuned", prefix, tgt);
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
        
        //Plotting stuff - Should move this elsewhere long-term
 
        { //Plotting cross sections
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
          //plotter.axis_maximum = 0.4;
          TObjArray* arr = new TObjArray();
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
          arr->Add(temp2p2h);
          arr->Add(tempdis);
          arr->Add(tempres);
          arr->Add(tempqe);
          arr->Add(tempother);

          std::string tmp; 
          std::string y_title;
          if (prefix=="nuke_pTmu")
          {
            tmp = "Muon p_{T}";
            y_title="d#sigma/dp_{t} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            tempxsec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix=="nuke_pZmu")
          {
            tmp = "Muon p_{Z}";
            y_title="d#sigma/dp_{Z} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            tempxsec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix=="nuke_BjorkenX")
          {
            tmp = "{X}";
            y_title="d#sigma/d{X} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
          }
          else if (prefix=="nuke_Erecoil")
          {
            tmp = "E_{recoil}";
            y_title="d#sigma/d{X} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
          } 
          else if (prefix=="nuke_Emu")
          {
            tmp = "E_{#mu}";
            y_title="d#sigma/dE_{#mu} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            tempxsec->GetXaxis()->SetRange(1, 11);
          } 

          std::string title = tmp + " " + tgt;
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;
          plotter.DrawDataStackedMC(tempxsec, arr, 1, "TR", "Data", 2, 1, 3001, "", y_title.c_str());
          plotter.AddHistoTitle(title.c_str(), 0.039);
          can.Update();
          can.Print((tgt+prefix+"_CrossSection.pdf").c_str());
          can.Print((tgt+prefix+"_CrossSection.png").c_str());
        }

        { //Plotting cross section ratios
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
          //plotter.axis_maximum = 0.4;
          std::string tmp; 
          if (prefix=="nuke_pTmu")
          {
            tmp = "Muon p_{T}";
            //tempxsec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix=="nuke_pZmu")
          {
            tmp = "Muon p_{Z}";
            //tempxsec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix=="nuke_BjorkenX") tmp = "{X}";
          else if (prefix=="nuke_Erecoil") tmp = "E_{recoil}";
          else if (prefix=="nuke_Emu")
          {
            tmp = "E_{#mu}";
            //tempxsec->GetXaxis()->SetRange(1, 11);
          } 

          std::string title = tmp + " " + tgt;
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;
          plotter.DrawDataMCRatio(crossSection, simulatedCrossSection, 1, "TR", "Data");
          plotter.AddHistoTitle(title.c_str(), 0.039);
          can.Update();
          can.Print((tgt+prefix+"_CrossSectionRatio.pdf").c_str());
          can.Print((tgt+prefix+"_CrossSectionRatio.png").c_str());
        }

        { //Plotting sidebands
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kDefaultStyle);
          //plotter.axis_maximum = 0.4;

          TObjArray* arr = new TObjArray();
          auto tempUSSignal = USSidebandSignal->Clone();
          auto tempUSUS = USSidebandUS->Clone();
          auto tempUSDS = USSidebandDS->Clone();
          auto tempUSOther = USSidebandOther->Clone();
          arr->Add(tempUSSignal);
          arr->Add(tempUSUS);
          arr->Add(tempUSDS);
          arr->Add(tempUSOther);
          auto tmpUSData = DataUSSideband;
          tmpUSData->GetXaxis()->SetRange(1, 81);
          std::string title =  tgt + " Upstream Plastic Sideband";
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_size_x=0.03;
          plotter.axis_title_offset_x=1.8;
          plotter.axis_title_offset_y=1.8;
          plotter.DrawDataStackedMC(tmpUSData, arr, (dataPOT/mcPOT), "TR", "Data", 2, 1, 3001, "MINER#nuA Segment", "# of Events");
          plotter.AddHistoTitle(title.c_str(), 0.039);
          can.Update();
          can.Print((tgt+prefix+"_USSideband.pdf").c_str());
          can.Print((tgt+prefix+"_USSideband.png").c_str());


          can.Clear();
          //std::string title =  tgt;
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_size_x=0.03;
          plotter.axis_title_offset_x=1.8;
          plotter.axis_title_offset_y=1.8;
          plotter.DrawDataStackedMC(tmpUSData, arr, USScaleFactor*(dataPOT/mcPOT), "TR", "Data", 2, 1, 3001, "MINER#nuA Segment", "# of Events");
          title =  tgt + " Upstream Plastic Sideband Tuned";
          plotter.AddHistoTitle(title.c_str(), 0.039);
          can.Update();
          can.Print((tgt+prefix+"_USSidebandScaled.pdf").c_str());
          can.Print((tgt+prefix+"_USSidebandScaled.png").c_str());
        }
        { //Plotting sidebands
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kDefaultStyle);
          //plotter.axis_maximum = 0.4;

          TObjArray* arr = new TObjArray();
          auto tempUSSignal = DSSidebandSignal->Clone();
          auto tempUSUS = DSSidebandUS->Clone();
          auto tempUSDS = DSSidebandDS->Clone();
          auto tempUSOther = DSSidebandOther->Clone();
          arr->Add(tempUSSignal);
          arr->Add(tempUSUS);
          arr->Add(tempUSDS);
          arr->Add(tempUSOther);
          auto tmpUSData = DataDSSideband;
          tmpUSData->GetXaxis()->SetRange(1, 81);
          std::string title =  tgt + " Downstream Plastic Sideband";
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_size_x=0.03;
          plotter.axis_title_offset_x=1.8;
          plotter.axis_title_offset_y=1.8;
          plotter.DrawDataStackedMC(tmpUSData, arr, (dataPOT/mcPOT), "TR", "Data", 2, 1, 3001, "MINER#nuA Segment", "# of Events");
          plotter.AddHistoTitle(title.c_str(), 0.039);
          can.Update();
          can.Print((tgt+prefix+"_DSSideband.pdf").c_str());
          can.Print((tgt+prefix+"_DSSideband.png").c_str());


          can.Clear();
          //std::string title =  tgt;
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_size_x=0.03;
          plotter.axis_title_offset_x=1.8;
          plotter.axis_title_offset_y=1.8;
          plotter.DrawDataStackedMC(tmpUSData, arr, DSScaleFactor*(dataPOT/mcPOT), "TR", "Data", 2, 1, 3001, "MINER#nuA Segment", "# of Events");
          title =  tgt + " Downstream Plastic Sideband Tuned";
          plotter.AddHistoTitle(title.c_str(), 0.039);
          can.Update();
          can.Print((tgt+prefix+"_DSSidebandScaled.pdf").c_str());
          can.Print((tgt+prefix+"_DSSidebandScaled.png").c_str());
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
          plotter.DrawNormalizedMigrationHistogram(migration, false, false, true, true);
          //plotter.AddHistoTitle(title.c_str(), 0.039);
          can.Update();
          //can.Print((tgt+prefix+"_DSSidebandScaled.pdf").c_str());
          can.Print((tgt+prefix+"_migrationTest.png").c_str());
        }
        /*




          auto trackerFile = TFile::Open("./tracker_BjorkenX_Daisy_crossSection.root", "READ");
          if (!trackerFile)
          {
            std::cerr << "Failed to open tracker file\n";
            return 2;
          }
          auto trackerFileNoDaisy = TFile::Open("./tracker_Emu_crossSection.root", "READ");
          if (!trackerFileNoDaisy)
          {
            std::cerr << "Failed to open tracker no daisy file\n";
            return 2;
          }
          std::string titlebase;
          PlotUtils::MnvH1D* plasticxsec;
          if (tgt == "Iron")
          {
            titlebase = "Fe";
            plasticxsec=util::GetIngredient<PlotUtils::MnvH1D>(*trackerFile, "tracker_Emu_Fe_CrossSection");
          }
          if (tgt == "Lead")
          {
            titlebase="Pb";
            plasticxsec=util::GetIngredient<PlotUtils::MnvH1D>(*trackerFile, "tracker_Emu_Pb_CrossSection");
          }
          if (tgt == "Carbon")
          {
            titlebase="C";
            plasticxsec=util::GetIngredient<PlotUtils::MnvH1D>(*trackerFile, "tracker_Emu_C_CrossSection");
          }
          if (tgt == "6000")
          {
            titlebase="H_{2}O";
            plasticxsec=util::GetIngredient<PlotUtils::MnvH1D>(*trackerFileNoDaisy, "crossSection");
          }
          titlebase+="/CH Ratio in E_#mu";
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
          //plotter.axis_maximum = 0.4;
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;
          plotter.DrawDataMCRatio(crossSectionTuned, plasticxsec, 1, "TR", "Data");
          plotter.AddHistoTitle(titlebase.c_str(), 0.039);
          can.Update();
          can.Print((prefix+"_"+tgt+"CHRatio.pdf").c_str());
          can.Print((prefix+"_"+tgt+"CHRatio.png").c_str());
 */

        outFile->Close();
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
