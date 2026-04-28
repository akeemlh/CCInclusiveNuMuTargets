#define HELP \
"\n*** Help: ***\n"\
" File: PreExtractionTargets.cpp\n"\
" Brief: Given data and MC files produced by runEventLoopTargets.cpp, performs the neccessary pre-processing required to perform warping stuides\n"\
"        And later extract a 1D differential cross section. It basically per-performs and extracts everything needed for a cross section, and stops\n"\
"        Just short of unfolding\n\n"\
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
"        Performs plastic sideband scaling, subtracts backgrounds \n"\
"        acceptance correction, divides by flux and number of nucleons and if option is selected .\n"\
"        Writes a .root file with the cross section histograms\n"\
"        To run for a single playlist (for example 1A) simply pass /path/to/dirs/1A as the directory path\n\n"\
" Usage: PreExtractionTargets <directory> <target> <pdg>\n"\
"        e.g:   PreExtractionTargets /path/to/dirs 2026 14 -- to extract xsecs for target 2 Iron over neutrino-mode playlists with 5 iterations\n"\
"        e.g:   PreExtractionTargets /path/to/dirs 6000 -14 -- to extract xsecs for the water target over antineutrino-mode playlists with 10 iterations\n"\
"        e.g:   PreExtractionTargets /path/to/dirs iron -14 -- to extract xsecs for the combined iron targets over the antineutrino-mode playlists with 10 iterations\n"\
"        e.g:   PreExtractionTargets /path/to/dirs all 14 -- to extract xsecs for all targets over the neutrino-mode playlists with 5 iterations\n"\
"        e.g:   PreExtractionTargets /path/to/dirs WaterFull -14 -- to extract xsecs for the water target over the antineutrino-mode playlists with the water target filled with 10 iterations\n"\
"        e.g:   PreExtractionTargets /path/to/dirs WaterEmpty -14 -- to extract xsecs for the water target over the antineutrino-mode playlists with the water target empty with 10 iterations\n\n"\

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


template <class Hist>
struct Chi2Calculator {
    Hist* dataHist;
    Hist* mcFixedHist;
    Hist* mcFloatHist;

    Chi2Calculator(Hist* data, Hist* mc_fix, Hist* mc_float) : dataHist(data), mcFixedHist(mc_fix), mcFloatHist(mc_float) {}

    /* double operator()(const double* val) const {
      double scale  = val[0];
      int ndof;
      PlotUtils::MnvPlotter plotter;
      auto tmpHist = dataHist->Clone("tmp_minimization_hist");
      double chiSq = plotter.Chi2DataMC(tmpHist, mcHist, ndof, scale, true, true); //Is this implemented properly for TH2Ds
      delete tmpHist; 
      return chiSq;
    } */
    double operator()(const double* val) const {
        double scale = val[0];
        double chi2 = 0.0;
        
        int nBinsX = mcFloatHist->GetNbinsX();
        // ROOT Magic: For TH1D, GetNbinsY() and GetNbinsZ() automatically return 1
        int nBinsY = mcFloatHist->GetNbinsY(); 
        int nBinsZ = mcFloatHist->GetNbinsZ();
        
        // Loop over the maximum dimensionality
        for (int x = 1; x <= nBinsX; ++x) {
            for (int y = 1; y <= nBinsY; ++y) {
                for (int z = 1; z <= nBinsZ; ++z) {
                    
                    // Convert multi-dimensional X/Y/Z into ROOT's 1D global bin index
                    int globalBin = mcFloatHist->GetBin(x, y, z);
                    
                    double dataVal = dataHist->GetBinContent(globalBin);
                    double dataErr = dataHist->GetBinError(globalBin);

                    double mcFixedVal = mcFixedHist->GetBinContent(globalBin);
                    double mcFixedErr = mcFixedHist->GetBinError(globalBin);
                    
                    double mcFloatVal = mcFloatHist->GetBinContent(globalBin) * scale;
                    double mcFloatErr = mcFloatHist->GetBinError(globalBin) * scale;
                    
                    // Only calculate chi2 for bins that physically contain events.
                    // This natively masks empty kinematic regions in 2D/3D spaces.
                    double totalVariance = (dataErr * dataErr) + (mcFixedErr * mcFixedErr) + (mcFixedErr * mcFixedErr);
                    
                    if (totalVariance > 0) {
                        double diff = dataVal - (mcFixedVal + mcFloatVal);
                        chi2 += (diff * diff) / totalVariance;
                    }
                }
            }
        }
        
        return chi2;
    }

};

bool isNumber(const std::string& str) {
  for (char c : str) {
    if (!std::isdigit(c)) {
      return false;
    }
  }
  return true;
}

template <typename THist>
std::vector<double> RunMinuitFit(THist* data, THist* mc_fix, THist* mc_float) {
    Chi2Calculator<THist> chi2Calc(data, mc_fix, mc_float);
    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2");
    minimizer->SetMaxFunctionCalls(1000000);
    minimizer->SetMaxIterations(100000);
    minimizer->SetTolerance(0.01);
    
    // CRITICAL: Set to 0 to avoid massive log spam when looping hundreds of universes
    minimizer->SetPrintLevel(0); 
    
    ROOT::Math::Functor f(chi2Calc, 1);
    minimizer->SetFunction(f);
    minimizer->SetVariable(0, "scale", 1.0, 0.001);
    
    minimizer->Minimize();
    
    double scale = minimizer->X()[0];
    double statErr = minimizer->Errors()[0];
    double chi2 = minimizer->MinValue();
    delete minimizer;
    
    return {scale, statErr, chi2};
}

// Plot a step in cross section extraction.
void Plot(PlotUtils::MnvH1D &hist, const std::string &stepName, const std::string &prefix, const std::string &target)
{
  //return;
  bool plotpngs = true;
  TCanvas can(stepName.c_str());
  hist.GetCVHistoWithError().Clone()->Draw();
  if (plotpngs) can.Print((target + "_" + prefix + "_" + stepName + ".png").c_str());

  // Uncertainty summary
  PlotUtils::MnvPlotter plotter;
  plotter.ApplyStyle(PlotUtils::kCCQENuInclusiveStyle);
  plotter.axis_maximum = 0.4;

  plotter.DrawErrorSummary(&hist);
  if (plotpngs) can.Print((target + "_" + prefix + "_" + stepName + "_uncertaintySummary.png").c_str());

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "Other");
  if (plotpngs) can.Print((target + "_" + prefix + "_" + stepName + "_otherUncertainties.png").c_str());
}

void Plot(PlotUtils::MnvH2D &hist, const std::string &stepName, const std::string &prefix, const std::string &target)
{
  return;
}

// The final step of cross section extraction: normalize by flux, bin width, POT, and number of targets
template <class MnvHist> 
MnvHist *normalize(MnvHist *efficiencyCorrected, MnvHist *fluxIntegral, const double nNucleons, const double POT)
{
  efficiencyCorrected->Divide(efficiencyCorrected, fluxIntegral);

  efficiencyCorrected->Scale(1. / nNucleons / POT);
  efficiencyCorrected->Scale(1.e4); // Flux histogram is in m^-2, but convention is to report cm^2
  efficiencyCorrected->Scale(1., "width");

  return efficiencyCorrected;
}

template <class MnvHist>
void ProcessTargets(const std::string& prefix, 
                         const std::string& tgtname, 
                         int pdg, 
                         const std::vector<std::string>& dirs,
                         int sidebandSubtraction)
{
    MnvHist *fluxIntReweighted = new MnvHist();

    MnvHist* flux = new MnvHist();
    MnvHist* folded = new MnvHist();
    PlotUtils::MnvH2D* migration = new PlotUtils::MnvH2D();
    MnvHist* effNum = new MnvHist();
    MnvHist* effDenom = new MnvHist();
    MnvHist* effDenom2P2H = new MnvHist();
    MnvHist* effDenomDIS = new MnvHist();
    MnvHist* effDenomRES = new MnvHist();
    MnvHist* effDenomQE = new MnvHist();
    MnvHist* effDenomOther = new MnvHist();

    MnvHist* USSidebandSignal = new MnvHist();
    MnvHist* USSidebandDS = new MnvHist();
    MnvHist* USSidebandUS = new MnvHist();
    MnvHist* USSidebandOther = new MnvHist();
    MnvHist* DSSidebandSignal = new MnvHist();
    MnvHist* DSSidebandDS = new MnvHist();
    MnvHist* DSSidebandUS = new MnvHist();
    MnvHist* DSSidebandOther = new MnvHist();
    MnvHist* DataUSSideband = new MnvHist();
    MnvHist* DataDSSideband = new MnvHist();
    MnvHist* DataSignal = new MnvHist();

    MnvHist* SelectedSignalReco = new MnvHist();

    MnvHist* BackgroundWrongSign = new MnvHist();
    MnvHist* BackgroundNC = new MnvHist();
    MnvHist* BackgroundOther = new MnvHist();
    MnvHist* BackgroundWaterTank = new MnvHist();
    MnvHist* BackgroundUSPlastic = new MnvHist();
    MnvHist* BackgroundDSPlastic = new MnvHist();
    MnvHist* BackgroundOtherTarget = new MnvHist();
    MnvHist* BackgroundVtxElsewhere = new MnvHist();

    MnvHist* MCData = new MnvHist();
    MnvHist* MCData2p2h = new MnvHist();
    MnvHist* MCDataDIS = new MnvHist();
    MnvHist* MCDataRES = new MnvHist();
    MnvHist* MCDataQE = new MnvHist();
    MnvHist* MCDataOther = new MnvHist();

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
        if (tgtname=="Iron") targetsInTgt = {"2026", "3026", "5026"};
        else if (tgtname=="Carbon") targetsInTgt = {"3006"};
        else if (tgtname=="Lead") targetsInTgt = {"2082", "3082", "4082", "5082"};
        else if (tgtname=="WaterEmpty")
        {
          targetsInTgt = {"6000"};
          waterFilledEmpty = 2;
        }
        else if (tgtname=="WaterFull")
        {
          targetsInTgt = {"6000"};
          waterFilledEmpty = 1;
        }
        else if (tgtname=="Water")
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
              //return;
            }
            double data_pot = util::GetIngredientSmart<TParameter<double>>(*dataFile, "POTUsed")->GetVal();
            auto mcFile = TFile::Open(mcpath.c_str(), "READ");
            if (!mcFile)
            {
              std::cerr << "Failed to open MC file " << mcpath.c_str() << ".\n";
              continue;
              //return;
            }
            double mc_pot = util::GetIngredientSmart<TParameter<double>>(*mcFile, "POTUsed")->GetVal();
            std::string playlistUsed = util::GetIngredientSmart<TNamed>(*mcFile, "PlaylistUsed")->GetTitle();
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
        else targetsInTgt = {tgtname};
        std::cout<<"Investigating directory1\n";
        //Flux parameters
        int n_flux_universes = 100; // Is this right
        const bool use_nue_constraint = true;
        const std::string project_dir = "targets_2345_jointNueIMD";
        double min_energy = 0;
        double max_energy = 120;

        double mcPOT = 0;
        double dataPOT = 0;
        std::cout<<"Investigating directory2\n";

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

            auto dataFile = TFile::Open(datapath.c_str(), "READ");
            if (!dataFile)
            {
              std::cerr << "Failed to open data file " << datapath.c_str() << ".\n";
              continue;
              //return;
            }

            auto mcFile = TFile::Open(mcpath.c_str(), "READ");
            if (!mcFile)
            {
              std::cerr << "Failed to open MC file " << mcpath.c_str() << ".\n";
              continue;
              //return;
            }
            
            std::string playlistUsed = util::GetIngredientSmart<TNamed>(*mcFile, "PlaylistUsed")->GetTitle();
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


            //util::AddHist(*flux,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("reweightedflux_integrated")), prefix).get(), mcscale);
            util::AddHist(*folded,util::GetIngredientSmart<MnvHist>(*dataFile, (std::string("data")), prefix).get(), datascale);
            util::AddHist(*migration,util::GetIngredientSmart<PlotUtils::MnvH2D>(*mcFile, (std::string("migration")), prefix).get(), mcscale);
            util::AddHist(*effNum,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("efficiency_numerator")), prefix).get(), mcscale);
            util::AddHist(*effDenom,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("efficiency_denominator")), prefix).get(), mcscale);
            util::AddHist(*effDenom2P2H,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("efficiency_denominator_intChannels_2p2h")), prefix).get(), mcscale);
            util::AddHist(*effDenomDIS,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("efficiency_denominator_intChannels_DIS")), prefix).get(), mcscale);
            util::AddHist(*effDenomRES,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("efficiency_denominator_intChannels_RES")), prefix).get(), mcscale);
            util::AddHist(*effDenomQE,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("efficiency_denominator_intChannels_QE")), prefix).get(), mcscale);
            util::AddHist(*effDenomOther,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("efficiency_denominator_intChannels_Other")), prefix).get(), mcscale);
            
            //If doing scalar sideband tuning
            if (sidebandSubtraction == 1)
            {
              util::AddHist(*USSidebandSignal,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("segment_US_sideband_Signal"))).get(), mcscale);
              util::AddHist(*USSidebandDS,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("segment_US_sideband_DS"))).get(), mcscale);
              util::AddHist(*USSidebandUS,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("segment_US_sideband_US"))).get(), mcscale);
              util::AddHist(*USSidebandOther,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("segment_US_sideband_Other"))).get(), mcscale);
              util::AddHist(*DSSidebandSignal,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("segment_DS_sideband_Signal"))).get(), mcscale);
              util::AddHist(*DSSidebandDS,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("segment_DS_sideband_DS"))).get(), mcscale);
              util::AddHist(*DSSidebandUS,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("segment_DS_sideband_US"))).get(), mcscale);
              util::AddHist(*DSSidebandOther,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("segment_DS_sideband_Other"))).get(), mcscale);
              util::AddHist(*DataUSSideband,util::GetIngredientSmart<MnvHist>(*dataFile, (std::string("segment_US_Sideband"))).get(), datascale);
              util::AddHist(*DataDSSideband,util::GetIngredientSmart<MnvHist>(*dataFile, (std::string("segment_DS_Sideband"))).get(), datascale);
              util::AddHist(*DataSignal,util::GetIngredientSmart<MnvHist>(*dataFile, (std::string("segment_data"))).get(), datascale);
            }
            //If doing bin-by-bin sideband tuning
            // Typically don't use since it's poorly defined for BjorkenX (since BjorkenX is so smeared)
            // and since we have empty bins in the 2D analysis, so a bin by bin scaling is impossible for it
            if (sidebandSubtraction == 2)
            {
              util::AddHist(*USSidebandSignal,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("US_sideband_Signal")), prefix).get(), mcscale);
              util::AddHist(*USSidebandDS,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("US_sideband_DS")), prefix).get(), mcscale);
              util::AddHist(*USSidebandUS,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("US_sideband_US")), prefix).get(), mcscale);
              util::AddHist(*USSidebandOther,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("US_sideband_Other")), prefix).get(), mcscale);
              util::AddHist(*DSSidebandSignal,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("DS_sideband_Signal")), prefix).get(), mcscale);
              util::AddHist(*DSSidebandDS,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("DS_sideband_DS")), prefix).get(), mcscale);
              util::AddHist(*DSSidebandUS,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("DS_sideband_US")), prefix).get(), mcscale);
              util::AddHist(*DSSidebandOther,util::GetIngredientSmart<MnvHist>(*mcFile, (std::string("DS_sideband_Other")), prefix).get(), mcscale);
              util::AddHist(*DataUSSideband,util::GetIngredientSmart<MnvHist>(*dataFile, (std::string("US_Sideband")), prefix).get(), datascale);
              util::AddHist(*DataDSSideband,util::GetIngredientSmart<MnvHist>(*dataFile, (std::string("DS_Sideband")), prefix).get(), datascale);
              util::AddHist(*DataSignal,util::GetIngredientSmart<MnvHist>(*dataFile, (std::string("data")), prefix).get(), datascale);
            }


            util::AddHist(*SelectedSignalReco,util::GetIngredient<MnvHist>(*mcFile, (std::string("selected_signal_reco")), prefix), mcscale);
            util::AddHist(*BackgroundWrongSign,util::GetIngredient<MnvHist>(*mcFile, (std::string("background_Wrong_Sign_Bkg")), prefix), mcscale);
            util::AddHist(*BackgroundNC,util::GetIngredient<MnvHist>(*mcFile, (std::string("background_NC_Bkg")), prefix), mcscale);
            util::AddHist(*BackgroundWaterTank,util::GetIngredient<MnvHist>(*mcFile, (std::string("background_Water_Tank_Bkg")), prefix), mcscale);
            util::AddHist(*BackgroundDSPlastic,util::GetIngredient<MnvHist>(*mcFile, (std::string("background_Downstream_Plastic_Bkg")), prefix), mcscale);
            util::AddHist(*BackgroundUSPlastic,util::GetIngredient<MnvHist>(*mcFile, (std::string("background_Upstream_Plastic_Bkg")), prefix), mcscale);
            util::AddHist(*BackgroundOtherTarget,util::GetIngredient<MnvHist>(*mcFile, (std::string("background_True_In_Other_Target_Bkg")), prefix), mcscale);
            util::AddHist(*BackgroundVtxElsewhere,util::GetIngredient<MnvHist>(*mcFile, (std::string("background_True_Vtx_Elsewhere_Bkg")), prefix), mcscale);
            util::AddHist(*BackgroundOther,util::GetIngredient<MnvHist>(*mcFile, (std::string("background_Other")), prefix), mcscale);

            util::AddHist(*MCData,util::GetIngredient<MnvHist>(*mcFile, (std::string("data")), prefix), mcscale);
            util::AddHist(*MCData2p2h,util::GetIngredient<MnvHist>(*mcFile, (std::string("intType_2p2h")), prefix), mcscale);
            util::AddHist(*MCDataDIS,util::GetIngredient<MnvHist>(*mcFile, (std::string("intType_DIS")), prefix), mcscale);
            util::AddHist(*MCDataRES,util::GetIngredient<MnvHist>(*mcFile, (std::string("intType_RES")), prefix), mcscale);
            util::AddHist(*MCDataQE,util::GetIngredient<MnvHist>(*mcFile, (std::string("intType_QE")), prefix), mcscale);
            util::AddHist(*MCDataOther,util::GetIngredient<MnvHist>(*mcFile, (std::string("intType_Other")), prefix), mcscale);

            if (d==0) //Only get POT once per playlist. E.g we don't wanna double count the POT for the same runs just because we're looking at target 2026 and 3026
            {
              mcPOT += util::GetIngredientSmart<TParameter<double>>(*mcFile, "POTUsed")->GetVal();
              std::cout<<"mcPOT: " << mcPOT << std::endl;
              double tempDataPOT = util::GetIngredientSmart<TParameter<double>>(*dataFile, "POTUsed")->GetVal();
              dataPOT+=tempDataPOT;
              std::cout<<"dataPOT: " << dataPOT << std::endl;

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
        Plot(*folded, "folded", prefix, tgtname);

        std::string nutypestr;
        if (pdg>0) nutypestr = "Neutrino";
        else nutypestr = "Antineutrino";

        std::string outpath = nutypestr+"_"+tgtname +"_" + prefix + "_preProcessed.root";
        auto outFile = TFile::Open(outpath.c_str(), "RECREATE");
        if (!outFile)
        {
            std::cerr << "Could not create a file called " << outpath << ".\n";
            //return;
        }

        MnvHist* mc_reco_untuned = (MnvHist*)SelectedSignalReco->Clone("mc_reco_untuned");
        mc_reco_untuned->Add(BackgroundWrongSign);
        mc_reco_untuned->Add(BackgroundNC);
        mc_reco_untuned->Add(BackgroundWaterTank);
        mc_reco_untuned->Add(BackgroundOtherTarget);
        mc_reco_untuned->Add(BackgroundVtxElsewhere);
        mc_reco_untuned->Add(BackgroundOther);
        mc_reco_untuned->Add(BackgroundUSPlastic);
        mc_reco_untuned->Add(BackgroundDSPlastic);
        mc_reco_untuned->Write("mc_reco_untuned");
        BackgroundUSPlastic->Write("BackgroundUSPlasticUntuned");
        BackgroundDSPlastic->Write("BackgroundDSPlasticUntuned");


        if (isNumber(tgtname))
        {
          int tgtCode = std::stoi(tgtname);
          if (tgtCode>=7 && tgtCode <= 19) sidebandSubtraction = 0; // No plastic sideband tuning possible for the pseudotargets
        } 
        double USScaleFactor = 1;
        double DSScaleFactor = 1;
        //************************************************
        //Sideband chisq minimisation - single scale factor across a given sideband region
        //************************************************
        if (sidebandSubtraction==1)
        {
          MnvHist* combinedUSFixed = USSidebandSignal->Clone("combinedUSFixedMC");
          //combinedUSFixed->Add(USSidebandUS);
          combinedUSFixed->Add(USSidebandDS);
          combinedUSFixed->Add(USSidebandOther);
          //GetCVHistoWithStatError();
          combinedUSFixed->Scale(dataPOT/mcPOT);

          MnvHist* combinedDSFixed = DSSidebandSignal->Clone("combinedDSFixedMC");
          combinedDSFixed->Add(DSSidebandUS);
          //combinedDSFixed->Add(DSSidebandDS);
          combinedDSFixed->Add(DSSidebandOther);
          combinedDSFixed->Scale(dataPOT/mcPOT);
          
          MnvHist* DSSidebandDSScaled = (MnvHist*) DSSidebandDS->Clone("DSSidebandDS_Scaled");
          MnvHist* USSidebandUSScaled = (MnvHist*) USSidebandUS->Clone("USSidebandUS_Scaled");
          DSSidebandDSScaled->Scale(dataPOT/mcPOT);
          USSidebandUSScaled->Scale(dataPOT/mcPOT);


          // Helper to fill all bins of a histogram with the scale factor 
          auto fillScaleHisto = [](TH1* hist, double val, double err = 0.0) {
              for(int i = 0; i < hist->GetNcells(); ++i) {
                  hist->SetBinContent(i, val);
                  hist->SetBinError(i, err);
              }
          };

          MnvHist* scaleHistoUS = (MnvHist*)DataUSSideband->Clone("scaleHisto_US");
          scaleHistoUS->ClearAllErrorBands();
          scaleHistoUS->Reset();
          scaleHistoUS->AddMissingErrorBandsAndFillWithCV(*combinedUSFixed); // Ensure it has MC error bands

          MnvHist* scaleHistoDS = (MnvHist*)DataDSSideband->Clone("scaleHisto_DS");
          scaleHistoDS->ClearAllErrorBands();
          scaleHistoDS->Reset();
          scaleHistoDS->AddMissingErrorBandsAndFillWithCV(*combinedDSFixed);

          // -------------------------------------------------------------
          // 1. Fit the Central Value (CV) First
          // -------------------------------------------------------------
          TH1* cvDataUS = (TH1*)DataUSSideband->GetCVHistoWithStatError().Clone("cvDataUS");
          TH1* cvMCUSFixed   = (TH1*)combinedUSFixed->GetCVHistoWithStatError().Clone("cvMCUSFixed");
          TH1* cvMCUSFloat   = (TH1*)USSidebandUSScaled->GetCVHistoWithStatError().Clone("cvMCUSFloat");
          auto resultsUS = RunMinuitFit(cvDataUS, cvMCUSFixed, cvMCUSFloat);
          USScaleFactor = resultsUS[0];
          fillScaleHisto(scaleHistoUS, USScaleFactor, resultsUS[1]);
          //fillScaleHisto(scaleHistoUS->GetCVHistoWithStatError().GetHist(), USScaleFactor, resultsUS[1]);
          
          TH1* cvDataDS = (TH1*)DataDSSideband->GetCVHistoWithStatError().Clone("cvDataDS");
          TH1* cvMCDSFixed   = (TH1*)combinedDSFixed->GetCVHistoWithStatError().Clone("cvMCDS");
          TH1* cvMCDSFloat   = (TH1*)DSSidebandDSScaled->GetCVHistoWithStatError().Clone("cvMCDS");
          auto resultsDS = RunMinuitFit(cvDataDS, cvMCDSFixed, cvMCDSFloat);
          DSScaleFactor = resultsDS[0];
          fillScaleHisto(scaleHistoDS, DSScaleFactor, resultsDS[1]);
          //fillScaleHisto(scaleHistoDS->GetCVHistoWithStatError().GetHist(), DSScaleFactor, resultsDS[1]);

          std::cout << "CV US Scale: " << USScaleFactor << " +/- " << resultsUS[1] << " (Stat Err)" << std::endl;
          std::cout << "CV DS Scale: " << DSScaleFactor << " +/- " << resultsDS[1] << " (Stat Err)" << std::endl;

          // Scale the CV Backgrounds
          BackgroundUSPlastic->Scale(USScaleFactor);
          BackgroundDSPlastic->Scale(DSScaleFactor);
          
          delete cvDataUS; delete cvMCUSFixed; delete cvMCUSFloat; delete cvDataDS; delete cvMCDSFixed; delete cvMCDSFloat;

          // -------------------------------------------------------------
          // 2. Loop Over Vertical Systematics Universes (e.g., Flux, GENIE)
          // -------------------------------------------------------------
          std::vector<std::string> vertBandNames = combinedUSFixed->GetVertErrorBandNames();
          for (const std::string& bandName : vertBandNames) {
              auto bandUSFixed_MC = combinedUSFixed->GetVertErrorBand(bandName);
              auto bandDSFixed_MC = combinedDSFixed->GetVertErrorBand(bandName);
              auto bandUSFloat_MC = USSidebandUSScaled->GetVertErrorBand(bandName);
              auto bandDSFloat_MC = DSSidebandDSScaled->GetVertErrorBand(bandName);
              auto bandUS_Bkg = BackgroundUSPlastic->GetVertErrorBand(bandName);
              auto bandDS_Bkg = BackgroundDSPlastic->GetVertErrorBand(bandName);

              auto bandScale_US = scaleHistoUS->GetVertErrorBand(bandName);
              auto bandScale_DS = scaleHistoDS->GetVertErrorBand(bandName);
              for (unsigned int i = 0; i < bandUSFixed_MC->GetNHists(); ++i) {
                  TH1* univFixedMC_US = bandUSFixed_MC->GetHist(i);
                  TH1* univFixedMC_DS = bandDSFixed_MC->GetHist(i);
                  
                  TH1* univFloatMC_US = bandUSFloat_MC->GetHist(i);
                  TH1* univFloatMC_DS = bandDSFloat_MC->GetHist(i);
                  
                  // Data does not have vertical errors, so we fit against CV Data
                  TH1* univData_US = (TH1*)DataUSSideband->GetCVHistoWithStatError().Clone("univDataUS");
                  TH1* univData_DS = (TH1*)DataDSSideband->GetCVHistoWithStatError().Clone("univDataDS");

                  double univUSScale = RunMinuitFit(univData_US, univFixedMC_US, univFloatMC_US)[0];
                  double univDSScale = RunMinuitFit(univData_DS, univFixedMC_DS, univFloatMC_DS)[0];

                  bandUS_Bkg->GetHist(i)->Scale(univUSScale);
                  bandDS_Bkg->GetHist(i)->Scale(univDSScale);
                  fillScaleHisto(bandScale_US->GetHist(i), univUSScale);
                  fillScaleHisto(bandScale_DS->GetHist(i), univDSScale);
                  delete univData_US; delete univData_DS;
              }
          }

          // -------------------------------------------------------------
          // 3. Loop Over Lateral Systematics Universes (e.g., Muon Tracking)
          // -------------------------------------------------------------
          std::vector<std::string> latBandNames = combinedUSFixed->GetLatErrorBandNames();
          for (const std::string& bandName : latBandNames) {
              auto bandUSFixed_MC = combinedUSFixed->GetLatErrorBand(bandName);
              auto bandDSFixed_MC = combinedDSFixed->GetLatErrorBand(bandName);
              auto bandUSFloat_MC = USSidebandUSScaled->GetLatErrorBand(bandName);
              auto bandDSFloat_MC = DSSidebandDSScaled->GetLatErrorBand(bandName);
              auto bandUS_Bkg = BackgroundUSPlastic->GetLatErrorBand(bandName);
              auto bandDS_Bkg = BackgroundDSPlastic->GetLatErrorBand(bandName);

              auto bandScale_US = scaleHistoUS->GetLatErrorBand(bandName);
              auto bandScale_DS = scaleHistoDS->GetLatErrorBand(bandName);

              // Lateral bands sometimes exist in data (e.g. momentum shifts applied to Data instead of MC)
              bool dataHasBandUS = DataUSSideband->HasLatErrorBand(bandName);
              bool dataHasBandDS = DataDSSideband->HasLatErrorBand(bandName);

              for (unsigned int i = 0; i < bandUSFixed_MC->GetNHists(); ++i) {
                  TH1* univFixedMC_US = bandUSFixed_MC->GetHist(i);
                  TH1* univFixedMC_DS = bandDSFixed_MC->GetHist(i);
                  
                  TH1* univFloatMC_US = bandUSFloat_MC->GetHist(i);
                  TH1* univFloatMC_DS = bandDSFloat_MC->GetHist(i);

                  TH1* univData_US = dataHasBandUS ? DataUSSideband->GetLatErrorBand(bandName)->GetHist(i) : (TH1*)DataUSSideband->GetCVHistoWithStatError().Clone("univDataUS");
                  TH1* univData_DS = dataHasBandDS ? DataDSSideband->GetLatErrorBand(bandName)->GetHist(i) : (TH1*)DataDSSideband->GetCVHistoWithStatError().Clone("univDataDS");

                  double univUSScale = RunMinuitFit(univData_US, univFixedMC_US, univFloatMC_US)[0];
                  double univDSScale = RunMinuitFit(univData_DS, univFixedMC_DS, univFloatMC_DS)[0];

                  bandUS_Bkg->GetHist(i)->Scale(univUSScale);
                  bandDS_Bkg->GetHist(i)->Scale(univDSScale);

                  fillScaleHisto(bandScale_US->GetHist(i), univUSScale);
                  fillScaleHisto(bandScale_DS->GetHist(i), univDSScale);

                  if (!dataHasBandUS) delete univData_US;
                  if (!dataHasBandDS) delete univData_DS;
              }
          }
          scaleHistoUS->Write("scaleFactor_US");
          scaleHistoDS->Write("scaleFactor_DS");
          /* 
          //Plot(*combinedUS, "combinedUS", prefix, tgtname);
          //Plot(*DataUSSideband, "DataUSSideband", prefix, tgtname);
          double minValUS, minValDS;
          {          //US Sideband
            Chi2Calculator<MnvHist> chi2Calc(DataUSSideband, combinedUS);
            ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2");
            minimizer->SetMaxFunctionCalls(1000000);
            minimizer->SetMaxIterations(100000);
            minimizer->SetTolerance(0.01);
            minimizer->SetPrintLevel(1);
            ROOT::Math::Functor f(chi2Calc, 1);
            minimizer->SetFunction(f);
            minimizer->SetVariable(0, "scale", 1, 0.001);
            minimizer->Minimize();
            const double *bestfit = minimizer->X();
            double bestChiSq = minimizer->MinValue();
            USScaleFactor = bestfit[0];
            std::cout<<"Minised chisq for USSideband. Scale factor: " << USScaleFactor << " best chisq: " << bestChiSq << std::endl;
          }

          {          //DS Sideband
            Chi2Calculator<MnvHist> chi2Calc(DataUSSideband, combinedUS);
            ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2");
            minimizer->SetMaxFunctionCalls(1000000);
            minimizer->SetMaxIterations(100000);
            minimizer->SetTolerance(0.01);
            minimizer->SetPrintLevel(1);
            ROOT::Math::Functor f(chi2Calc, 1);
            minimizer->SetFunction(f);
            minimizer->SetVariable(0, "scale", 1, 0.001);
            minimizer->Minimize();
            const double *bestfit = minimizer->X();
            double bestChiSq = minimizer->MinValue();
            DSScaleFactor = bestfit[0];
            std::cout<<"Minised chisq for DSSideband. Scale factor: " << DSScaleFactor << " best chisq: " << bestChiSq << std::endl;
          }
          BackgroundUSPlastic->Scale(USScaleFactor);
          BackgroundDSPlastic->Scale(DSScaleFactor); */
        }
        //************************************************
        //End - Sideband chisq minimisation - single scale factor across a given sideband region
        //************************************************
        else if (sidebandSubtraction==2)
        {
          MnvHist* combinedUS = USSidebandSignal->Clone("CombinedUSMC");
          combinedUS->Add(USSidebandUS);
          combinedUS->Add(USSidebandDS);
          combinedUS->Add(USSidebandOther);
          //GetCVHistoWithStatError();
          combinedUS->Scale(dataPOT/mcPOT);

          MnvHist* combinedDS = DSSidebandSignal->Clone("CombinedUSMC");
          combinedDS->Add(DSSidebandUS);
          combinedDS->Add(DSSidebandDS);
          combinedDS->Add(DSSidebandOther);
          combinedDS->Scale(dataPOT/mcPOT);

          MnvHist* USScale = DataUSSideband->Clone("UpstreamScaling");
          MnvHist* DSScale = DataDSSideband->Clone("DownstreamScaling");

          USScale->Divide(USScale, combinedUS);
          DSScale->Divide(DSScale, combinedDS);

          BackgroundUSPlastic->Multiply(BackgroundUSPlastic, USScale);
          BackgroundDSPlastic->Multiply(BackgroundDSPlastic, DSScale);
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

        // TODO: Remove these debugging plots when done


        Plot(*BackgroundUSPlastic, "BackgroundUSPlastic", prefix, tgtname);
        Plot(*BackgroundDSPlastic, "BackgroundDSPlastic", prefix, tgtname);

        std::vector<MnvHist*> sidebandScaled = {BackgroundUSPlastic, BackgroundDSPlastic};
        std::vector<MnvHist*> backgrounds = { BackgroundWrongSign, BackgroundNC, BackgroundOther, BackgroundWaterTank, BackgroundOtherTarget, BackgroundVtxElsewhere};
        auto toSubtract = std::accumulate(std::next(backgrounds.begin()), backgrounds.end(), (*backgrounds.begin())->Clone(),
                                          [](auto sum, const auto hist)
                                          {
                                            sum->Add(hist);
                                            return sum;
                                          });
        Plot(*toSubtract, "BackgroundSum", prefix, tgtname);
        auto toSubtractSideband = std::accumulate(std::next(sidebandScaled.begin()), sidebandScaled.end(), (*sidebandScaled.begin())->Clone(),
                                          [](auto sum, const auto hist)
                                          {
                                            sum->Add(hist);
                                            return sum;
                                          });
        Plot(*toSubtractSideband, "sidebandScaledSum", prefix, tgtname);

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

        //sum up backgrounds without sideband scaling
        auto MCsidebandScaledSubtracted = std::accumulate(backgrounds.begin(), backgrounds.end(), MCData->Clone(),
                                              [](auto sum, const auto hist)
                                              {
                                                sum->Add(hist, -1);
                                                return sum;
                                              });
        //Both background subtracted and sideband scaled
        auto MCbkgScaledSubtracted = std::accumulate(sidebandScaled.begin(), sidebandScaled.end(), MCsidebandScaledSubtracted->Clone(),
                                              [](auto sum, const auto hist)
                                              {
                                                sum->Add(hist, -1);
                                                return sum;
                                              });

        Plot(*bkgSubtracted, "backgroundSubtracted", prefix, tgtname);
        Plot(*bkgScaledSubtracted, "backgroundSubtractedScaled", prefix, tgtname);
        Plot(*sidebandScaledSubtracted, "sidebandScaledSubtracted", prefix, tgtname);


        toSubtractSideband->Clone()->Write("toSubtractSideband"); // TODO: Seg fault first appears when I uncomment this line
        toSubtract->Clone()->Write("toSubtract"); // TODO: Seg fault first appears when I uncomment this line

        bkgSubtracted->Clone()->Write("bkgSubtracted"); // TODO: Seg fault first appears when I uncomment this line
        bkgScaledSubtracted->Clone()->Write("bkgScaledSubtracted"); // TODO: Seg fault first appears when I uncomment this line
        sidebandScaledSubtracted->Clone()->Write("sidebandScaledSubtracted"); // TODO: Seg fault first appears when I uncomment this line

        folded->Clone()->Write("folded"); // TODO: Seg fault first appears when I uncomment this line

        bkgSubtracted->Write("backgroundSubtracted");
        bkgScaledSubtracted->Write("backgroundSubtractedSidebandScaled");

        MCbkgScaledSubtracted->Write("MCbackgroundSubtractedSidebandScaled");

        //This is the selected signal WITH the inflated/tuned backgrounds - This is important to pass into transwarp since it has an internal mechanism
        //for background signal sepration using the migration matrix but that only works of the background passed to it matches that which we use in our
        //analysis
        MnvHist* mc_reco_tuned = (MnvHist*)SelectedSignalReco->Clone("mc_reco_tuned");
        mc_reco_tuned->Add(BackgroundWrongSign);
        mc_reco_tuned->Add(BackgroundNC);
        mc_reco_tuned->Add(BackgroundWaterTank);
        mc_reco_tuned->Add(BackgroundOtherTarget);
        mc_reco_tuned->Add(BackgroundVtxElsewhere);
        mc_reco_tuned->Add(BackgroundOther);
        mc_reco_tuned->Add(BackgroundUSPlastic);
        mc_reco_tuned->Add(BackgroundDSPlastic);
        mc_reco_tuned->Write("mc_reco_tuned");
        BackgroundUSPlastic->Write("BackgroundUSPlasticTuned");
        BackgroundDSPlastic->Write("BackgroundDSPlasticTuned");
        fluxIntReweighted->Write("fluxIntReweighted");
        flux->Write("flux");
        folded->Write("folded");
        migration->Write("migration");
        effNum->Write("effNum");
        effDenom->Write("effDenom");
        effDenom2P2H->Write("effDenom2P2H");
        effDenomDIS->Write("effDenomDIS");
        effDenomRES->Write("effDenomRES");
        effDenomQE->Write("effDenomQE");
        effDenomOther->Write("effDenomOther");
        USSidebandSignal->Write("USSidebandSignal");
        USSidebandDS->Write("USSidebandDS");
        USSidebandUS->Write("USSidebandUS");
        USSidebandOther->Write("USSidebandOther");
        DSSidebandSignal->Write("DSSidebandSignal");
        DSSidebandDS->Write("DSSidebandDS");
        DSSidebandUS->Write("DSSidebandUS");
        DSSidebandOther->Write("DSSidebandOther");
        DataUSSideband->Write("DataUSSideband");
        DataDSSideband->Write("DataDSSideband");
        DataSignal->Write("DataSignal");
        SelectedSignalReco->Write("SelectedSignalReco");
        BackgroundWrongSign->Write("BackgroundWrongSign");
        BackgroundNC->Write("BackgroundNC");
        BackgroundOther->Write("BackgroundOther");
        BackgroundWaterTank->Write("BackgroundWaterTank");
        BackgroundUSPlastic->Write("BackgroundUSPlastic");
        BackgroundDSPlastic->Write("BackgroundDSPlastic");
        MCData->Write("MCData");
        MCData2p2h->Write("MCData2p2h");
        MCDataDIS->Write("MCDataDIS");
        MCDataRES->Write("MCDataRES");
        MCDataQE->Write("MCDataQE");
        MCDataOther->Write("MCDataOther");
        outFile->Close();

    }
    catch (const std::runtime_error &e)
    {
    std::cerr << "Failed to extract a cross section for prefix " << prefix << " and target " << tgtname << " : " << e.what() << "\n";
    //return;
    // break;
    }

    delete fluxIntReweighted;
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
    delete USSidebandSignal;
    delete USSidebandDS;
    delete USSidebandUS;
    delete USSidebandOther;
    delete DSSidebandSignal;
    delete DSSidebandDS;
    delete DSSidebandUS;
    delete DSSidebandOther;
    delete DataUSSideband;
    delete DataDSSideband;
    delete DataSignal;
    delete SelectedSignalReco;
    delete BackgroundWrongSign;
    delete BackgroundNC;
    delete BackgroundOther;
    delete BackgroundWaterTank;
    delete BackgroundUSPlastic;
    delete BackgroundDSPlastic;
    delete MCData;
    delete MCData2p2h;
    delete MCDataDIS;
    delete MCDataRES;
    delete MCDataQE;
    delete MCDataOther;
}


int main(const int argc, const char **argv)
{
#ifndef NCINTEX
    ROOT::Cintex::Cintex::Enable(); // Needed to look up dictionaries for PlotUtils classes like MnvH1D
#endif

    TH1::AddDirectory(kFALSE); // Needed so that MnvH1D gets to clean up its own MnvLatErrorBands (which are TH1Ds).
    if (argc < 4)
    {
        std::cerr << "Expected 3 or more arguments, but I got " << argc - 1 << ".\n" << HELP << std::endl;
        return 1;
    }
    std::string indir = std::string(argv[1]);
    int pdg = std::stoi(argv[2]);
    std::string intgt, inptl;
    for (int i = 3; i < argc; i++)
    {
        std::cout<<"argv["<< i <<"]: " << argv[i] << std::endl;
        if (std::string(argv[i]) == "--target")
        {
            std::cout<<"Here123\n";
            if (argc < i+1) //Worth doing a more comprehensive check
            {
                std::cerr<<"Did you pass the --target flag and then not provide a target?\n";
                return -1;
            }
            else
            {
                intgt = std::string(argv[i+1]);
            }
        }
        else if (std::string(argv[i]) == "--petal")
        {
            if (argc < i+1) //Worth doing a more comprehensive check
            {
                std::cerr<<"Did you pass the --petal flag and then not provide a petal?\n";
                return -1;
            }
            else
            {
                inptl = std::string(argv[i+1]);
            }
        }

    }
    std::cout<<"intgt: " << intgt << std::endl;


    //Set sideband subtraction type
    int sidebandSubtraction = 1; //0 = No subtraction, 1 = scalar subraction, 2 = bin-by-bin subtraction (Avoid for now since it seems more fragile for 2D where some bins have a very low occupancy)
    
    std::vector<std::string> target_dirs = util::findContainingDirectories(indir, "Targets", true, true); //Should 
    if (target_dirs.size()==0) target_dirs = {"./"};
    std::vector<std::string> tracker_dirs = util::findContainingDirectories(indir, "Tracker", true, true); //Should 
    if (target_dirs.size()==0) tracker_dirs = {"./"};
    if (target_dirs.size() != tracker_dirs.size()) std::cout<<"Just letting you know that I found a different number of playlists containing Target event loops than tracker event loops\n";

    std::vector<std::string> targets;
    if (intgt == "ALL" || intgt == "all")
    {
        targets = {"1026", "1082", "2026", "2082", "3006", "3026", "3082", "4082", "5026", "5082", "6000", "Iron", "Carbon", "Lead", "WaterFull", "WaterEmpty"};
    }
    else
    {
        targets = {intgt};
    }
    std::vector<std::string> petals;
    if (inptl == "ALL" || inptl == "all")
    {
        petals = {"0","1","2","3","4","5","6","7","8","9","10","11"};
    }
    else
    {
        petals = {inptl};
    }
    std::vector<std::string> crossSectionPrefixes = {"pTmu", "pZmu", "BjorkenX", "Erecoil", "Emu"};
    //std::vector<std::string> crossSectionPrefixes = {"pTmu", "pZmu", "BjorkenX", "BjorkenY", "Erecoil", "Emu" , "beamAngle", "segment", "ANNConfidence"};
    //std::vector<std::string> crossSectionPrefixes = {"BjorkenX"};
    //std::vector<std::string> crossSectionPrefixes = {"beamAngle"};
    //std::vector<std::string> crossSectionPrefixes = {"pTmu", "pZmu"};
    //std::vector<std::string> crossSectionPrefixes = {"ANNConfidence"};

    for (const auto &prefix : crossSectionPrefixes)
    {
        bool is2D = false;
        if (prefix == "pTmu_pZmu") is2D = true;
        for (std::string &tgtname : targets)
        {
            if (is2D)
            {
                ProcessTargets<PlotUtils::MnvH2D>(prefix, tgtname, pdg, target_dirs, sidebandSubtraction);
            }
            else
            {
                ProcessTargets<PlotUtils::MnvH1D>(prefix, tgtname, pdg, target_dirs, sidebandSubtraction);
            }
        }
        //Tracker not implemented yet
        /* for (std::string &ptl : petals)
        {
            if (is2D)
            {
                ProcessTracker<PlotUtils::MnvH2D>(prefix, ptl, pdg, tracker_dirs);
            }
            else
            {
                ProcessTracker<PlotUtils::MnvH1D>(prefix, ptl, pdg, tracker_dirs);
            }
        } */
    }
    return 0;
}
