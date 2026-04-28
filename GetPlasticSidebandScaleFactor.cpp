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
    
    //Sgettin to 0 to avoid massive log spam when looping hundreds of universes, comment out if debugging
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


void CalcScaleFactor(const std::string& tgtname, 
                         int pdg, 
                         const std::vector<std::string>& dirs)
{
    PlotUtils::MnvH1D USSidebandSignal;
    PlotUtils::MnvH1D USSidebandDS;
    PlotUtils::MnvH1D USSidebandUS;
    PlotUtils::MnvH1D USSidebandOther;
    PlotUtils::MnvH1D DSSidebandSignal;
    PlotUtils::MnvH1D DSSidebandDS;
    PlotUtils::MnvH1D DSSidebandUS;
    PlotUtils::MnvH1D DSSidebandOther;
    PlotUtils::MnvH1D DataUSSideband;
    PlotUtils::MnvH1D DataDSSideband;
    PlotUtils::MnvH1D BackgroundDSPlastic;
    PlotUtils::MnvH1D BackgroundUSPlastic;

    std::cout << "Currently working on target: " << tgtname << std::endl;
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
            
            util::AddHist(USSidebandSignal,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_US_sideband_Signal"))).get(), mcscale);
            util::AddHist(USSidebandDS,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_US_sideband_DS"))).get(), mcscale);
            util::AddHist(USSidebandUS,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_US_sideband_US"))).get(), mcscale);
            util::AddHist(USSidebandOther,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_US_sideband_Other"))).get(), mcscale);
            util::AddHist(DSSidebandSignal,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_DS_sideband_Signal"))).get(), mcscale);
            util::AddHist(DSSidebandDS,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_DS_sideband_DS"))).get(), mcscale);
            util::AddHist(DSSidebandUS,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_DS_sideband_US"))).get(), mcscale);
            util::AddHist(DSSidebandOther,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_DS_sideband_Other"))).get(), mcscale);
            util::AddHist(DataUSSideband,util::GetIngredientSmart<PlotUtils::MnvH1D>(*dataFile, (std::string("segment_US_Sideband"))).get(), datascale);
            util::AddHist(DataDSSideband,util::GetIngredientSmart<PlotUtils::MnvH1D>(*dataFile, (std::string("segment_DS_Sideband"))).get(), datascale);
            util::AddHist(BackgroundDSPlastic,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_background_Downstream_Plastic_Bkg"))), mcscale);
            util::AddHist(BackgroundUSPlastic,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_background_Upstream_Plastic_Bkg"))), mcscale);

            if (d==0) //Only get POT once per playlist. E.g we don't wanna double count the POT for the same runs just because we're looking at target 2026 and 3026
            {
              mcPOT += util::GetIngredientSmart<TParameter<double>>(*mcFile, "POTUsed")->GetVal();
              std::cout<<"mcPOT current sum " << mcPOT << std::endl;
              dataPOT+= util::GetIngredientSmart<TParameter<double>>(*dataFile, "POTUsed")->GetVal();
              std::cout<<"dataPOT current sum " << dataPOT << std::endl;
            }
            std::cout<<"tmp2" << std::endl;

            dataFile->Close();
            mcFile->Close();
          }
        }

        std::cout<<"mcPOT: " << mcPOT << std::endl;
        std::cout<<"dataPOT: " << dataPOT << std::endl;

        std::string nutypestr;
        if (pdg>0) nutypestr = "Neutrino";
        else nutypestr = "Antineutrino";

        std::string outpath = nutypestr+"_"+tgtname + "_PlasticScaleFactors.root";
        auto outFile = TFile::Open(outpath.c_str(), "RECREATE");
        if (!outFile)
        {
            std::cerr << "Could not create a file called " << outpath << ".\n";
            //return;
        }

        if (isNumber(tgtname))
        {
          int tgtCode = std::stoi(tgtname);
          if (tgtCode>=7 && tgtCode <= 19)
          std::cerr << "We don't do a plastic sideband scaling or subtraction for target: " << tgtCode << ".\n";
          // No plastic sideband tuning possible for the pseudotargets
          return;
        } 
        double USScaleFactor = 1;
        double DSScaleFactor = 1;
        //************************************************
        //Sideband chisq minimisation - single scale factor across a given sideband region
        //************************************************
        PlotUtils::MnvH1D* combinedUSFixed = (PlotUtils::MnvH1D*)USSidebandSignal.Clone("combinedUSFixedMC");
        combinedUSFixed->Add(&USSidebandDS);
        combinedUSFixed->Add(&USSidebandOther);
        combinedUSFixed->Scale(dataPOT/mcPOT);

        PlotUtils::MnvH1D* combinedDSFixed = (PlotUtils::MnvH1D*)DSSidebandSignal.Clone("combinedDSFixedMC");
        combinedDSFixed->Add(&DSSidebandUS);
        combinedDSFixed->Add(&DSSidebandOther);
        combinedDSFixed->Scale(dataPOT/mcPOT);
        
        PlotUtils::MnvH1D* DSSidebandDSScaled = (PlotUtils::MnvH1D*) DSSidebandDS.Clone("DSSidebandDS_Scaled");
        PlotUtils::MnvH1D* USSidebandUSScaled = (PlotUtils::MnvH1D*) USSidebandUS.Clone("USSidebandUS_Scaled");
        DSSidebandDSScaled->Scale(dataPOT/mcPOT);
        USSidebandUSScaled->Scale(dataPOT/mcPOT);


        // Helper to fill all bins of a histogram with the scale factor 
        auto fillScaleHisto = [](TH1* hist, double val, double err = 0.0) {
            for(int i = 0; i < hist->GetNcells(); ++i) {
                hist->SetBinContent(i, val);
                hist->SetBinError(i, err);
            }
        };

        //Only need 1 bin since we're just saving a factor, we propogate these across different variable bins in Extract?DCrossSectionTargets later
        //Using MnvH1D to keep track of systematics
        PlotUtils::MnvH1D scaleFactorUS("USScaleFactor", "USScaleFactor", 1, 0.0, 1.0); 
        PlotUtils::MnvH1D scaleFactorDS("DSScaleFactor", "DSScaleFactor", 1, 0.0, 1.0); //Only need 1 bin since we're just saving a factor
        scaleFactorUS.ClearAllErrorBands();
        scaleFactorUS.Reset();
        scaleFactorUS.AddMissingErrorBandsAndFillWithCV(*combinedUSFixed); // Ensure it has MC error bands
        scaleFactorDS.ClearAllErrorBands();
        scaleFactorDS.Reset();
        scaleFactorDS.AddMissingErrorBandsAndFillWithCV(*combinedDSFixed); // Ensure it has MC error bands

        //Fitting CV
        TH1* cvDataUS = (TH1*)DataUSSideband.GetCVHistoWithStatError().Clone("cvDataUS");
        TH1* cvMCUSFixed   = (TH1*)combinedUSFixed->GetCVHistoWithStatError().Clone("cvMCUSFixed");
        TH1* cvMCUSFloat   = (TH1*)USSidebandUSScaled->GetCVHistoWithStatError().Clone("cvMCUSFloat");
        auto resultsUS = RunMinuitFit(cvDataUS, cvMCUSFixed, cvMCUSFloat);
        USScaleFactor = resultsUS[0];
        fillScaleHisto(&scaleFactorUS, USScaleFactor, resultsUS[1]);
        
        TH1* cvDataDS = (TH1*)DataDSSideband.GetCVHistoWithStatError().Clone("cvDataDS");
        TH1* cvMCDSFixed   = (TH1*)combinedDSFixed->GetCVHistoWithStatError().Clone("cvMCDSFixed");
        TH1* cvMCDSFloat   = (TH1*)DSSidebandDSScaled->GetCVHistoWithStatError().Clone("cvMCDSFloat");
        auto resultsDS = RunMinuitFit(cvDataDS, cvMCDSFixed, cvMCDSFloat);
        DSScaleFactor = resultsDS[0];
        fillScaleHisto(&scaleFactorDS, DSScaleFactor, resultsDS[1]);
        
        delete cvDataUS; delete cvMCUSFixed; delete cvMCUSFloat; delete cvDataDS; delete cvMCDSFixed; delete cvMCDSFloat;

        //Loop Over Vertical Bands
        std::vector<std::string> vertBandNames = combinedUSFixed->GetVertErrorBandNames();
        for (const std::string& bandName : vertBandNames) {
            auto bandUSFixed_MC = combinedUSFixed->GetVertErrorBand(bandName);
            auto bandDSFixed_MC = combinedDSFixed->GetVertErrorBand(bandName);
            auto bandUSFloat_MC = USSidebandUSScaled->GetVertErrorBand(bandName);
            auto bandDSFloat_MC = DSSidebandDSScaled->GetVertErrorBand(bandName);

            auto bandScale_US = scaleFactorUS.GetVertErrorBand(bandName);
            auto bandScale_DS = scaleFactorDS.GetVertErrorBand(bandName);
            for (unsigned int i = 0; i < bandUSFixed_MC->GetNHists(); ++i) {
                TH1* univFixedMC_US = bandUSFixed_MC->GetHist(i);
                TH1* univFixedMC_DS = bandDSFixed_MC->GetHist(i);
                
                TH1* univFloatMC_US = bandUSFloat_MC->GetHist(i);
                TH1* univFloatMC_DS = bandDSFloat_MC->GetHist(i);
                
                // Data does not have vertical errors, so we fit against CV Data
                TH1* univData_US = (TH1*)DataUSSideband.GetCVHistoWithStatError().Clone("univDataUS");
                TH1* univData_DS = (TH1*)DataDSSideband.GetCVHistoWithStatError().Clone("univDataDS");

                double univUSScale = RunMinuitFit(univData_US, univFixedMC_US, univFloatMC_US)[0];
                double univDSScale = RunMinuitFit(univData_DS, univFixedMC_DS, univFloatMC_DS)[0];

                fillScaleHisto(bandScale_US->GetHist(i), univUSScale);
                fillScaleHisto(bandScale_DS->GetHist(i), univDSScale);
                delete univData_US; delete univData_DS;
            }
        }

        //Loop Over Lateral Bands (Do I even have any?)
        std::vector<std::string> latBandNames = combinedUSFixed->GetLatErrorBandNames();
        for (const std::string& bandName : latBandNames) {
            auto bandUSFixed_MC = combinedUSFixed->GetLatErrorBand(bandName);
            auto bandDSFixed_MC = combinedDSFixed->GetLatErrorBand(bandName);
            auto bandUSFloat_MC = USSidebandUSScaled->GetLatErrorBand(bandName);
            auto bandDSFloat_MC = DSSidebandDSScaled->GetLatErrorBand(bandName);

            auto bandScale_US = scaleFactorUS.GetLatErrorBand(bandName);
            auto bandScale_DS = scaleFactorDS.GetLatErrorBand(bandName);

            // Lateral bands sometimes exist in data (e.g. momentum shifts applied to Data instead of MC)
            bool dataHasBandUS = DataUSSideband.HasLatErrorBand(bandName);
            bool dataHasBandDS = DataDSSideband.HasLatErrorBand(bandName);

            for (unsigned int i = 0; i < bandUSFixed_MC->GetNHists(); ++i) {
                TH1* univFixedMC_US = bandUSFixed_MC->GetHist(i);
                TH1* univFixedMC_DS = bandDSFixed_MC->GetHist(i);
                
                TH1* univFloatMC_US = bandUSFloat_MC->GetHist(i);
                TH1* univFloatMC_DS = bandDSFloat_MC->GetHist(i);

                TH1* univData_US = dataHasBandUS ? DataUSSideband.GetLatErrorBand(bandName)->GetHist(i) : (TH1*)DataUSSideband.GetCVHistoWithStatError().Clone("univDataUS");
                TH1* univData_DS = dataHasBandDS ? DataDSSideband.GetLatErrorBand(bandName)->GetHist(i) : (TH1*)DataDSSideband.GetCVHistoWithStatError().Clone("univDataDS");

                double univUSScale = RunMinuitFit(univData_US, univFixedMC_US, univFloatMC_US)[0];
                double univDSScale = RunMinuitFit(univData_DS, univFixedMC_DS, univFloatMC_DS)[0];

                fillScaleHisto(bandScale_US->GetHist(i), univUSScale);
                fillScaleHisto(bandScale_DS->GetHist(i), univDSScale);

                if (!dataHasBandUS) delete univData_US;
                if (!dataHasBandDS) delete univData_DS;
            }
        }
        TH1D h_syst_errUS = scaleFactorUS.GetTotalError(false);
        TH1D h_syst_errDS = scaleFactorDS.GetTotalError(false);
        std::cout << "US Scale: " << scaleFactorUS.GetBinContent(1) << " +/- " << scaleFactorUS.GetBinError(1) /* resultsUS[1] */ << " (Stat Err) +/- " << h_syst_errUS.GetBinContent(1) << " (Syst Err) with Minimised ChiSq: " << resultsUS[2] << std::endl;
        std::cout << "DS Scale: " << scaleFactorDS.GetBinContent(1) << " +/- " << scaleFactorDS.GetBinError(1) /* resultsDS[1] */ << " (Stat Err) +/- " << h_syst_errDS.GetBinContent(1) << " (Syst Err) with Minimised ChiSq: " << resultsDS[2] << std::endl;
        scaleFactorUS.Write("scaleFactor_US");
        scaleFactorDS.Write("scaleFactor_DS");
        //************************************************
        //End - Sideband chisq minimisation - single scale factor across a given sideband region
        //************************************************
        outFile->Close();

    }
    catch (const std::runtime_error &e)
    {
    std::cerr << "Failed to extract a cross section for target " << tgtname << " : " << e.what() << "\n";
    //return;
    // break;
    }
}


int main(const int argc, const char **argv)
{
#ifndef NCINTEX
    ROOT::Cintex::Cintex::Enable(); // Needed to look up dictionaries for PlotUtils classes like MnvH1D
#endif

    TH1::AddDirectory(kFALSE); // Needed so that MnvH1D gets to clean up its own MnvLatErrorBands (which are TH1Ds).
    if (argc < 2)
    {
        std::cerr << "Expected 3 or more arguments, but I got " << argc - 1 << ".\n" << HELP << std::endl;
        return 1;
    }
    std::string indir = std::string(argv[1]);
    int pdg = std::stoi(argv[2]);
    std::string intgt;
    for (int i = 3; i < argc; i++)
    {
        std::cout<<"argv["<< i <<"]: " << argv[i] << std::endl;
        if (std::string(argv[i]) == "--target")
        {
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

    }
    std::cout<<"intgt: " << intgt << std::endl;


    //Set sideband subtraction type
    int sidebandSubtraction = 1; //0 = No subtraction, 1 = scalar subraction, 2 = bin-by-bin subtraction (Avoid for now since it seems more fragile for 2D where some bins have a very low occupancy)
    
    std::vector<std::string> target_dirs = util::findContainingDirectories(indir, "Targets", true, true); //Should 
    if (target_dirs.size()==0) target_dirs = {"./"};

    std::vector<std::string> targets;
    if (intgt == "ALL" || intgt == "all")
    {
        targets = {"1026", "1082", "2026", "2082", "3006", "3026", "3082", "4082", "5026", "5082", "6000", "Iron", "Carbon", "Lead", "WaterFull", "WaterEmpty"};
    }
    else if (intgt.size()==0) targets = {"Iron", "Lead", "Carbon", "Water"};
    else
    {
        targets = {intgt};
    }
    for (std::string &tgtname : targets)
    {
      CalcScaleFactor(tgtname, pdg, target_dirs);
    }
    return 0;
}
