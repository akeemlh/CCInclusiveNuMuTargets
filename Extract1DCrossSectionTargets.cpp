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

// Helper to map bins of a histogram with the scale factor 
PlotUtils::MnvH1D* mapSFtoVarHist(PlotUtils::MnvH1D& h_ref, PlotUtils::MnvH1D& scaleFactor_hist)
{
  PlotUtils::MnvH1D* h_new = (PlotUtils::MnvH1D*)h_ref.Clone("h_new");
  h_new->SetTitle("Expanded MnvH1D");
  h_new->Reset(); 

  int nBins = h_new->GetNbinsX();

  // 2. Fill the Central Value (CV)
  double cv_val = scaleFactor_hist.GetBinContent(1);
  double cv_err = scaleFactor_hist.GetBinError(1);
  for (int i = 1; i <= nBins; ++i) {
      h_new->SetBinContent(i, cv_val);
      h_new->SetBinError(i, cv_err);
  }

  //Loop Over Vertical Bands
  std::vector<std::string> vertBandNames = scaleFactor_hist.GetVertErrorBandNames();
  for (const std::string& bandName : vertBandNames) {
    PlotUtils::MnvVertErrorBand*  refBand = h_new->GetVertErrorBand(bandName);
    PlotUtils::MnvVertErrorBand* sfBand = scaleFactor_hist.GetVertErrorBand(bandName);
    for (unsigned int u = 0; u < refBand->GetNHists(); ++u) {
      double uni_val = sfBand->GetHist(u)->GetBinContent(1);
      double uni_err = sfBand->GetHist(u)->GetBinError(1);
      for (int i = 1; i <= nBins; ++i) {
          refBand->GetHist(u)->SetBinContent(i, uni_val);
          refBand->GetHist(u)->SetBinError(i, uni_err);
      }
    }
  }

  //Loop Over Lateral Bands
  std::vector<std::string> latBandNames = scaleFactor_hist.GetLatErrorBandNames();
  for (const std::string& bandName : latBandNames) {
    PlotUtils::MnvLatErrorBand* refBand = h_new->GetLatErrorBand(bandName);
    PlotUtils::MnvLatErrorBand* sfBand = scaleFactor_hist.GetLatErrorBand(bandName);
    for (unsigned int u = 0; u < refBand->GetNHists(); ++u) {
      double uni_val = sfBand->GetHist(u)->GetBinContent(1);
      double uni_err = sfBand->GetHist(u)->GetBinError(1);
      for (int i = 1; i <= nBins; ++i) {
          refBand->GetHist(u)->SetBinContent(i, uni_val);
          refBand->GetHist(u)->SetBinError(i, uni_err);
      }
    }
  }
  return h_new;
};


struct Chi2Calculator {
    TH1D* dataHist;
    TH1D* mcFixedHist;
    TH1D* mcFloatHist;

    Chi2Calculator(TH1D* data, TH1D* mc_fix, TH1D* mc_float) : dataHist(data), mcFixedHist(mc_fix), mcFloatHist(mc_float) {}

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

std::vector<double> RunMinuitFit(TH1D* data, TH1D* mc_fix, TH1D* mc_float) {
    Chi2Calculator chi2Calc(data, mc_fix, mc_float);
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

/* 
template <class MnvH>
struct Chi2Calculator {
    MnvH* dataHist;
    MnvH* mcHist;

    Chi2Calculator(MnvH* data, MnvH* mc) : dataHist(data), mcHist(mc) {}

    double operator()(const double* val) const {
      double scale  = val[0];
      int ndof;
      PlotUtils::MnvPlotter plotter;
      plotter.ApplyStyle(PlotUtils::kCCQENuInclusiveStyle);
      auto tmpHist = dataHist->Clone("tmp_minimization_hist");
      double chiSq = plotter.Chi2DataMC(tmpHist, mcHist, ndof, scale, true, true);
      delete tmpHist; 
      return chiSq;
    }
}; */


bool isNumber(const std::string& str) {
  for (char c : str) {
    if (!std::isdigit(c)) {
      return false;
    }
  }
  return true;
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
  if (plotpngs) plotter.MultiPrint(&can, target + "_" + prefix + "_" + stepName + "_uncertaintySummary", "png,pdf,C");

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "Other");
  if (plotpngs) plotter.MultiPrint(&can, target + "_" + prefix + "_" + stepName + "_otherUncertainties", "png,pdf,C");
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
  if (argc != 6)
  {
    std::cerr << "Expected 4 arguments, but I got " << argc - 1 << ".\n" << HELP << std::endl;
    return 1;
  }
  const int nIterations = std::stoi(argv[1]);
  std::string indir = std::string(argv[2]);
  std::string intgt = std::string(argv[3]);
  int pdg = std::stoi(argv[4]);
  std::string var = std::string(argv[5]);

  //Set sideband subtraction type
  int sidebandSubtraction = 1; //0 = No subtraction, 1 = scalar subraction, 2 = bin-by-bin subtraction
  
  std::vector<std::string> dirs = util::findContainingDirectories(indir, "Targets", true, true);

  if (dirs.size()==0) dirs = {"./"};

  PlotUtils::FluxReweighter *frw = nullptr;

  std::vector<std::string> targets;
  if (intgt == "ALL" || intgt == "all")
  {
    targets = {"1026", "1082", "2026", "2082", "3006", "3026", "3082", "4082", "5026", "5082", "6000", "Iron", "Carbon", "Lead", "WaterFull", "WaterEmpty"};
  }
  else
  {
    targets = {intgt};
  }
  for (std::string &tgtname : targets)
  {

    std::string material = "tracker";
    if(tgtname == "3006") material = "carbon";
    else if(tgtname == "1026" || tgtname == "2026" ||tgtname == "3026" || tgtname == "5026" ) material = "iron";
    else if(tgtname == "1082" || tgtname == "2082" ||tgtname == "3082" || tgtname == "4082" || tgtname == "5082" ) material = "lead";
    else if (tgtname == "Lead") material = "lead";
    else if (tgtname == "Iron") material = "iron";
    else if (tgtname == "Carbon") material = "carbon";
    //Is it ok to just let water use the tracker flux?
    std::cout<<"material: " << material << std::endl;

    //std::vector<std::string> crossSectionPrefixes = {"pTmu", "pZmu", "BjorkenX", "Erecoil", "Emu"};
    
    //std::vector<std::string> crossSectionPrefixes = {"pTmu", "pZmu", "BjorkenX", "BjorkenY", "Erecoil", "Emu" , "beamAngle", "segment", "ANNConfidence"};
    //std::vector<std::string> crossSectionPrefixes = {"BjorkenX"};
    //std::vector<std::string> crossSectionPrefixes = {"beamAngle"};
    //std::vector<std::string> crossSectionPrefixes = {"pTmu", "pZmu"};
    //std::vector<std::string> crossSectionPrefixes = {"ANNConfidence"};
    std::vector<std::string> crossSectionPrefixes = {var};

    for (const auto &prefix : crossSectionPrefixes)
    {

      PlotUtils::MnvH1D *fluxIntReweighted = new PlotUtils::MnvH1D();

      PlotUtils::MnvH1D* flux = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* folded = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* bkgScaledSubtracted = new PlotUtils::MnvH1D();
      PlotUtils::MnvH2D* migration = new PlotUtils::MnvH2D();
      PlotUtils::MnvH1D* effCorrectedFoldedSummed = new PlotUtils::MnvH1D();
      //PlotUtils::MnvH1D* effNum = new PlotUtils::MnvH1D();
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

      PlotUtils::MnvH1D* SelectedSignalReco = new PlotUtils::MnvH1D();

      PlotUtils::MnvH1D* BackgroundWrongSign = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* BackgroundNC = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* BackgroundOther = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* BackgroundWaterTank = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* BackgroundUSPlastic = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* BackgroundDSPlastic = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* BackgroundOtherTarget = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* BackgroundVtxElsewhere = new PlotUtils::MnvH1D();

      PlotUtils::MnvH1D* BackgroundWrongSignTmp = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* BackgroundNCTmp = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* BackgroundOtherTmp = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* BackgroundWaterTankTmp = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* BackgroundUSPlasticTmp = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* BackgroundDSPlasticTmp = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* BackgroundOtherTargetTmp = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* BackgroundVtxElsewhereTmp = new PlotUtils::MnvH1D();

      PlotUtils::MnvH1D* BackgroundUSPlasticScaled = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* BackgroundDSPlasticScaled = new PlotUtils::MnvH1D();

      PlotUtils::MnvH1D* MCData = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* MCData2p2h = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* MCDataDIS = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* MCDataRES = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* MCDataQE = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* MCDataOther = new PlotUtils::MnvH1D();


      PlotUtils::MnvH1D* unfolded = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D* unfolded_tuned = new PlotUtils::MnvH1D();

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
              return 2;
            }
            double data_pot = util::GetIngredientSmart<TParameter<double>>(*dataFile, "POTUsed")->GetVal();
            auto mcFile = TFile::Open(mcpath.c_str(), "READ");
            if (!mcFile)
            {
              std::cerr << "Failed to open MC file " << mcpath.c_str() << ".\n";
              continue;
              return 3;
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

        std::string nutypestr;
        if (pdg>0) nutypestr = "Neutrino";
        else nutypestr = "Antineutrino";

        double totalMCPOTFull = 0; //Only needed for water subtractive analysis
        double totalDataPOTFull = 0;
        double totalMCPOTEmpty = 0; //Only needed for water subtractive analysis
        double totalDataPOTEMpt = 0;
        double culmulativeMCPOT = 0; //Only needed for water subtractive analysis
        double culmulativeDataPOT = 0;

        std::string CapitalisedTgt = "";
        if(tgtname == "3006") CapitalisedTgt = "Carbon";
        else if(tgtname == "1026" || tgtname == "2026" ||tgtname == "3026" || tgtname == "5026" ) CapitalisedTgt = "Iron";
        else if(tgtname == "1082" || tgtname == "2082" ||tgtname == "3082" || tgtname == "4082" || tgtname == "5082" ) CapitalisedTgt = "Lead";
        else if (tgtname == "Lead") CapitalisedTgt = "Lead";
        else if (tgtname == "Iron") CapitalisedTgt = "Iron";
        else if (tgtname == "Carbon") CapitalisedTgt = "Carbon";

        std::string scalefactorfile = indir + "/ScaleFactors/"+nutypestr+"_"+CapitalisedTgt+"_PlasticScaleFactors.root";
        auto sfFile = TFile::Open(scalefactorfile.c_str(), "READ");
        PlotUtils::MnvH1D* US_scaleFactor_tmp = util::GetIngredient<PlotUtils::MnvH1D>(*sfFile, "scaleFactor_US");
        PlotUtils::MnvH1D* DS_scaleFactor_tmp = util::GetIngredient<PlotUtils::MnvH1D>(*sfFile, "scaleFactor_DS");

        for (int d = 0; d<targetsInTgt.size(); d++)
        {
          //PlotUtils::MnvH1D* folded = new PlotUtils::MnvH1D();
          PlotUtils::MnvH1D* effNum = new PlotUtils::MnvH1D();
          PlotUtils::MnvH1D* effDenomTmp = new PlotUtils::MnvH1D();
          PlotUtils::MnvH1D* foldedTmp = new PlotUtils::MnvH1D();
          PlotUtils::MnvH2D* migrationTmp = new PlotUtils::MnvH2D();


          PlotUtils::MnvH1D* BackgroundUSPlasticTmp = new PlotUtils::MnvH1D();
          PlotUtils::MnvH1D* BackgroundDSPlasticTmp = new PlotUtils::MnvH1D();
          std::cout<<"Investigating target " << targetsInTgt[d] << " which is " << d+1 <<"/"<<targetsInTgt.size() <<" targets selected" <<std::endl;
          std::cout<<"dirs " << dirs.size()<<"\n";
          for (int c = 0; c<dirs.size(); c++)
          {
            std::cout<<"Investigating directory " << dirs[c] << " which is " << c+1 <<"/"<<dirs.size() <<" playlists identified" <<std::endl;

            std::string datapath = dirs[c] + "/runEventLoopTargetsData" + targetsInTgt[d] + ".root";
            std::string mcpath = dirs[c] + "/runEventLoopTargetsMC" + targetsInTgt[d] + ".root";

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


            //util::AddHist(*flux,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("reweightedflux_integrated")), prefix).get(), mcscale);
            util::AddHist(*foldedTmp,util::GetIngredientSmart<PlotUtils::MnvH1D>(*dataFile, (std::string("data")), prefix).get(), datascale);
            util::AddHist(*migrationTmp,util::GetIngredientSmart<PlotUtils::MnvH2D>(*mcFile, (std::string("migration")), prefix).get(), mcscale);
            util::AddHist(*effNum,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_numerator")), prefix).get(), mcscale);
            util::AddHist(*effDenomTmp,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator")), prefix).get(), mcscale);
            util::AddHist(*effDenom2P2H,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator_intChannels_2p2h")), prefix).get(), mcscale);
            util::AddHist(*effDenomDIS,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator_intChannels_DIS")), prefix).get(), mcscale);
            util::AddHist(*effDenomRES,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator_intChannels_RES")), prefix).get(), mcscale);
            util::AddHist(*effDenomQE,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator_intChannels_QE")), prefix).get(), mcscale);
            util::AddHist(*effDenomOther,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator_intChannels_Other")), prefix).get(), mcscale);
            //If doing scalar sideband tuning
            if (sidebandSubtraction == 1)
            {
              util::AddHist(*USSidebandSignal,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_US_sideband_Signal"))).get(), mcscale);
              util::AddHist(*USSidebandDS,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_US_sideband_DS"))).get(), mcscale);
              util::AddHist(*USSidebandUS,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_US_sideband_US"))).get(), mcscale);
              util::AddHist(*USSidebandOther,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_US_sideband_Other"))).get(), mcscale);
              util::AddHist(*DSSidebandSignal,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_DS_sideband_Signal"))).get(), mcscale);
              util::AddHist(*DSSidebandDS,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_DS_sideband_DS"))).get(), mcscale);
              util::AddHist(*DSSidebandUS,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_DS_sideband_US"))).get(), mcscale);
              util::AddHist(*DSSidebandOther,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_DS_sideband_Other"))).get(), mcscale);
              util::AddHist(*DataUSSideband,util::GetIngredientSmart<PlotUtils::MnvH1D>(*dataFile, (std::string("segment_US_Sideband"))).get(), datascale);
              util::AddHist(*DataDSSideband,util::GetIngredientSmart<PlotUtils::MnvH1D>(*dataFile, (std::string("segment_DS_Sideband"))).get(), datascale);
              util::AddHist(*DataSignal,util::GetIngredientSmart<PlotUtils::MnvH1D>(*dataFile, (std::string("segment_data"))).get(), datascale);
            }
            //If doing bin-by-bin sideband tuning
            if (sidebandSubtraction == 2)
            {
              util::AddHist(*USSidebandSignal,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("US_sideband_Signal")), prefix).get(), mcscale);
              util::AddHist(*USSidebandDS,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("US_sideband_DS")), prefix).get(), mcscale);
              util::AddHist(*USSidebandUS,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("US_sideband_US")), prefix).get(), mcscale);
              util::AddHist(*USSidebandOther,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("US_sideband_Other")), prefix).get(), mcscale);
              util::AddHist(*DSSidebandSignal,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("DS_sideband_Signal")), prefix).get(), mcscale);
              util::AddHist(*DSSidebandDS,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("DS_sideband_DS")), prefix).get(), mcscale);
              util::AddHist(*DSSidebandUS,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("DS_sideband_US")), prefix).get(), mcscale);
              util::AddHist(*DSSidebandOther,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("DS_sideband_Other")), prefix).get(), mcscale);
              util::AddHist(*DataUSSideband,util::GetIngredientSmart<PlotUtils::MnvH1D>(*dataFile, (std::string("US_Sideband")), prefix).get(), datascale);
              util::AddHist(*DataDSSideband,util::GetIngredientSmart<PlotUtils::MnvH1D>(*dataFile, (std::string("DS_Sideband")), prefix).get(), datascale);
              util::AddHist(*DataSignal,util::GetIngredientSmart<PlotUtils::MnvH1D>(*dataFile, (std::string("data")), prefix).get(), datascale);
            }


            util::AddHist(*SelectedSignalReco,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("selected_signal_reco")), prefix), mcscale);
            util::AddHist(*BackgroundWrongSignTmp,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_Wrong_Sign_Bkg")), prefix), mcscale);
            util::AddHist(*BackgroundNCTmp,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_NC_Bkg")), prefix), mcscale);
            util::AddHist(*BackgroundWaterTankTmp,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_Water_Tank_Bkg")), prefix), mcscale);
            util::AddHist(*BackgroundDSPlasticTmp,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_Downstream_Plastic_Bkg")), prefix), mcscale);
            util::AddHist(*BackgroundUSPlasticTmp,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_Upstream_Plastic_Bkg")), prefix), mcscale);
            util::AddHist(*BackgroundOtherTargetTmp,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_True_In_Other_Target_Bkg")), prefix), mcscale);
            util::AddHist(*BackgroundVtxElsewhereTmp,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_True_Vtx_Elsewhere_Bkg")), prefix), mcscale);
            util::AddHist(*BackgroundOtherTmp,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_Other")), prefix), mcscale);

            util::AddHist(*BackgroundWrongSign, BackgroundWrongSignTmp);
            util::AddHist(*BackgroundNC, BackgroundNCTmp);
            util::AddHist(*BackgroundWaterTank, BackgroundWaterTankTmp);
            util::AddHist(*BackgroundDSPlastic, BackgroundDSPlasticTmp);
            util::AddHist(*BackgroundUSPlastic, BackgroundUSPlasticTmp);
            util::AddHist(*BackgroundOtherTarget, BackgroundOtherTargetTmp);
            util::AddHist(*BackgroundVtxElsewhere, BackgroundVtxElsewhereTmp);
            util::AddHist(*BackgroundOther, BackgroundOtherTmp);




            util::AddHist(*MCData,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("data")), prefix), mcscale);
            util::AddHist(*MCData2p2h,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("intType_2p2h")), prefix), mcscale);
            util::AddHist(*MCDataDIS,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("intType_DIS")), prefix), mcscale);
            util::AddHist(*MCDataRES,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("intType_RES")), prefix), mcscale);
            util::AddHist(*MCDataQE,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("intType_QE")), prefix), mcscale);
            util::AddHist(*MCDataOther,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("intType_Other")), prefix), mcscale);

            if (d==0) //Only get POT once per playlist. E.g we don't wanna double count the POT for the same runs just because we're looking at target 2026 and 3026
            {
              mcPOT += util::GetIngredientSmart<TParameter<double>>(*mcFile, "POTUsed")->GetVal();
              std::cout<<"mcPOT: " << mcPOT << std::endl;
              double tempDataPOT = util::GetIngredientSmart<TParameter<double>>(*dataFile, "POTUsed")->GetVal();
              dataPOT+=tempDataPOT;
              std::cout<<"dataPOT: " << dataPOT << std::endl;

              //bool useMuCorr = true;
              //if (nupdg<14) useMuCorr == false;
              //PlotUtils::FluxReweighter frw_temp = PlotUtils::flux_reweighter(playlistUsed, nupdg, use_nue_constraint, n_flux_universes);
              //auto tempIntFlux = frw_temp.GetIntegratedFluxReweighted(nupdg, effDenom, min_energy, max_energy, useMuCorr)->Clone();

              //Scale the integrated flux for this playlist by playlist data POT for appropriate scaling when applying across different playlists
              //PlotUtils::FluxReweighter frw = PlotUtils::FluxReweighter( nupdg, use_nue_constraint, playlistUsed, PlotUtils::FluxReweighter::gen2thin, PlotUtils::FluxReweighter::g4numiv6, n_flux_universes );                    
              //auto tempIntFlux = frw.GetIntegratedFluxReweighted(nupdg, effDenomTmp, min_energy, max_energy, true)->Clone();
              //tempIntFlux->Scale(tempDataPOT);

              //PlotUtils::FluxReweighter frw = PlotUtils::FluxReweighter( nupdg, use_nue_constraint, "minervame5A", PlotUtils::FluxReweighter::gen2thin, PlotUtils::FluxReweighter::g4numiv6, n_flux_universes );                    
              //PlotUtils::FluxReweighter  frw = PlotUtils::flux_reweighter("minervame5A", nupdg, use_nue_constraint, n_flux_universes);
              //flux = frw.GetTargetFluxMnvH1D(nupdg, material, project_dir);
              //auto tempIntFlux = frw.GetIntegratedTargetFlux(nupdg, material, effDenomTmp, min_energy, max_energy, project_dir);

              frw = &PlotUtils::flux_reweighter(playlistUsed, nupdg, use_nue_constraint,n_flux_universes);
              auto tempIntFlux = frw->GetIntegratedTargetFlux(nupdg, material, effDenomTmp, min_energy, max_energy, project_dir);
              //This wont work for water will it?

              util::AddHist(*fluxIntReweighted,tempIntFlux,tempDataPOT);
            }

            dataFile->Close();
            mcFile->Close();
          }


          //Applying background scaling and subtraction
          PlotUtils::MnvH1D* DS_scaleFactor = mapSFtoVarHist(*effDenomTmp, *DS_scaleFactor_tmp);
          PlotUtils::MnvH1D* US_scaleFactor = mapSFtoVarHist(*effDenomTmp, *US_scaleFactor_tmp);
          //Plot(*DS_scaleFactor, "DS_scaleFactor_"+targetsInTgt[d], prefix, tgtname);
          //Plot(*US_scaleFactor, "US_scaleFactor_"+targetsInTgt[d], prefix, tgtname);

          //Plot(*BackgroundUSPlasticTmp, "BackgroundUSPlasticTmp_Prescaling_"+targetsInTgt[d], prefix, tgtname);
          //Plot(*BackgroundDSPlasticTmp, "BackgroundDSPlasticTmp_Prescaling_"+targetsInTgt[d], prefix, tgtname);
          BackgroundUSPlasticTmp->Multiply(BackgroundUSPlasticTmp, US_scaleFactor);
          BackgroundDSPlasticTmp->Multiply(BackgroundDSPlasticTmp, DS_scaleFactor);
          //Plot(*BackgroundUSPlasticTmp, "BackgroundUSPlasticTmp_"+targetsInTgt[d], prefix, tgtname);
          //Plot(*BackgroundDSPlasticTmp, "BackgroundDSPlasticTmp_"+targetsInTgt[d], prefix, tgtname);

          util::AddHist(*BackgroundDSPlasticScaled, BackgroundDSPlasticTmp);
          util::AddHist(*BackgroundUSPlasticScaled, BackgroundUSPlasticTmp);
          
          //Plot(*BackgroundUSPlastic, "BackgroundUSPlastic", prefix, tgtname);
          //Plot(*BackgroundDSPlastic, "BackgroundDSPlastic", prefix, tgtname);
          Plot(*foldedTmp, "foldedTmp_"+targetsInTgt[d], prefix, tgtname);

          std::vector<PlotUtils::MnvH1D*> sidebandScaled = {BackgroundUSPlasticTmp, BackgroundDSPlasticTmp};
          std::vector<PlotUtils::MnvH1D*> backgrounds = { BackgroundWrongSignTmp, BackgroundNCTmp, BackgroundOtherTmp, BackgroundWaterTankTmp, BackgroundOtherTargetTmp, BackgroundVtxElsewhereTmp};
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
          auto bkgSubtracted = std::accumulate(backgrounds.begin(), backgrounds.end(), foldedTmp->Clone(),
                                                [mcPOT, dataPOT](auto sum, const auto hist)
                                                {
                                                  std::cout << "Subtracting " << hist->GetName() << " scaled by " << -dataPOT / mcPOT << " from " << sum->GetName() << "\n";
                                                  sum->Add(hist, -dataPOT / mcPOT);
                                                  return sum;
                                                });
          Plot(*bkgSubtracted, "bkgSubtracted_"+targetsInTgt[d], prefix, tgtname);

          //remove scaled sideband
          auto sidebandScaledSubtracted = std::accumulate(sidebandScaled.begin(), sidebandScaled.end(), foldedTmp->Clone(),
                                                [mcPOT, dataPOT](auto sum, const auto hist)
                                                {
                                                  std::cout << "Subtracting " << hist->GetName() << " scaled by " << -dataPOT / mcPOT << " from " << sum->GetName() << "\n";
                                                  sum->Add(hist, -dataPOT / mcPOT);
                                                  return sum;
                                                });
          //Both background subtracted and sideband scaled
          auto bkgScaledSubtractedTmp = std::accumulate(sidebandScaled.begin(), sidebandScaled.end(), bkgSubtracted->Clone(),
                                                [mcPOT, dataPOT](auto sum, const auto hist)
                                                {
                                                  std::cout << "Subtracting " << hist->GetName() << " scaled by " << -dataPOT / mcPOT << " from " << sum->GetName() << "\n";
                                                  sum->Add(hist, -dataPOT / mcPOT);
                                                  return sum;
                                                });
          util::AddHist(*bkgScaledSubtracted, bkgScaledSubtractedTmp);
          //End - Applying background scaling and subtraction
          // Basing my unfolding procedure for a differential cross section on Alex's MINERvA 101 talk at https://minerva-docdb.fnal.gov/cgi-bin/private/RetrieveFile?docid=27438&filename=whatsACrossSection.pdf&version=1
          // d'Aogstini unfolding
          auto unfoldedTmp = UnfoldHist(bkgSubtracted, migrationTmp, nIterations);
          if (!unfoldedTmp)
            throw std::runtime_error(std::string("Failed to unfold ") + bkgSubtracted->GetName() + " using " + migration->GetName());
          unfoldedTmp->Clone()->Write("unfolded"); // TODO: Seg fault first appears when I uncomment this line
          Plot(*bkgScaledSubtractedTmp, "bkgScaledSubtractedTmp_"+targetsInTgt[d], prefix, tgtname);
          auto unfolded_tunedTmp = UnfoldHist(bkgScaledSubtractedTmp, migrationTmp, nIterations);
          if (!unfolded_tunedTmp)
            throw std::runtime_error(std::string("Failed to unfold ") + bkgScaledSubtracted->GetName() + " using " + migration->GetName());
          unfolded_tunedTmp->Clone()->Write("unfolded_sidebandTuned"); // TODO: Seg fault first appears when I uncomment this line
          std::cout << "Survived writing the unfolded histogram.\n"
                    << std::flush; // This is evidence that the problem is on the final file Write() and not unfolded->Clone()->Write().
          // End - d'Aogstini unfolding

          //Efficiency correcting
          effNum->Divide(effNum, effDenomTmp, 1.0, 1.0, "B"); // Only the 2 parameter version of MnvH1D::Divide()
          unfoldedTmp->Divide(unfoldedTmp, effNum);
          Plot(*unfolded_tunedTmp, "unfolded_tunedTmp_precorrection_"+targetsInTgt[d], prefix, tgtname);
          unfolded_tunedTmp->Divide(unfolded_tunedTmp, effNum);  //Only the 2 parameter version of MnvH1D::Divide() //handles systematics correctly.
          
          Plot(*unfolded_tunedTmp, "unfolded_tunedTmp_"+targetsInTgt[d], prefix, tgtname);
          Plot(*unfolded, "unfolded_"+targetsInTgt[d], prefix, tgtname);

          util::AddHist(*unfolded, unfoldedTmp);
          util::AddHist(*unfolded_tuned, unfolded_tunedTmp);
          util::AddHist(*effDenom, effDenomTmp);
          util::AddHist(*folded, foldedTmp);
          util::AddHist(*migration, migrationTmp);
          //End - Efficiency correcting
          delete unfoldedTmp, unfolded_tunedTmp, effNum, effDenomTmp, BackgroundUSPlasticTmp, BackgroundDSPlasticTmp;
          //delete folded;
        }

        std::cout<<"mcPOT: " << mcPOT << std::endl;
        std::cout<<"dataPOT: " << dataPOT << std::endl;

        //Normalising the integrated flux across different playlists by dataPOT
        fluxIntReweighted->Scale(1/dataPOT);
        Plot(*effCorrectedFoldedSummed, "EfficiencyCorrectedFoldedSummed", prefix, tgtname);
        
        auto simEventRate = effDenom->Clone(); // Make a copy for later
        auto simEventRate2P2H = effDenom2P2H->Clone(); // Make a copy for later
        auto simEventRateDIS = effDenomDIS->Clone(); // Make a copy for later
        auto simEventRateRES = effDenomRES->Clone(); // Make a copy for later
        auto simEventRateQE = effDenomQE->Clone(); // Make a copy for later
        auto simEventRateOther = effDenomOther->Clone(); // Make a copy for later
        // There are no error bands in the data, but I need somewhere to put error bands on the results I derive from it.
        effCorrectedFoldedSummed->AddMissingErrorBandsAndFillWithCV(*migration);


        //Plot(*bkgSubtracted, "backgroundSubtracted", prefix, tgtname);
        //Plot(*bkgScaledSubtracted, "backgroundSubtractedScaled", prefix, tgtname);
        //Plot(*sidebandScaledSubtracted, "sidebandScaledSubtracted", prefix, tgtname);

        auto outFile = TFile::Open((nutypestr+"_"+tgtname +"_" + prefix + "_crossSection.root").c_str(), "RECREATE");
        if (!outFile)
        {
          std::cerr << "Could not create a file called " << (nutypestr+"_"+tgtname +"_"+ prefix + "_crossSection.root") << ".  Does it already exist?\n";
          return 5;
        }

        //toSubtractSideband->Clone()->Write("toSubtractSideband"); // TODO: Seg fault first appears when I uncomment this line
        //toSubtract->Clone()->Write("toSubtract"); // TODO: Seg fault first appears when I uncomment this line

        //bkgSubtracted->Clone()->Write("bkgSubtracted"); // TODO: Seg fault first appears when I uncomment this line
        //bkgScaledSubtracted->Clone()->Write("bkgScaledSubtracted"); // TODO: Seg fault first appears when I uncomment this line
        //sidebandScaledSubtracted->Clone()->Write("sidebandScaledSubtracted"); // TODO: Seg fault first appears when I uncomment this line

        //effCorrectedFoldedSummed->Clone()->Write("effCorrectedFoldedSummed"); // TODO: Seg fault first appears when I uncomment this line

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
        migration->Write("Migration");
        //bkgSubtracted->Write("backgroundSubtracted");
        //bkgScaledSubtracted->Write("backgroundSubtractedSidebandScaled");


        /* Plot(*effNum, "efficiencyNumerator", prefix, tgtname);
        Plot(*effDenom, "efficiencyDenominator", prefix, tgtname);
        effNum->Divide(effNum, effDenom); // Only the 2 parameter version of MnvH1D::Divide()
                                          // handles systematics correctly.
        Plot(*effNum, "efficiency", prefix, tgtname);

        unfolded->Divide(unfolded, effNum);
        Plot(*unfolded, "efficiencyCorrected", prefix, tgtname);

        unfolded_tuned->Divide(unfolded_tuned, effNum); */


        Plot(*folded, "folded", prefix, tgtname);
        Plot(*effDenom, "efficiencyDenominator", prefix, tgtname);
        Plot(*effDenom2P2H, "efficiencyDenominator2P2H", prefix, tgtname);
        Plot(*simEventRate, "simEventRate", prefix, tgtname);


        Plot(*unfolded_tuned, "efficiencyCorrected_sidebandTuned", prefix, tgtname);

        // double nnucleons = nNucleons->GetVal()/numMergedPlaylists;
        double nnucleons = 0;

        const size_t base = tgtname.find("Scaled");
        if (base != std::string::npos)
        {
          tgtname= tgtname.substr(0, base);
        }

        double nnucleonsData;
        PlotUtils::TargetUtils targetInfo;
        targetInfo.SetDistToDivCut(25);

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
        else if (isNumber(tgtname))
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
          if (tgtNum<7 && tgtCode>1000)
          {
            nnucleons = targetInfo.GetPassiveTargetNNucleons(tgtNum, tgtMat, true);
            nnucleonsData = targetInfo.GetPassiveTargetNNucleons(tgtNum, tgtMat, false);
          }
          else if (tgtCode>=7 && tgtCode <= 19)
          {
            std::cout<<"tgtCode: " << tgtCode <<std::endl;
            if (tgtCode==7) nnucleons = targetInfo.GetTrackerNNucleons( 7, true); 
            if (tgtCode==8) nnucleons = targetInfo.GetTrackerNNucleons( 6, true); 
            if (tgtCode==9) nnucleons = targetInfo.GetTrackerNNucleons( 6, true); 
            if (tgtCode==10) nnucleons = targetInfo.GetTrackerNNucleons( 6, true); 
            if (tgtCode==11) nnucleons = targetInfo.GetTrackerNNucleons( 6, true); 
            if (tgtCode==12) nnucleons = targetInfo.GetTrackerNNucleons( 2, true); 
            if (tgtCode==7) nnucleonsData = targetInfo.GetTrackerNNucleons( 7, false); 
            if (tgtCode==8) nnucleonsData = targetInfo.GetTrackerNNucleons( 6, false); 
            if (tgtCode==9) nnucleonsData = targetInfo.GetTrackerNNucleons( 6, false); 
            if (tgtCode==10) nnucleonsData = targetInfo.GetTrackerNNucleons( 6, false); 
            if (tgtCode==11) nnucleonsData = targetInfo.GetTrackerNNucleons( 6, false); 
            if (tgtCode==12) nnucleonsData = targetInfo.GetTrackerNNucleons( 2, false); 
            if (tgtCode==13) nnucleons = targetInfo.GetTrackerNNucleons( 8, true); 
            if (tgtCode==13) nnucleonsData = targetInfo.GetTrackerNNucleons( 8, false); 
            if (tgtCode==14) nnucleons = targetInfo.GetTrackerNNucleons( 24, true); 
            if (tgtCode==14) nnucleonsData = targetInfo.GetTrackerNNucleons( 24, false); 
            if (tgtCode==15) nnucleons = targetInfo.GetTrackerNNucleons( 24, true); 
            if (tgtCode==15) nnucleonsData = targetInfo.GetTrackerNNucleons( 24, false); 
            if (tgtCode==16) nnucleons = targetInfo.GetTrackerNNucleons( 24, true); 
            if (tgtCode==16) nnucleonsData = targetInfo.GetTrackerNNucleons( 24, false); 
            if (tgtCode==17) nnucleons = targetInfo.GetTrackerNNucleons( 24, true); 
            if (tgtCode==17) nnucleonsData = targetInfo.GetTrackerNNucleons( 24, false); 
            if (tgtCode==18) nnucleons = targetInfo.GetTrackerNNucleons( 24, true); 
            if (tgtCode==18) nnucleonsData = targetInfo.GetTrackerNNucleons( 24, false); 
            if (tgtCode==19) nnucleons = targetInfo.GetTrackerNNucleons( 18, true); 
            if (tgtCode==19) nnucleonsData = targetInfo.GetTrackerNNucleons( 18, false); 
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


        if (tgtname == "Water")
        {
          dataPOT = waterFilledPOTData;
          mcPOT = waterFilledPOTMC;
        }

        std::string tgtString = "";
        //Commented the below lines out to make this script more generalisable to the ANu case as well
        //if (tgtname == "WaterEmpty") tgtString = "Water (Empty Playlists - Approx 63% of ME POT)";
        //else if (tgtname == "WaterFull") tgtString = "Water (Filled Playlists - Approx 37% of ME POT)";
        if (tgtname == "WaterEmpty") tgtString = "Water (Empty Playlists)";
        else if (tgtname == "WaterFull") tgtString = "Water (Filled Playlists)";
        else if (tgtname == "Water") tgtString = "Water (Filled - Empty)";
        else if (tgtname == "Iron") tgtString = "Iron";
        else if (tgtname == "Carbon") tgtString = "Carbon";
        else if (tgtname == "Lead") tgtString = "Lead";
        else
        {
          int tgtInt = std::stoi(tgtname);
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
        fluxIntReweighted->Write("fluxIntReweighted");
        auto crossSection = normalize(unfolded, fluxIntReweighted, nnucleonsData, dataPOT);
        Plot(*crossSection, "crossSection", prefix, tgtname);
        crossSection->Write("crossSection");
        auto crossSectionTuned = normalize(unfolded_tuned, fluxIntReweighted, nnucleonsData, dataPOT);
        Plot(*crossSectionTuned, "crossSectionTuned", prefix, tgtname);
        crossSectionTuned->Write("crossSectionTuned");
        simEventRate->Write("simulatedEventRate");
        fluxIntReweighted->Write("fluxIntReweighted");

        MCData->Write("MCData");
        MCData2p2h->Write("MCData2p2h");
        MCDataDIS->Write("MCDataDIS");
        MCDataRES->Write("MCDataRES");
        MCDataQE->Write("MCDataQE");
        MCDataOther->Write("MCDataOther");
        bkgScaledSubtracted->Write("EventRateBackgroundSubtracted");
        folded->Write("EventRate");

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



        { //Plotting event rates
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kDefaultStyle);
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

          auto tempUS = BackgroundUSPlastic->Clone();
          auto tempDS = BackgroundDSPlastic->Clone();
          auto tempUSScaled = BackgroundUSPlasticScaled->Clone();
          auto tempDSScaled = BackgroundDSPlasticScaled->Clone();
          auto tempTank = BackgroundWaterTank->Clone();
          auto tempWS = BackgroundWrongSign->Clone();
          auto tempNC = BackgroundNC->Clone();
          auto tempBackgroundOther = BackgroundOther->Clone();
          auto tempOtherTgt = BackgroundOtherTarget->Clone();
          auto tempVtxElsewhere = BackgroundVtxElsewhere->Clone();
          tempBackgroundOther->Add(BackgroundWaterTank);
          auto tempsignal = SelectedSignalReco->Clone();

          auto tempevrate = folded->Clone();
          auto tempevrateScaled = bkgScaledSubtracted->Clone();
          
          auto tempsim = MCData->Clone();

          Plot(*tempUS, "tempUS", prefix, tgtname);
          Plot(*tempDS, "tempDS", prefix, tgtname);

          tempUS->SetTitle("Upstream Plastic");
          tempDS->SetTitle("Downstream Plastic");
          //tempTank->SetTitle("Water Tank");

          tempWS->SetTitle("Wrong Sign");
          tempNC->SetTitle("Neutral Current");
          tempBackgroundOther->SetTitle("Other");
          tempOtherTgt->SetTitle("Different Target");
          tempVtxElsewhere->SetTitle("True elsewhere");
          tempBackgroundOther->Add(BackgroundWaterTank);
          tempsignal->SetTitle("Signal");


          temp2p2h->Scale(1e-3);
          tempdis->Scale(1e-3);
          tempres->Scale(1e-3);
          tempqe->Scale(1e-3);
          tempChannelOther->Scale(1e-3);

          tempUS->Scale(1e-3);
          tempDS->Scale(1e-3);
          tempUSScaled->Scale(1e-3);
          tempDSScaled->Scale(1e-3);
          //tempTank->Scale(1e-3);
          tempWS->Scale(1e-3);
          tempNC->Scale(1e-3);
          tempBackgroundOther->Scale(1e-3);
          tempOtherTgt->Scale(1e-3);
          tempVtxElsewhere->Scale(1e-3);

          Plot(*tempUS, "tempUS0", prefix, tgtname);
          Plot(*tempDS, "tempDS0", prefix, tgtname);

          tempsim->Scale(1e-3);
          tempsignal->Scale(1e-3);


          tempevrate->Scale(1e-3);
          tempevrateScaled->Scale(1e-3);

          if (tgtname == "Water")
          {
            dataPOT = waterFilledPOTData;
            mcPOT = waterFilledPOTMC;
          }

          temp2p2h->Scale(dataPOT/mcPOT);
          tempdis->Scale(dataPOT/mcPOT);
          tempres->Scale(dataPOT/mcPOT);
          tempqe->Scale(dataPOT/mcPOT);
          tempChannelOther->Scale(dataPOT/mcPOT);

          Plot(*tempUS, "tempUS1", prefix, tgtname);
          Plot(*tempDS, "tempDS1", prefix, tgtname);

          std::cout<<"dataPOT: " << dataPOT << " mcPOT " << mcPOT << std::endl;
          tempUS->Scale(dataPOT/mcPOT);
          tempDS->Scale(dataPOT/mcPOT);
          //tempTank->Scale(dataPOT/mcPOT);

          tempWS->Scale(dataPOT/mcPOT);
          tempNC->Scale(dataPOT/mcPOT);
          tempBackgroundOther->Scale(dataPOT/mcPOT);
          tempOtherTgt->Scale(dataPOT/mcPOT);
          tempVtxElsewhere->Scale(dataPOT/mcPOT);
          tempsignal->Scale(dataPOT/mcPOT);

          tempsim->Scale(dataPOT/mcPOT);
          tempevrate->Scale(dataPOT/dataPOT);
          //tempevrateScaled->Scale(dataPOT/dataPOT);

          TObjArray* channelArr = new TObjArray();
          channelArr->Add(temp2p2h);
          channelArr->Add(tempdis);
          channelArr->Add(tempres);
          channelArr->Add(tempqe);
          channelArr->Add(tempChannelOther);

          Plot(*tempUS, "tempUS2", prefix, tgtname);
          Plot(*tempDS, "tempDS2", prefix, tgtname);

          TObjArray* bkgArr = new TObjArray();
          bkgArr->Add(tempUS);
          bkgArr->Add(tempDS);
          //bkgArr->Add(tempTank);
          bkgArr->Add(tempWS);
          bkgArr->Add(tempNC);
          bkgArr->Add(tempBackgroundOther);
          bkgArr->Add(tempOtherTgt);
          bkgArr->Add(tempVtxElsewhere);
          bkgArr->Add(tempsignal);

          std::string tmp; 
          std::string y_title="Events#times10^{3}";
          double xmax = -1111; //Default value
          std::string xaxislabel = "";
          if (prefix=="pTmu")
          {
            tmp = "Muon p_{T}";
            xaxislabel = "p_{T, #mu} [GeV/c]";
            y_title += " /(GeV/c)";
            xmax = 13;
            tempsim->GetXaxis()->SetRange(1, xmax);
            tempsignal->GetXaxis()->SetRange(1, xmax);
          }
          else if (prefix=="pZmu")
          {
            tmp = "Muon p_{Z}";
            xaxislabel = "p_{Z, #mu} [GeV/c]";
            y_title += " /(GeV/c)";
            xmax = 13;
            tempsim->GetXaxis()->SetRange(1, xmax);
            tempsignal->GetXaxis()->SetRange(1, xmax);
          }
          else if (prefix=="BjorkenX")
          {
            tmp = "Bjorken X";
            xaxislabel = "Bjorken X";
            y_title += " per unit X";
            xmax = 7;
            tempsim->GetXaxis()->SetRange(1, xmax);
            tempsignal->GetXaxis()->SetRange(1, xmax);
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
            tempsim->GetXaxis()->SetRange(1, xmax);
            tempsignal->GetXaxis()->SetRange(1, xmax);
          } 

          std::string title = tmp + " " + tgtString;
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;

          plotter.DrawErrorSummary(tempsim);
          plotter.MultiPrint(&can, tgtname+"_"+prefix+"_EvRateErrors", "png,pdf,C");
          can.Clear();


          // Adjust pad heights manually (top bigger, bottom smaller)
          TPad *pad1 = new TPad("pad1","pad1",0,0.35,1,1.0); // top pad
          TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.35); // bottom pad
          pad1->SetBottomMargin(0); // remove x-axis label space for top
          pad2->SetTopMargin(0);
          pad2->SetBottomMargin(0.25);
          pad1->Draw();
          pad2->Draw();

          if (prefix=="BjorkenX")
          {
            pad1->SetLogx();
            pad2->SetLogx();
          }

          double defaultmin = plotter.axis_minimum;
          plotter.legend_offset_x = 0.1325;
          plotter.legend_border_size = 0;
          plotter.legend_text_size = 0.03;

          plotter.mc_line_width=0; //Why??
          if (tgtname == "Water")
          {
            std::stringstream POTstrFilled, POTstrEmpty;
            POTstrFilled << std::fixed << std::scientific << std::setprecision(2) << waterFilledPOTData;
            POTstrEmpty << std::fixed << std::scientific << std::setprecision(2) << waterEmptyPOTData;
            std::cout << "waterFilledPOTData: "<< waterFilledPOTData <<std::endl;
            std::cout << "waterEmptyPOTData: "<< waterEmptyPOTData <<std::endl;
            std::string potstrFilled = std::string("Water Filled POT ")+POTstrFilled.str();
            std::string potstrEmpty = std::string("Water Empty POT: ")+POTstrEmpty.str();// + " (Scaled to match filled)";
            plotter.AddPlotLabel(potstrFilled.c_str() , 0.18, 0.93, 0.025, 1, 42, 13, 0);
            plotter.AddPlotLabel(potstrEmpty.c_str() , 0.18, 0.90, 0.025, 1, 42, 13, 0);
          }
          else
          {
            std::stringstream POTstr;
            POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOT;
            std::string potstr = std::string("POT Normalised: ")+POTstr.str();
            plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.93, 0.025, 1, 42, 13, 0);
          }
          plotter.WritePreliminary( 0.5, 0.15, 0.035, true);
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 25, 1, 43, 22, 0);
          plotter.AddPlotLabel(y_title.c_str() , 0.08, 0.9, 0.03, 1, 42, 33, 90);
          plotter.AddPlotLabel("Data/MC" , 0.08, 0.3, 0.03, 1, 42, 33, 90);
          plotter.AddPlotLabel(xaxislabel.c_str() , 0.9, 0.08, 0.03, 1, 42, 33, 0);

          //Interaction Type
          pad1->cd(); 
          plotter.axis_minimum=0.001;

          int channelColArr[5] = {TColor::GetColor("#d2c271"), TColor::GetColor("#CA5454"), TColor::GetColor("#40B0A6"), TColor::GetColor("#9285d5"), TColor::GetColor("#c2c2c2")};
          plotter.DrawDataStackedMC(tempevrate, channelArr, channelColArr, 1, "TR", "Data (Stat. Only)", 1001, "", "", false, 1, xmax);
          pad2->cd();
          plotter.axis_minimum=defaultmin;
          plotter.DrawDataMCRatio(tempevrate, tempsim, 1, true, true, 0.5, 1.5, "");
          can.SetLogx();
          pad1->cd();
          tempevrate->Draw("AXIS SAME");
          can.cd();

          plotter.MultiPrint(&can, tgtname+"_"+prefix+"_EvRateByIntType", "png,pdf,C");
          //Are these duplicates??
          pad1->cd(); 
          plotter.axis_minimum=0.001;
          plotter.DrawDataStackedMC(tempevrate, channelArr, channelColArr, 1, "TR", "Data (Stat. Only)", 1001, "", "", false, 1, xmax);
          pad2->cd();
          plotter.axis_minimum=defaultmin;
          plotter.DrawDataMCRatio(tempevrate, tempsim, 1, true, true, 0.5, 1.5, "");
          can.SetLogx();
          pad1->cd();
          tempevrate->Draw("AXIS SAME"); //!!??
          can.cd();
          plotter.MultiPrint(&can, tgtname+"_"+prefix+"_EvRateByIntType", "png,pdf,C");

          //Backgrounds
          pad1->cd(); 
          plotter.axis_minimum=0.001;
          int bkgColArr[9] = {TColor::GetColor("#ee6677"), TColor::GetColor("#228833"),/*  TColor::GetColor("#4477aa"), */ TColor::GetColor("#ccbb44"), TColor::GetColor("#66ccee"), TColor::GetColor("#aa3377"), TColor::GetColor("#bbbbbb"), TColor::GetColor("#CA5454"), TColor::GetColor("#40B0A6")};
          plotter.DrawDataStackedMC(tempevrate, bkgArr, bkgColArr, 1, "TR", "Data (Stat. Only)", 1001, "", "", false, 1, xmax);
          pad2->cd();
          plotter.axis_minimum=defaultmin;
          plotter.DrawDataMCRatio(tempevrate, tempsim, 1, true, true, 0.5, 1.5, "");
          can.SetLogx();
          pad1->cd();
          tempevrate->Draw("AXIS SAME");
          can.cd();
          plotter.MultiPrint(&can, tgtname+"_"+prefix+"_EvRateByBackgroundsSidebandScaled", "png,pdf,C");
          can.SetLogx(0);
          can.SetLogy(0);

          //Backgrounds With Sideband Scale Applied
          TObjArray* bkgArrScaled = new TObjArray();
          bkgArrScaled->Add(tempUSScaled);
          bkgArrScaled->Add(tempDSScaled);
          //bkgArrScaled->Add(tempTank);
          bkgArrScaled->Add(tempWS);
          bkgArrScaled->Add(tempNC);
          bkgArrScaled->Add(tempBackgroundOther);
          bkgArrScaled->Add(tempOtherTgt);
          bkgArrScaled->Add(tempVtxElsewhere);
          bkgArrScaled->Add(tempsignal);

          auto tmpScaledMCSummed = tempsignal->Clone();
          tmpScaledMCSummed->Add(tempUSScaled);
          tmpScaledMCSummed->Add(tempDSScaled);
          tmpScaledMCSummed->Add(tempWS);
          tmpScaledMCSummed->Add(tempNC);
          tmpScaledMCSummed->Add(tempBackgroundOther);
          tmpScaledMCSummed->Add(tempOtherTgt);
          tmpScaledMCSummed->Add(tempVtxElsewhere);

          pad1->cd(); 
          plotter.axis_minimum=0.001;
          plotter.DrawDataStackedMC(tempevrate, bkgArrScaled, bkgColArr, 1, "TR", "Data (Stat. Only)", 1001, "", "", false, 1, xmax);
          pad2->cd();
          plotter.axis_minimum=defaultmin;
          plotter.DrawDataMCRatio(tempevrate, tmpScaledMCSummed, 1, true, true, 0.5, 1.5, "");
          can.SetLogx();
          pad1->cd();
          tempevrate->Draw("AXIS SAME");
          can.cd();
          plotter.MultiPrint(&can, tgtname+"_"+prefix+"_EvRateByBackgrounds", "png,pdf,C");
          can.SetLogx(0);
          can.SetLogy(0);

          //Background subtracted
          
          //Interaction Type
          pad1->cd(); 
          plotter.axis_minimum=0.001;

          plotter.DrawDataStackedMC(tempevrateScaled, channelArr, channelColArr, 1, "TR", "Data (Stat. Only)", 1001, "", "", false, 1, xmax);
          pad2->cd();
          plotter.axis_minimum=defaultmin;
          plotter.DrawDataMCRatio(tempevrateScaled, tempsim, 1, true, true, 0.5, 1.5, "");
          can.SetLogx();
          pad1->cd();
          tempevrateScaled->Draw("AXIS SAME");
          can.cd();
          plotter.MultiPrint(&can, tgtname+"_"+prefix+"_EvRateByIntType_BackgroundSubtracted", "png,pdf,C");

          pad1->cd(); 
          plotter.axis_minimum=0.001;
          plotter.DrawDataStackedMC(tempevrateScaled, channelArr, channelColArr, 1, "TR", "Data (Stat. Only)", 1001, "", "", false, 1, xmax);
          pad2->cd();
          plotter.axis_minimum=defaultmin;
          plotter.DrawDataMCRatio(tempevrateScaled, tempsim, 1, true, true, 0.5, 1.5, "");
          can.SetLogx();
          pad1->cd();
          tempevrateScaled->Draw("AXIS SAME"); //!!??
          can.cd();
          plotter.MultiPrint(&can, tgtname+"_"+prefix+"_EvRateByIntTypeSidebandScaled_BackgroundSubtracted", "png,pdf,C");

          //Backgrounds
          TObjArray* bkgArr2 = new TObjArray();
          bkgArr2->Add(tempsignal);

          pad1->cd(); 
          plotter.axis_minimum=0.001;
          int bkgColArr2[7] = {TColor::GetColor("#ee6677"), TColor::GetColor("#228833")};
          plotter.DrawDataStackedMC(tempevrateScaled, bkgArr2, bkgColArr2, 1, "TR", "Data (Stat. Only)", 1001, "", "", false, 1, xmax);
          pad2->cd();
          plotter.axis_minimum=defaultmin;
          plotter.DrawDataMCRatio(tempevrateScaled, tempsignal, 1, true, true, 0.5, 1.5, "");
          can.SetLogx();
          pad1->cd();
          tempevrateScaled->Draw("AXIS SAME");
          can.cd();
          plotter.MultiPrint(&can, tgtname+"_"+prefix+"_EvRateByBackgrounds_BackgroundSubtracted", "png,pdf,C");
          can.SetLogx(0);
          can.SetLogy(0);



          TObjArray* bkgArr3 = new TObjArray();
          bkgArr3->Add(tempWS);
          bkgArr3->Add(tempNC);
          bkgArr3->Add(tempBackgroundOther);
          bkgArr3->Add(tempOtherTgt);
          bkgArr3->Add(tempVtxElsewhere);
          bkgArr3->Add(tempsignal);

          //bkgArr3->Add(tempsignal);

          pad1->cd(); 
          plotter.axis_minimum=0.001;
          int bkgColArr3[7] = {TColor::GetColor("#ee6677"), TColor::GetColor("#228833"),/*  TColor::GetColor("#4477aa"), */ TColor::GetColor("#ccbb44"), TColor::GetColor("#66ccee"), TColor::GetColor("#aa3377"), TColor::GetColor("#bbbbbb")};
          plotter.DrawDataStackedMC(tempevrateScaled, bkgArr3, bkgColArr3, 1, "TR", "Data (Stat. Only)", 1001, "", "", false, 1, xmax);
          pad2->cd();
          plotter.axis_minimum=defaultmin;
          plotter.DrawDataMCRatio(tempevrateScaled, tempsignal, 1, true, true, 0.5, 1.5, "");
          can.SetLogx();
          pad1->cd();
          tempWS->Draw("AXIS SAME");
          can.cd();
          plotter.MultiPrint(&can, tgtname+"_"+prefix+"_EvRateByBackgrounds_SidebandScaledAndSubtracted", "png,pdf,C");
          can.SetLogx(0);
          can.SetLogy(0);

          
        } //End plotting event rates

        auto simulatedCrossSection = normalize(simEventRate, fluxIntReweighted, nnucleons, mcPOT);
        Plot(*simulatedCrossSection, "simulatedCrossSection", prefix, tgtname);
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
        TParameter<double> EmptyPOTData = TParameter<double>("dataPOTWaterEmpty", waterEmptyPOTData);
        TParameter<double> FilledPOTData = TParameter<double>("dataPOTWaterFilled", waterFilledPOTData);
        potMC.Write();
        potData.Write();
        EmptyPOTData.Write();
        FilledPOTData.Write();
        
        //Plotting stuff - Should move this elsewhere long-term

        /* { //Plotting cross sections
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
          if (prefix=="pTmu")
          {
            tmp = "Muon p_{T}";
            y_title="d#sigma/dp_{t} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            tempxsec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix=="pZmu")
          {
            tmp = "Muon p_{Z}";
            y_title="d#sigma/dp_{Z} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            tempxsec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix=="BjorkenX")
          {
            tmp = "{X}";
            y_title="d#sigma/d{X} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
          }
          else if (prefix=="Erecoil")
          {
            tmp = "E_{recoil}";
            y_title="d#sigma/d{X} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
          } 
          else if (prefix=="Emu")
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

          std::string title = tmp + " " + tgt;
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;
          plotter.DrawDataMCRatio(crossSection, simulatedCrossSection, 1, "TR", "Data");
          plotter.AddHistoTitle(title.c_str(), 0.039);
          can.Update();
          can.Print((tgt+prefix+"_CrossSectionRatio.pdf").c_str());
          can.Print((tgt+prefix+"_CrossSectionRatio.png").c_str());
        } */
        std::cout<<"Here1\n";
        { //Plotting migration
          TCanvas can("Temp");
          std::cout<<"Here2\n";
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kDefaultStyle);
          //plotter.axis_maximum = 0.4;

          std::cout<<"Here3\n";
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

          std::cout<<"Here4\n";
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
          else if (prefix.find("BjorkenX") != std::string::npos)
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

          std::cout<<"Here5\n";
          //Int_t colors[] = {0, 1};
          //gStyle->SetPalette(2, colors);
          plotter.axis_title_offset_z = 1.25;

          std::cout<<"Here6\n";
          plotter.DrawNormalizedMigrationHistogram(migration, false, false, true, false, 0.02);

          std::cout<<"Here7\n";
          plotter.WritePreliminary( 0.65, 0.2, 0.035, true);
          std::string title = tgtString + " " + tmp;
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 25, 1, 43, 22, 0);
          plotter.AddPlotLabel(yaxislabel.c_str() , 0.08, 0.9, 0.03, 1, 62, 33, 90);
          plotter.AddPlotLabel(xaxislabel.c_str() , 0.9, 0.08, 0.03, 1, 62, 33, 0);

          std::cout<<"Here8\n";
          //plotter.AddHistoTitle(title.c_str(), 0.039);
          can.Update();
          plotter.MultiPrint(&can, tgtname+"_"+prefix+"_migration", "png,pdf,C");

          std::cout<<"Here9\n";
          can.SetLogx(0);
          can.SetLogy(0);

          std::cout<<"Here10\n";
        }
        { //Plotting sidebands
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kDefaultStyle);
          plotter.draw_normalized_to_bin_width=false; //Already width normlized by normalise()
          plotter.axis_label_font=43;
          plotter.axis_label_size = 25;
          plotter.headroom = 1.2;
          plotter.mc_line_width=0;  //Why??
          plotter.ApplyAxisStyle();
          // Adjust pad heights manually (top bigger, bottom smaller)
          TPad *pad1 = new TPad("pad1","pad1",0,0.35,1,1.0); // top pad
          TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.35); // bottom pad
          pad1->SetBottomMargin(0); // remove x-axis label space for top
          pad2->SetTopMargin(0);
          pad2->SetBottomMargin(0.25);
          pad1->Draw();
          pad2->Draw();
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

          PlotUtils::MnvH1D* tmpBase = tempUSSignal->Clone();
          tmpBase->Add(tempUSOther);

          TObjArray* arrScaled = new TObjArray();
          PlotUtils::MnvH1D* US_scaleFactor = mapSFtoVarHist(*tempUSUS, *US_scaleFactor_tmp);
          PlotUtils::MnvH1D* tempUSScaled = tempUSUS->Clone();
          tempUSScaled->Multiply(tempUSScaled, US_scaleFactor);
          arrScaled->Add(tempUSSignal);
          arrScaled->Add(tempUSScaled);
          arrScaled->Add(tempUSDS);
          arrScaled->Add(tempUSOther);

          PlotUtils::MnvH1D* tmpUSSummed = tmpBase->Clone();
          tmpUSSummed->Add(tempUSUS);
          tmpUSSummed->Add(tempUSDS);
          PlotUtils::MnvH1D* tmpUSSummedScaled = tmpBase->Clone();
          tmpUSSummedScaled->Add(tempUSScaled);
          tmpUSSummedScaled->Add(tempUSDS);

          auto tmpUSData = DataUSSideband;
          tmpUSData->Scale(1./10000.);
          tmpUSData->GetXaxis()->SetRange(1, 81);
          tmpUSSummed->GetXaxis()->SetRange(1, 81);
          tmpUSSummedScaled->GetXaxis()->SetRange(1, 81);
          std::string title =  tgtname + " Upstream Plastic Sideband - Untuned";
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_size_x=0.03;
          plotter.axis_title_offset_x=1.8;
          plotter.axis_title_offset_y=1.8;
          //int channelColArr[5] = {TColor::GetColor("#d2c271"), TColor::GetColor("#CA5454"), TColor::GetColor("#40B0A6"), TColor::GetColor("#9285d5"), TColor::GetColor("#c2c2c2")};
          int channelColArr[4] = {TColor::GetColor("#d2c271"), TColor::GetColor("#CA5454"), TColor::GetColor("#40B0A6"), TColor::GetColor("#9285d5")};

          pad1->cd();

          plotter.DrawDataStackedMC(tmpUSData, arr, channelColArr, (dataPOT/(mcPOT*10000)), "TR", "Data", 1001, "MINER#nuA Segment", "Events per segment (#times 10^4)");
          plotter.AddHistoTitle(title.c_str(), 0.039);
          pad2->cd();
          plotter.DrawDataMCRatio(tmpUSData, tmpUSSummed, (dataPOT/(mcPOT*10000)), true, true, 0.5, 1.5, "");
          can.Update();
          plotter.MultiPrint(&can, tgtname+"_"+prefix+"_Sideband_US", "png,pdf,C");

          //std::string title =  tgt;
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_size_x=0.03;
          plotter.axis_title_offset_x=1.8;
          plotter.axis_title_offset_y=1.8;
          plotter.DrawDataStackedMC(tmpUSData, arrScaled, channelColArr, (dataPOT/(mcPOT*10000)), "TR", "Data", 1001, "MINER#nuA Segment", "Events per segment (#times 10^4)");
          title =  tgtname + " Upstream Plastic Sideband - Tuned";
          plotter.AddHistoTitle(title.c_str(), 0.039);
          pad2->cd();
          plotter.DrawDataMCRatio(tmpUSData, tmpUSSummedScaled, (dataPOT/(mcPOT*10000)), true, true, 0.5, 1.5, "");
          can.Update();
          plotter.MultiPrint(&can, tgtname+"_"+prefix+"_Sideband_US_Tuned", "png,pdf,C");

        }
        { //Plotting sidebands
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kDefaultStyle);
          //plotter.axis_maximum = 0.4;

          TObjArray* arr = new TObjArray();
          auto tempDSSignal = DSSidebandSignal->Clone();
          auto tempDSUS = DSSidebandUS->Clone();
          auto tempDSDS = DSSidebandDS->Clone();
          auto tempDSOther = DSSidebandOther->Clone();
          arr->Add(tempDSSignal);
          arr->Add(tempDSUS);
          arr->Add(tempDSDS);
          arr->Add(tempDSOther);

          TObjArray* arrScaled = new TObjArray();
          PlotUtils::MnvH1D* DS_scaleFactor = mapSFtoVarHist(*tempDSDS, *DS_scaleFactor_tmp);
          tempDSDS->Multiply(tempDSDS, DS_scaleFactor);
          arrScaled->Add(tempDSSignal);
          arrScaled->Add(tempDSUS);
          arrScaled->Add(tempDSDS);
          arrScaled->Add(tempDSOther);

          auto tmpDSData = DataDSSideband;
          tmpDSData->Scale(1./10000.);
          tmpDSData->GetXaxis()->SetRange(1, 81);
          std::string title =  tgtname + " Downstream Plastic Sideband - Untuned";
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_size_x=0.03;
          plotter.axis_title_offset_x=1.8;
          plotter.axis_title_offset_y=1.8;
          plotter.DrawDataStackedMC(tmpDSData, arr, (dataPOT/(mcPOT*10000)), "TR", "Data", 2, 1, 3001, "MINER#nuA Segment", "Events per segment (#times 10^4)");
          plotter.AddHistoTitle(title.c_str(), 0.039);
          can.Update();
          plotter.MultiPrint(&can, tgtname+"_"+prefix+"_Sideband_DS", "png,pdf,C");


          can.Clear();
          //std::string title =  tgt;
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_size_x=0.03;
          plotter.axis_title_offset_x=1.8;
          plotter.axis_title_offset_y=1.8;
          plotter.DrawDataStackedMC(tmpDSData, arrScaled, (dataPOT/(mcPOT*10000)), "TR", "Data", 2, 1, 3001, "MINER#nuA Segment", "Events per segment (#times 10^4)");
          title =  tgtname + " Downstream Plastic Sideband - Tuned";
          plotter.AddHistoTitle(title.c_str(), 0.039);
          can.Update();
          plotter.MultiPrint(&can, tgtname+"_"+prefix+"_Sideband_DS_Tuned", "png,pdf,C");
        }

        { //Plotting cross-sections
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


          TObjArray* arr = new TObjArray();
          auto temp2p2h = simulatedCrossSection2P2H->Clone();
          auto tempdis = simulatedCrossSectionDIS->Clone();
          auto tempres = simulatedCrossSectionRES->Clone();
          auto tempqe = simulatedCrossSectionQE->Clone();
          auto tempother = simulatedCrossSectionOther->Clone();
          auto tempsimxsec = simulatedCrossSection->Clone();
          auto tempxsec = crossSectionTuned->Clone();
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
          plotter.legend_offset_x = 0.0825;
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

          std::string title = nutypestr+ " " +tmp + " " + tgtString;
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;

          plotter.DrawErrorSummary(tempsimxsec, "TR", true, true,  0.00001, false, "", true, "", true);
          plotter.MultiPrint(&can, tgtname+"_"+prefix+"_"+nutypestr+"_XSecErrors", "png,pdf,C");
          can.Clear();


          // Adjust pad heights manually (top bigger, bottom smaller)
          TPad *pad1 = new TPad("pad1","pad1",0,0.35,1,1.0); // top pad
          TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.35); // bottom pad
          pad1->SetBottomMargin(0); // remove x-axis label space for top
          pad2->SetTopMargin(0);
          pad2->SetBottomMargin(0.25);
          pad1->Draw();
          pad2->Draw();

          if (prefix=="BjorkenX")
          {
            pad1->SetLogx();
            pad2->SetLogx();
          }


          if (tgtname == "Water")
          {
            std::stringstream POTstrFilled, POTstrEmpty;
            POTstrFilled << std::fixed << std::scientific << std::setprecision(2) << waterFilledPOTData;
            POTstrEmpty << std::fixed << std::scientific << std::setprecision(2) << waterEmptyPOTData;
            std::cout << "waterFilledPOTData: "<< waterFilledPOTData <<std::endl;
            std::cout << "waterEmptyPOTData: "<< waterEmptyPOTData <<std::endl;
            std::string potstrFilled = std::string("Water Filled POT ")+POTstrFilled.str();
            std::string potstrEmpty = std::string("Water Empty POT: ")+POTstrEmpty.str();// + " (Scaled to match filled)";
            plotter.AddPlotLabel(potstrFilled.c_str() , 0.18, 0.93, 0.025, 1, 42, 13, 0);
            plotter.AddPlotLabel(potstrEmpty.c_str() , 0.18, 0.90, 0.025, 1, 42, 13, 0);
          }
          else
          {
            std::stringstream POTstr;
            POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOT;
            std::string potstr = std::string("POT Used: ")+POTstr.str();
            plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.93, 0.025, 1, 42, 13, 0);
          }
          plotter.WritePreliminary( 0.5, 0.15, 0.035, true);
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 25, 1, 43, 22, 0);
          plotter.AddPlotLabel(y_title.c_str() , 0.08, 0.9, 0.03, 1, 42, 33, 90);
          plotter.AddPlotLabel("Data/MC" , 0.08, 0.3, 0.03, 1, 42, 33, 90);
          plotter.AddPlotLabel(xaxislabel.c_str() , 0.9, 0.08, 0.03, 1, 42, 33, 0);

          //Interaction Type
          pad1->cd(); 
          plotter.axis_minimum=0.001;

          int channelColArr[5] = {TColor::GetColor("#d2c271"), TColor::GetColor("#CA5454"), TColor::GetColor("#40B0A6"), TColor::GetColor("#9285d5"), TColor::GetColor("#c2c2c2")};
          plotter.DrawDataStackedMC(tempxsec, arr, channelColArr, 1, "TR", "Data", 1001, "", "", false, 1, xmax);
          pad2->cd();
          plotter.axis_minimum=defaultmin;
          plotter.DrawDataMCRatio(tempxsec, tempsimxsec, 1, true, true, 0.5, 1.5, "");
          can.SetLogx();
          pad1->cd();
          tempxsec->Draw("AXIS SAME");
          can.cd();
          plotter.MultiPrint(&can, tgtname+"_"+prefix+"_"+nutypestr+"_CrossSection", "png,pdf,C");
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
          
          auto tempsimxsec = simulatedCrossSection->Clone();
          auto tempxsec = crossSectionTuned->Clone();
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

          std::string title = nutypestr+ " " +tmp + " " + tgtString;
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;

          plotter.DrawErrorSummary(tempsimxsec, "TR", true, true,  0.00001, false, "", true, "", true);
          plotter.MultiPrint(&can, tgtname+"_"+prefix+"_"+nutypestr+"_XSecErrors", "png,pdf,C");
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


          if (tgtname == "Water")
          {
            std::stringstream POTstrFilled, POTstrEmpty;
            POTstrFilled << std::fixed << std::scientific << std::setprecision(2) << waterFilledPOTData;
            POTstrEmpty << std::fixed << std::scientific << std::setprecision(2) << waterEmptyPOTData;
            std::cout << "waterFilledPOTData: "<< waterFilledPOTData <<std::endl;
            std::cout << "waterEmptyPOTData: "<< waterEmptyPOTData <<std::endl;
            std::string potstrFilled = std::string("Water Filled POT ")+POTstrFilled.str();
            std::string potstrEmpty = std::string("Water Empty POT: ")+POTstrEmpty.str();// + " (Scaled to match filled)";
            plotter.AddPlotLabel(potstrFilled.c_str() , 0.18, 0.93, 0.025, 1, 42, 13, 0);
            plotter.AddPlotLabel(potstrEmpty.c_str() , 0.18, 0.90, 0.025, 1, 42, 13, 0);
          }
          else
          {
            std::stringstream POTstr;
            POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOT;
            std::string potstr = std::string("POT Used: ")+POTstr.str();
            plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.93, 0.025, 1, 42, 13, 0);
          }
          plotter.WritePreliminary( 0.5, 0.15, 0.035, true);
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 25, 1, 43, 22, 0);
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

          AnezkasResults->SetBinContent(15,0);
          AnezkasResults->SetBinError(15,0);
          AnezkasResultsWithErr->SetBinContent(15,0);
          AnezkasResultsWithErr->SetBinError(15,0);
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
          if (tgtname == "Iron")
          {
            AnezkasResults->SetBinContent(1,0.8049324);
            AnezkasResults->SetBinContent(2,3.07727);
            AnezkasResults->SetBinContent(3,7.291774);
            AnezkasResults->SetBinContent(4,12.37992);
            AnezkasResults->SetBinContent(5,16.24553);
            AnezkasResults->SetBinContent(6,19.54535);
            AnezkasResults->SetBinContent(7,21.76995);
            AnezkasResults->SetBinContent(8,21.48726);
            AnezkasResults->SetBinContent(9,17.36142);
            AnezkasResults->SetBinContent(10,11.5971);
            AnezkasResults->SetBinContent(11,5.383072);
            AnezkasResults->SetBinContent(12,1.516578);
            AnezkasResults->SetBinContent(13,0.2127459);
            AnezkasResults->SetBinContent(14,0.007011034);
            AnezkasResults->SetBinError(1,0.03904017);
            AnezkasResults->SetBinError(2,0.07088147);
            AnezkasResults->SetBinError(3,0.09911176);
            AnezkasResults->SetBinError(4,0.1442425);
            AnezkasResults->SetBinError(5,0.1686212);
            AnezkasResults->SetBinError(6,0.18823);
            AnezkasResults->SetBinError(7,0.207742);
            AnezkasResults->SetBinError(8,0.1723641);
            AnezkasResults->SetBinError(9,0.1652355);
            AnezkasResults->SetBinError(10,0.1409085);
            AnezkasResults->SetBinError(11,0.08428957);
            AnezkasResults->SetBinError(12,0.04414538);
            AnezkasResults->SetBinError(13,0.009251238);
            AnezkasResults->SetBinError(14,0.001713926);

            AnezkasResults->SetEntries(83535.27);
            AnezkasResults->SetStats(0);
            AnezkasResults->SetMarkerStyle(20);
            AnezkasResults->SetMarkerStyle(20);
            AnezkasResults->SetLineColor(kGreen);
            AnezkasResults->SetMarkerColor(kGreen);
            AnezkasResults->Draw("SAME HIST P E1 X0");

            AnezkasResultsWithErr->SetBinContent(1,0.8049324);
            AnezkasResultsWithErr->SetBinContent(2,3.07727);
            AnezkasResultsWithErr->SetBinContent(3,7.291774);
            AnezkasResultsWithErr->SetBinContent(4,12.37992);
            AnezkasResultsWithErr->SetBinContent(5,16.24553);
            AnezkasResultsWithErr->SetBinContent(6,19.54535);
            AnezkasResultsWithErr->SetBinContent(7,21.76995);
            AnezkasResultsWithErr->SetBinContent(8,21.48726);
            AnezkasResultsWithErr->SetBinContent(9,17.36142);
            AnezkasResultsWithErr->SetBinContent(10,11.5971);
            AnezkasResultsWithErr->SetBinContent(11,5.383072);
            AnezkasResultsWithErr->SetBinContent(12,1.516578);
            AnezkasResultsWithErr->SetBinContent(13,0.2127459);
            AnezkasResultsWithErr->SetBinContent(14,0.007011034);
            AnezkasResultsWithErr->SetBinError(1,0.08640005);
            AnezkasResultsWithErr->SetBinError(2,0.2700852);
            AnezkasResultsWithErr->SetBinError(3,0.5637781);
            AnezkasResultsWithErr->SetBinError(4,0.8657621);
            AnezkasResultsWithErr->SetBinError(5,1.055752);
            AnezkasResultsWithErr->SetBinError(6,1.191148);
            AnezkasResultsWithErr->SetBinError(7,1.270389);
            AnezkasResultsWithErr->SetBinError(8,1.265996);
            AnezkasResultsWithErr->SetBinError(9,1.139308);
            AnezkasResultsWithErr->SetBinError(10,0.9240369);
            AnezkasResultsWithErr->SetBinError(11,0.4985791);
            AnezkasResultsWithErr->SetBinError(12,0.1648089);
            AnezkasResultsWithErr->SetBinError(13,0.02620666);
            AnezkasResultsWithErr->SetBinError(14,0.001713926);

            AnezkasResultsWithErr->SetEntries(83535.27);
            AnezkasResultsWithErr->SetStats(0);
            AnezkasResultsWithErr->SetMarkerStyle(20);
            AnezkasResultsWithErr->SetMarkerStyle(20);
            AnezkasResultsWithErr->SetLineColor(kGreen);
            AnezkasResultsWithErr->SetMarkerColor(kGreen);
            AnezkasResultsWithErr->Draw("HIST E1 P SAME X0");

            AnezkasMCOther->SetBinContent(1,0.2332567);
            AnezkasMCOther->SetBinContent(2,0.5624879);
            AnezkasMCOther->SetBinContent(3,0.633126);
            AnezkasMCOther->SetBinContent(4,0.496955);
            AnezkasMCOther->SetBinContent(5,0.3215828);
            AnezkasMCOther->SetBinContent(6,0.2034721);
            AnezkasMCOther->SetBinContent(7,0.1084601);
            AnezkasMCOther->SetBinContent(8,0.0344629);
            AnezkasMCOther->SetBinContent(9,0.006515664);
            AnezkasMCOther->SetBinContent(10,0.0008418757);
            AnezkasMCOther->SetBinContent(11,0.0001914333);
            AnezkasMCOther->SetBinError(1,0.0068298);
            AnezkasMCOther->SetBinError(2,0.01010822);
            AnezkasMCOther->SetBinError(3,0.008868003);
            AnezkasMCOther->SetBinError(4,0.008983267);
            AnezkasMCOther->SetBinError(5,0.007298487);
            AnezkasMCOther->SetBinError(6,0.005908919);
            AnezkasMCOther->SetBinError(7,0.004356724);
            AnezkasMCOther->SetBinError(8,0.001781609);
            AnezkasMCOther->SetBinError(9,0.0008155981);
            AnezkasMCOther->SetBinError(10,0.0003013329);
            AnezkasMCOther->SetBinError(11,0.0001107458);
            AnezkasMCOther->SetEntries(16140.05);
            
            AnezkasMC2P2H->SetBinContent(1,0.2145634);
            AnezkasMC2P2H->SetBinContent(2,0.7113748);
            AnezkasMC2P2H->SetBinContent(3,1.433497);
            AnezkasMC2P2H->SetBinContent(4,2.353725);
            AnezkasMC2P2H->SetBinContent(5,3.159733);
            AnezkasMC2P2H->SetBinContent(6,3.52234);
            AnezkasMC2P2H->SetBinContent(7,2.870219);
            AnezkasMC2P2H->SetBinContent(8,1.630484);
            AnezkasMC2P2H->SetBinContent(9,0.6937129);
            AnezkasMC2P2H->SetBinContent(10,0.197009);
            AnezkasMC2P2H->SetBinContent(11,0.01291229);
            AnezkasMC2P2H->SetBinError(1,0.005516488);
            AnezkasMC2P2H->SetBinError(2,0.01027323);
            AnezkasMC2P2H->SetBinError(3,0.01481279);
            AnezkasMC2P2H->SetBinError(4,0.02967102);
            AnezkasMC2P2H->SetBinError(5,0.04384085);
            AnezkasMC2P2H->SetBinError(6,0.05003351);
            AnezkasMC2P2H->SetBinError(7,0.0391064);
            AnezkasMC2P2H->SetBinError(8,0.01417049);
            AnezkasMC2P2H->SetBinError(9,0.007013355);
            AnezkasMC2P2H->SetBinError(10,0.003702207);
            AnezkasMC2P2H->SetBinError(11,0.0007369864);
            AnezkasMC2P2H->SetEntries(37858.84);

            AnezkasMCQE->SetBinContent(1,0.08927685);
            AnezkasMCQE->SetBinContent(2,0.4338516);
            AnezkasMCQE->SetBinContent(3,1.176086);
            AnezkasMCQE->SetBinContent(4,2.10123);
            AnezkasMCQE->SetBinContent(5,2.797872);
            AnezkasMCQE->SetBinContent(6,3.326905);
            AnezkasMCQE->SetBinContent(7,3.617254);
            AnezkasMCQE->SetBinContent(8,3.526771);
            AnezkasMCQE->SetBinContent(9,2.679006);
            AnezkasMCQE->SetBinContent(10,1.674482);
            AnezkasMCQE->SetBinContent(11,0.6908877);
            AnezkasMCQE->SetBinContent(12,0.1624833);
            AnezkasMCQE->SetBinContent(13,0.008358044);
            AnezkasMCQE->SetBinContent(14,4.193384e-05);
            AnezkasMCQE->SetBinError(1,0.002824692);
            AnezkasMCQE->SetBinError(2,0.006069272);
            AnezkasMCQE->SetBinError(3,0.008855276);
            AnezkasMCQE->SetBinError(4,0.01429327);
            AnezkasMCQE->SetBinError(5,0.01685143);
            AnezkasMCQE->SetBinError(6,0.0189547);
            AnezkasMCQE->SetBinError(7,0.02058241);
            AnezkasMCQE->SetBinError(8,0.01534648);
            AnezkasMCQE->SetBinError(9,0.01414756);
            AnezkasMCQE->SetBinError(10,0.01141024);
            AnezkasMCQE->SetBinError(11,0.005619911);
            AnezkasMCQE->SetBinError(12,0.002678962);
            AnezkasMCQE->SetBinError(13,0.0003065252);
            AnezkasMCQE->SetBinError(14,2.144101e-05);
            AnezkasMCQE->SetEntries(248406.6);
            
            AnezkasMCDIS->SetBinContent(4,0.0009924486);
            AnezkasMCDIS->SetBinContent(5,0.006213854);
            AnezkasMCDIS->SetBinContent(6,0.02194196);
            AnezkasMCDIS->SetBinContent(7,0.07976177);
            AnezkasMCDIS->SetBinContent(8,0.5111634);
            AnezkasMCDIS->SetBinContent(9,1.89432);
            AnezkasMCDIS->SetBinContent(10,2.162233);
            AnezkasMCDIS->SetBinContent(11,1.146156);
            AnezkasMCDIS->SetBinContent(12,0.4910884);
            AnezkasMCDIS->SetBinContent(13,0.1109456);
            AnezkasMCDIS->SetBinContent(14,0.004940402);
            AnezkasMCDIS->SetBinError(4,0.0004231666);
            AnezkasMCDIS->SetBinError(5,0.001013966);
            AnezkasMCDIS->SetBinError(6,0.001877995);
            AnezkasMCDIS->SetBinError(7,0.003499099);
            AnezkasMCDIS->SetBinError(8,0.006022901);
            AnezkasMCDIS->SetBinError(9,0.01148782);
            AnezkasMCDIS->SetBinError(10,0.01242399);
            AnezkasMCDIS->SetBinError(11,0.00716388);
            AnezkasMCDIS->SetBinError(12,0.004825225);
            AnezkasMCDIS->SetBinError(13,0.001156979);
            AnezkasMCDIS->SetBinError(14,0.0002393437);
            AnezkasMCDIS->SetEntries(99340.93);

            AnezkasMCSoftDIS->SetBinContent(1,0.03248427);
            AnezkasMCSoftDIS->SetBinContent(2,0.2509049);
            AnezkasMCSoftDIS->SetBinContent(3,1.041699);
            AnezkasMCSoftDIS->SetBinContent(4,2.469559);
            AnezkasMCSoftDIS->SetBinContent(5,4.150664);
            AnezkasMCSoftDIS->SetBinContent(6,5.615304);
            AnezkasMCSoftDIS->SetBinContent(7,6.904437);
            AnezkasMCSoftDIS->SetBinContent(8,7.197054);
            AnezkasMCSoftDIS->SetBinContent(9,4.797907);
            AnezkasMCSoftDIS->SetBinContent(10,2.538382);
            AnezkasMCSoftDIS->SetBinContent(11,1.188943);
            AnezkasMCSoftDIS->SetBinContent(12,0.3268911);
            AnezkasMCSoftDIS->SetBinContent(13,0.02829924);
            AnezkasMCSoftDIS->SetBinContent(14,0.0006698488);
            AnezkasMCSoftDIS->SetBinError(1,0.002041519);
            AnezkasMCSoftDIS->SetBinError(2,0.00568662);
            AnezkasMCSoftDIS->SetBinError(3,0.009966142);
            AnezkasMCSoftDIS->SetBinError(4,0.01772028);
            AnezkasMCSoftDIS->SetBinError(5,0.02302679);
            AnezkasMCSoftDIS->SetBinError(6,0.02687907);
            AnezkasMCSoftDIS->SetBinError(7,0.02984525);
            AnezkasMCSoftDIS->SetBinError(8,0.02155445);
            AnezkasMCSoftDIS->SetBinError(9,0.01754847);
            AnezkasMCSoftDIS->SetBinError(10,0.01265984);
            AnezkasMCSoftDIS->SetBinError(11,0.006730957);
            AnezkasMCSoftDIS->SetBinError(12,0.003602525);
            AnezkasMCSoftDIS->SetBinError(13,0.0005555895);
            AnezkasMCSoftDIS->SetBinError(14,8.816347e-05);
            AnezkasMCSoftDIS->SetEntries(372515.3);

            AnezkasMCRES->SetBinContent(1,0.6699978);
            AnezkasMCRES->SetBinContent(2,2.016599);
            AnezkasMCRES->SetBinContent(3,3.980916);
            AnezkasMCRES->SetBinContent(4,5.731707);
            AnezkasMCRES->SetBinContent(5,6.868337);
            AnezkasMCRES->SetBinContent(6,7.393548);
            AnezkasMCRES->SetBinContent(7,7.440685);
            AnezkasMCRES->SetBinContent(8,6.488466);
            AnezkasMCRES->SetBinContent(9,4.536648);
            AnezkasMCRES->SetBinContent(10,2.703977);
            AnezkasMCRES->SetBinContent(11,1.170722);
            AnezkasMCRES->SetBinContent(12,0.28281);
            AnezkasMCRES->SetBinContent(13,0.01850774);
            AnezkasMCRES->SetBinContent(14,0.0001079154);
            AnezkasMCRES->SetBinError(1,0.009618373);
            AnezkasMCRES->SetBinError(2,0.01657863);
            AnezkasMCRES->SetBinError(3,0.02055881);
            AnezkasMCRES->SetBinError(4,0.02907097);
            AnezkasMCRES->SetBinError(5,0.0321493);
            AnezkasMCRES->SetBinError(6,0.03357979);
            AnezkasMCRES->SetBinError(7,0.03381782);
            AnezkasMCRES->SetBinError(8,0.02236415);
            AnezkasMCRES->SetBinError(9,0.01876228);
            AnezkasMCRES->SetBinError(10,0.01456123);
            AnezkasMCRES->SetBinError(11,0.007539753);
            AnezkasMCRES->SetBinError(12,0.00380736);
            AnezkasMCRES->SetBinError(13,0.0005096632);
            AnezkasMCRES->SetBinError(14,3.967509e-05);
            AnezkasMCRES->SetEntries(400076.3);



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

          }
          else if (tgtname == "Lead")
          {
            AnezkasResults->SetBinContent(1,0.6561926);
            AnezkasResults->SetBinContent(2,2.441443);
            AnezkasResults->SetBinContent(3,6.129783);
            AnezkasResults->SetBinContent(4,10.60599);
            AnezkasResults->SetBinContent(5,14.75021);
            AnezkasResults->SetBinContent(6,18.55033);
            AnezkasResults->SetBinContent(7,20.93295);
            AnezkasResults->SetBinContent(8,20.4205);
            AnezkasResults->SetBinContent(9,16.04993);
            AnezkasResults->SetBinContent(10,10.51605);
            AnezkasResults->SetBinContent(11,4.826657);
            AnezkasResults->SetBinContent(12,1.519215);
            AnezkasResults->SetBinContent(13,0.2152156);
            AnezkasResults->SetBinContent(14,0.009208849);
            AnezkasResults->SetBinError(1,0.0338608);
            AnezkasResults->SetBinError(2,0.06054109);
            AnezkasResults->SetBinError(3,0.08586898);
            AnezkasResults->SetBinError(4,0.1243123);
            AnezkasResults->SetBinError(5,0.1469578);
            AnezkasResults->SetBinError(6,0.1658356);
            AnezkasResults->SetBinError(7,0.18216);
            AnezkasResults->SetBinError(8,0.1516043);
            AnezkasResults->SetBinError(9,0.1440136);
            AnezkasResults->SetBinError(10,0.1211898);
            AnezkasResults->SetBinError(11,0.0732774);
            AnezkasResults->SetBinError(12,0.04079451);
            AnezkasResults->SetBinError(13,0.008292306);
            AnezkasResults->SetBinError(14,0.001820963);
            AnezkasResults->SetEntries(92830.36);
            AnezkasResults->SetStats(0);
            AnezkasResults->SetMarkerStyle(20);
            AnezkasResults->SetMarkerStyle(20);
            AnezkasResults->SetLineColor(kGreen);
            AnezkasResults->SetMarkerColor(kGreen);
            AnezkasResults->Draw("SAME HIST P E1 X0");

            AnezkasResultsWithErr->SetBinContent(1,0.6561926);
            AnezkasResultsWithErr->SetBinContent(2,2.441443);
            AnezkasResultsWithErr->SetBinContent(3,6.129783);
            AnezkasResultsWithErr->SetBinContent(4,10.60599);
            AnezkasResultsWithErr->SetBinContent(5,14.75021);
            AnezkasResultsWithErr->SetBinContent(6,18.55033);
            AnezkasResultsWithErr->SetBinContent(7,20.93295);
            AnezkasResultsWithErr->SetBinContent(8,20.4205);
            AnezkasResultsWithErr->SetBinContent(9,16.04993);
            AnezkasResultsWithErr->SetBinContent(10,10.51605);
            AnezkasResultsWithErr->SetBinContent(11,4.826657);
            AnezkasResultsWithErr->SetBinContent(12,1.519215);
            AnezkasResultsWithErr->SetBinContent(13,0.2152156);
            AnezkasResultsWithErr->SetBinContent(14,0.009208849);
            AnezkasResultsWithErr->SetBinError(1,0.07677816);
            AnezkasResultsWithErr->SetBinError(2,0.2313552);
            AnezkasResultsWithErr->SetBinError(3,0.4977336);
            AnezkasResultsWithErr->SetBinError(4,0.7779536);
            AnezkasResultsWithErr->SetBinError(5,0.9743027);
            AnezkasResultsWithErr->SetBinError(6,1.142635);
            AnezkasResultsWithErr->SetBinError(7,1.230503);
            AnezkasResultsWithErr->SetBinError(8,1.212708);
            AnezkasResultsWithErr->SetBinError(9,1.090595);
            AnezkasResultsWithErr->SetBinError(10,0.8416227);
            AnezkasResultsWithErr->SetBinError(11,0.4436641);
            AnezkasResultsWithErr->SetBinError(12,0.1553531);
            AnezkasResultsWithErr->SetBinError(13,0.02462683);
            AnezkasResultsWithErr->SetBinError(14,0.001820963);
            AnezkasResultsWithErr->SetEntries(92830.36);
            AnezkasResultsWithErr->SetStats(0);
            AnezkasResultsWithErr->SetMarkerStyle(20);
            AnezkasResultsWithErr->SetMarkerStyle(20);
            AnezkasResultsWithErr->SetLineColor(kGreen);
            AnezkasResultsWithErr->SetMarkerColor(kGreen);
            AnezkasResultsWithErr->Draw("HIST E1 p SAME X0");


            AnezkasMCOther->SetBinContent(1,0.1561364);
            AnezkasMCOther->SetBinContent(2,0.3256761);
            AnezkasMCOther->SetBinContent(3,0.2938892);
            AnezkasMCOther->SetBinContent(4,0.1913129);
            AnezkasMCOther->SetBinContent(5,0.1094337);
            AnezkasMCOther->SetBinContent(6,0.04710294);
            AnezkasMCOther->SetBinContent(7,0.01823111);
            AnezkasMCOther->SetBinContent(8,0.004850307);
            AnezkasMCOther->SetBinContent(9,0.0004005317);
            AnezkasMCOther->SetBinError(1,0.005368682);
            AnezkasMCOther->SetBinError(2,0.007193071);
            AnezkasMCOther->SetBinError(3,0.005704332);
            AnezkasMCOther->SetBinError(4,0.005383224);
            AnezkasMCOther->SetBinError(5,0.004109163);
            AnezkasMCOther->SetBinError(6,0.002756873);
            AnezkasMCOther->SetBinError(7,0.001745703);
            AnezkasMCOther->SetBinError(8,0.0006742341);
            AnezkasMCOther->SetBinError(9,0.0001811571);
            AnezkasMCOther->SetEntries(7734.685);

            AnezkasMC2P2H->SetBinContent(1,0.2037768);
            AnezkasMC2P2H->SetBinContent(2,0.6652423);
            AnezkasMC2P2H->SetBinContent(3,1.343499);
            AnezkasMC2P2H->SetBinContent(4,2.190119);
            AnezkasMC2P2H->SetBinContent(5,2.940682);
            AnezkasMC2P2H->SetBinContent(6,3.093817);
            AnezkasMC2P2H->SetBinContent(7,2.655823);
            AnezkasMC2P2H->SetBinContent(8,1.592262);
            AnezkasMC2P2H->SetBinContent(9,0.6774053);
            AnezkasMC2P2H->SetBinContent(10,0.1839579);
            AnezkasMC2P2H->SetBinContent(11,0.009885401);
            AnezkasMC2P2H->SetBinError(1,0.005199156);
            AnezkasMC2P2H->SetBinError(2,0.009613387);
            AnezkasMC2P2H->SetBinError(3,0.01371744);
            AnezkasMC2P2H->SetBinError(4,0.02776389);
            AnezkasMC2P2H->SetBinError(5,0.04100238);
            AnezkasMC2P2H->SetBinError(6,0.04412376);
            AnezkasMC2P2H->SetBinError(7,0.03627244);
            AnezkasMC2P2H->SetBinError(8,0.01353946);
            AnezkasMC2P2H->SetBinError(9,0.006668482);
            AnezkasMC2P2H->SetBinError(10,0.003447541);
            AnezkasMC2P2H->SetBinError(11,0.0006229185);
            AnezkasMC2P2H->SetEntries(38644.54);

            AnezkasMCQE->SetBinContent(1,0.02843428);
            AnezkasMCQE->SetBinContent(2,0.2436697);
            AnezkasMCQE->SetBinContent(3,0.8477265);
            AnezkasMCQE->SetBinContent(4,1.669573);
            AnezkasMCQE->SetBinContent(5,2.356901);
            AnezkasMCQE->SetBinContent(6,2.926238);
            AnezkasMCQE->SetBinContent(7,3.249881);
            AnezkasMCQE->SetBinContent(8,3.150217);
            AnezkasMCQE->SetBinContent(9,2.437865);
            AnezkasMCQE->SetBinContent(10,1.498064);
            AnezkasMCQE->SetBinContent(11,0.6306986);
            AnezkasMCQE->SetBinContent(12,0.1507706);
            AnezkasMCQE->SetBinContent(13,0.007408994);
            AnezkasMCQE->SetBinContent(14,1.200343e-05);
            AnezkasMCQE->SetBinError(1,0.00151774);
            AnezkasMCQE->SetBinError(2,0.004409623);
            AnezkasMCQE->SetBinError(3,0.007107976);
            AnezkasMCQE->SetBinError(4,0.01186107);
            AnezkasMCQE->SetBinError(5,0.01453164);
            AnezkasMCQE->SetBinError(6,0.01683605);
            AnezkasMCQE->SetBinError(7,0.01876005);
            AnezkasMCQE->SetBinError(8,0.01403453);
            AnezkasMCQE->SetBinError(9,0.01309157);
            AnezkasMCQE->SetBinError(10,0.01043286);
            AnezkasMCQE->SetBinError(11,0.005182219);
            AnezkasMCQE->SetBinError(12,0.002486239);
            AnezkasMCQE->SetBinError(13,0.0002760264);
            AnezkasMCQE->SetBinError(14,9.724431e-06);
            AnezkasMCQE->SetEntries(234765.3);

            AnezkasMCDIS->SetBinContent(4,0.0005425822);
            AnezkasMCDIS->SetBinContent(5,0.004619062);
            AnezkasMCDIS->SetBinContent(6,0.01810059);
            AnezkasMCDIS->SetBinContent(7,0.07383283);
            AnezkasMCDIS->SetBinContent(8,0.4788902);
            AnezkasMCDIS->SetBinContent(9,1.767908);
            AnezkasMCDIS->SetBinContent(10,2.000405);
            AnezkasMCDIS->SetBinContent(11,1.045321);
            AnezkasMCDIS->SetBinContent(12,0.4329718);
            AnezkasMCDIS->SetBinContent(13,0.09741804);
            AnezkasMCDIS->SetBinContent(14,0.004538899);
            AnezkasMCDIS->SetBinError(4,0.0002563791);
            AnezkasMCDIS->SetBinError(5,0.0008462986);
            AnezkasMCDIS->SetBinError(6,0.001614761);
            AnezkasMCDIS->SetBinError(7,0.003222523);
            AnezkasMCDIS->SetBinError(8,0.005577529);
            AnezkasMCDIS->SetBinError(9,0.01058429);
            AnezkasMCDIS->SetBinError(10,0.01142514);
            AnezkasMCDIS->SetBinError(11,0.006521924);
            AnezkasMCDIS->SetBinError(12,0.004331015);
            AnezkasMCDIS->SetBinError(13,0.001040553);
            AnezkasMCDIS->SetBinError(14,0.0002194781);
            AnezkasMCDIS->SetEntries(100184.5);

            AnezkasMCSoftDIS->SetBinContent(1,0.03372605);
            AnezkasMCSoftDIS->SetBinContent(2,0.2495663);
            AnezkasMCSoftDIS->SetBinContent(3,1.049851);
            AnezkasMCSoftDIS->SetBinContent(4,2.461786);
            AnezkasMCSoftDIS->SetBinContent(5,4.062952);
            AnezkasMCSoftDIS->SetBinContent(6,5.550518);
            AnezkasMCSoftDIS->SetBinContent(7,6.79058);
            AnezkasMCSoftDIS->SetBinContent(8,7.010942);
            AnezkasMCSoftDIS->SetBinContent(9,4.685369);
            AnezkasMCSoftDIS->SetBinContent(10,2.463943);
            AnezkasMCSoftDIS->SetBinContent(11,1.159185);
            AnezkasMCSoftDIS->SetBinContent(12,0.3209256);
            AnezkasMCSoftDIS->SetBinContent(13,0.02734468);
            AnezkasMCSoftDIS->SetBinContent(14,0.0005226055);
            AnezkasMCSoftDIS->SetBinError(1,0.002026829);
            AnezkasMCSoftDIS->SetBinError(2,0.005462103);
            AnezkasMCSoftDIS->SetBinError(3,0.009608662);
            AnezkasMCSoftDIS->SetBinError(4,0.01706577);
            AnezkasMCSoftDIS->SetBinError(5,0.02192089);
            AnezkasMCSoftDIS->SetBinError(6,0.02562423);
            AnezkasMCSoftDIS->SetBinError(7,0.02849947);
            AnezkasMCSoftDIS->SetBinError(8,0.02047803);
            AnezkasMCSoftDIS->SetBinError(9,0.0167475);
            AnezkasMCSoftDIS->SetBinError(10,0.01202345);
            AnezkasMCSoftDIS->SetBinError(11,0.006425167);
            AnezkasMCSoftDIS->SetBinError(12,0.003445254);
            AnezkasMCSoftDIS->SetBinError(13,0.0005212152);
            AnezkasMCSoftDIS->SetBinError(14,7.395283e-05);
            AnezkasMCSoftDIS->SetEntries(394042.5);

            AnezkasMCRES->SetBinContent(1,0.7004375);
            AnezkasMCRES->SetBinContent(2,2.102252);
            AnezkasMCRES->SetBinContent(3,4.21529);
            AnezkasMCRES->SetBinContent(4,6.092341);
            AnezkasMCRES->SetBinContent(5,7.286666);
            AnezkasMCRES->SetBinContent(6,7.859736);
            AnezkasMCRES->SetBinContent(7,7.87675);
            AnezkasMCRES->SetBinContent(8,6.817545);
            AnezkasMCRES->SetBinContent(9,4.705162);
            AnezkasMCRES->SetBinContent(10,2.823745);
            AnezkasMCRES->SetBinContent(11,1.22242);
            AnezkasMCRES->SetBinContent(12,0.3025925);
            AnezkasMCRES->SetBinContent(13,0.01880532);
            AnezkasMCRES->SetBinContent(14,0.0001239769);
            AnezkasMCRES->SetBinError(1,0.009418198);
            AnezkasMCRES->SetBinError(2,0.01624629);
            AnezkasMCRES->SetBinError(3,0.02036675);
            AnezkasMCRES->SetBinError(4,0.02893163);
            AnezkasMCRES->SetBinError(5,0.03205692);
            AnezkasMCRES->SetBinError(6,0.03351968);
            AnezkasMCRES->SetBinError(7,0.03377675);
            AnezkasMCRES->SetBinError(8,0.02225132);
            AnezkasMCRES->SetBinError(9,0.01850397);
            AnezkasMCRES->SetBinError(10,0.01442289);
            AnezkasMCRES->SetBinError(11,0.007475384);
            AnezkasMCRES->SetBinError(12,0.003845522);
            AnezkasMCRES->SetBinError(13,0.0004990566);
            AnezkasMCRES->SetBinError(14,4.454069e-05);
            AnezkasMCRES->SetEntries(450100.9);

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
          }
          else if (tgtname == "Carbon")
          {
            AnezkasResults->SetBinContent(1,1.076004);
            AnezkasResults->SetBinContent(2,3.947597);
            AnezkasResults->SetBinContent(3,8.576456);
            AnezkasResults->SetBinContent(4,13.91362);
            AnezkasResults->SetBinContent(5,17.69146);
            AnezkasResults->SetBinContent(6,20.40097);
            AnezkasResults->SetBinContent(7,22.94916);
            AnezkasResults->SetBinContent(8,22.67642);
            AnezkasResults->SetBinContent(9,17.50646);
            AnezkasResults->SetBinContent(10,12.64096);
            AnezkasResults->SetBinContent(11,5.889888);
            AnezkasResults->SetBinContent(12,1.67176);
            AnezkasResults->SetBinContent(13,0.2231268);
            AnezkasResults->SetBinContent(14,0.009600314);
            AnezkasResults->SetBinError(1,0.07895694);
            AnezkasResults->SetBinError(2,0.1421401);
            AnezkasResults->SetBinError(3,0.190154);
            AnezkasResults->SetBinError(4,0.2708475);
            AnezkasResults->SetBinError(5,0.3068982);
            AnezkasResults->SetBinError(6,0.3334209);
            AnezkasResults->SetBinError(7,0.3643222);
            AnezkasResults->SetBinError(8,0.3089035);
            AnezkasResults->SetBinError(9,0.2933375);
            AnezkasResults->SetBinError(10,0.2667108);
            AnezkasResults->SetBinError(11,0.1633153);
            AnezkasResults->SetBinError(12,0.08528932);
            AnezkasResults->SetBinError(13,0.01706879);
            AnezkasResults->SetBinError(14,0.004045316);
            AnezkasResults->SetEntries(29242.21);
            AnezkasResults->SetStats(0);
            AnezkasResults->SetMarkerStyle(20);
            AnezkasResults->SetMarkerStyle(20);
            AnezkasResults->SetLineColor(kGreen);
            AnezkasResults->SetMarkerColor(kGreen);
            AnezkasResults->Draw("SAME HIST P E1 X0");

            AnezkasResultsWithErr->SetBinContent(1,1.076004);
            AnezkasResultsWithErr->SetBinContent(2,3.947597);
            AnezkasResultsWithErr->SetBinContent(3,8.576456);
            AnezkasResultsWithErr->SetBinContent(4,13.91362);
            AnezkasResultsWithErr->SetBinContent(5,17.69146);
            AnezkasResultsWithErr->SetBinContent(6,20.40097);
            AnezkasResultsWithErr->SetBinContent(7,22.94916);
            AnezkasResultsWithErr->SetBinContent(8,22.67642);
            AnezkasResultsWithErr->SetBinContent(9,17.50646);
            AnezkasResultsWithErr->SetBinContent(10,12.64096);
            AnezkasResultsWithErr->SetBinContent(11,5.889888);
            AnezkasResultsWithErr->SetBinContent(12,1.67176);
            AnezkasResultsWithErr->SetBinContent(13,0.2231268);
            AnezkasResultsWithErr->SetBinContent(14,0.009600314);
            AnezkasResultsWithErr->SetBinError(1,0.1232578);
            AnezkasResultsWithErr->SetBinError(2,0.342858);
            AnezkasResultsWithErr->SetBinError(3,0.6556305);
            AnezkasResultsWithErr->SetBinError(4,0.9964636);
            AnezkasResultsWithErr->SetBinError(5,1.16301);
            AnezkasResultsWithErr->SetBinError(6,1.256069);
            AnezkasResultsWithErr->SetBinError(7,1.370369);
            AnezkasResultsWithErr->SetBinError(8,1.380239);
            AnezkasResultsWithErr->SetBinError(9,1.249442);
            AnezkasResultsWithErr->SetBinError(10,1.010363);
            AnezkasResultsWithErr->SetBinError(11,0.5789273);
            AnezkasResultsWithErr->SetBinError(12,0.1972076);
            AnezkasResultsWithErr->SetBinError(13,0.03212316);
            AnezkasResultsWithErr->SetBinError(14,0.004045316);
            AnezkasResultsWithErr->SetEntries(29242.21);
            AnezkasResultsWithErr->SetStats(0);
            AnezkasResultsWithErr->SetStats(0);
            AnezkasResultsWithErr->SetMarkerStyle(20);
            AnezkasResultsWithErr->SetMarkerStyle(20);
            AnezkasResultsWithErr->SetLineColor(kGreen);
            AnezkasResultsWithErr->SetMarkerColor(kGreen);
            AnezkasResultsWithErr->Draw("HIST E1 p SAME X0");


            AnezkasMCOther->SetBinContent(1,0.2347042);
            AnezkasMCOther->SetBinContent(2,0.6614324);
            AnezkasMCOther->SetBinContent(3,0.8574958);
            AnezkasMCOther->SetBinContent(4,0.7850597);
            AnezkasMCOther->SetBinContent(5,0.6588785);
            AnezkasMCOther->SetBinContent(6,0.4457613);
            AnezkasMCOther->SetBinContent(7,0.352723);
            AnezkasMCOther->SetBinContent(8,0.1736719);
            AnezkasMCOther->SetBinContent(9,0.04092633);
            AnezkasMCOther->SetBinContent(10,0.01039743);
            AnezkasMCOther->SetBinContent(11,0.002213694);
            AnezkasMCOther->SetBinContent(12,0.0002075966);
            AnezkasMCOther->SetBinError(1,0.01275057);
            AnezkasMCOther->SetBinError(2,0.02126815);
            AnezkasMCOther->SetBinError(3,0.02001447);
            AnezkasMCOther->SetBinError(4,0.02160056);
            AnezkasMCOther->SetBinError(5,0.01964789);
            AnezkasMCOther->SetBinError(6,0.01617358);
            AnezkasMCOther->SetBinError(7,0.01441865);
            AnezkasMCOther->SetBinError(8,0.00727417);
            AnezkasMCOther->SetBinError(9,0.003656266);
            AnezkasMCOther->SetBinError(10,0.001866059);
            AnezkasMCOther->SetBinError(11,0.0007118107);
            AnezkasMCOther->SetBinError(12,0.0002075966);
            AnezkasMCOther->SetEntries(7407.965);

            AnezkasMC2P2H->SetBinContent(1,0.2290859);
            AnezkasMC2P2H->SetBinContent(2,0.6662635);
            AnezkasMC2P2H->SetBinContent(3,1.384212);
            AnezkasMC2P2H->SetBinContent(4,2.339865);
            AnezkasMC2P2H->SetBinContent(5,3.233154);
            AnezkasMC2P2H->SetBinContent(6,3.223004);
            AnezkasMC2P2H->SetBinContent(7,2.950755);
            AnezkasMC2P2H->SetBinContent(8,1.67929);
            AnezkasMC2P2H->SetBinContent(9,0.6975054);
            AnezkasMC2P2H->SetBinContent(10,0.2013269);
            AnezkasMC2P2H->SetBinContent(11,0.0159375);
            AnezkasMC2P2H->SetBinError(1,0.01100673);
            AnezkasMC2P2H->SetBinError(2,0.01926947);
            AnezkasMC2P2H->SetBinError(3,0.02857379);
            AnezkasMC2P2H->SetBinError(4,0.0580588);
            AnezkasMC2P2H->SetBinError(5,0.08719479);
            AnezkasMC2P2H->SetBinError(6,0.09103463);
            AnezkasMC2P2H->SetBinError(7,0.07771348);
            AnezkasMC2P2H->SetBinError(8,0.02856847);
            AnezkasMC2P2H->SetBinError(9,0.01354608);
            AnezkasMC2P2H->SetBinError(10,0.007150189);
            AnezkasMC2P2H->SetBinError(11,0.001564069);
            AnezkasMC2P2H->SetEntries(9985.966);

            AnezkasMCQE->SetBinContent(1,0.132169);
            AnezkasMCQE->SetBinContent(2,0.5534028);
            AnezkasMCQE->SetBinContent(3,1.406729);
            AnezkasMCQE->SetBinContent(4,2.380878);
            AnezkasMCQE->SetBinContent(5,3.130233);
            AnezkasMCQE->SetBinContent(6,3.632809);
            AnezkasMCQE->SetBinContent(7,3.731393);
            AnezkasMCQE->SetBinContent(8,3.564532);
            AnezkasMCQE->SetBinContent(9,2.766282);
            AnezkasMCQE->SetBinContent(10,1.706072);
            AnezkasMCQE->SetBinContent(11,0.7637538);
            AnezkasMCQE->SetBinContent(12,0.1856746);
            AnezkasMCQE->SetBinContent(13,0.01007721);
            AnezkasMCQE->SetBinContent(14,4.557366e-05);
            AnezkasMCQE->SetBinError(1,0.006464369);
            AnezkasMCQE->SetBinError(2,0.01302211);
            AnezkasMCQE->SetBinError(3,0.01853815);
            AnezkasMCQE->SetBinError(4,0.02918155);
            AnezkasMCQE->SetBinError(5,0.03448448);
            AnezkasMCQE->SetBinError(6,0.03851421);
            AnezkasMCQE->SetBinError(7,0.0404669);
            AnezkasMCQE->SetBinError(8,0.02956611);
            AnezkasMCQE->SetBinError(9,0.02745944);
            AnezkasMCQE->SetBinError(10,0.02206261);
            AnezkasMCQE->SetBinError(11,0.01139503);
            AnezkasMCQE->SetBinError(12,0.005576629);
            AnezkasMCQE->SetBinError(13,0.0006528596);
            AnezkasMCQE->SetBinError(14,4.557366e-05);
            AnezkasMCQE->SetEntries(71847.15);

          AnezkasMCDIS->SetBinContent(4,0.0009873169);
          AnezkasMCDIS->SetBinContent(5,0.006217955);
          AnezkasMCDIS->SetBinContent(6,0.01498795);
          AnezkasMCDIS->SetBinContent(7,0.07820934);
          AnezkasMCDIS->SetBinContent(8,0.5671901);
          AnezkasMCDIS->SetBinContent(9,2.031549);
          AnezkasMCDIS->SetBinContent(10,2.389959);
          AnezkasMCDIS->SetBinContent(11,1.223928);
          AnezkasMCDIS->SetBinContent(12,0.5134881);
          AnezkasMCDIS->SetBinContent(13,0.1157349);
          AnezkasMCDIS->SetBinContent(14,0.004344186);
          AnezkasMCDIS->SetBinError(4,0.0006981403);
          AnezkasMCDIS->SetBinError(5,0.002033867);
          AnezkasMCDIS->SetBinError(6,0.002957023);
          AnezkasMCDIS->SetBinError(7,0.006639578);
          AnezkasMCDIS->SetBinError(8,0.01235666);
          AnezkasMCDIS->SetBinError(9,0.02298054);
          AnezkasMCDIS->SetBinError(10,0.02520459);
          AnezkasMCDIS->SetBinError(11,0.01429955);
          AnezkasMCDIS->SetBinError(12,0.009537408);
          AnezkasMCDIS->SetBinError(13,0.002286753);
          AnezkasMCDIS->SetBinError(14,0.0004359495);
          AnezkasMCDIS->SetEntries(28787.06);

          AnezkasMCSoftDIS->SetBinContent(1,0.03620217);
          AnezkasMCSoftDIS->SetBinContent(2,0.2458733);
          AnezkasMCSoftDIS->SetBinContent(3,1.007919);
          AnezkasMCSoftDIS->SetBinContent(4,2.523109);
          AnezkasMCSoftDIS->SetBinContent(5,4.1756);
          AnezkasMCSoftDIS->SetBinContent(6,5.866563);
          AnezkasMCSoftDIS->SetBinContent(7,7.126496);
          AnezkasMCSoftDIS->SetBinContent(8,7.542857);
          AnezkasMCSoftDIS->SetBinContent(9,4.928868);
          AnezkasMCSoftDIS->SetBinContent(10,2.624999);
          AnezkasMCSoftDIS->SetBinContent(11,1.265241);
          AnezkasMCSoftDIS->SetBinContent(12,0.3495936);
          AnezkasMCSoftDIS->SetBinContent(13,0.02728756);
          AnezkasMCSoftDIS->SetBinContent(14,0.0004640762);
          AnezkasMCSoftDIS->SetBinError(1,0.004279062);
          AnezkasMCSoftDIS->SetBinError(2,0.01081395);
          AnezkasMCSoftDIS->SetBinError(3,0.01878313);
          AnezkasMCSoftDIS->SetBinError(4,0.03454548);
          AnezkasMCSoftDIS->SetBinError(5,0.04460561);
          AnezkasMCSoftDIS->SetBinError(6,0.05308649);
          AnezkasMCSoftDIS->SetBinError(7,0.05855016);
          AnezkasMCSoftDIS->SetBinError(8,0.04263792);
          AnezkasMCSoftDIS->SetBinError(9,0.0342632);
          AnezkasMCSoftDIS->SetBinError(10,0.0247152);
          AnezkasMCSoftDIS->SetBinError(11,0.01339069);
          AnezkasMCSoftDIS->SetBinError(12,0.007196595);
          AnezkasMCSoftDIS->SetBinError(13,0.001068296);
          AnezkasMCSoftDIS->SetBinError(14,0.0001407476);
          AnezkasMCSoftDIS->SetEntries(103461.2);

          AnezkasMCRES->SetBinContent(1,0.6330794);
          AnezkasMCRES->SetBinContent(2,1.92406);
          AnezkasMCRES->SetBinContent(3,3.878374);
          AnezkasMCRES->SetBinContent(4,5.545959);
          AnezkasMCRES->SetBinContent(5,6.749652);
          AnezkasMCRES->SetBinContent(6,7.306854);
          AnezkasMCRES->SetBinContent(7,7.275363);
          AnezkasMCRES->SetBinContent(8,6.427497);
          AnezkasMCRES->SetBinContent(9,4.658767);
          AnezkasMCRES->SetBinContent(10,2.8251);
          AnezkasMCRES->SetBinContent(11,1.220214);
          AnezkasMCRES->SetBinContent(12,0.3143608);
          AnezkasMCRES->SetBinContent(13,0.02013586);
          AnezkasMCRES->SetBinContent(14,0.0001953734);
          AnezkasMCRES->SetBinError(1,0.0178625);
          AnezkasMCRES->SetBinError(2,0.03100077);
          AnezkasMCRES->SetBinError(3,0.03894752);
          AnezkasMCRES->SetBinError(4,0.05501879);
          AnezkasMCRES->SetBinError(5,0.0614835);
          AnezkasMCRES->SetBinError(6,0.0645396);
          AnezkasMCRES->SetBinError(7,0.06455953);
          AnezkasMCRES->SetBinError(8,0.04309249);
          AnezkasMCRES->SetBinError(9,0.03679076);
          AnezkasMCRES->SetBinError(10,0.02884715);
          AnezkasMCRES->SetBinError(11,0.01486384);
          AnezkasMCRES->SetBinError(12,0.007841854);
          AnezkasMCRES->SetBinError(13,0.001033812);
          AnezkasMCRES->SetBinError(14,0.0001220556);
          AnezkasMCRES->SetEntries(106874.5);

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
          TH1D *tempTH12p2h = (TH1D*)temp2p2h;
          TH1D *tempTH1dis = (TH1D*)tempdis;
          TH1D *tempTH1res = (TH1D*)tempres;
          TH1D *tempTH1qe = (TH1D*)tempqe;
          TH1D *tempTH1other = (TH1D*)tempother;
          AnezkasMCDISSum->Add(AnezkasMCSoftDIS);
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
          //tempxsec->Draw("SAME HIST P E1 X0"); //Draw it back on top
          can.cd();
          plotter.MultiPrint(&can, tgtname+"_"+prefix+"_"+nutypestr+"_CrossSectionWithAnezkas", "png,pdf,C");
        }

        std::string daisyFilepath = "./"+nutypestr+"_"+prefix+"_Daisy_crossSection.root";
        std::string trackerNoDaisyFilepath = "./"+nutypestr+"_"+prefix+"_crossSection.root";
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

          std::unique_ptr<PlotUtils::MnvH1D> plasticXsecSmart;
          std::unique_ptr<PlotUtils::MnvH1D> plasticXsecMCSmart;
          std::cout<<"Here111\n";
          std::cout<<"tgtname"<<tgtname<<"\n";
          if (tgtname == "Iron" || tgtname == "2026" || tgtname == "3026")
          {
            titlebase = "Fe";
            plasticXsecSmart=util::GetIngredientSmart<PlotUtils::MnvH1D>(*trackerFile, prefix+"_Fe_CrossSection");
            plasticXsecMCSmart=util::GetIngredientSmart<PlotUtils::MnvH1D>(*trackerFile, prefix+"_Fe_CrossSectionSimulated");
          }
          if (tgtname == "Lead")
          {
            titlebase="Pb";
            plasticXsecSmart=util::GetIngredientSmart<PlotUtils::MnvH1D>(*trackerFile, prefix+"_Pb_CrossSection");
            plasticXsecMCSmart=util::GetIngredientSmart<PlotUtils::MnvH1D>(*trackerFile, prefix+"_Pb_CrossSectionSimulated");
          }
          if (tgtname == "Carbon")
          {
            titlebase="C";
            plasticXsecSmart=util::GetIngredientSmart<PlotUtils::MnvH1D>(*trackerFile, prefix+"_C_CrossSection");
            plasticXsecMCSmart=util::GetIngredientSmart<PlotUtils::MnvH1D>(*trackerFile, prefix+"_C_CrossSectionSimulated");
          }
           if (tgtname == "Water")
          {
            titlebase="H_{2}O";
            plasticXsecSmart=util::GetIngredientSmart<PlotUtils::MnvH1D>(*trackerFileNoDaisy, "crossSection");
            plasticXsecMCSmart=util::GetIngredientSmart<PlotUtils::MnvH1D>(*trackerFileNoDaisy, "simulatedCrossSection");
          } 
          std::string ratioStr = titlebase+"/CH Ratio";
          titlebase = ratioStr + " in "+prefix;
          std::cout<<"Here112\n";
          PlotUtils::MnvH1D* plasticxsec = plasticXsecSmart.get();
          PlotUtils::MnvH1D* plasticxsecMC = plasticXsecMCSmart.get();
          if (prefix=="pTmu")
          {
            plasticxsec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix=="pZmu")
          {
            plasticxsec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix=="BjorkenX")
          {
            plasticxsec->GetXaxis()->SetRange(1, 6);
          }
          else if (prefix=="Erecoil")
          {
          } 
          else if (prefix=="Emu")
          {
            plasticxsec->GetXaxis()->SetRange(1, 11);
          } 
          std::cout<<"Here113\n";

          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
          //plotter.axis_maximum = 0.4;
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;
          PlotUtils::MnvH1D* dataRatio = crossSectionTuned->Clone("DataRatio");
          PlotUtils::MnvH1D* mcRatio = simulatedCrossSection->Clone("MCRatio");
          dataRatio->Divide(dataRatio, plasticxsec);
          mcRatio->Divide(mcRatio, plasticxsecMC);

          dataRatio->SetMaximum( 1.5 );
          dataRatio->SetMinimum( 0.5 );

          dataRatio->GetYaxis()->SetTitle(ratioStr.c_str());
          dataRatio->GetXaxis()->SetTitleFont(plotter.axis_title_font_x);
          dataRatio->GetYaxis()->SetTitleFont(plotter.axis_title_font_y);
          dataRatio->GetXaxis()->SetTitleSize(plotter.axis_title_size_x);
          dataRatio->GetYaxis()->SetTitleSize(plotter.axis_title_size_y);
          dataRatio->GetXaxis()->SetLabelFont(plotter.axis_label_font);
          dataRatio->GetYaxis()->SetLabelFont(plotter.axis_label_font);
          dataRatio->GetXaxis()->SetLabelSize(plotter.axis_label_size);
          dataRatio->GetYaxis()->SetLabelSize(plotter.axis_label_size);
          dataRatio->GetXaxis()->CenterTitle(kTRUE);
          dataRatio->SetMarkerStyle(plotter.ratio_marker);
          dataRatio->SetMarkerSize(plotter.ratio_marker_size);
          dataRatio->SetLineWidth(plotter.ratio_line_width);
          dataRatio->SetLineColor(plotter.ratio_color);
          dataRatio->GetXaxis()->SetNdivisions(509);
          dataRatio->GetYaxis()->SetNdivisions(5, 0, 0);
          //dataRatio->DrawCopy("E1 X0");
          TH1D dataRatioTotalErr = dataRatio->GetCVHistoWithError();

          dataRatioTotalErr.SetMaximum( 1.5 );
          dataRatioTotalErr.SetMinimum( 0.5 );
          dataRatioTotalErr.GetYaxis()->SetTitle(ratioStr.c_str());
          dataRatioTotalErr.GetXaxis()->SetTitleFont(plotter.axis_title_font_x);
          dataRatioTotalErr.GetYaxis()->SetTitleFont(plotter.axis_title_font_y);
          dataRatioTotalErr.GetXaxis()->SetTitleSize(plotter.axis_title_size_x);
          dataRatioTotalErr.GetYaxis()->SetTitleSize(plotter.axis_title_size_y);
          dataRatioTotalErr.GetXaxis()->SetLabelFont(plotter.axis_label_font);
          dataRatioTotalErr.GetYaxis()->SetLabelFont(plotter.axis_label_font);
          dataRatioTotalErr.GetXaxis()->SetLabelSize(plotter.axis_label_size);
          dataRatioTotalErr.GetYaxis()->SetLabelSize(plotter.axis_label_size);
          dataRatioTotalErr.GetXaxis()->CenterTitle(kTRUE);
          dataRatioTotalErr.SetMarkerStyle(plotter.ratio_marker);
          dataRatioTotalErr.SetMarkerSize(plotter.ratio_marker_size);
          dataRatioTotalErr.SetLineWidth(plotter.ratio_line_width);
          dataRatioTotalErr.SetLineColor(plotter.ratio_color);
          dataRatioTotalErr.GetXaxis()->SetNdivisions(509);
          dataRatioTotalErr.GetYaxis()->SetNdivisions(5, 0, 0);
          dataRatioTotalErr.DrawCopy("E1 X0");

          TH1D dataRatioStatErr = dataRatio->GetCVHistoWithStatError();
          dataRatioStatErr.Draw("E1 X0 SAME");
          mcRatio->DrawCopy("HIST SAME");
          const TAxis *axis = dataRatio->GetXaxis();
          double lowX  = axis->GetBinLowEdge( axis->GetFirst() );
          double highX = axis->GetBinUpEdge(  axis->GetLast() );
          TLine line;
          line.SetLineStyle(2);
          line.SetLineWidth(3);
          line.SetLineColor(36);
          line.DrawLine(lowX, 1., highX, 1.); //creates a new line which is owned by gPad

          TLegend *leg = new TLegend(0.6,0.7,0.75,0.9);
          leg->SetBorderSize(0);   // no border
          leg->SetFillColorAlpha(0, 0);    // transparent background
          leg->AddEntry(dataRatio, "Data", "ple");
          leg->AddEntry(mcRatio,"Tune 4.3.0 (CV)","l");
          leg->Draw();
          //plotter.DrawDataMCRatio(crossSectionTuned, plasticxsec, 1, true, true, 0.5, 1.5);
          plotter.AddHistoTitle(titlebase.c_str(), 0.039);
          can.Update();
          plotter.MultiPrint(&can, tgtname+"CHRatio_"+prefix, "png,pdf,C");
        }

        if (prefix!="segment" && prefix!="ANNConfidence")
        { //Plotting cross-sections
          std::string modeltgt;
          if (tgtname == "Iron" || tgtname == "2026" || tgtname == "3026" || tgtname == "5026") modeltgt = "Iron";
          else if (tgtname == "Lead" || tgtname == "2082" || tgtname == "3082" || tgtname == "4082" || tgtname == "5082") modeltgt = "Lead";
          else if (tgtname == "Carbon" || tgtname == "3006") modeltgt = "Carbon";
          else if (tgtname == "Water" || tgtname == "WaterEmpty" || tgtname == "WaterFull" || tgtname == "6000") modeltgt = "Water";
          bool isLead = false;
          if (modeltgt=="Lead") isLead = true;
          
          //GENIE3 G18_02a
          std::string G18_02a_02_11a_string = std::string("/exp/minerva/app/users/alhart/MAT_AL9/MINERvA-101-Cross-Section/modelComparions/")+nutypestr+"/"+modeltgt+"_GENIE_G18_02a_02_11a_cross_sections.root";
          auto G18_02a_02_11a_file = TFile::Open(G18_02a_02_11a_string.c_str());
          if (!G18_02a_02_11a_file)
          {
            std::cerr << "Failed to open model comparison file " << G18_02a_02_11a_string << "\n";
            return 2;
          }
          //GENIE3 G18_02a
          std::string G18_02b_02_11a_string = std::string("/exp/minerva/app/users/alhart/MAT_AL9/MINERvA-101-Cross-Section/modelComparions/")+nutypestr+"/"+modeltgt+"_GENIE_G18_02b_02_11a_cross_sections.root";
          auto G18_02b_02_11a_file = TFile::Open(G18_02b_02_11a_string.c_str());
          if (!G18_02b_02_11a_file)
          {
            std::cerr << "Failed to open model comparison file " << G18_02b_02_11a_string << "\n";
            return 2;
          }
          //GENIE3 G18_10a
          std::string G18_10a_02_11a_string = std::string("/exp/minerva/app/users/alhart/MAT_AL9/MINERvA-101-Cross-Section/modelComparions/")+nutypestr+"/"+modeltgt+"_GENIE_G18_10a_02_11a_cross_sections.root";
          auto G18_10a_02_11a_file = TFile::Open(G18_10a_02_11a_string.c_str());
          if (!G18_10a_02_11a_file)
          {
            std::cerr << "Failed to open model comparison file " << G18_10a_02_11a_string << "\n";
            return 2;
          }
          //GENIE3 G18_10b
          std::string G18_10b_02_11a_string = std::string("/exp/minerva/app/users/alhart/MAT_AL9/MINERvA-101-Cross-Section/modelComparions/")+nutypestr+"/"+modeltgt+"_GENIE_G18_10b_02_11a_cross_sections.root";
          auto G18_10b_02_11a_file = TFile::Open(G18_10b_02_11a_string.c_str());
          if (!G18_10b_02_11a_file)
          {
            std::cerr << "Failed to open model comparison file " << G18_10b_02_11a_string << "\n";
            return 2;
          }
          //NEUT LFG ma105
          std::string NEUT_tune_LFG_maqe1_05_string = std::string("/exp/minerva/app/users/alhart/MAT_AL9/MINERvA-101-Cross-Section/modelComparions/")+nutypestr+"/"+modeltgt+"_NEUT_SF_ma103_cross_sections.root";
          auto NEUT_tune_LFG_maqe1_05_file = TFile::Open(NEUT_tune_LFG_maqe1_05_string.c_str());
          if (!NEUT_tune_LFG_maqe1_05_file)
          {
            std::cerr << "Failed to open model comparison file " << NEUT_tune_LFG_maqe1_05_string << "\n";
            return 2;
          }
          //NEUT SFG ma103
          std::string NEUT_tune_SFG_maqe1_03_string = std::string("/exp/minerva/app/users/alhart/MAT_AL9/MINERvA-101-Cross-Section/modelComparions/")+nutypestr+"/"+modeltgt+"_NEUT_LFG_ma105_cross_sections.root";
          auto NEUT_tune_SFG_maqe1_03_file = TFile::Open(NEUT_tune_SFG_maqe1_03_string.c_str());
          if (!NEUT_tune_SFG_maqe1_03_file && !isLead)
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


          TObjArray* arr = new TObjArray();
          auto tempsimxsec = simulatedCrossSection->Clone();
          auto tempxsec = crossSectionTuned->Clone();

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


          //Used in ratio drawing later
          PlotUtils::MnvH1D * cv_hist_ratio= (PlotUtils::MnvH1D*)tempxsec->Clone("MnvTune 4.3.0 (CV)");
          PlotUtils::MnvH1D * hist_G18_02a_ratio= (PlotUtils::MnvH1D*)tempxsec->Clone("GENIE3 G18_02a");
          PlotUtils::MnvH1D * hist_G18_02b_ratio= (PlotUtils::MnvH1D*)tempxsec->Clone("GENIE3 G18_02b");
          PlotUtils::MnvH1D * hist_G18_10a_ratio= (PlotUtils::MnvH1D*)tempxsec->Clone("GENIE3 G18_10a");
          PlotUtils::MnvH1D * hist_G18_10b_ratio= (PlotUtils::MnvH1D*)tempxsec->Clone("GENIE3 G18_10b");
          PlotUtils::MnvH1D * hist_NEUT_LFG_ratio= (PlotUtils::MnvH1D*)tempxsec->Clone("NEUT LFG ma105");

          cv_hist_ratio->Divide(cv_hist_ratio, tempsimxsec);
          hist_G18_02a_ratio->Divide(hist_G18_02a_ratio, &hist_G18_02a);
          hist_G18_02b_ratio->Divide(hist_G18_02b_ratio, &hist_G18_02b);
          hist_G18_10a_ratio->Divide(hist_G18_10a_ratio, &hist_G18_10a);
          hist_G18_10b_ratio->Divide(hist_G18_10b_ratio, &hist_G18_10b);
          hist_NEUT_LFG_ratio->Divide(hist_NEUT_LFG_ratio, &hist_NEUT_LFG);

          cv_hist_ratio->SetLineColor(46);
          hist_G18_02a_ratio->SetLineColor(TColor::GetColor("#332288")); 
          hist_G18_02b_ratio->SetLineColor(TColor::GetColor("#44AA99")); 
          hist_G18_10a_ratio->SetLineColor(TColor::GetColor("#88CCEE")); 
          hist_G18_10b_ratio->SetLineColor(TColor::GetColor("#DDCC77")); 
          hist_NEUT_LFG_ratio->SetLineColor(TColor::GetColor("#AA4499")); 

          cv_hist_ratio->SetMarkerColor(46);
          hist_G18_02a_ratio->SetMarkerColor(TColor::GetColor("#332288")); 
          hist_G18_02b_ratio->SetMarkerColor(TColor::GetColor("#44AA99")); 
          hist_G18_10a_ratio->SetMarkerColor(TColor::GetColor("#88CCEE")); 
          hist_G18_10b_ratio->SetMarkerColor(TColor::GetColor("#DDCC77")); 
          hist_NEUT_LFG_ratio->SetMarkerColor(TColor::GetColor("#AA4499")); 

          cv_hist_ratio->SetMinimum(0.0001);
          cv_hist_ratio->SetMaximum(1.9999);

          plotter.ApplyNextLineStyle(cv_hist_ratio, true, false);
          plotter.ApplyNextLineStyle(cv_hist_ratio, false, false);
          plotter.ApplyNextLineStyle(hist_G18_02a_ratio, false, false);
          plotter.ApplyNextLineStyle(hist_G18_02b_ratio, false, false);
          plotter.ApplyNextLineStyle(hist_G18_10a_ratio, false, false);
          plotter.ApplyNextLineStyle(hist_G18_10b_ratio, false, false);
          plotter.ApplyNextLineStyle(hist_NEUT_LFG_ratio, false, false);

          PlotUtils::MnvH1D * hist_NEUT_SFG_ratio = nullptr;
          PlotUtils::MnvH1D hist_NEUT_SFG;
          if (!isLead)
          {
            TH1D* NEUT_tune_SFG_maqe1_03_xsec=(TH1D*)util::GetIngredientSmart<TH1D>(*NEUT_tune_SFG_maqe1_03_file, prefix+"_rw").get()->Clone("NEUT SF ma=103");
            NEUT_tune_SFG_maqe1_03_xsec->Scale(1., "width");
            hist_NEUT_SFG = PlotUtils::MnvH1D(*NEUT_tune_SFG_maqe1_03_xsec);
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
          }


          

          std::string title = tmp + " " + tgtString + " Model Comparisons";
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

          if (prefix=="BjorkenX")
          {
            pad1->SetLogx();
            pad2->SetLogx();
          }


          if (tgtname == "Water")
          {
            std::stringstream POTstrFilled, POTstrEmpty;
            POTstrFilled << std::fixed << std::scientific << std::setprecision(2) << waterFilledPOTData;
            POTstrEmpty << std::fixed << std::scientific << std::setprecision(2) << waterEmptyPOTData;
            std::cout << "waterFilledPOTData: "<< waterFilledPOTData <<std::endl;
            std::cout << "waterEmptyPOTData: "<< waterEmptyPOTData <<std::endl;
            std::string potstrFilled = std::string("Water Filled POT ")+POTstrFilled.str();
            std::string potstrEmpty = std::string("Water Empty POT: ")+POTstrEmpty.str();// + " (Scaled to match filled)";
            plotter.AddPlotLabel(potstrFilled.c_str() , 0.18, 0.93, 0.025, 1, 42, 13, 0);
            plotter.AddPlotLabel(potstrEmpty.c_str() , 0.18, 0.90, 0.025, 1, 42, 13, 0);
          }
          else
          {
            std::stringstream POTstr;
            POTstr << std::fixed << std::scientific << std::setprecision(2) << dataPOT;
            std::string potstr = std::string("POT Used: ")+POTstr.str();
            plotter.AddPlotLabel(potstr.c_str() , 0.18, 0.93, 0.025, 1, 42, 13, 0);
          }
          plotter.WritePreliminary( 0.5, 0.15, 0.035, true);
          plotter.AddPlotLabel(title.c_str() , 0.5, 0.97, 25, 1, 43, 22, 0);
          plotter.AddPlotLabel(y_title.c_str() , 0.08, 0.9, 0.03, 1, 42, 33, 90);
          plotter.AddPlotLabel("Data/MC" , 0.08, 0.3, 0.03, 1, 42, 33, 90);
          plotter.AddPlotLabel(xaxislabel.c_str() , 0.9, 0.08, 0.03, 1, 42, 33, 0);

          //Interaction Type
          pad1->cd(); 
          plotter.axis_minimum=0.001;

          //tempxsec->Draw("HIST");
          //tempsimxsec->Draw("HIST SAME");
          //plotter.DrawDataMC(tempxsec, tempsimxsec);
          tempxsec->GetXaxis()->SetLabelFont(43);
          tempxsec->GetXaxis()->SetLabelSize(25);
          std::string legPos = (prefix == "BjorkenX") ? "TL" : "TR";
          plotter.DrawDataMCVariations(tempxsec, arr, 1, legPos, true, true, false, false, true, false);
          pad2->cd();
          plotter.axis_minimum=defaultmin;

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

          //plotter.DrawDataMCRatio(tempxsec, tempsimxsec, 1, false, true, 0.5, 1.5, "");
          //plotter.DrawDataMCRatio(tempxsec, &hist_G18_02a, 1, false, true, 0.5, 1.5, "");
          //can.SetLogx();
          pad1->cd();
          can.cd();
          plotter.MultiPrint(&can, tgtname+"_"+prefix+"_ModelComparisonsCrossSection", "png,pdf,C");
        }

        std::string anuFilepath = "./Antineutrino_"+tgtname+"_"+prefix+"_crossSection.root";
        if (pdg>0 &&std::filesystem::exists(anuFilepath))
        {
          auto anuFile = TFile::Open(anuFilepath.c_str());
          if (!anuFile)
          {
            std::cerr << "Failed to open antineutrino file\n";
            return 2;
          }
          auto anuXSecSmart =util::GetIngredientSmart<PlotUtils::MnvH1D>(*anuFile, "crossSectionTuned");
          auto anuXSec =anuXSecSmart.get();

          std::cout<<"Antineutrino extracted cross sections found - plotting ratios\n";
          std::string titlebase;
          std::cout<<"here2221\n";
          if (tgtname == "Iron" || tgtname == "2026" || tgtname == "3026")
          {
            titlebase = "Fe";
          }
          if (tgtname == "Lead")
          {
            titlebase="Pb";
          }
          if (tgtname == "Carbon")
          {
            titlebase="C";
          }
           if (tgtname == "Water")
          {
            titlebase="H_{2}O";
          } 
          titlebase+=" Neutrino/Antineutrino Ratio in "+prefix;
          std::cout<<"here2222\n";
          if (prefix=="pTmu")
          {
            anuXSec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix=="pZmu")
          {
            anuXSec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix=="BjorkenX")
          {
            anuXSec->GetXaxis()->SetRange(1, 6);
          }
          else if (prefix=="Erecoil")
          {
          } 
          else if (prefix=="Emu")
          {
            anuXSec->GetXaxis()->SetRange(1, 11);
          } 
          std::cout<<"here2223\n";

          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
          //plotter.axis_maximum = 0.4;
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;
          plotter.DrawDataMCRatio(crossSectionTuned, anuXSec, 1, true, true);
          plotter.AddHistoTitle(titlebase.c_str(), 0.039);
          can.Update();
          std::cout<<(tgtname+"_NuAnuRatio_"+prefix+".pdf\n");
          plotter.MultiPrint(&can, tgtname+"_NuAnuRatio_"+prefix, "png,pdf,C");
        }

        //Plotting nu/anu mat/ch double ratio
        //std::string daisyFilepath = "./"+nutypestr+"_"+prefix+"_Daisy_crossSection.root";
        std::string anudaisyFilepath = "./Antineutrino_"+prefix+"_Daisy_crossSection.root";
        //std::string anuFilepath = "./Antineutrino_"+tgtname+"_"+prefix+"_crossSection.root";
        //std::string trackerNoDaisyFilepath = "./"+nutypestr+"_"+prefix+"_crossSection.root";
        std::string anutrackerNoDaisyFilepath = "./Antineutrino_"+prefix+"_crossSection.root";
        if (pdg>0 &&std::filesystem::exists(anuFilepath)&&std::filesystem::exists(daisyFilepath)&&std::filesystem::exists(anudaisyFilepath))
        {
          auto anuFile = TFile::Open(anuFilepath.c_str());
          if (!anuFile)
          {
            std::cerr << "Failed to open antineutrino file\n";
            return 2;
          }
          auto trackerFile = TFile::Open(daisyFilepath.c_str());
          if (!trackerFile)
          {
            std::cerr << "Failed to open tracker file\n";
            return 2;
          }
          auto anutrackerFile = TFile::Open(anudaisyFilepath.c_str());
          if (!anutrackerFile)
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
          auto anuTrackerFileNoDaisy = TFile::Open(anutrackerNoDaisyFilepath.c_str());
          if (!anuTrackerFileNoDaisy)
          {
            std::cerr << "Failed to open anu tracker file\n";
            return 2;
          }

          auto anuXSecSmart =util::GetIngredientSmart<PlotUtils::MnvH1D>(*anuFile, "crossSectionTuned");
          auto anuXSec =anuXSecSmart.get();

          std::string titlebase;

          std::unique_ptr<PlotUtils::MnvH1D> plasticXsecSmart, anuplasticXsecSmart;
          if (tgtname == "Iron" || tgtname == "2026" || tgtname == "3026")
          {
            titlebase = "Fe";
            plasticXsecSmart=util::GetIngredientSmart<PlotUtils::MnvH1D>(*trackerFile, prefix+"_Fe_CrossSection");
            anuplasticXsecSmart=util::GetIngredientSmart<PlotUtils::MnvH1D>(*anutrackerFile, prefix+"_Fe_CrossSection");
          }
          if (tgtname == "Lead")
          {
            titlebase="Pb";
            plasticXsecSmart=util::GetIngredientSmart<PlotUtils::MnvH1D>(*trackerFile, prefix+"_Pb_CrossSection");
            anuplasticXsecSmart=util::GetIngredientSmart<PlotUtils::MnvH1D>(*anutrackerFile, prefix+"_Pb_CrossSection");
          }
          if (tgtname == "Carbon")
          {
            titlebase="C";
            plasticXsecSmart=util::GetIngredientSmart<PlotUtils::MnvH1D>(*trackerFile, prefix+"_C_CrossSection");
            anuplasticXsecSmart=util::GetIngredientSmart<PlotUtils::MnvH1D>(*anutrackerFile, prefix+"_C_CrossSection");
          }
          if (tgtname == "Water")
          {
            titlebase="H_{2}O";
            plasticXsecSmart=util::GetIngredientSmart<PlotUtils::MnvH1D>(*trackerFileNoDaisy, "crossSection");
            anuplasticXsecSmart=util::GetIngredientSmart<PlotUtils::MnvH1D>(*anuTrackerFileNoDaisy, "crossSection");
          } 
          PlotUtils::MnvH1D* plasticxsec = plasticXsecSmart.get();
          PlotUtils::MnvH1D* anuplasticxsec = anuplasticXsecSmart.get();

          std::cout<<"Antineutrino extracted cross sections found - plotting ratios\n";

          std::string title=" #frac{"+titlebase+"}{Tracker}#nu  / #frac{"+titlebase+"}{Tracker}#bar{nu} Ratio in "+prefix;
          std::cout<<"here2222\n";
          if (prefix=="pTmu")
          {
            anuXSec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix=="pZmu")
          {
            anuXSec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix=="BjorkenX")
          {
            anuXSec->GetXaxis()->SetRange(1, 6);
          }
          else if (prefix=="Erecoil")
          {
          } 
          else if (prefix=="Emu")
          {
            anuXSec->GetXaxis()->SetRange(1, 11);
          } 
          std::cout<<"here2223\n";

          auto nuRatio = crossSectionTuned->Clone();
          nuRatio->Divide(nuRatio, plasticxsec);
          auto anuRatio = anuXSec->Clone();
          anuRatio->Divide(anuRatio, anuplasticxsec);
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
          //plotter.axis_maximum = 0.4;
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;
          plotter.DrawDataMCRatio(nuRatio, anuRatio, 1, true, true);
          plotter.AddHistoTitle(title.c_str(), 0.039);
          can.Update();
          std::cout<<(tgtname+"_NuAnuDoubleRatio_"+prefix+".pdf\n");

          plotter.MultiPrint(&can, tgtname+"_NuAnuDoubleRatio_"+prefix, "png,pdf,C");
          plotter.DrawDataMCRatio(crossSectionTuned, plasticxsec, 1, true, true);
          plotter.AddHistoTitle(title.c_str(), 0.039);
          can.Update();
          std::cout<<(tgtname+"_NuAnuDoubleRatio1_"+prefix+".pdf\n");
          plotter.MultiPrint(&can, tgtname+"_NuAnuDoubleRatio1_"+prefix, "png,pdf,C");
          plotter.DrawDataMCRatio(anuXSec, anuplasticxsec, 1, true, true);
          plotter.AddHistoTitle(title.c_str(), 0.039);
          can.Update();
          std::cout<<(tgtname+"_NuAnuDoubleRatio2_"+prefix+".pdf\n");
          plotter.MultiPrint(&can, tgtname+"_NuAnuDoubleRatio2_"+prefix, "png,pdf,C");
        }


        outFile->Close();
      }
      catch (const std::runtime_error &e)
      {
        std::cerr << "Failed to extract a cross section for prefix " << prefix << " and target " << tgtname << " : " << e.what() << "\n";
        return 4;
        // break;
      }

      delete fluxIntReweighted;

      delete flux;
      delete folded;
      delete migration;
      //delete effNum;
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
  }
  return 0;
}
