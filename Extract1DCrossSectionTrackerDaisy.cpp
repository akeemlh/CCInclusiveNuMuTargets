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

  std::vector<std::string> crossSectionPrefixes = {"tracker_pTmu", "tracker_pZmu", "tracker_BjorkenX", "tracker_Erecoil", "tracker_Emu"};

  double mcPOT = 0;
  double dataPOT = 0;

  std::vector<std::pair<std::string, double>> playlistPOTpair; //Used for daisy reweight POT normalising fluxes

  for (const auto &prefix : crossSectionPrefixes)
  {
    //if (!(prefix == "tracker_Erecoil" || prefix == "tracker_pTmu")) continue; //Used for testing with only subset of prefixes
    std::cout << "Currently working on variable: " << prefix << std::endl;
    try
    {
      //Flux parameters
      int n_flux_universes = 100; // Is this right
      const bool use_nue_constraint = true;
      const std::string project_dir = "targets_2345_jointNueIMD";
      double min_energy = 0;
      double max_energy = 100;
      
      PlotUtils::MnvH1D *fluxIntReweighted = new PlotUtils::MnvH1D();
      PlotUtils::MnvH1D *DaisyEffNum[12], *DaisyFolded[12];
      PlotUtils::MnvH2D *DaisyMigration[12];
      std::vector<PlotUtils::MnvH1D*> DaisyBackgrounds[12];
      PlotUtils::MnvH1D *DaisyEffDenom[12], *DaisyEffDenom2P2H[12], *DaisyEffDenomDIS[12], *DaisyEffDenomRES[12], *DaisyEffDenomQE[12], *DaisyEffDenomOther[12];

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

      std::vector<PlotUtils::MnvH1D *> backgrounds;

      for (int c = 0; c<dirs.size(); c++)
      {
        std::cout<<"Investigating directory " << dirs[c] << " which is " << c+1 <<"/"<<dirs.size() <<" playlists identified" <<std::endl;
        std::string datapath = dirs[c] + "/runEventLoopTrackerData.root";
        std::string mcpath = dirs[c] + "/runEventLoopTrackerMC.root";

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

        if (pdg!=nupdg)
        {
          std::cout<<"Skipping this set of files because this playlist pdg is " << nupdg << " but you specified " << pdg << std::endl;
          continue;
        }
        util::AddHist(*folded,util::GetIngredient<PlotUtils::MnvH1D>(*dataFile, (std::string("data")), prefix));
        util::AddHist(*migration,util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("migration")), prefix));
        util::AddHist(*effNum,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_numerator")), prefix));
        util::AddHist(*effDenom,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator")), prefix));
        util::AddHist(*effDenom2P2H,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator_intChannels_2p2h")), prefix));
        util::AddHist(*effDenomDIS,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator_intChannels_DIS")), prefix));
        util::AddHist(*effDenomRES,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator_intChannels_RES")), prefix));
        util::AddHist(*effDenomQE,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator_intChannels_QE")), prefix));
        util::AddHist(*effDenomOther,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator_intChannels_Other")), prefix));

        util::AddHist(*BackgroundWrongSign,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_Wrong_Sign")), prefix));
        util::AddHist(*BackgroundNC,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_NC")), prefix));
        util::AddHist(*BackgroundOther,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_Other")), prefix));

        backgrounds.push_back(BackgroundWrongSign);
        backgrounds.push_back(BackgroundNC);
        backgrounds.push_back(BackgroundOther);        

        mcPOT += util::GetIngredient<TParameter<double>>(*mcFile, "POTUsed")->GetVal();
        double tempDataPOT = util::GetIngredient<TParameter<double>>(*dataFile, "POTUsed")->GetVal();
        dataPOT+=tempDataPOT;

        playlistPOTpair.push_back(std::make_pair(playlistUsed, tempDataPOT));

        PlotUtils::FluxReweighter frw_temp = PlotUtils::FluxReweighter( nupdg, use_nue_constraint, playlistUsed, PlotUtils::FluxReweighter::gen2thin, PlotUtils::FluxReweighter::g4numiv6, n_flux_universes );
        auto tempIntFlux = frw_temp.GetIntegratedFluxReweighted(nupdg, effDenom, min_energy, max_energy, true)->Clone();
        tempIntFlux->Scale(tempDataPOT);
        util::AddHist(*fluxIntReweighted,tempIntFlux);

        dataFile->Close();
        mcFile->Close();

        //***********************************************
        //Daisy reweight
        //Getting ingredients
        //***********************************************
        for (int petal=0; petal<12; petal++){
          std::string datapathdaisy = dirs[c] + "/runEventLoopTrackerData_petal_"+std::to_string(petal)+".root";
          std::string mcpathdaisy = dirs[c] + "/runEventLoopTrackerMC_petal_"+std::to_string(petal)+".root";
          bool petalFilesExist = std::filesystem::exists(datapathdaisy) && std::filesystem::exists(mcpathdaisy);
          doDaisy = doDaisy && petalFilesExist;
        }
        doDaisy=false;
        if (!doDaisy) std::cout<<"Will not do daisy petal reweight in this analysis\n";
        for (int petal=0; petal<12 && doDaisy; petal++){
          std::string datapathdaisy = dirs[c] + "/runEventLoopTrackerData_petal_"+std::to_string(petal)+".root";
          std::string mcpathdaisy = dirs[c] + "/runEventLoopTrackerMC_petal_"+std::to_string(petal)+".root";
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


          util::AddHist(*DaisyEffNum[petal], util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("efficiency_numerator")), prefix));
          util::AddHist(*DaisyEffDenom[petal], util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("efficiency_denominator")), prefix));
          util::AddHist(*DaisyMigration[petal], util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, (std::string("migration")), prefix));
          util::AddHist(*DaisyFolded[petal], util::GetIngredient<PlotUtils::MnvH1D>(*dataDaisyFile, (std::string("data")), prefix));
          DaisyBackgrounds[petal].push_back(util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("background_Wrong_Sign")), prefix));
          DaisyBackgrounds[petal].push_back(util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("background_NC")), prefix));
          DaisyBackgrounds[petal].push_back(util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("background_Other")), prefix));
          //Backgrounds
          util::AddHist(*DaisyEffDenom2P2H[petal], util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("efficiency_denominator_intChannels_2p2h")), prefix));
          util::AddHist(*DaisyEffDenomDIS[petal], util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("efficiency_denominator_intChannels_DIS")), prefix));
          util::AddHist(*DaisyEffDenomRES[petal], util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("efficiency_denominator_intChannels_RES")), prefix));
          util::AddHist(*DaisyEffDenomQE[petal], util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("efficiency_denominator_intChannels_QE")), prefix));
          util::AddHist(*DaisyEffDenomOther[petal], util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, (std::string("efficiency_denominator_intChannels_Other")), prefix));

          dataDaisyFile->Close();
          mcDaisyFile->Close();

        }
      }

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
      for (int petal=0; petal<12 && doDaisy; petal++) DaisyFolded[petal]->AddMissingErrorBandsAndFillWithCV(*(DaisyMigration[petal]));

      // Basing my unfolding procedure for a differential cross section on Alex's MINERvA 101 talk at https://minerva-docdb.fnal.gov/cgi-bin/private/RetrieveFile?docid=27438&filename=whatsACrossSection.pdf&version=1

      // TODO: Remove these debugging plots when done
      auto toSubtract = std::accumulate(std::next(backgrounds.begin()), backgrounds.end(), (*backgrounds.begin())->Clone(),
                                        [](auto sum, const auto hist)
                                        {
                                          sum->Add(hist);
                                          return sum;
                                        });
      Plot(*toSubtract, "BackgroundSum", prefix);

      auto bkgSubtracted = std::accumulate(backgrounds.begin(), backgrounds.end(), folded->Clone(),
                                            [mcPOT, dataPOT](auto sum, const auto hist)
                                            {
                                              std::cout << "Subtracting " << hist->GetName() << " scaled by " << -dataPOT / mcPOT << " from " << sum->GetName() << "\n";
                                              sum->Add(hist, -dataPOT / mcPOT);
                                              return sum;
                                            });
      Plot(*bkgSubtracted, "backgroundSubtracted", prefix);
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
          auto unfoldedDaisy = UnfoldHist(bkgSubtractedDaisy, DaisyMigration[petal], nIterations);
          outFileDaisy->cd();
          unfoldedDaisy->Write((prefix+"_unfoldedDaisy_"+petal));
          if(!unfoldedDaisy) throw std::runtime_error(std::string("Failed to unfold ") + DaisyFolded[petal]->GetName() + " using " + DaisyMigration[petal]->GetName());

          DaisyEffNum[petal]->Divide(DaisyEffNum[petal],DaisyEffDenom[petal]);
          unfoldedDaisy->Divide(unfoldedDaisy, DaisyEffNum[petal]);

          daisy_petal_hists[petal]=unfoldedDaisy->Clone();

          daisy_petal_hists_eff_denom[petal]= DaisyEffDenom[petal];
          daisy_petal_hists_eff_denom_2p2h[petal]= DaisyEffDenom2P2H[petal];
          daisy_petal_hists_eff_denom_dis[petal]= DaisyEffDenomDIS[petal];
          daisy_petal_hists_eff_denom_res[petal]= DaisyEffDenomRES[petal];
          daisy_petal_hists_eff_denom_qe[petal]= DaisyEffDenomQE[petal];
          daisy_petal_hists_eff_denom_other[petal]= DaisyEffDenomOther[petal];

        }
      
        PlotUtils::MnvH1D *fluxIntegralC, *fluxIntegralFe, *fluxIntegralPb; 
        PlotUtils::MnvH1D *DaisyCorrectedC, *DaisyCorrectedFe, *DaisyCorrectedPb; 
        PlotUtils::MnvH1D *DaisyCorrectedSimC, *DaisyCorrectedSimFe, *DaisyCorrectedSimPb; 
        PlotUtils::MnvH1D *DaisyCorrectedSim2p2hC, *DaisyCorrectedSim2p2hFe, *DaisyCorrectedSim2p2hPb; 
        PlotUtils::MnvH1D *DaisyCorrectedSimDISC, *DaisyCorrectedSimDISFe, *DaisyCorrectedSimDISPb; 
        PlotUtils::MnvH1D *DaisyCorrectedSimRESC, *DaisyCorrectedSimRESFe, *DaisyCorrectedSimRESPb; 
        PlotUtils::MnvH1D *DaisyCorrectedSimQEC, *DaisyCorrectedSimQEFe, *DaisyCorrectedSimQEPb; 
        PlotUtils::MnvH1D *DaisyCorrectedSimOtherC, *DaisyCorrectedSimOtherFe, *DaisyCorrectedSimOtherPb; 
        for (auto pp : playlistPOTpair)
        {
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

        fluxIntegralC->Scale(1/dataPOT);
        fluxIntegralFe->Scale(1/dataPOT);
        fluxIntegralPb->Scale(1/dataPOT);
        DaisyCorrectedC->Scale(1/dataPOT);
        DaisyCorrectedFe->Scale(1/dataPOT);
        DaisyCorrectedPb->Scale(1/dataPOT);
        
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

        //Carbon
        {
          DaisyCorrectedC->Write((prefix+"DaisyCorrectedC").c_str());
          crossSectionC = normalize(DaisyCorrectedC, fluxIntegralC, nnucleonsData, dataPOT);
          outFileDaisy->cd();
          crossSectionC->Write((prefix+"_C_CrossSection").c_str());
          crossSectionSimC = normalize(DaisyCorrectedSimC, fluxIntegralC, nnucleons, mcPOT);
          crossSectionSimC->Write((prefix+"_C_CrossSectionSimulated").c_str());
          crossSectionSimC2p2h = normalize(DaisyCorrectedSim2p2hC, fluxIntegral, nnucleons, mcPOT);
          crossSectionSimC2p2h->Write((prefix+"_C_CrossSectionSimulated2p2h").c_str());
          crossSectionSimCDIS = normalize(DaisyCorrectedSimDISC, fluxIntegral, nnucleons, mcPOT);
          crossSectionSimCDIS->Write((prefix+"_C_CrossSectionSimulatedDIS").c_str());
          crossSectionSimCRes = normalize(DaisyCorrectedSimRESC, fluxIntegral, nnucleons, mcPOT);
          crossSectionSimCRes->Write((prefix+"_C_CrossSectionSimulatedRES").c_str());
          crossSectionSimCQE = normalize(DaisyCorrectedSimQEC, fluxIntegral, nnucleons, mcPOT);
          crossSectionSimCQE->Write((prefix+"_C_CrossSectionSimulatedQE").c_str());
          crossSectionSimCOther= normalize(DaisyCorrectedSimOtherC, fluxIntegral, nnucleons, mcPOT);
          crossSectionSimCOther->Write((prefix+"_C_CrossSectionSimulatedOther").c_str());
        }

        //Iron
        {
          DaisyCorrectedFe->Write((prefix+"DaisyCorrectedFe").c_str());
          crossSectionFe = normalize(DaisyCorrectedFe, fluxIntegralFe, nnucleonsData, dataPOT);
          outFileDaisy->cd();
          crossSectionFe->Write((prefix+"_Fe_CrossSection").c_str());
          crossSectionSimFe = normalize(DaisyCorrectedSimFe, fluxIntegralFe, nnucleons, mcPOT);
          crossSectionSimFe->Write((prefix+"_Fe_CrossSectionSimulated").c_str());
          /* crossSectionSimFe2p2h = normalize(DaisyCorrectedSim2p2hFe, fluxIntegral, nnucleons, mcPOT);
          crossSectionSimFe2p2h->Write((prefix+"_Fe_CrossSectionSimulated2p2h").c_str());
          crossSectionSimFeDIS = normalize(DaisyCorrectedSimDISFe, fluxIntegral, nnucleons, mcPOT);
          crossSectionSimFeDIS->Write((prefix+"_Fe_CrossSectionSimulatedDIS").c_str());
          crossSectionSimFeRes = normalize(DaisyCorrectedSimRESFe, fluxIntegral, nnucleons, mcPOT);
          crossSectionSimFeRes->Write((prefix+"_Fe_CrossSectionSimulatedRES").c_str());
          crossSectionSimFeQE = normalize(DaisyCorrectedSimQEFe, fluxIntegral, nnucleons, mcPOT);
          crossSectionSimFeQE->Write((prefix+"_Fe_CrossSectionSimulatedQE").c_str());
          crossSectionSimFeOther= normalize(DaisyCorrectedSimOtherFe, fluxIntegral, nnucleons, mcPOT);
          crossSectionSimFeOther->Write((prefix+"_Fe_CrossSectionSimulatedOther").c_str()); */
        }

        //Lead
        {
          DaisyCorrectedPb->Write((prefix+"DaisyCorrectedPb").c_str());
          crossSectionPb = normalize(DaisyCorrectedPb, fluxIntegralPb, nnucleonsData, dataPOT);
          outFileDaisy->cd();
          crossSectionPb->Write((prefix+"_Pb_CrossSection").c_str());
          crossSectionSimPb = normalize(DaisyCorrectedSimPb, fluxIntegralPb, nnucleons, mcPOT);
          crossSectionSimPb->Write((prefix+"_Pb_CrossSectionSimulated").c_str());
          /* crossSectionSimPb2p2h = normalize(DaisyCorrectedSim2p2hPb, fluxIntegral, nnucleons, mcPOT);
          crossSectionSimPb2p2h->Write((prefix+"_Pb_CrossSectionSimulated2p2h").c_str());
          crossSectionSimPbDIS = normalize(DaisyCorrectedSimDISPb, fluxIntegral, nnucleons, mcPOT);
          crossSectionSimPbDIS->Write((prefix+"_Pb_CrossSectionSimulatedDIS").c_str());
          crossSectionSimPbRes = normalize(DaisyCorrectedSimRESPb, fluxIntegral, nnucleons, mcPOT);
          crossSectionSimPbRes->Write((prefix+"_Pb_CrossSectionSimulatedRES").c_str());
          crossSectionSimPbQE = normalize(DaisyCorrectedSimQEPb, fluxIntegral, nnucleons, mcPOT);
          crossSectionSimPbQE->Write((prefix+"_Pb_CrossSectionSimulatedQE").c_str());
          crossSectionSimPbOther= normalize(DaisyCorrectedSimOtherPb, fluxIntegral, nnucleons, mcPOT);
          crossSectionSimPbOther->Write((prefix+"_Pb_CrossSectionSimulatedOther").c_str()); */
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
          if (prefix=="tracker_pTmu")
          {
            tmp = "Muon p_{T}";
            //tempxsec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix=="tracker_pZmu")
          {
            tmp = "Muon p_{Z}";
            //tempxsec->GetXaxis()->SetRange(1, 13);
          }
          else if (prefix=="tracker_BjorkenX") tmp = "{X}";
          else if (prefix=="tracker_Erecoil") tmp = "E_{recoil}";
          else if (prefix=="tracker_Emu")
          {
            tmp = "E_{#mu}";
            //tempxsec->GetXaxis()->SetRange(1, 11);
          } 
          plotter.axis_title_size_y=0.025;
          plotter.axis_title_offset_y=1.8;

          std::string title;

          title = tmp + " " + "Tracker Cross Section to Daisy Reweighted Cross Section ratio - Lead";
          plotter.DrawDataMCRatio(crossSection, crossSectionPb, 1, "TR", "Data");
          plotter.AddHistoTitle(title.c_str(), 0.028);
          can.Update();
          can.Print(("LeadDaisy"+prefix+"_CrossSectionRatio.pdf").c_str());
          can.Print(("LeadDaisy"+prefix+"_CrossSectionRatio.png").c_str());

          title = tmp + " " + "Tracker Cross Section to Daisy Reweighted Cross Section ratio - Iron";
          plotter.DrawDataMCRatio(crossSection, crossSectionFe, 1, "TR", "Data");
          plotter.AddHistoTitle(title.c_str(), 0.028);
          can.Update();
          can.Print(("IronDaisy"+prefix+"_CrossSectionRatio.pdf").c_str());
          can.Print(("IronDaisy"+prefix+"_CrossSectionRatio.png").c_str());

          title = tmp + " " + "Tracker Cross Section to Daisy Reweighted Cross Section ratio - Carbon";
          plotter.DrawDataMCRatio(crossSection, crossSectionC, 1, "TR", "Data");
          plotter.AddHistoTitle(title.c_str(), 0.028);
          can.Update();
          can.Print(("CarbonDaisy"+prefix+"_CrossSectionRatio.pdf").c_str());
          can.Print(("CarbonDaisy"+prefix+"_CrossSectionRatio.png").c_str());
        }
        { //Plotting cross sections
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
          //plotter.axis_maximum = 0.4;
          TObjArray* arr = new TObjArray();
          auto temp2p2h = crossSectionSimPb2p2h->Clone();
          auto tempdis = crossSectionSimPbDIS->Clone();
          auto tempres = crossSectionSimPbRes->Clone();
          auto tempqe = crossSectionSimPbQE->Clone();
          auto tempother = crossSectionSimPbOther->Clone();
          auto tempxsec = crossSectionPb->Clone();
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
            tmp = "{X}";
            y_title="d#sigma/d{X} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
          }
          else if (prefix.find("Erecoil")!=std::string::npos)
          {
            tmp = "E_{recoil}";
            y_title="d#sigma/d{X} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
          } 
          else if (prefix.find("Emu")!=std::string::npos)
          {
            tmp = "E_{#mu}";
            y_title="d#sigma/dE_{#mu} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            tempxsec->GetXaxis()->SetRange(1, 11);
          } 

          std::string title = tmp;
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;
          plotter.DrawDataStackedMC(tempxsec, arr, 1, "TR", "Data", 2, 1, 3001, "", y_title.c_str());
          plotter.AddHistoTitle(title.c_str(), 0.039);
          can.Update();
          can.Print((prefix+"_PbDaisyCrossSection.pdf").c_str());
          can.Print((prefix+"_PbDaisyCrossSection.png").c_str());
        }
        { //Plotting cross sections
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
          //plotter.axis_maximum = 0.4;
          TObjArray* arr = new TObjArray();
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
          arr->Add(temp2p2h);
          arr->Add(tempdis);
          arr->Add(tempres);
          arr->Add(tempqe);
          arr->Add(tempother);

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
            tmp = "{X}";
            y_title="d#sigma/d{X} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
          }
          else if (prefix.find("Erecoil")!=std::string::npos)
          {
            tmp = "E_{recoil}";
            y_title="d#sigma/d{X} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
          } 
          else if (prefix.find("Emu")!=std::string::npos)
          {
            tmp = "E_{#mu}";
            y_title="d#sigma/dE_{#mu} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            tempxsec->GetXaxis()->SetRange(1, 11);
          } 

          std::string title = tmp;
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;
          plotter.DrawDataStackedMC(tempxsec, arr, 1, "TR", "Data", 2, 1, 3001, "", y_title.c_str());
          plotter.AddHistoTitle(title.c_str(), 0.039);
          can.Update();
          can.Print((prefix+"_FeDaisyCrossSection.pdf").c_str());
          can.Print((prefix+"_FeDaisyCrossSection.png").c_str());
        }
        { //Plotting cross sections
          TCanvas can("Temp");
          // Uncertainty summary
          PlotUtils::MnvPlotter plotter;
          plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
          //plotter.axis_maximum = 0.4;
          TObjArray* arr = new TObjArray();
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
          arr->Add(temp2p2h);
          arr->Add(tempdis);
          arr->Add(tempres);
          arr->Add(tempqe);
          arr->Add(tempother);

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
            tmp = "{X}";
            y_title="d#sigma/d{X} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
          }
          else if (prefix.find("Erecoil")!=std::string::npos)
          {
            tmp = "E_{recoil}";
            y_title="d#sigma/d{X} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
          } 
          else if (prefix.find("Emu")!=std::string::npos)
          {
            tmp = "E_{#mu}";
            y_title="d#sigma/dE_{#mu} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            tempxsec->GetXaxis()->SetRange(1, 11);
          } 

          std::string title = tmp;
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;
          plotter.DrawDataStackedMC(tempxsec, arr, 1, "TR", "Data", 2, 1, 3001, "", y_title.c_str());
          plotter.AddHistoTitle(title.c_str(), 0.039);
          can.Update();
          can.Print((prefix+"_CDaisyCrossSection.pdf").c_str());
          can.Print((prefix+"_CDaisyCrossSection.png").c_str());
        }
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
            tmp = "{X}";
            y_title="d#sigma/d{X} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
          }
          else if (prefix.find("Erecoil")!=std::string::npos)
          {
            tmp = "E_{recoil}";
            y_title="d#sigma/d{X} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
          } 
          else if (prefix.find("Emu")!=std::string::npos)
          {
            tmp = "E_{#mu}";
            y_title="d#sigma/dE_{#mu} (x10^{-39}) (cm^{2}/(GeV/c)^{2}/Nucleon)";
            tempxsec->GetXaxis()->SetRange(1, 11);
          } 

          std::string title = tmp;
          plotter.axis_title_size_y=0.03;
          plotter.axis_title_offset_y=1.8;
          plotter.DrawDataStackedMC(tempxsec, arr, 1, "TR", "Data", 2, 1, 3001, "", y_title.c_str());
          plotter.AddHistoTitle(title.c_str(), 0.039);
          can.Update();
          can.Print((prefix+"_NoDaisyCrossSection.pdf").c_str());
          can.Print((prefix+"_NoDaisyCrossSection.png").c_str());
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
          plotter.DrawNormalizedMigrationHistogram(DaisyMigration[i], false, false, true, false);
          //plotter.AddHistoTitle(title.c_str(), 0.039);
          can.Update();
          //can.Print((tgt+prefix+"_DSSidebandScaled.pdf").c_str());
          can.Print(std::string(prefix+"_petal_"+i+"_migrationTest.png").c_str());
        }

        
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
      can.Print(std::string(prefix+"_migrationTest.png").c_str());
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
