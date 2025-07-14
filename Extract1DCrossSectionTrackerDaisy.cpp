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

  if (argc != 5)
  {
    std::cerr << "Expected 4 arguments, but I got " << argc - 1 << ".\n"
              << "USAGE: ExtractCrossSection <unfolding iterations> <directory> <playlistname-for-flux> <pdg>\n"
              << "e.g: ExtractCrossSection 5 ./ minervame1A 14\n";
    return 1;
  }

  const int nIterations = std::stoi(argv[1]);
  std::string directory = std::string(argv[2]);
  std::string playlistname = std::string(argv[3]);
  int pdg = std::stoi(argv[4]);

  std::string datapath = directory + "/runEventLoopTrackerData.root";
  std::string mcpath = directory + "/runEventLoopTrackerMC.root";

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

  std::vector<std::string> crossSectionPrefixes = {"tracker_pTmu", "tracker_pZmu", "tracker_BjorkenX", "tracker_Erecoil", "tracker_Emu"};

  //Automatically figure out prefixes
  /* for (auto key : *dataFile->GetListOfKeys())
  {
    const std::string keyName = key->GetName();
    if (keyName == "POTUsed")
      continue;
    std::cout << "keyName " << keyName << std::endl;
    const size_t endOfPrefix = keyName.find("_data");
    std::string prefix = keyName.substr(0, endOfPrefix);
    std::cout << "prefix " << prefix << std::endl;
    bool twoDimension = (keyName == "tracker_pTmu_pZmu_data");

    bool alreadyInVector = std::find(crossSectionPrefixes.begin(), crossSectionPrefixes.end(), prefix) != crossSectionPrefixes.end();
    std::cout << "twoDimension " << twoDimension << std::endl;
    std::cout << "alreadyInVector " << alreadyInVector << std::endl;
    if (endOfPrefix != std::string::npos && !alreadyInVector && !twoDimension)
      crossSectionPrefixes.push_back(prefix);
  } */
  const double mcPOT = util::GetIngredient<TParameter<double>>(*mcFile, "POTUsed")->GetVal(),
                dataPOT = util::GetIngredient<TParameter<double>>(*dataFile, "POTUsed")->GetVal();
  std::cout << "DataPOT: " << dataPOT << " mcPOT " << mcPOT << std::endl;
  for (const auto &prefix : crossSectionPrefixes)
  {
    //if (!(prefix == "tracker_Erecoil" || prefix == "tracker_pTmu")) continue; //Used for testing with only subset of prefixes
    std::cout << "Current working on prefix: " << prefix << std::endl;
    try
    {

      PlotUtils::MnvH1D *DaisyEffNum[12], *DaisyEffDenom[12], *DaisyFolded[12];
      PlotUtils::MnvH2D *DaisyMigration[12];
      std::vector<PlotUtils::MnvH1D*> DaisyBackgrounds[12];


      //auto flux = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("reweightedflux_integrated")), prefix);
      //flux->Scale(1.0);
      auto folded = util::GetIngredient<PlotUtils::MnvH1D>(*dataFile, (std::string("data")), prefix);
      //Plot(*folded, "data", prefix);
      auto migration = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("migration")), prefix);
      auto effNum = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_numerator")), prefix);
      auto effDenom = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator")), prefix);
      auto effDenom2P2H = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator_intChannels_2p2h")), prefix);
      auto effDenomDIS = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator_intChannels_DIS")), prefix);
      auto effDenomRES = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator_intChannels_RES")), prefix);
      auto effDenomQE = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator_intChannels_QE")), prefix);
      auto effDenomOther = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator_intChannels_Other")), prefix);

      //auto nNucleons = util::GetIngredient<TParameter<double>>(*mcFile, (std::string("fiducial_nucleons")), prefix); // Dan: Use the same truth fiducial volume for all extractions.  The acceptance correction corrects data back to this fiducial even if the reco fiducial cut is different.
      //double nNucleonsVal = nNucleons->GetVal();
      std::vector<PlotUtils::MnvH1D *> backgrounds;
      for (auto key : *mcFile->GetListOfKeys())
      {
        if (std::string(key->GetName()).find(prefix + "_background_") != std::string::npos)
        {
          std::cout << "Found and adding background: " << key->GetName() << std::endl;
          backgrounds.push_back(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, key->GetName()));
        }
      }

      //***********************************************
      //Daisy reweight
      //Getting ingredients
      //***********************************************
      bool doDaisy = true;
      for (int petal=0; petal<12; petal++){
        std::string datapathdaisy = directory + "/runEventLoopTrackerData_petal_"+std::to_string(petal);
        std::string mcpathdaisy = directory + "/runEventLoopTrackerMCDaisy_petal_"+std::to_string(petal);
        bool petalFilesExist = std::filesystem::exists(datapathdaisy) && std::filesystem::exists(mcpathdaisy);
        doDaisy = doDaisy && petalFilesExist;
      }
      for (int petal=0; petal<12 && doDaisy; petal++){
        std::string datapathdaisy = directory + "/runEventLoopTrackerData_petal_"+std::to_string(petal);
        std::string mcpathdaisy = directory + "/runEventLoopTrackerMCDaisy_petal_"+std::to_string(petal);
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


        DaisyEffNum[petal] = util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, ("Daisy_EffNum_"+std::to_string(petal)).c_str(), prefix);
        DaisyEffDenom[petal] = util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, ("Daisy_EffDenom_"+std::to_string(petal)).c_str(), prefix);
        DaisyMigration[petal] = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, ("Daisy_Migration_"+std::to_string(petal)).c_str(), prefix);
        DaisyFolded[petal] = util::GetIngredient<PlotUtils::MnvH1D>(*dataDaisyFile, ("Daisy_Data_"+std::to_string(petal)).c_str(), prefix);
        
        //Backgrounds
        for(auto key: *mcFile->GetListOfKeys())
        {
          //std::string temp = prefix + "_Daisy_Background_"+std::to_string(petal);
          //std::cout<<"Looking for background: " <<  temp << std::endl;
          if(std::string(key->GetName()).find(prefix + "_Daisy_Background_"+std::to_string(petal)) != std::string::npos)
          {
            std::cout << "Found and adding Daisy background: " << key->GetName() << std::endl;
            DaisyBackgrounds[petal].push_back(util::GetIngredient<PlotUtils::MnvH1D>(*mcDaisyFile, key->GetName()));
          }
        }
        DaisyFolded[petal]->AddMissingErrorBandsAndFillWithCV(*(DaisyMigration[petal]));

        dataDaisyFile->Close();
        mcDaisyFile->Close();

      }





      std::cout << "Test " <<std::endl;
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
      std::cout << "Test 2 " <<std::endl;
      auto outFile = TFile::Open((prefix + "_crossSection.root").c_str(), "RECREATE");
      if (!outFile)
      {
        std::cerr << "Could not create a file called " << prefix + "_crossSection.root" << ".  Does it already exist?\n";
        return 5;
      }

      bkgSubtracted->Write("backgroundSubtracted");

      // d'Aogstini unfolding
      std::cout << "Test 3 " <<std::endl;
      auto unfolded = UnfoldHist(bkgSubtracted, migration, nIterations);
      std::cout << "Test 4 " <<std::endl;
      if (!unfolded)
        throw std::runtime_error(std::string("Failed to unfold ") + folded->GetName() + " using " + migration->GetName());
      std::cout << "Test 5 " <<std::endl;
      Plot(*unfolded, "unfolded", prefix);
      unfolded->Clone()->Write("unfolded"); // TODO: Seg fault first appears when I uncomment this line
      std::cout << "Survived writing the unfolded histogram.\n"
                << std::flush; // This is evidence that the problem is on the final file Write() and not unfolded->Clone()->Write().

      effNum->Divide(effNum, effDenom); // Only the 2 parameter version of MnvH1D::Divide()
                                        // handles systematics correctly.
      Plot(*effNum, "efficiency", prefix);
      unfolded->Divide(unfolded, effNum);
      Plot(*unfolded, "efficiencyCorrected", prefix);

      // double nnucleons = nNucleons->GetVal()/numMergedPlaylists;
      int n_flux_universes = 100; // Is this right
      const bool use_nue_constraint = true;
      const std::string project_dir = "targets_2345_jointNueIMD";
      double min_energy = 0;
      double max_energy = 100;

      PlotUtils::TargetUtils targetInfo;
      double nnucleons = targetInfo.GetTrackerNNucleons(PlotUtils::TargetProp::Tracker::Face, PlotUtils::TargetProp::Tracker::Back, true);
      double nnucleonsData = targetInfo.GetTrackerNNucleons(PlotUtils::TargetProp::Tracker::Face, PlotUtils::TargetProp::Tracker::Back, true);

      PlotUtils::MnvH1D *flux2;
      PlotUtils::MnvH1D *fluxIntegral;
      PlotUtils::MnvH1D *fluxRebinned;

      std::cout<<"pdg: " << pdg << std::endl;
      std::cout<<"use_nue_constraint: " << use_nue_constraint << std::endl;
      std::cout<<"n_flux_universes: " << n_flux_universes << std::endl;
      std::cout<<"min_energy: " << min_energy << std::endl;
      std::cout<<"max_energy: " << max_energy << std::endl;
      std::cout<<"nnucleons: " << nnucleons << std::endl;
      std::cout<<"mcPOT: " << mcPOT << std::endl;
      //PlotUtils::FluxReweighter frw = PlotUtils::flux_reweighter(playlistname, pdg, use_nue_constraint, n_flux_universes);
      PlotUtils::FluxReweighter *frw = new PlotUtils::FluxReweighter( pdg, use_nue_constraint, playlistname, PlotUtils::FluxReweighter::gen2thin, PlotUtils::FluxReweighter::g4numiv6, n_flux_universes );
      //std::cout << "ABC123-1 " << std::endl;
      PlotUtils::MnvH1D *fluxIntReweighted = frw->GetIntegratedFluxReweighted(pdg, simEventRate, min_energy, max_energy, true);
      //fluxIntegral = frw->GetIntegratedTargetFlux(pdg, material, simEventRate, min_energy, max_energy, project_dir);
      //auto &frw2 = PlotUtils::flux_reweighter("minervame1A", pdg, use_nue_constraint, n_flux_universes);
      //flux2 = frw2.GetTargetFluxMnvH1D(pdg, material, project_dir);

      //fluxRebinned = frw->GetRebinnedFluxReweighted_FromInputFlux(flux2, simEventRate); // issue here
      //PlotUtils::MnvH1D *Integrated_fluxGenie = frw->GetIntegratedFluxReweighted_FromInputFlux(flux2, simEventRate, min_energy, max_energy);
      //PlotUtils::MnvH1D *Integrated_fluxGenie2 = frw->GetIntegratedFluxReweighted_FromInputFlux(flux, simEventRate, min_energy, max_energy);
      std::cout << "Test 6 " <<std::endl;
      outFile->cd();
      std::cout << "Test 7 " <<std::endl;
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
      std::cout << "Test 8 " <<std::endl;
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
      std::cout << "Test 9 " <<std::endl;




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
        std::cout << "Test 10 " <<std::endl;
        //auto& frw = PlotUtils::flux_reweighter(playlistname, pdg, true, 100);
        std::map<int, PlotUtils::MnvH1D*> daisy_petal_hists;

        for (int petal=0; petal<12; petal++){
          std::cout << "Test 11 " <<std::endl;
          std::vector<PlotUtils::MnvH1D*> DaisyBackgroundsTemp = DaisyBackgrounds[petal];

          auto toSubtractDaisy = std::accumulate(std::next(DaisyBackgroundsTemp.begin()), DaisyBackgroundsTemp.end(), (*DaisyBackgroundsTemp.begin())->Clone(),
                                          [](auto sum, const auto hist)
                                          {
                                            sum->Add(hist);
                                            return sum;
                                          });
          //Plot(*toSubtract, "BackgroundSum", prefix);
          std::cout << "Test 11.1 " <<std::endl;
          outFileDaisy->cd();
          std::cout << "Test 11.2 " <<std::endl;
          toSubtractDaisy->Write((prefix+"_toSubtractDaisy_"+petal));
          std::cout << "Test 11.3 " <<std::endl;
          auto bkgSubtractedDaisy = std::accumulate(DaisyBackgrounds[petal].begin(), DaisyBackgrounds[petal].end(), DaisyFolded[petal]->Clone(),
                                              [mcPOT, dataPOT](auto sum, const auto hist)
                                              {
                                                std::cout << "Subtracting " << hist->GetName() << " scaled by " << -dataPOT/mcPOT << " from " << sum->GetName() << "\n";
                                                sum->Add(hist, -dataPOT/mcPOT);
                                                return sum;
                                              });
          std::cout << "Test 11.4 " <<std::endl;
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
          std::cout << "Test 12 " <<std::endl;
          if(!unfoldedDaisy) throw std::runtime_error(std::string("Failed to unfold ") + DaisyFolded[petal]->GetName() + " using " + DaisyMigration[petal]->GetName());

          DaisyEffNum[petal]->Divide(DaisyEffNum[petal],DaisyEffDenom[petal]);
          unfoldedDaisy->Divide(unfoldedDaisy, DaisyEffNum[petal]);

          daisy_petal_hists[petal]=unfoldedDaisy->Clone();
          std::cout << "Test 13 " <<std::endl;
        }
        PlotUtils::MnvH1D* DaisyCorrectedC = frw->GetReweightedDaisySum(14, "carbon", daisy_petal_hists, project_dir );
        PlotUtils::MnvH1D* DaisyCorrectedFe = frw->GetReweightedDaisySum(14, "iron", daisy_petal_hists, project_dir );
        PlotUtils::MnvH1D* DaisyCorrectePb = frw->GetReweightedDaisySum(14, "lead", daisy_petal_hists, project_dir );
        //Carbon
        {
          auto fluxIntegral = frw->GetIntegratedTargetFlux(14, "carbon", DaisyCorrectedC, min_energy, max_energy, project_dir);
          auto fluxC = frw->GetTargetFluxMnvH1D(14, "carbon", project_dir);
          PlotUtils::MnvH1D* crossSectionC = normalize(DaisyCorrectedC, fluxC, nnucleons, dataPOT);
          outFileDaisy->cd();
          crossSectionC->Write((prefix+"_C_CrossSection").c_str());
        }
        //Iron
        {
          auto fluxIntegral = frw->GetIntegratedTargetFlux(14, "iron", DaisyCorrectedFe, min_energy, max_energy, project_dir);
          auto fluxC = frw->GetTargetFluxMnvH1D(14, "iron", project_dir);
          PlotUtils::MnvH1D* crossSectionC = normalize(DaisyCorrectedC, fluxC, nnucleons, dataPOT);
          outFileDaisy->cd();
          crossSectionC->Write((prefix+"_Fe_CrossSection").c_str());
        }
        //Lead
        {
          auto fluxIntegral = frw->GetIntegratedTargetFlux(14, "lead", DaisyCorrectePb, min_energy, max_energy, project_dir);
          auto fluxC = frw->GetTargetFluxMnvH1D(14, "lead", project_dir);
          PlotUtils::MnvH1D* crossSectionC = normalize(DaisyCorrectedC, fluxC, nnucleons, dataPOT);
          outFileDaisy->cd();
          crossSectionC->Write((prefix+"_Pb_CrossSection").c_str());
        }

        outFileDaisy->Close();
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
  dataFile->Close();
  mcFile->Close();
  return 0;
}
