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
void Plot(PlotUtils::MnvH2D &hist, const std::string &stepName, const std::string &prefix)
{
 /*  bool plotpngs = false;
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
  if (plotpngs) can.Print((prefix + "_" + stepName + "_otherUncertainties.png").c_str()); */
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
  std::string migrationpath = directory + "/runEventLoopTracker2DMigration.root";

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
  std::vector<std::string> crossSectionPrefixes = {"tracker_pTmu_pZmu"};

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

      PlotUtils::MnvH2D *DaisyEffNum[12], *DaisyFolded[12];
      PlotUtils::MnvH2D *DaisyMigration[12], *DaisyReco[12], *DaisyTruth[12];;
      std::vector<PlotUtils::MnvH2D*> DaisyBackgrounds[12];
      PlotUtils::MnvH2D *DaisyEffDenom[12], *DaisyEffDenom2P2H[12], *DaisyEffDenomDIS[12], *DaisyEffDenomRES[12], *DaisyEffDenomQE[12], *DaisyEffDenomOther[12];

      //auto flux = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("reweightedflux_integrated")), prefix);
      //flux->Scale(1.0);data
      auto folded = util::GetIngredient<PlotUtils::MnvH2D>(*dataFile, "_data_tracker_pTmu_pZmu");
      //Plot(*folded, "data", prefix);
      auto migration = util::GetIngredient<PlotUtils::MnvH2D>(*migrationFile, (std::string("migration;1")), prefix);
      auto reco = util::GetIngredient<PlotUtils::MnvH2D>(*migrationFile, (std::string("migration;2")), prefix);
      auto truth = util::GetIngredient<PlotUtils::MnvH2D>(*migrationFile, (std::string("migration;3")), prefix);
      auto effNum = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, "tracker_pZmu_tracker_pTmu_efficiency_numerator");
      auto effDenom = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, "tracker_pZmu_tracker_pTmu_efficiency_denominator");
      auto effDenom2P2H = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_2p2h")), prefix);
      auto effDenomDIS = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_DIS")), prefix);
      auto effDenomRES = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_RES")), prefix);
      auto effDenomQE = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_QE")), prefix);
      auto effDenomOther = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator_intChannels_Other")), prefix);

      //auto nNucleons = util::GetIngredient<TParameter<double>>(*mcFile, (std::string("fiducial_nucleons")), prefix); // Dan: Use the same truth fiducial volume for all extractions.  The acceptance correction corrects data back to this fiducial even if the reco fiducial cut is different.
      //double nNucleonsVal = nNucleons->GetVal();
      std::vector<PlotUtils::MnvH2D *> backgrounds;
      for (auto key : *mcFile->GetListOfKeys())
      {
        if (std::string(key->GetName()).find(prefix + "_by_BKG_") != std::string::npos)
        {
          std::cout << "Found and adding background: " << key->GetName() << std::endl;
          backgrounds.push_back(util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, key->GetName()));
        }
      }

      //***********************************************
      //Daisy reweight
      //Getting ingredients
      //***********************************************
      bool doDaisy = true;
      for (int petal=0; petal<12; petal++){
        std::string datapathdaisy = directory + "/runEventLoopTrackerData_petal_"+std::to_string(petal)+".root";
        std::string mcpathdaisy = directory + "/runEventLoopTrackerMC_petal_"+std::to_string(petal)+".root";
        bool petalFilesExist = std::filesystem::exists(datapathdaisy) && std::filesystem::exists(mcpathdaisy);
        doDaisy = doDaisy && petalFilesExist;
      }
      if (doDaisy) std::cout<<"Will do daisy petal reweight in this analysis\n";
      else std::cout<<"Will not do daisy petal reweight in this analysis\n";
      for (int petal=0; petal<12 && doDaisy; petal++){
        std::string datapathdaisy = directory + "/runEventLoopTrackerData_petal_"+std::to_string(petal)+".root";
        std::string mcpathdaisy = directory + "/runEventLoopTrackerMC_petal_"+std::to_string(petal)+".root";
        std::string migrationpathdaisy = directory + "/runEventLoopTracker2DMigration_petal_"+std::to_string(petal)+".root";

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

        auto migDaisyFile = TFile::Open(migrationpathdaisy.c_str(), "READ");
        if (!migDaisyFile)
        {
          std::cerr << "Failed to open MC file " << migrationpathdaisy.c_str() << ".\n";
          return 3;
        }


        DaisyEffNum[petal] = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, "tracker_pZmu_tracker_pTmu_efficiency_numerator");
        DaisyEffDenom[petal] = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, "tracker_pZmu_tracker_pTmu_efficiency_denominator");
        DaisyMigration[petal] = util::GetIngredient<PlotUtils::MnvH2D>(*migrationFile, (std::string("migration;1")), prefix);
        DaisyReco[petal] = util::GetIngredient<PlotUtils::MnvH2D>(*migrationFile, (std::string("migration;2")), prefix);
        DaisyTruth[petal] = util::GetIngredient<PlotUtils::MnvH2D>(*migrationFile, (std::string("migration;3")), prefix);

        DaisyEffDenom2P2H[petal] = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, "tracker_pTmu_pZmu_efficiency_denominator_intChannels_2p2h");
        DaisyEffDenomDIS[petal] = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, "tracker_pTmu_pZmu_efficiency_denominator_intChannels_DIS");
        DaisyEffDenomRES[petal] = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, "tracker_pTmu_pZmu_efficiency_denominator_intChannels_RES");
        DaisyEffDenomQE[petal] = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, "tracker_pTmu_pZmu_efficiency_denominator_intChannels_QE");
        DaisyEffDenomOther[petal] = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, "tracker_pTmu_pZmu_efficiency_denominator_intChannels_Other");

        DaisyFolded[petal] = util::GetIngredient<PlotUtils::MnvH2D>(*dataDaisyFile, "_data_tracker_pTmu_pZmu");
        auto WSBkg = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, (std::string("by_BKG_Label_Bkg_Wrong_Sign")), prefix);
        auto NCBkg =  util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, (std::string("by_BKG_Label_NC_Bkg")), prefix);
        auto OtherBkg = util::GetIngredient<PlotUtils::MnvH2D>(*mcDaisyFile, (std::string("by_BKG_Label_Other")), prefix);
        DaisyBackgrounds[petal] = {WSBkg, NCBkg, OtherBkg};
        //Backgrounds
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
      auto unfolded = UnfoldHist(bkgSubtracted, migration, reco, truth, nIterations);
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
      double nnucleonsData = targetInfo.GetTrackerNNucleons(PlotUtils::TargetProp::Tracker::Face, PlotUtils::TargetProp::Tracker::Back, false);

      PlotUtils::MnvH2D *flux2;
      PlotUtils::MnvH2D *fluxIntegral;
      PlotUtils::MnvH2D *fluxRebinned;

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
      PlotUtils::MnvH2D *fluxIntReweighted = frw->GetIntegratedFluxReweighted(pdg, simEventRate, min_energy, max_energy, true);
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
        std::map<int, PlotUtils::MnvH2D*> daisy_petal_hists;
        std::map<int, PlotUtils::MnvH2D*> daisy_petal_hists_eff_denom;
        std::map<int, PlotUtils::MnvH2D*> daisy_petal_hists_eff_denom_2p2h;
        std::map<int, PlotUtils::MnvH2D*> daisy_petal_hists_eff_denom_dis;
        std::map<int, PlotUtils::MnvH2D*> daisy_petal_hists_eff_denom_res;
        std::map<int, PlotUtils::MnvH2D*> daisy_petal_hists_eff_denom_qe;
        std::map<int, PlotUtils::MnvH2D*> daisy_petal_hists_eff_denom_other;


        for (int petal=0; petal<12; petal++){
          std::cout << "Test 11 " <<std::endl;
          std::vector<PlotUtils::MnvH2D*> DaisyBackgroundsTemp = DaisyBackgrounds[petal];
          std::cout << "Petal " <<petal <<std::endl;
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
          auto unfoldedDaisy = UnfoldHist(bkgSubtractedDaisy, DaisyMigration[petal], DaisyReco[petal], DaisyTruth[petal], nIterations);
          outFileDaisy->cd();
          unfoldedDaisy->Write((prefix+"_unfoldedDaisy_"+petal));
          std::cout << "Test 12 " <<std::endl;
          if(!unfoldedDaisy) throw std::runtime_error(std::string("Failed to unfold ") + DaisyFolded[petal]->GetName() + " using " + DaisyMigration[petal]->GetName());

          DaisyEffNum[petal]->Divide(DaisyEffNum[petal],DaisyEffDenom[petal]);
          unfoldedDaisy->Divide(unfoldedDaisy, DaisyEffNum[petal]);

          daisy_petal_hists[petal]=unfoldedDaisy->Clone();


          //1!!!!!!!!!!!!

          daisy_petal_hists_eff_denom[petal]= DaisyEffDenom[petal];
          daisy_petal_hists_eff_denom_2p2h[petal]= DaisyEffDenom2P2H[petal];
          daisy_petal_hists_eff_denom_dis[petal]= DaisyEffDenomDIS[petal];
          daisy_petal_hists_eff_denom_res[petal]= DaisyEffDenomRES[petal];
          daisy_petal_hists_eff_denom_qe[petal]= DaisyEffDenomQE[petal];
          daisy_petal_hists_eff_denom_other[petal]= DaisyEffDenomOther[petal];

          std::cout << "Test 13 " <<std::endl;
        }
        
        PlotUtils::MnvH2D* DaisyCorrectedC = frw->GetReweightedDaisySum(14, "carbon", daisy_petal_hists, project_dir );
        PlotUtils::MnvH2D* DaisyCorrectedFe = frw->GetReweightedDaisySum(14, "iron", daisy_petal_hists, project_dir );
        PlotUtils::MnvH2D* DaisyCorrectedPb = frw->GetReweightedDaisySum(14, "lead", daisy_petal_hists, project_dir );

        PlotUtils::MnvH2D* DaisyCorrectedSimC = frw->GetReweightedDaisySum(14, "carbon", daisy_petal_hists_eff_denom, project_dir );
        PlotUtils::MnvH2D* DaisyCorrectedSimFe = frw->GetReweightedDaisySum(14, "iron", daisy_petal_hists_eff_denom, project_dir );
        PlotUtils::MnvH2D* DaisyCorrectedSimPb = frw->GetReweightedDaisySum(14, "lead", daisy_petal_hists_eff_denom, project_dir );
        PlotUtils::MnvH2D* DaisyCorrectedSim2p2hC = frw->GetReweightedDaisySum(14, "carbon", daisy_petal_hists_eff_denom_2p2h, project_dir );
        PlotUtils::MnvH2D* DaisyCorrectedSim2p2hFe = frw->GetReweightedDaisySum(14, "iron", daisy_petal_hists_eff_denom_2p2h, project_dir );
        PlotUtils::MnvH2D* DaisyCorrectedSim2p2hPb = frw->GetReweightedDaisySum(14, "lead", daisy_petal_hists_eff_denom_2p2h, project_dir );
        PlotUtils::MnvH2D* DaisyCorrectedSimDISC = frw->GetReweightedDaisySum(14, "carbon", daisy_petal_hists_eff_denom_dis, project_dir );
        PlotUtils::MnvH2D* DaisyCorrectedSimDISFe = frw->GetReweightedDaisySum(14, "iron", daisy_petal_hists_eff_denom_dis, project_dir );
        PlotUtils::MnvH2D* DaisyCorrectedSimDISPb = frw->GetReweightedDaisySum(14, "lead", daisy_petal_hists_eff_denom_dis, project_dir );
        PlotUtils::MnvH2D* DaisyCorrectedSimRESC = frw->GetReweightedDaisySum(14, "carbon", daisy_petal_hists_eff_denom_res, project_dir );
        PlotUtils::MnvH2D* DaisyCorrectedSimRESFe = frw->GetReweightedDaisySum(14, "iron", daisy_petal_hists_eff_denom_res, project_dir );
        PlotUtils::MnvH2D* DaisyCorrectedSimRESPb = frw->GetReweightedDaisySum(14, "lead", daisy_petal_hists_eff_denom_res, project_dir );
        PlotUtils::MnvH2D* DaisyCorrectedSimQE = frw->GetReweightedDaisySum(14, "carbon", daisy_petal_hists_eff_denom_qe, project_dir );
        PlotUtils::MnvH2D* DaisyCorrectedSimQEFe = frw->GetReweightedDaisySum(14, "iron", daisy_petal_hists_eff_denom_qe, project_dir );
        PlotUtils::MnvH2D* DaisyCorrectedSimQEPb = frw->GetReweightedDaisySum(14, "lead", daisy_petal_hists_eff_denom_qe, project_dir );
        PlotUtils::MnvH2D* DaisyCorrectedSimOtherC = frw->GetReweightedDaisySum(14, "carbon", daisy_petal_hists_eff_denom_other, project_dir );
        PlotUtils::MnvH2D* DaisyCorrectedSimOtherFe = frw->GetReweightedDaisySum(14, "iron", daisy_petal_hists_eff_denom_other, project_dir );
        PlotUtils::MnvH2D* DaisyCorrectedSimOtherPb = frw->GetReweightedDaisySum(14, "lead", daisy_petal_hists_eff_denom_other, project_dir );

        outFileDaisy->cd();
        DaisyCorrectedSimC->Write("DaisyCorrectedSimC");
        PlotUtils::MnvH2D *crossSectionC, *crossSectionFe, *crossSectionPb;
        outFileDaisy->cd();

        PlotUtils::MnvH2D *crossSectionSimC, *crossSectionSimFe, *crossSectionSimPb;
        PlotUtils::MnvH2D *crossSectionSimC2p2h, *crossSectionSimFe2p2h, *crossSectionSimPb2p2h;
        PlotUtils::MnvH2D *crossSectionSimCDIS, *crossSectionSimFeDIS, *crossSectionSimPbDIS;
        PlotUtils::MnvH2D *crossSectionSimCRes, *crossSectionSimFeRes, *crossSectionSimPbRes;
        PlotUtils::MnvH2D *crossSectionSimCQE, *crossSectionSimFeQE, *crossSectionSimPbQE;
        PlotUtils::MnvH2D *crossSectionSimCOther, *crossSectionSimFeOther, *crossSectionSimPbOther;
        //Carbon
        {
          outFileDaisy->cd();
          DaisyCorrectedC->Write((prefix+"DaisyCorrectedC").c_str());
          auto fluxIntegral = frw->GetIntegratedTargetFlux(14, "carbon", DaisyCorrectedC, min_energy, max_energy, project_dir);
          auto fluxC = frw->GetTargetFluxMnvH1D(14, "carbon", project_dir);
          crossSectionC = normalize(DaisyCorrectedC, fluxIntegral, nnucleonsData, dataPOT);
          outFileDaisy->cd();
          crossSectionC->Write((prefix+"_C_CrossSection").c_str());
          crossSectionSimC = normalize(DaisyCorrectedSimC, fluxIntegral, nnucleons, mcPOT);
          outFileDaisy->cd();
          crossSectionSimC->Write((prefix+"_C_CrossSectionSimulated").c_str());
          crossSectionSimC2p2h = normalize(DaisyCorrectedSim2p2hC, fluxIntegral, nnucleons, mcPOT);
          outFileDaisy->cd();
          crossSectionSimC2p2h->Write((prefix+"_C_CrossSectionSimulated2p2h").c_str());
          crossSectionSimCDIS = normalize(DaisyCorrectedSimDISC, fluxIntegral, nnucleons, mcPOT);
          outFileDaisy->cd();
          crossSectionSimCDIS->Write((prefix+"_C_CrossSectionSimulatedDIS").c_str());
          crossSectionSimCRes = normalize(DaisyCorrectedSimRESC, fluxIntegral, nnucleons, mcPOT);
          outFileDaisy->cd();
          crossSectionSimCRes->Write((prefix+"_C_CrossSectionSimulatedRES").c_str());
          crossSectionSimCQE = normalize(DaisyCorrectedSimQE, fluxIntegral, nnucleons, mcPOT);
          outFileDaisy->cd();
          crossSectionSimCQE->Write((prefix+"_C_CrossSectionSimulatedQE").c_str());
          crossSectionSimCOther= normalize(DaisyCorrectedSimOtherC, fluxIntegral, nnucleons, mcPOT);
          outFileDaisy->cd();
          crossSectionSimCOther->Write((prefix+"_C_CrossSectionSimulatedOther").c_str());
        }
        //Iron
        {
          outFileDaisy->cd();
          DaisyCorrectedFe->Write((prefix+"DaisyCorrectedFe").c_str());
          auto fluxIntegral = frw->GetIntegratedTargetFlux(14, "iron", DaisyCorrectedFe, min_energy, max_energy, project_dir);
          auto fluxFe = frw->GetTargetFluxMnvH1D(14, "iron", project_dir);
          //auto fluxRebinned = frw->GetRebinnedFluxReweighted_FromInputFlux(fluxFe, DaisyCorrectedFe);
          crossSectionFe = normalize(DaisyCorrectedFe, fluxIntegral, nnucleonsData, dataPOT);
          outFileDaisy->cd();
          crossSectionFe->Write((prefix+"_Fe_CrossSection").c_str());
          outFileDaisy->cd();
          crossSectionSimFe = normalize(DaisyCorrectedSimFe, fluxIntegral, nnucleons, mcPOT);
          outFileDaisy->cd();
          crossSectionSimFe->Write((prefix+"_Fe_CrossSectionSimulated").c_str());
          crossSectionSimFe2p2h = normalize(DaisyCorrectedSim2p2hFe, fluxIntegral, nnucleons, mcPOT);
          outFileDaisy->cd();
          crossSectionSimFe2p2h->Write((prefix+"_Fe_CrossSectionSimulated2p2h").c_str());
          crossSectionSimFeDIS = normalize(DaisyCorrectedSimDISFe, fluxIntegral, nnucleons, mcPOT);
          outFileDaisy->cd();
          crossSectionSimFeDIS->Write((prefix+"_Fe_CrossSectionSimulatedDIS").c_str());
          crossSectionSimFeRes = normalize(DaisyCorrectedSimRESFe, fluxIntegral, nnucleons, mcPOT);
          outFileDaisy->cd();
          crossSectionSimFeRes->Write((prefix+"_Fe_CrossSectionSimulatedRES").c_str());
          crossSectionSimFeQE = normalize(DaisyCorrectedSimQEFe, fluxIntegral, nnucleons, mcPOT);
          outFileDaisy->cd();
          crossSectionSimFeQE->Write((prefix+"_Fe_CrossSectionSimulatedQE").c_str());
          crossSectionSimFeOther= normalize(DaisyCorrectedSimOtherFe, fluxIntegral, nnucleons, mcPOT);
          outFileDaisy->cd();
          crossSectionSimFeOther->Write((prefix+"_Fe_CrossSectionSimulatedOther").c_str());
        }
        //Lead
        {
          outFileDaisy->cd();
          DaisyCorrectedPb->Write((prefix+"DaisyCorrectedPb").c_str());
          auto fluxIntegral = frw->GetIntegratedTargetFlux(14, "lead", DaisyCorrectedPb, min_energy, max_energy, project_dir);
          auto fluxPb = frw->GetTargetFluxMnvH1D(14, "lead", project_dir);
          crossSectionPb = normalize(DaisyCorrectedPb, fluxIntegral, nnucleonsData, dataPOT);
          outFileDaisy->cd();
          crossSectionPb->Write((prefix+"_Pb_CrossSection").c_str());
          crossSectionSimPb = normalize(DaisyCorrectedSimPb, fluxIntegral, nnucleons, mcPOT);
          outFileDaisy->cd();
          crossSectionSimPb->Write((prefix+"_Pb_CrossSectionSimulated").c_str());
          crossSectionSimPb2p2h = normalize(DaisyCorrectedSim2p2hPb, fluxIntegral, nnucleons, mcPOT);
          outFileDaisy->cd();
          crossSectionSimPb2p2h->Write((prefix+"_Pb_CrossSectionSimulated2p2h").c_str());
          crossSectionSimPbDIS = normalize(DaisyCorrectedSimDISPb, fluxIntegral, nnucleons, mcPOT);
          outFileDaisy->cd();
          crossSectionSimPbDIS->Write((prefix+"_Pb_CrossSectionSimulatedDIS").c_str());
          crossSectionSimPbRes = normalize(DaisyCorrectedSimRESPb, fluxIntegral, nnucleons, mcPOT);
          outFileDaisy->cd();
          crossSectionSimPbRes->Write((prefix+"_Pb_CrossSectionSimulatedRES").c_str());
          crossSectionSimPbQE = normalize(DaisyCorrectedSimQEPb, fluxIntegral, nnucleons, mcPOT);
          outFileDaisy->cd();
          crossSectionSimPbQE->Write((prefix+"_Pb_CrossSectionSimulatedQE").c_str());
          crossSectionSimPbOther= normalize(DaisyCorrectedSimOtherPb, fluxIntegral, nnucleons, mcPOT);
          outFileDaisy->cd();
          crossSectionSimPbOther->Write((prefix+"_Pb_CrossSectionSimulatedOther").c_str());
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
