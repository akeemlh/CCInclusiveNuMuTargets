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

  if (argc != 5)
  {
    std::cerr << "Expected 4 arguments, but I got " << argc - 1 << ".\n"
              << "USAGE: ExtractCrossSection <unfolding iterations> <directory> <target> <playlistname-for-flux>\n"
              << "Where <target> is the name of a target or \"all\" to loop over all targets found\n"
              << "or the name of the material(+'_tuned', if applicable), if looking at combined playlists\n"
              << "e.g: ExtractCrossSection 5 ./ 2026 minervame1A  -- to just extract xsecs for tgt2 iron\n"
              << "e.g: ExtractCrossSection 5 ./ all minervame1A  -- to just extract xsecs for all targets\n"
              << "e.g: ExtractCrossSection 5 ./ iron minervame1A  -- to just extract xsecs for combined iron\n"
              << "e.g: ExtractCrossSection 5 ./ iron_tuned minervame1A  -- to just extract xsecs for combined iron after sideband tune\n";
    return 1;
  }

  const int nIterations = std::stoi(argv[1]);
  std::string directory = std::string(argv[2]);
  std::string intgt = std::string(argv[3]);
  std::string platlistname = std::string(argv[4]);

  std::vector<std::string> targets;
  if (intgt == "ALL" || intgt == "all")
  {
    for (const auto &entry : std::filesystem::directory_iterator(directory))
    {
      std::string path = entry.path();
      std::cout << "path " << path << std::endl;
      const size_t base = path.find("runEventLoopTargetsData");
      if (base != std::string::npos)
      {
        const size_t ext = path.find(".root");
        std::cout << "base " << base << std::endl;
        std::cout << "ext " << ext << std::endl;
        int lenTgtString = ext - base - 23;
        std::cout << "lenTgtString: " << lenTgtString << std::endl;
        std::string tgtstring = path.substr(base + 23, lenTgtString);
        std::cout << "tgtstring: " << tgtstring << std::endl;
        targets.push_back(tgtstring);
      }
    }
  }
  else targets = {intgt};
  for (std::string &tgt : targets)
  {
    //if (tgt != "3006") continue; //Used for testing with only one target
    std::string datapath = directory + "/runEventLoopTargetsData" + tgt + ".root";
    std::string mcpath = directory + "/runEventLoopTargetsMC" + tgt + ".root";
    std::string migrationpath = directory + "/runEventLoopTargets2DMigration" + tgt + ".root";

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

    std::vector<std::string> crossSectionPrefixes = {"nuke_pTmu_pZmu"};

    for (auto key : *dataFile->GetListOfKeys())
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
    }
    const double mcPOT = util::GetIngredient<TParameter<double>>(*mcFile, "POTUsed")->GetVal(),
                  dataPOT = util::GetIngredient<TParameter<double>>(*dataFile, "POTUsed")->GetVal();
    std::cout << "Targets DataPOT: " << dataPOT << " mcPOT " << mcPOT << std::endl;
    for (const auto &prefix : crossSectionPrefixes)
    {
      //if (!(prefix == "nuke_Erecoil" || prefix == "nuke_pTmu")) continue; //Used for testing with only subset of prefixes
      std::cout << "Current working on prefix: " << prefix << std::endl;
      try
      {
        //auto flux = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("reweightedflux_integrated")), prefix);
        //flux->Scale(1.0);
        auto folded = util::GetIngredient<PlotUtils::MnvH2D>(*dataFile, (std::string("data")), prefix);
        Plot(*folded, "data", prefix, tgt);
        auto migration = util::GetIngredient<PlotUtils::MnvH2D>(*migrationFile, (std::string("migration")), prefix);
        auto migration_reco = util::GetIngredient<PlotUtils::MnvH2D>(*migrationFile, (std::string("reco")), prefix);
        auto migration_truth = util::GetIngredient<PlotUtils::MnvH2D>(*migrationFile, (std::string("truth")), prefix);
        auto effNum = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_numerator")), "nuke_pZmu_nuke_pTmu");
        auto effDenom = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("efficiency_denominator")), "nuke_pZmu_nuke_pTmu");

        //auto nNucleons = util::GetIngredient<TParameter<double>>(*mcFile, (std::string("fiducial_nucleons")), prefix); // Dan: Use the same truth fiducial volume for all extractions.  The acceptance correction corrects data back to this fiducial even if the reco fiducial cut is different.
        //double nNucleonsVal = nNucleons->GetVal();
        std::vector<PlotUtils::MnvH2D *> backgrounds;
        for (auto key : *mcFile->GetListOfKeys())
        {
          if (std::string(key->GetName()).find(prefix + "_by_BKG_Label_") != std::string::npos)
          {
            std::cout << "Found and adding background: " << key->GetName() << std::endl;
            backgrounds.push_back(util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, key->GetName()));
          }
        }

        auto simEventRate = effDenom->Clone(); // Make a copy for later
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
        Plot(*toSubtract, "BackgroundSum", prefix, tgt);

        auto bkgSubtracted = std::accumulate(backgrounds.begin(), backgrounds.end(), folded->Clone(),
                                              [mcPOT, dataPOT](auto sum, const auto hist)
                                              {
                                                std::cout << "Subtracting " << hist->GetName() << " scaled by " << -dataPOT / mcPOT << " from " << sum->GetName() << "\n";
                                                sum->Add(hist, -dataPOT / mcPOT);
                                                return sum;
                                              });
        Plot(*bkgSubtracted, "backgroundSubtracted", prefix, tgt);

        auto outFile = TFile::Open((tgt + prefix + "_crossSection.root").c_str(), "RECREATE");
        if (!outFile)
        {
          std::cerr << "Could not create a file called " << prefix + "_crossSection.root" << ".  Does it already exist?\n";
          return 5;
        }

        bkgSubtracted->Write("backgroundSubtracted");

        // d'Aogstini unfolding
        auto unfolded = UnfoldHist(bkgSubtracted, migration, migration_reco, migration_truth, nIterations);
        if (!unfolded)
          throw std::runtime_error(std::string("Failed to unfold ") + folded->GetName() + " using " + migration->GetName());
        Plot(*unfolded, "unfolded", prefix, tgt);
        unfolded->Clone()->Write("unfolded"); // TODO: Seg fault first appears when I uncomment this line
        std::cout << "Survived writing the unfolded histogram.\n"
                  << std::flush; // This is evidence that the problem is on the final file Write() and not unfolded->Clone()->Write().

        effNum->Divide(effNum, effDenom); // Only the 2 parameter version of MnvH1D::Divide()
                                          // handles systematics correctly.
        Plot(*effNum, "efficiency", prefix, tgt);

        unfolded->Divide(unfolded, effNum);
        Plot(*unfolded, "efficiencyCorrected", prefix, tgt);

        double nnucleons = 0;
        double nnucleonsData = 0;


        int n_flux_universes = 100; // Is this right
        int nu_pdg = 14;
        const bool use_nue_constraint = true;
        const std::string project_dir = "targets_2345_jointNueIMD";
        double min_energy = 0;
        double max_energy = 100;



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
          if (tgtNum<7) nnucleons = targetInfo.GetPassiveTargetNNucleons(tgtNum, tgtMat, true);
          else
          {
            if (tgtNum==7) nnucleons = targetInfo.GetTrackerNNucleons( 7, true); 
            if (tgtNum==8) nnucleons = targetInfo.GetTrackerNNucleons( 6, true); 
            if (tgtNum==9) nnucleons = targetInfo.GetTrackerNNucleons( 6, true); 
            if (tgtNum==10) nnucleons = targetInfo.GetTrackerNNucleons( 6, true); 
            if (tgtNum==11) nnucleons = targetInfo.GetTrackerNNucleons( 6, true); 
            if (tgtNum==12) nnucleons = targetInfo.GetTrackerNNucleons( 2, true); 
          }
          if (tgtNum<7) nnucleonsData = targetInfo.GetPassiveTargetNNucleons(tgtNum, tgtMat, false);
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

        std::cout<<"nu_pdg: " << nu_pdg << std::endl;
        std::cout<<"use_nue_constraint: " << use_nue_constraint << std::endl;
        std::cout<<"n_flux_universes: " << n_flux_universes << std::endl;
        std::cout<<"min_energy: " << min_energy << std::endl;
        std::cout<<"max_energy: " << max_energy << std::endl;
        std::cout<<"nnucleons: " << nnucleons << std::endl;
        std::cout<<"nnucleonsData: " << nnucleonsData << std::endl;
        std::cout<<"mcPOT: " << mcPOT << std::endl;
        std::cout<<"dataPOT: " << dataPOT << std::endl;
        //PlotUtils::FluxReweighter frw = PlotUtils::flux_reweighter(platlistname, nu_pdg, use_nue_constraint, n_flux_universes);
        PlotUtils::FluxReweighter *frw = new PlotUtils::FluxReweighter( nu_pdg, use_nue_constraint, platlistname, PlotUtils::FluxReweighter::gen2thin, PlotUtils::FluxReweighter::g4numiv6, n_flux_universes );
        //std::cout << "ABC123-1 " << std::endl;
        PlotUtils::MnvH2D *fluxIntReweighted = frw->GetIntegratedFluxReweighted(nu_pdg, simEventRate, min_energy, max_energy, true);
        //fluxIntegral = frw->GetIntegratedTargetFlux(nu_pdg, material, simEventRate, min_energy, max_energy, project_dir);
        //auto &frw2 = PlotUtils::flux_reweighter("minervame1A", nu_pdg, use_nue_constraint, n_flux_universes);
        //flux2 = frw2.GetTargetFluxMnvH1D(nu_pdg, material, project_dir);

        //fluxRebinned = frw->GetRebinnedFluxReweighted_FromInputFlux(flux2, simEventRate); // issue here
        //PlotUtils::MnvH1D *Integrated_fluxGenie = frw->GetIntegratedFluxReweighted_FromInputFlux(flux2, simEventRate, min_energy, max_energy);
        //PlotUtils::MnvH1D *Integrated_fluxGenie2 = frw->GetIntegratedFluxReweighted_FromInputFlux(flux, simEventRate, min_energy, max_energy);
        outFile->cd();
        unfolded->Write("unfolded2");
        fluxIntReweighted->Write("fluxIntReweighted2");
        auto crossSection = normalize(unfolded, fluxIntReweighted, nnucleonsData, dataPOT);
        Plot(*crossSection, "crossSection", prefix, tgt);
        crossSection->Clone()->Write("crossSection");
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
        auto crossSection2 = normalize(simEventRate, fluxIntReweighted, nnucleons, mcPOT);
        Plot(*crossSection2, "simulatedCrossSection", prefix, tgt);
        crossSection2->Write("simulatedCrossSection");
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
    dataFile->Close();
    mcFile->Close();
  }
  return 0;
}
