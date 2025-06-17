//File: ExtractCrossSection.cpp
//Brief: Given data and MC files from analyses/studies/CrossSection.h, extract a 1D differential cross section.
//       Subtracts backgrounds, performs unfolding, applies efficiency x acceptance correction, and 
//       divides by flux and number of nucleons.  Writes a .root file with the cross section histogram.
//
//Usage: ExtractCrossSection <unfolding iterations> <data.root> <mc.root>
//
//Author: Andrew Olivier aolivier@ur.rochester.edu

//util includes
#include "util/GetIngredient.h"

//UnfoldUtils includes
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#include "MinervaUnfold/MnvUnfold.h"

//PlotUtils includes
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvPlotter.h"
#include "util/NukeUtils.h"
#include "PlotUtils/TargetUtils.h"
#include "PlotUtils/FluxReweighter.h"
#pragma GCC diagnostic pop

//ROOT includes
#include "TH1D.h"
#include "TFile.h"
#include "TKey.h"
#include "TParameter.h"
#include "TCanvas.h"
#include "string.h"

//Cintex is only needed for older ROOT versions like the GPVMs.
////Let CMake decide whether it's needed.
#ifndef NCINTEX
#include "Cintex/Cintex.h"
#endif

//c++ includes
#include <iostream>
#include <exception>
#include <algorithm>
#include <numeric>

//Convince the STL to talk to TIter so I can use std::find_if()
namespace std
{
  template <>
  struct iterator_traits<TIter>
  {
    using value_type = TObject;
    using pointer = TObject*;
    using reference = TObject&;
    using iterator_category = forward_iterator_tag;
  };
}


double GetTotalScatteringCenters(int targetZ, bool isMC)
{
  // TARGET INFO
  PlotUtils::TargetUtils targetInfo;
  double Nucleons;

  // Target 1 is generally excluded due to rock muon contamination (in the inclusive analysis)
  if(targetZ == 6){
    Nucleons = targetInfo.GetPassiveTargetNNucleons( 3, targetZ, isMC ); // Target 3
  }
  
  if(targetZ == 26){
    Nucleons = targetInfo.GetPassiveTargetNNucleons( 2, targetZ, isMC ) // Target 2
             + targetInfo.GetPassiveTargetNNucleons( 3, targetZ, isMC ) // Target 3
             + targetInfo.GetPassiveTargetNNucleons( 5, targetZ, isMC );// Target 5
  }

  if(targetZ == 82){
    Nucleons = targetInfo.GetPassiveTargetNNucleons( 2, targetZ, isMC ) // Target 2
             + targetInfo.GetPassiveTargetNNucleons( 3, targetZ, isMC ) // Target 3
             + targetInfo.GetPassiveTargetNNucleons( 4, targetZ, isMC ) // Target 4
             + targetInfo.GetPassiveTargetNNucleons( 5, targetZ, isMC );// Target 5
  }
  if(targetZ > 90 ){
    Nucleons = targetInfo.GetTrackerNNucleons(5980, 8422, isMC, 850);
    //double TargetUtils::GetTrackerNNucleons( double minZ, double maxZ, bool isMC, double apothem /* = 850. */ ) const

  }

  return Nucleons;     
}



//Plot a step in cross section extraction.
void Plot(PlotUtils::MnvH1D& hist, const std::string& stepName, const std::string& prefix)
{
  TCanvas can(stepName.c_str());
  hist.GetCVHistoWithError().Clone()->Draw();
  can.Print((prefix + "_" + stepName + ".png").c_str());

  //Uncertainty summary
  PlotUtils::MnvPlotter plotter;
  plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
  plotter.axis_maximum = 0.4;

  plotter.DrawErrorSummary(&hist);
  can.Print((prefix + "_" + stepName + "_uncertaintySummary.png").c_str());

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "Other");
  can.Print((prefix + "_" + stepName + "_otherUncertainties.png").c_str());
}

//Unfolding function from Aaron Bercelle
//TODO: Trim it down a little?  Remove that static?
PlotUtils::MnvH1D* UnfoldHist( PlotUtils::MnvH1D* h_folded, PlotUtils::MnvH2D* h_migration, int num_iter )
{
  static MinervaUnfold::MnvUnfold unfold;
  PlotUtils::MnvH1D* h_unfolded = nullptr;

  //bool bUnfolded = false;

  TMatrixD dummyCovMatrix;
  if(!unfold.UnfoldHisto( h_unfolded, dummyCovMatrix, h_migration, h_folded, RooUnfold::kBayes, num_iter, true, false ))
    return nullptr;

  /////////////////////////////////////////////////////////////////////////////////////////  
  //No idea if this is still needed
  //Probably.  This gets your stat unfolding covariance matrix
  TMatrixD unfoldingCovMatrixOrig; 
  int correctNbins;
  int matrixRows;  
  TH1D* hUnfoldedDummy  = new TH1D(h_unfolded->GetCVHistoWithStatError());
  TH1D* hRecoDummy      = new TH1D(h_migration->ProjectionX()->GetCVHistoWithStatError());
  TH1D* hTruthDummy     = new TH1D(h_migration->ProjectionY()->GetCVHistoWithStatError());
  TH1D* hBGSubDataDummy = new TH1D(h_folded->GetCVHistoWithStatError());
  TH2D* hMigrationDummy = new TH2D(h_migration->GetCVHistoWithStatError());
  unfold.UnfoldHisto(hUnfoldedDummy, unfoldingCovMatrixOrig, hMigrationDummy, hRecoDummy, hTruthDummy, hBGSubDataDummy,RooUnfold::kBayes, num_iter);//Stupid RooUnfold.  This is dummy, we don't need iterations

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

//The final step of cross section extraction: normalize by flux, bin width, POT, and number of targets
PlotUtils::MnvH1D* normalize(PlotUtils::MnvH1D* efficiencyCorrected, PlotUtils::MnvH1D* fluxIntegral, const double nNucleons, const double POT)
{
  efficiencyCorrected->Divide(efficiencyCorrected, fluxIntegral);

  efficiencyCorrected->Scale(1./nNucleons/POT);
  efficiencyCorrected->Scale(1.e4); //Flux histogram is in m^-2, but convention is to report cm^2
  efficiencyCorrected->Scale(1., "width");

  return efficiencyCorrected;
}

int main(const int argc, const char** argv)
{
  #ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable(); //Needed to look up dictionaries for PlotUtils classes like MnvH1D
  #endif

  TH1::AddDirectory(kFALSE); //Needed so that MnvH1D gets to clean up its own MnvLatErrorBands (which are TH1Ds).

  if(argc != 4 && argc != 5 )
  {
    std::cerr << "Expected 3 or 4 arguments, but I got " << argc-1 << ".\n"
              << "USAGE: ExtractCrossSection <unfolding iterations> <data.root> <mc.root> <OPTIONAL:numPlaylists>\n"
              << "Where <numPlaylists> is the number of playlists merged to make <data.root> and <mc.root>. If empty, it will be assumed to be 1.\n";
    return 1;
  }
  int numMergedPlaylists = 1;

  const int nIterations = std::stoi(argv[1]);
  auto dataFile = TFile::Open(argv[2], "READ");
  if(!dataFile)
  {
    std::cerr << "Failed to open data file " << argv[2] << ".\n";
    return 2;
  }

  auto mcFile = TFile::Open(argv[3], "READ");
  if(!mcFile)
  {
    std::cerr << "Failed to open MC file " << argv[3] << ".\n";
    return 3;
  }

  if (argc == 5)
  {
    numMergedPlaylists = std::stoi(argv[4]);
  }

  //std::vector<std::string> crossSectionPrefixes = {"nuke_pTmu" , "nuke_pZmu", "nuke_BjorkenX", "nuke_Erecoil", "nuke_Emu"};
  std::vector<std::string> crossSectionPrefixes = {"nuke_pTmu"};

  //std::vector<std::string> targets = { "Target7", "Target8", "Target9", "Target10", "Target11", "Target12", "1026", "1082", "2026", "2082", "3006", "3026", "3082", "4082", "5026", "5082",  "Water"}; //Is there any benefit to getting this programatically like above for the 1D prefixes?
  //std::vector<std::string> targets = {"2026", "2082", "3006", "3026", "3082", "4082", "5026", "5082"}; 
  std::vector<std::string> targets = {"2026"};

  for(auto key: *dataFile->GetListOfKeys())
  {
    const std::string keyName = key->GetName();
    if (keyName == "POTUsed") continue;
    std::cout << "keyName " << keyName <<std::endl;
    const size_t endOfPrefix = keyName.find("_data");
    std::string prefix = keyName.substr(0, endOfPrefix);
    std::cout << "prefix " << prefix <<std::endl;
    //This counts the number of '_' in the prefix, which should match the dimension, eg 2D would be nuke_pTmu_pZmu, which has 2 underscores
    bool twoDimension = (keyName == "_data_nuke_pTmu_pZmu");

    bool alreadyInVector = std::find(crossSectionPrefixes.begin(), crossSectionPrefixes.end(), prefix) != crossSectionPrefixes.end();
    std::cout << "twoDimension " << twoDimension <<std::endl;
    std::cout << "alreadyInVector " << alreadyInVector <<std::endl;
    if(endOfPrefix != std::string::npos && !alreadyInVector && !twoDimension) crossSectionPrefixes.push_back(prefix);
  }
  const double mcPOT = util::GetIngredient<TParameter<double>>(*mcFile, "POTUsed")->GetVal(),
               dataPOT = util::GetIngredient<TParameter<double>>(*dataFile, "POTUsed")->GetVal();
  std::cout<<"Targets DataPOT: " << dataPOT << " mcPOT " << mcPOT << std::endl;
  for(const auto& prefix: crossSectionPrefixes)
  {
    for (std::string& tgt : targets)
    {
      std::cout<< "Current working on prefix: " << prefix << std::endl;
      try
      {
        std::string temptgtstr = tgt;
        if (tgt == "Water") temptgtstr = "TargetWater"; //Very hacky, should just change the leaf name in the event loop so the water branch naming schema is consistent with everything else
        auto flux = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (temptgtstr+std::string("_reweightedflux_integrated")), prefix);
        flux->Scale(1.0/numMergedPlaylists);
        auto folded = util::GetIngredient<PlotUtils::MnvH1D>(*dataFile, (std::string("by_TargetCode_Data_")+tgt), prefix);
        Plot(*folded, "data", prefix);
        auto migration = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, (std::string("migration_")+tgt), prefix);
        auto effNum = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_numerator_")+tgt), prefix);
        auto effDenom = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("efficiency_denominator_")+tgt), prefix);


        auto nNucleons = util::GetIngredient<TParameter<double>>(*mcFile, (temptgtstr+std::string("_fiducial_nucleons")), prefix); //Dan: Use the same truth fiducial volume for all extractions.  The acceptance correction corrects data back to this fiducial even if the reco fiducial cut is different.
        double nNucleonsVal = nNucleons->GetVal();
        std::vector<PlotUtils::MnvH1D*> backgrounds;
        backgrounds.push_back(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, tgt+std::string("_Wrong_Material_Bkg"), prefix));
        backgrounds.push_back(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, tgt+std::string("_Wrong_Sign_Bkg"), prefix));
        backgrounds.push_back(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, tgt+std::string("_NC_Bkg"), prefix));

        for(auto key: *mcFile->GetListOfKeys())
        {
          if(std::string(key->GetName()).find(prefix + "_background_") != std::string::npos)
          {
            backgrounds.push_back(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, key->GetName()));
          }
        }

        auto simEventRate = effDenom->Clone(); //Make a copy for later
        //There are no error bands in the data, but I need somewhere to put error bands on the results I derive from it.
        folded->AddMissingErrorBandsAndFillWithCV(*migration);

        //Basing my unfolding procedure for a differential cross section on Alex's MINERvA 101 talk at https://minerva-docdb.fnal.gov/cgi-bin/private/RetrieveFile?docid=27438&filename=whatsACrossSection.pdf&version=1

        //TODO: Remove these debugging plots when done
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
                                              std::cout << "Subtracting " << hist->GetName() << " scaled by " << -dataPOT/mcPOT << " from " << sum->GetName() << "\n";
                                              sum->Add(hist, -dataPOT/mcPOT);
                                              return sum;
                                            });
        Plot(*bkgSubtracted, "backgroundSubtracted", prefix);

        auto outFile = TFile::Open((tgt+prefix + "_crossSection.root").c_str(), "RECREATE");
        if(!outFile)
        {
          std::cerr << "Could not create a file called " << prefix + "_crossSection.root" << ".  Does it already exist?\n";
          return 5;
        }

        bkgSubtracted->Write("backgroundSubtracted");

        //d'Aogstini unfolding
        auto unfolded = UnfoldHist(bkgSubtracted, migration, nIterations);
        if(!unfolded) throw std::runtime_error(std::string("Failed to unfold ") + folded->GetName() + " using " + migration->GetName());
        Plot(*unfolded, "unfolded", prefix);
        unfolded->Clone()->Write("unfolded"); //TODO: Seg fault first appears when I uncomment this line
        std::cout << "Survived writing the unfolded histogram.\n" << std::flush; //This is evidence that the problem is on the final file Write() and not unfolded->Clone()->Write().

        effNum->Divide(effNum, effDenom); //Only the 2 parameter version of MnvH1D::Divide()
                                          //handles systematics correctly.
        Plot(*effNum, "efficiency", prefix);

        unfolded->Divide(unfolded, effNum);
        Plot(*unfolded, "efficiencyCorrected", prefix);

        double nnucleons = nNucleons->GetVal()/numMergedPlaylists;
        std::cout<<"nnucleons: " << nnucleons <<std::endl;
        std::cout<<"tgt: " << tgt <<std::endl;
        //if (tgt == "1026") nnucleons = PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(1, 26, true, 850);
        //else if (tgt == "1082") nnucleons = PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(1, 82, true, 850);
        if (tgt == "2026") nnucleons = PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(2, 26, true);
        else if (tgt == "2082") nnucleons = PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(2, 82, true);
        else if (tgt == "3006") nnucleons = PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(3, 6, true);
        else if (tgt == "3026") nnucleons = PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(3, 26, true);
        else if (tgt == "3082") nnucleons = PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(3, 82, true);
        else if (tgt == "4082") nnucleons = PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(4, 82, true);
        else if (tgt == "5026") nnucleons = PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(5, 26, true);
        else if (tgt == "5082") nnucleons = PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(5, 82, true);
        else if (tgt == "Water") nnucleons = PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(6, 1, true);
        std::cout<<"nnucleons: " << nnucleons <<std::endl;




        int n_flux_universes = 100;
        int nu_pdg = 14;
        const bool use_nue_constraint = true;
        const std::string project_dir = "targets_2345_jointNueIMD";
        double min_energy = 0;
        double max_energy = 100;
        //PlotUtils::MnvH1D *fluxIntegral;
        std::string material;

        PlotUtils::MnvH1D* flux2;
        PlotUtils::MnvH1D* fluxIntegral;
        PlotUtils::MnvH1D* fluxRebinned;
        std::cout<< "Test123.1: " <<std::endl;

        //FluxReweighter* frw = new FluxReweighter( 14, use_nue_constraint, "minervame1A", FluxReweighter::gen2thin, FluxReweighter::g4numiv6, 100);

        auto& frw = PlotUtils::flux_reweighter("minervame1A", nu_pdg, use_nue_constraint, n_flux_universes);
        auto& frw2 = PlotUtils::flux_reweighter("minervame1A", nu_pdg, use_nue_constraint, n_flux_universes);
        std::cout<< "Test123.2: " <<std::endl;  
        if(tgt == "3006") material = "carbon";
        else if(tgt == "1026" || tgt == "2026" ||tgt == "3026" || tgt == "5026" ) material = "iron";
        else if(tgt == "1082" || tgt == "2082" ||tgt == "3082" || tgt == "4082" || tgt == "5082" ) material = "lead";
        else material = "tracker";
        std::cout<< "Test123.3: " <<std::endl;
        //auto &frw = PlotUtils::flux_reweighter("minervame6A", nu_pdg, use_nue_constraint, n_flux_universes);
        
        //MnvH1D* Integrated_flux = fluxReweighter->GetIntegratedFluxReweighted_FromInputFlux(flux, xsecHist, Emin, Emax);
        //xsec.integratedFlux = Integrated_flux->Clone("IntegratedFlux");


        fluxIntegral = frw.GetIntegratedTargetFlux(nu_pdg, material, unfolded, min_energy, max_energy, project_dir);
        //fluxIntegral = frw.GetIntegratedFluxReweighted(14, simEventRate, 0, 100, false);
        std::cout<< "Test123.4: " <<std::endl;
        flux2 = frw2.GetTargetFluxMnvH1D(nu_pdg, material, project_dir);

        fluxRebinned = frw2.GetRebinnedFluxReweighted_FromInputFlux(flux2, simEventRate); // issue here
        std::cout << "Using target generated flux!" << std::endl;
  







        outFile->cd();
        auto crossSection = normalize(unfolded, fluxIntegral, nnucleons, dataPOT);
                std::cout<< "Test123.5: " <<std::endl;

        Plot(*crossSection, "crossSection", prefix);
                std::cout<< "Test123.6: " <<std::endl;

        crossSection->Clone()->Write("crossSection");
                std::cout<< "Test123.7: " <<std::endl;

        simEventRate->Write("simulatedEventRate");
                std::cout<< "Test123.8: " <<std::endl;

        flux->Write("flux");
        std::cout<< "Test123.9: " <<std::endl;
        flux2->Write("flux2");
        fluxRebinned->Write("fluxRebinned");
        fluxIntegral->Write("fluxIntegral");
        //Write a "simulated cross section" to compare to the data I just extracted.
        //If this analysis passed its closure test, this should be the same cross section as
        //what GENIEXSecExtract would produce.

        int tgtCode = std::stoi(tgt);
        int targetZ = tgtCode % 1000;
        int targetID = (tgtCode - targetZ) / 1000;
        double passiveNucleons = 0;
        if (targetID<7) passiveNucleons = PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(targetID, targetZ, true); 
        else
        {
          if (targetID==7) passiveNucleons = PlotUtils::TargetUtils::Get().GetTrackerNNucleons( 7, true); 
          if (targetID==8) passiveNucleons = PlotUtils::TargetUtils::Get().GetTrackerNNucleons( 6, true); 
          if (targetID==9) passiveNucleons = PlotUtils::TargetUtils::Get().GetTrackerNNucleons( 6, true); 
          if (targetID==10) passiveNucleons = PlotUtils::TargetUtils::Get().GetTrackerNNucleons( 6, true); 
          if (targetID==11) passiveNucleons = PlotUtils::TargetUtils::Get().GetTrackerNNucleons( 6, true); 
          if (targetID==12) passiveNucleons = PlotUtils::TargetUtils::Get().GetTrackerNNucleons( 2, true); 

        }

        /* if (tgt == "1026") nnucleons =  GetTotalScatteringCenters(26, true);
        else if (tgt == "1082") nnucleons =  GetTotalScatteringCenters(82, true);
        else if (tgt == "2026") nnucleons =  GetTotalScatteringCenters(26, true);
        else if (tgt == "2082") nnucleons =  GetTotalScatteringCenters(82, true);
        else if (tgt == "3006") nnucleons =  GetTotalScatteringCenters(6, true);
        else if (tgt == "3026") nnucleons =  GetTotalScatteringCenters(26, true);
        else if (tgt == "3082") nnucleons =  GetTotalScatteringCenters(82, true);
        else if (tgt == "4082") nnucleons =  GetTotalScatteringCenters(82, true);
        else if (tgt == "5026") nnucleons =  GetTotalScatteringCenters(26, true);
        else if (tgt == "5082") nnucleons =  GetTotalScatteringCenters(82, true); */
        std::cout<<"nnucleons: " << nnucleons <<std::endl;
        std::cout<<"passiveNucleons: " << passiveNucleons <<std::endl;
        //fluxIntegral->Scale(1.103);
        auto crossSection2 = normalize(simEventRate, flux, passiveNucleons, mcPOT);

        ///!Temporary to undo the issue induced in runXseclooper where i erroneusly include the plastic before tgt2 in norm value calculations

        //double trackerAtomsC = PlotUtils::TargetUtils::Get().GetTrackerElementNAtoms( 6, 108, true );
        //double trackerAtomsC = TargetUtils::Get().GetTrackerElementNAtoms( 6, 4358, PlotUtils::TargetProp::NukeRegion::Back, true, 850.0);
        //double trackerAtomsC = TargetUtils::Get().GetTrackerElementNAtoms( 6, 5970, 8450, true, 850.0);
        //double trackerAtomsCInclTgt1 = TargetUtils::Get().GetTrackerElementNAtoms( 6, PlotUtils::TargetProp::NukeRegion::Face, PlotUtils::TargetProp::NukeRegion::Back, true, 850.0);
        
        //std::cout<< "trackerAtomsC: " << trackerAtomsC <<std::endl;

        //std::cout<< "trackerAtomsCInclTgt1: " << trackerAtomsCInclTgt1 <<std::endl;
        //double scale =  (trackerAtomsC / trackerAtomsCInclTgt1);
        //std::cout<<"scale: " <<scale<<std::endl;

        //crossSection2->Scale(scale);
        
        Plot(*crossSection2, "simulatedCrossSection", prefix);
        crossSection2->Write("simulatedCrossSection");
        std::cout<< "Test123.10: " <<std::endl;
        outFile->Close();
      }
      catch(const std::runtime_error& e)
      {
        std::cerr << "Failed to extract a cross section for prefix " << prefix << ": " << e.what() << "\n";
        return 4;
        //break;
      }
    }
  }
  dataFile->Close();
  mcFile->Close();
  return 0;
}
