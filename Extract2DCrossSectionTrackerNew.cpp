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
#pragma GCC diagnostic pop

//ROOT includes
#include "TH1D.h"
#include "TFile.h"
#include "TKey.h"
#include "TParameter.h"
#include "TCanvas.h"

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

  if(argc < 5 )
  {
    std::cerr << "Expected 4 or more arguments, but I got " << argc-1 << ".\n"
              << "USAGE: ExtractCrossSection <unfolding iterations> <data_directory> <mc_directory> <playlists>\n"
              << "EG: ExtractCrossSection 5 /data/files/here/ /mc/files/here/ 1A 1B 1C\n"
              << "Note 23/Jun/2025: DONT TRY TO RUN MULTIPLE PLAYLISTS TOGETHER, JUST DO ONE AT A TIME I.E ONLY PASS ONE PLAYLIST AS <playlists>\n";
              
    return 1;
  }
  int numMergedPlaylists = argc - 4;

  const int nIterations = std::stoi(argv[1]);

  std::vector<std::string> variables = {"pTmu", "pZmu", "BjorkenX", "Erecoil", "Emu"};

  double mcPOT = 0;
  double dataPOT = 0;
  //std::cout<<"Data POT: " << dataPOT << " mcPOT " << mcPOT << std::endl;
  for(const auto& var: variables)
  {
    std::string prefix = "tracker_" + var;
    try
    {
      
      PlotUtils::MnvH1D *flux, *folded, *effNum, *effDenom, *DaisyEffNum[12], *DaisyEffDenom[12];
      PlotUtils::MnvH2D* migration, *DaisyMigration[12];
      std::vector<PlotUtils::MnvH1D*> backgrounds, DaisyBackgrounds[12];
      double nNucleonsVal = 0;

      for (int n = 0; n < numMergedPlaylists; n++)
      {
        std::string dataFileName = std::string(argv[2]) + "/" + std::string(argv[4+n])+"-runEventLoopDataTracker.root";
        std::string mcFileName = std::string(argv[3]) + "/" + std::string(argv[4+n])+"-runEventLoopMCTracker.root";
        auto dataFile = TFile::Open(dataFileName.c_str(), "READ");
        if(!dataFile)
        {
          std::cerr << "Failed to open data file " << dataFileName << ".\n";
          return 2;
        }

        auto mcFile = TFile::Open(mcFileName.c_str(), "READ");
        if(!mcFile)
        {
          std::cerr << "Failed to open MC file " << mcFileName << ".\n";
          return 3;
        }

        double fileMCPOT= util::GetIngredient<TParameter<double>>(*mcFile, "POTUsed")->GetVal();
        double fileDataPOT= util::GetIngredient<TParameter<double>>(*dataFile, "POTUsed")->GetVal();
        mcPOT += fileMCPOT;
        dataPOT += fileDataPOT;
        std::cout<<"Data POT: " << dataPOT << " mcPOT " << mcPOT << std::endl;

        if ( n == 0)
        {
          flux = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "reweightedflux_integrated", prefix);
          flux->Scale(fileMCPOT);
          folded = util::GetIngredient<PlotUtils::MnvH1D>(*dataFile, "data", prefix);
          migration = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, "migration", prefix);
          effNum = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "efficiency_numerator", prefix);
          effDenom = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "efficiency_denominator", prefix);

          //Daisy reweight
          for (int petal=0; petal<12; petal++){
            DaisyEffNum[petal] = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, ("Daisy_EffNum_"+std::to_string(petal)).c_str(), prefix);
            DaisyEffDenom[petal] = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, ("Daisy_EffDenom_"+std::to_string(petal)).c_str(), prefix);
            DaisyEffDenom[petal] = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, ("Daisy_Migration_"+std::to_string(petal)).c_str(), prefix);
            //Backgrounds
            for(auto key: *mcFile->GetListOfKeys())
            {
              if(std::string(key->GetName()).find(prefix + "_Daisy_Background_"+std::to_string(petal)) != std::string::npos)
              {
                DaisyBackgrounds[petal].push_back(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, key->GetName()));
              }
            }
          }
        }
        else //Not functional yet/ever -- DONT TRY TO RUN MULTIPLE PLAYLISTS TOGETHER, JUST DO ONE AT A TIME
        {
          PlotUtils::MnvH1D* tempFlux = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "reweightedflux_integrated", prefix);
          tempFlux->Scale(fileMCPOT);
          flux->Add(tempFlux);
          folded->Add(util::GetIngredient<PlotUtils::MnvH1D>(*dataFile, "data", prefix));
          migration->Add(util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, "migration", prefix));
          effNum->Add(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "efficiency_numerator", prefix));
          effDenom->Add(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "efficiency_denominator", prefix));
        }
        auto nNucleons = util::GetIngredient<TParameter<double>>(*mcFile, ("fiducial_nucleons"), prefix); //Dan: Use the same truth fiducial volume for all extractions.  The acceptance correction corrects data back to this fiducial even if the reco fiducial cut is different.
        nNucleonsVal += nNucleons->GetVal();
        for(auto key: *mcFile->GetListOfKeys())
        {
          if(std::string(key->GetName()).find(prefix + "_background_") != std::string::npos)
          {
            backgrounds.push_back(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, key->GetName()));
          }
        }
        dataFile->Close();
        mcFile->Close();

      }
      flux->Scale(1.0/mcPOT); //POT Weighted average flux -- Is this the right way to do it?  

      std::cout<< "Current working on prefix: " << prefix << std::endl;
      //Plot(*folded, "data", prefix);
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
      //Plot(*toSubtract, "BackgroundSum", prefix);

      auto bkgSubtracted = std::accumulate(backgrounds.begin(), backgrounds.end(), folded->Clone(),
                                           [mcPOT, dataPOT](auto sum, const auto hist)
                                           {
                                             std::cout << "Subtracting " << hist->GetName() << " scaled by " << -dataPOT/mcPOT << " from " << sum->GetName() << "\n";
                                             sum->Add(hist, -dataPOT/mcPOT);
                                             return sum;
                                           });
      //Plot(*bkgSubtracted, "backgroundSubtracted", prefix);

      auto outFile = TFile::Open("ExtractedCrossSection.root", "UPDATE");
      if(!outFile)
      {
        std::cerr << "Could not create a file called ExtractedCrossSection.root\n";
        return 5;
      }

      bkgSubtracted->Write((var+"_backgroundSubtracted").c_str());

      //d'Aogstini unfolding
      auto unfolded = UnfoldHist(bkgSubtracted, migration, nIterations);
      if(!unfolded) throw std::runtime_error(std::string("Failed to unfold ") + folded->GetName() + " using " + migration->GetName());
      //Plot(*unfolded, "unfolded", prefix);
      unfolded->Clone()->Write((var+"_unfolded").c_str()); //TODO: Seg fault first appears when I uncomment this line
      std::cout << "Survived writing the unfolded histogram.\n" << std::flush; //This is evidence that the problem is on the final file Write() and not unfolded->Clone()->Write().

      effNum->Divide(effNum, effDenom); //Only the 2 parameter version of MnvH1D::Divide()
                                        //handles systematics correctly.
      //Plot(*effNum, "efficiency", prefix);

      unfolded->Divide(unfolded, effNum);
      Plot(*unfolded, "efficiencyCorrected", prefix);

      auto crossSection = normalize(unfolded, flux, nNucleonsVal/numMergedPlaylists, dataPOT);
      //Plot(*crossSection, "crossSection", prefix);
      crossSection->Clone()->Write((var+"_crossSection").c_str());

      simEventRate->Write((var+"_simulatedEventRate").c_str());
      flux->Write((var+"_flux_reweighted").c_str());
      //Write a "simulated cross section" to compare to the data I just extracted.
      //If this analysis passed its closure test, this should be the same cross section as
      //what GENIEXSecExtract would produce.
      normalize(simEventRate, flux, nNucleonsVal/numMergedPlaylists, mcPOT);
      
      //Plot(*simEventRate, "simulatedCrossSection", prefix);
      simEventRate->Write((var+"_simulatedCrossSection").c_str());
      outFile->Close();
    }
    catch(const std::runtime_error& e)
    {
      std::cerr << "Failed to extract a cross section for prefix " << prefix << ": " << e.what() << "\n";
      return 4;
    }
  }

  return 0;
}
