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
#include "TSystemDirectory.h"

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
#include <utility>

#include <util/NukeUtils.h>

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
void Plot(PlotUtils::MnvH1D& hist, const std::string& stepName, const std::string& prefix, const std::string& tgt)
{
  TCanvas can(stepName.c_str());
  hist.GetCVHistoWithError().Clone()->Draw();
  can.Print((tgt+"_"+prefix + "_" + stepName + ".png").c_str());

  //Uncertainty summary
  PlotUtils::MnvPlotter plotter;
  plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
  plotter.axis_maximum = 0.4;

  plotter.DrawErrorSummary(&hist);
  can.Print((tgt+"_"+prefix + "_" + stepName + "_uncertaintySummary.png").c_str());

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "Other");
  can.Print((tgt+"_"+prefix + "_" + stepName + "_otherUncertainties.png").c_str());
}

std::pair<double, double> getCombinedPOTs(std::string dir, std::vector<std::string> playlists)
{
  double mcPOT = 0;
  double dataPOT = 0;
  for (auto p : playlists)
  {
    std::string datafilepath = dir+p+"-runEventLoopDataTargets.root";
    std::string mcfilepath = dir+p+"-runEventLoopMCTargets.root";
    TFile *data = TFile::Open((datafilepath).c_str());
    dataPOT += util::GetIngredient<TParameter<double>>(*data, "POTUsed")->GetVal();
    std::cout << "Data POT: " << dataPOT <<std::endl;
    TFile *mc = TFile::Open((mcfilepath).c_str());
    mcPOT += util::GetIngredient<TParameter<double>>(*mc, "POTUsed")->GetVal();
    std::cout << "MC POT: " << mcPOT <<std::endl;
  }
  return std::make_pair(dataPOT, mcPOT);
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

  if (!(argc>3))
  {
    std::cerr << "Expected 4 or more arguments, but I got " << argc-1 << ".\n"
              << "USAGE: ExtractCrossSection <unfolding iterations> <root file directory> <playlists>\n"
              << "eg: ExtractCrossSection 5 /files/in/here/ 1A 1C 1D\n";
    return 1;
  }

  const int nIterations = std::stoi(argv[1]);
  std::string dir = std::string(argv[2]);
  if (dir.back()!='/') dir+="/";
  std::vector<std::string> playlists;
  for (int i = 3; i<argc; i++)
  {
    std::string datafilepath = dir+std::string(argv[i])+"-runEventLoopDataTargets.root";
    auto dataFile = TFile::Open(datafilepath.c_str(), "READ");
    if(!dataFile)
    {
      std::cerr << "Failed to open  " << datafilepath << ".\n";
      return 2;
    };
    dataFile->Close();
    std::string mcfilepath = dir+std::string(argv[i])+"-runEventLoopMCTargets.root";
    auto mcFile = TFile::Open(mcfilepath.c_str(), "READ");
    if(!mcFile)
    {
      std::cerr << "Failed to open  " << mcfilepath << ".\n";
      return 3;
    };
    mcFile->Close();
    playlists.push_back(std::string(argv[i]));
  }

  std::string dataFilePath = dir+std::string(playlists[0])+"-runEventLoopDataTargets.root";
  auto dataFile = TFile::Open(dataFilePath.c_str(), "READ");

  std::cout<<"Here-11\n";

  //Getting prefixes
  std::vector<std::string> crossSectionPrefixes;
  for(auto key: *dataFile->GetListOfKeys())
  {
    const std::string keyName = key->GetName();
    const size_t endOfPrefix = keyName.find("_by_TargetCode_Data_");
    std::string prefix = keyName.substr(0, endOfPrefix);
    //+4 below because "nuke", which prefaces all of the nuclear target stuff is 4 letters long
    //This counts the number of '_' in the prefix, which should match the dimension, eg 2D would be nuke_pTmu_pZmu, which has 2 underscores
    bool oneDimension = std::count_if( prefix.begin()+4, prefix.end(), []( char c ){return c =='_';})==1;
    bool alreadyInVector = std::find(crossSectionPrefixes.begin(), crossSectionPrefixes.end(), prefix) != crossSectionPrefixes.end();
    if(endOfPrefix != std::string::npos && !alreadyInVector && oneDimension) crossSectionPrefixes.push_back(prefix);
  }
  dataFile->Close();
  std::cout<<"Here-2\n";
  //const double mcPOT = util::GetIngredient<TParameter<double>>(*mcFile, "POTUsed")->GetVal(),
  //             dataPOT = util::GetIngredient<TParameter<double>>(*dataFile, "POTUsed")->GetVal();
  std::pair<double, double> POTs = getCombinedPOTs(dir, playlists); //Note !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! need to exclude the combined/madd-ed files
  const double dataPOT = POTs.first;
  const double mcPOT = POTs.second;
  std::cout<<"Data POT: " << dataPOT << " mcPOT " << mcPOT << std::endl;
  std::cout<<"Here-3\n";


  std::cout<<"NumTgts " << util::TgtCodeLabelsNuke.size();
  for(const auto& prefix: crossSectionPrefixes)
  {
    //std::map<int, std::string> testTgt = {{3026, "3026"}, {1026, "1026"}};
    //for(const auto& TgtCode: testTgt) 
    for(const auto& TgtCode: util::TgtCodeLabelsNuke) 
    {
      PlotUtils::MnvH1D* flux = nullptr;
      PlotUtils::MnvH1D* folded = nullptr;
      PlotUtils::MnvH2D* migration = nullptr;
      PlotUtils::MnvH1D* effNum = nullptr;
      PlotUtils::MnvH1D* effDenom = nullptr;
      PlotUtils::MnvH1D* simEventRate = nullptr;
      std::map<std::string, PlotUtils::MnvH1D*> backgrounds;
      TParameter<double>* nNucleons;

      try
      {
        for (auto playlist : playlists)
        {
          std::string dataFilePath = dir + std::string(playlist)+"-runEventLoopDataTargets.root";
          dataFile = TFile::Open(dataFilePath.c_str(), "READ");

          std::string mcFilePath = dir + std::string(playlist)+"-runEventLoopMCTargets.root";
          auto mcFile = TFile::Open(mcFilePath.c_str(), "READ");

          std::cout<<"Here1\n";
          std::cout<<"TgtCode " << TgtCode.second <<std::endl;
          std::string tgtVal;
          if (TgtCode.first>=7 && TgtCode.first <=12) tgtVal=std::to_string(TgtCode.first);
          else tgtVal=TgtCode.second;
          if (flux==nullptr) flux = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "target"+tgtVal+"_reweightedflux_integrated", prefix);
          else flux->Add(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "target"+tgtVal+"_reweightedflux_integrated", prefix));
          if (folded==nullptr) folded = util::GetIngredient<PlotUtils::MnvH1D>(*dataFile,  "by_TargetCode_Data_"+TgtCode.second, prefix);
          else folded->Add(util::GetIngredient<PlotUtils::MnvH1D>(*dataFile,  "by_TargetCode_Data_"+TgtCode.second, prefix));
          Plot(*folded, "data", prefix, TgtCode.second);
          if (migration==nullptr) migration = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, "migration_"+TgtCode.second, prefix);
          else migration->Add(util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, "migration_"+TgtCode.second, prefix));
          if (effNum==nullptr) effNum = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "efficiency_numerator_"+TgtCode.second, prefix);
          else effNum->Add(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "efficiency_numerator_"+TgtCode.second, prefix));
          if (effDenom==nullptr) effDenom = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "efficiency_denominator_"+TgtCode.second, prefix);
          else effDenom->Add(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "efficiency_denominator_"+TgtCode.second, prefix));
          simEventRate = effDenom->Clone(); //Make a copy for later
          /* const auto fiducialFound = std::find_if(mcFile->GetListOfKeys()->begin(), mcFile->GetListOfKeys()->end(),
                                                  [&prefix](const auto key)
                                                  {
                                                    const std::string keyName = key->GetName();
                                                    const size_t fiducialEnd = keyName.find("_fiducial_nucleons");
                                                    return (fiducialEnd != std::string::npos) && (prefix.find(keyName.substr(0, fiducialEnd)) != std::string::npos);
                                                  });
          if(fiducialFound == mcFile->GetListOfKeys()->end()) throw std::runtime_error("Failed to find a number of nucleons that matches prefix " + prefix);
          */
          //auto nNucleons = util::GetIngredient<TParameter<double>>(*mcFile, (*fiducialFound)->GetName()); //Dan: Use the same truth fiducial volume for all extractions.  The acceptance correction corrects data back to this fiducial even if the reco fiducial cut is different.
          nNucleons = util::GetIngredient<TParameter<double>>(*mcFile, "target"+tgtVal+"_fiducial_nucleons", prefix); //Dan: Use the same truth fiducial volume for all extractions.  The acceptance correction corrects data back to this fiducial even if the reco fiducial cut is different.

          //Look for backgrounds with <prefix>_<analysis>_Background_<name>
          /* for(auto key: *mcFile->GetListOfKeys())
          {
            if(std::string(key->GetName()).find(prefix + "_background_") != std::string::npos)
            {
              backgrounds.push_back(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, key->GetName()));
            }
          } */

          if (backgrounds.find("NC") == backgrounds.end()) {
            backgrounds["MC"]=util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "tgt"+TgtCode.second+"_NC_Bkg", prefix);
          } else {
            backgrounds["MC"]->Add(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "tgt"+TgtCode.second+"_NC_Bkg", prefix));
          }
          if (backgrounds.find("WrongMaterial") == backgrounds.end()) {
            backgrounds["WrongMaterial"]=util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "tgt"+TgtCode.second+"_Wrong_Material_Bkg", prefix);
          } else {
            backgrounds["WrongMaterial"]->Add(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "tgt"+TgtCode.second+"_Wrong_Material_Bkg", prefix));
          }
          dataFile->Close();
          mcFile->Close();
        }
        //There are no error bands in the data, but I need somewhere to put error bands on the results I derive from it.
        folded->AddMissingErrorBandsAndFillWithCV(*migration);

        //Basing my unfolding procedure for a differential cross section on Alex's MINERvA 101 talk at https://minerva-docdb.fnal.gov/cgi-bin/private/RetrieveFile?docid=27438&filename=whatsACrossSection.pdf&version=1

        //TODO: Remove these debugging plots when done
        auto toSubtract = std::accumulate(std::next(backgrounds.begin()), backgrounds.end(), (*backgrounds.begin()).second->Clone(),
                                          [](auto sum, const auto hist)
                                          {
                                            sum->Add(hist.second);
                                            return sum;
                                          });
        Plot(*toSubtract, "BackgroundSum", prefix, TgtCode.second);

        auto bkgSubtracted = std::accumulate(backgrounds.begin(), backgrounds.end(), folded->Clone(),
                                            [mcPOT, dataPOT](auto sum, const auto hist)
                                            {
                                              std::cout << "Subtracting " << (hist.second)->GetName() << " scaled by " << -dataPOT/mcPOT << " from " << sum->GetName() << "\n";
                                              sum->Add(hist.second, -dataPOT/mcPOT);
                                              return sum;
                                            });
        Plot(*bkgSubtracted, "backgroundSubtracted", prefix, TgtCode.second);
        auto outFile = TFile::Open((prefix + TgtCode.second+"_crossSection.root").c_str(), "RECREATE");
        if(!outFile)
        {
          std::cerr << "Could not create a file called " << prefix + "_crossSection.root" << ".  Does it already exist?\n";
          return 5;
        }

        bkgSubtracted->Write("backgroundSubtracted");

        //d'Aogstini unfolding
        auto unfolded = UnfoldHist(bkgSubtracted, migration, nIterations);
        if(!unfolded) throw std::runtime_error(std::string("Failed to unfold ") + folded->GetName() + " using " + migration->GetName());
        Plot(*unfolded, "unfolded", prefix, TgtCode.second);
        unfolded->Clone()->Write("unfolded"); //TODO: Seg fault first appears when I uncomment this line
        std::cout << "Survived writing the unfolded histogram.\n" << std::flush; //This is evidence that the problem is on the final file Write() and not unfolded->Clone()->Write().

        effNum->Divide(effNum, effDenom); //Only the 2 parameter version of MnvH1D::Divide()
                                          //handles systematics correctly.
        Plot(*effNum, "efficiency", prefix, TgtCode.second);

        unfolded->Divide(unfolded, effNum);
        Plot(*unfolded, "efficiencyCorrected", prefix, TgtCode.second);

        auto crossSection = normalize(unfolded, flux, nNucleons->GetVal(), dataPOT);
        Plot(*crossSection, "crossSection", prefix, TgtCode.second);
        crossSection->Clone()->Write("crossSection");

        simEventRate->Write("simulatedEventRate");
        flux->Write("flux");
        //Write a "simulated cross section" to compare to the data I just extracted.
        //If this analysis passed its closure test, this should be the same cross section as
        //what GENIEXSecExtract would produce.
        normalize(simEventRate, flux, nNucleons->GetVal(), mcPOT);
        
        Plot(*simEventRate, "simulatedCrossSection", prefix, TgtCode.second);
        simEventRate->Write("simulatedCrossSection");

        outFile->Close();
      }
      catch(const std::runtime_error& e)
      {
        std::cerr << "Failed to extract a cross section for prefix " << prefix << ": " << e.what() << "\n";
        return 4;
      }
    }
  }

  return 0;
}
