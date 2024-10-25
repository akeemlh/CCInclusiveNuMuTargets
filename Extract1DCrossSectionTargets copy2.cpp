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

std::pair<double, double> getCombinedPOTs(std::string datadir, std::string mcdir)
{
  double mcPOT = 0;
  double dataPOT = 0;
  TSystemDirectory dir_data(datadir.c_str(), datadir.c_str());
  TList *data_files = dir_data.GetListOfFiles();
  if (data_files) {
      TSystemFile *file;
      TString fname;
      TIter next(data_files);
      while ((file=(TSystemFile*)next())) {
          fname = file->GetName();
          if (!file->IsDirectory()) {
              std::cout << (dir_data.GetName()+fname).Data() << std::endl;
              if(fname.EndsWith("Data.root")) 
              {
                TFile *data = TFile::Open((std::string(dir_data.GetName())+"/"+std::string(fname)).c_str());
                dataPOT += util::GetIngredient<TParameter<double>>(*data, "POTUsed")->GetVal();
                std::cout << "Data POT: " << dataPOT <<std::endl;
              }
          }
      }
  }
  TSystemDirectory dir_mc(mcdir.c_str(), mcdir.c_str());
  TList *mc_files = dir_mc.GetListOfFiles();
  if (mc_files) {
      TSystemFile *file;
      TString fname;
      TIter next(mc_files);
      while ((file=(TSystemFile*)next())) {
          fname = file->GetName();
          if (!file->IsDirectory()) {
              std::cout << (dir_mc.GetName()+fname).Data() << std::endl;
              if(fname.EndsWith("combinedMC.root"))
              {
                TFile *mc = TFile::Open((std::string(dir_mc.GetName())+"/"+std::string(fname)).c_str());
                mcPOT += util::GetIngredient<TParameter<double>>(*mc, "POTUsed")->GetVal();
                std::cout << "MC POT: " << mcPOT <<std::endl;
              }
          }
      }
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

  if(argc != 4)
  {
    std::cerr << "Expected 3 arguments, but I got " << argc-1 << ".\n"
              << "USAGE: ExtractCrossSection <unfolding iterations> <data directory> <mc directory>\n";
    return 1;
  }
  std::cout<<"Here11\n";
  const int nIterations = std::stoi(argv[1]);

  TFile *sample_data_file;
  TFile *sample_mc_file;
  std::string datadir = std::string(argv[2]);
  std::string mcdir = std::string(argv[3]);

  TSystemDirectory dirdata(datadir.c_str(), datadir.c_str());
  TList *datafiles = dirdata.GetListOfFiles();
  std::cout<<"Here12\n";
  if (datafiles) {
      TSystemFile *file;
      TString fname;
      TIter next(datafiles);
      while ((file=(TSystemFile*)next())) {
          fname = file->GetName();
          if (!file->IsDirectory()) {
              if(fname.EndsWith("Data.root"))
                {
                  sample_data_file = TFile::Open((std::string(dirdata.GetName())+"/"+std::string(fname)).c_str());
                  break;
                }
            }
      }
  }
  std::cout<<"Here13\n";

  TSystemDirectory dirmc(mcdir.c_str(), mcdir.c_str());
  TList *mcfiles = dirmc.GetListOfFiles();
  if (mcfiles) {
      TSystemFile *file;
      TString fname;
      TIter next(mcfiles);
      while ((file=(TSystemFile*)next())) {
          fname = file->GetName();
          if (!file->IsDirectory()) {
              if(fname.EndsWith("MC.root"))
                {
                  TFile *sample_mc_file = TFile::Open((std::string(dirmc.GetName())+"/"+std::string(fname)).c_str());
                  std::cout<<"mc file: " << (std::string(dirmc.GetName())+"/"+std::string(fname)).c_str() <<std::endl;
                  break;
                }
            }
      }
  }

  if(!sample_data_file)
  {
    std::cerr << "Failed to open sample data file " << argv[2] << ". Weird\n";
    return 2;
  }
  if(!sample_mc_file)
  {
    std::cerr << "Failed to open sample MC file " << argv[3] << ". Weird\n";
    return 3;
  };
  std::cout<<"Here-11\n";
  std::vector<std::string> crossSectionPrefixes;
  for(auto key: *sample_data_file->GetListOfKeys())
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
  std::cout<<"Here-2\n";
  //const double mcPOT = util::GetIngredient<TParameter<double>>(*mcFile, "POTUsed")->GetVal(),
  //             dataPOT = util::GetIngredient<TParameter<double>>(*dataFile, "POTUsed")->GetVal();
  std::pair<double, double> POTs = getCombinedPOTs(argv[2], argv[3]); //Note !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! need to exclude the combined/madd-ed files
  const double dataPOT = POTs.first;
  const double mcPOT = POTs.second;
  std::cout<<"Data POT: " << dataPOT << " mcPOT " << mcPOT << std::endl;
  std::cout<<"Here-3\n";
  for(const auto& prefix: crossSectionPrefixes)
  {
    std::map<int, std::string> testTgt = {{3026, "3026"}};

    for(const auto& TgtCode: testTgt) 
    //for(const auto& TgtCode: util::TgtCodeLabelsNuke) 
    {
      try
      {

        PlotUtils::MnvH1D *flux;
        PlotUtils::MnvH1D *folded;
        PlotUtils::MnvH2D *migration;
        PlotUtils::MnvH1D *effNum;
        PlotUtils::MnvH1D *effDenom;
        bool t_flux = false;
        bool t_folded = false;
        bool t_migration = false; 
        bool t_effNum;
        bool t_effDenom;
        TParameter<double>* nNucleons;
        std::vector<PlotUtils::MnvH1D*> backgrounds;
        //Combine stuff
        TSystemDirectory dir_data(datadir.c_str(), datadir.c_str());
        TList *files_data = dir_data.GetListOfFiles();
        if (files_data) {
            TSystemFile *file;
            TString fname;
            TIter next(files_data);
            while ((file=(TSystemFile*)next())) {
                fname = file->GetName();
                if (!file->IsDirectory()) {
                    std::cout << (dir_data.GetName()+fname).Data() << std::endl;
                    if(fname.EndsWith("Data.root")) 
                    {
                      TFile *dataFile = TFile::Open((std::string(dir_data.GetName())+"/"+std::string(fname)).c_str());
                      folded->Add(util::GetIngredient<PlotUtils::MnvH1D>(*dataFile,  "by_TargetCode_Data_"+TgtCode.second, prefix));
                    }
                }
            }
        }
        TSystemDirectory dir_mc(mcdir.c_str(), mcdir.c_str());
        TList *files_mc = dir_mc.GetListOfFiles();
        if (files_mc) {
            TSystemFile *file;
            TString fname;
            TIter next(files_mc);
            while ((file=(TSystemFile*)next())) {
                fname = file->GetName();
                if (!file->IsDirectory()) {
                    std::cout << (dir_mc.GetName()+fname).Data() << std::endl;
                    if(fname.EndsWith("combinedMC.root"))
                    {
                      TFile *mcFile = TFile::Open((std::string(dir_mc.GetName())+"/"+std::string(fname)).c_str());
                      if (!t_flux) flux->Add(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "target"+TgtCode.second+"_reweightedflux_integrated", prefix));
                      else flux->Add(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "target"+TgtCode.second+"_reweightedflux_integrated", prefix));
                      if (!t_migration) migration->Add(util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, "migration_"+TgtCode.second, prefix));
                      else migration->Add(util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, "migration_"+TgtCode.second, prefix));
                      if (!t_effNum) effNum->Add(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "efficiency_numerator_"+TgtCode.second, prefix));
                      else effNum->Add(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "efficiency_numerator_"+TgtCode.second, prefix));
                      if (!t_effDenom) effDenom->Add(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "efficiency_denominator_"+TgtCode.second, prefix));
                      else effDenom->Add(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "efficiency_denominator_"+TgtCode.second, prefix));
                      nNucleons = (util::GetIngredient<TParameter<double>>(*mcFile, "target"+TgtCode.second+"_fiducial_nucleons", prefix)); //Dan: Use the same truth fiducial volume for all extractions.  The acceptance correction corrects data back to this fiducial even if the reco fiducial cut is different.
                      backgrounds.push_back(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "tgt"+TgtCode.second+"_NC_Bkg", prefix));
                      backgrounds.push_back(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "tgt"+TgtCode.second+"_Wrong_Material_Bkg", prefix));
                    }
                }
            }
        }


        std::cout<<"Here1\n";
        Plot(*folded, "data", prefix);
        auto simEventRate = effDenom->Clone(); //Make a copy for later

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

        //Look for backgrounds with <prefix>_<analysis>_Background_<name>
        /* for(auto key: *mcFile->GetListOfKeys())
        {
          if(std::string(key->GetName()).find(prefix + "_background_") != std::string::npos)
          {
            backgrounds.push_back(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, key->GetName()));
          }
        } */

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

        auto outFile = TFile::Open((prefix + TgtCode.second+"_crossSection.root").c_str(), "CREATE");
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

        auto crossSection = normalize(unfolded, flux, nNucleons->GetVal(), dataPOT);
        Plot(*crossSection, "crossSection", prefix);
        crossSection->Clone()->Write("crossSection");

        //Write a "simulated cross section" to compare to the data I just extracted.
        //If this analysis passed its closure test, this should be the same cross section as
        //what GENIEXSecExtract would produce.
        normalize(simEventRate, flux, nNucleons->GetVal(), mcPOT);
        
        Plot(*simEventRate, "simulatedCrossSection", prefix);
        simEventRate->Write("simulatedCrossSection");
      }
      catch(const std::runtime_error& e)
      {
        std::cerr << "Failed to extra a cross section for prefix " << prefix << ": " << e.what() << "\n";
        return 4;
      }
    }
  }

  return 0;
}
