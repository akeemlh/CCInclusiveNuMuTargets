#ifndef VARIABLE1DNUKE_H
#define VARIABLE1DNUKE_H

//Includes from this package
#include "event/CVUniverse.h"
#include "util/SafeROOTName.h"
#include "util/Categorized.h"

//PlotUtils includes
#include "PlotUtils/VariableBase.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "util/NukeUtils.h"

class Variable1DNuke: public PlotUtils::VariableBase<CVUniverse>
{
  private:
    typedef PlotUtils::HistWrapper<CVUniverse> Hist;
    typedef PlotUtils::Hist2DWrapper<CVUniverse> Hist2D;
  public:
    template <class ...ARGS>
    Variable1DNuke(ARGS... args): PlotUtils::VariableBase<CVUniverse>(args...){}

    //Map of target codes to investigate, default, all of them
    std::map<int, std::string> TgtCodeLabels = util::TgtCodeLabelsNuke; // = {{-1, "Tracker"}, {1026, "1026"}, {1082, "1082"}, {2026, "2026"}, {2082, "2082"}, {3006, "3006"}, {3026, "3026"}, {3082, "3082"}, {4082, "4082"}, {5026, "5026"}, {5082, "5082"}, {-999, "Water"}, {7, "Target7"}, {8, "Target8"}, {9, "Target9"}, {10, "Target10"}, {11, "Target11"}, {12, "Target12"}};

    //TODO: It's really silly to have to make 2 sets of error bands just because they point to different trees.
    //      I'd rather the physics of the error bands remain the same and just change which tree they point to.
    void InitializeMCHists(std::map<std::string, std::vector<CVUniverse*>>& mc_error_bands,
                           std::map<std::string, std::vector<CVUniverse*>>& truth_error_bands)
    {

      m_HistsByTgtCodeMC = new util::Categorized<Hist, int>((GetName() + "_by_TargetCode_MC").c_str(),
        GetName().c_str(), TgtCodeLabels,
        GetBinVec(), mc_error_bands);

      for(auto& target: util::TargetNums)
      {
        //For each target set the categorised sets of histograms to store the US MC sideband distrubtions
        m_sidebandHistSetUSMC.insert({target.first, 
              new util::Categorized<Hist, int>((GetName() + std::string("_tgt") + target.second + std::string("_US_sideband")).c_str(),
              GetName().c_str(), util::SidebandCategories,
              GetBinVec(), mc_error_bands)
          });
        //For each target set the categorised sets of histograms to store the DS MC sideband distrubtions
        m_sidebandHistSetDSMC.insert({target.first, 
              new util::Categorized<Hist, int>((GetName() + std::string("_tgt") + target.second + std::string("_DS_sideband")).c_str(),
              GetName().c_str(), util::SidebandCategories,
              GetBinVec(), mc_error_bands)
          });
      }

      for(auto& tgtCode: TgtCodeLabels)
      {
        //For each target set the histogram to store the interaction channel
        m_intChannelsByTgtCode.insert({tgtCode.first, 
              new util::Categorized<Hist, int>((GetName() + std::string("_tgt") + tgtCode.second).c_str(),
              GetName().c_str(), util::GENIELabels,
              GetBinVec(), mc_error_bands)
          });

        //For each target set the histogram to store the backgrounds
        m_bkgsByTgtCode.insert({tgtCode.first, 
              new util::Categorized<Hist, int>((GetName() + std::string("_tgt") + tgtCode.second).c_str(),
              GetName().c_str(), util::BKGLabels,
              GetBinVec(), mc_error_bands)
          });
      }

      m_HistsByTgtCodeEfficiencyNumerator  = new util::Categorized<Hist, int>((GetName() + "_efficiency_numerator").c_str(),
        GetName().c_str(), TgtCodeLabels,
        GetBinVec(), mc_error_bands);
      m_HistsByTgtCodeEfficiencyDenominator  = new util::Categorized<Hist, int>((GetName() + "_efficiency_denominator").c_str(),
        GetName().c_str(), TgtCodeLabels,
        GetBinVec(), truth_error_bands);
      m_HistsByTgtCodeMigration  = new util::Categorized<PlotUtils::Hist2DWrapper<CVUniverse>, int>((GetName() + "_migration").c_str(),
        GetName().c_str(), TgtCodeLabels, GetBinVec(), GetBinVec(), mc_error_bands);

      //efficiencyNumerator = new Hist((GetNameX() + "_" + GetNameY() + "_efficiency_numerator").c_str(), GetName().c_str(), GetBinVecX(), GetBinVecY(), mc_error_bands);
      //efficiencyDenominator = new Hist((GetNameX() + "_" + GetNameY() + "_efficiency_denominator").c_str(), GetName().c_str(), GetBinVecX(), GetBinVecY(), truth_error_bands);
      //selectedMCReco = new Hist((GetName() + "_selected_mc_reco").c_str(), GetName().c_str(), GetBinVecX(), GetBinVecY(), mc_error_bands);
    }

    //Histograms to be filled
    //These histograms plot the events that we reconstruct as being WITHIN a nuclear target
    util::Categorized<Hist, int>* m_HistsByTgtCodeMC; ////-
    util::Categorized<Hist, int>* m_HistsByTgtCodeData; ////-

    //These histograms plot the distrubution of interaction channels per target code 
    std::map<int, util::Categorized<Hist, int>* > m_intChannelsByTgtCode; ////-

    //These histograms plot the distrubution of background channels per target code 
    std::map<int, util::Categorized<Hist, int>* > m_bkgsByTgtCode; ////-

    //These histograms plot the events that we reconstruct as being WITHIN a nuclear target
    //For each US or DS plane we want a set of hists to store where it really came from
    //For each event reconstructed within an US plane we store the real event vertex 
    std::map<int, util::Categorized<Hist, int>* > m_sidebandHistSetUSMC; ////-
    //For each event reconstructed within an DS plane we store the real event vertex 
    std::map<int, util::Categorized<Hist, int>* > m_sidebandHistSetDSMC; ////-

    //For each US or DS plane in reco/data we want to save just the events we see
    //These histograms plot the events that we reconstruct as being UPSTREAM of a nuclear target
    util::Categorized<Hist, int>* m_sidebandHistsUSData; ////-
    //These histograms plot the events that we reconstruct as being DOWNSTREAM of a nuclear target
    util::Categorized<Hist, int>* m_sidebandHistsDSData; ////-
    //No equivalent for MC since we can simply get all the MC upstream and downstream events by summing the m_sidebandHistSetUSMC and m_sidebandHistSetDSMC 


    util::Categorized<Hist, int>* m_HistsByTgtCodeEfficiencyNumerator;
    util::Categorized<Hist, int>* m_HistsByTgtCodeEfficiencyDenominator;
    util::Categorized<Hist2D, int>* m_HistsByTgtCodeMigration;
    //Hist* dataHist;  
    //Hist* efficiencyNumerator;
    //Hist* efficiencyDenominator;
    //Hist* selectedMCReco; //Treat the MC CV just like data for the closure test

    void InitializeDATAHists(std::vector<CVUniverse*>& data_error_bands)
    {

      std::string strName = GetName();
      const char* name = strName.c_str();
      std::cout << "Hist base name: " <<name<< std::endl;

      m_HistsByTgtCodeData = new util::Categorized<Hist, int>((GetName() + "_by_TargetCode_Data").c_str(),
        GetName().c_str(), TgtCodeLabels,
        GetBinVec(), data_error_bands);

      m_sidebandHistsUSData = new util::Categorized<Hist, int>((GetName() + "_US_sideband_by_Target_Data").c_str(),
        GetName().c_str(), util::TargetNums,
        GetBinVec(), data_error_bands);

      m_sidebandHistsDSData = new util::Categorized<Hist, int>((GetName() + "_DS_sideband_by_Target_Data").c_str(),
        GetName().c_str(), util::TargetNums,
        GetBinVec(), data_error_bands);
    }

    void WriteData(TFile& file)
    {

      m_sidebandHistsUSData->visit([&file](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(&file);
                                      categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                    });
      m_sidebandHistsDSData->visit([&file](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(&file);
                                      categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                    });
      m_HistsByTgtCodeData->visit([&file](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(&file);
                                      categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                    });
    }

    void WriteMC(TFile& file)
    {
      SyncCVHistos();
      file.cd();

      for(auto& histSet: m_sidebandHistSetUSMC)
      {
        histSet.second->visit([&file](Hist& categ)
                                      {
                                        categ.hist->SetDirectory(&file);
                                        categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                      });
      }
      for(auto& histSet: m_sidebandHistSetDSMC)
      {
        histSet.second->visit([&file](Hist& categ)
                                      {
                                        categ.hist->SetDirectory(&file);
                                        categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                      });
      }
      m_HistsByTgtCodeMC->visit([&file](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(&file);
                                      categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                    });

      for(auto& histSet: m_intChannelsByTgtCode)
      {
        histSet.second->visit([&file](Hist& categ)
                                      {
                                        categ.hist->SetDirectory(&file);
                                        categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                      });
      }
      for(auto& histSet: m_bkgsByTgtCode)
      {
        histSet.second->visit([&file](Hist& categ)
                                      {
                                        categ.hist->SetDirectory(&file);
                                        categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                      });
      }
      m_HistsByTgtCodeEfficiencyNumerator->visit([&file](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(&file);
                                      categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                    });
      m_HistsByTgtCodeEfficiencyDenominator->visit([&file](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(&file);
                                      categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                    });
      m_HistsByTgtCodeMigration->visit([&file](Hist2D& categ)
                                    {
                                      categ.hist->SetDirectory(&file);
                                      categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                    });
    }
    //Only call this manually if you Draw(), Add(), or Divide() plots in this
    //program.
    //Makes sure that all error bands know about the CV.  In the Old Systematics
    //Framework, this was implicitly done by the event loop.
    void SyncCVHistos()
    {
      /*m_backgroundHists->visit([](Hist& categ) { categ.SyncCVHistos(); });
      if(dataHist) dataHist->SyncCVHistos();
      if(efficiencyNumerator) efficiencyNumerator->SyncCVHistos();
      if(efficiencyDenominator) efficiencyDenominator->SyncCVHistos();
      if(selectedSignalReco) selectedSignalReco->SyncCVHistos();
      if(selectedMCReco) selectedMCReco->SyncCVHistos();
      if(migration) migration->SyncCVHistos();*/
    }
};

#endif //VARIABLE1DNUKE_H