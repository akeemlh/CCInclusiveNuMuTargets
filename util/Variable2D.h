#ifndef VARIABLE2D_H
#define VARIABLE2D_H

#include "util/SafeROOTName.h"
#include "PlotUtils/Variable2DBase.h"
#include "util/Categorized.h"

class Variable2D: public PlotUtils::Variable2DBase<CVUniverse>
{
  private:
    typedef PlotUtils::Hist2DWrapper<CVUniverse> Hist;
  public:
    template <class ...ARGS>
    Variable2D(ARGS... args): PlotUtils::Variable2DBase<CVUniverse>(args...)
    {
    }

    //TODO: It's really silly to have to make 2 sets of error bands just because they point to different trees.
    //      I'd rather the physics of the error bands remain the same and just change which tree they point to.
    void InitializeMCHists(std::map<std::string, std::vector<CVUniverse*>>& mc_error_bands,
                           std::map<std::string, std::vector<CVUniverse*>>& truth_error_bands)
    {

      std::map<int, std::string> BKGLabels = {{0, "NC_Bkg"},
					       {1, "Bkg_Wrong_Sign"}};
      
      m_backgroundHists = new util::Categorized<Hist, int>((GetName() + "_by_BKG_Label").c_str(),
							   GetName().c_str(), BKGLabels,
							   GetBinVecX(), GetBinVecY(), mc_error_bands);
      
      m_HistsByTgtCodeMC = new util::Categorized<Hist, int>((GetName() + "_by_TargetCode_MC").c_str(),
        GetName().c_str(), TgtCodeLabels,
        GetBinVecX(), GetBinVecY(), mc_error_bands);

      for(auto& target: TargetNums)
      {
        //For each target set the categorised sets of histograms to store the US MC sideband distrubtions
        m_sidebandHistSetUSMC.insert({target.first, 
              new util::Categorized<Hist, int>((GetName() + std::string("_tgt") + target.second + std::string("_US_sideband")).c_str(),
              GetName().c_str(), SidebandCategories,
              GetBinVecX(), GetBinVecY(), mc_error_bands)
          });
        //For each target set the categorised sets of histograms to store the DS MC sideband distrubtions
        m_sidebandHistSetDSMC.insert({target.first, 
              new util::Categorized<Hist, int>((GetName() + std::string("_tgt") + target.second + std::string("_DS_sideband")).c_str(),
              GetName().c_str(), SidebandCategories,
              GetBinVecX(), GetBinVecY(), mc_error_bands)
          });
      }

      for(auto& tgtCode: TgtCodeLabels)
      {
        //For each target set the histogram to store the interaction channel
        m_intChannelsByTgtCode.insert({tgtCode.first, 
              new util::Categorized<Hist, int>((GetName() + std::string("_tgt") + tgtCode.second).c_str(),
              GetName().c_str(), GENIELabels,
              GetBinVecX(), GetBinVecY(), mc_error_bands)
          });
      }

      for(auto& tgtCode: TgtCodeLabels)
      {
        //For each target set the histogram to store the interaction channel
        m_bkgsByTgtCode.insert({tgtCode.first, 
              new util::Categorized<Hist, int>((GetName() + std::string("_bkg_tgt") + tgtCode.second).c_str(),
              GetName().c_str(), BKGLabels,
              GetBinVecX(), GetBinVecY(), mc_error_bands)
          });
      }

      efficiencyNumerator = new Hist((GetNameX() + "_" + GetNameY() + "_efficiency_numerator").c_str(), GetName().c_str(), GetBinVecX(), GetBinVecY(), mc_error_bands);
      efficiencyDenominator = new Hist((GetNameX() + "_" + GetNameY() + "_efficiency_denominator").c_str(), GetName().c_str(), GetBinVecX(), GetBinVecY(), truth_error_bands);
      selectedMCReco = new Hist((GetName() + "_selected_mc_reco").c_str(), GetName().c_str(), GetBinVecX(), GetBinVecY(), mc_error_bands);
    }

    //Histograms to be filled
    util::Categorized<Hist, int>* m_backgroundHists; //m_sidebandHistSetDSMC

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


    Hist* dataHist;  
    Hist* efficiencyNumerator;
    Hist* efficiencyDenominator;
    Hist* selectedMCReco; //Treat the MC CV just like data for the closure test

    //Format: "(Material*1000)ID" ie 
    std::map<int, std::string> TargetNums = {{1, "1"}, {2, "2"}, {3, "3"}, {4, "4"}, {5, "5"}, {6, "6"}};
    std::map<int, std::string> SidebandCategories = {{0, "US"}, {1, "DS"}, {2, "Signal"}};
    std::map<int, std::string> TgtCodeLabels = {{1026, "1026"}, {1082, "1082"}, {2026, "2026"}, {2082, "2082"}, {3006, "3006"}, {3026, "3026"}, {3082, "3082"}, {4082, "4082"}, {5026, "5026"}, {5082, "5082"}, {-999, "Water"}};
    std::map<int, std::string> GENIELabels = {{1, "QE"},
                                              {8, "2p2h"},
                                              {2, "RES"},
                                              {3, "DIS"}};

    void InitializeDATAHists(std::vector<CVUniverse*>& data_error_bands)
    {
      std::string strName = GetName();
      const char* name = strName.c_str();
      std::cout << "Hist base name: " <<name<< std::endl;
      dataHist = new Hist(Form("_data_%s", name), name, GetBinVecX(), GetBinVecY(), data_error_bands);

      m_HistsByTgtCodeData = new util::Categorized<Hist, int>((GetName() + "_by_TargetCode_Data").c_str(),
        GetName().c_str(), TgtCodeLabels,
        GetBinVecX(), GetBinVecY(), data_error_bands);

      m_sidebandHistsUSData = new util::Categorized<Hist, int>((GetName() + "_US_sideband_by_Target_Data").c_str(),
        GetName().c_str(), TargetNums,
        GetBinVecX(), GetBinVecY(), data_error_bands);

      m_sidebandHistsDSData = new util::Categorized<Hist, int>((GetName() + "_DS_sideband_by_Target_Data").c_str(),
        GetName().c_str(), TargetNums,
        GetBinVecX(), GetBinVecY(), data_error_bands);
    }

    void WriteData(TFile& file)
    {
      if (dataHist->hist) {
        dataHist->hist->SetDirectory(&file);
        dataHist->hist->Write();
      }

      std::size_t found = GetName().find("tracker");
      if (found==std::string::npos) //If this isn't a tracker variable
      {
        m_HistsByTgtCodeData->visit([&file](Hist& categ)
                                      {
                                        categ.hist->SetDirectory(&file);
                                        categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                      });
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
      }
    }

    void WriteMC(TFile& file)
    {
      SyncCVHistos();
      file.cd();

      m_backgroundHists->visit([&file](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(&file);
                                      categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                    });

      std::size_t found = GetName().find("tracker");
      if (found==std::string::npos) //If this isn't a tracker variable
      {
        m_HistsByTgtCodeMC->visit([&file](Hist& categ)
                                      {
                                        categ.hist->SetDirectory(&file);
                                        categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                      });
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
      }

      if(efficiencyNumerator)
      {
        efficiencyNumerator->hist->SetDirectory(&file); //TODO: Can I get around having to call SetDirectory() this many times somehow?
        efficiencyNumerator->hist->Write();
      }

      if(efficiencyDenominator)
      {
        efficiencyDenominator->hist->SetDirectory(&file);
        efficiencyDenominator->hist->Write();
      }
      if(selectedMCReco)
      {
        selectedMCReco->hist->SetDirectory(&file);
        selectedMCReco->hist->Write((GetName() + "_data").c_str()); //Make this histogram look just like the data for closure tests
      }
    }

    //Only call this manually if you Draw(), Add(), or Divide() plots in this
    //program.
    //Makes sure that all error bands know about the CV.  In the Old Systematics
    //Framework, this was implicitly done by the event loop.
    void SyncCVHistos()
    {
      m_backgroundHists->visit([](Hist& categ) { categ.SyncCVHistos(); });
      if(dataHist) dataHist->SyncCVHistos();
      if(efficiencyNumerator) efficiencyNumerator->SyncCVHistos();
      if(efficiencyDenominator) efficiencyDenominator->SyncCVHistos();
    }
};

#endif //VARIABLE2D_H
