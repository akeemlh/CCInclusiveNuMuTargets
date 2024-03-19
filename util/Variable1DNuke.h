#ifndef VARIABLE1DNUKE
#define VARIABLE1DNUKE

//Includes from this package
#include "event/CVUniverse.h"
#include "util/SafeROOTName.h"
#include "util/Categorized.h"

//PlotUtils includes
#include "PlotUtils/VariableBase.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "util/Variable.h"


class Variable1DNuke: public Variable
{
  private:
    typedef PlotUtils::HistWrapper<CVUniverse> Hist;
  public:
    template <class ...ARGS>
    Variable1DNuke(ARGS... args): Variable(args...)
    {
    }

    std::map<int, std::string> BKGLabels = {{0, "NC"},
                {1, "Wrong_Sign"}};
    std::map<int, std::string> TargetNums = {{1, "1"}, {2, "2"}, {3, "3"}, {4, "4"}, {5, "5"}, {6, "6"}};
    std::map<int, std::string> SidebandCategories = {{0, "US"}, {1, "DS"}, {2, "Signal"}};
    std::map<int, std::string> TgtCodeLabels = {{1026, "1026"}, {1082, "1082"}, {2026, "2026"}, {2082, "2082"}, {3006, "3006"}, {3026, "3026"}, {3082, "3082"}, {4082, "4082"}, {5026, "5026"}, {5082, "5082"}, {-999, "Water"}};
    std::map<int, std::string> GENIELabels = {{1, "QE"},
                                              {8, "2p2h"},
                                              {2, "RES"},
                                              {3, "DIS"}};

    //TODO: It's really silly to have to make 2 sets of error bands just because they point to different trees.
    //      I'd rather the physics of the error bands remain the same and just change which tree they point to.
    void InitializeMCHists(std::map<std::string, std::vector<CVUniverse*>>& mc_error_bands,
                           std::map<std::string, std::vector<CVUniverse*>>& truth_error_bands)
    {

      //For each target set the histogram to store the interaction channel
      m_bkgsByTgtCode.insert({tgtCode.first, 
            new util::Categorized<Hist, int>((GetName() + std::string("_bkg_tgt") + tgtCode.second).c_str(),
            GetName().c_str(), BKGLabels,
            GetBinVec(), mc_error_bands)
      });


      efficiencyNumerator = new Hist((GetName() + "_efficiency_numerator").c_str(), GetName().c_str(), GetBinVec(), mc_error_bands);
      efficiencyDenominator = new Hist((GetName() + "_efficiency_denominator").c_str(), GetName().c_str(), GetBinVec(), truth_error_bands);
      selectedSignalReco = new Hist((GetName() + "_selected_signal_reco").c_str(), GetName().c_str(), GetBinVec(), mc_error_bands);
      selectedMCReco = new Hist((GetName() + "_selected_mc_reco").c_str(), GetName().c_str(), GetBinVec(), mc_error_bands);
      migration = new PlotUtils::Hist2DWrapper<CVUniverse>((GetName() + "_migration").c_str(), GetName().c_str(), GetBinVec(), GetBinVec(), mc_error_bands);
    
      for(auto& target: TargetNums)
      {
        //For each target set the categorised sets of histograms to store the US MC sideband distrubtions
        m_sidebandHistSetUSMC.insert({target.first, 
              new util::Categorized<Hist, int>((GetName() + std::string("_tgt") + target.second + std::string("_US_sideband")).c_str(),
              GetName().c_str(), SidebandCategories,
              GetBinVec(), mc_error_bands)
          });
        //For each target set the categorised sets of histograms to store the DS MC sideband distrubtions
        m_sidebandHistSetDSMC.insert({target.first, 
              new util::Categorized<Hist, int>((GetName() + std::string("_tgt") + target.second + std::string("_DS_sideband")).c_str(),
              GetName().c_str(), SidebandCategories,
              GetBinVec(), mc_error_bands)
          });
      }

      for(auto& tgtCode: TgtCodeLabels)
      {
        //For each target set the histogram to store the interaction channel
        m_intChannelsByTgtCode.insert({tgtCode.first, 
              new util::Categorized<Hist, int>((GetName() + std::string("_tgt") + tgtCode.second).c_str(),
              GetName().c_str(), GENIELabels,
              GetBinVec(), mc_error_bands)
          });
      }

      for(auto& tgtCode: TgtCodeLabels)
      {
        //For each target set the histogram to store the interaction channel
        m_bkgsByTgtCode.insert({tgtCode.first, 
              new util::Categorized<Hist, int>((GetName() + std::string("_bkg_tgt") + tgtCode.second).c_str(),
              GetName().c_str(), BKGLabels,
              GetBinVec(), mc_error_bands)
          });
      }

    }

    //Histograms to be filled
    //These histograms plot the distrubution of background channels per target code 
    std::map<int, util::Categorized<Hist, int>* > m_backgroundHistsByTgtCode; ////-

    Hist* dataHist;
    Hist* efficiencyNumerator;
    Hist* efficiencyDenominator;
    Hist* selectedSignalReco; //Effectively "true background subtracted" distribution for warping studies.
                              //Also useful for a bakground breakdown plot that you'd use to start background subtraction studies.
    Hist* selectedMCReco; //Treat the MC CV just like data for the closure test
    PlotUtils::Hist2DWrapper<CVUniverse>* migration;




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




    void InitializeDATAHists(std::vector<CVUniverse*>& data_error_bands)
    {
      dataHist = new Hist((GetName() + "_data").c_str(), GetName().c_str(), GetBinVec(), data_error_bands);
          m_HistsByTgtCodeData = new util::Categorized<Hist, int>((GetName() + "_by_TargetCode_Data").c_str(),
        GetName().c_str(), TgtCodeLabels,
        GetBinVec(), data_error_bands);

      m_sidebandHistsUSData = new util::Categorized<Hist, int>((GetName() + "_US_sideband_by_Target_Data").c_str(),
        GetName().c_str(), TargetNums,
        GetBinVec(), data_error_bands);

      m_sidebandHistsDSData = new util::Categorized<Hist, int>((GetName() + "_DS_sideband_by_Target_Data").c_str(),
        GetName().c_str(), TargetNums,
        GetBinVec(), data_error_bands);

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

      if(migration)
      {
        migration->hist->SetDirectory(&file); 
        migration->hist->Write();
      }

      if(selectedSignalReco)
      {
        selectedSignalReco->hist->SetDirectory(&file);
        selectedSignalReco->hist->Write();
      }

      if(selectedMCReco)
      {
        selectedMCReco->hist->SetDirectory(&file);
        selectedMCReco->hist->Write((GetName() + "_data").c_str()); //Make this histogram look just like the data for closure tests
      }
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
    }

    //Only call this manually if you Draw(), Add(), or Divide() plots in this
    //program.
    //Makes sure that all error bands know about the CV.  In the Old Systematics
    //Framework, this was implicitly done by the event loop.
    void SyncCVHistos()
    {
      if(dataHist) dataHist->SyncCVHistos();
      if(efficiencyNumerator) efficiencyNumerator->SyncCVHistos();
      if(efficiencyDenominator) efficiencyDenominator->SyncCVHistos();
      if(selectedSignalReco) selectedSignalReco->SyncCVHistos();
      if(selectedMCReco) selectedMCReco->SyncCVHistos();
      if(migration) migration->SyncCVHistos();
    }
};

#endif //VARIABLE_H




#ifndef VARIABLE_H
#define VARIABLE_H

//Includes from this package
#include "event/CVUniverse.h"
#include "util/SafeROOTName.h"
#include "util/Categorized.h"

//PlotUtils includes
#include "PlotUtils/VariableBase.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/Hist2DWrapper.h"

class Variable1DNuke: public Variable<CVUniverse>
{
  private:
    typedef PlotUtils::HistWrapper<CVUniverse> Hist;
  public:
    template <class ...ARGS>
    Variable1DNuke(ARGS... args): PlotUtils::Variable<CVUniverse>(args...)
    {
    }

    std::map<int, std::string> BKGLabels = {{0, "NC"},
                {1, "Wrong_Sign"}};
    std::map<int, std::string> TargetNums = {{1, "1"}, {2, "2"}, {3, "3"}, {4, "4"}, {5, "5"}, {6, "6"}};
    std::map<int, std::string> SidebandCategories = {{0, "US"}, {1, "DS"}, {2, "Signal"}};
    std::map<int, std::string> TgtCodeLabels = {{1026, "1026"}, {1082, "1082"}, {2026, "2026"}, {2082, "2082"}, {3006, "3006"}, {3026, "3026"}, {3082, "3082"}, {4082, "4082"}, {5026, "5026"}, {5082, "5082"}, {-999, "Water"}};
    std::map<int, std::string> GENIELabels = {{1, "QE"},
                                              {8, "2p2h"},
                                              {2, "RES"},
                                              {3, "DIS"}};

    //TODO: It's really silly to have to make 2 sets of error bands just because they point to different trees.
    //      I'd rather the physics of the error bands remain the same and just change which tree they point to.
    void InitializeMCHists(std::map<std::string, std::vector<CVUniverse*>>& mc_error_bands,
                           std::map<std::string, std::vector<CVUniverse*>>& truth_error_bands)
    {

      efficiencyNumerator = new Hist((GetName() + "_efficiency_numerator").c_str(), GetName().c_str(), GetBinVec(), mc_error_bands);
      efficiencyDenominator = new Hist((GetName() + "_efficiency_denominator").c_str(), GetName().c_str(), GetBinVec(), truth_error_bands);
      selectedSignalReco = new Hist((GetName() + "_selected_signal_reco").c_str(), GetName().c_str(), GetBinVec(), mc_error_bands);
      selectedMCReco = new Hist((GetName() + "_selected_mc_reco").c_str(), GetName().c_str(), GetBinVec(), mc_error_bands);
      migration = new PlotUtils::Hist2DWrapper<CVUniverse>((GetName() + "_migration").c_str(), GetName().c_str(), GetBinVec(), GetBinVec(), mc_error_bands);
    }

    //Histograms to be filled
    Hist* dataHist;
    Hist* efficiencyNumerator;
    Hist* efficiencyDenominator;
    Hist* selectedSignalReco; //Effectively "true background subtracted" distribution for warping studies.
                              //Also useful for a bakground breakdown plot that you'd use to start background subtraction studies.
    Hist* selectedMCReco; //Treat the MC CV just like data for the closure test
    PlotUtils::Hist2DWrapper<CVUniverse>* migration;

    void InitializeDATAHists(std::vector<CVUniverse*>& data_error_bands)
    {
      dataHist = new Hist((GetName() + "_data").c_str(), GetName().c_str(), GetBinVec(), data_error_bands);
    }

    void WriteData(TFile& file)
    {
      if (dataHist->hist) {
                dataHist->hist->SetDirectory(&file);
                dataHist->hist->Write();
      }
    }

    void WriteMC(TFile& file)
    {
      SyncCVHistos();
      file.cd();


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

      if(migration)
      {
        migration->hist->SetDirectory(&file); 
        migration->hist->Write();
      }

      if(selectedSignalReco)
      {
        selectedSignalReco->hist->SetDirectory(&file);
        selectedSignalReco->hist->Write();
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
      if(selectedSignalReco) selectedSignalReco->SyncCVHistos();
      if(selectedMCReco) selectedMCReco->SyncCVHistos();
      if(migration) migration->SyncCVHistos();
    }
};
#endif //VARIABLE1DNUKE
