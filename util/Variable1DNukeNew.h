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
  public:
    template <class ...ARGS>
    Variable1DNuke(ARGS... args): PlotUtils::VariableBase<CVUniverse>(args...)
    {
    }

    //TODO: It's really silly to have to make 2 sets of error bands just because they point to different trees.
    //      I'd rather the physics of the error bands remain the same and just change which tree they point to.
    void InitializeMCHists(std::map<std::string, std::vector<CVUniverse*>>& mc_error_bands,
                           std::map<std::string, std::vector<CVUniverse*>>& truth_error_bands)
    {
      
      m_backgroundHists = new util::Categorized<Hist, int>((GetName() + "_background").c_str(),
							   GetName().c_str(), util::BKGLabels,
							   GetBinVec(), mc_error_bands);

      m_intChannelsEffDenom = new util::Categorized<Hist, int>((GetName() + "_efficiency_denominator_intChannels"),
              GetName().c_str(), util::GENIELabels,
              GetBinVec(), truth_error_bands);

      m_interactionTypeHists = new util::Categorized<Hist, int>((GetName() + "_intType").c_str(),
        GetName().c_str(), util::GENIELabels,
        GetBinVec(), mc_error_bands);

      efficiencyNumerator = new Hist((GetName() + "_efficiency_numerator").c_str(), GetName().c_str(), GetBinVec(), mc_error_bands);
      efficiencyDenominator = new Hist((GetName() + "_efficiency_denominator").c_str(), GetName().c_str(), GetBinVec(), truth_error_bands);
      selectedSignalReco = new Hist((GetName() + "_selected_signal_reco").c_str(), GetName().c_str(), GetBinVec(), mc_error_bands);
      selectedMCReco = new Hist((GetName() + "_selected_mc_reco").c_str(), GetName().c_str(), GetBinVec(), mc_error_bands);
      migration = new PlotUtils::Hist2DWrapper<CVUniverse>((GetName() + "_migration").c_str(), GetName().c_str(), GetBinVec(), GetBinVec(), mc_error_bands);
    
      m_sidebandHistSetUSMC = new util::Categorized<Hist, int>((GetName() + "_US_sideband").c_str(),
              GetName().c_str(), util::SidebandCategories,
              GetBinVec(), mc_error_bands);
      m_sidebandHistSetDSMC = new util::Categorized<Hist, int>((GetName() + "_DS_sideband").c_str(),
              GetName().c_str(), util::SidebandCategories,
              GetBinVec(), mc_error_bands);
    }

    //Histograms to be filled
    util::Categorized<Hist, int>* m_backgroundHists;
    Hist* dataHist;
    Hist* efficiencyNumerator;
    Hist* efficiencyDenominator;
    Hist* selectedSignalReco; //Effectively "true background subtracted" distribution for warping studies.
                              //Also useful for a bakground breakdown plot that you'd use to start background subtraction studies.
    Hist* selectedMCReco; //Treat the MC CV just like data for the closure test
    PlotUtils::Hist2DWrapper<CVUniverse>* migration;

    //These histograms plot the events that we reconstruct as being WITHIN a nuclear target
    //For each US or DS plane we want a set of hists to store where it really came from
    //For each event reconstructed within an US plane we store the real event vertex 
    util::Categorized<Hist, int>* m_sidebandHistSetUSMC; ////-
    //For each event reconstructed within an DS plane we store the real event vertex 
    util::Categorized<Hist, int>* m_sidebandHistSetDSMC; ////-

    //For each US or DS plane in reco/data we want to save just the events we see
    //These histograms plot the events that we reconstruct as being UPSTREAM of a nuclear target
    Hist* m_US_Sideband_Data; ////-
    //These histograms plot the events that we reconstruct as being DOWNSTREAM of a nuclear target
    Hist* m_DS_Sideband_Data; ////-
    //No equivalent for MC since we can simply get all the MC upstream and downstream events by summing the m_sidebandHistSetUSMC and m_sidebandHistSetDSMC 
    
    //These histograms plot the distrubution of interaction channels
    util::Categorized<Hist, int>*  m_interactionTypeHists;
    util::Categorized<Hist, int>* m_intChannelsEffDenom;

    void InitializeDATAHists(std::vector<CVUniverse*>& data_error_bands)
    {
      dataHist = new Hist((GetName() + "_data").c_str(), GetName().c_str(), GetBinVec(), data_error_bands);
      m_US_Sideband_Data = new Hist((GetName() + "_US_Sideband").c_str(), GetName().c_str(), GetBinVec(), data_error_bands);
      m_DS_Sideband_Data = new Hist((GetName() + "_DS_Sideband").c_str(), GetName().c_str(), GetBinVec(), data_error_bands);
    }

    void WriteData(TFile& file)
    {
      if (dataHist->hist) {
                dataHist->hist->SetDirectory(&file);
                dataHist->hist->Write();
      }
      if (m_US_Sideband_Data->hist) {
                m_US_Sideband_Data->hist->SetDirectory(&file);
                m_US_Sideband_Data->hist->Write();
      }
      if (m_DS_Sideband_Data->hist) {
                m_DS_Sideband_Data->hist->SetDirectory(&file);
                m_DS_Sideband_Data->hist->Write();
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

      m_sidebandHistSetUSMC->visit([&file](Hist& categ)
                              {
                                categ.hist->SetDirectory(&file);
                                categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                              });

      m_sidebandHistSetDSMC->visit([&file](Hist& categ)
                              {
                                categ.hist->SetDirectory(&file);
                                categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                              });

      m_intChannelsEffDenom->visit([&file](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(&file);
                                      categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                    });

      m_interactionTypeHists->visit([&file](Hist& categ)
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
      m_backgroundHists->visit([](Hist& categ) { categ.SyncCVHistos(); });
      m_sidebandHistSetUSMC->visit([](Hist& categ) { categ.SyncCVHistos(); });
      m_sidebandHistSetDSMC->visit([](Hist& categ) { categ.SyncCVHistos(); });
      m_interactionTypeHists->visit([](Hist& categ) { categ.SyncCVHistos(); });
      m_intChannelsEffDenom->visit([](Hist& categ) { categ.SyncCVHistos(); });
      if(dataHist) dataHist->SyncCVHistos();
      if(m_US_Sideband_Data) m_US_Sideband_Data->SyncCVHistos();
      if(m_DS_Sideband_Data) m_DS_Sideband_Data->SyncCVHistos();
      if(efficiencyNumerator) efficiencyNumerator->SyncCVHistos();
      if(efficiencyDenominator) efficiencyDenominator->SyncCVHistos();
      if(selectedSignalReco) selectedSignalReco->SyncCVHistos();
      if(selectedMCReco) selectedMCReco->SyncCVHistos();
      if(migration) migration->SyncCVHistos();
    }
};

#endif //VARIABLE1DNUKE_H