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

class Variable: public PlotUtils::VariableBase<CVUniverse>
{
  private:
    typedef PlotUtils::HistWrapper<CVUniverse> Hist;
  public:
    template <class ...ARGS>
    Variable(ARGS... args): PlotUtils::VariableBase<CVUniverse>(args...)
    {
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

    Hist* EffNumDaisy[12];
    Hist* EffDenomDaisy[12];
    Hist* selectedMCRecoDaisy[12], *selectedSignalRecoDaisy[12];
    Hist* dataDaisy[12];
    PlotUtils::Hist2DWrapper<CVUniverse>* MigrationDaisy[12];
    util::Categorized<Hist, int>* BackgroundsDaisy[12];

    //TODO: It's really silly to have to make 2 sets of error bands just because they point to different trees.
    //      I'd rather the physics of the error bands remain the same and just change which tree they point to.
    void InitializeMCHists(std::map<std::string, std::vector<CVUniverse*>>& mc_error_bands,
                           std::map<std::string, std::vector<CVUniverse*>>& truth_error_bands)
    {

      std::map<int, std::string> BKGLabels = {{0, "NC"},
					       {1, "Wrong_Sign"}};
      
      m_backgroundHists = new util::Categorized<Hist, int>((GetName() + "_background").c_str(),
							   GetName().c_str(), BKGLabels,
							   GetBinVec(), mc_error_bands);

      efficiencyNumerator = new Hist((GetName() + "_efficiency_numerator").c_str(), GetName().c_str(), GetBinVec(), mc_error_bands);
      efficiencyDenominator = new Hist((GetName() + "_efficiency_denominator").c_str(), GetName().c_str(), GetBinVec(), truth_error_bands);
      selectedSignalReco = new Hist((GetName() + "_selected_signal_reco").c_str(), GetName().c_str(), GetBinVec(), mc_error_bands);
      selectedMCReco = new Hist((GetName() + "_selected_mc_reco").c_str(), GetName().c_str(), GetBinVec(), mc_error_bands);
      migration = new PlotUtils::Hist2DWrapper<CVUniverse>((GetName() + "_migration").c_str(), GetName().c_str(), GetBinVec(), GetBinVec(), mc_error_bands);

      for (int petal = 0; petal<12; petal++)
      {
        EffNumDaisy[petal] = new Hist((GetName() + "_Daisy_EffNum_"+petal), GetName().c_str(), GetBinVec(), mc_error_bands);
        EffDenomDaisy[petal] = new Hist((GetName() + "_Daisy_EffDenom_"+petal), GetName().c_str(), GetBinVec(), truth_error_bands);
        MigrationDaisy[petal] = new PlotUtils::Hist2DWrapper<CVUniverse>(std::string(GetName() + "_Daisy_Migration_"+petal).c_str(), GetName().c_str(), GetBinVec(), GetBinVec(), mc_error_bands);
        BackgroundsDaisy[petal] = new util::Categorized<Hist, int>((std::string(GetName() + "_Daisy_Background_"+petal).c_str()),
							   GetName().c_str(), BKGLabels,
							   GetBinVec(), mc_error_bands);
        selectedMCRecoDaisy[petal]  = new Hist((GetName() + "_Daisy_selected_mc_reco_"+petal), GetName().c_str(), GetBinVec(), mc_error_bands);
        selectedSignalRecoDaisy[petal]  = new Hist((GetName() + "_Daisy_selected_signal_reco_"+petal), GetName().c_str(), GetBinVec(), mc_error_bands);
      }
    }
    
    void InitializeDATAHists(std::vector<CVUniverse*>& data_error_bands)
    {
      dataHist = new Hist((GetName() + "_data").c_str(), GetName().c_str(), GetBinVec(), data_error_bands);
      for (int petal = 0; petal<12; petal++)
      {
        dataDaisy[petal] = new Hist((GetName() + "_Daisy_Data_"+petal), GetName().c_str(), GetBinVec(), data_error_bands);
      }
    }

    void WriteData(TFile& file)
    {
      if (dataHist->hist) {
                dataHist->hist->SetDirectory(&file);
                dataHist->hist->Write();
      }
      for (int daisy = 0; daisy<12; daisy++)
      {
        if (dataDaisy[daisy])
        {
          dataDaisy[daisy]->hist->SetDirectory(&file);
          dataDaisy[daisy]->hist->Write();
        }
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

      for (int petal = 0; petal<12; petal++)
      {
        if (EffNumDaisy[petal])
        {
          EffNumDaisy[petal]->hist->SetDirectory(&file);
          EffNumDaisy[petal]->hist->Write();
        }
        if (EffDenomDaisy[petal])
        {
          EffDenomDaisy[petal]->hist->SetDirectory(&file);
          EffDenomDaisy[petal]->hist->Write();
        }
        if (MigrationDaisy[petal])
        {
          MigrationDaisy[petal]->hist->SetDirectory(&file);
          MigrationDaisy[petal]->hist->Write();
        }
        if (selectedMCRecoDaisy[petal])
        {
          selectedMCRecoDaisy[petal]->hist->SetDirectory(&file);
          selectedMCRecoDaisy[petal]->hist->Write((GetName() + "_Daisy_Data_" + petal)); //Make this histogram look just like the data for closure tests
        }
        if (selectedSignalRecoDaisy[petal])
        {
          selectedSignalRecoDaisy[petal]->hist->SetDirectory(&file);
          selectedSignalRecoDaisy[petal]->hist->Write();
        }
        BackgroundsDaisy[petal]->visit([&file](Hist& categ)
                                      {
                                        categ.hist->SetDirectory(&file);
                                        categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                      });
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

      //What about daisy?
    }
};

#endif //VARIABLE_H
