#ifndef VARIABLE2D_H
#define VARIABLE2D_H

#include "util/SafeROOTName.h"
#include "PlotUtils/Variable2DBase.h"
#include "util/Categorized.h"

#include "MinervaUnfold/MnvResponse.h"
#include "util/NukeUtils.h"

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

      efficiencyNumerator = new Hist((GetNameX() + "_" + GetNameY() + "_efficiency_numerator").c_str(), GetName().c_str(), GetBinVecX(), GetBinVecY(), mc_error_bands);
      efficiencyDenominator = new Hist((GetNameX() + "_" + GetNameY() + "_efficiency_denominator").c_str(), GetName().c_str(), GetBinVecX(), GetBinVecY(), truth_error_bands);
      selectedMCReco = new Hist((GetName() + "_selected_mc_reco").c_str(), GetName().c_str(), GetBinVecX(), GetBinVecY(), mc_error_bands);
      selectedSignalReco = new Hist((GetName() + "_selected_signal_reco").c_str(), GetName().c_str(), GetBinVecX(), GetBinVecY(), mc_error_bands);
      
      // taken from Noah
      // make a temp universe map to make Response happy
      std::map<const std::string, int> response_bands; // necessary?
      for (auto band : mc_error_bands){
        std::string name = band.first;
        const std::string realname = (band.second)[0]->ShortName();
        int nuniv = band.second.size();

        response_bands[realname] = nuniv;
      }

      migration  = new MinervaUnfold::MnvResponse(GetName().c_str(),
        GetName().c_str(),
        GetBinVecX(), GetBinVecY(), response_bands);    


      m_intChannelsEffDenom = new util::Categorized<Hist, int>((GetName() + "_efficiency_denominator_intChannels"),
              GetName().c_str(), util::GENIELabels,
              GetBinVecX(), GetBinVecY(), truth_error_bands);

      m_intChannels = new util::Categorized<Hist, int>((GetName() + "_intChannels"),
              GetName().c_str(), util::GENIELabels,
              GetBinVecX(), GetBinVecY(), mc_error_bands);

    }

    //Histograms to be filled
    util::Categorized<Hist, int>* m_backgroundHists;
    Hist* dataHist;  
    Hist* efficiencyNumerator;
    Hist* efficiencyDenominator;
    Hist* selectedSignalReco; //Effectively "true background subtracted" distribution for warping studies.
                              //Also useful for a bakground breakdown plot that you'd use to start background subtraction studies.

    Hist* selectedMCReco; //Treat the MC CV just like data for the closure test

    MinervaUnfold::MnvResponse* migration;

    //These histograms plot the distrubution of interaction channels
    util::Categorized<Hist, int>* m_intChannels; ////-
    util::Categorized<Hist, int>* m_intChannelsEffDenom;

    void InitializeDATAHists(std::vector<CVUniverse*>& data_error_bands)
    {
  	  dataHist = new Hist((std::string("_data_") + GetName()).c_str(), GetName().c_str(), GetBinVecX(), GetBinVecY(), data_error_bands);
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

      m_backgroundHists->visit([&file](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(&file);
                                      categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                    });
      m_intChannels->visit([&file](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(&file);
                                      categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                    });

      m_intChannelsEffDenom->visit([&file](Hist& categ)
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

    void WriteMigration(TFile& file)
    {
      //SyncCVHistos();
      std::cout<<"Writing 2D migration matrices\n";
      //file->cd();
      {
        MnvH2D* migration_hist = NULL;
        MnvH2D* reco_hist = NULL;
        MnvH2D* truth_hist = NULL;
        migration->GetMigrationObjects(migration_hist, reco_hist, truth_hist);
        migration_hist->SetDirectory(&file); 
        migration_hist->Write();
        reco_hist->SetDirectory(&file); 
        reco_hist->Write();
        truth_hist->SetDirectory(&file); 
        truth_hist->Write();
      }
    }


    //Only call this manually if you Draw(), Add(), or Divide() plots in this
    //program.
    //Makes sure that all error bands know about the CV.  In the Old Systematics
    //Framework, this was implicitly done by the event loop.
    void SyncCVHistos()
    {
      m_backgroundHists->visit([](Hist& categ) { categ.SyncCVHistos(); });
      m_intChannels->visit([](Hist& categ) { categ.SyncCVHistos(); });
      m_intChannelsEffDenom->visit([](Hist& categ) { categ.SyncCVHistos(); });
      if(dataHist) dataHist->SyncCVHistos();
      if(efficiencyNumerator) efficiencyNumerator->SyncCVHistos();
      if(efficiencyDenominator) efficiencyDenominator->SyncCVHistos();
      if(selectedMCReco) selectedMCReco->SyncCVHistos();
    }
};

#endif //VARIABLE2D_H
