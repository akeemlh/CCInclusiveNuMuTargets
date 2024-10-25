#ifndef VARIABLE2DTRACKER_H
#define VARIABLE2DTRACKER_H

#include "util/SafeROOTName.h"
#include "PlotUtils/Variable2DBase.h"
#include "util/Categorized.h"
#include "util/NukeUtils.h"
#include "MinervaUnfold/MnvResponse.h"

class Variable2DTracker: public PlotUtils::Variable2DBase<CVUniverse>
{
  private:
    typedef PlotUtils::Hist2DWrapper<CVUniverse> Hist;
  public:
    template <class ...ARGS>
    Variable2DTracker(ARGS... args): PlotUtils::Variable2DBase<CVUniverse>(args...)
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

      m_interactionChannels = new util::Categorized<Hist, int>((GetName() + "_by_int_channel").c_str(),
							   GetName().c_str(), util::GENIELabels,
							   GetBinVecX(), GetBinVecY(), mc_error_bands);

      efficiencyNumerator = new Hist((GetNameX() + "_" + GetNameY() + "_efficiency_numerator").c_str(), GetName().c_str(), GetBinVecX(), GetBinVecY(), mc_error_bands);
      efficiencyDenominator = new Hist((GetNameX() + "_" + GetNameY() + "_efficiency_denominator").c_str(), GetName().c_str(), GetBinVecX(), GetBinVecY(), truth_error_bands);
      selectedMCReco = new Hist((GetName() + "_selected_mc_reco").c_str(), GetName().c_str(), GetBinVecX(), GetBinVecY(), mc_error_bands);
      
      // taken from Noah
      // make a temp universe map to make Response happy
      std::map<const std::string, int> response_bands; // necessary?
      for (auto band : mc_error_bands){
        std::string name = band.first;
        const std::string realname = (band.second)[0]->ShortName();
        int nuniv = band.second.size();

        response_bands[realname] = nuniv;
      }

      migration = new MinervaUnfold::MnvResponse((GetName() + "_migration").c_str(), GetName().c_str(), GetBinVecX(), GetBinVecY() , response_bands);
      
    
    }

    //Histograms to be filled
    util::Categorized<Hist, int>* m_backgroundHists;
    util::Categorized<Hist, int>* m_interactionChannels;
    Hist* dataHist;  
    Hist* efficiencyNumerator;
    Hist* efficiencyDenominator;
    Hist* selectedMCReco; //Treat the MC CV just like data for the closure test
    MinervaUnfold::MnvResponse* migration;
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

    TH1D* xsamplehist = new TH1D("x", "x", GetBinVecX().size()-1, &GetBinVecX()[0]);
    TH1D* ysamplehist = new TH1D("y", "y", GetBinVecY().size()-1, &GetBinVecY()[0]);

    void WriteMC(TFile& file)
    {
      SyncCVHistos();
      file.cd();

      m_backgroundHists->visit([&file](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(&file);
                                      categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                    });

      m_interactionChannels->visit([&file](Hist& categ)
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

      if(selectedMCReco)
      {
        selectedMCReco->hist->SetDirectory(&file);
        selectedMCReco->hist->Write((GetName() + "_data").c_str()); //Make this histogram look just like the data for closure tests
      }
      if(migration)
      {
        MnvH2D* migration_hist;
        MnvH2D* reco_hist;
        MnvH2D* truth_hist;
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
      m_interactionChannels->visit([](Hist& categ) { categ.SyncCVHistos(); });
      if(dataHist) dataHist->SyncCVHistos();
      if(efficiencyNumerator) efficiencyNumerator->SyncCVHistos();
      if(efficiencyDenominator) efficiencyDenominator->SyncCVHistos();
      if(selectedMCReco) selectedMCReco->SyncCVHistos();
    }
};

#endif //VARIABLE2DTRACKER_H
