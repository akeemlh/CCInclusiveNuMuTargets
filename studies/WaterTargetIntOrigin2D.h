//studies includes
#include "studies/Study.h"

//Mehreen's includes
#include "event/MichelEvent.h"
#include "util/Categorized.h"
#include "event/CVUniverse.h"

//c++ includes
#include <functional> //for std::function

//We noticed that we were seeing some interesting behaviour from playlist 1C
//We were seeing events in the what was supposed to be an empty water target indiciative of it being not really empty
//This study is to investigate that in 2D

int getVtxCode(const CVUniverse& univ)
{
  int truthMod = univ.GetTruthVtxModule();
  if (truthMod==6)
  //if (mc_z>5153.77 && mc_z<5456.74)
  {
    //std::cout<<"In water\n";
    return 0;
  }
  int truthTgtCode = univ.GetTruthTargetCode();
  if (truthTgtCode>0)
  {
    //std::cout<<"In target " << truthTgtCode <<std::endl;
    return truthTgtCode;
  }
  else
  { 
    
    if ((truthMod >= 11) && (truthMod <= 14))
    {
      //std::cout<<"In US\n";
      return 7;
    }
    else if ((truthMod >= 15) && (truthMod <= 18)) 
    {
      //std::cout<<"In DS\n";
      return 8;
    }
    else
    {
      //std::cout<<"In Other\n";
      return 9;
    }
  }
}

class WaterTargetIntOrigin2D: public Study
{
  public:
    using reco_t = std::function<double(const CVUniverse&, const MichelEvent&)>;

    WaterTargetIntOrigin2D(reco_t reco_x, reco_t reco_y, const std::string& varName, const std::string& varUnits, const std::vector<double> xBins, const std::vector<double> yBins, const std::map<std::string, std::vector<CVUniverse*>>& univs): Study(), fReco_x(reco_x), fReco_y(reco_y)
    {
      std::map<int, std::string> TargetZ = {{1, "Hydrogen"},
                                                {8, "Oxygen"}};
      std::map<int, std::string> TargetVtx = {
                {0, "Water"},
                {1, "Target1"},
                {2, "Target2"},
                {3, "Target3"},
                {4, "Target4"},
                {5, "Target5"} ,
                {7, "US"},
                {8, "DS"}};

      m_VarToTargetZ = new util::Categorized<HIST, int>(varName+"_tgtZ", varName + " [" + varUnits + "]", TargetZ, xBins, yBins, univs);
      m_VarToTargetVtx = new util::Categorized<HIST, int>(varName, varName + " [" + varUnits + "]", TargetVtx, xBins, yBins, univs);
    }

    void SaveOrDraw(TDirectory& outDir)
    {
      std::cout<<"Saving Study Histograms\n";
       outDir.cd();
       m_VarToTargetZ->visit([](HIST& wrapper)
                                {
                                  wrapper.SyncCVHistos();
                                  wrapper.hist->Write();
                                });
       m_VarToTargetVtx->visit([](HIST& wrapper)
                                {
                                  wrapper.SyncCVHistos();
                                  wrapper.hist->Write();
                                });
       //TODO: You could do plotting here
    }

  private:
    using HIST = PlotUtils::Hist2DWrapper<CVUniverse>;

    reco_t fReco_x, fReco_y;

    util::Categorized<HIST, int>* m_VarToTargetZ;
    util::Categorized<HIST, int>* m_VarToTargetVtx;

    //Overriding base class functions
    //Do nothing for now...  Good place for data comparisons in the future. 
    void fillSelected(const CVUniverse& univ, const MichelEvent& evt, const double weight){}

    //All of your plots happen here so far.
    void fillSelectedSignal(const CVUniverse& univ, const MichelEvent& evt, const double weight)
    {
        //if reconstructed in water
        if (univ.GetANNSegment()==36 || univ.GetANNVtxModule()==6)
        {
          int vtxCode = getVtxCode(univ);
          /* std::cout<<"getVtxCode(univ): " << vtxCode <<std::endl;
          std::cout<<"univ.GetTruthTargetID() " << univ.GetTruthTargetID() <<std::endl;
          std::cout<<"univ.GetTruthTargetCode() " << univ.GetTruthTargetCode() <<std::endl;
          std::cout<<"univ.GetTruthVtxModule() " << univ.GetTruthVtxModule() <<std::endl;
          std::cout<<"univ.GetTrueVertex().Z()() " << univ.GetTrueVertex().Z() <<std::endl; */

          (*m_VarToTargetZ)[univ.GetMCTargetZ()].FillUniverse(&univ, fReco_x(univ, evt), fReco_y(univ, evt), weight);
          (*m_VarToTargetVtx)[vtxCode].FillUniverse(&univ, fReco_x(univ, evt), fReco_y(univ, evt), weight);
        }
    }

    //Do nothing for now...  Good place for efficiency denominators in the future.
    void fillTruthSignal(const CVUniverse& univ, const MichelEvent& evt, const double weight) {}
};
