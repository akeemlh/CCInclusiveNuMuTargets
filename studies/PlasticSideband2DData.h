//studies includes
#include "studies/Study.h"

//Mehreen's includes
#include "event/MichelEvent.h"
#include "util/Categorized.h"
#include "event/CVUniverse.h"
#include "util/PlasticSidebands.h"

//c++ includes
#include <functional> //for std::function

class PlasticSideband2DData: public Study
{
  public:
    //PerMichelVarByGENIELabel fills a histogram with 1 entry per Michel with some variable calculated from that Michel.  Your function will get to see the CVUniverse, the MichelEvent (= reconstructed Michels), and which Michel it's looping over.
    using reco_t = std::function<double(const CVUniverse&, const MichelEvent&)>;

    PlasticSideband2DData(reco_t reco_x, reco_t reco_y, const std::string& varName, const std::string& varUnits, const std::vector<double> xBins, const std::vector<double> yBins, const std::map<std::string, std::vector<CVUniverse*>>& univs): Study(), fReco_x(reco_x), fReco_y(reco_y)
    {
      m_HistsByTgtCodeData = new util::Categorized<Hist, int>((varName + "_by_TargetCode_Data").c_str(),
        varName.c_str(), TgtCodeLabels,
        GetBinVecX(), GetBinVecY(), data_error_bands);

      m_sidebandHistsUSData = new util::Categorized<Hist, int>((varName + "_US_sideband_by_Target_Data").c_str(),
        varName.c_str(), TargetNums,
        GetBinVecX(), GetBinVecY(), data_error_bands);

      m_sidebandHistsDSData = new util::Categorized<Hist, int>((varName + "_DS_sideband_by_Target_Data").c_str(),
        varName.c_str(), TargetNums,
        GetBinVecX(), GetBinVecY(), data_error_bands);
    }

    void SaveOrDraw(TDirectory& outDir)
    {
      outDir.cd();
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
      //TODO: You could do plotting here
    }

    //These histograms plot the events that we reconstruct as being WITHIN a nuclear target
    util::Categorized<Hist, int>* m_HistsByTgtCodeMC; ////-


    //These histograms plot the events that we reconstruct as being WITHIN a nuclear target
    //For each US or DS plane we want a set of hists to store where it really came from
    //For each event reconstructed within an US plane we store the real event vertex 
    std::map<int, util::Categorized<Hist, int>* > m_sidebandHistSetUSMC; ////-
    //For each event reconstructed within an DS plane we store the real event vertex 
    std::map<int, util::Categorized<Hist, int>* > m_sidebandHistSetDSMC; ////-

  private:
    using HIST = PlotUtils::Hist2DWrapper<CVUniverse>;

    reco_t fReco_x, fReco_y;

    //Overriding base class functions
    //Do nothing for now...  Good place for data comparisons in the future. 
    void fillSelected(const CVUniverse& univ, const MichelEvent& evt, const double weight) 
    {
      //Nuke Target Study
      int annTgtCode = univ->GetANNTargetCode();
      //If this has a segment num 36 it came from water target
      bool inWaterSegment = (univ->GetANNSegment()==36);
      if (annTgtCode>00 || inWaterSegment) //If this event occurs inside a nuclear target
      {
        int code = inWaterSegment ? -999 : annTgtCode;
        //Plot events that occur within the nuclear targets grouped by which target they occur in
        (m_HistsByTgtCodeData)[code].FillUniverse(univ, fReco_x(univ, evt), fReco_y(univ, evt), 1);
      }
      else
      {
        int tmpModCode = (univ->GetANNVtxModule()*10)+univ->GetANNVtxPlane();
        int USModNum = util::isUSPlane(tmpModCode);
        int DSModNum = util::isDSPlane(tmpModCode);
        if (USModNum>0) //Is event reconstructed immediately upstream of a nuclear target
        {
          (m_sidebandHistsUSData)[USModNum].FillUniverse(univ, fReco_x(univ, evt), fReco_y(univ, evt), 1);
        }
        else if (DSModNum>0) //Or is event reconstructed immediately downstream of a nuclear target
        {
          (m_sidebandHistsDSData)[DSModNum].FillUniverse(univ, fReco_x(univ, evt), fReco_y(univ, evt), 1);
        }
      }
      //End - Nuke Target Study
    }
    //All of your plots happen here so far.
    void fillSelectedSignal(const CVUniverse& univ, const MichelEvent& evt, const double weight){}

    //Do nothing for now...  Good place for efficiency denominators in the future.
    void fillTruthSignal(const CVUniverse& univ, const MichelEvent& evt, const double weight) {}
};
