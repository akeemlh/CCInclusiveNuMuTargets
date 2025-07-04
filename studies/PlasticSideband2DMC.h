//studies includes
#include "studies/Study.h"

//Mehreen's includes
#include "event/MichelEvent.h"
#include "util/Categorized.h"
#include "event/CVUniverse.h"
#include "util/PlasticSidebands.h"

//c++ includes
#include <functional> //for std::function

class PlasticSideband2DMC: public Study
{
  public:
    //PerMichelVarByGENIELabel fills a histogram with 1 entry per Michel with some variable calculated from that Michel.  Your function will get to see the CVUniverse, the MichelEvent (= reconstructed Michels), and which Michel it's looping over.
    using reco_t = std::function<double(const CVUniverse&, const MichelEvent&)>;

    PlasticSideband2DMC(reco_t reco_x, reco_t reco_y, const std::string& varName, const std::string& varUnits, const std::vector<double> xBins, const std::vector<double> yBins, const std::map<std::string, std::vector<CVUniverse*>>& univs): Study(), fReco_x(reco_x), fReco_y(reco_y)
    {
      m_SelectedMCRecoByTgtCode = new util::Categorized<Hist, int>((varName + "_by_TargetCode_MC").c_str(),
        varName.c_str(), util::TgtCodeLabels,
        GetBinVecX(), GetBinVecY(), mc_error_bands);

      for(auto& target: util::TargetNums)
      {
        //For each target set the categorised sets of histograms to store the US MC sideband distrubtions
        m_sidebandHistSetUSMC.insert({target.first, 
              new util::Categorized<Hist, int>((varName + std::string("_tgt") + target.second + std::string("_US_sideband")).c_str(),
              varName.c_str(), util::SidebandCategories,
              GetBinVecX(), GetBinVecY(), mc_error_bands)
          });
        //For each target set the categorised sets of histograms to store the DS MC sideband distrubtions
        m_sidebandHistSetDSMC.insert({target.first, 
              new util::Categorized<Hist, int>((varName + std::string("_tgt") + target.second + std::string("_DS_sideband")).c_str(),
              varName.c_str(), util::SidebandCategories,
              GetBinVecX(), GetBinVecY(), mc_error_bands)
          });
      }

    }

    void SaveOrDraw(TDirectory& outDir)
    {
      outDir.cd();

      m_SelectedMCRecoByTgtCode->visit([&file](Hist& categ)
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

       //TODO: You could do plotting here
    }

    std::map<int, std::string> SidebandCategories;
    //These histograms plot the events that we reconstruct as being WITHIN a nuclear target
    util::Categorized<Hist, int>* m_SelectedMCRecoByTgtCode; ////-


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
    void fillSelected(const CVUniverse& univ, const MichelEvent& evt, const double weight) {}

    //All of your plots happen here so far.
    void fillSelectedSignal(const CVUniverse& univ, const MichelEvent& evt, const double weight)
    {
      //To do: use universe->hasMLPred()
      //Nuke Target Study
      int annTgtCode = univ->GetANNTargetCode();
      //If this has a segment num 36 it came from water target
      bool inWaterSegment = (univ->GetANNSegment()==36);
      if (annTgtCode>00 || inWaterSegment) //If this event occurs inside a nuclear target
      {
        int code = inWaterSegment ? -999 : annTgtCode;
        //Plot events that occur within the nuclear targets grouped by which target they occur in
        (m_SelectedMCRecoByTgtCode)[code].FillUniverse(univ, univfReco_x(univ, evt), fReco_y(univ, evt), 1);
      }
      else
      {
        int tmpModCode = (univ->GetANNVtxModule()*10)+univ->GetANNVtxPlane();
        auto USTgtID = util::USModPlaCodeToTgtId.find(tmpModCode);
        auto DSTgtID = util::DSModPlaCodeToTgtId.find(tmpModCode);
        if (USTgtID != util::USModPlaCodeToTgtId.end()) //Is event reconstructed immediately upstream of a nuclear target
        {
          //std::cout<<"Found in US target "<< USTgtID->second<< "\n";
          //Check where it really interacted
          int tmpTruthModCode = (univ->GetTruthVtxModule()*10)+univ->GetTruthVtxPlane();
          //Checking if this interaction truthfully originated in an upstream or downstream sideband
          //Checking if this interaction truthfully originated in the neigbouring nuclear target
          //If this has a non-zero target code then it actually originated in a nuke target
          //If this has a segment num 36 it came from water target
          int truthTgtID = univ->GetTruthTargetID();
          if (truthTgtID > 0 || inWaterSegment) 
          {
            (m_sidebandHistSetUSMC[USTgtID->second])[2].FillUniverse(univ, univfReco_x(univ, evt), fReco_y(univ, evt), 1);
          }
          else if (util::isUSPlane(tmpModCode)>0) //If originated immediately upstream of nuke target
          {
            (m_sidebandHistSetUSMC[USTgtID->second])[0].FillUniverse(univ, univfReco_x(univ, evt), fReco_y(univ, evt), 1);
          }
          else if (util::isDSPlane(tmpModCode)>0) //If originated immediately downstream of nuke target
          {
            (m_sidebandHistSetUSMC[USTgtID->second])[1].FillUniverse(univ, univfReco_x(univ, evt), fReco_y(univ, evt), 1);
          }
          else //Fill "Other" histogram if this event didn't really have a vtx in the nuke target or sideband
          {
            (m_sidebandHistSetUSMC[USTgtID->second])[-1].FillUniverse(univ, univfReco_x(univ, evt), fReco_y(univ, evt), 1);
          }
        }
        else if (DSTgtID != util::DSModPlaCodeToTgtId.end()) //Or is event reconstructed immediately downstream of a nuclear target
        {
          //Check where it really interacted
          int tmpTruthModCode = (univ->GetTruthVtxModule()*10)+univ->GetTruthVtxPlane();
          //Checking if this interaction truthfully originated in an upstream or downstream sideband
          //Checking if this interaction truthfully originated in the neigbouring nuclear target
          //If this has a non-zero target code then it actually originated in a nuke target
          //If this has a segment num 36 it came from water target
          int truthTgtID = univ->GetTruthTargetID();
          if (truthTgtID > 0 || inWaterSegment) 
          {
            (m_sidebandHistSetDSMC[DSTgtID->second])[2].FillUniverse(univ, univfReco_x(univ, evt), fReco_y(univ, evt), 1);
          }
          else if (util::isUSPlane(tmpModCode)>0) //Is event reconstructed immediately upstream of a nuclear target
          {
            (m_sidebandHistSetDSMC[DSTgtID->second])[0].FillUniverse(univ, univfReco_x(univ, evt), fReco_y(univ, evt), 1);
          }
          else if (util::isDSPlane(tmpModCode)>0) //Is event reconstructed immediately upstream of a nuclear target
          {
            (m_sidebandHistSetDSMC[DSTgtID->second])[1].FillUniverse(univ, univfReco_x(univ, evt), fReco_y(univ, evt), 1);
          }
          else //Fill "Other" histogram if this event didn't really have a vtx in the nuke target or sideband
          {
            (m_sidebandHistSetDSMC[DSTgtID->second])[-1].FillUniverse(univ, univfReco_x(univ, evt), fReco_y(univ, evt), 1);
          }
        }
      }
      //End - Nuke Target Study
    }

    //Do nothing for now...  Good place for efficiency denominators in the future.
    void fillTruthSignal(const CVUniverse& univ, const MichelEvent& evt, const double weight) {}
};
