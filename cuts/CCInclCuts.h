#ifndef CCINCLCUTS_H
#define CCINCLCUTS_H
//PlotUtils includes
#include "PlotUtils/Cut.h"
#include "util/PlasticSidebands.h"

//Package includes
#include "event/CVUniverse.h"
namespace reco
{

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class ANNConfidenceCut: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
      ANNConfidenceCut(const double conf): PlotUtils::Cut<UNIVERSE, EVENT>(std::string("ANN confidence > ") + std::to_string(conf)), fConf(conf) {}

    private:
      bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
      {
        return univ.GetANNProb() > fConf; 
      }
      const double fConf;
  };

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class HasMLPredCut: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
      HasMLPredCut(): PlotUtils::Cut<UNIVERSE, EVENT>(std::string("Has ML Prediction")){}

    private:
      bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
      {
        return univ.hasMLPred(); 
      }
  };

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  using MuonEnergyMax = PlotUtils::Maximum<UNIVERSE, double, &UNIVERSE::GetEmu, EVENT>;

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class ANNDistToTgtDivCut: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
      ANNDistToTgtDivCut(const double dist): PlotUtils::Cut<UNIVERSE, EVENT>(std::string("Interaction Distance from target division > ") + std::to_string(dist)), fDistCut(dist) {}

    private:
      bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
      {
        int tgtID = univ.GetANNTargetID();
        if (tgtID == 1 || tgtID == 2 || tgtID == 3 || tgtID == 5) return (univ.GetANNTargetDistToDiv() >= fDistCut); //Only perform chech for targets that have this division
        else return true;
      }
      const double fDistCut;
  };
}
#endif //CCINCLCUTS_H