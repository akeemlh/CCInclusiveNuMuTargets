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
  using MuonEnergyMax = PlotUtils::Maximum<UNIVERSE, double, &UNIVERSE::GetEmu, EVENT>;

}
#endif //CCINCLCUTS_H