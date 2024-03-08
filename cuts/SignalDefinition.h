#ifndef SIGNALDEFINITION_H
#define SIGNALDEFINITION_H
//This is what we'll do for reco
/*template <class UNIVERSE>
using q3Signal = PlotUtils::Maximum<UNIVERSE, double, &UNIVERSE::GetTrueq3>;*/ //TODO: Except I didn't actually write this for SignalConstraints :(

//PlotUtils includes
#include "PlotUtils/Cut.h"
#include "util/PlasticSidebands.h"

//Package includes
#include "event/CVUniverse.h"

template <class UNIVERSE>
class Q3Limit: public PlotUtils::SignalConstraint<UNIVERSE>
{
  public:
    Q3Limit(const double q3Max): PlotUtils::SignalConstraint<UNIVERSE>("q3 < " + std::to_string(q3Max) + " GeV"),
                                 fQ3Max(q3Max)
    {
    }

  private:
    double fQ3Max; //Maximum q3 allowed in GeV/c

    bool checkConstraint(const UNIVERSE& univ) const //override
    {
      double trueq3 = univ.Getq3True()/pow(10,3);
      return trueq3 < fQ3Max;
    }
};


namespace truth
{
    
  template <class UNIVERSE>
  class MuonEnergyMin: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      MuonEnergyMin(const double val, const std::string& name): PlotUtils::SignalConstraint<UNIVERSE>(name), fVal(val) {}

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        return univ.GetTruthMuE() > fVal;
      }
      double fVal;
  };

  template <class UNIVERSE>
  class MuonEnergyMax: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      MuonEnergyMax(const double val, const std::string& name): PlotUtils::SignalConstraint<UNIVERSE>(name), fVal(val) {}

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        return univ.GetTruthMuE() < fVal;
      }
      double fVal;
  };
}

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

  //Is this vertex upstream of a nuclear target
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class USScintillator: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
      USScintillator(): PlotUtils::Cut<UNIVERSE, EVENT>(std::string("Upstream scintillator sideband")) {}

    private:
      bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
      {
        //To do
      }
  };

  //Is this vertex downstream of a nuclear target
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class DSScintillator: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
      DSScintillator(): PlotUtils::Cut<UNIVERSE, EVENT>(std::string("Downstream scintillator sideband")) {}

    private:
      bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
      {
        //To do 
      }
  };
}
#endif //SIGNALDEFINITION_H