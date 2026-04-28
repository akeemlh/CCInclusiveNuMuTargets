//File: Q2TpiReweighter.h
//Brief: A Reweighter changes the CV model into a different model using just a multiplicative
//       constant.  All vertical systematics are implemented by taking ratios to such weights.
//       Some Reweighters are mutually exclusive, and others are only needed for specific systematics.

#ifndef Q2Tpi_REWEIGHTER_H
#define Q2Tpi_REWEIGHTER_H

#if __cplusplus < 201103L
  #define override
#endif

//c++ includes
#include <string>
#include <vector>
#include "PlotUtils/Reweighter.h"
#include "PlotUtils/WeightFunctions.h"
namespace PlotUtils
{

  template <class UNIVERSE, class EVENT = MichelEvent>
  class Q2TpiReweighter:public Reweighter<UNIVERSE, EVENT>
  {
    public:
      Q2TpiReweighter() = default;
      virtual ~Q2TpiReweighter() = default;

      virtual double GetWeight(const UNIVERSE& univ, const EVENT& myevent /*event*/) const override{
        return univ.GetCCOnePiWeight();
      };
      virtual std::string GetName() const {return "Q2TpiReweighterReweighter"; }

      virtual bool DependsReco() const {return false;}
      //virtual bool DependsTruth() const {return true;}; //Not needed as of time of writing.
      //PlotUtils::PionReweighter& PionReweighter();
      //virtual bool IsCompatible(const PionReweighter& /*other*/) const { return true; }
      //virtual std::vector<UNIVERSE*> GetRequiredUniverses() const { return std::vector<UNIVERSE*>{}; }
  };
}

#endif //DIFFRACTIVE_REWEIGHTER_H