//File: AsitsRobFineFluxReweighter.h
//Brief: A Reweighter changes the CV model into a different model using just a multiplicative
//       constant.  All vertical systematics are implemented by taking ratios to such weights.
//       Some Reweighters are mutually exclusive, and others are only needed for specific systematics.


#ifndef ASITROBFINE_REWEIGHTER_H
#define ASITROBFINE_REWEIGHTER_H

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
  class AsitsRobFineFluxReweighter:public Reweighter<UNIVERSE, EVENT>
  {
    public:
      AsitsRobFineFluxReweighter() = default;
      virtual ~AsitsRobFineFluxReweighter() = default;
      // from MINERvA Document 33762-v3
      virtual double GetWeight(const UNIVERSE& univ, const EVENT& myevent /*event*/) const override{
            double Enu = univ.GetEnuTrueGeV();
            double weight = 1.0;
            // now only want Asit's weight above 7.5 GeV only
            if (Enu > 7.5){
                weight = univ.GetAsitsRobFineFluxWeight(); // comment back in for low nu flux!
            }
            return weight;
        }


      
      virtual std::string GetName() const {return "AsitsRobFineFluxReweighter"; }

      virtual bool DependsReco() const {return false;}
      //virtual bool DependsTruth() const {return true;}; //Not needed as of time of writing.
      //PlotUtils::PionReweighter& PionReweighter();
      //virtual bool IsCompatible(const PionReweighter& /*other*/) const { return true; }
      //virtual std::vector<UNIVERSE*> GetRequiredUniverses() const { return std::vector<UNIVERSE*>{}; }
  };
}

#endif //ASITROBFINE_REWEIGHTER_H