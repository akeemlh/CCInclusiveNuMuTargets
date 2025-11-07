//File: COHPionReweighter.h
//Brief: A Reweighter changes the CV model into a different model using just a multiplicative
//       constant.  All vertical systematics are implemented by taking ratios to such weights.
//       Some Reweighters are mutually exclusive, and others are only needed for specific systematics.


#ifndef COHPion_REWEIGHTER_H
#define COHPion_REWEIGHTER_H

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
  class COHPionReweighter:public Reweighter<UNIVERSE, EVENT>
  {
    public:
      COHPionReweighter() = default;
      virtual ~COHPionReweighter() = default;

      virtual double GetWeight(const UNIVERSE& univ, const EVENT& myevent /*event*/) const override{
	if( univ.GetInt("mc_intType") == 4){
	   
 	    double weight = univ.GetCOHPionWeight();//GetCoherentPiWeight(angle, KE); //PlotUtils::weightCoherentPi().get_combined_weight(angle, KE);
            return weight;
	 }
         else return 1.0;

                
	
      };
      virtual std::string GetName() const {return "COHPionReweighter"; }

      virtual bool DependsReco() const {return false;}
      //virtual bool DependsTruth() const {return true;}; //Not needed as of time of writing.
      //PlotUtils::PionReweighter& PionReweighter();
      //virtual bool IsCompatible(const PionReweighter& /*other*/) const { return true; }
      //virtual std::vector<UNIVERSE*> GetRequiredUniverses() const { return std::vector<UNIVERSE*>{}; }
  };
}

#endif //DIFFRACTIVE_REWEIGHTER_H