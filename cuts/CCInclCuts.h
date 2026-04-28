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

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class Helicity: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
      Helicity(const int pdg): PlotUtils::Cut<UNIVERSE, EVENT>(getCutName(pdg)), fPDG(pdg) {}

    private:
      bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
      {
        if (fPDG == 14)
        { 
          return univ.GetHelicity() == 1;
        }
        else if (fPDG == -14)
        {
          return univ.GetHelicity() == 2;
        }
      }
      std::string getCutName(int pdg)
      {
        if (fPDG == 14)
        { 
          return std::string("Is Neutrino");
        }
        else if (fPDG == -14)
        {
          return std::string("Is Antineutrino");
        }
      }
      const double fPDG;
  };
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class MinMuonAngle: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
      MinMuonAngle(const double max): PlotUtils::Cut<UNIVERSE, EVENT>(std::string("Muon Theta > ") + std::to_string(max)), fMax(max*M_PI/180.) {}

    private:
      bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
      {
        return univ.GetThetamu() >= fMax;
      }

      const double fMax; //radians
  };
}
#endif //CCINCLCUTS_H