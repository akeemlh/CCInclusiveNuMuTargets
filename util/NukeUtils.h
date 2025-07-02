#ifndef UTIL_NUKEUTILS_H
#define UTIL_NUKEUTILS_H

#include <map>
#include <algorithm>
#include <string>

#include "PlotUtils/Cut.h"

namespace util
{
    //Just used to store objects useful for nuclear target study
    //From https://cdcvs.fnal.gov/redmine/projects/minerva/wiki/Z_Positions_of_Planes_in_the_Full_MINERvA_Detector
    //And https://cdcvs.fnal.gov/redmine/projects/minerva-sw/wiki/MINERvA_Detector
    //And https://cdcvs.fnal.gov/redmine/projects/minerva-ops/wiki/MINERvA_Detector
    //And https://cdcvs.fnal.gov/redmine/projects/minerva/wiki/Z_Positions_of_Planes_in_the_Downstream_Detector -- Not used
    
    std::map<int, std::string> TargetNums = {{1, "1"}, {2, "2"}, {3, "3"}, {4, "4"}, {5, "5"}, {6, "6"}};
    std::map<int, std::string> SidebandCategories = {{0, "US"}, {1, "DS"}, {2, "Signal"}};
    std::map<int, std::string> TgtCodeLabelsTracker = {{-1, "Tracker"}};
    std::map<int, std::string> TgtCodeLabelsNuke = {{1026, "1026"}, {1082, "1082"}, {2026, "2026"}, {2082, "2082"}, {3006, "3006"}, {3026, "3026"}, {3082, "3082"}, {4082, "4082"}, {5026, "5026"}, {5082, "5082"}, {6000, "Water"}, {7, "Target7"}, {8, "Target8"}, {9, "Target9"}, {10, "Target10"}, {11, "Target11"}, {12, "Target12"}};

    //Without pseudotargets
    std::map<int, std::string> TgtCodeLabelsNukeNoPS = {{1026, "1026"}, {1082, "1082"}, {2026, "2026"}, {2082, "2082"}, {3006, "3006"}, {3026, "3026"}, {3082, "3082"}, {4082, "4082"}, {5026, "5026"}, {5082, "5082"}, {6000, "Water"}};

    std::map<int, std::string> GENIELabels = {{1, "QE"},
                                              {8, "2p2h"},
                                              {2, "RES"},
                                              {3, "DIS"}}; //kNoInteraction=0, kQEL=1, kRES=2, kDIS=3, kCOHPI=4, kAMNUGAMMA=5, kIMD=6, kNUEEL=7, k2P2H=8

    std::map<int, std::string> BKGLabels = {{0, "NC_Bkg"},
					       {1, "Wrong_Sign_Bkg"}, {2, "US_Plastic"}, {3, "DS_Plastic"}};

    int nuOrAntiNuMode(std::string playlist)
    {
        std::vector<std::string> nuVector = {"minervame1A", "minervame1B", "minervame1C", "minervame1D", "minervame1E", "minervame1F", "minervame1G", "minervame1L", "minervame1M", "minervame1N", "minervame1O", "minervame1P"};
        std::vector<std::string> anuVector = {"minervame5A", "minervame6A", "minervame6B", "minervame6C", "minervame6D", "minervame6E", "minervame6F", "minervame6G", "minervame6H", "minervame6I", "minervame6J"};
        if (std::find(std::begin(nuVector), std::end(nuVector), playlist) != std::end(nuVector)) return 1;
        else if (std::find(std::begin(anuVector), std::end(anuVector), playlist) != std::end(anuVector)) return 2;
        //Is this the most efficient way? Probably not
        return 0;
    }
}
namespace reco
{
    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class ZRangeANN: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
        public:
        ZRangeANN(const std::string& name, const double zMin, const double zMax): PlotUtils::Cut<UNIVERSE, EVENT>(name), fMin(zMin), fMax(zMax)
        {
        }

        private:
        bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
        {
            return univ.GetANNVertex().Z() >= fMin && univ.GetANNVertex().Z() <= fMax;
        } 

        const double fMin;
        const double fMax;
    };

    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class MuonCurveSignificance: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
        public:
        MuonCurveSignificance( const double zMin): PlotUtils::Cut<UNIVERSE, EVENT>("Muon Curvature significance in MINOS"), fMin(zMin)
        {
        }

        private:
        bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
        {
            //Don't make a significance cut if we reconstructed by range
            if(univ.GetInt((univ.GetAnaToolName() + "_minos_used_curvature").c_str()) != 1) return true;
            double relativeErr = univ.GetMuonQPErr();
            if (univ.GetAnalysisNuPDG()>0) return ( relativeErr <= -fMin );
            else if (univ.GetAnalysisNuPDG()<0) return ( relativeErr >= fMin );
        } 

        const double fMin;
    };
}
#endif //UTIL_NUKEUTILS_H