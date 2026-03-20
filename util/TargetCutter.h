//A slight extension/redefintion for the MAT cutter class to accomodate the plastic sidebands in a nuclear target analysis

//The existing cutter class sets all sidebands to false in checkSelection if any of the reco cuts fail. In a nuclear target
//analysis you are likely to have a selection cut for your nuclear targets, since our sideband region (the plastic inbetween targets)
//would clearly fail such a cut in the existing checkSelection there any capacity for working with these sidebands
#ifndef TARGETCUTTER_H
#define TARGETCUTTER_H

#include "PlotUtils/Cutter.h"


namespace PlotUtils
{
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class TargetCutter : PlotUtils::Cutter
  {
    public:

        TargetCutter(reco_t&& recoPre, reco_t&& recoSideband, truth_t&& truthSignal, truth_t&& truthPhaseSpace) : PlotUtils::Cutter(recoPre, recoSideband, truthSignal, truthPhaseSpace) {};
        template <class UNIVERSE, class EVENT>
        std::bitset<64> Cutter<UNIVERSE, EVENT>::checkSelection(const UNIVERSE& univ, EVENT& event, double weight, const bool isSignalForCuts)
        {
            if(!std::all_of(fRecoPreCuts.begin(), fRecoPreCuts.end(), [&univ, &event, weight, isSignalForCuts](std::unique_ptr<Cut<UNIVERSE, EVENT> >& cut)
                                                                    { return cut->passesCut(univ, event, weight, isSignalForCuts); }))
            {
            return 0; //All sideband cuts set to false
            }

            constexpr int nCutsMax = 64;
            assert(1 + fRecoSidebandCuts.size() < nCutsMax && "Cutter needs to be upgraded to support more than 63 sideband cuts.");
            std::bitset<nCutsMax> result;
            result.set(); //Set all Cuts to passed by default so user can just check isSelected().all()

            int whichCut = 0;
            for(auto& cut: fRecoSidebandCuts)
            {
            if(!cut->passesCut(univ, event, weight, isSignalForCuts))
            {
                result.set(whichCut, false);
                weight = 0; //Keep checking for sideband, but this event has not been selected.
            }
            ++whichCut;
            }

            return result;
        }
  }
}

#endif //TARGETCUTTER_H
