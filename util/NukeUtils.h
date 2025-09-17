#ifndef UTIL_NUKEUTILS_H
#define UTIL_NUKEUTILS_H

#include <map>
#include <algorithm>
#include <string>
#include <filesystem>

#include "PlotUtils/Cut.h"
#include "event/MichelEvent.h"
#include "PlotUtils/Cutter.h"
#include "PlotUtils/TargetUtils.h"
#include "cuts/SignalDefinition.h"
#include "PlotUtils/CCInclusiveCuts.h"
#include "PlotUtils/CCInclusiveSignal.h"
#include "PlotUtils/Reweighter.h"
#include "PlotUtils/FluxAndCVReweighter.h"
#include "PlotUtils/GENIEReweighter.h"
#include "PlotUtils/LowRecoil2p2hReweighter.h"
#include "PlotUtils/RPAReweighter.h"
#include "PlotUtils/MINOSEfficiencyReweighter.h"
#include "PlotUtils/LowQ2PiReweighter.h"
#include "PlotUtils/AMUDISReweighter.h"
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

    int filledOrEmptyMEPlaylist(std::string playlist)
    {
        std::vector<std::string> filled = {"minervame1L", "minervame1M", "minervame1N", "minervame1O", "minervame1P"};
        std::vector<std::string> empty = {"minervame1A", "minervame1B", "minervame1C", "minervame1D", "minervame1E", "minervame1F", "minervame1G"/*,  "minervame5A", "minervame6A", "minervame6B", "minervame6C", "minervame6D", "minervame6E", "minervame6F", "minervame6G", "minervame6H", "minervame6I", "minervame6J" */};
        if (std::find(std::begin(filled), std::end(filled), playlist) != std::end(filled)) return 1;
        else if (std::find(std::begin(empty), std::end(empty), playlist) != std::end(empty)) return 2;
        //Is this the most efficient way? Probably not
        return 0;
    }

    //Taken from Oscar's code (/exp/minerva/app/users/omorenop/cmtuser/git-Mat/Personal/Test/InclusiveUtils.h)
    int planecode(int mdl, int plane)
    {

        if (mdl == -5) {
            if (plane == 1) {return 1;}
            if (plane == 2) {return 2;}}
        if (mdl == -4) {
            if (plane == 1) {return 3;}
            if (plane == 2) {return 4;}}
        if (mdl == -3) {
            if (plane == 1) {return 5;}
            if (plane == 2) {return 6;}}
        if (mdl == -2) {
            if (plane == 1) {return 7;}
            if (plane == 2) {return 8;}}
        if (mdl == -1) {
            if (plane == 1) {return 9;} //target 1
            if (plane == 2) {return 9;}}
        if (mdl == 0) {
            if (plane == 1) {return 10;}
            if (plane == 2) {return 11;}}
        if (mdl == 1) {
            if (plane == 1) {return 12;}
            if (plane == 2) {return 13;}}
        if (mdl == 2) {
            if (plane == 1) {return 14;}
            if (plane == 2) {return 15;}}
        if (mdl == 3) {
            if (plane == 1) {return 16;}
            if (plane == 2) {return 17;}}
        if (mdl == 4) {
            if (plane == 1) {return 18;}} //Target 2!}
        if (mdl == 5) {
            if (plane == 1) {return 19;}
            if (plane == 2) {return 20;}}
        if (mdl == 6) {
            if (plane == 1) {return 21;}
            if (plane == 2) {return 22;}
        }
        if (mdl == 7) {
            if (plane == 1) {return 23;}
            if (plane == 2) {return 24;}}
        if (mdl == 8) {
            if (plane == 1) {return 25;}
            if (plane == 2) {return 26;}}
        if (mdl == 9) {
            if (plane == 1) {return 27;} //target 3!
            if (plane == 2) {return 27;}}
        //water?
        //if (mdl == 10) {
        //    if (plane == 1) {
        //        return 1; //weird!!
        //    }
        //    if (plane == 2) {
        //        return 0; //weird!!
        //    }
        //}

        if (mdl == 11) {
            if (plane == 1) {return 28;}
            if (plane == 2) {return 29;}}
        if (mdl == 12) {
            if (plane == 1) {return 30;}
            if (plane == 2) {return 31;}}
        if (mdl == 13) {
            if (plane == 1) {return 32;}
            if (plane == 2) {return 33;}}
        if (mdl == 14) {
            if (plane == 1) {return 34;}
            if (plane == 2) {return 35;}}
        //Here goes Water!!!
        if (mdl == -999) {
            if(plane == -999) {return 36;}}
        if (mdl == 15) {
            if (plane == 1) {return 37;}
            if (plane == 2) {return 38;}}
        if (mdl == 16) {
            if (plane == 1) {return 39;}
            if (plane == 2) {return 40;}}
        if (mdl == 17) {
            if (plane == 1) {return 41;}
            if (plane == 2) {return 42;}}
        if (mdl == 18) {
            if (plane == 1) {return 43;}
            if (plane == 2) {return 44;}}
        if (mdl == 19) {
            if (plane == 1) {return 45; }//Target 4!!}
            if (plane == 2) {return 45;}} //Target 4!!}}
        if (mdl == 20) {
            if (plane == 1) {return 46;}
            if (plane == 2) {return 47;}}
        if (mdl == 21) {
            if (plane == 1) {return 48;}
            if (plane == 2) {return 49;}  }
        if (mdl == 22) {
            if (plane == 1) {return 50; }//target 5!!}
            if (plane == 2) {return 50;}} //target 5!!}  }
        if (mdl == 23) {
            if (plane == 1) {return 51;}
            if (plane == 2) {return 52;}}
        if (mdl == 24) {
            if (plane == 1) {return 53;}
            if (plane == 2) {return 54;}  }
        if (mdl == 25) {
            if (plane == 1) {return 55;}
            if (plane == 2) {return 56;}  }
        if (mdl == 26) {
            if (plane == 1) {return 57;}
            if (plane == 2) {return 58;}  }
        if (mdl == 27) {
            if (plane == 1) {return 59;}
            if (plane == 2) {return 60;}  }
        if (mdl == 28) {
            if (plane == 1) {return 61;}
            if (plane == 2) {return 62;}  }
        if (mdl == 29) {
            if (plane == 1) {return 63;}
            if (plane == 2) {return 64;}  }
        if (mdl == 30) {
            if (plane == 1) {return 65;}
            if (plane == 2) {return 66;}  }
        if (mdl == 31) {
            if (plane == 1) {return 67;}
            if (plane == 2) {return 68;}  }
        if (mdl == 32) {
            if (plane == 1) {return 69;}
            if (plane == 2) {return 70;}  }
        if (mdl == 33) {
            if (plane == 1) {return 71;}
            if (plane == 2) {return 72;}  }
        if (mdl == 34) {
            if (plane == 1) {return 73;}
            if (plane == 2) {return 74;}  }
        if (mdl == 35) {
            if (plane == 1) {return 75;}
            if (plane == 2) {return 76;}  }
        if (mdl == 36) {
            if (plane == 1) {return 77;}
            if (plane == 2) {return 78;}  }
        if (mdl == 37) {
            if (plane == 1) {return 79;}
            if (plane == 2) {return 80;}  }
        if (mdl == 38) {
            if (plane == 1) {return 81;}
            if (plane == 2) {return 82;}  }
        if (mdl == 39) {
            if (plane == 1) {return 83;}
            if (plane == 2) {return 84;}  }
        if (mdl == 40) {
            if (plane == 1) {return 85;}
            if (plane == 2) {return 86;}  }
        if (mdl == 41) {
            if (plane == 1) {return 87;}
            if (plane == 2) {return 88;}  }
        if (mdl == 42) {
            if (plane == 1) {return 89;}
            if (plane == 2) {return 90;}  }
        if (mdl == 43) {
            if (plane == 1) {return 91;}
            if (plane == 2) {return 92;}  }
        if (mdl == 44) {
            if (plane == 1) {return 93;}
            if (plane == 2) {return 94;}  }
        if (mdl == 45) {
            if (plane == 1) {return 95;}
            if (plane == 2) {return 96;}  }
        if (mdl == 46) {
            if (plane == 1) {return 97;}
            if (plane == 2) {return 98;}  }
        if (mdl == 47) {
            if (plane == 1) {return 99;}
            if (plane == 2) {return 100;}  }
        if (mdl == 48) {
            if (plane == 1) {return 101;}
            if (plane == 2) {return 102;}  }
        if (mdl == 49) {
            if (plane == 1) {return 103;}
            if (plane == 2) {return 104;}  }
        if (mdl == 50) {
            if (plane == 1) {return 105;}
            if (plane == 2) {return 106;}  }
        if (mdl == 51) {
            if (plane == 1) {return 107;}
            if (plane == 2) {return 108;}  }
        if (mdl == 52) {
            if (plane == 1) {return 109;}
            if (plane == 2) {return 110;}  }
        if (mdl == 53) {
            if (plane == 1) {return 111;}
            if (plane == 2) {return 112;}  }
        if (mdl == 54) {
            if (plane == 1) {return 113;}
            if (plane == 2) {return 114;}  }
        if (mdl == 55) {
            if (plane == 1) {return 115;}
            if (plane == 2) {return 116;}  }
        if (mdl == 56) {
            if (plane == 1) {return 117;}
            if (plane == 2) {return 118;}  }
        if (mdl == 57) {
            if (plane == 1) {return 119;}
            if (plane == 2) {return 120;}  }
        if (mdl == 58) {
            if (plane == 1) {return 121;}
            if (plane == 2) {return 122;}  }
        if (mdl == 59) {
            if (plane == 1) {return 123;}
            if (plane == 2) {return 124;}  }
        if (mdl == 60) {
            if (plane == 1) {return 125;}
            if (plane == 2) {return 126;}  }
        if (mdl == 61) {
            if (plane == 1) {return 127;}
            if (plane == 2) {return 128;}  }
        if (mdl == 62) {
            if (plane == 1) {return 129;}
            if (plane == 2) {return 130;}  }
        if (mdl == 63) {
            if (plane == 1) {return 131;}
            if (plane == 2) {return 132;}  }
        if (mdl == 64) {
            if (plane == 1) {return 133;}
            if (plane == 2) {return 134;}  }
        if (mdl == 65) {
            if (plane == 1) {return 135;}
            if (plane == 2) {return 136;}  }
        if (mdl == 66) {
            if (plane == 1) {return 1;}
            if (plane == 2) {return 138;}  }
        if (mdl == 67) {
            if (plane == 1) {return 139;}
            if (plane == 2) {return 140;}  }
        if (mdl == 68) {
            if (plane == 1) {return 141;}
            if (plane == 2) {return 142;}  }
        if (mdl == 69) {
            if (plane == 1) {return 143;}
            if (plane == 2) {return 144;}  }
        if (mdl == 70) {
            if (plane == 1) {return 145;}
            if (plane == 2) {return 146;}  }
        if (mdl == 71) {
            if (plane == 1) {return 147;}
            if (plane == 2) {return 148;}  }
        if (mdl == 72) {
            if (plane == 1) {return 149;}
            if (plane == 2) {return 150;}  }
        if (mdl == 73) {
            if (plane == 1) {return 151;}
            if (plane == 2) {return 152;}  }
        if (mdl == 74) {
            if (plane == 1) {return 153;}
            if (plane == 2) {return 154;}  }
        if (mdl == 75) {
            if (plane == 1) {return 155;}
            if (plane == 2) {return 156;}
        }
        if (mdl == 76) {
            if (plane == 1) {return 157;}
            if (plane == 2) {return 158;}  }
        if (mdl == 77) {
            if (plane == 1) {return 159;}
            if (plane == 2) {return 160;}  }
        if (mdl == 78) {
            if (plane == 1) {return 161;}
            if (plane == 2) {return 162;}  }
        if (mdl == 79) {
            if (plane == 1) {return 163;}
            if (plane == 2) {return 164;}  }
        if (mdl == 80) {
            if (plane == 1) {return 165;}
            if (plane == 2) {return 166;}  }
        if (mdl == 81) {
            if (plane == 1) {return 167;}
            if (plane == 2) {return 168;}  }
        if (mdl == 82) {
            if (plane == 1) {return 169;}
            if (plane == 2) {return 170;}  }
        if (mdl == 83) {
            if (plane == 1) {return 171;}
            if (plane == 2) {return 172;}  }
        if (mdl == 84) {
            if (plane == 1) {return 173;}
            if (plane == 2) {return 174;}  }
        if (mdl == 85) {
            if (plane == 1) {return 175;}
            if (plane == 2) {return 176;}  }
        if (mdl == 86) {
            if (plane == 1) {return 177;}
            if (plane == 2) {return 178;}  }
        if (mdl == 87) {
            if (plane == 1) {return 179;}
            if (plane == 2) {return 180;}  }
        if (mdl == 88) {
            if (plane == 1) {return 181;}
            if (plane == 2) {return 182;}  }
        if (mdl == 89) {
            if (plane == 1) {return 183;}
            if (plane == 2) {return 184;}  }
        if (mdl == 90) {
            if (plane == 1) {return 185;}
            if (plane == 2) {return 186;}  }
        if (mdl == 91) {
            if (plane == 1) {return 187;}
            if (plane == 2) {return 188;}  }
        if (mdl == 92) {
            if (plane == 1) {return 189;}
            if (plane == 2) {return 190;}  }
        if (mdl == 93) {
            if (plane == 1) {return 191;}
            if (plane == 2) {return 192;}  }
        if (mdl == 94) {
            if (plane == 1) {return 193;}
            if (plane == 2) {return 194;}  }
        if (mdl == 95) {
            if (plane == 1) {return 195;}
            if (plane == 2) {return 195;}  }
        if (mdl == 96) {
            if (plane == 1) {return 196;}
            if (plane == 2) {return 196;}  }
        if (mdl == 97) {
            if (plane == 1) {return 197;}
            if (plane == 2) {return 197;}  }
        if (mdl == 98) {
            if (plane == 1) {return 198;} //not needed}
            if (plane == 2) {return 198;}  }
        if (mdl == 99) {
            if (plane == 1) {return 199;} //not needed}
            if (plane == 2) {return 199;}  }
        if (mdl == 100) {
            if (plane == 2) {return 200;}  }
        if (mdl == 101) {
            if (plane == 2) {return 201;}  }
        if (mdl == 102) {
            if (plane == 2) {return 202;}  }
        if (mdl == 103) {
            if (plane == 2) {return 203;}  }
        if (mdl == 104) {
            if (plane == 2) {return 204;}  }
        if (mdl == 105) {
            if (plane == 2) {return 205;}  }
        if (mdl == 106) {
            if (plane == 2) {return 206;}  }
        if (mdl == 107) {
            if (plane == 2) {return 207;}  }
        if (mdl == 108) {
            if (plane == 2) {return 208;}  }
        if (mdl == 109) {
            if (plane == 2) {return 209;}  }
        if (mdl == 110) {
            if (plane == 2) {return 210;}  }
        if (mdl == 111) {
            if (plane == 2) {return 211;}  }
        if (mdl == 112) {
            if (plane == 2) {return 212;}  }
        if (mdl == 113) {
            if (plane == 2) {return 213;}  }
        if (mdl == 114) {
            if (plane == 2) {return 214;}  }
        if (mdl == 115) {
            if (plane == 2) {return 2;}  }
        if (mdl == 116) {
            if (plane == 2) {return 0;}  }
        if (mdl == 117) {
            if (plane == 2) {return 1;}  }
        if (mdl == 118) {
            if (plane == 2) {return 0;}  }
        if (mdl == 119) {
            if (plane == 2) {return 2;}  }

    }



    int getTargetCodeFromVtxInfo(double vtx_x, double vtx_y, double vtx_z, int mod, int plane)
    {
        PlotUtils::TargetUtils tgtUtil;
        bool distanceToDivCut = true;
        //tgtUtil.SetDistToDivCut( ??? ); //What was Anezka's + what is it in GenieXSecExtract/ Should this be tweaked?
        //2/Jul/2025 Answer to above question, Anezka used 25mm which is also the default so I will leave it at that
        if (tgtUtil.InIron1VolMC( vtx_x, vtx_y, vtx_z, 850., distanceToDivCut ) ) return 1026;
        else if (tgtUtil.InLead1VolMC( vtx_x, vtx_y, vtx_z, 850., distanceToDivCut ) ) return 1082;
        else if (tgtUtil.InIron2VolMC( vtx_x, vtx_y, vtx_z, 850., distanceToDivCut ) ) return 2026;
        else if (tgtUtil.InLead2VolMC( vtx_x, vtx_y, vtx_z, 850., distanceToDivCut ) ) return 2082;
        else if (tgtUtil.InCarbon3VolMC( vtx_x, vtx_y, vtx_z, 850., distanceToDivCut ) ) return 3006;
        else if (tgtUtil.InIron3VolMC( vtx_x, vtx_y, vtx_z, 850., distanceToDivCut ) ) return 3026;
        else if (tgtUtil.InLead3VolMC( vtx_x, vtx_y, vtx_z, 850., distanceToDivCut ) ) return 3082;
        else if (tgtUtil.InLead4VolMC( vtx_x, vtx_y, vtx_z, 850. ) ) return 4082;
        else if (tgtUtil.InIron5VolMC( vtx_x, vtx_y, vtx_z, 850., distanceToDivCut ) ) return 5026;
        else if (tgtUtil.InLead5VolMC( vtx_x, vtx_y, vtx_z, 850., distanceToDivCut ) ) return 5082;
        else if (tgtUtil.InWaterTargetVolMC( vtx_x, vtx_y, vtx_z, 850. ) ) return 6000;
        
        //Including the planes immediate up/downstream
        bool includingUSDSplanes = true; //true; //This is Anezka's fix
        //Utilising methods to get the z center position of each target to determine what material the x and y positions of interactions in the US/DS planes would correspond to, eg GetTarget1CenterZMC
        if (includingUSDSplanes)
        {
            //Upstream
            if (mod == -2 && plane == 2) //US of tgt 1
            {
            if (tgtUtil.InIron1VolMC( vtx_x, vtx_y, tgtUtil.GetTarget1CenterZMC(), 850., true )) return 1026;
            else if (tgtUtil.InLead1VolMC( vtx_x, vtx_y, tgtUtil.GetTarget1CenterZMC(), 850., true )) return 1082;
            }
            else if (mod == 3 && plane == 2) //US of tgt 2
            {
            if (tgtUtil.InIron2VolMC( vtx_x, vtx_y, tgtUtil.GetTarget2CenterZMC(), 850., true )) return 2026;
            else if (tgtUtil.InLead2VolMC( vtx_x, vtx_y, tgtUtil.GetTarget2CenterZMC(), 850., true )) return 2082;
            }
            else if (mod == 8 && plane == 2) //US of tgt 3
            {
            if (tgtUtil.InCarbon3VolMC( vtx_x, vtx_y, tgtUtil.GetTarget3CenterZMC(), 850., true )) return 3006;
            else if (tgtUtil.InIron3VolMC( vtx_x, vtx_y, tgtUtil.GetTarget3CenterZMC(), 850., true )) return 3026;
            else if (tgtUtil.InLead3VolMC( vtx_x, vtx_y, tgtUtil.GetTarget3CenterZMC(), 850., true )) return 3082;
            }
            else if (mod == 18 && plane == 2) return 4082; //US of tgt 4 //Should I still check the x and y vertex are within the target?
            else if (mod == 21 && plane == 2) //US of tgt 5
            {
            if (tgtUtil.InIron5VolMC( vtx_x, vtx_y, tgtUtil.GetTarget5CenterZMC(), 850., true )) return 5026;
            else if (tgtUtil.InLead5VolMC( vtx_x, vtx_y, tgtUtil.GetTarget5CenterZMC(), 850., true )) return 5082;
            }
            else if (mod == 14 && plane == 2) return 6000; //US of water target //Should I still check the x and y vertex are within the target?

            //Downstream
            if (mod == 0 && plane == 1) //US of tgt 1
            {
            if (tgtUtil.InIron1VolMC( vtx_x, vtx_y, tgtUtil.GetTarget1CenterZMC(), 850., true )) return 1026;
            else if (tgtUtil.InLead1VolMC( vtx_x, vtx_y, tgtUtil.GetTarget1CenterZMC(), 850., true )) return 1082;
            }
            else if (mod == 5 && plane == 1) //US of tgt 2
            {
            if (tgtUtil.InIron2VolMC( vtx_x, vtx_y, tgtUtil.GetTarget2CenterZMC(), 850., true )) return 2026;
            else if (tgtUtil.InLead2VolMC( vtx_x, vtx_y, tgtUtil.GetTarget2CenterZMC(), 850., true )) return 2082;
            }
            else if (mod == 11 && plane == 1) //US of tgt 3
            {
            if (tgtUtil.InCarbon3VolMC( vtx_x, vtx_y, tgtUtil.GetTarget3CenterZMC(), 850., true )) return 3006;
            else if (tgtUtil.InIron3VolMC( vtx_x, vtx_y, tgtUtil.GetTarget3CenterZMC(), 850., true )) return 3026;
            else if (tgtUtil.InLead3VolMC( vtx_x, vtx_y, tgtUtil.GetTarget3CenterZMC(), 850., true )) return 3082;
            }
            else if (mod == 20 && plane == 1) return 4082; //US of tgt 4 //Should I still check the x and y vertex are within the target?
            else if (mod == 23 && plane == 1) //US of tgt 5
            {
            if (tgtUtil.InIron5VolMC( vtx_x, vtx_y, tgtUtil.GetTarget5CenterZMC(), 850., true )) return 5026;
            else if (tgtUtil.InLead5VolMC( vtx_x, vtx_y, tgtUtil.GetTarget5CenterZMC(), 850., true )) return 5082;
            }
            else if (mod == 15 && plane == 1) return 6000; //US of water target //Should I still check the x and y vertex are within the target?
        }
        return -1;
    }

    template<class HIST>
    void AddHist(HIST& hist1, HIST* hist2, double scale = 1)
    {
        if (hist1.GetEntries() > 0) //If not null
        {
            hist1.Add(hist2, scale);
        }
        else //If null
        {
            hist1 = *hist2;
            hist1.Scale(scale);
        }
    }

    //Helper function used to find directories containing root files
    std::vector<std::string> findContainingDirectories(std::string dir, std::string type, bool recursiveSearch = true, bool skipTest = true, int curdepth = 0, int maxdepth = 2)
    {
        std::vector<std::string> dirpaths;
        for (const auto &entry : std::filesystem::directory_iterator(dir))
        {
            std::string path = entry.path();
            std::filesystem::file_type ft = std::filesystem::status(path).type();
            if (ft == std::filesystem::file_type::regular)
            {
                const size_t base = path.find("runEventLoop"+type);
                if ((base != std::string::npos) && (std::find(dirpaths.begin(), dirpaths.end(), dir) == dirpaths.end()))
                {
                    if (!(skipTest && (dir.find("/Test-") != std::string::npos))) dirpaths.push_back(dir); //If skipTest is false or if skipTest is true and the test string isn't found in the path, add to vector
                }
            }
            else if ((ft == std::filesystem::file_type::directory || ft == std::filesystem::file_type::symlink) && recursiveSearch && curdepth<=maxdepth)
            {
            std::vector<std::string> subdirs = findContainingDirectories(path, type, recursiveSearch, skipTest, curdepth+1);
            for (auto subdirpath : subdirs)
            {
                if (std::find(dirpaths.begin(), dirpaths.end(), subdirpath) == dirpaths.end()) dirpaths.push_back(subdirpath);
            }
            }
        }
        std::sort(dirpaths.begin(), dirpaths.end());
        return dirpaths;
        }

};

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
            double relativeErr = 1/univ.GetMuonQPErr(); //Need to do 1/err because of how MINOS reports significance
            if (univ.GetAnalysisNuPDG()>0) return ( relativeErr <= -fMin );
            else if (univ.GetAnalysisNuPDG()<0) return ( relativeErr >= fMin );
        } 

        const double fMin;
    };

    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class RockMuonCut: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
        public:
        RockMuonCut(): PlotUtils::Cut<UNIVERSE, EVENT>("Rock Muon Cut")
        {
        }

        private:
        bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
        {
            return (univ.GetInt("rock_muons_removed") == 1);
        }
    };

    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class VetoWall: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
        public:
        VetoWall(): PlotUtils::Cut<UNIVERSE, EVENT>("Veto Wall Cut")
        {
        }

        private:
        bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
        {
            return (univ.GetInt("VetoWall_event_IsVeto") == 1);
        }
    };
};

namespace util
{

  std::vector<double> PTBins = {0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1, 1.25, 1.5, 2.5, 4.5},
                      PzBins = {1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20, 40, 60},
                      EmuBins = {0,1,2,3,4,5,7,9,12,15,18,22,36,50,75,100,120},
                      bjorkenXbins = {0.001, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.6, 1 , 2.2},
                      Erecoilbins = {0, 0.2, 0.75, 1.5, 2.5, 3.5, 5, 8, 12, 14.5, 20}; //GeV
                      //Erecoilbins = {0, 200, 750, 1500, 2500, 3500, 5000, 8000, 12000, 14500, 20000}; //MeV


  //Boiler plate precuts for CC inclusive NuMu ME Analysis
  const double apothem = 850; //All in mm
  PlotUtils::Cutter<CVUniverse, MichelEvent>::reco_t GetAnalysisCuts(int pdg)
  {
    PlotUtils::Cutter<CVUniverse, MichelEvent>::reco_t cuts;
    cuts.emplace_back(new reco::Apothem<CVUniverse, MichelEvent>(apothem));
    cuts.emplace_back(new reco::MaxMuonAngle<CVUniverse, MichelEvent>(17.));
    cuts.emplace_back(new reco::HasMINOSMatch<CVUniverse, MichelEvent>());
    cuts.emplace_back(new reco::NoDeadtime<CVUniverse, MichelEvent>(1, "Deadtime"));
    if (pdg>0)  cuts.emplace_back(new reco::IsNeutrino<CVUniverse, MichelEvent>()); //Used minos curvature
    else if (pdg<0)  cuts.emplace_back(new reco::IsAntiNeutrino<CVUniverse, MichelEvent>()); //Used minos curvature
    cuts.emplace_back(new reco::MuonCurveSignificance<CVUniverse, MichelEvent>(5));
    cuts.emplace_back(new reco::MuonEnergyMin<CVUniverse, MichelEvent>(2000.0, "EMu Min"));
    cuts.emplace_back(new reco::MuonEnergyMax<CVUniverse, MichelEvent>(20000.0, "EMu Max"));
    cuts.emplace_back(new reco::ANNConfidenceCut<CVUniverse, MichelEvent>(0.40)); //Reccomended at 0.4 for P6 ML vertexing and 0.2 for P4 vertexing
    //cuts.emplace_back(new reco::RockMuonCut<CVUniverse, MichelEvent>()); //Reccomended for P6 ML vertexing
    //cuts.emplace_back(new reco::VetoWall<CVUniverse, MichelEvent>()); //Reccomended for P6 ML vertexing
    return cuts;
  }

  PlotUtils::Cutter<CVUniverse, MichelEvent>::truth_t GetPhaseSpace()
  {
    PlotUtils::Cutter<CVUniverse, MichelEvent>::truth_t phasespace;
    phasespace.emplace_back(new truth::Apothem<CVUniverse>(apothem));
    phasespace.emplace_back(new truth::MuonAngle<CVUniverse>(17.));
    phasespace.emplace_back(new truth::MuonEnergyMin<CVUniverse>(2000.0, "EMu Min"));
    phasespace.emplace_back(new truth::MuonEnergyMax<CVUniverse>(20000.0, "EMu Max"));
    return phasespace;
  }
};
#endif //UTIL_NUKEUTILS_H