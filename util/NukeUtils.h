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
#include "cuts/CCInclCuts.h"
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
    std::map<int, std::string> SidebandCategoriesExtWater = {{0, "US"}, {1, "DS"}, {2, "Signal"}, {3, "WaterTank"}};
    std::map<int, std::string> TgtCodeLabelsTracker = {{-1, "Tracker"}};
    std::map<int, std::string> TgtCodeLabelsNuke = {{1026, "1026"}, {1082, "1082"}, {2026, "2026"}, {2082, "2082"}, {3006, "3006"}, {3026, "3026"}, {3082, "3082"}, {4082, "4082"}, {5026, "5026"}, {5082, "5082"}, {6000, "Water"}, {7, "Target7"}, {8, "Target8"}, {9, "Target9"}, {10, "Target10"}, {11, "Target11"}, {12, "Target12"}};

    //Without pseudotargets
    std::map<int, std::string> TgtCodeLabelsNukeNoPS = {{1026, "1026"}, {1082, "1082"}, {2026, "2026"}, {2082, "2082"}, {3006, "3006"}, {3026, "3026"}, {3082, "3082"}, {4082, "4082"}, {5026, "5026"}, {5082, "5082"}, {6000, "Water"}};

    std::map<int, std::string> GENIELabels = {{1, "QE"},
                                              {8, "2p2h"},
                                              {2, "RES"},
                                              {3, "DIS"}}; //kNoInteraction=0, kQEL=1, kRES=2, kDIS=3, kCOHPI=4, kAMNUGAMMA=5, kIMD=6, kNUEEL=7, k2P2H=8

    std::map<int, std::string> BKGLabels = {{0, "NC_Bkg"},
					       {1, "Wrong_Sign_Bkg"}};

    std::map<int, std::string> BKGLabelsWithPlasticSidebands = {{0, "NC_Bkg"},
					       {1, "Wrong_Sign_Bkg"}, {2, "Upstream_Plastic_Bkg"}, {3, "Downstream_Plastic_Bkg"}, {4, "Water_Tank_Bkg"}, {5, "True_In_Other_Target_Bkg"}, {6, "True_Vtx_Elsewhere_Bkg"}};

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
        std::vector<std::string> filled = {"minervame1L", "minervame1M", "minervame1N", "minervame1O", "minervame1P", "1L", "1M", "1N", "1O", "1P"};
        std::vector<std::string> empty = {"minervame1A", "minervame1B", "minervame1C", "minervame1D", "minervame1E", "minervame1F", "minervame1G", "1A", "1B", "1C", "1D", "1E", "1F", "1G"/*,  "minervame5A", "minervame6A", "minervame6B", "minervame6C", "minervame6D", "minervame6E", "minervame6F", "minervame6G", "minervame6H", "minervame6I", "minervame6J" */};
        if (std::find(std::begin(filled), std::end(filled), playlist) != std::end(filled)) return 1;
        else if (std::find(std::begin(empty), std::end(empty), playlist) != std::end(empty)) return 2;
        //Is this the most efficient way? Probably not
        return 0;
    }


    double getZPosFromSegment(int segment)
    {
        if (segment ==1) return 4293.04;
        if (segment == 2) return 4313.68;
        if (segment == 3) return 4337.25;
        if (segment == 4) return 4357.9;
        if (segment == 5) return 4381.47;
        if (segment == 6) return 4402.11;
        if (segment == 7) return 4425.68;
        if (segment == 8) return 4446.33;
        if (segment == 9) return 4481.21; //Target 1
        if (segment == 10) return 4514.11;
        if (segment == 11) return 4534.76;
        if (segment == 12) return 4558.33;
        if (segment == 13) return 4578.97;
        if (segment == 14) return 4602.54;
        if (segment == 15) return 4623.19;
        if (segment == 16) return 4646.76;
        if (segment == 17) return 4667.4;
        if (segment == 18) return 4702.29; //Target 2
        if (segment == 19) return 4735.19;
        if (segment == 20) return 4755.83;
        if (segment == 21) return 4779.4;
        if (segment == 22) return 4800.05;
        if (segment == 23) return 4823.62;
        if (segment == 24) return 4844.26;
        if (segment == 25) return 4867.83;
        if (segment == 26) return 4888.48;
        if (segment == 27) return 4923.36; //Target 3
        if (segment == 28) return 5000.48;
        if (segment == 29) return 5021.12;
        if (segment == 30) return 5044.69;
        if (segment == 31) return 5065.34;
        if (segment == 32) return 5088.91;
        if (segment == 33) return 5109.55;
        if (segment == 34) return 5133.12;
        if (segment == 35) return 5153.77;
        if (segment == 36) return 5310; //Water Target -- Find better number, this is a guesstimate
        if (segment == 37) return 5456.74;
        if (segment == 38) return 5477.38;
        if (segment == 39) return 5500.95;
        if (segment == 40) return 5521.6;
        if (segment == 41) return 5545.17;
        if (segment == 42) return 5565.81;
        if (segment == 43) return 5589.38;
        if (segment == 44) return 5610.02;
        if (segment == 45) return 5644.91; //Target 4
        if (segment == 46) return 5677.81;
        if (segment == 47) return 5698.45;
        if (segment == 48) return 5722.03;
        if (segment == 49) return 5742.67;
        if (segment == 50) return 5777.55; //Target 5
        if (segment == 51) return 5810.45;
        if (segment == 52) return 5831.1;
        if (segment == 53) return 5855.68;
        if (segment == 54) return 5876.33;
        if (segment == 55) return 5900.91;
        if (segment == 56) return 5921.56;
        if (segment == 57) return 5946.14;
        if (segment == 58) return 5966.79;
        if (segment == 59) return 5991.37;
        if (segment == 60) return 6012.01;
        if (segment == 61) return 6036.6;
        if (segment == 62) return 6057.24;
        if (segment == 63) return 6081.83;
        if (segment == 64) return 6102.47;
        if (segment == 65) return 6127.06;
        if (segment == 66) return 6147.7;
        if (segment == 67) return 6172.29;
        if (segment == 68) return 6192.93;
        if (segment == 69) return 6217.52;
        if (segment == 70) return 6238.16;
        if (segment == 71) return 6262.74;
        if (segment == 72) return 6283.39;
        if (segment == 73) return 6307.97;
        if (segment == 74) return 6328.62;
        if (segment == 75) return 6353.2;
        if (segment == 76) return 6373.85;
        if (segment == 77) return 6398.43;
        if (segment == 78) return 6419.08;
        if (segment == 79) return 6443.66;
        if (segment == 80) return 6464.3;
        if (segment == 81) return 6488.89;
        if (segment == 82) return 6509.53;
        if (segment == 83) return 6534.12;
        if (segment == 84) return 6554.76;
        if (segment == 85) return 6579.35;
        if (segment == 86) return 6599.99;
        if (segment == 87) return 6624.58;
        if (segment == 88) return 6645.22;
        if (segment == 89) return 6669.81;
        if (segment == 90) return 6690.45;
        if (segment == 91) return 6715.03;
        if (segment == 92) return 6735.68;
        if (segment == 93) return 6760.26;
        if (segment == 94) return 6780.91;
        if (segment == 95) return 6805.49;
        if (segment == 96) return 6826.14;
        if (segment == 97) return 6850.72;
        if (segment == 98) return 6871.37;
        if (segment == 99) return 6895.95;
        if (segment == 100) return 6916.59;
        if (segment == 101) return 6941.18;
        if (segment == 102) return 6961.82;
        if (segment == 103) return 6986.41;
        if (segment == 104) return 7007.05;
        if (segment == 105) return 7031.64;
        if (segment == 106) return 7052.28;
        if (segment == 107) return 7076.87;
        if (segment == 108) return 7097.51;
        if (segment == 109) return 7122.1;
        if (segment == 110) return 7142.74;
        if (segment == 111) return 7167.32;
        if (segment == 112) return 7187.97;
        if (segment == 113) return 7212.55;
        if (segment == 114) return 7233.2;
        if (segment == 115) return 7257.78;
        if (segment == 116) return 7278.43;
        if (segment == 117) return 7303.01;
        if (segment == 118) return 7323.66;
        if (segment == 119) return 7348.24;
        if (segment == 120) return 7368.88;
        if (segment == 121) return 7393.47;
        if (segment == 122) return 7414.11;
        if (segment == 123) return 7438.7;
        if (segment == 124) return 7459.34;
        if (segment == 125) return 7483.93;
        if (segment == 126) return 7504.57;
        if (segment == 127) return 7529.16;
        if (segment == 128) return 7549.8;
        if (segment == 129) return 7574.39;
        if (segment == 130) return 7595.03;
        if (segment == 131) return 7619.61;
        if (segment == 132) return 7640.26;
        if (segment == 133) return 7664.84;
        if (segment == 134) return 7685.49;
        if (segment == 135) return 7710.07;
        if (segment == 136) return 7730.72;
        if (segment == 137) return 7755.3;
        if (segment == 138) return 7775.95;
        if (segment == 139) return 7800.53;
        if (segment == 140) return 7821.17;
        if (segment == 141) return 7845.76;
        if (segment == 142) return 7866.4;
        if (segment == 143) return 7890.99;
        if (segment == 144) return 7911.63;
        if (segment == 145) return 7936.22;
        if (segment == 146) return 7956.86;
        if (segment == 147) return 7981.45;
        if (segment == 148) return 8002.09;
        if (segment == 149) return 8026.68;
        if (segment == 150) return 8047.32;
        if (segment == 151) return 8071.9;
        if (segment == 152) return 8092.55;
        if (segment == 153) return 8117.13;
        if (segment == 154) return 8137.78;
        if (segment == 155) return 8162.36;
        if (segment == 156) return 8183.01;
        if (segment == 157) return 8207.59;
        if (segment == 158) return 8228.24;
        if (segment == 159) return 8252.82;
        if (segment == 160) return 8273.46;
        if (segment == 161) return 8298.05;
        if (segment == 162) return 8318.69;
        if (segment == 163) return 8343.28;
        if (segment == 164) return 8363.92;
        if (segment == 165) return 8388.51;
        if (segment == 166) return 8409.15;
        if (segment == 167) return 8433.74;
        if (segment == 168) return 8454.38;
        if (segment == 169) return 8478.97;
        if (segment == 170) return 8499.61;
        if (segment == 171) return 8524.19;
        if (segment == 172) return 8544.84;
        if (segment == 173) return 8569.42;
        if (segment == 174) return 8590.07;
        if (segment == 175) return 8614.65;
        if (segment == 176) return 8635.3;
        if (segment == 177) return 8659.46;
        if (segment == 178) return 8680.1;
        if (segment == 179) return 8704.26;
        if (segment == 180) return 8724.9;
        if (segment == 181) return 8749.06;
        if (segment == 182) return 8769.71;
        if (segment == 183) return 8793.86;
        if (segment == 184) return 8814.51;
        if (segment == 185) return 8838.67;
        if (segment == 186) return 8859.31;
        if (segment == 187) return 8883.47;
        if (segment == 188) return 8904.11;
        if (segment == 189) return 8928.27;
        if (segment == 190) return 8948.92;
        if (segment == 191) return 8973.08;
        if (segment == 192) return 8993.72;
        if (segment == 193) return 9017.88;
        if (segment == 194) return 9038.52;
        if (segment == 195) return 9088.08;
        if (segment == 196) return 9135.41;
        if (segment == 197) return 9182.75;
        if (segment == 198) return 9230.08;
        if (segment == 199) return 9277.41;
        if (segment == 200) return 9324.74;
        if (segment == 201) return 9372.08;
        if (segment == 202) return 9419.41;
        if (segment == 203) return 9466.74;
        if (segment == 204) return 9514.07;
        if (segment == 205) return 9561.41;
        if (segment == 206) return 9608.74;
        if (segment == 207) return 9656.07;
        if (segment == 208) return 9703.4;
        if (segment == 209) return 9750.74;
        if (segment == 210) return 9798.07;
        if (segment == 211) return 9845.4;
        if (segment == 212) return 9892.73;
        if (segment == 213) return 9940.07;
        if (segment == 214) return 9987.4;
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

    int getExtendedTarget(int mod, int plane, double vtx_x, double vtx_y) //Get which target this module and plane corresponds to if using extended target definition
    {
        PlotUtils::TargetUtils tgtUtil;
        bool distanceToDivCut = true;
        //Upstream
        if ((mod == -2 && plane == 2) || (mod == 0 && plane == 1)) //US/DS of tgt 1
        {
        if (tgtUtil.InIron1VolMC( vtx_x, vtx_y, tgtUtil.GetTarget1CenterZMC(), 850., distanceToDivCut )) return 1026;
        else if (tgtUtil.InLead1VolMC( vtx_x, vtx_y, tgtUtil.GetTarget1CenterZMC(), 850., distanceToDivCut )) return 1082;
        }
        else if ((mod == 3 && plane == 2) || (mod == 5 && plane == 1)) //US/DS of tgt 2
        {
        if (tgtUtil.InIron2VolMC( vtx_x, vtx_y, tgtUtil.GetTarget2CenterZMC(), 850., distanceToDivCut )) return 2026;
        else if (tgtUtil.InLead2VolMC( vtx_x, vtx_y, tgtUtil.GetTarget2CenterZMC(), 850., distanceToDivCut )) return 2082;
        }
        else if ((mod == 8 && plane == 2) || (mod == 11 && plane == 1)) //US/DS of tgt 3
        {
        if (tgtUtil.InCarbon3VolMC( vtx_x, vtx_y, tgtUtil.GetTarget3CenterZMC(), 850., distanceToDivCut )) return 3006;
        else if (tgtUtil.InIron3VolMC( vtx_x, vtx_y, tgtUtil.GetTarget3CenterZMC(), 850., distanceToDivCut )) return 3026;
        else if (tgtUtil.InLead3VolMC( vtx_x, vtx_y, tgtUtil.GetTarget3CenterZMC(), 850., distanceToDivCut )) return 3082;
        }
        else if ((mod == 18 && plane == 2) || (mod == 20 && plane == 1)) //US/DS of tgt 4 //still check the x and y vertex are within the target
        {
            if (tgtUtil.InLead4VolMC( vtx_x, vtx_y, tgtUtil.GetTarget4CenterZMC(), 850. )) return 4082;
        }
        else if ((mod == 21 && plane == 2) || (mod == 23 && plane == 1)) //US/DS of tgt 5
        {
        if (tgtUtil.InIron5VolMC( vtx_x, vtx_y, tgtUtil.GetTarget5CenterZMC(), 850., distanceToDivCut )) return 5026;
        else if (tgtUtil.InLead5VolMC( vtx_x, vtx_y, tgtUtil.GetTarget5CenterZMC(), 850., distanceToDivCut )) return 5082;
        }
        else if ((mod == 14 && plane == 2) || (mod == 15 && plane == 1)) //US/DS of water target //Still check the x and y vertex are within the target?
        {
            if (tgtUtil.InWaterTargetVolMC( vtx_x, vtx_y, (TargetProp::WaterTarget::Face+TargetProp::WaterTarget::Back)/2, 850. )) return 6000;
        }
        return -1;
    }

    int getTargetCodeFromVtxInfo(double vtx_x, double vtx_y, double vtx_z, int mod, int plane, bool extendedTargetDefinition = true)
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
        if (extendedTargetDefinition) return getExtendedTarget(mod, plane, vtx_x, vtx_y); //Get which target this module and plane corresponds to if using extended target definition
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
            hist1 = *(hist2->Clone());
            hist1.Scale(scale);
        }
    }

    /* template<class HIST>
    void AddHist(HIST* hist1, HIST* hist2, double scale = 1)
    {
        if (hist1->GetEntries() > 0) //If not null
        {
            hist1->Add(hist2, scale);
        }
        else //If null
        {
            hist1 = hist2;
            hist1->Scale(scale);
        }
    } */

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
                    if (!(skipTest && (dir.find("/Test") != std::string::npos))) dirpaths.push_back(dir); //If skipTest is false or if skipTest is true and the test string isn't found in the path, add to vector
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

    // Treat the plastic between planes as though they were nuclear targets themselves
    //
    int getPlasticPseudoTargetCode(int mod, int plane, bool extTarget = true)
    {
    if (extTarget) //If using extended target definition, exclude planes that would be inlcuded
    {
        if (mod == -2 && plane == 2)
        return -1;
        if (mod == 0 && plane == 1)
        return -1;
        if (mod == 3 && plane == 2)
        return -1;
        if (mod == 5 && plane == 1)
        return -1;
        if (mod == 8 && plane == 2)
        return -1;
        if (mod == 11 && plane == 1)
        return -1;
        if (mod == 14 && plane == 2)
        return -1;
        if (mod == 15 && plane == 1)
        return -1;
        if (mod == 18 && plane == 2)
        return -1;
        if (mod == 20 && plane == 1)
        return -1;
        if (mod == 21 && plane == 2)
        return -1;
    }
    if (mod >= -5 && mod <= -2)
        return 7;
    else if (mod >= 0 && mod <= 3)
        return 8;
    else if (mod >= 5 && mod <= 8)
        return 9;
    else if (mod >= 11 && mod <= 14)
        return 10;
    else if (mod >= 15 && mod <= 18)
        return 11;
    else if (mod >= 20 && mod <= 21)
        return 12;

    //Segmenting the tracker region for the sake of CH-CH comparisons
    if (mod>=23 && mod <27) return 13;
    else if (mod>=27 && mod <81)
    {   
        //int num = 13+std::floor((mod - 23)/3);
        int num = 14+std::floor((mod - 24)/12);
        //std::cout<<"NUM: "<<num<<std::endl;
        return (num);
    }
    else if (mod>=81 && mod <=84) return 19;
    else
        return -1;
    }

    int getTgtCode(const CVUniverse *universe, bool truth, bool useExtendedTarget )
    {
    double vtx_x, vtx_y, vtx_z;
    int mod, plane;
    if (truth) // Truth
    {
        //First step, what target does the tuple say it's in
        int truthTgtCode = universe->GetTruthTargetCode();
        if (truthTgtCode > 0) return truthTgtCode; //Maybe do a more explicit check for if it's a valid target code?

        mod = universe->GetTruthVtxModule();
        plane = universe->GetTruthVtxPlane();
        ROOT::Math::XYZTVector Vtx = universe->GetTrueVertex();
        vtx_x = Vtx.X();
        vtx_y = Vtx.Y();
        vtx_z = Vtx.Z();
        //Water is a special case, due to (I think?) a bug in the MAT tuples. To figure out if an event happened in water we cant use the tuple branch, instead compare it's vertex position
        //if (universe->GetTruthVtxModule() == -999 && universe->GetTruthVtxPlane() == -999  && (mcTgtZ = 6 || mcTgtZ == 1) ) return 6000;
        //if (universe->GetTruthVtxModule() == 15 && universe->GetTruthVtxPlane() == 1  && (mcTgtZ = 6 || mcTgtZ == 1) ) return 6000;
        PlotUtils::TargetUtils tgtUtil;
        if (tgtUtil.InWaterTargetVolMC( vtx_x, vtx_y, vtx_z, 850. )) return 6000;
    }
    else // ANN
    {
        //First step, what target does the tuple say it's in
        int ANNTgtCode = universe->GetANNTargetCode();
        if (ANNTgtCode > 0) return ANNTgtCode; //Maybe do a more explicit check for if it's a valid target code?

        //Water is a special case, due to (I think?) a bug in the MAT tuples. To figure out if an event happened in water we cant use the regular tuple branch, instead we check if the ANN reconstruction put in the water segment (36)
        //if (universe->GetANNSegment() == 36 ) return 6000;

        mod = universe->GetANNVtxModule();
        plane = universe->GetANNVtxPlane();
        ROOT::Math::XYZVector Vtx = universe->GetANNVertex();
        vtx_x = Vtx.X();
        vtx_y = Vtx.Y();
        vtx_z = Vtx.Z();
        PlotUtils::TargetUtils tgtUtil;
        if (tgtUtil.InWaterTargetVolMC( vtx_x, vtx_y, vtx_z, 850. )) return 6000;
    }
    //Next determine if the interaction vertex is within the extended target definition (if applicable)
    if (useExtendedTarget && mod <24)
    {
        int extendedTarget = util::getExtendedTarget(mod, plane, vtx_x, vtx_y);
        if (extendedTarget!=-1) return extendedTarget;
    }
    //Lastly, check if it's in a pseudotarget
    return getPlasticPseudoTargetCode(mod, plane); 

    //Return value for function will be -1 if no target or pseudotarget is determined for this event
    }


    bool isTargetSideband(CVUniverse *universe, int mode /*mc = 0, ANN = 1, TB = 2*/, int targetCode, int USorDS /*Some planes are both the US planes of one target but the DS planes of another and so need to be treated twice*/, bool removeNeighbors = true)
    {
        PlotUtils::TargetUtils tgtUtil;
        double vtx_x, vtx_y, vtx_z;
        int mod, plane;
        if (mode == 0) // Truth
        {
            mod = universe->GetTruthVtxModule();
            plane = universe->GetTruthVtxPlane();
            ROOT::Math::XYZTVector Vtx = universe->GetTrueVertex();
            vtx_x = Vtx.X();
            vtx_y = Vtx.Y();
            vtx_z = Vtx.Z();
        }
        else if (mode == 1) // ANN
        {
            mod = universe->GetANNVtxModule();
            plane = universe->GetANNVtxPlane();
            ROOT::Math::XYZVector Vtx = universe->GetANNVertex();
            vtx_x = Vtx.X();
            vtx_y = Vtx.Y();
            vtx_z = Vtx.Z();
        }
        else if (mode == 2) // TB
        {
            mod = universe->GetMADVtxModule();
            plane = universe->GetMADVtxPlane();
            ROOT::Math::XYZTVector Vtx = universe->GetVertex();
            vtx_x = Vtx.X();
            vtx_y = Vtx.Y();
            vtx_z = Vtx.Z();
        }

        //Sideband for psuedotargets - experimental
        if (targetCode < 1000) return false; //No sideband for tracker pseudotargets
        if (targetCode<12)
        {
            if (targetCode==7 && !USorDS) return false;
            if (targetCode==7 && USorDS && mod == -2 && plane == 2) return true;
            if (targetCode==8 && USorDS && mod == 3 && plane == 2) return true;
            if (targetCode==8 && !USorDS && mod == 0 && plane == 1) return true;
            if (targetCode==9 && USorDS && mod == 8 && plane == 2) return true;
            if (targetCode==9 && !USorDS && mod == 5 && plane == 1) return true;
            if (targetCode==10 && USorDS && mod == 14 && plane == 2) return true;
            if (targetCode==10 && !USorDS && mod == 11 && plane == 1) return true;
            if (targetCode==11 && USorDS && mod == 18 && plane == 2) return true;
            if (targetCode==11 && !USorDS && mod == 15 && plane == 1) return true;
            if (targetCode==12 && USorDS && mod == 21 && plane == 2) return true;
            if (targetCode==12 && !USorDS && mod == 20 && plane == 1) return true;
        }

        // Removing planes immediately up and downstream of targets if we're looking at a nuclear target
        if (removeNeighbors && targetCode > 1000)
        {
            if (mod == -2 && plane == 2)
            return false;
            if (mod == 0 && plane == 1)
            return false;
            if (mod == 3 && plane == 2)
            return false;
            if (mod == 5 && plane == 1)
            return false;
            if (mod == 8 && plane == 2)
            return false;
            if (mod == 11 && plane == 1)
            return false;
            if (mod == 14 && plane == 2)
            return false;
            if (mod == 15 && plane == 1)
            return false;
            if (mod == 18 && plane == 2)
            return false;
            if (mod == 20 && plane == 1)
            return false;
            if (mod == 21 && plane == 2)
            return false;
        }
        if (mod >= 0 && mod <= 3) // DS of target 1 (Iron and Lead) and US of target 2 (Iron and Lead)
        {
            if (USorDS) // If true we're looking for US Planes therefore we're looking at US of tgt2
            {
                if (targetCode == 2026 && tgtUtil.InIron2VolMC(vtx_x, vtx_y, tgtUtil.GetTarget2CenterZMC(), 850., true)) return true;
                else if (targetCode == 2082 && tgtUtil.InLead2VolMC(vtx_x, vtx_y, tgtUtil.GetTarget2CenterZMC(), 850., true)) return true;
            }
            else // If false we're looking for DS Planes therefore we're looking at DS of tgt1
            {
                if (targetCode == 1026 && tgtUtil.InIron1VolMC(vtx_x, vtx_y, tgtUtil.GetTarget1CenterZMC(), 850., true)) return true;
                else if (targetCode == 1082 && tgtUtil.InLead1VolMC(vtx_x, vtx_y, tgtUtil.GetTarget1CenterZMC(), 850., true)) return true;
            }
        }
        if (mod >= 5 && mod <= 8) // DS of target 2 (Iron and Lead) and US of target 3 (Carbon, iron and Lead)
        {
            if (USorDS) // If true we're looking for US Planes therefore we're looking at US of tgt3
            {
                if (targetCode == 3006 && tgtUtil.InCarbon3VolMC(vtx_x, vtx_y, tgtUtil.GetTarget3CenterZMC(), 850., true)) return true;
                else if (targetCode == 3026 && tgtUtil.InIron3VolMC(vtx_x, vtx_y, tgtUtil.GetTarget3CenterZMC(), 850., true)) return true;
                else if (targetCode == 3082 && tgtUtil.InLead3VolMC(vtx_x, vtx_y, tgtUtil.GetTarget3CenterZMC(), 850., true)) return true;
            }
            else // If false we're looking for DS Planes therefore we're looking at DS of tgt2
            {
                if (targetCode == 2026 && tgtUtil.InIron2VolMC(vtx_x, vtx_y, tgtUtil.GetTarget2CenterZMC(), 850., true)) return true;
                else if (targetCode == 2082 && tgtUtil.InLead2VolMC(vtx_x, vtx_y, tgtUtil.GetTarget2CenterZMC(), 850., true)) return true;
            }
        }
        if (mod >= 11 && mod <= 14) // DS of target 3 (Carbon, iron and Lead) and US of water target
        {
            if (USorDS) // If true we're looking for US Planes therefore we're looking at US of water
            {
                if (targetCode == 6000 && tgtUtil.InWaterTargetVolMC(vtx_x, vtx_y, (PlotUtils::TargetProp::WaterTarget::Face + PlotUtils::TargetProp::WaterTarget::Back) / 2, 850.)) return true;
            }
            else // If false we're looking for DS Planes therefore we're looking at DS of tgt3
            {
                if (targetCode == 3006 && tgtUtil.InCarbon3VolMC(vtx_x, vtx_y, tgtUtil.GetTarget3CenterZMC(), 850., true)) return true;
                else if (targetCode == 3026 && tgtUtil.InIron3VolMC(vtx_x, vtx_y, tgtUtil.GetTarget3CenterZMC(), 850., true)) return true;
                else if (targetCode == 3082 && tgtUtil.InLead3VolMC(vtx_x, vtx_y, tgtUtil.GetTarget3CenterZMC(), 850., true)) return true;
            }
        }
        if (mod >= 15 && mod <= 18) // DS of water target and US of target 4 ( Lead)
        {
            if (USorDS) // If true we're looking for US Planes therefore we're looking at US of target 4
            {
                if (targetCode == 4082 && tgtUtil.InLead4VolMC(vtx_x, vtx_y, tgtUtil.GetTarget4CenterZMC(), 850.)) return true;
            }
            else // If false we're looking for DS Planes therefore we're looking at DS of water target
            {
                if (targetCode == 6000 && tgtUtil.InWaterTargetVolMC(vtx_x, vtx_y, (PlotUtils::TargetProp::WaterTarget::Face + PlotUtils::TargetProp::WaterTarget::Back) / 2, 850.)) return true;
            }
        }
        if (mod >= 20 && mod <= 21) // DS of target 4 (Lead) and US of target 5 (Iron and Lead)
        {
            if (USorDS) // If true we're looking for US Planes therefore we're looking at US of tgt3
            {
                if (targetCode == 5026 && tgtUtil.InIron5VolMC(vtx_x, vtx_y, tgtUtil.GetTarget5CenterZMC(), 850., true)) return true;
                else if (targetCode == 5082 && tgtUtil.InLead5VolMC(vtx_x, vtx_y, tgtUtil.GetTarget5CenterZMC(), 850., true)) return true;
            }
            else // If false we're looking for DS Planes therefore we're looking at DS of tgt4
            {
                if (targetCode == 4082 && tgtUtil.InLead4VolMC(vtx_x, vtx_y, tgtUtil.GetTarget4CenterZMC(), 850.)) return true;
            }
        }
        if (mod >= 20 && mod <= 21) // DS of target 4 (Lead) and US of target 5 (Iron and Lead)
        {
            if (USorDS) // If true we're looking for US Planes therefore we're looking at US of tgt3
            {
                if (targetCode == 5026 && tgtUtil.InIron5VolMC(vtx_x, vtx_y, tgtUtil.GetTarget5CenterZMC(), 850., true)) return true;
                else if (targetCode == 5082 && tgtUtil.InLead5VolMC(vtx_x, vtx_y, tgtUtil.GetTarget5CenterZMC(), 850., true)) return true;
            }
            else // If false we're looking for DS Planes therefore we're looking at DS of tgt4
            {
                if (targetCode == 4082 && tgtUtil.InLead4VolMC(vtx_x, vtx_y, tgtUtil.GetTarget4CenterZMC(), 850.)) return true;
            }
        }
        if (mod >= 23 && mod <= 26) // DS of target 5 (Iron and Lead)
        {
            if (USorDS) // If true we're looking for US Planes therefore we're looking at US of tgt3
            {
                return false; // This is immedaitely after target 5 and before the tracker region, so these planes are not the upstream region for any targets
            }
            else // If false we're looking for DS Planes therefore we're looking at DS of tgt4
            {
                if (targetCode == 5026 && tgtUtil.InIron5VolMC(vtx_x, vtx_y, tgtUtil.GetTarget5CenterZMC(), 850., true)) return true;
                else if (targetCode == 5082 && tgtUtil.InLead5VolMC(vtx_x, vtx_y, tgtUtil.GetTarget5CenterZMC(), 850., true)) return true;
            }
        }
        return false;
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
            return (univ.GetInt("rock_muons_removed") != 1);
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
            return (univ.GetInt("VetoWall_event_IsVeto") != 1);
        }
    };

    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class IsInTarget: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
        public:
        IsInTarget(int targetCode, bool useExtendedTarget): PlotUtils::Cut<UNIVERSE, EVENT>("Is In Target " + std::to_string(targetCode)), fTargetCode(targetCode), fUseExtendedTarget(useExtendedTarget)
        {
        }

        private:
        bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
        {   
            int tgt = util::getTgtCode(&univ, false, fUseExtendedTarget);
            //std::cout<< "Target code: " << fTargetCode << std::endl;
            //std::cout<< "Target code1: " << tgt << std::endl;
            return tgt == fTargetCode;
        }
        int fTargetCode;
        bool fUseExtendedTarget;
    };

    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    using ANNMuonEnergyMinGeV = PlotUtils::Minimum<UNIVERSE, double, &UNIVERSE::GetANNEmuGeV, EVENT>;
    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    using ANNMuonEnergyMaxGeV = PlotUtils::Maximum<UNIVERSE, double, &UNIVERSE::GetANNEmuGeV, EVENT>;
};


namespace truth
{
  template <class UNIVERSE>
  class IsInTarget: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      IsInTarget(const int targetCode, const bool useExtendedTarget): PlotUtils::SignalConstraint<UNIVERSE>("Is In Target " + std::to_string(targetCode)), fTargetCode(targetCode), fUseExtendedTarget(useExtendedTarget)
      {
      }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        return util::getTgtCode(&univ, true, fUseExtendedTarget) == fTargetCode;
      }
    int fTargetCode;
    bool fUseExtendedTarget;
  };
};

namespace util
{

  std::vector<double> PTBins = {0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1, 1.25, 1.5, 2.5, 4.5},
                      PzBins = {1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20, 40, 60},
                      EmuBins = {/* 0,1, */2,3,4,5,7,9,12,15,20},
                      //bjorkenXbins = {0.001, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.6, 1 , 2.2},
                      bjorkenXbins = {0.001, 0.02, 0.05, 0.1, 0.2, 0.4, 1 , 2.2},
                      bjorkenYbins = {0.001, 0.02, 0.05, 0.1, 0.2, 0.4, 1 , 2.2},
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
    cuts.emplace_back(new reco::ANNMuonEnergyMinGeV<CVUniverse, MichelEvent>(2, "EMu Min"));
    cuts.emplace_back(new reco::ANNMuonEnergyMaxGeV<CVUniverse, MichelEvent>(20, "EMu Max"));
    cuts.emplace_back(new reco::ANNConfidenceCut<CVUniverse, MichelEvent>(0.40)); //Reccomended at 0.4 for P6 ML vertexing and 0.2 for P4 vertexing
    cuts.emplace_back(new reco::RockMuonCut<CVUniverse, MichelEvent>()); //Reccomended for P6 ML vertexing
    cuts.emplace_back(new reco::VetoWall<CVUniverse, MichelEvent>()); //Reccomended for P6 ML vertexing
    return cuts;
  }

  PlotUtils::Cutter<CVUniverse, MichelEvent>::truth_t GetPhaseSpace()
  {
    PlotUtils::Cutter<CVUniverse, MichelEvent>::truth_t phasespace;
    phasespace.emplace_back(new truth::Apothem<CVUniverse>(apothem));
    phasespace.emplace_back(new truth::MuonAngle<CVUniverse>(17.));
    phasespace.emplace_back(new truth::MuonEnergyMinGeV<CVUniverse>(2, "EMu Min"));
    phasespace.emplace_back(new truth::MuonEnergyMaxGeV<CVUniverse>(20, "EMu Max"));
    return phasespace;
  }
};

/* namespace truth
{

  template <class UNIVERSE>
  class PZMuMin: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      PZMuMin(const double min): PlotUtils::SignalConstraint<UNIVERSE>(std::string("PzMu > ") + std::to_string(min)), fMin(min)
      {
      }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        return univ.GetPlepTrue() * cos(univ.GetThetalepTrue()) >= fMin;
      }

      const double fMin;
  };
}; */
#endif //UTIL_NUKEUTILS_H