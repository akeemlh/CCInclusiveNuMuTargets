#ifndef UTIL_NUKEUTILS_H
#define UTIL_NUKEUTILS_H

#include <map>
#include <algorithm>

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
    std::map<int, std::string> TgtCodeLabelsNuke = {{1026, "1026"}, {1082, "1082"}, {2026, "2026"}, {2082, "2082"}, {3006, "3006"}, {3026, "3026"}, {3082, "3082"}, {4082, "4082"}, {5026, "5026"}, {5082, "5082"}, {-999, "Water"}, {7, "Target7"}, {8, "Target8"}, {9, "Target9"}, {10, "Target10"}, {11, "Target11"}, {12, "Target12"}};

    std::map<int, std::string> TgtCodeLabelsNuke2 = {{1026, "1026"}, {1082, "1082"}, {2026, "2026"}, {2082, "2082"}, {3006, "3006"}, {3026, "3026"}, {3082, "3082"}, {4082, "4082"}, {5026, "5026"}, {5082, "5082"}, {-999, "Water"}, {7, "7"}, {8, "8"}, {9, "9"}, {10, "10"}, {11, "11"}, {12, "12"}};

    std::map<int, std::string> GENIELabels = {{1, "QE"},
                                              {8, "2p2h"},
                                              {2, "RES"},
                                              {3, "DIS"}};

    std::map<int, std::string> BKGLabels = {{0, "NC_Bkg"},
					       {1, "Wrong_Sign_Bkg"}, {2, "Wrong_Material_Bkg"}};
}
#endif //UTIL_NUKEUTILS_H