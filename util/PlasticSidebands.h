#ifndef UTIL_PLASTICSIDEBANDS_H
#define UTIL_PLASTICSIDEBANDS_H

#include <map>
#include <algorithm>

namespace util
{
    //From https://cdcvs.fnal.gov/redmine/projects/minerva/wiki/Z_Positions_of_Planes_in_the_Full_MINERvA_Detector
    //And https://cdcvs.fnal.gov/redmine/projects/minerva-sw/wiki/MINERvA_Detector
    //And https://cdcvs.fnal.gov/redmine/projects/minerva-ops/wiki/MINERvA_Detector
    //And https://cdcvs.fnal.gov/redmine/projects/minerva/wiki/Z_Positions_of_Planes_in_the_Downstream_Detector -- Not used
    std::map<int, int> TargetUSPlanes = {{1, 1339555840}, {2, 1211629568}, {3, 1216872448}, {4, 1227358208}, {5, 1230503936}};
    std::map<int, int> TargetDSPlanes = {{1, 1208221696}, {2, 1213464576}, {3, 1219756032}, {4, 1229193216}, {5, 1500774400}};
    
    //Is this plane upstream of one of the nuclear targets?
    bool isUSPlane(int planeID)
    {
        std::map<int, int>::iterator result = std::find_if(
            TargetUSPlanes.begin(),
            TargetUSPlanes.end(),
        [planeID](const auto& obj) {return obj.second == planeID; });

        if(result == TargetUSPlanes.end())
        {
            return false;
        }
        else
        {
            return true;
        }
    }

    //Is this plane downstream of one of the nuclear targets?
    bool isDSPlane(int planeID)
    {
        std::map<int, int>::iterator result = std::find_if(
            TargetDSPlanes.begin(),
            TargetDSPlanes.end(),
        [planeID](const auto& obj) {return obj.second == planeID; });

        if(result == TargetDSPlanes.end())
        {
            return false;
        }
        else
        {
            return true;
        }
    }
}
#endif