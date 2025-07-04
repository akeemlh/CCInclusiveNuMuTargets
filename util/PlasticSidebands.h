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
    static std::map<int, int> TargetUSPlanes;
    static std::map<int, int> TargetDSPlanes;
    
    static std::map<int, int> TargetUSPlanesFlipped;
    static std::map<int, int> TargetDSPlanesFlipped;
    //Is this plane upstream of one of the nuclear targets?
    static std::map<int, int> plane1ModIDMap;

    //Plane2
    static std::map<int, int> plane2ModIDMap;

    //Verify US and DS planes for water target, i.e tgt 6.
    //To do: Change from ModPla code to segment system? See MasterAnaDev::GetTargetFromSegment()
    static std::map<int, int> USModPlaCodeToTgtId {{-18, 1}, {32, 2}, {82, 3}, {182, 4}, {212, 5}, {142, 6}};
    static std::map<int, int> DSModPlaCodeToTgtId {{1, 1}, {51, 2}, {111, 3}, {201, 4}, {231, 5}, {151, 6}};

    //Is this plane upstream of one of the nuclear targets?
    //modPlaCode is just modnum*10+plane ; modnum is in range 1-5 and plane is either 1, 2
    static int isUSPlane(int modPlaCode)
    {
        auto USTgtID = USModPlaCodeToTgtId.find(modPlaCode);
        if (USTgtID != util::USModPlaCodeToTgtId.end()) //Is event reconstructed immediately upstream of a nuclear target
        {
            return USTgtID->second;
        }
        else
        {
            return 0;
        }
    }

    static int isDSPlane(int modPlaCode)
    {
        auto DSTgtID = DSModPlaCodeToTgtId.find(modPlaCode);
        if (DSTgtID != util::DSModPlaCodeToTgtId.end()) //Is event reconstructed immediately upstream of a nuclear target
        {
            return DSTgtID->second;
        }
        else
        {
            return 0;
        }
    }

    static int isUSSideband(int modPlaCode)
    {
        auto USTgtID = USModPlaCodeToTgtId.find(modPlaCode);
        if (USTgtID != util::USModPlaCodeToTgtId.end()) //Is event reconstructed immediately upstream of a nuclear target
        {
            return USTgtID->second;
        }
        else
        {
            return 0;
        }
    }

    static int isDSSideband(int modPlaCode)
    {
        auto DSTgtID = DSModPlaCodeToTgtId.find(modPlaCode);
        if (DSTgtID != util::DSModPlaCodeToTgtId.end()) //Is event reconstructed immediately upstream of a nuclear target
        {
            return DSTgtID->second;
        }
        else
        {
            return 0;
        }
    }
}
#endif //UTIL_PLASTICSIDEBANDS_H