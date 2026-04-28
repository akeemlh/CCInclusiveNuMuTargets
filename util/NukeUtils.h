#ifndef UTIL_NUKEUTILS_H
#define UTIL_NUKEUTILS_H

#include <map>
#include <algorithm>
#include <string>
#include <filesystem>

#include "event/CVUniverse.h"
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
#include "utilities/SMSCalibrator.h"

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
                                              {2, "RES"},
                                              {3, "DIS"},
                                              {8, "2p2h"},
                                              {30, "SoftDIS"} //This is not a GENIE label, I'm adding this, I hope there isn't already a GENIE enum that has value 30 because then we'll have a clash
                                             }; //kNoInteraction=0, kQEL=1, kRES=2, kDIS=3, kCOHPI=4, kAMNUGAMMA=5, kIMD=6, kNUEEL=7, k2P2H=8

    std::map<int, std::string> BKGLabels = {{0, "NC_Bkg"},
					       {1, "Wrong_Sign_Bkg"}};

    std::map<int, std::string> BKGLabelsWithPlasticSidebands = {{0, "NC_Bkg"},
					       {1, "Wrong_Sign_Bkg"}, {2, "Upstream_Plastic_Bkg"}, {3, "Downstream_Plastic_Bkg"}, {4, "Water_Tank_Bkg"}, {5, "True_In_Other_Target_Bkg"}, {6, "True_Vtx_Elsewhere_Bkg"}};

    int nuOrAntiNuMode(std::string playlist)
    {
        static std::vector<std::string> nuVector = {"minervame1A", "minervame1B", "minervame1C", "minervame1D", "minervame1E", "minervame1F", "minervame1G", "minervame1L", "minervame1M", "minervame1N", "minervame1O", "minervame1P"};
        static std::vector<std::string> anuVector = {"minervame5A", "minervame6A", "minervame6B", "minervame6C", "minervame6D", "minervame6E", "minervame6F", "minervame6G", "minervame6H", "minervame6I", "minervame6J"};
        if (std::find(std::begin(nuVector), std::end(nuVector), playlist) != std::end(nuVector)) return 1;
        else if (std::find(std::begin(anuVector), std::end(anuVector), playlist) != std::end(anuVector)) return 2;
        //Is this the most efficient way? Probably not
        return 0;
    }

    int filledOrEmptyMEPlaylist(std::string playlist)
    {
        static std::vector<std::string> filled = {"minervame1L", "minervame1M", "minervame1N", "minervame1O", "minervame1P", "1L", "1M", "1N", "1O", "1P",  "minervame5A", "minervame6A", "minervame6B", "minervame6C", "minervame6D", "minervame6E", "minervame6F", "minervame6G", "minervame6H", "5A", "6A", "6B", "6C", "6D", "6E", "6F", "6G", "6H"};
        static std::vector<std::string> empty = {"minervame1A", "minervame1B", "minervame1C", "minervame1D", "minervame1E", "minervame1F", "minervame1G", "1A", "1B", "1C", "1D", "1E", "1F", "1G", "minervame6J", "6J"/*,  "minervame5A", "minervame6A", "minervame6B", "minervame6C", "minervame6D", "minervame6E", "minervame6F", "minervame6G", "minervame6H", "minervame6I", "minervame6J" */};
        if (std::find(std::begin(filled), std::end(filled), playlist) != std::end(filled)) return 1;
        else if (std::find(std::begin(empty), std::end(empty), playlist) != std::end(empty)) return 2;
        //Is this the most efficient way? Probably not
        return 0;
    }



    //Taken from Oscar's code (/exp/minerva/app/users/omorenop/cmtuser/git-Mat/Personal/Test/InclusiveUtils.h)
    int planecode(int mdl, int plane)
    {
        if (mdl == -999 && plane == -999) return 36;

        if (plane != 1 && plane != 2) return -1;

        static constexpr std::array<std::array<int, 2>, 125> planeMap = {{
            /* -5  */ {{1, 2}},    {{3, 4}},    {{5, 6}},    {{7, 8}},    {{9, 9}},    // -5 to -1 (Target 1 at -1)
            /* 0  */ {{10, 11}},  {{12, 13}},  {{14, 15}},  {{16, 17}},  {{18, -1}},  //  0 to 4  (Target 2 at 4)
            /* 5  */ {{19, 20}},  {{21, 22}},  {{23, 24}},  {{25, 26}},  {{27, 27}},  //  5 to 9  (Target 3 at 9)
            /* 10  */ {{-1, -1}},  {{28, 29}},  {{30, 31}},  {{32, 33}},  {{34, 35}},  // 10 to 14 (10 unused)
            /* 15  */ {{37, 38}},  {{39, 40}},  {{41, 42}},  {{43, 44}},  {{45, 45}},  // 15 to 19 (Target 4 at 19)
            /* 20  */ {{46, 47}},  {{48, 49}},  {{50, 50}},  {{51, 52}},  {{53, 54}},  // 20 to 24 (Target 5 at 22)
            /* 25  */ {{55, 56}},  {{57, 58}},  {{59, 60}},  {{61, 62}},  {{63, 64}},  // 25 to 29
            /* 30  */ {{65, 66}},  {{67, 68}},  {{69, 70}},  {{71, 72}},  {{73, 74}},  // 30 to 34
            /* 35  */ {{75, 76}},  {{77, 78}},  {{79, 80}},  {{81, 82}},  {{83, 84}},  // 35 to 39
            /* 40  */ {{85, 86}},  {{87, 88}},  {{89, 90}},  {{91, 92}},  {{93, 94}},  // 40 to 44
            /* 45  */ {{95, 96}},  {{97, 98}},  {{99, 100}}, {{101, 102}},{{103, 104}}, // 45 to 49
            /* 50  */ {{105, 106}},{{107, 108}},{{109, 110}},{{111, 112}},{{113, 114}}, // 50 to 54
            /* 55  */ {{115, 116}},{{117, 118}},{{119, 120}},{{121, 122}},{{123, 124}}, // 55 to 59
            /* 60  */ {{125, 126}},{{127, 128}},{{129, 130}},{{131, 132}},{{133, 134}}, // 60 to 64
            /* 65  */ {{135, 136}},{{137, 138}},{{139, 140}},{{141, 142}},{{143, 144}}, // 65 to 69
            /* 70  */ {{145, 146}},{{147, 148}},{{149, 150}},{{151, 152}},{{153, 154}}, // 70 to 74
            /* 75  */ {{155, 156}},{{157, 158}},{{159, 160}},{{161, 162}},{{163, 164}}, // 75 to 79
            /* 80  */ {{165, 166}},{{167, 168}},{{169, 170}},{{171, 172}},{{173, 174}}, // 80 to 84
            /* 85  */ {{175, 176}},{{177, 178}},{{179, 180}},{{181, 182}},{{183, 184}}, // 85 to 89
            /* 90  */ {{185, 186}},{{187, 188}},{{189, 190}},{{191, 192}},{{193, 194}}, // 90 to 94
            /* 95  */ {{195, 195}},{{196, 196}},{{197, 197}},{{198, 198}},{{199, 199}}, // 95 to 99
            /* 100 */ {{-1, 200}}, {{-1, 201}}, {{-1, 202}}, {{-1, 203}}, {{-1, 204}},  // 100 to 104
            /* 105 */ {{-1, 205}}, {{-1, 206}}, {{-1, 207}}, {{-1, 208}}, {{-1, 209}},  // 105 to 109
            /* 110 */ {{-1, 210}}, {{-1, 211}}, {{-1, 212}}, {{-1, 213}}, {{-1, 214}},  // 110 to 114
            /* 115 */ {{-1, 2}},   {{-1, 0}},   {{-1, 1}},   {{-1, 0}},   {{-1, 2}}     // 115 to 119
        }};

        if (mdl >= -5 && mdl <= 119) {
            return planeMap[mdl + 5][plane - 1];
        }
        return -1; 
    }

    static PlotUtils::TargetUtils tgtUtil; 

    int getExtendedTarget(int mod, int plane, double vtx_x, double vtx_y) //Get which target this module and plane corresponds to if using extended target definition
    {
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
            if (scale!=1) hist1.Scale(scale);
        }
    }

    template<class HIST>
    void AddHist(HIST& hist1, HIST& hist2, double scale = 1)
    {
        if (hist1.GetEntries() > 0) //If not null
        {
            hist1.Add(*hist2, scale);
        }
        else //If null
        {
            hist1 = *(hist2.Clone());
            if (scale!=1) hist1.Scale(scale);
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
    int getPlasticPseudoTargetCode(int mod, int plane, bool extTarget = true)
    {
        if (extTarget) //If using extended target definition, exclude planes that would be included
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
        else return -1;
    }

    int getTgtCode(const CVUniverse *universe, bool truth, bool useExtendedTarget, bool useCalibrator = true)
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
            if (tgtUtil.InWaterTargetVolMC( vtx_x, vtx_y, vtx_z, 850. )) return 6000;
        }
        else // ANN
        {
            //First step, what target does the tuple say it's in
            int ANNTgtCode = universe->GetANNTargetCode(useCalibrator);
            if ((ANNTgtCode > 1000) && (ANNTgtCode < 7000)) //In the range that the nuclear target codes are in
            {
                return ANNTgtCode;
            }
            //Water is a special case, due to (I think?) a bug in the MAT tuples. To figure out if an event happened in water we cant use the regular tuple branch, instead we check if the ANN reconstruction put in the water segment (36)
            //if (universe->GetANNSegment() == 36 ) return 6000;

            mod = universe->GetANNVtxModule(useCalibrator);
            plane = universe->GetANNVtxPlane(useCalibrator);
            ROOT::Math::XYZVector Vtx = universe->GetANNVertex(useCalibrator);
            vtx_x = Vtx.X();
            vtx_y = Vtx.Y();
            vtx_z = Vtx.Z();
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
        if (targetCode < 1000) return false; //No such thing as a plastic sideband for the pseudo targets

        int mod, plane;
        double vtx_x = 0, vtx_y = 0;

        if (mode == 0) {
            mod = universe->GetTruthVtxModule();
            plane = universe->GetTruthVtxPlane();
        } else if (mode == 1) {
            mod = universe->GetANNVtxModule();
            plane = universe->GetANNVtxPlane();
        } else {
            mod = universe->GetMADVtxModule();
            plane = universe->GetMADVtxPlane();
        }

        // Dont cound plane immediately neighboring the nuclear targets
        if (removeNeighbors) {
            switch (mod) {
                case -2: if (plane == 2) return false; break;
                case  0: if (plane == 1) return false; break;
                case  3: if (plane == 2) return false; break;
                case  5: if (plane == 1) return false; break;
                case  8: if (plane == 2) return false; break;
                case 11: if (plane == 1) return false; break;
                case 14: if (plane == 2) return false; break;
                case 15: if (plane == 1) return false; break;
                case 18: if (plane == 2) return false; break;
                case 20: if (plane == 1) return false; break;
                case 21: if (plane == 2) return false; break;
            }
        }

        // Early exit if mod is out of all target ranges
        if (!((mod >= 0  && mod <= 3)  ||
            (mod >= 5  && mod <= 8)  ||
            (mod >= 11 && mod <= 14) ||
            (mod >= 15 && mod <= 18) ||
            (mod >= 20 && mod <= 26))) {
            return false;
        }

        if (mode == 0) {
            ROOT::Math::XYZTVector Vtx = universe->GetTrueVertex();
            vtx_x = Vtx.X();
            vtx_y = Vtx.Y();
        } else if (mode == 1) {
            ROOT::Math::XYZVector Vtx = universe->GetANNVertex();
            vtx_x = Vtx.X();
            vtx_y = Vtx.Y();
        } else {
            ROOT::Math::XYZTVector Vtx = universe->GetVertex();
            vtx_x = Vtx.X();
            vtx_y = Vtx.Y();
        }

        static const double zTgt1 = tgtUtil.GetTarget1CenterZMC();
        static const double zTgt2 = tgtUtil.GetTarget2CenterZMC();
        static const double zTgt3 = tgtUtil.GetTarget3CenterZMC();
        static const double zTgt4 = tgtUtil.GetTarget4CenterZMC();
        static const double zTgt5 = tgtUtil.GetTarget5CenterZMC();
        static const double zWater = (PlotUtils::TargetProp::WaterTarget::Face + 
                                    PlotUtils::TargetProp::WaterTarget::Back) / 2;

        if (mod >= 0 && mod <= 3) {
            if (USorDS) {
                if (targetCode == 2026) return tgtUtil.InIron2VolMC(vtx_x, vtx_y, zTgt2, 850., true);
                if (targetCode == 2082) return tgtUtil.InLead2VolMC(vtx_x, vtx_y, zTgt2, 850., true);
            } else {
                if (targetCode == 1026) return tgtUtil.InIron1VolMC(vtx_x, vtx_y, zTgt1, 850., true);
                if (targetCode == 1082) return tgtUtil.InLead1VolMC(vtx_x, vtx_y, zTgt1, 850., true);
            }
            return false;
        }
        if (mod >= 5 && mod <= 8) {
            if (USorDS) {
                if (targetCode == 3006) return tgtUtil.InCarbon3VolMC(vtx_x, vtx_y, zTgt3, 850., true);
                if (targetCode == 3026) return tgtUtil.InIron3VolMC(vtx_x, vtx_y, zTgt3, 850., true);
                if (targetCode == 3082) return tgtUtil.InLead3VolMC(vtx_x, vtx_y, zTgt3, 850., true);
            } else {
                if (targetCode == 2026) return tgtUtil.InIron2VolMC(vtx_x, vtx_y, zTgt2, 850., true);
                if (targetCode == 2082) return tgtUtil.InLead2VolMC(vtx_x, vtx_y, zTgt2, 850., true);
            }
            return false;
        }
        if (mod >= 11 && mod <= 14) {
            if (USorDS) {
                if (targetCode == 6000) return tgtUtil.InWaterTargetVolMC(vtx_x, vtx_y, zWater, 850.);
            } else {
                if (targetCode == 3006) return tgtUtil.InCarbon3VolMC(vtx_x, vtx_y, zTgt3, 850., true);
                if (targetCode == 3026) return tgtUtil.InIron3VolMC(vtx_x, vtx_y, zTgt3, 850., true);
                if (targetCode == 3082) return tgtUtil.InLead3VolMC(vtx_x, vtx_y, zTgt3, 850., true);
            }
            return false;
        }
        if (mod >= 15 && mod <= 18) {
            if (USorDS) {
                if (targetCode == 4082) return tgtUtil.InLead4VolMC(vtx_x, vtx_y, zTgt4, 850.);
            } else {
                if (targetCode == 6000) return tgtUtil.InWaterTargetVolMC(vtx_x, vtx_y, zWater, 850.);
            }
            return false;
        }
        if (mod >= 20 && mod <= 21) {
            if (USorDS) {
                if (targetCode == 5026) return tgtUtil.InIron5VolMC(vtx_x, vtx_y, zTgt5, 850., true);
                if (targetCode == 5082) return tgtUtil.InLead5VolMC(vtx_x, vtx_y, zTgt5, 850., true);
            } else {
                if (targetCode == 4082) return tgtUtil.InLead4VolMC(vtx_x, vtx_y, zTgt4, 850.);
            }
            return false;
        }
        if (mod >= 23 && mod <= 26) {
            if (!USorDS) {
                if (targetCode == 5026) return tgtUtil.InIron5VolMC(vtx_x, vtx_y, zTgt5, 850., true);
                if (targetCode == 5082) return tgtUtil.InLead5VolMC(vtx_x, vtx_y, zTgt5, 850., true);
            }
            return false;
        }

        return false;
    }

/*     {
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
    } */


    SMSCalibrator loadCalibrator(bool nuMode /* True = FHC, False = RHC */, bool filled /* True = Filled, False = Empty */)
    {
        std::string calibratorEnv = std::string(std::getenv("SMSCALIBRATORFILES"));
        std::string nuModeStr = nuMode ? "FHC" : "RHC";
        std::string filledStr = filled ? "Filled" :  "Empty";
        std::string calibPath = std::string(calibratorEnv)+"/"+nuModeStr+"/"+filledStr+"/sms_calibrator.json";
        SMSCalibrator calibrator(calibPath);
        if (calibrator.num_classes == -1) 
        {
            std::string errStr = "Failed to load calibrator from file " + calibPath;
            throw std::runtime_error(errStr);
        }
            else
        {
            std::cout<<"Successfully loaded calibrator "<< calibPath <<"\n";
            return calibrator;
        }
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
            //if(univ.GetInt((univ.GetAnaToolName() + "_minos_used_curvature").c_str()) != 1) return true;
            double relativeErr = 1/univ.GetMuonQPErr(); //Need to do 1/err because of how MINOS reports significance
            if (univ.GetAnalysisNuPDG()>0) return ( relativeErr <= -fMin );
            else if (univ.GetAnalysisNuPDG()<0) return ( relativeErr >= fMin );
            return false;
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
            //std::cout<< "Target code: " << fTargetCode << " Reco target code: " << tgt << " universe: " << univ.ShortName() << std::endl;
            return tgt == fTargetCode;
        }
        int fTargetCode;
        bool fUseExtendedTarget;
    };

    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class IsInPetal: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
        public:
        IsInPetal(int petal): PlotUtils::Cut<UNIVERSE, EVENT>("Is In Petal " + std::to_string(petal)), fPetal(petal)
        {
        }

        private:
        bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
        {   
            PlotUtils::TargetUtils tgtUtil;
            int petal = tgtUtil.GetDaisyPetal(univ.GetANNVertex().X(), univ.GetANNVertex().Y() );
            return petal==fPetal;
        }
        int fPetal;
    };

    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class InANNFiducial: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
        public:
        InANNFiducial(): PlotUtils::Cut<UNIVERSE, EVENT>("Is In ANN Fiducial Area")
        {
        }

        private:
        bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
        {   
            return univ.InANNFiducial();
        }
    };

    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class MuCoilCut: public PlotUtils::Cut<UNIVERSE, EVENT> //Unused in this analysis
    {
        public:
        MuCoilCut(): PlotUtils::Cut<UNIVERSE, EVENT>("MINOS Muon Coil Cut")
        {
        }

        private:
        bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
        {   
            //Code taken from Anezka
            const double coilXPos = 1219.0;
            const double coilYPos = 393.0;
            const double minos_x = univ.GetDouble((univ.GetAnaToolName() + "_minos_trk_end_x").c_str()) + coilXPos;
            const double minos_y = univ.GetDouble((univ.GetAnaToolName() + "_minos_trk_end_y").c_str()) + coilYPos;
            double minosR = sqrt(pow(minos_x,2) + pow(minos_y,2) );
            /*if( !((pow(minos_x,2) + pow(minos_y,2) )>= pow(coilR, 2)) )
            cout << minos_x << " " << minos_y << " " << coilR << endl;
            */
            return (minosR > 210 && minosR < 2500 );
        }
    };

    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class MuCorrectedPBug: public PlotUtils::Cut<UNIVERSE, EVENT> //Unused in this analysis
    {
        public:
        MuCorrectedPBug(): PlotUtils::Cut<UNIVERSE, EVENT>("Filtering muon_corrected_p bug events")
        {
        }

        private:
        bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
        {   
            return univ.GetInt("muon_corrected_p") != -999;
        }
    };

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class ANNApothem: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    ANNApothem(const double apothem): PlotUtils::Cut<UNIVERSE, EVENT>(std::string("Apothem ") + std::to_string(apothem)), fApothem(apothem), fSlope(-1./sqrt(3.))//A regular hexagon has angles of 2*M_PI/3, so I can find this is 1/tan(M_PI/3.)
      {
      }

    private:
      bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
      {
        const ROOT::Math::XYZVector vertex = univ.GetANNVertex();
        return (fabs(vertex.y()) < fSlope*fabs(vertex.x()) + 2.*fApothem/sqrt(3.))
               && (fabs(vertex.x()) < fApothem);
      }

      const double fApothem;
      const double fSlope; 
  };

    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    using ANNMuonEnergyMinGeV = PlotUtils::Minimum<UNIVERSE, double, &UNIVERSE::GetANNEmuGeV, EVENT>;
    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    using ANNMuonEnergyMaxGeV = PlotUtils::Maximum<UNIVERSE, double, &UNIVERSE::GetANNEmuGeV, EVENT>;
    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    using ANNMaxEnuGeV = PlotUtils::Maximum<UNIVERSE, double, &UNIVERSE::GetANNEnuGeV, EVENT>;
};


namespace truth
{
  static PlotUtils::TargetUtils tgtUtil; 
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

  template <class UNIVERSE>
  class IsInPetal: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      IsInPetal(const int petal): PlotUtils::SignalConstraint<UNIVERSE>("Is In Petal " + std::to_string(petal)), fPetal(petal)
      {
      }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        int petal = tgtUtil.GetDaisyPetal(univ.GetTrueVertex().X(), univ.GetTrueVertex().Y() );
        return petal==fPetal;
      }
    int fPetal;
  };

  template <class UNIVERSE>
  class InHexagonTrue: public PlotUtils::SignalConstraint<UNIVERSE> //Redundant with the truth::Apothem cut
  {
    public:
      InHexagonTrue(): PlotUtils::SignalConstraint<UNIVERSE>("Is In Fiducial Hexagon")
      {
      }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        return tgtUtil.IsInHexagon(univ.GetTrueVertex().X(), univ.GetTrueVertex().Y(), 850.); //Is there a mc "_in_fiducial_area" branch and if so should I use that?
      }
  };

};

namespace util
{

  std::vector<double> PTBins = {0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1, 1.25, 1.5, 2.5, 4.5},
                      PzBins = {1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20},//, 40, 60},
                      EmuBins = {/* 0,1, */2,3,4,5,7,9,12,15,20},
                      EnuBins = { 2.0, 3.0, 4.0, 5.0, 6.25, 7.5, 8.75, 10., 12.5, 15., 20.},
                      //bjorkenXbins = {0.001, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.6, 1 , 2.2},
                      //bjorkenXbins = {0.001, 0.02, 0.05, 0.1, 0.2, 0.4, 1 , 2.2},
                      bjorkenXbins = {0.001, 0.02, 0.05, 0.1, 0.2, 0.4, 2.2},
                      bjorkenXbins2 = {0.001, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8, 2.2},
                      bjorkenXbins3 = {0.001, 0.02, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 2.2},
                      bjorkenXbins4 = {0.001, 0.02, 0.05, 0.1, 0.2, 0.6, 2.2},
                      bjorkenXbinsBrian = {0.0, 0.1, 0.3, 0.7, 0.9, 1.1, 1.5},
                      bjorkenXbinsBrianMod = {0.0, 0.1, 0.3, 0.7, 0.9, 1.5},
                      bjorkenXbinsTest = {0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2},
                      //bjorkenXbins = {0.001, 0.02, 0.05, 0.1, 0.25, 1 , 2.2},
                      bjorkenYbins = {0.001, 0.02, 0.05, 0.1, 0.2, 0.4, 1 },
                      Erecoilbins = {0, 0.2, 0.75, 1.5, 2.5, 3.5, 5, 8, 11, 14.5, 20}; //GeV
                      //Erecoilbins = {0, 200, 750, 1500, 2500, 3500, 5000, 8000, 12000, 14500, 20000}; //MeV
  
  const double apothem = 850; //All in mm

  //Boiler plate precuts for CC inclusive NuMu ME Analysis
  //Note, the way I have these cuts organised is not the most computationally efficient
  //For maximum efficiency we would have the cuts that are quicker/cheaper to evaluate early on
  //So we don't run the more taxing evaluations only for those events to be cut later on.
  //Additionally, you'd also want to put the cuts that are most likely to fail earlier on so that 
  //the short-circuit evaluation logic we have triggers a fast fail in failing cases
  //Instead, I've used this order for readability/interpretability, grouping cuts by function
  PlotUtils::Cutter<CVUniverse, MichelEvent>::reco_t GetAnalysisCuts(int pdg)
  {
    PlotUtils::Cutter<CVUniverse, MichelEvent>::reco_t cuts;
    //Data Quality Cuts
    cuts.emplace_back(new reco::MuCorrectedPBug<CVUniverse, MichelEvent>());
    cuts.emplace_back(new reco::HasMINOSMatch<CVUniverse, MichelEvent>()); //Not in Anezkas!
    cuts.emplace_back(new reco::NoDeadtime<CVUniverse, MichelEvent>(1, "Deadtime"));
    cuts.emplace_back(new reco::RockMuonCut<CVUniverse, MichelEvent>()); //Reccomended for P6 ML vertexing
    cuts.emplace_back(new reco::VetoWall<CVUniverse, MichelEvent>()); //Reccomended for P6 ML vertexing
    cuts.emplace_back(new reco::MuonCurveSignificance<CVUniverse, MichelEvent>(5));
    //End - Data Quality Cuts
    //ANN Reco Cuts
    cuts.emplace_back(new reco::HasMLPredCut<CVUniverse, MichelEvent>());  //Not in Anezkas!
    cuts.emplace_back(new reco::ANNConfidenceCut<CVUniverse, MichelEvent>(0.40)); //Recomended at 0.4 for P6 ML vertexing and 0.2 for P4 vertexing //Anezkas is at 0.2
    //End - ANN Reco Cuts
    //Geometry/Kinematic Cuts
    cuts.emplace_back(new reco::ANNApothem<CVUniverse, MichelEvent>(apothem)); //May be redundant with the cut below //Not used in Anezkas reco
    cuts.emplace_back(new reco::InANNFiducial<CVUniverse, MichelEvent>());
    cuts.emplace_back(new reco::MinMuonAngle<CVUniverse, MichelEvent>(0.)); //May be unnecessary
    cuts.emplace_back(new reco::MaxMuonAngle<CVUniverse, MichelEvent>(17.));
    //if (pdg>0)  cuts.emplace_back(new reco::IsNeutrino<CVUniverse, MichelEvent>()); //Uses minos curvature
    //else if (pdg<0)  cuts.emplace_back(new reco::IsAntiNeutrino<CVUniverse, MichelEvent>()); //Uses minos curvature
    cuts.emplace_back(new reco::Helicity<CVUniverse, MichelEvent>(pdg));
    cuts.emplace_back(new reco::ANNMuonEnergyMinGeV<CVUniverse, MichelEvent>(2, "Emu Min (GeV)")); //!!! Anezka uses the old Emu!!!?!
    cuts.emplace_back(new reco::ANNMuonEnergyMaxGeV<CVUniverse, MichelEvent>(20, "Emu Max (GeV)")); //!!! Anezka uses the old Emu!!!?!
    cuts.emplace_back(new reco::ANNMaxEnuGeV<CVUniverse, MichelEvent>(120., "Enu Max (GeV)")); //!!! Anezka uses the old Enu!!!?!
    //End - Geometry Cuts
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

  bool IsTrueDIS(const CVUniverse& univ)
  {
    // if(q2 >= 1.0 && 1.3 <=W &&  W < 4.0){
    if(univ.GetTrueExperimentersQ2() >= 1 && univ.GetTrueExperimentersW() >= 2)
    {
        return true;
    }
    else
    {
        return false;
    }
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