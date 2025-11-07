/* 
To whoever stumbles across this in the future:

This was intially intended to be a one-file-contains-all special event loop just to
produce validation plots or, more specifically, plots that aren't directly utilised
in extracting a cross section but are important to check or investigate something.
This was initially derived from runEventLoop.cpp but has since significantly changed

It has the same cuts and signal definition as the cross section event loops but this 
allows me to keep those event loops, which I run more often, less bloated so they run
faster, are easier to maintain/modify and have a lower memory footprint when running,
so I can make lower grid memory requests and get higher priority.

Currently main function of this macro is to produce z vertex plots initially intended
to determine how our vertexing fares with different cuts or models in different 
detector regions and between different vertexing approaches (ML vs Track Based). Now
also includes an ERecoil segmentation into 4 quartiles and 

Consequently, because this event loop is where I just throw *everything* that isn't
fed into the cross section pipeline, it's possible that it may in the future become
bloated and require fragmentation. */

/*-----------------------------------------------
Akeem Hart 
This script is used to generate plots of the true vertex posistions of interactions by whether they were reconstructed as being inside or outside the water target
by different reconstruction approaches

Build with:
g++ -g runEventLoopTracker.cpp util/**.cpp -o runEventLoopTracker -I $CCQENU_ANA/ -I$CCQENU_INCLUDE_PATH -ggdb -Wall -fPIC -I$PLOTUTILSROOT/../include/ -I/exp/minerva/app/users/alhart/MAT_AL9/opt/include -I. `root-config --cflags --glibs --libs` -L$PLOTUTILSROOT -lMAT -lMAT-MINERvA -lUnfoldUtils -L.

Or run with root -b -l TargetRegionValidation.cpp++g

root -b -l "TargetRegionValidation.cpp++g('/exp/minerva/app/users/alhart/MAT_AL9/MINERvA-101-Cross-Section/PlaylistFiles/me-playlists/DataP6/Test1A-Data.txt', '/exp/minerva/app/users/alhart/MAT_AL9/MINERvA-101-Cross-Section/PlaylistFiles/me-playlists/MCP6/Test1A-MC.txt')"
-----------------------------------------------*/

#define USAGE \
"\n*** USAGE ***\n"\
"runEventLoop <dataPlaylist.txt> <mcPlaylist.txt>\n\n"\
"*** Explanation ***\n"\
"Reduce MasterAnaDev AnaTuples to event selection histograms to extract a\n"\
"single-differential inclusive cross section for the 2021 MINERvA 101 tutorial.\n\n"\
"*** The Input Files ***\n"\
"Playlist files are plaintext files with 1 file name per line.  Filenames may be\n"\
"xrootd URLs or refer to the local filesystem.  The first playlist file's\n"\
"entries will be treated like data, and the second playlist's entries must\n"\
"have the \"Truth\" tree to use for calculating the efficiency denominator.\n\n"\
"*** Output ***\n"\
"Produces a two files, !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!, with\n"\
"all histograms needed for the ExtractCrossSection program also built by this\n"\
"package.  You'll need a .rootlogon.C that loads ROOT object definitions from\n"\
"PlotUtils to access systematics information from these files.\n\n"\
"*** Environment Variables ***\n"\
"Setting up this package appends to PATH and LD_LIBRARY_PATH.  PLOTUTILSROOT,\n"\
"MPARAMFILESROOT, and MPARAMFILES must be set according to the setup scripts in\n"\
"those packages for systematics and flux reweighters to function.\n"\
"If MNV101_SKIP_SYST is defined at all, output histograms will have no error bands.\n"\
"This is useful for debugging the CV and running warping studies.\n\n"\
"*** Return Codes ***\n"\
"0 indicates success.  All histograms are valid only in this case.  Any other\n"\
"return code indicates that histograms should not be used.  Error messages\n"\
"about what went wrong will be printed to stderr.  So, they'll end up in your\n"\
"terminal, but you can separate them from everything else with something like:\n"\
"\"runEventLoop data.txt mc.txt 2> errors.txt\"\n"

enum ErrorCodes
{
  success = 0,
  badCmdLine = 1,
  badInputFile = 2,
  badFileRead = 3,
  badOutputFile = 4
};

//PlotUtils includes
//No junk from PlotUtils please!  I already
//know that MnvH1D does horrible horrible things.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"

//Includes from this package
#include "event/CVUniverse.h"
#include "event/MichelEvent.h"
#include "systematics/Systematics.h"
#include "util/GetFluxIntegral.h"
#include "util/GetPlaylist.h"
#include "cuts/SignalDefinition.h"
//#include "Binning.h" //TODO: Fix me

//PlotUtils includes
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/MacroUtil.h"
#include "PlotUtils/MnvPlotter.h"
#include "PlotUtils/CCInclusiveCuts.h"
#include "PlotUtils/CCInclusiveSignal.h"
//#include "PlotUtils/CrashOnROOTMessage.h" //Sets up ROOT's debug callbacks by itself
#include "PlotUtils/Cutter.h"
#include "PlotUtils/Model.h"
#include "PlotUtils/FluxAndCVReweighter.h"
#include "PlotUtils/GENIEReweighter.h"
#include "PlotUtils/LowRecoil2p2hReweighter.h"
#include "PlotUtils/RPAReweighter.h"
#include "PlotUtils/MINOSEfficiencyReweighter.h"
#include "PlotUtils/LowQ2PiReweighter.h"
#include "PlotUtils/AMUDISReweighter.h"
#include "PlotUtils/SuSAFromValencia2p2hReweighter.h"
#include "PlotUtils/FSIReweighter.h"
#include "PlotUtils/TargetUtils.h"
#include "util/NukeUtils.h"

#include "util/COHPionReweighter.h"
#include "util/DiffractiveReweighter.h"
#include "PlotUtils/BodekRitchieReweighter.h"

#pragma GCC diagnostic pop

//ROOT includes
#include "TParameter.h"

#include "Math/Vector3D.h"
#include "TH3D.h"
#include "TTree.h"
#include "TLatex.h"
#include "TFrame.h"
#include "TMath.h"
#include <chrono>

//c++ includes
#include <iostream>
#include <cstdlib> //getenv()
#include <fstream>

//These 2 variables are used when doing broken down runs to speed up running on the grid
int nSubruns = 0;
int nProcess = 0;
bool verbose = false;

//std::vector<double> vertexBins = {4293.04, 4337.25, 4381.47, 4425.68, 4514.11, 4558.33, 4602.54, 4646.76, 4735.19, 4779.4, 4823.62, 4867.83, 5000.48, 5044.69, 5088.91, 5133.12, 5456.74, 5500.95, 5545.17, 5589.38, 5677.81, 5722.03, 5810.45, 5855.68, 5900.91, 5946.14, 5991.37, 6036.6, 6081.83, 6127.06, 6172.29, 6217.52, 6262.74, 6307.97, 6353.2, 6398.43, 6443.66, 6488.89, 6534.12, 6579.35, 6624.58, 6669.81, 6715.03, 6760.26, 6805.49, 6850.72, 6895.95, 6941.18, 6986.41, 7031.64, 7076.87, 7122.1, 7167.32, 7212.55, 7257.78, 7303.01, 7348.24, 7393.47, 7438.7, 7483.93, 7529.16, 7574.39, 7619.61, 7664.84, 7710.07, 7755.3, 7800.53, 7845.76, 7890.99, 7936.22, 7981.45, 8026.68, 8071.9, 8117.13, 8162.36, 8207.59, 8252.82, 8298.05, 8343.28, 8388.51, 8433.74, 8478.97, 8524.19, 8569.42, 8614.65};

std::vector<double> vertexBins = {4293.04, 4313.68, 4337.25, 4357.9, 4381.47, 4402.11, 4425.68, 4446.33, 4514.11, 4534.76, 4558.33, 4578.97, 4602.54, 4623.19, 4646.76, 4667.4, 4735.19, 4755.83, 4779.4, 4800.05, 4823.62, 4844.26, 4867.83, 4888.48, 5000.48, 5021.12, 5044.69, 5065.34, 5088.91, 5109.55, 5133.12, 5153.77, 5456.74, 5477.38, 5500.95, 5521.6, 5545.17, 5565.81, 5589.38, 5610.02, 5677.81, 5698.45, 5722.03, 5742.67, 5810.45, 5831.1, 5855.68, 5876.33, 5900.91, 5921.56, 5946.14, 5966.79, 5991.37, 6012.01, 6036.6, 6057.24, 6081.83, 6102.47, 6127.06, 6147.7, 6172.29, 6192.93, 6217.52, 6238.16, 6262.74, 6283.39, 6307.97, 6328.62, 6353.2, 6373.85, 6398.43, 6419.08, 6443.66, 6464.3, 6488.89, 6509.53, 6534.12, 6554.76, 6579.35, 6599.99, 6624.58, 6645.22, 6669.81, 6690.45, 6715.03, 6735.68, 6760.26, 6780.91, 6805.49, 6826.14, 6850.72, 6871.37, 6895.95, 6916.59, 6941.18, 6961.82, 6986.41, 7007.05, 7031.64, 7052.28, 7076.87, 7097.51, 7122.1, 7142.74, 7167.32, 7187.97, 7212.55, 7233.2, 7257.78, 7278.43, 7303.01, 7323.66, 7348.24, 7368.88, 7393.47, 7414.11, 7438.7, 7459.34, 7483.93, 7504.57, 7529.16, 7549.8, 7574.39, 7595.03, 7619.61, 7640.26, 7664.84, 7685.49, 7710.07, 7730.72, 7755.3, 7775.95, 7800.53, 7821.17, 7845.76, 7866.4, 7890.99, 7911.63, 7936.22, 7956.86, 7981.45, 8002.09, 8026.68, 8047.32, 8071.9, 8092.55, 8117.13, 8137.78, 8162.36, 8183.01, 8207.59, 8228.24, 8252.82, 8273.46, 8298.05, 8318.69, 8343.28, 8363.92, 8388.51, 8409.15, 8433.74, 8454.38, 8478.97, 8499.61, 8524.19, 8544.84, 8569.42, 8590.07, 8614.65, 8635.3, 8659.46, 8680.1, 8704.26, 8724.9, 8749.06, 8769.71, 8793.86, 8814.51, 8838.67, 8859.31, 8883.47, 8904.11, 8928.27, 8948.92, 8973.08, 8993.72, 9017.88, 9038.52, 9088.08, 9135.41, 9182.75, 9230.08, 9277.41, 9324.74, 9372.08, 9419.41, 9466.74, 9514.07, 9561.41, 9608.74, 9656.07, 9703.4, 9750.74, 9798.07, 9845.4, 9892.73, 9940.07, 9987.4};

double rebinNum = 10;
double nbins = 3400;
double xlow = 4200;
double xhigh = 5900;


//std::unique_ptr<PlotUtils::Cut<UNIVERSE, EVENT> >
auto ANNConfCut = new reco::ANNConfidenceCut<CVUniverse, MichelEvent>(0.20);


//==============================================================================
// Global - Declaring histograms
//==============================================================================

TH1D *ANNVerticesMC_ByModule = new TH1D ("ANNVerticesMC_ByModule", "ANNVerticesMC_ByModule", vertexBins.size()-1, &vertexBins[0]);
TH1D *TBVerticesMC_ByModule = new TH1D ("TBVerticesMC_ByModule", "TBVerticesMC_ByModule", vertexBins.size()-1, &vertexBins[0]);
TH1D *ANNVerticesData_ByModule = new TH1D ("ANNVerticesData_ByModule", "ANNVerticesData_ByModule", vertexBins.size()-1, &vertexBins[0]);
TH1D *TBVerticesData_ByModule = new TH1D ("TBVerticesData_ByModule", "TBVerticesData_ByModule", vertexBins.size()-1, &vertexBins[0]);
TH1D *TruthVerticesMC_ByModule = new TH1D ("TruthVerticesMC_ByModule", "TruthVerticesMC_ByModule", vertexBins.size()-1, &vertexBins[0]);
TH1D *ANNVerticesMC_ByZPos = new TH1D ("ANNVerticesMC_ByZPos", "ANNVerticesMC_ByZPos", 9000, 4200, 8700);
TH1D *TBVerticesMC_ByZPos = new TH1D ("TBVerticesMC_ByZPos", "TBVerticesMC_ByModuleTruthVerticesMC_ByZPosHighRes", 9000, 4200, 8700);
TH1D *ANNVerticesData_ByZPos = new TH1D ("ANNVerticesData_ByZPos", "ANNVerticesData_ByZPos", 9000, 4200, 8700);
TH1D *TBVerticesData_ByZPos = new TH1D ("TBVerticesData_ByZPos", "TBVerticesData_ByZPos", 9000, 4200, 8700);
TH1D *TruthVerticesMC_ByZPos = new TH1D ("TruthVerticesMC_ByZPos", "TruthVerticesMC_ByZPos", 9000, 4200, 8700);
TH1D *ANNVerticesData_ByMod = new TH1D ("ANNVerticesData_ByMod", "ANNVerticesData_ByMod", 125, -5, 120);
TH1D *TBVerticesData_ByMod = new TH1D ("TBVerticesData_ByMod", "TBVerticesData_ByMod", 125, -5, 120);
TH1D *TruthVerticesMC_ByMod = new TH1D ("TruthVerticesMC_ByMod", "TruthVerticesMC_ByMod", 125, -5, 120);

TH2D *ANNZWeightedVerticesData_VsCutoff = new TH2D ("ANNZWeightedVerticesData_VsCutoff", "ANNZWeightedVerticesData_VsCutoff", 9000, 4200, 8700, 100, 0, 100);
TH2D *ANNZWeightedVerticesMC_VsCutoff = new TH2D ("ANNZWeightedVerticesMC_VsCutoff", "ANNZWeightedVerticesMC_VsCutoff", 9000, 4200, 8700, 100, 0 , 100);

TH2D *ANNZWeightedVerticesData_VsEhad = new TH2D ("ANNZWeightedVerticesData_VsEhad", "ANNZWeightedVerticesData_VsEhad", 9000, 4200, 8700, 100, 0, 20);
TH2D *ANNZWeightedVerticesMC_VsEhad = new TH2D ("ANNZWeightedVerticesMC_VsEhad", "ANNZWeightedVerticesMC_VsEhad", 9000, 4200, 8700, 100, 0 , 20);

TH1D *ANNVerticesMC_BySegment = new TH1D ("ANNVerticesMC_BySegment", "ANNVerticesMC_BySegment", 220, 0, 220);
TH1D *TruthVerticesMC_BySegment = new TH1D ("TruthVerticesMC_BySegment", "TruthVerticesMC_BySegment", 220, 0, 220);

TH1D *ANNVerticesData_BySegment = new TH1D ("ANNVerticesData_BySegment", "ANNVerticesData_BySegment", 220, 0, 220);

TH1D *ANNZVertexResidual = new TH1D ("ANNZVertexResidual", "ANNZVertexResidual", 800, -400, 400);
TH1D *ANNZVertexResidualTracker = new TH1D ("ANNZVertexResidualTracker", "ANNZVertexResidualTracker", 800, -400, 400);
TH1D *ANNZVertexResidualTgt2Iron = new TH1D ("ANNZVertexResidualTgt2Iron", "ANNZVertexResidualTgt2Iron", 800, -400, 400);
TH1D *ANNZVertexResidualTgt2Lead = new TH1D ("ANNZVertexResidualTgt2Lead", "ANNZVertexResidualTgt2Lead", 800, -400, 400);
TH1D *ANNZVertexResidualTgt3Iron = new TH1D ("ANNZVertexResidualTgt3Iron", "ANNZVertexResidualTgt3Iron", 800, -400, 400);
TH1D *ANNZVertexResidualTgt3Lead = new TH1D ("ANNZVertexResidualTgt3Lead", "ANNZVertexResidualTgt3Lead", 800, -400, 400);
TH1D *ANNZVertexResidualTgt3Carbon = new TH1D ("ANNZVertexResidualTgt3Carbon", "ANNZVertexResidualTgt3Carbon", 800, -400, 400);
TH1D *ANNZVertexResidualTgt4Lead = new TH1D ("ANNZVertexResidualTgt4Lead", "ANNZVertexResidualTgt4Lead", 800, -400, 400);
TH1D *ANNZVertexResidualTgt5Iron = new TH1D ("ANNZVertexResidualTgt5Iron", "ANNZVertexResidualTgt5Iron", 800, -400, 400);
TH1D *ANNZVertexResidualTgt5Lead = new TH1D ("ANNZVertexResidualTgt5Lead", "ANNZVertexResidualTgt5Lead", 800, -400, 400);
TH1D *ANNZVertexResidualWater = new TH1D ("ANNZVertexResidualWater", "ANNZVertexResidualWater", 800, -400, 400);

TH2D *ANNZVertexResidualVsConf = new TH2D ("ANNZVertexResidualVsConf", "ANNZVertexResidualVsConf", 800, -400, 400, 10, 0, 1);
TH2D *ANNZVertexResidualVsConfTracker = new TH2D ("ANNZVertexResidualVsConfTracker", "ANNZVertexResidualVsConfTracker", 800, -400, 400, 10, 0, 1);
TH2D *ANNZVertexResidualVsConfTgt2Iron = new TH2D ("ANNZVertexResidualVsConfTgt2Iron", "ANNZVertexResidualVsConfTgt2Iron", 800, -400, 400, 10, 0, 1);
TH2D *ANNZVertexResidualVsConfTgt2Lead = new TH2D ("ANNZVertexResidualVsConfTgt2Lead", "ANNZVertexResidualVsConfTgt2Lead", 800, -400, 400, 10, 0, 1);
TH2D *ANNZVertexResidualVsConfTgt3Iron = new TH2D ("ANNZVertexResidualVsConfTgt3Iron", "ANNZVertexResidualVsConfTgt3Iron", 800, -400, 400, 10, 0, 1);
TH2D *ANNZVertexResidualVsConfTgt3Lead = new TH2D ("ANNZVertexResidualVsConfTgt3Lead", "ANNZVertexResidualVsConfTgt3Lead", 800, -400, 400, 10, 0, 1);
TH2D *ANNZVertexResidualVsConfTgt3Carbon = new TH2D ("ANNZVertexResidualVsConfTgt3Carbon", "ANNZVertexResidualVsConfTgt3Carbon", 800, -400, 400, 10, 0, 1);
TH2D *ANNZVertexResidualVsConfTgt4Lead = new TH2D ("ANNZVertexResidualVsConfTgt4Lead", "ANNZVertexResidualVsConfTgt4Lead", 800, -400, 400, 10, 0, 1);
TH2D *ANNZVertexResidualVsConfTgt5Iron = new TH2D ("ANNZVertexResidualVsConfTgt5Iron", "ANNZVertexResidualVsConfTgt5Iron", 800, -400, 400, 10, 0, 1);
TH2D *ANNZVertexResidualVsConfTgt5Lead = new TH2D ("ANNZVertexResidualVsConfTgt5Lead", "ANNZVertexResidualVsConfTgt5Lead", 800, -400, 400, 10, 0, 1);
TH2D *ANNZVertexResidualVsConfWater = new TH2D ("ANNZVertexResidualVsConfWater", "ANNZVertexResidualVsConfWater", 800, -400, 400, 10, 0, 1);

TH2D *ANNWeightedZVertexResidualVsConf = new TH2D ("ANNWeightedZVertexResidualVsConf", "ANNWeightedZVertexResidualVsConf", 800, -400, 400, 10, 0, 1);
TH2D *ANNWeightedZVertexResidualVsConfTracker = new TH2D ("ANNWeightedZVertexResidualVsConfTracker", "ANNWeightedZVertexResidualVsConfTracker", 800, -400, 400, 10, 0, 1);
TH2D *ANNWeightedZVertexResidualVsConfTgt2Iron = new TH2D ("ANNWeightedZVertexResidualVsConfTgt2Iron", "ANNWeightedZVertexResidualVsConfTgt2Iron", 800, -400, 400, 10, 0, 1);
TH2D *ANNWeightedZVertexResidualVsConfTgt2Lead = new TH2D ("ANNWeightedZVertexResidualVsConfTgt2Lead", "ANNWeightedZVertexResidualVsConfTgt2Lead", 800, -400, 400, 10, 0, 1);
TH2D *ANNWeightedZVertexResidualVsConfTgt3Iron = new TH2D ("ANNWeightedZVertexResidualVsConfTgt3Iron", "ANNWeightedZVertexResidualVsConfTgt3Iron", 800, -400, 400, 10, 0, 1);
TH2D *ANNWeightedZVertexResidualVsConfTgt3Lead = new TH2D ("ANNWeightedZVertexResidualVsConfTgt3Lead", "ANNWeightedZVertexResidualVsConfTgt3Lead", 800, -400, 400, 10, 0, 1);
TH2D *ANNWeightedZVertexResidualVsConfTgt3Carbon = new TH2D ("ANNWeightedZVertexResidualVsConfTgt3Carbon", "ANNWeightedZVertexResidualVsConfTgt3Carbon", 800, -400, 400, 10, 0, 1);
TH2D *ANNWeightedZVertexResidualVsConfTgt4Lead = new TH2D ("ANNWeightedZVertexResidualVsConfTgt4Lead", "ANNWeightedZVertexResidualVsConfTgt4Lead", 800, -400, 400, 10, 0, 1);
TH2D *ANNWeightedZVertexResidualVsConfTgt5Iron = new TH2D ("ANNWeightedZVertexResidualVsConfTgt5Iron", "ANNWeightedZVertexResidualVsConfTgt5Iron", 800, -400, 400, 10, 0, 1);
TH2D *ANNWeightedZVertexResidualVsConfTgt5Lead = new TH2D ("ANNWeightedZVertexResidualVsConfTgt5Lead", "ANNWeightedZVertexResidualVsConfTgt5Lead", 800, -400, 400, 10, 0, 1);
TH2D *ANNWeightedZVertexResidualVsConfWater = new TH2D ("ANNWeightedZVertexResidualVsConfWater", "ANNWeightedZVertexResidualVsConfWater", 800, -400, 400, 10, 0, 1);

TH2D *ANNWeightedZVertexResidual = new TH2D ("ANNWeightedZVertexResidual", "ANNWeightedZVertexResidual", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedZVertexResidualTracker = new TH2D ("ANNWeightedZVertexResidualTracker", "ANNWeightedZVertexResidualTracker", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedZVertexResidualTgt2Iron = new TH2D ("ANNWeightedZVertexResidualTgt2Iron", "ANNWeightedZVertexResidualTgt2Iron", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedZVertexResidualTgt2Lead = new TH2D ("ANNWeightedZVertexResidualTgt2Lead", "ANNWeightedZVertexResidualTgt2Lead", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedZVertexResidualTgt3Iron = new TH2D ("ANNWeightedZVertexResidualTgt3Iron", "ANNWeightedZVertexResidualTgt3Iron", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedZVertexResidualTgt3Lead = new TH2D ("ANNWeightedZVertexResidualTgt3Lead", "ANNWeightedZVertexResidualTgt3Lead", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedZVertexResidualTgt3Carbon = new TH2D ("ANNWeightedZVertexResidualTgt3Carbon", "ANNWeightedZVertexResidualTgt3Carbon", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedZVertexResidualTgt4Lead = new TH2D ("ANNWeightedZVertexResidualTgt4Lead", "ANNWeightedZVertexResidualTgt4Lead", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedZVertexResidualTgt5Iron = new TH2D ("ANNWeightedZVertexResidualTgt5Iron", "ANNWeightedZVertexResidualTgt5Iron", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedZVertexResidualTgt5Lead = new TH2D ("ANNWeightedZVertexResidualTgt5Lead", "ANNWeightedZVertexResidualTgt5Lead", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedZVertexResidualWater = new TH2D ("ANNWeightedZVertexResidualWater", "ANNWeightedZVertexResidualWater", 800, -400, 400, 101, 0, 100);

TH2D *ANNWeightedVsUnweightedZVertexDifferenceMC = new TH2D ("ANNWeightedVsUnweightedZVertexDifferenceMC", "ANNWeightedVsUnweightedZVertexDifferenceMC", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedVsUnweightedZVertexDifferenceTrackerMC = new TH2D ("ANNWeightedVsUnweightedZVertexDifferenceTrackerMC", "ANNWeightedVsUnweightedZVertexDifferenceTrackerMC", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedVsUnweightedZVertexDifferenceTgt2IronMC = new TH2D ("ANNWeightedVsUnweightedZVertexDifferenceTgt2IronMC", "ANNWeightedVsUnweightedZVertexDifferenceTgt2IronMC", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedVsUnweightedZVertexDifferenceTgt2LeadMC = new TH2D ("ANNWeightedVsUnweightedZVertexDifferenceTgt2LeadMC", "ANNWeightedVsUnweightedZVertexDifferenceTgt2LeadMC", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedVsUnweightedZVertexDifferenceTgt3IronMC = new TH2D ("ANNWeightedVsUnweightedZVertexDifferenceTgt3IronMC", "ANNWeightedVsUnweightedZVertexDifferenceTgt3IronMC", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedVsUnweightedZVertexDifferenceTgt3LeadMC = new TH2D ("ANNWeightedVsUnweightedZVertexDifferenceTgt3LeadMC", "ANNWeightedVsUnweightedZVertexDifferenceTgt3LeadMC", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedVsUnweightedZVertexDifferenceTgt3CarbonMC = new TH2D ("ANNWeightedVsUnweightedZVertexDifferenceTgt3CarbonMC", "ANNWeightedVsUnweightedZVertexDifferenceTgt3CarbonMC", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedVsUnweightedZVertexDifferenceTgt4LeadMC = new TH2D ("ANNWeightedVsUnweightedZVertexDifferenceTgt4LeadMC", "ANNWeightedVsUnweightedZVertexDifferenceTgt4LeadMC", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedVsUnweightedZVertexDifferenceTgt5IronMC = new TH2D ("ANNWeightedVsUnweightedZVertexDifferenceTgt5IronMC", "ANNWeightedVsUnweightedZVertexDifferenceTgt5IronMC", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedVsUnweightedZVertexDifferenceTgt5LeadMC = new TH2D ("ANNWeightedVsUnweightedZVertexDifferenceTgt5LeadMC", "ANNWeightedVsUnweightedZVertexDifferenceTgt5LeadMC", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedVsUnweightedZVertexDifferenceWaterMC = new TH2D ("ANNWeightedVsUnweightedZVertexDifferenceWaterMC", "ANNWeightedVsUnweightedZVertexDifferenceWaterMC", 800, -400, 400, 101, 0, 100);

TH2D *ANNWeightedVsUnweightedZVertexDifferenceData = new TH2D ("ANNWeightedVsUnweightedZVertexDifferenceData", "ANNWeightedVsUnweightedZVertexDifferenceData", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedVsUnweightedZVertexDifferenceTrackerData = new TH2D ("ANNWeightedVsUnweightedZVertexDifferenceTrackerData", "ANNWeightedVsUnweightedZVertexDifferenceTrackerData", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedVsUnweightedZVertexDifferenceTgt2IronData = new TH2D ("ANNWeightedVsUnweightedZVertexDifferenceTgt2IronData", "ANNWeightedVsUnweightedZVertexDifferenceTgt2IronData", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedVsUnweightedZVertexDifferenceTgt2LeadData = new TH2D ("ANNWeightedVsUnweightedZVertexDifferenceTgt2LeadData", "ANNWeightedVsUnweightedZVertexDifferenceTgt2LeadData", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedVsUnweightedZVertexDifferenceTgt3IronData = new TH2D ("ANNWeightedVsUnweightedZVertexDifferenceTgt3IronData", "ANNWeightedVsUnweightedZVertexDifferenceTgt3IronData", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedVsUnweightedZVertexDifferenceTgt3LeadData = new TH2D ("ANNWeightedVsUnweightedZVertexDifferenceTgt3LeadData", "ANNWeightedVsUnweightedZVertexDifferenceTgt3LeadData", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedVsUnweightedZVertexDifferenceTgt3CarbonData = new TH2D ("ANNWeightedVsUnweightedZVertexDifferenceTgt3CarbonData", "ANNWeightedVsUnweightedZVertexDifferenceTgt3CarbonData", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedVsUnweightedZVertexDifferenceTgt4LeadData = new TH2D ("ANNWeightedVsUnweightedZVertexDifferenceTgt4LeadData", "ANNWeightedVsUnweightedZVertexDifferenceTgt4LeadData", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedVsUnweightedZVertexDifferenceTgt5IronData = new TH2D ("ANNWeightedVsUnweightedZVertexDifferenceTgt5IronData", "ANNWeightedVsUnweightedZVertexDifferenceTgt5IronData", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedVsUnweightedZVertexDifferenceTgt5LeadData = new TH2D ("ANNWeightedVsUnweightedZVertexDifferenceTgt5LeadData", "ANNWeightedVsUnweightedZVertexDifferenceTgt5LeadData", 800, -400, 400, 101, 0, 100);
TH2D *ANNWeightedVsUnweightedZVertexDifferenceWaterData = new TH2D ("ANNWeightedVsUnweightedZVertexDifferenceWaterData", "ANNWeightedVsUnweightedZVertexDifferenceWaterData", 800, -400, 400, 101, 0, 100);

TH2D *ANNZResidualVsConfDifference = new TH2D ("ANNZResidualVsConfDifference", "ANNZResidualVsConfDifference", 100, -400, 400, 100, 0, 1);
TH2D *ANNZResidualVsConfDifferenceTracker = new TH2D ("ANNZResidualVsConfDifferenceTracker", "ANNZResidualVsConfDifferenceTracker", 100, -400, 400, 100, 0, 1);
TH2D *ANNZResidualVsConfDifferenceTgt2Iron = new TH2D ("ANNZResidualVsConfDifferenceTgt2Iron", "ANNZResidualVsConfDifferenceTgt2Iron", 100, -400, 400, 100, 0, 1);
TH2D *ANNZResidualVsConfDifferenceTgt2Lead = new TH2D ("ANNZResidualVsConfDifferenceTgt2Lead", "ANNZResidualVsConfDifferenceTgt2Lead", 100, -400, 400, 100, 0, 1);
TH2D *ANNZResidualVsConfDifferenceTgt3Iron = new TH2D ("ANNZResidualVsConfDifferenceTgt3Iron", "ANNZResidualVsConfDifferenceTgt3Iron", 100, -400, 400, 100, 0, 1);
TH2D *ANNZResidualVsConfDifferenceTgt3Lead = new TH2D ("ANNZResidualVsConfDifferenceTgt3Lead", "ANNZResidualVsConfDifferenceTgt3Lead", 100, -400, 400, 100, 0, 1);
TH2D *ANNZResidualVsConfDifferenceTgt3Carbon = new TH2D ("ANNZResidualVsConfDifferenceTgt3Carbon", "ANNZResidualVsConfDifferenceTgt3Carbon", 100, -400, 400, 100, 0, 1);
TH2D *ANNZResidualVsConfDifferenceTgt4Lead = new TH2D ("ANNZResidualVsConfDifferenceTgt4Lead", "ANNZResidualVsConfDifferenceTgt4Lead", 100, -400, 400, 100, 0, 1);
TH2D *ANNZResidualVsConfDifferenceTgt5Iron = new TH2D ("ANNZResidualVsConfDifferenceTgt5Iron", "ANNZResidualVsConfDifferenceTgt5Iron", 100, -400, 400, 100, 0, 1);
TH2D *ANNZResidualVsConfDifferenceTgt5Lead = new TH2D ("ANNZResidualVsConfDifferenceTgt5Lead", "ANNZResidualVsConfDifferenceTgt5Lead", 100, -400, 400, 100, 0, 1);
TH2D *ANNZResidualVsConfDifferenceWater = new TH2D ("ANNZResidualVsConfDifferenceWater", "ANNZResidualVsConfDifferenceWater", 100, -400, 400, 100, 0, 1);

TH2D *ANNPlaneProbabilityVsEhadData = new TH2D ("ANNPlaneProbabilityVsEhadData", "ANNPlaneProbabilityVsEhadData", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTrackerVsEhadData = new TH2D ("ANNPlaneProbabilityTrackerVsEhadData", "ANNPlaneProbabilityTrackerVsEhadData", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt2IronVsEhadData = new TH2D ("ANNPlaneProbabilityTgt2IronVsEhadData", "ANNPlaneProbabilityTgt2IronVsEhadData", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt2LeadVsEhadData = new TH2D ("ANNPlaneProbabilityTgt2LeadVsEhadData", "ANNPlaneProbabilityTgt2LeadVsEhadData", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt3IronVsEhadData = new TH2D ("ANNPlaneProbabilityTgt3IronVsEhadData", "ANNPlaneProbabilityTgt3IronVsEhadData", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt3LeadVsEhadData = new TH2D ("ANNPlaneProbabilityTgt3LeadVsEhadData", "ANNPlaneProbabilityTgt3LeadVsEhadData", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt3CarbonVsEhadData = new TH2D ("ANNPlaneProbabilityTgt3CarbonVsEhadData", "ANNPlaneProbabilityTgt3CarbonVsEhadData", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt4LeadVsEhadData = new TH2D ("ANNPlaneProbabilityTgt4LeadVsEhadData", "ANNPlaneProbabilityTgt4LeadVsEhadData", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt5IronVsEhadData = new TH2D ("ANNPlaneProbabilityTgt5IronVsEhadData", "ANNPlaneProbabilityTgt5IronVsEhadData", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt5LeadVsEhadData = new TH2D ("ANNPlaneProbabilityTgt5LeadVsEhadData", "ANNPlaneProbabilityTgt5LeadVsEhadData", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityWaterVsEhadData = new TH2D ("ANNPlaneProbabilityWaterVsEhadData", "ANNPlaneProbabilityWaterVsEhadData", 100, 0, 1, 100, 0, 20);

TH2D *ANNPlaneProbabilityVsEhadMC = new TH2D ("ANNPlaneProbabilityVsEhadMC", "ANNPlaneProbabilityVsEhadMC", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTrackerVsEhadMC = new TH2D ("ANNPlaneProbabilityTrackerVsEhadMC", "ANNPlaneProbabilityTrackerVsEhadMC", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt2IronVsEhadMC = new TH2D ("ANNPlaneProbabilityTgt2IronVsEhadMC", "ANNPlaneProbabilityTgt2IronVsEhadMC", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt2LeadVsEhadMC = new TH2D ("ANNPlaneProbabilityTgt2LeadVsEhadMC", "ANNPlaneProbabilityTgt2LeadVsEhadMC", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt3IronVsEhadMC = new TH2D ("ANNPlaneProbabilityTgt3IronVsEhadMC", "ANNPlaneProbabilityTgt3IronVsEhadMC", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt3LeadVsEhadMC = new TH2D ("ANNPlaneProbabilityTgt3LeadVsEhadMC", "ANNPlaneProbabilityTgt3LeadVsEhadMC", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt3CarbonVsEhadMC = new TH2D ("ANNPlaneProbabilityTgt3CarbonVsEhadMC", "ANNPlaneProbabilityTgt3CarbonVsEhadMC", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt4LeadVsEhadMC = new TH2D ("ANNPlaneProbabilityTgt4LeadVsEhadMC", "ANNPlaneProbabilityTgt4LeadVsEhadMC", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt5IronVsEhadMC = new TH2D ("ANNPlaneProbabilityTgt5IronVsEhadMC", "ANNPlaneProbabilityTgt5IronVsEhadMC", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt5LeadVsEhadMC = new TH2D ("ANNPlaneProbabilityTgt5LeadVsEhadMC", "ANNPlaneProbabilityTgt5LeadVsEhadMC", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityWaterVsEhadMC = new TH2D ("ANNPlaneProbabilityWaterVsEhadMC", "ANNPlaneProbabilityWaterVsEhadMC", 100, 0, 1, 100, 0, 20);

//ERecoil
TH1D *ErecoilMC = new TH1D ("ErecoilMC", "ErecoilMC", 1000, 0, 20);
TH1D *ErecoilData = new TH1D ("ErecoilData", "ErecoilData", 1000, 0, 20);
TH2D *TruthVerticesMCERecoil_ByModule = new TH2D ("TruthVerticesMCERecoil_ByModule", "TruthVerticesMCERecoil_ByModule", vertexBins.size()-1, &vertexBins[0], 1000, 0, 20);
TH2D *ANNVerticesMCERecoil_ByModule = new TH2D ("ANNVerticesMCERecoil_ByModule", "ANNVerticesMCERecoil_ByModule", vertexBins.size()-1, &vertexBins[0], 1000, 0, 20);
TH2D *TBVerticesMCERecoil_ByModule = new TH2D ("TBVerticesMCERecoil_ByModule", "TBVerticesMCERecoil_ByModule", vertexBins.size()-1, &vertexBins[0], 1000, 0, 20);
TH2D *ANNVerticesDataERecoil_ByModule = new TH2D ("ANNVerticesDataERecoil_ByModule", "ANNVerticesDataERecoil_ByModule", vertexBins.size()-1, &vertexBins[0], 1000, 0, 20);
TH2D *TBVerticesDataERecoil_ByModule = new TH2D ("TBVerticesDataERecoil_ByModule", "TBVerticesDataERecoil_ByModule", vertexBins.size()-1, &vertexBins[0], 1000, 0, 20);


TH2D *ANNVerticesData_ByZPosVsERecoil = new TH2D ("ANNVerticesData_ByZPosVsERecoil", "ANNVerticesData_ByZPosVsERecoil", 9000, 4200, 8700, 100, 0, 20);
TH2D *ANNVerticesMC_ByZPosVsERecoil = new TH2D ("ANNVerticesMC_ByZPosVsERecoil", "ANNVerticesMC_ByZPosVsERecoil", 9000, 4200, 8700, 100, 0, 20);


TH2D *ANNVerticesMC_BySegmentVsERecoil = new TH2D ("ANNVerticesMC_BySegmentVsERecoil", "ANNVerticesMC_BySegmentVsERecoil", 220, 0, 220, 1000, 0, 20);
TH2D *TruthVerticesMC_BySegmentVsERecoil = new TH2D ("TruthVerticesMC_BySegmentVsERecoil", "TruthVerticesMC_BySegmentVsERecoil", 220, 0, 220, 1000, 0, 20);
TH2D *ANNVerticesData_BySegmentVsERecoil = new TH2D ("ANNVerticesData_BySegmentVsERecoil", "ANNVerticesData_BySegmentVsERecoil", 220, 0, 220, 1000, 0, 20);


//Multiplicity
TH2D *TruthVerticesMCMultiplicity_ByModule = new TH2D ("TruthVerticesMCMultiplicity_ByModule", "TruthVerticesMCMultiplicity_ByModule", vertexBins.size()-1, &vertexBins[0], 10, 0, 10);
TH2D *ANNVerticesMCMultiplicity_ByModule = new TH2D ("ANNVerticesMCMultiplicity_ByModule", "ANNVerticesMCMultiplicity_ByModule", vertexBins.size()-1, &vertexBins[0], 10, 0, 10);
TH2D *TBVerticesMCMultiplicity_ByModule = new TH2D ("TBVerticesMCMultiplicity_ByModule", "TBVerticesMCMultiplicity_ByModule", vertexBins.size()-1, &vertexBins[0], 10, 0, 10);
TH2D *ANNVerticesDataMultiplicity_ByModule = new TH2D ("ANNVerticesDataMultiplicity_ByModule", "ANNVerticesDataMultiplicity_ByModule", vertexBins.size()-1, &vertexBins[0], 10, 0, 10);
TH2D *TBVerticesDataMultiplicity_ByModule = new TH2D ("TBVerticesDataMultiplicity_ByModule", "TBVerticesDataMultiplicity_ByModule", vertexBins.size()-1, &vertexBins[0], 10, 0, 10);

//ANN Confidence
TH2D *TruthVerticesMCANNConf_ByModule = new TH2D ("TruthVerticesMCANNConf_ByModule", "TruthVerticesMCANNConf_ByModule", vertexBins.size()-1, &vertexBins[0], 100, 0, 1);
TH2D *ANNVerticesMCANNConf_ByModule = new TH2D ("ANNVerticesMCANNConf_ByModule", "ANNVerticesMCANNConf_ByModule", vertexBins.size()-1, &vertexBins[0], 100, 0, 1);
TH2D *ANNVerticesDataANNConf_ByModule = new TH2D ("ANNVerticesDataANNConf_ByModule", "ANNVerticesDataANNConf_ByModule", vertexBins.size()-1, &vertexBins[0], 100, 0, 1);

TH2D *ANNVerticesData_ByZPosVsANNConf = new TH2D ("ANNVerticesData_ByZPosVsANNConf", "ANNVerticesData_ByZPosVsANNConf", 9000, 4200, 8700, 100, 0, 1);
TH2D *ANNVerticesMC_ByZPosVsANNConf = new TH2D ("ANNVerticesMC_ByZPosVsANNConf", "ANNVerticesMC_ByZPosVsANNConf", 9000, 4200, 8700, 100, 0, 1);

//For events misreconstructed in each target, what's the probability distribution
TH1D *ANNConfMisrecoInTgt2Iron = new TH1D ("ANNConfMisrecoInTgt2Iron", "ANNConfMisrecoInTgt2Iron", 100, 0, 1);
TH1D *ANNConfMisrecoInTgt2Lead = new TH1D ("ANNConfMisrecoInTgt2Lead", "ANNConfMisrecoInTgt2Lead", 100, 0, 1);
TH1D *ANNConfMisrecoInTgt3Iron = new TH1D ("ANNConfMisrecoInTgt3Iron", "ANNConfMisrecoInTgt3Iron", 100, 0, 1);
TH1D *ANNConfMisrecoInTgt3Lead = new TH1D ("ANNConfMisrecoInTgt3Lead", "ANNConfMisrecoInTgt3Lead", 100, 0, 1);
TH1D *ANNConfMisrecoInTgt3Carbon = new TH1D ("ANNConfMisrecoInTgt3Carbon", "ANNConfMisrecoInTgt3Carbon", 100, 0, 1);
TH1D *ANNConfMisrecoInTgt4Lead = new TH1D ("ANNConfMisrecoInTgt4Lead", "ANNConfMisrecoInTgt4Lead", 100, 0, 1);
TH1D *ANNConfMisrecoInTgt5Iron = new TH1D ("ANNConfMisrecoInTgt5Iron", "ANNConfMisrecoInTgt5Iron", 100, 0, 1);
TH1D *ANNConfMisrecoInTgt5Lead = new TH1D ("ANNConfMisrecoInTgt5Lead", "ANNConfMisrecoInTgt5Lead", 100, 0, 1);
TH1D *ANNConfMisrecoInWater = new TH1D ("ANNConfMisrecoInWater", "ANNConfMisrecoInWater", 100, 0, 1);

TH2D *ANNPlaneProbabilityVsPmuMC = new TH2D ("ANNPlaneProbabilityVsPmuMC", "ANNPlaneProbabilityVsPmuMC", 100, 0, 1, 50, 0, 50);
TH2D *ANNPlaneProbabilityVsPmuData = new TH2D ("ANNPlaneProbabilityVsPmuData", "ANNPlaneProbabilityVsPmuData", 100, 0, 1, 50, 0, 50);

//Efficiency as a function of ANNConfidence and EHad bin

TH2D *ANNPlaneProbabilityVsEhadNumerator = new TH2D ("ANNPlaneProbabilityVsEhadNumerator", "ANNPlaneProbabilityVsEhadNumerator", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTrackerVsEhadNumerator = new TH2D ("ANNPlaneProbabilityTrackerVsEhadNumerator", "ANNPlaneProbabilityTrackerVsEhadNumerator", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt2IronVsEhadNumerator = new TH2D ("ANNPlaneProbabilityTgt2IronVsEhadNumerator", "ANNPlaneProbabilityTgt2IronVsEhadNumerator", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt2LeadVsEhadNumerator = new TH2D ("ANNPlaneProbabilityTgt2LeadVsEhadNumerator", "ANNPlaneProbabilityTgt2LeadVsEhadNumerator", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt3IronVsEhadNumerator = new TH2D ("ANNPlaneProbabilityTgt3IronVsEhadNumerator", "ANNPlaneProbabilityTgt3IronVsEhadNumerator", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt3LeadVsEhadNumerator = new TH2D ("ANNPlaneProbabilityTgt3LeadVsEhadNumerator", "ANNPlaneProbabilityTgt3LeadVsEhadNumerator", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt3CarbonVsEhadNumerator = new TH2D ("ANNPlaneProbabilityTgt3CarbonVsEhadNumerator", "ANNPlaneProbabilityTgt3CarbonVsEhadNumerator", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt4LeadVsEhadNumerator = new TH2D ("ANNPlaneProbabilityTgt4LeadVsEhadNumerator", "ANNPlaneProbabilityTgt4LeadVsEhadNumerator", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt5IronVsEhadNumerator = new TH2D ("ANNPlaneProbabilityTgt5IronVsEhadNumerator", "ANNPlaneProbabilityTgt5IronVsEhadNumerator", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt5LeadVsEhadNumerator = new TH2D ("ANNPlaneProbabilityTgt5LeadVsEhadNumerator", "ANNPlaneProbabilityTgt5LeadVsEhadNumerator", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityWaterVsEhadNumerator = new TH2D ("ANNPlaneProbabilityWaterVsEhadNumerator", "ANNPlaneProbabilityWaterVsEhadNumerator", 100, 0, 1, 100, 0, 20);

TH2D *ANNPlaneProbabilityVsEhadDenominator = new TH2D ("ANNPlaneProbabilityVsEhadDenominator", "ANNPlaneProbabilityVsEhadDenominator", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTrackerVsEhadDenominator = new TH2D ("ANNPlaneProbabilityTrackerVsEhadDenominator", "ANNPlaneProbabilityTrackerVsEhadDenominator", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt2IronVsEhadDenominator = new TH2D ("ANNPlaneProbabilityTgt2IronVsEhadDenominator", "ANNPlaneProbabilityTgt2IronVsEhadDenominator", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt2LeadVsEhadDenominator = new TH2D ("ANNPlaneProbabilityTgt2LeadVsEhadDenominator", "ANNPlaneProbabilityTgt2LeadVsEhadDenominator", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt3IronVsEhadDenominator = new TH2D ("ANNPlaneProbabilityTgt3IronVsEhadDenominator", "ANNPlaneProbabilityTgt3IronVsEhadDenominator", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt3LeadVsEhadDenominator = new TH2D ("ANNPlaneProbabilityTgt3LeadVsEhadDenominator", "ANNPlaneProbabilityTgt3LeadVsEhadDenominator", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt3CarbonVsEhadDenominator = new TH2D ("ANNPlaneProbabilityTgt3CarbonVsEhadDenominator", "ANNPlaneProbabilityTgt3CarbonVsEhadDenominator", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt4LeadVsEhadDenominator = new TH2D ("ANNPlaneProbabilityTgt4LeadVsEhadDenominator", "ANNPlaneProbabilityTgt4LeadVsEhadDenominator", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt5IronVsEhadDenominator = new TH2D ("ANNPlaneProbabilityTgt5IronVsEhadDenominator", "ANNPlaneProbabilityTgt5IronVsEhadDenominator", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityTgt5LeadVsEhadDenominator = new TH2D ("ANNPlaneProbabilityTgt5LeadVsEhadDenominator", "ANNPlaneProbabilityTgt5LeadVsEhadDenominator", 100, 0, 1, 100, 0, 20);
TH2D *ANNPlaneProbabilityWaterVsEhadDenominator = new TH2D ("ANNPlaneProbabilityWaterVsEhadDenominator", "ANNPlaneProbabilityWaterVsEhadDenominator", 100, 0, 1, 100, 0, 20);


//Curvature signifiance
TH2D *TruthVerticesMCCurvSig_ByModule = new TH2D ("TruthVerticesMCCurvSig_ByModule", "TruthVerticesMCCurvSig_ByModule", vertexBins.size()-1, &vertexBins[0], 100, 0, 10);
TH2D *ANNVerticesMCCurvSig_ByModule = new TH2D ("ANNVerticesMCCurvSig_ByModule", "ANNVerticesMCCurvSig_ByModule", vertexBins.size()-1, &vertexBins[0], 100, 0, 10);
TH2D *TBVerticesMCCurvSig_ByModule = new TH2D ("TBVerticesMCCurvSig_ByModule", "TBVerticesMCCurvSig_ByModule", vertexBins.size()-1, &vertexBins[0], 100, 0, 10);
TH2D *ANNVerticesDataCurvSig_ByModule = new TH2D ("ANNVerticesDataCurvSig_ByModule", "ANNVerticesDataCurvSig_ByModule", vertexBins.size()-1, &vertexBins[0], 100, 0, 10);
TH2D *TBVerticesDataCurvSig_ByModule = new TH2D ("TBVerticesDataCurvSig_ByModule", "TBVerticesDataCurvSig_ByModule", vertexBins.size()-1, &vertexBins[0], 100, 0, 10);

//Migration ANN only
TH2D *ANNVerticesConfusion_ByModule = new TH2D ("ANNVerticesConfusion_ByModule", "ANNVerticesConfusion_ByModule", vertexBins.size()-1, &vertexBins[0], vertexBins.size()-1, &vertexBins[0]);
TH2D *ANNVerticesConfusion_ByZPos = new TH2D ("ANNVerticesConfusion_ByZPos", "ANNVerticesConfusion_ByZPos", 4500, 4200, 8700, 4500, 4200, 8700);

//==============================================================================
// End - Declaring histograms
//==============================================================================

int nuOrAntiNuMode(std::string playlist)
{
    std::vector<std::string> nuVector = {"minervame1A", "minervame1B", "minervame1C", "minervame1D", "minervame1E", "minervame1F", "minervame1G", "minervame1L", "minervame1M", "minervame1N", "minervame1O", "minervame1P"};
    std::vector<std::string> anuVector = {"minervame5A", "minervame6A", "minervame6B", "minervame6C", "minervame6D", "minervame6E", "minervame6F", "minervame6G", "minervame6H", "minervame6I", "minervame6J"};
    if (std::find(std::begin(nuVector), std::end(nuVector), playlist) != std::end(nuVector)) return 1;
    else if (std::find(std::begin(anuVector), std::end(anuVector), playlist) != std::end(anuVector)) return 2;
    //Is this the most efficient way? Probably not
    return 0;
}

//==============================================================================
// Plotting Functions
//==============================================================================
void PlotRegions(double ymin, double ymax, double mcpot_input, double datapot_input)
//Add option to draw by module
//Add option to draw immediate US and DS scintillator
{
        //These were all taken from https://cdcvs.fnal.gov/redmine/projects/minerva/wiki/Z_Positions_of_Planes_in_the_Full_MINERvA_Detector
    //For target 3 the centre was taken as the mid point between the US and DS planes, the stated centre gives weird looking results I suspect because it's a 3 material target the reported centre may not be actual centre
    double tgt1centre = 4481; 
    double tgt1thickness = 25.78; 
    double tgt2centre = 4702; 
    double tgt2thickness = 25.81; 
    double tgt3centre = 4944.48; 
    double tgt3thickness = 76.2; 
    double tgt4centre = 5641.68; 
    double tgt4thickness = 7.95; 
    double tgt5centre = 5778; 
    double tgt5thickness = 13.17; 

    TBox *tgt1 = new TBox(tgt1centre-tgt1thickness/2, ymin, tgt1centre+tgt1thickness/2, ymax);
    TBox *tgt2 = new TBox(tgt2centre-tgt2thickness/2, ymin, tgt2centre+tgt2thickness/2, ymax);
    TBox *tgt3 = new TBox(tgt3centre-tgt3thickness/2, ymin, tgt3centre+tgt3thickness/2, ymax);
    TBox *tgt4 = new TBox(tgt4centre-tgt4thickness/2, ymin, tgt4centre+tgt4thickness/2, ymax);
    TBox *tgt5 = new TBox(tgt5centre-tgt5thickness/2, ymin, tgt5centre+tgt5thickness/2, ymax);
    TBox *tgt6 = new TBox(5200, ymin, 5420, ymax);

    std::stringstream mcpot, datapot;
    mcpot << std::fixed << std::scientific << std::setprecision(3) << mcpot_input;
    datapot << std::fixed << std::scientific << std::setprecision(3) << datapot_input;
    std::string mcpotstr = std::string("MC POT: ")+mcpot.str();
    std::string datapotstr = std::string("Data POT: ")+datapot.str();

    std::string potstr = std::string("Normalised POT: ")+datapot.str(); 
    /* TLatex filledpotlabel = TLatex( 4500 ,1000,filledpotstr.c_str());
    filledpotlabel.SetTextAlign(22);
    filledpotlabel.SetTextFont(43);
    filledpotlabel.SetTextSize(20); 
    //xaxlabel.SetTextFont(22);
    TLatex emptypotlabel = TLatex( 4500 ,1000,emptypotstr.c_str());
    emptypotlabel.SetTextAlign(22);
    emptypotlabel.SetTextFont(43);
    emptypotlabel.SetTextSize(20); 
    //xaxlabel.SetTextFont(22); */

    TLatex potlabel = TLatex( xlow , ymax*1.04,potstr.c_str());
    potlabel.SetTextAlign(12);
    potlabel.SetTextFont(43);
    potlabel.SetTextSize(20);

    double resolution = ((xhigh - xlow )/ nbins) * rebinNum;
    std::stringstream resolutionStr;
    resolutionStr << std::fixed << std::setprecision(1) << resolution;
    std::string resStr = std::string("Resolution/Bin Width: ")+resolutionStr.str()+ std::string(" mm");
    TLatex resLabel = TLatex( xhigh , ymax*1.04,resStr.c_str());
    resLabel.SetTextAlign(32);
    resLabel.SetTextFont(43);
    resLabel.SetTextSize(20); 
    resLabel.Draw("SAME");
    
    TLatex titlelabel = TLatex( 5000 , ymax*1.08, (std::string("True Vertex Positions of Interactions reconstructed inside the ")).c_str());

    titlelabel.SetTextAlign(22);
    titlelabel.SetTextFont(43);
    titlelabel.SetTextSize(35); 

    //filledpotlabel.Draw("SAME");
    //emptypotlabel.Draw("SAME");
    potlabel.Draw("SAME");
    titlelabel.Draw("same");

    tgt1->SetFillColorAlpha(kBlack,0.5);
    tgt1->SetFillStyle(3004);
    tgt1->SetLineWidth(0);
    tgt1->SetLineColor(1);
    tgt1->Draw("same");
    tgt2->SetFillColorAlpha(kBlack,0.5);
    tgt2->SetFillStyle(3004);
    tgt2->SetLineWidth(0);
    tgt2->SetLineColor(1);
    tgt2->Draw("same");
    tgt3->SetFillColorAlpha(kBlack,0.5);
    tgt3->SetFillStyle(3004);
    tgt3->SetLineWidth(0);
    tgt3->SetLineColor(1);
    tgt3->Draw("same");
    tgt4->SetFillColorAlpha(kBlack,0.5);
    tgt4->SetFillStyle(3004);
    tgt4->SetLineWidth(0);
    tgt4->SetLineColor(1);
    tgt4->Draw("same");
    tgt5->SetFillColorAlpha(kBlack,0.5);
    tgt5->SetFillStyle(3004);
    tgt5->SetLineWidth(0);
    tgt5->SetLineColor(1);
    tgt5->Draw("same");
    tgt6->SetFillColorAlpha(kBlue,0.5);
    tgt6->SetFillStyle(3004);
    tgt6->SetLineWidth(0);
    tgt6->SetLineColor(1);
    tgt6->Draw("same");
}


void DrawDataMC(TH1D* DataIn, TH1D* MCIn, std::string name, double mcpot, double datapot)
{
        TCanvas *canvas = new TCanvas("canvas", "canvas", 1500, 1000);
        TH1D *DataHist = (TH1D*)DataIn->Clone();
        TH1D *MCHist = (TH1D*)MCIn->Clone();
        MCHist->SetTitle("MC");
        DataHist->SetTitle("Data");
        
        //std::cout<<"name: " <<name  <<std::endl;
        name.erase(std::remove(name.begin(), name.end(), ' '), name.end()); //remove spaces
        //std::cout<<"name after spaces: " <<name  <<std::endl;
        name.erase(std::remove(name.begin(), name.end(), '-'), name.end()); //remove dashes
        //std::cout<<"name after dashes: " <<name  <<std::endl;

        PlotUtils::MnvPlotter plotter;
        plotter.ApplyStyle(PlotUtils::kDefaultStyle);
        //plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
        plotter.DrawDataMC(DataHist, MCHist, 1.0, "TL", true);
        double ymin = gPad->GetFrame()->GetY1();
        double ymax = gPad->GetFrame()->GetY2();
        PlotRegions(ymin, ymax, mcpot, datapot);
        std::cout<<"Test0\n";
        canvas->SaveAs((std::string("/pnfs/minerva/persistent/users/alhart/31JulTemp/")+name+std::string(".png")).c_str());
        std::cout<<"Test1\n";
        canvas->Clear();
        std::cout<<"Test2\n";
        //canvas->SetLogy(0);
        plotter.DrawDataMCRatio(DataHist, MCHist, 1.0, true, -1.0, -1.0, "Data MC Ratio");
        std::cout<<"Test3\n";
        ymin = gPad->GetFrame()->GetY1();
        std::cout<<"Test4\n";
        ymax = gPad->GetFrame()->GetY2();
        std::cout<<"Test5\n";
        PlotRegions(ymin, ymax, mcpot, datapot);
        std::cout<<"Test6\n";
        canvas->SaveAs((std::string("/pnfs/minerva/persistent/users/alhart/31JulTemp/")+name+std::string("Ratio.png")).c_str());

        canvas->Clear();
        canvas->SetLogy();
        DataHist->SetMinimum(10);
        MCHist->SetMinimum(10);
        std::cout<< "Literally just set min to 1\n";
        plotter.DrawDataMC(DataHist, MCHist, 1.0, "TL", true);
        ymin = gPad->GetFrame()->GetY1();
        ymax = gPad->GetFrame()->GetY2();
        PlotRegions(std::pow(10,ymin), std::pow(10,ymax), mcpot, datapot);
        std::cout<<"LogY set: " << gPad->GetLogy() << std::endl;
        TH1D* tmpMC = (TH1D*)canvas->GetPrimitive("tmpMC"); //The name of the MC histogram from DrawDataMC, which is drawn first
        std::cout << "minimum: " << tmpMC->GetMinimum() << std::endl; //For log plot we start at 1 not 0
        canvas->SaveAs((std::string("/exp/minerva/app/users/alhart/MAT_AL9/AnalysisPlotScripts/")+name+std::string("Log.png")).c_str());
}



void DrawFullAndEmptyDataMC(TH1D* DataInFull, TH1D* MCInFull, TH1D* DataInEmpty, TH1D* MCInEmpty, std::string name, double mcpotfull, double datapotfull, double mcpotempty, double datapotempty)
{
        TCanvas *canvas = new TCanvas("canvas", "canvas", 1500, 1000);

        //TH1D *DataHist = (TH1D*)DataIn->Clone();
        //TH1D *MCHist = (TH1D*)MCIn->Clone();
        MCInFull->SetTitle("MC");
        DataInFull->SetTitle("Data");
        MCInEmpty->SetTitle("MC");
        DataInEmpty->SetTitle("Data");
        
        PlotUtils::MnvPlotter plotter;
        plotter.ApplyStyle(PlotUtils::kCCQENuStyle);

        //Draw filled linear plot
        plotter.DrawDataMC(DataInFull, MCInFull, 1.0, "TL", true);
        double ymin = gPad->GetFrame()->GetY1();
        double ymax = gPad->GetFrame()->GetY2();
        PlotRegions(ymin, ymax, 0, 0);
        canvas->SaveAs((std::string("/exp/minerva/app/users/alhart/MAT_AL9/AnalysisPlotScripts/TargetRegionValidationPlotsQ4/")+name+std::string("Filled.png")).c_str());
        
        //Draw empty linear plot
        canvas->Clear();
        plotter.DrawDataMC(DataInEmpty, MCInEmpty, 1.0, "TL", true);
        ymin = gPad->GetFrame()->GetY1();
        ymax = gPad->GetFrame()->GetY2();
        PlotRegions(ymin, ymax, 0, 0);
        canvas->SaveAs((std::string("/exp/minerva/app/users/alhart/MAT_AL9/AnalysisPlotScripts/TargetRegionValidationPlotsQ4/")+name+std::string("Empty.png")).c_str());


        //Draw filled ratio plot
        canvas->Clear();
        plotter.DrawDataMCRatio(DataInFull, MCInFull, 1.0, true, -1.0, -1.0, "Data MC Ratio");
        ymin = gPad->GetFrame()->GetY1();
        ymax = gPad->GetFrame()->GetY2();
        PlotRegions(ymin, ymax, 0, 0);
        canvas->SaveAs((std::string("/exp/minerva/app/users/alhart/MAT_AL9/AnalysisPlotScripts/TargetRegionValidationPlotsQ4/")+name+std::string("FilledRatio.png")).c_str());


        //Draw empty ratio plot
        canvas->Clear();
        plotter.DrawDataMCRatio(DataInEmpty, MCInEmpty, 1.0, true, -1.0, -1.0, "Data MC Ratio");
        ymin = gPad->GetFrame()->GetY1();
        ymax = gPad->GetFrame()->GetY2();
        PlotRegions(ymin, ymax, 0, 0);
        canvas->SaveAs((std::string("/exp/minerva/app/users/alhart/MAT_AL9/AnalysisPlotScripts/TargetRegionValidationPlotsQ4/")+name+std::string("FilledRatio.png")).c_str());

        //Draw filled log plot
        canvas->Clear();
        canvas->SetLogy();
        DataInFull->SetMinimum(1);
        MCInFull->SetMinimum(1);
        plotter.DrawDataMC(DataInFull, MCInFull, 1.0, "TL", true);
        ymin = gPad->GetFrame()->GetY1();
        ymax = gPad->GetFrame()->GetY2();
        PlotRegions(std::pow(10,ymin), std::pow(10,ymax), 0, 0);
        canvas->SaveAs((std::string("/exp/minerva/app/users/alhart/MAT_AL9/AnalysisPlotScripts/TargetRegionValidationPlotsQ4/")+name+std::string("FilledLog.png")).c_str());

        //Draw empty log plot
        canvas->Clear();
        canvas->SetLogy();
        DataInEmpty->SetMinimum(1);
        MCInEmpty->SetMinimum(1);
        plotter.DrawDataMC(DataInEmpty, MCInEmpty, 1.0, "TL", true);
        ymin = gPad->GetFrame()->GetY1();
        ymax = gPad->GetFrame()->GetY2();
        PlotRegions(std::pow(10,ymin), std::pow(10,ymax), 0, 0);
        canvas->SaveAs((std::string("/exp/minerva/app/users/alhart/MAT_AL9/AnalysisPlotScripts/TargetRegionValidationPlotsQ4/")+name+std::string("EmptyLog.png")).c_str());


        //Drawing split plots, with both filled and empty in the same plot
        //Preparing canvas
        canvas->Clear();
        canvas->SetLogy(false); //Unsetting log axis
        TPad *p1 = new TPad("p1","top",0,0.45,1,0.9); // xlow,yloy, xup,yup
        TPad *p2 = new TPad("p2","bottom",0,0,1,0.45);
        auto tex3 = new TLatex(.5,.95,"Title2");
        tex3->SetTextAlign(22);
        tex3->SetTextSize(0.025);
        std::string title = name + " Z vertex distributions - water filled (top) and empty (bottom) for 1.57<ERecoil<50";
        tex3->DrawLatexNDC(.5,.95, title.c_str());
        //Draw linear plots
        p1->Draw();
        p1->cd();
        //Draw filled linear plot
        plotter.DrawDataMC(DataInFull, MCInFull, 1.0, "TL", true);
        ymin = gPad->GetFrame()->GetY1();
        ymax = gPad->GetFrame()->GetY2();
        PlotRegions(ymin, ymax, 0, 0);
        canvas->cd();
        p2->Draw();
        p2->cd();
        //Draw empty linear plot
        plotter.DrawDataMC(DataInEmpty, MCInEmpty, 1.0, "TL", true);
        ymin = gPad->GetFrame()->GetY1();
        ymax = gPad->GetFrame()->GetY2();
        PlotRegions(ymin, ymax, 0, 0);
        canvas->cd();
        canvas->SaveAs((std::string("/exp/minerva/app/users/alhart/MAT_AL9/AnalysisPlotScripts/TargetRegionValidationPlotsQ4/")+name+std::string("Split.png")).c_str());


        //Draw ratio plots
        canvas->Clear();
        p1 = new TPad("p1","top",0,0.45,1,0.9); // xlow,yloy, xup,yup
        p2 = new TPad("p2","bottom",0,0,1,0.45);

        tex3->DrawLatexNDC(.5,.95, title.c_str());
        //Draw filled ratio plot
        p1->Draw();
        p1->cd();
        plotter.DrawDataMCRatio(DataInFull, MCInFull, 1.0, true, -1.0, -1.0, "Data MC Ratio");
        ymin = gPad->GetFrame()->GetY1();
        ymax = gPad->GetFrame()->GetY2();
        PlotRegions(ymin, ymax, 0, 0);
        auto tex = new TLatex(0.5,0.5,"Title");
        tex->SetTextAlign(22);
        tex->SetTextSize(0.25);
        tex->Draw("SAME");
        //Draw empty ratio plot
        canvas->cd();
        p2->Draw();
        p2->cd();
        plotter.DrawDataMCRatio(DataInEmpty, MCInEmpty, 1.0, true, -1.0, -1.0, "Data MC Ratio");
        ymin = gPad->GetFrame()->GetY1();
        ymax = gPad->GetFrame()->GetY2();
        PlotRegions(ymin, ymax, 0, 0);
        auto tex2 = new TLatex(0.5,0.5,"Title1");
        tex2->SetTextAlign(22);
        tex2->SetTextSize(0.25);
        tex2->Draw("SAME");
        canvas->cd();

        canvas->SaveAs((std::string("/exp/minerva/app/users/alhart/MAT_AL9/AnalysisPlotScripts/TargetRegionValidationPlotsQ4/")+name+std::string("RatioSplit.png")).c_str());




        //Draw log plots
        canvas->Clear();
        p1 = new TPad("p1","top",0,0.45,1,0.9); // xlow,yloy, xup,yup
        p2 = new TPad("p2","bottom",0,0,1,0.45);
        tex3->DrawLatexNDC(.5,.95, title.c_str());

        //Draw filled log plot
        p1->Draw();
        p1->cd();
        p1->SetLogy();
        DataInFull->SetMinimum(1);
        MCInFull->SetMinimum(1);
        plotter.DrawDataMC(DataInFull, MCInFull, 1.0, "TL", true);
        ymin = gPad->GetFrame()->GetY1();
        ymax = gPad->GetFrame()->GetY2();
        std::cout << "Filled ymin: " << ymin << " ymax: " << ymax << std::endl;
        PlotRegions(std::pow(10,ymin), std::pow(10,ymax), 0, 0);
        //Draw empty log plot
        canvas->cd();
        p2->Draw();
        p2->cd();
        p2->SetLogy();
        //canvas->SetLogy();
        DataInEmpty->SetMinimum(1);
        MCInEmpty->SetMinimum(1);
        plotter.DrawDataMC(DataInEmpty, MCInEmpty, 1.0, "TL", true);
        ymin = gPad->GetFrame()->GetY1();
        ymax = gPad->GetFrame()->GetY2();
        std::cout << "Empty ymin: " << ymin << " ymax: " << ymax << std::endl;
        PlotRegions(std::pow(10,ymin), std::pow(10,ymax), 0, 0);
        canvas->cd();
        canvas->SaveAs((std::string("/exp/minerva/app/users/alhart/MAT_AL9/AnalysisPlotScripts/TargetRegionValidationPlotsQ4/")+name+std::string("LogSplit.png")).c_str());
}

//==============================================================================
// End - Plotting Functions
//==============================================================================



//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFillMC(
    PlotUtils::ChainWrapper* chain,
    std::map<std::string, std::vector<CVUniverse*> > error_bands,
    PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts,
    PlotUtils::Model<CVUniverse, MichelEvent>& model
    )
{
  std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
  assert(!error_bands["cv"].empty() && "\"cv\" error band is empty!  Can't set Model weight.");
  auto& cvUniv = error_bands["cv"].front();

  std::cout << "Starting MC reco loop...\n";

  const int nEntries = chain->GetEntries();
  //const int nEntries = 10000;
  int  startNum = 0;
  int endNum = nEntries;
  if (nSubruns != 0 )
  {
    int chunkSize = std::floor(nEntries/nSubruns);
    startNum = nProcess*chunkSize;
    if (nProcess == (nSubruns-1)) //The last one, run to the end
    {
      endNum = nEntries; //probably redundant with the same line just above
    }
    else
    {
      endNum = (nProcess+1)*chunkSize;
    }
  }
  // const int nEntries = 10000;
  for (int i = startNum; i < endNum; ++i)
  {
    if(i%1000==0)
    { 
      std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
      double mselapsed = std::chrono::duration_cast<std::chrono::milliseconds> (now - startTime).count();
      double msperevent = mselapsed/i;
      int sleft = (nEntries-i)*msperevent/1000;
      double selapsed = ((int)mselapsed)/1000;
      std::cout << i << " / " << nEntries << ". Time Elapsed: " << selapsed <<" Est Remaining: " << sleft <<"\r" <<std::flush;
    }
    MichelEvent cvEvent;
    cvUniv->SetEntry(i);
    model.SetEntry(*cvUniv, cvEvent);

    const double cvWeight = model.GetWeight(*cvUniv, cvEvent);

    //=========================================
    // Systematics loop(s)
    //=========================================
    for (auto band : error_bands)
    {
      std::vector<CVUniverse*> error_band_universes = band.second;
      int univCount = 0;
      for (auto universe : error_band_universes)
      {
        univCount++; // Put the iterator right at the start so it's executed even in paths that lead to a continue, don't forget to subtract by 1 when we use it
        MichelEvent myevent; // make sure your event is inside the error band loop. 
        // Tell the Event which entry in the TChain it's looking at
        universe->SetEntry(i);
        ROOT::Math::XYZVector ANNVtx = universe->GetANNVertex();
        ROOT::Math::XYZTVector TrackBasedVtx = universe->GetVertex();
        ROOT::Math::XYZTVector TrueVtx = universe->GetTrueVertex();


        PlotUtils::TargetUtils* m_TargetUtils=new PlotUtils::TargetUtils();
        const double weight = model.GetWeight(*universe, myevent); //Only calculate the per-universe weight for events that will actually use it.
        double ANNProb = universe->GetANNProb();
        double erecoil = universe->GetRecoilE()/pow(10,3);
        double curvsig = 1/universe->GetMuonQPErr();
        double pmu = universe->GetPmu()/1000;

        // This is where you would Access/create a Michel

        //weight is ignored in isMCSelected() for all but the CV Universe.
        if (!michelcuts.isMCSelected(*universe, myevent, cvWeight).all()) continue;   //All cuts except ANN confidence cut

        if(ANNVtx.X()!=-1.0)
        {
          ANNPlaneProbabilityVsEhadMC->Fill(ANNProb, erecoil, cvWeight);
          if (m_TargetUtils->InIron2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNPlaneProbabilityTgt2IronVsEhadMC->Fill(ANNProb, erecoil, cvWeight);
          if (m_TargetUtils->InLead2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNPlaneProbabilityTgt2LeadVsEhadMC->Fill(ANNProb, erecoil, cvWeight);
          if (m_TargetUtils->InIron3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNPlaneProbabilityTgt3IronVsEhadMC->Fill(ANNProb, erecoil, cvWeight);
          if (m_TargetUtils->InLead3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNPlaneProbabilityTgt3LeadVsEhadMC->Fill(ANNProb, erecoil, cvWeight);
          if (m_TargetUtils->InCarbon3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNPlaneProbabilityTgt3CarbonVsEhadMC->Fill(ANNProb, erecoil, cvWeight);
          if (m_TargetUtils->InLead4VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNPlaneProbabilityTgt4LeadVsEhadMC->Fill(ANNProb, erecoil, cvWeight);
          if (m_TargetUtils->InIron5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNPlaneProbabilityTgt5IronVsEhadMC->Fill(ANNProb, erecoil, cvWeight);
          if (m_TargetUtils->InLead5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNPlaneProbabilityTgt5LeadVsEhadMC->Fill(ANNProb, erecoil, cvWeight);
          if (m_TargetUtils->InWaterTargetVolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNPlaneProbabilityWaterVsEhadMC->Fill(ANNProb, erecoil, cvWeight);
          if (m_TargetUtils->InTracker(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNPlaneProbabilityTrackerVsEhadMC->Fill(ANNProb, erecoil, cvWeight);
          ANNPlaneProbabilityVsPmuMC->Fill(ANNProb, pmu, cvWeight);
          ANNPlaneProbabilityVsEhadMC->Fill(ANNProb, erecoil, cvWeight);

          ANNVerticesMC_ByZPosVsANNConf->Fill(ANNVtx.Z(), ANNProb, cvWeight);
          ANNVerticesMCANNConf_ByModule->Fill( ANNVtx.Z(), ANNProb, cvWeight);

          //For events misreconstructed in each target, what's the probability distribution
          if (m_TargetUtils->InIron2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true) && !m_TargetUtils->InIron2VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true)) ANNConfMisrecoInTgt2Iron->Fill(ANNProb);
          if (m_TargetUtils->InLead2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true) && !m_TargetUtils->InLead2VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true)) ANNConfMisrecoInTgt2Lead->Fill(ANNProb);
          if (m_TargetUtils->InIron3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true) && !m_TargetUtils->InIron3VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true)) ANNConfMisrecoInTgt3Iron->Fill(ANNProb);
          if (m_TargetUtils->InLead3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true) && !m_TargetUtils->InLead3VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true)) ANNConfMisrecoInTgt3Lead->Fill(ANNProb);
          if (m_TargetUtils->InCarbon3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true) && !m_TargetUtils->InCarbon3VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true)) ANNConfMisrecoInTgt3Carbon->Fill(ANNProb);
          if (m_TargetUtils->InLead4VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()) && !m_TargetUtils->InLead4VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z())) ANNConfMisrecoInTgt4Lead->Fill(ANNProb);
          if (m_TargetUtils->InIron5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true) && !m_TargetUtils->InIron5VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true)) ANNConfMisrecoInTgt5Iron->Fill(ANNProb);
          if (m_TargetUtils->InLead5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true) && !m_TargetUtils->InLead5VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true)) ANNConfMisrecoInTgt5Lead->Fill(ANNProb);
          if (m_TargetUtils->InWaterTargetVolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()) && !m_TargetUtils->InWaterTargetVolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z())) ANNConfMisrecoInWater->Fill(ANNProb);

        }
        TBVerticesMCCurvSig_ByModule->Fill( TrackBasedVtx.Z(), curvsig, cvWeight);
        TruthVerticesMCANNConf_ByModule->Fill( TrackBasedVtx.Z(), ANNProb, cvWeight);
        TruthVerticesMCCurvSig_ByModule->Fill( TrackBasedVtx.Z(), curvsig, cvWeight);

        if (!ANNConfCut->passesCut(*universe, myevent, cvWeight)) continue;   //ANN confidence cut

        //Performing vtx validation check Deborah suggested
        double batchPOT = cvUniv->GetBatchPOT();
        double efficiency = 0.5563 - (0.01353*batchPOT); //Based on MINERvA-doc-21436
        //cvWeight/=efficiency; //Not doing efficiency correction right now

        //Hadron Energy Spectrum plot
        int multiplicity = cvUniv->GetMultiplicity();
        //-----------------------------------------------------------
        // Making a note of events that may be interesting to view in Arachne
        //-----------------------------------------------------------
        //if (m_TargetUtils->InWaterTargetVolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()) && (m_TargetUtils->InWaterTargetVolMC(TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z())))  std::cout<<"Check in arachne: \n" << " ev_run: " << cvUniv->GetInt("ev_run") << " ev_subrun: " << cvUniv->GetInt("ev_subrun") << " ev_gate: " << cvUniv->GetInt("ev_gate") << std::endl;
        if ( std::abs(TrueVtx.Z() - ANNVtx.Z()) > 2000 && ANNProb>0.4 && verbose)
        {
          std::cout<<"Difference between true and ann reco z larger than 2000: \tANN Segment: " << universe->GetANNSegment() << "\tTrue Segment: " << universe->GetTruthSegment() << " mc_run: " << universe->GetInt("mc_run") << " mc_subrun: " << universe->GetInt("mc_subrun")  << " mc_nthEvtInSpill: " << universe->GetInt("mc_nthEvtInSpill") << " mc_nthEvtInFile: " << universe->GetInt("mc_nthEvtInFile") << " ev_global_gate: " << universe->GetInt("ev_global_gate")  << std::endl;
          std::string ArachneLink = "https://mnvevdgpvm02.fnal.gov/Arachne/?det=SIM_minerva&recoVer=v22r1p1&run="+std::to_string(universe->GetInt("mc_run"))+"&subrun="+std::to_string(universe->GetInt("mc_subrun"))+"&gate="+std::to_string(universe->GetInt("mc_nthEvtInFile")+1)+"&slice=-1";
          std::cout<<"Arachne Link: " << ArachneLink << std::endl;
        }

        if(ANNVtx.X()!=-1.0)
        {
          ANNVerticesMC_ByModule->Fill( ANNVtx.Z(), cvWeight);
          ANNVerticesMC_ByZPos->Fill( ANNVtx.Z(), cvWeight);
          ANNVerticesMC_ByZPosVsERecoil->Fill(ANNVtx.Z(), erecoil, cvWeight);

          ANNVerticesMCERecoil_ByModule->Fill( ANNVtx.Z(), erecoil, cvWeight);
          ANNVerticesMCMultiplicity_ByModule->Fill( ANNVtx.Z(), multiplicity, cvWeight);
          ANNVerticesMCCurvSig_ByModule->Fill( ANNVtx.Z(), curvsig, cvWeight);

          ANNVerticesConfusion_ByModule->Fill(ANNVtx.Z(), TrueVtx.Z(), cvWeight);
          ANNVerticesConfusion_ByZPos->Fill(ANNVtx.Z(), TrueVtx.Z(), cvWeight);

          ANNVerticesMC_BySegment->Fill(universe->GetANNSegment(), cvWeight);
          ANNVerticesMC_BySegmentVsERecoil->Fill(universe->GetANNSegment(), erecoil, cvWeight);
          //if (universe->GetANNSegment() == 45) std::cout<<"Found plane in ANN mod 45, truth seg: " << universe->GetTruthSegment() << " truth mod: " << universe->GetTruthVtxModule()<< " truth pla: " << universe->GetTruthVtxPlane()<< " truth tgtid: " << universe->GetTruthTargetID()<< " truth tgtz: " << universe->GetTruthTargetZ()<<  std::endl;

          ANNZVertexResidual->Fill((ANNVtx.Z() - TrueVtx.Z()), cvWeight);
          if (m_TargetUtils->InIron2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZVertexResidualTgt2Iron->Fill((ANNVtx.Z() - TrueVtx.Z()), cvWeight);
          if (m_TargetUtils->InLead2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZVertexResidualTgt2Lead->Fill((ANNVtx.Z() - TrueVtx.Z()), cvWeight);
          if (m_TargetUtils->InIron3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZVertexResidualTgt3Iron->Fill((ANNVtx.Z() - TrueVtx.Z()), cvWeight);
          if (m_TargetUtils->InLead3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZVertexResidualTgt3Lead->Fill((ANNVtx.Z() - TrueVtx.Z()), cvWeight);
          if (m_TargetUtils->InCarbon3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZVertexResidualTgt3Carbon->Fill((ANNVtx.Z() - TrueVtx.Z()), cvWeight);
          if (m_TargetUtils->InLead4VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNZVertexResidualTgt4Lead->Fill((ANNVtx.Z() - TrueVtx.Z()), cvWeight);
          if (m_TargetUtils->InIron5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZVertexResidualTgt5Iron->Fill((ANNVtx.Z() - TrueVtx.Z()), cvWeight);
          if (m_TargetUtils->InLead5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZVertexResidualTgt5Lead->Fill((ANNVtx.Z() - TrueVtx.Z()), cvWeight);
          if (m_TargetUtils->InWaterTargetVolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNZVertexResidualWater->Fill((ANNVtx.Z() - TrueVtx.Z()), cvWeight);
          if (m_TargetUtils->InTracker(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNZVertexResidualTracker->Fill((ANNVtx.Z() - TrueVtx.Z()), cvWeight);

          ANNZVertexResidualVsConf->Fill((ANNVtx.Z() - TrueVtx.Z()), ANNProb, cvWeight);
          if (m_TargetUtils->InIron2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZVertexResidualVsConfTgt2Iron->Fill((ANNVtx.Z() - TrueVtx.Z()), ANNProb, cvWeight);
          if (m_TargetUtils->InLead2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZVertexResidualVsConfTgt2Lead->Fill((ANNVtx.Z() - TrueVtx.Z()), ANNProb, cvWeight);
          if (m_TargetUtils->InIron3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZVertexResidualVsConfTgt3Iron->Fill((ANNVtx.Z() - TrueVtx.Z()), ANNProb, cvWeight);
          if (m_TargetUtils->InLead3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZVertexResidualVsConfTgt3Lead->Fill((ANNVtx.Z() - TrueVtx.Z()), ANNProb, cvWeight);
          if (m_TargetUtils->InCarbon3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZVertexResidualVsConfTgt3Carbon->Fill((ANNVtx.Z() - TrueVtx.Z()), ANNProb, cvWeight);
          if (m_TargetUtils->InLead4VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNZVertexResidualVsConfTgt4Lead->Fill((ANNVtx.Z() - TrueVtx.Z()), ANNProb, cvWeight);
          if (m_TargetUtils->InIron5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZVertexResidualVsConfTgt5Iron->Fill((ANNVtx.Z() - TrueVtx.Z()), ANNProb, cvWeight);
          if (m_TargetUtils->InLead5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZVertexResidualVsConfTgt5Lead->Fill((ANNVtx.Z() - TrueVtx.Z()), ANNProb, cvWeight);
          if (m_TargetUtils->InWaterTargetVolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNZVertexResidualVsConfWater->Fill((ANNVtx.Z() - TrueVtx.Z()), ANNProb, cvWeight);
          if (m_TargetUtils->InTracker(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNZVertexResidualVsConfTracker->Fill((ANNVtx.Z() - TrueVtx.Z()), ANNProb, cvWeight);

          double WeightedANNZ = universe->GetANNSegmentsZPosWeighted2(1000000); //Just some large number
          ANNWeightedZVertexResidualVsConf->Fill((ANNVtx.Z() - TrueVtx.Z()), ANNProb, cvWeight);
          if (m_TargetUtils->InIron2VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZVertexResidualVsConfTgt2Iron->Fill((ANNVtx.Z() - TrueVtx.Z()), ANNProb, cvWeight);
          if (m_TargetUtils->InLead2VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZVertexResidualVsConfTgt2Lead->Fill((ANNVtx.Z() - TrueVtx.Z()), ANNProb, cvWeight);
          if (m_TargetUtils->InIron3VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZVertexResidualVsConfTgt3Iron->Fill((ANNVtx.Z() - TrueVtx.Z()), ANNProb, cvWeight);
          if (m_TargetUtils->InLead3VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZVertexResidualVsConfTgt3Lead->Fill((ANNVtx.Z() - TrueVtx.Z()), ANNProb, cvWeight);
          if (m_TargetUtils->InCarbon3VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZVertexResidualVsConfTgt3Carbon->Fill((ANNVtx.Z() - TrueVtx.Z()), ANNProb, cvWeight);
          if (m_TargetUtils->InLead4VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ))  ANNWeightedZVertexResidualVsConfTgt4Lead->Fill((ANNVtx.Z() - TrueVtx.Z()), ANNProb, cvWeight);
          if (m_TargetUtils->InIron5VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZVertexResidualVsConfTgt5Iron->Fill((ANNVtx.Z() - TrueVtx.Z()), ANNProb, cvWeight);
          if (m_TargetUtils->InLead5VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZVertexResidualVsConfTgt5Lead->Fill((ANNVtx.Z() - TrueVtx.Z()), ANNProb, cvWeight);
          if (m_TargetUtils->InWaterTargetVolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ))  ANNWeightedZVertexResidualVsConfWater->Fill((ANNVtx.Z() - TrueVtx.Z()), ANNProb, cvWeight);
          if (m_TargetUtils->InTracker(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ))  ANNWeightedZVertexResidualVsConfTracker->Fill((ANNVtx.Z() - TrueVtx.Z()), ANNProb, cvWeight);

          ANNZWeightedVerticesMC_VsEhad->Fill(WeightedANNZ, erecoil, cvWeight); //Some arbitrarily large cutoff
          //for (double cutoff = 0; cutoff <=0.5; cutoff+=0.05)
          for (int cutoff = 0; cutoff <100; cutoff++)
          {
            double WeightedANNZPos = universe->GetANNSegmentsZPosWeighted2(cutoff);
            ANNZWeightedVerticesMC_VsCutoff->Fill(WeightedANNZPos, cutoff, cvWeight);

            ANNWeightedVsUnweightedZVertexDifferenceMC->Fill((WeightedANNZPos - ANNVtx.Z()), cutoff, cvWeight);
            if (m_TargetUtils->InIron2VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZPos, 850, true))  ANNWeightedVsUnweightedZVertexDifferenceTgt2IronMC->Fill((WeightedANNZPos - ANNVtx.Z()), cutoff, cvWeight);
            if (m_TargetUtils->InLead2VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZPos, 850, true))  ANNWeightedVsUnweightedZVertexDifferenceTgt2LeadMC->Fill((WeightedANNZPos - ANNVtx.Z()), cutoff, cvWeight);
            if (m_TargetUtils->InIron3VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZPos, 850, true))  ANNWeightedVsUnweightedZVertexDifferenceTgt3IronMC->Fill((WeightedANNZPos - ANNVtx.Z()), cutoff, cvWeight);
            if (m_TargetUtils->InLead3VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZPos, 850, true))  ANNWeightedVsUnweightedZVertexDifferenceTgt3LeadMC->Fill((WeightedANNZPos - ANNVtx.Z()), cutoff, cvWeight);
            if (m_TargetUtils->InCarbon3VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZPos, 850, true))  ANNWeightedVsUnweightedZVertexDifferenceTgt3CarbonMC->Fill((WeightedANNZPos - ANNVtx.Z()), cutoff, cvWeight);
            if (m_TargetUtils->InLead4VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZPos))  ANNWeightedVsUnweightedZVertexDifferenceTgt4LeadMC->Fill((WeightedANNZPos - ANNVtx.Z()), cutoff, cvWeight);
            if (m_TargetUtils->InIron5VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZPos, 850, true))  ANNWeightedVsUnweightedZVertexDifferenceTgt5IronMC->Fill((WeightedANNZPos - ANNVtx.Z()), cutoff, cvWeight);
            if (m_TargetUtils->InLead5VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZPos, 850, true))  ANNWeightedVsUnweightedZVertexDifferenceTgt5LeadMC->Fill((WeightedANNZPos - ANNVtx.Z()), cutoff, cvWeight);
            if (m_TargetUtils->InWaterTargetVolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZPos))  ANNWeightedVsUnweightedZVertexDifferenceWaterMC->Fill((WeightedANNZPos - ANNVtx.Z()), cutoff, cvWeight);
            if (m_TargetUtils->InTracker(ANNVtx.X(), ANNVtx.Y(), WeightedANNZPos))  ANNWeightedVsUnweightedZVertexDifferenceTrackerMC->Fill((WeightedANNZPos - ANNVtx.Z()), cutoff, cvWeight);

            ANNWeightedZVertexResidual->Fill((WeightedANNZPos - TrueVtx.Z()), cutoff, cvWeight);
            if (m_TargetUtils->InIron2VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZPos, 850, true))  ANNWeightedZVertexResidualTgt2Iron->Fill((WeightedANNZPos - TrueVtx.Z()), cutoff, cvWeight);
            if (m_TargetUtils->InLead2VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZPos, 850, true))  ANNWeightedZVertexResidualTgt2Lead->Fill((WeightedANNZPos - TrueVtx.Z()), cutoff, cvWeight);
            if (m_TargetUtils->InIron3VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZPos, 850, true))  ANNWeightedZVertexResidualTgt3Iron->Fill((WeightedANNZPos - TrueVtx.Z()), cutoff, cvWeight);
            if (m_TargetUtils->InLead3VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZPos, 850, true))  ANNWeightedZVertexResidualTgt3Lead->Fill((WeightedANNZPos - TrueVtx.Z()), cutoff, cvWeight);
            if (m_TargetUtils->InCarbon3VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZPos, 850, true))  ANNWeightedZVertexResidualTgt3Carbon->Fill((WeightedANNZPos - TrueVtx.Z()), cutoff, cvWeight);
            if (m_TargetUtils->InLead4VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZPos))  ANNWeightedZVertexResidualTgt4Lead->Fill((WeightedANNZPos - TrueVtx.Z()), cutoff, cvWeight);
            if (m_TargetUtils->InIron5VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZPos, 850, true))  ANNWeightedZVertexResidualTgt5Iron->Fill((WeightedANNZPos - TrueVtx.Z()), cutoff, cvWeight);
            if (m_TargetUtils->InLead5VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZPos, 850, true))  ANNWeightedZVertexResidualTgt5Lead->Fill((WeightedANNZPos - TrueVtx.Z()), cutoff, cvWeight);
            if (m_TargetUtils->InWaterTargetVolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZPos))  ANNWeightedZVertexResidualWater->Fill((WeightedANNZPos - TrueVtx.Z()), cutoff, cvWeight);
            if (m_TargetUtils->InTracker(ANNVtx.X(), ANNVtx.Y(), WeightedANNZPos))  ANNWeightedZVertexResidualTracker->Fill((WeightedANNZPos - TrueVtx.Z()), cutoff, cvWeight);

          }

          ANNZResidualVsConfDifference->Fill((ANNVtx.Z() - TrueVtx.Z()), (universe->GetVecElem("ANN_plane_probs", 0) - universe->GetVecElem("ANN_plane_probs", 1)), cvWeight);
          if (m_TargetUtils->InIron2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZResidualVsConfDifferenceTgt2Iron->Fill((ANNVtx.Z() - TrueVtx.Z()), (universe->GetVecElem("ANN_plane_probs", 0) - universe->GetVecElem("ANN_plane_probs", 1)), cvWeight);
          if (m_TargetUtils->InLead2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZResidualVsConfDifferenceTgt2Lead->Fill((ANNVtx.Z() - TrueVtx.Z()), (universe->GetVecElem("ANN_plane_probs", 0) - universe->GetVecElem("ANN_plane_probs", 1)), cvWeight);
          if (m_TargetUtils->InIron3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZResidualVsConfDifferenceTgt3Iron->Fill((ANNVtx.Z() - TrueVtx.Z()), (universe->GetVecElem("ANN_plane_probs", 0) - universe->GetVecElem("ANN_plane_probs", 1)), cvWeight);
          if (m_TargetUtils->InLead3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZResidualVsConfDifferenceTgt3Lead->Fill((ANNVtx.Z() - TrueVtx.Z()), (universe->GetVecElem("ANN_plane_probs", 0) - universe->GetVecElem("ANN_plane_probs", 1)), cvWeight);
          if (m_TargetUtils->InCarbon3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZResidualVsConfDifferenceTgt3Carbon->Fill((ANNVtx.Z() - TrueVtx.Z()), (universe->GetVecElem("ANN_plane_probs", 0) - universe->GetVecElem("ANN_plane_probs", 1)), cvWeight);
          if (m_TargetUtils->InLead4VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNZResidualVsConfDifferenceTgt4Lead->Fill((ANNVtx.Z() - TrueVtx.Z()), (universe->GetVecElem("ANN_plane_probs", 0) - universe->GetVecElem("ANN_plane_probs", 1)), cvWeight);
          if (m_TargetUtils->InIron5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZResidualVsConfDifferenceTgt5Iron->Fill((ANNVtx.Z() - TrueVtx.Z()), (universe->GetVecElem("ANN_plane_probs", 0) - universe->GetVecElem("ANN_plane_probs", 1)), cvWeight);
          if (m_TargetUtils->InLead5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZResidualVsConfDifferenceTgt5Lead->Fill((ANNVtx.Z() - TrueVtx.Z()), (universe->GetVecElem("ANN_plane_probs", 0) - universe->GetVecElem("ANN_plane_probs", 1)), cvWeight);
          if (m_TargetUtils->InWaterTargetVolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNZResidualVsConfDifferenceWater->Fill((ANNVtx.Z() - TrueVtx.Z()), (universe->GetVecElem("ANN_plane_probs", 0) - universe->GetVecElem("ANN_plane_probs", 1)), cvWeight);
          if (m_TargetUtils->InTracker(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNZResidualVsConfDifferenceTracker->Fill((ANNVtx.Z() - TrueVtx.Z()), (universe->GetVecElem("ANN_plane_probs", 0) - universe->GetVecElem("ANN_plane_probs", 1)), cvWeight);

        }
        TBVerticesMC_ByModule->Fill(TrackBasedVtx.Z(), cvWeight);
        TBVerticesMC_ByZPos->Fill( TrackBasedVtx.Z(), cvWeight);
        TruthVerticesMC_ByModule->Fill( TrueVtx.Z(), cvWeight);
        TruthVerticesMC_ByZPos->Fill( TrueVtx.Z(), cvWeight);

        TruthVerticesMC_ByMod->Fill(universe->GetTruthVtxModule(), cvWeight);

        TBVerticesMCERecoil_ByModule->Fill( TrackBasedVtx.Z(), erecoil, cvWeight);
        TBVerticesMCMultiplicity_ByModule->Fill( TrackBasedVtx.Z(), multiplicity, cvWeight);
        TruthVerticesMCERecoil_ByModule->Fill( TrueVtx.Z(), erecoil, cvWeight);
        TruthVerticesMCMultiplicity_ByModule->Fill( TrueVtx.Z(), multiplicity, cvWeight);

        ErecoilMC->Fill( erecoil, cvWeight);

        TruthVerticesMC_BySegment->Fill(universe->GetTruthSegment(), cvWeight);
        TruthVerticesMC_BySegmentVsERecoil->Fill(universe->GetTruthSegment(), erecoil, cvWeight);


        const bool isSignal = michelcuts.isSignal(*universe, weight);
        if (isSignal) // If it is signal
        {
          if(ANNVtx.X()!=-1.0)
          {
            ANNPlaneProbabilityVsEhadNumerator->Fill(ANNProb, erecoil, cvWeight);
            if (m_TargetUtils->InTracker(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()) && m_TargetUtils->InTracker(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z())) ANNPlaneProbabilityTrackerVsEhadNumerator->Fill(ANNProb, erecoil, cvWeight);
            if (m_TargetUtils->InIron2VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true) && m_TargetUtils->InIron2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true)) ANNPlaneProbabilityTgt2IronVsEhadNumerator->Fill(ANNProb, erecoil, cvWeight);
            if (m_TargetUtils->InLead2VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true) && m_TargetUtils->InLead2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true)) ANNPlaneProbabilityTgt2LeadVsEhadNumerator->Fill(ANNProb, erecoil, cvWeight);
            if (m_TargetUtils->InIron3VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true) && m_TargetUtils->InIron3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true)) ANNPlaneProbabilityTgt3IronVsEhadNumerator->Fill(ANNProb, erecoil, cvWeight);
            if (m_TargetUtils->InLead3VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true) && m_TargetUtils->InLead3VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true)) ANNPlaneProbabilityTgt3LeadVsEhadNumerator->Fill(ANNProb, erecoil, cvWeight);
            if (m_TargetUtils->InCarbon3VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true) && m_TargetUtils->InCarbon3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true)) ANNPlaneProbabilityTgt3CarbonVsEhadNumerator->Fill(ANNProb, erecoil, cvWeight);
            if (m_TargetUtils->InLead4VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()) && m_TargetUtils->InLead4VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z())) ANNPlaneProbabilityTgt4LeadVsEhadNumerator->Fill(ANNProb, erecoil, cvWeight);
            if (m_TargetUtils->InIron5VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true) && m_TargetUtils->InIron5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true)) ANNPlaneProbabilityTgt5IronVsEhadNumerator->Fill(ANNProb, erecoil, cvWeight);
            if (m_TargetUtils->InLead5VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true) && m_TargetUtils->InLead5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true)) ANNPlaneProbabilityTgt5LeadVsEhadNumerator->Fill(ANNProb, erecoil, cvWeight);
            if (m_TargetUtils->InWaterTargetVolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()) && m_TargetUtils->InWaterTargetVolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z())) ANNPlaneProbabilityWaterVsEhadNumerator->Fill(ANNProb, erecoil, cvWeight);
          }
        }
      }
    }
  } //End entries loop
  std::cout << "Finished MC reco loop.\n"; 
}

void LoopAndFillData( PlotUtils::ChainWrapper* data,
			        std::vector<CVUniverse*> data_band,
				PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts
        )

{
  std::cout << "Starting data loop...\n";
  std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
  const int nEntries = data->GetEntries();
  //const int nEntries = 10000;
  int  startNum = 0;
  int endNum = nEntries;
  if (nSubruns != 0 )
  {
    int chunkSize = std::floor(nEntries/nSubruns);
    startNum = nProcess*chunkSize;
    if (nProcess == (nSubruns-1)) //The last one, run to the end
    {
      endNum = nEntries; //probably redundant with the same line just above
    }
    else
    {
      endNum = (nProcess+1)*chunkSize;
    }
  }
  // const int nEntries = 10000;
  for (int i = startNum; i < endNum; ++i)
  {
    for (auto universe : data_band) {
      universe->SetEntry(i);
      if(i%1000==0)
      { 
        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
        double mselapsed = std::chrono::duration_cast<std::chrono::milliseconds> (now - startTime).count();
        double msperevent = mselapsed/i;
        int sleft = (nEntries-i)*msperevent/1000;
        double selapsed = ((int)mselapsed)/1000;
        std::cout << i << " / " << nEntries << ". Time Elapsed: " << selapsed <<" Est Remaining: " << sleft <<"\r" <<std::flush;
      }
      MichelEvent myevent; 
      ROOT::Math::XYZVector ANNVtx = universe->GetANNVertex();
      ROOT::Math::XYZTVector TrackBasedVtx = universe->GetVertex();
      double batchPOT = universe->GetBatchPOT();
      //Incorporate batch POT efficiency scaling as applied above
      double erecoil = universe->GetRecoilE()/pow(10,3);
      double ANNProb = universe->GetANNProb();
      double curvsig = 1/universe->GetMuonQPErr();
      int multiplicity = universe->GetMultiplicity();
      double pmu = universe->GetPmu()/1000;
      PlotUtils::TargetUtils* m_TargetUtils=new PlotUtils::TargetUtils();

      if (!michelcuts.isDataSelected(*universe, myevent).all()) continue;   //All cuts except ANN confidence cut
  
      if(ANNVtx.X()!=-1.0)
      {
        ANNPlaneProbabilityVsPmuData->Fill(ANNProb, pmu);
        ANNPlaneProbabilityVsEhadData->Fill(ANNProb, erecoil);
        if (m_TargetUtils->InIron2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNPlaneProbabilityTgt2IronVsEhadData->Fill(ANNProb, erecoil);
        if (m_TargetUtils->InLead2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNPlaneProbabilityTgt2LeadVsEhadData->Fill(ANNProb, erecoil);
        if (m_TargetUtils->InIron3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNPlaneProbabilityTgt3IronVsEhadData->Fill(ANNProb, erecoil);
        if (m_TargetUtils->InLead3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNPlaneProbabilityTgt3LeadVsEhadData->Fill(ANNProb, erecoil);
        if (m_TargetUtils->InCarbon3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNPlaneProbabilityTgt3CarbonVsEhadData->Fill(ANNProb, erecoil);
        if (m_TargetUtils->InLead4VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNPlaneProbabilityTgt4LeadVsEhadData->Fill(ANNProb, erecoil);
        if (m_TargetUtils->InIron5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNPlaneProbabilityTgt5IronVsEhadData->Fill(ANNProb, erecoil);
        if (m_TargetUtils->InLead5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNPlaneProbabilityTgt5LeadVsEhadData->Fill(ANNProb, erecoil);
        if (m_TargetUtils->InWaterTargetVolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNPlaneProbabilityWaterVsEhadData->Fill(ANNProb, erecoil);
        if (m_TargetUtils->InTracker(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNPlaneProbabilityTrackerVsEhadData->Fill(ANNProb, erecoil);

        ANNVerticesDataANNConf_ByModule->Fill( ANNVtx.Z(), ANNProb);
        ANNVerticesDataCurvSig_ByModule->Fill( ANNVtx.Z(), curvsig);

        ANNVerticesData_ByZPosVsANNConf->Fill( ANNVtx.Z(), ANNProb);
      }

      if (!ANNConfCut->passesCut(*universe, myevent)) continue;   //ANN confidence cut
      
      if(ANNVtx.X()!=-1.0)
      {
        ANNVerticesData_ByModule->Fill(ANNVtx.Z());
        ANNVerticesData_ByZPos->Fill(ANNVtx.Z());
        //for (double cutoff = 0; cutoff <=0.5; cutoff+=0.05)
        ANNZWeightedVerticesData_VsEhad->Fill(universe->GetANNSegmentsZPosWeighted2(1000000), erecoil); //Some arbitrarily large cutoff

        for (int cutoff = 0; cutoff <100; cutoff++)
        {
          double WeightedANNZPos = universe->GetANNSegmentsZPosWeighted2(cutoff);
          ANNZWeightedVerticesData_VsCutoff->Fill(WeightedANNZPos, cutoff);
          ANNWeightedVsUnweightedZVertexDifferenceData->Fill((WeightedANNZPos - ANNVtx.Z()), cutoff);
          if (m_TargetUtils->InIron2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNWeightedVsUnweightedZVertexDifferenceTgt2IronData->Fill((WeightedANNZPos - ANNVtx.Z()), cutoff);
          if (m_TargetUtils->InLead2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNWeightedVsUnweightedZVertexDifferenceTgt2LeadData->Fill((WeightedANNZPos - ANNVtx.Z()), cutoff);
          if (m_TargetUtils->InIron3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNWeightedVsUnweightedZVertexDifferenceTgt3IronData->Fill((WeightedANNZPos - ANNVtx.Z()), cutoff);
          if (m_TargetUtils->InLead3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNWeightedVsUnweightedZVertexDifferenceTgt3LeadData->Fill((WeightedANNZPos - ANNVtx.Z()), cutoff);
          if (m_TargetUtils->InCarbon3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNWeightedVsUnweightedZVertexDifferenceTgt3CarbonData->Fill((WeightedANNZPos - ANNVtx.Z()), cutoff);
          if (m_TargetUtils->InLead4VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNWeightedVsUnweightedZVertexDifferenceTgt4LeadData->Fill((WeightedANNZPos - ANNVtx.Z()), cutoff);
          if (m_TargetUtils->InIron5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNWeightedVsUnweightedZVertexDifferenceTgt5IronData->Fill((WeightedANNZPos - ANNVtx.Z()), cutoff);
          if (m_TargetUtils->InLead5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNWeightedVsUnweightedZVertexDifferenceTgt5LeadData->Fill((WeightedANNZPos - ANNVtx.Z()), cutoff);
          if (m_TargetUtils->InWaterTargetVolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNWeightedVsUnweightedZVertexDifferenceWaterData->Fill((WeightedANNZPos - ANNVtx.Z()), cutoff);
          if (m_TargetUtils->InTracker(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNWeightedVsUnweightedZVertexDifferenceTrackerData->Fill((WeightedANNZPos - ANNVtx.Z()), cutoff);

        }

        ANNVerticesData_ByZPosVsERecoil->Fill(ANNVtx.Z(), erecoil);

        ANNVerticesData_ByMod->Fill(universe->GetANNVtxModule());

        ANNVerticesDataERecoil_ByModule->Fill( ANNVtx.Z(), erecoil);
        ANNVerticesDataMultiplicity_ByModule->Fill( ANNVtx.Z(), multiplicity);

        ANNVerticesDataCurvSig_ByModule->Fill( ANNVtx.Z(), curvsig);

        ANNVerticesData_BySegment->Fill(universe->GetANNSegment());
        ANNVerticesData_BySegmentVsERecoil->Fill(universe->GetANNSegment(), erecoil);
      }
      TBVerticesData_ByZPos->Fill(TrackBasedVtx.Z());
      TBVerticesData_ByModule->Fill(TrackBasedVtx.Z());
      TBVerticesData_ByMod->Fill(universe->GetMADVtxModule());

      TBVerticesDataERecoil_ByModule->Fill( TrackBasedVtx.Z(), erecoil);
      TBVerticesDataMultiplicity_ByModule->Fill( TrackBasedVtx.Z(), multiplicity);
      TBVerticesDataCurvSig_ByModule->Fill( TrackBasedVtx.Z(), curvsig);
      ErecoilData->Fill( erecoil);
    }
  }
  std::cout << "Finished data loop.\n";
}


void LoopAndFillEffDenom(PlotUtils::ChainWrapper *truth,
                         std::map<std::string, std::vector<CVUniverse *>> truth_bands,
                         PlotUtils::Cutter<CVUniverse, MichelEvent> &michelcuts,
                         PlotUtils::Model<CVUniverse, MichelEvent> &model)
{
  assert(!truth_bands["cv"].empty() && "\"cv\" error band is empty!  Could not set Model entry.");
  auto &cvUniv = truth_bands["cv"].front();

  std::cout << "Starting efficiency denominator loop...\n";
  const int nEntries = truth->GetEntries();
  //const int nEntries = 10000;
  for (int i = 0; i < nEntries; ++i)
  {
    if (i % 1000 == 0)
      std::cout << i << " / " << nEntries << "\r" << std::flush;

    MichelEvent cvEvent;
    cvUniv->SetEntry(i);
    model.SetEntry(*cvUniv, cvEvent);
    const double cvWeight = model.GetWeight(*cvUniv, cvEvent);

    //=========================================
    // Systematics loop(s)
    //=========================================
    for (auto band : truth_bands)
    {
      std::vector<CVUniverse *> truth_band_universes = band.second;
      for (auto universe : truth_band_universes)
      {

        MichelEvent myevent; // Only used to keep the Model happy

        // Tell the Event which entry in the TChain it's looking at
        universe->SetEntry(i);
        if (!michelcuts.isEfficiencyDenom(*universe, cvWeight))
          continue;                                                // Weight is ignored for isEfficiencyDenom() in all but the CV universe
        const double weight = model.GetWeight(*universe, myevent); // Only calculate the weight for events that will use it
        ROOT::Math::XYZTVector TrueVtx = universe->GetTrueVertex();
        double ANNProb = universe->GetANNProb();
        double erecoil = universe->GetRecoilE()/pow(10,3);
        PlotUtils::TargetUtils* m_TargetUtils=new PlotUtils::TargetUtils();
    
        ANNPlaneProbabilityVsEhadDenominator->Fill(ANNProb, erecoil, cvWeight);
        if (m_TargetUtils->InTracker(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z())) ANNPlaneProbabilityTrackerVsEhadDenominator->Fill(ANNProb, erecoil, weight);
        if (m_TargetUtils->InIron2VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true)) ANNPlaneProbabilityTgt2IronVsEhadDenominator->Fill(ANNProb, erecoil, weight);
        if (m_TargetUtils->InLead2VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true)) ANNPlaneProbabilityTgt2LeadVsEhadDenominator->Fill(ANNProb, erecoil, weight);
        if (m_TargetUtils->InIron3VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true)) ANNPlaneProbabilityTgt3IronVsEhadDenominator->Fill(ANNProb, erecoil, weight);
        if (m_TargetUtils->InLead3VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true)) ANNPlaneProbabilityTgt3LeadVsEhadDenominator->Fill(ANNProb, erecoil, weight);
        if (m_TargetUtils->InCarbon3VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true)) ANNPlaneProbabilityTgt3CarbonVsEhadDenominator->Fill(ANNProb, erecoil, weight);
        if (m_TargetUtils->InLead4VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z())) ANNPlaneProbabilityTgt4LeadVsEhadDenominator->Fill(ANNProb, erecoil, weight);
        if (m_TargetUtils->InIron5VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true)) ANNPlaneProbabilityTgt5IronVsEhadDenominator->Fill(ANNProb, erecoil, weight);
        if (m_TargetUtils->InLead5VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true)) ANNPlaneProbabilityTgt5LeadVsEhadDenominator->Fill(ANNProb, erecoil, weight);
        if (m_TargetUtils->InWaterTargetVolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z())) ANNPlaneProbabilityWaterVsEhadDenominator->Fill(ANNProb, erecoil, weight);



      }
    }
  }
  std::cout << "Finished efficiency denominator loop.\n";
}


//Returns false if recoTreeName could not be inferred
bool inferRecoTreeNameAndCheckTreeNames(const std::string& mcPlaylistName, const std::string& dataPlaylistName, std::string& recoTreeName)
{
  const std::vector<std::string> knownTreeNames = {"Truth", "Meta"};
  bool areFilesOK = false;

  std::ifstream playlist(mcPlaylistName);
  std::string firstFile = "";
  playlist >> firstFile;
  auto testFile = TFile::Open(firstFile.c_str());
  if(!testFile)
  {
    std::cerr << "Failed to open the first MC file at " << firstFile << "\n";
    return false;
  }

  //Does the MC playlist have the Truth tree?  This is needed for the efficiency denominator.
  const auto truthTree = testFile->Get("Truth");
  if(truthTree == nullptr || !truthTree->IsA()->InheritsFrom(TClass::GetClass("TTree")))
  {
    std::cerr << "Could not find the \"Truth\" tree in MC file named " << firstFile << "\n";
    return false;
  }

  //Figure out what the reco tree name is
  for(auto key: *testFile->GetListOfKeys())
  {
    if(static_cast<TKey*>(key)->ReadObj()->IsA()->InheritsFrom(TClass::GetClass("TTree"))
       && std::find(knownTreeNames.begin(), knownTreeNames.end(), key->GetName()) == knownTreeNames.end())
    {
      recoTreeName = key->GetName();
      areFilesOK = true;
    }
  }
  delete testFile;
  testFile = nullptr;

  //Make sure the data playlist's first file has the same reco tree
  playlist.open(dataPlaylistName);
  playlist >> firstFile;
  testFile = TFile::Open(firstFile.c_str());
  if(!testFile)
  {
    std::cerr << "Failed to open the first data file at " << firstFile << "\n";
    return false;
  }

  const auto recoTree = testFile->Get(recoTreeName.c_str());
  if(recoTree == nullptr || !recoTree->IsA()->InheritsFrom(TClass::GetClass("TTree")))
  {
    std::cerr << "Could not find the \"" << recoTreeName << "\" tree in data file named " << firstFile << "\n";
    return false;
  }

  return areFilesOK;
}

std::string breakUpInputFileList(std::string inputFile, int numTotal, int thisNum) {
    std::ifstream in(inputFile);
    if (!in.is_open()) {
        std::cerr << "Error opening input file: " << inputFile << std::endl;
        return "";
    }

    // Read all lines into memory
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        lines.push_back(line);
    }
    in.close();

    size_t totalLines = lines.size();
    if (totalLines == 0) {
        std::cerr << "Input file is empty!" << std::endl;
        return "";
    }

    // Split into 10 chunks
    size_t chunkSize = std::floor(totalLines / numTotal);
    std::string outfile = inputFile+"tmp"+std::to_string(thisNum)+".txt";
    std::ofstream out(outfile);
    if (!out.is_open()) {
        std::cerr << "Error opening output file: " << outfile << std::endl;
        return "";
    }
    int start = thisNum * chunkSize;
    int end = (thisNum +1)* chunkSize;
    if (thisNum == numTotal-1 ) end = lines.size();
    for (size_t j = start; j < end; ++j) {
        out << lines[j] << "\n";
    }
    out.close();

    return outfile;
}

//==============================================================================
// Main
//==============================================================================
int main(const int argc, char** argv)
{
  TH1::AddDirectory(false);

  //Validate input.
  //I expect a data playlist file name and an MC playlist file name which is exactly 2 arguments.
  const int nArgsExpected = 3;
  if(argc > nArgsExpected + 1) //argc is the size of argv.  I check for number of arguments + 1 because
                                //argv[0] is always the path to the executable.
  {
    std::cerr << "Expected up to" << nArgsExpected << " arguments, but got " << argc - 1 << "\n" << USAGE << "\n";
    return badCmdLine;
  }
  if (argc==4 and std::string(argv[3])=="-v") verbose = true; 

  //One playlist must contain only MC files, and the other must contain only data files.
  //Only checking the first file in each playlist because opening each file an extra time
  //remotely (e.g. through xrootd) can get expensive.
  //TODO: Look in INSTALL_DIR if files not found?
  std::string  mc_file_list = argv[2],
                    data_file_list = argv[1];

  //Check that necessary TTrees exist in the first file of mc_file_list and data_file_list
  std::string reco_tree_name;
  if(!inferRecoTreeNameAndCheckTreeNames(mc_file_list, data_file_list, reco_tree_name))
  {
    std::cerr << "Failed to find required trees in MC playlist " << mc_file_list << " and/or data playlist " << data_file_list << ".\n" << USAGE << "\n";
    return badInputFile;
  }

  const bool doCCQENuValidation = (reco_tree_name == "CCQENu"); //Enables extra histograms and might influence which systematics I use.


  char* numSubruns = getenv("NumGridSubruns");
  char* numProcess = getenv("PROCESS");

  if (numSubruns != nullptr && numProcess != nullptr)
  {
    nSubruns = std::stoi(numSubruns);
    nProcess = std::stoi(numProcess);
    //mc_file_list= breakUpInputFileList(mc_file_list, nSubruns , nProcess);
    //data_file_list= breakUpInputFileList(data_file_list, nSubruns , nProcess);
  }

  //const bool is_grid = false; //TODO: Are we going to put this back?  Gonzalo needs it iirc.
  std::cout<<"Here\n";
  PlotUtils::MacroUtil options(reco_tree_name, mc_file_list, data_file_list, "minervame1A", true); //minervame1A is just a placeholder, it gets overwritted immediately below
  std::cout<<"Here1\n";
  options.m_plist_string = util::GetPlaylist(*options.m_mc, true); //TODO: Put GetPlaylist into PlotUtils::MacroUtil
  std::cout<<"Here2\n";

  // You're required to make some decisions
  PlotUtils::MinervaUniverse::SetNuEConstraint(true);
  PlotUtils::MinervaUniverse::SetPlaylist(options.m_plist_string); //TODO: Infer this from the files somehow?
  int nuoranu = nuOrAntiNuMode(options.m_plist_string);
  int nupdg;
  if (nuoranu==1) nupdg = 14;
  else if (nuoranu==2) nupdg = -14;
  PlotUtils::MinervaUniverse::SetAnalysisNuPDG(nupdg);
  //PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
  PlotUtils::MinervaUniverse::SetZExpansionFaReweight(false);

  PlotUtils::MinervaUniverse::RPAMaterials(true); 

  //Now that we've defined what a cross section is, decide which sample and model
  //we're extracting a cross section for.
  PlotUtils::Cutter<CVUniverse, MichelEvent>::reco_t sidebands, preCuts;
  PlotUtils::Cutter<CVUniverse, MichelEvent>::truth_t signalDefinition, phaseSpace;

  //preCuts = util::GetAnalysisCuts(nupdg);
  preCuts.emplace_back(new reco::Apothem<CVUniverse, MichelEvent>(util::apothem));
  preCuts.emplace_back(new reco::MaxMuonAngle<CVUniverse, MichelEvent>(17.));
  preCuts.emplace_back(new reco::HasMINOSMatch<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::NoDeadtime<CVUniverse, MichelEvent>(1, "Deadtime"));
  if (nupdg>0)  preCuts.emplace_back(new reco::IsNeutrino<CVUniverse, MichelEvent>()); //Used minos curvature
  else if (nupdg<0)  preCuts.emplace_back(new reco::IsAntiNeutrino<CVUniverse, MichelEvent>()); //Used minos curvature
  preCuts.emplace_back(new reco::MuonCurveSignificance<CVUniverse, MichelEvent>(5));
  preCuts.emplace_back(new reco::MuonEnergyMin<CVUniverse, MichelEvent>(2000.0, "EMu Min"));
  preCuts.emplace_back(new reco::MuonEnergyMax<CVUniverse, MichelEvent>(20000.0, "EMu Max"));
  ///////
  preCuts.emplace_back(new reco::ZRangeANN<CVUniverse, MichelEvent>("Z pos", PlotUtils::TargetProp::NukeRegion::Face, PlotUtils::TargetProp::Tracker::Back));

  //preCuts.emplace_back(new reco::ANNConfidenceCut<CVUniverse, MichelEvent>(0.40)); //Reccomended at 0.4 for P6 ML vertexing and 0.2 for P4 vertexing
  //preCuts.emplace_back(new reco::RockMuonCut<CVUniverse, MichelEvent>()); //Reccomended for P6 ML vertexing
  //preCuts.emplace_back(new reco::VetoWall<CVUniverse, MichelEvent>()); //Reccomended for P6 ML vertexing

                                                                                                                                   
  signalDefinition.emplace_back(new truth::IsNeutrino<CVUniverse>());
  signalDefinition.emplace_back(new truth::IsCC<CVUniverse>());
                                                                                                                                                   
  phaseSpace = util::GetPhaseSpace();
  phaseSpace.emplace_back(new truth::ZRange<CVUniverse>("Z pos", PlotUtils::TargetProp::NukeRegion::Face, PlotUtils::TargetProp::Tracker::Back));


  PlotUtils::Cutter<CVUniverse, MichelEvent> mycuts(std::move(preCuts), std::move(sidebands) , std::move(signalDefinition),std::move(phaseSpace));

  //Tune version vA.B.C
  int tuneA = 1;
  int tuneB = 0;
  int tuneC = 0;
  const char* mnvTuneIn =  getenv("MnvTune");
  if (mnvTuneIn != nullptr)
  {
    std::string mnvTuneStr = std::string(mnvTuneIn);
    if (mnvTuneStr.size()!=3 || !std::isdigit(mnvTuneStr[0]) || !std::isdigit(mnvTuneStr[1]) || !std::isdigit(mnvTuneStr[2])) "Unrecognised tune, using default";
    else
    { 
      int tune = std::stoi(mnvTuneStr);
      tuneC=tune%10;
      tuneB = ((tune-tuneC)/10)%10;
      tuneA = (tune -(tuneC + tuneB*10))/100;
    }
  }
  std::cout<< "Using minerva tune v" << tuneA << "."<< tuneB << "."<< tuneC << "\n";

  std::vector<std::unique_ptr<PlotUtils::Reweighter<CVUniverse, MichelEvent>>> MnvTune;
  //Setting the A component of mnvtune vA.B.C
  if (tuneA == 1)
  {
    MnvTune.emplace_back(new PlotUtils::FluxAndCVReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::GENIEReweighter<CVUniverse, MichelEvent>(true, false));
    MnvTune.emplace_back(new PlotUtils::MINOSEfficiencyReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::RPAReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::LowRecoil2p2hReweighter<CVUniverse, MichelEvent>());
  }
  if (tuneA == 2)
  {
    MnvTune.emplace_back(new PlotUtils::FluxAndCVReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::GENIEReweighter<CVUniverse, MichelEvent>(true, false));
    MnvTune.emplace_back(new PlotUtils::MINOSEfficiencyReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::RPAReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::LowRecoil2p2hReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::LowQ2PiReweighter<CVUniverse, MichelEvent>("JOINT")); //Is JOINT the correct option for mnvtune2?
  }
  if (tuneA == 3)
  {
    MnvTune.emplace_back(new PlotUtils::FluxAndCVReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::GENIEReweighter<CVUniverse, MichelEvent>(true, false));
    MnvTune.emplace_back(new PlotUtils::MINOSEfficiencyReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::RPAReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::SuSAFromValencia2p2hReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::BodekRitchieReweighter<CVUniverse, MichelEvent>(2)); //Is 2 the right mode?
  }
  if (tuneA == 4)
  {
    PlotUtils::MinervaUniverse::SetReadoutVolume("Nuke");
    PlotUtils::MinervaUniverse::SetMHRWeightNeutronCVReweight( true );
    PlotUtils::MinervaUniverse::SetMHRWeightElastics( true );
    MnvTune.emplace_back(new PlotUtils::FluxAndCVReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::GENIEReweighter<CVUniverse, MichelEvent>(true, true));
    MnvTune.emplace_back(new PlotUtils::MINOSEfficiencyReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::RPAReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::LowRecoil2p2hReweighter<CVUniverse, MichelEvent>());
    //Other decisions to add for MnvTunev4.3.1
  }
  //Setting the B component of mnvtune vA.B.C
  if (tuneB == 3)
  {
    MnvTune.emplace_back(new PlotUtils::LowQ2PiReweighter<CVUniverse, MichelEvent>("MENU1PI"));
    MnvTune.emplace_back(new PlotUtils::DiffractiveReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::COHPionReweighter<CVUniverse, MichelEvent>());
  }
  //Setting the C component of mnvtune vA.B.C
  if (tuneC == 1)
  {
    MnvTune.emplace_back(new PlotUtils::FSIReweighter<CVUniverse, MichelEvent>(true, true));
  }

  std::cout<<"Tune components applied:\n";
  for (auto&& t : MnvTune) std::cout<< "\t"<< t->GetName() <<std::endl;
  //Do we need all this for v 4.3.1? I found it somewhere else but idk if I need it here
  //https://github.com/MinervaExpt/LowRecoilPions/blob/902f51bd72e1dff74d26e0df7158f27750947521/studies2DEventLoop.cpp
  //Could also wrap all v431 cuts in one reweighter like https://github.com/MinervaExpt/LowRecoilPions/blob/902f51bd72e1dff74d26e0df7158f27750947521/twoDEventLoopSide.cpp
  //MnvTunev4.emplace_back(new PlotUtils::GeantNeutronCVReweighter<CVUniverse, MichelEvent>());
  //MnvTunev4.emplace_back(new PlotUtils::TargetMassReweighter<CVUniverse, MichelEvent>());

  //What about this?
  //MnvTunev4.emplace_back(new PlotUtils::PionReweighter<CVUniverse,MichelEvent>());

  PlotUtils::Model<CVUniverse, MichelEvent> model(std::move(MnvTune));
  

  PlotUtils::MinervaUniverse::SetNFluxUniverses(2); //Necessary to get Flux integral later...  Doesn't work with just 1 flux universe though because _that_ triggers "spread errors".

  // Loop entries and fill
  //try
  {
    std::map< std::string, std::vector<CVUniverse*> > error_bands;
    std::map<std::string, std::vector<CVUniverse*> > band_flux = PlotUtils::GetFluxSystematicsMap<CVUniverse>(options.m_mc, CVUniverse::GetNFluxUniverses());
    //error_bands.insert(band_flux.begin(), band_flux.end()); //Necessary to get flux integral later...
    error_bands["cv"] = {new CVUniverse(options.m_mc)};
    std::map< std::string, std::vector<CVUniverse*> > truth_bands;
    truth_bands["cv"] = {new CVUniverse(options.m_truth)};
    CVUniverse* data_univers = new CVUniverse(options.m_data);
    std::vector<CVUniverse*> data_band = {data_univers};
    CVUniverse::SetTruth(false);
    LoopAndFillMC(options.m_mc, error_bands, mycuts, model);
    CVUniverse::SetTruth(true);
    //LoopAndFillEffDenom(options.m_truth, truth_bands, mycuts, model);
    options.PrintMacroConfiguration(argv[0]);
    std::cout << "MC cut summary:\n" << mycuts << "\n";
    mycuts.resetStats();
    CVUniverse::SetTruth(false);
    LoopAndFillData(options.m_data, data_band, mycuts);
    //std::cout << "Data cut summary:\n" << mycuts << "\n";

    std::string outFileName = "VertexValidations.root";
    if (nSubruns != 0 ) outFileName = "VertexValidations_n"+std::to_string(nProcess)+ ".root";
    TFile* OutDir = TFile::Open(outFileName.c_str(), "RECREATE");

    ANNVerticesMC_ByModule->Write();
    TBVerticesMC_ByModule->Write();
    ANNVerticesData_ByModule->Write();
    TBVerticesData_ByModule->Write();
    TruthVerticesMC_ByModule->Write();
    ANNVerticesMC_ByZPos->Write();
    TBVerticesMC_ByZPos->Write();
    ANNVerticesData_ByZPos->Write();
    TBVerticesData_ByZPos->Write();
    TruthVerticesMC_ByZPos->Write();
    ANNVerticesData_ByMod->Write();
    TBVerticesData_ByMod->Write();
    TruthVerticesMC_ByMod->Write();

    ANNZWeightedVerticesData_VsCutoff->Write();
    ANNZWeightedVerticesMC_VsCutoff->Write();

    ANNZWeightedVerticesData_VsEhad->Write();
    ANNZWeightedVerticesMC_VsEhad->Write();

    ANNVerticesMC_BySegment->Write();
    TruthVerticesMC_BySegment->Write();

    ANNVerticesData_BySegment->Write();

    ANNZVertexResidual->Write();
    ANNZVertexResidualTracker->Write();
    ANNZVertexResidualTgt2Iron->Write();
    ANNZVertexResidualTgt2Lead->Write();
    ANNZVertexResidualTgt3Iron->Write();
    ANNZVertexResidualTgt3Lead->Write();
    ANNZVertexResidualTgt3Carbon->Write();
    ANNZVertexResidualTgt4Lead->Write();
    ANNZVertexResidualTgt5Iron->Write();
    ANNZVertexResidualTgt5Lead->Write();
    ANNZVertexResidualWater->Write();

    ANNZVertexResidualVsConf->Write();
    ANNZVertexResidualVsConfTracker->Write();
    ANNZVertexResidualVsConfTgt2Iron->Write();
    ANNZVertexResidualVsConfTgt2Lead->Write();
    ANNZVertexResidualVsConfTgt3Iron->Write();
    ANNZVertexResidualVsConfTgt3Lead->Write();
    ANNZVertexResidualVsConfTgt3Carbon->Write();
    ANNZVertexResidualVsConfTgt4Lead->Write();
    ANNZVertexResidualVsConfTgt5Iron->Write();
    ANNZVertexResidualVsConfTgt5Lead->Write();
    ANNZVertexResidualVsConfWater->Write();

    ANNWeightedZVertexResidualVsConf->Write();
    ANNWeightedZVertexResidualVsConfTracker->Write();
    ANNWeightedZVertexResidualVsConfTgt2Iron->Write();
    ANNWeightedZVertexResidualVsConfTgt2Lead->Write();
    ANNWeightedZVertexResidualVsConfTgt3Iron->Write();
    ANNWeightedZVertexResidualVsConfTgt3Lead->Write();
    ANNWeightedZVertexResidualVsConfTgt3Carbon->Write();
    ANNWeightedZVertexResidualVsConfTgt4Lead->Write();
    ANNWeightedZVertexResidualVsConfTgt5Iron->Write();
    ANNWeightedZVertexResidualVsConfTgt5Lead->Write();
    ANNWeightedZVertexResidualVsConfWater->Write();

    ANNWeightedZVertexResidual->Write();
    ANNWeightedZVertexResidualTracker->Write();
    ANNWeightedZVertexResidualTgt2Iron->Write();
    ANNWeightedZVertexResidualTgt2Lead->Write();
    ANNWeightedZVertexResidualTgt3Iron->Write();
    ANNWeightedZVertexResidualTgt3Lead->Write();
    ANNWeightedZVertexResidualTgt3Carbon->Write();
    ANNWeightedZVertexResidualTgt4Lead->Write();
    ANNWeightedZVertexResidualTgt5Iron->Write();
    ANNWeightedZVertexResidualTgt5Lead->Write();
    ANNWeightedZVertexResidualWater->Write();

    ANNWeightedVsUnweightedZVertexDifferenceMC->Write();
    ANNWeightedVsUnweightedZVertexDifferenceTrackerMC->Write();
    ANNWeightedVsUnweightedZVertexDifferenceTgt2IronMC->Write();
    ANNWeightedVsUnweightedZVertexDifferenceTgt2LeadMC->Write();
    ANNWeightedVsUnweightedZVertexDifferenceTgt3IronMC->Write();
    ANNWeightedVsUnweightedZVertexDifferenceTgt3LeadMC->Write();
    ANNWeightedVsUnweightedZVertexDifferenceTgt3CarbonMC->Write();
    ANNWeightedVsUnweightedZVertexDifferenceTgt4LeadMC->Write();
    ANNWeightedVsUnweightedZVertexDifferenceTgt5IronMC->Write();
    ANNWeightedVsUnweightedZVertexDifferenceTgt5LeadMC->Write();
    ANNWeightedVsUnweightedZVertexDifferenceWaterMC->Write();

    ANNWeightedVsUnweightedZVertexDifferenceData->Write();
    ANNWeightedVsUnweightedZVertexDifferenceTrackerData->Write();
    ANNWeightedVsUnweightedZVertexDifferenceTgt2IronData->Write();
    ANNWeightedVsUnweightedZVertexDifferenceTgt2LeadData->Write();
    ANNWeightedVsUnweightedZVertexDifferenceTgt3IronData->Write();
    ANNWeightedVsUnweightedZVertexDifferenceTgt3LeadData->Write();
    ANNWeightedVsUnweightedZVertexDifferenceTgt3CarbonData->Write();
    ANNWeightedVsUnweightedZVertexDifferenceTgt4LeadData->Write();
    ANNWeightedVsUnweightedZVertexDifferenceTgt5IronData->Write();
    ANNWeightedVsUnweightedZVertexDifferenceTgt5LeadData->Write();
    ANNWeightedVsUnweightedZVertexDifferenceWaterData->Write();

    ANNZResidualVsConfDifference->Write();
    ANNZResidualVsConfDifferenceTracker->Write();
    ANNZResidualVsConfDifferenceTgt2Iron->Write();
    ANNZResidualVsConfDifferenceTgt2Lead->Write();
    ANNZResidualVsConfDifferenceTgt3Iron->Write();
    ANNZResidualVsConfDifferenceTgt3Lead->Write();
    ANNZResidualVsConfDifferenceTgt3Carbon->Write();
    ANNZResidualVsConfDifferenceTgt4Lead->Write();
    ANNZResidualVsConfDifferenceTgt5Iron->Write();
    ANNZResidualVsConfDifferenceTgt5Lead->Write();
    ANNZResidualVsConfDifferenceWater->Write();

    ANNPlaneProbabilityVsEhadData->Write();
    ANNPlaneProbabilityTrackerVsEhadData->Write();
    ANNPlaneProbabilityTgt2IronVsEhadData->Write();
    ANNPlaneProbabilityTgt2LeadVsEhadData->Write();
    ANNPlaneProbabilityTgt3IronVsEhadData->Write();
    ANNPlaneProbabilityTgt3LeadVsEhadData->Write();
    ANNPlaneProbabilityTgt3CarbonVsEhadData->Write();
    ANNPlaneProbabilityTgt4LeadVsEhadData->Write();
    ANNPlaneProbabilityTgt5IronVsEhadData->Write();
    ANNPlaneProbabilityTgt5LeadVsEhadData->Write();
    ANNPlaneProbabilityWaterVsEhadData->Write();

    ANNPlaneProbabilityVsEhadMC->Write();
    ANNPlaneProbabilityTrackerVsEhadMC->Write();
    ANNPlaneProbabilityTgt2IronVsEhadMC->Write();
    ANNPlaneProbabilityTgt2LeadVsEhadMC->Write();
    ANNPlaneProbabilityTgt3IronVsEhadMC->Write();
    ANNPlaneProbabilityTgt3LeadVsEhadMC->Write();
    ANNPlaneProbabilityTgt3CarbonVsEhadMC->Write();
    ANNPlaneProbabilityTgt4LeadVsEhadMC->Write();
    ANNPlaneProbabilityTgt5IronVsEhadMC->Write();
    ANNPlaneProbabilityTgt5LeadVsEhadMC->Write();
    ANNPlaneProbabilityWaterVsEhadMC->Write();

    //ERecoil
    ErecoilMC->Write();
    ErecoilData->Write();
    TruthVerticesMCERecoil_ByModule->Write();
    ANNVerticesMCERecoil_ByModule->Write();
    TBVerticesMCERecoil_ByModule->Write();
    ANNVerticesDataERecoil_ByModule->Write();
    TBVerticesDataERecoil_ByModule->Write();


    ANNVerticesData_ByZPosVsERecoil->Write();
    ANNVerticesMC_ByZPosVsERecoil->Write();


    ANNVerticesMC_BySegmentVsERecoil->Write();
    TruthVerticesMC_BySegmentVsERecoil->Write();
    ANNVerticesData_BySegmentVsERecoil->Write();


    //Multiplicity
    TruthVerticesMCMultiplicity_ByModule->Write();
    ANNVerticesMCMultiplicity_ByModule->Write();
    TBVerticesMCMultiplicity_ByModule->Write();
    ANNVerticesDataMultiplicity_ByModule->Write();
    TBVerticesDataMultiplicity_ByModule->Write();

    //ANN Confidence
    TruthVerticesMCANNConf_ByModule->Write();
    ANNVerticesMCANNConf_ByModule->Write();
    ANNVerticesDataANNConf_ByModule->Write();

    ANNVerticesData_ByZPosVsANNConf->Write();
    ANNVerticesMC_ByZPosVsANNConf->Write();

    //For events misreconstructed in each target, what's the probability distribution
    ANNConfMisrecoInTgt2Iron->Write();
    ANNConfMisrecoInTgt2Lead->Write();
    ANNConfMisrecoInTgt3Iron->Write();
    ANNConfMisrecoInTgt3Lead->Write();
    ANNConfMisrecoInTgt3Carbon->Write();
    ANNConfMisrecoInTgt4Lead->Write();
    ANNConfMisrecoInTgt5Iron->Write();
    ANNConfMisrecoInTgt5Lead->Write();
    ANNConfMisrecoInWater->Write();

    ANNPlaneProbabilityVsPmuMC->Write();
    ANNPlaneProbabilityVsPmuData->Write();

    //Efficiency as a function of ANNConfidence and EHad bin

    ANNPlaneProbabilityVsEhadNumerator->Write();
    ANNPlaneProbabilityTrackerVsEhadNumerator->Write();
    ANNPlaneProbabilityTgt2IronVsEhadNumerator->Write();
    ANNPlaneProbabilityTgt2LeadVsEhadNumerator->Write();
    ANNPlaneProbabilityTgt3IronVsEhadNumerator->Write();
    ANNPlaneProbabilityTgt3LeadVsEhadNumerator->Write();
    ANNPlaneProbabilityTgt3CarbonVsEhadNumerator->Write();
    ANNPlaneProbabilityTgt4LeadVsEhadNumerator->Write();
    ANNPlaneProbabilityTgt5IronVsEhadNumerator->Write();
    ANNPlaneProbabilityTgt5LeadVsEhadNumerator->Write();
    ANNPlaneProbabilityWaterVsEhadNumerator->Write();

    ANNPlaneProbabilityVsEhadDenominator->Write();
    ANNPlaneProbabilityTrackerVsEhadDenominator->Write();
    ANNPlaneProbabilityTgt2IronVsEhadDenominator->Write();
    ANNPlaneProbabilityTgt2LeadVsEhadDenominator->Write();
    ANNPlaneProbabilityTgt3IronVsEhadDenominator->Write();
    ANNPlaneProbabilityTgt3LeadVsEhadDenominator->Write();
    ANNPlaneProbabilityTgt3CarbonVsEhadDenominator->Write();
    ANNPlaneProbabilityTgt4LeadVsEhadDenominator->Write();
    ANNPlaneProbabilityTgt5IronVsEhadDenominator->Write();
    ANNPlaneProbabilityTgt5LeadVsEhadDenominator->Write();
    ANNPlaneProbabilityWaterVsEhadDenominator->Write();


    //Curvature signifiance
    TruthVerticesMCCurvSig_ByModule->Write();
    ANNVerticesMCCurvSig_ByModule->Write();
    TBVerticesMCCurvSig_ByModule->Write();
    ANNVerticesDataCurvSig_ByModule->Write();
    TBVerticesDataCurvSig_ByModule->Write();

    //Migration ANN only
    ANNVerticesConfusion_ByModule->Write();
    ANNVerticesConfusion_ByZPos->Write();


    double potMC = options.m_mc_pot ;
    double potData = options.m_data_pot ;
    if (nSubruns>0)
    {
      potMC/=nSubruns;
      potData/=nSubruns;
    }
    auto MCPOT = new TParameter<double>("MCPOT", (potMC));
    MCPOT->Write();
    auto DataPOT = new TParameter<double>("DataPOT", (potData));
    DataPOT->Write();
    OutDir->Close();


    //==============================================================================
    // Drawing Plots
    //==============================================================================

    //DrawDataMC(ANNVerticesFilledMC_ByModule, TruthVerticesFilledMC_ByModule, "Test", optionsFilled.m_mc_pot, optionsFilled.m_data_pot);
    //DrawFullAndEmptyDataMC(ANNVerticesFilledData_ByZPos, ANNVerticesFilledMC_ByZPos, ANNVerticesEmptyData_ByZPos, ANNVerticesEmptyMC_ByZPos, "Test", 0, 0, 0, 0);
    //==============================================================================
    // End - Drawing Plots
    //==============================================================================
    std::cout << "Success" << std::endl;
  }
  /* catch(const ROOT::exception& e)
  {
    std::cerr << "Ending on a ROOT error message.  No histograms will be produced.\n"
              << "If the message talks about \"TNetXNGFile\", this could be a problem with dCache.  The message is:\n"
              << e.what() << "\n" << USAGE << "\n";
    return badFileRead;
  } */

  return success;
}

void TargetRegionValidationNew(std::string datafile, std::string mcfile)
{

    std::vector<std::string> strings {"one", "two", "three"};
    std::vector<char*> cstrings;
    cstrings.reserve(strings.size());
    cstrings.push_back(const_cast<char*>(datafile.c_str()));
    cstrings.push_back(const_cast<char*>(mcfile.c_str()));
    char ** instr = &cstrings[0];
  main(2, instr);
}