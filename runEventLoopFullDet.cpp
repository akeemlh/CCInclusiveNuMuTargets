#define MC_OUT_FILE_NAME "runEventLoopFullDetMC.root"
#define DATA_OUT_FILE_NAME "runEventLoopFullDetData.root"

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
"Produces a two files, " MC_OUT_FILE_NAME " and " DATA_OUT_FILE_NAME ", with\n"\
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
#include "cuts/MaxPzMu.h"
#include "util/Variable.h"
#include "util/Variable2D.h"
#include "util/GetFluxIntegral.h"
#include "util/GetPlaylist.h"
#include "cuts/SignalDefinition.h"
#include "cuts/q3RecoCut.h"
#include "studies/Study.h"
//#include "Binning.h" //TODO: Fix me

//PlotUtils includes
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/MacroUtil.h"
#include "PlotUtils/MnvPlotter.h"
#include "PlotUtils/CCInclusiveCuts.h"
#include "PlotUtils/CCInclusiveSignal.h"
#include "PlotUtils/CrashOnROOTMessage.h" //Sets up ROOT's debug callbacks by itself
#include "PlotUtils/Cutter.h"
#include "PlotUtils/Model.h"
#include "PlotUtils/FluxAndCVReweighter.h"
#include "PlotUtils/GENIEReweighter.h"
#include "PlotUtils/LowRecoil2p2hReweighter.h"
#include "PlotUtils/RPAReweighter.h"
#include "PlotUtils/MINOSEfficiencyReweighter.h"
#include "PlotUtils/TargetUtils.h"
#pragma GCC diagnostic pop

//ROOT includes
#include "TParameter.h"

#include "Math/Vector3D.h"
#include "TH3D.h"

//c++ includes
#include <iostream>
#include <cstdlib> //getenv()

std::vector<double> vertexBins = {4293.04, 4337.25, 4381.47, 4425.68, 4514.11, 4558.33, 4602.54, 4646.76, 4735.19, 4779.4, 4823.62, 4867.83, 5000.48, 5044.69, 5088.91, 5133.12, 5456.74, 5500.95, 5545.17, 5589.38, 5677.81, 5722.03, 5810.45, 5855.68, 5900.91, 5946.14, 5991.37, 6036.6, 6081.83, 6127.06, 6172.29, 6217.52, 6262.74, 6307.97, 6353.2, 6398.43, 6443.66, 6488.89, 6534.12, 6579.35, 6624.58, 6669.81, 6715.03, 6760.26, 6805.49, 6850.72, 6895.95, 6941.18, 6986.41, 7031.64, 7076.87, 7122.1, 7167.32, 7212.55, 7257.78, 7303.01, 7348.24, 7393.47, 7438.7, 7483.93, 7529.16, 7574.39, 7619.61, 7664.84, 7710.07, 7755.3, 7800.53, 7845.76, 7890.99, 7936.22, 7981.45, 8026.68, 8071.9, 8117.13, 8162.36, 8207.59, 8252.82, 8298.05, 8343.28, 8388.51, 8433.74, 8478.97, 8524.19, 8569.42, 8614.65};

TH1D *ANNVerticesMC = new TH1D ("ANNVerticesMC", "ANNVerticesMC", vertexBins.size()-1, &vertexBins[0]);
TH1D *TBVerticesMC = new TH1D ("TBVerticesMC", "TBVerticesMC", vertexBins.size()-1, &vertexBins[0]);
TH1D *ANNVerticesData = new TH1D ("ANNVerticesData", "ANNVerticesData", vertexBins.size()-1, &vertexBins[0]);
TH1D *TBVerticesData = new TH1D ("TBVerticesData", "TBVerticesData", vertexBins.size()-1, &vertexBins[0]);

TH3D *TBVerticesGranularData = new TH3D ("TBVerticesGranularData", "TBVerticesGranularData", 100, -1000, 1000, 100, -1000, 1000, 3400, 4200, 5900);
TH3D *TBVerticesGranularMC = new TH3D ("TBVerticesGranularMC", "TBVerticesGranularMC", 100, -1000, 1000, 100, -1000, 1000, 3400, 4200, 5900);
TH3D *ANNVerticesGranularData = new TH3D ("ANNVerticesGranularData", "ANNVerticesGranularData", 100, -1000, 1000, 100, -1000, 1000, 3400, 4200, 5900);
TH3D *ANNVerticesGranularMC = new TH3D ("ANNVerticesGranularMC", "ANNVerticesGranularMC", 100, -1000, 1000, 100, -1000, 1000, 3400, 4200, 5900);
TH3D *TrueVerticesGranular = new TH3D ("TrueVerticesGranular", "TrueVerticesGranular", 100, -1000, 1000, 100, -1000, 1000, 3400, 4200, 5900);

TH3D *TrueVtxANNRecoInWater = new TH3D ("TrueVtxANNRecoInWater", "TrueVtxANNRecoInWater", 100, -1000, 1000, 100, -1000, 1000, 3400, 4200, 5900);
TH3D *TrueVtxTBRecoInWater = new TH3D ("TrueVtxTBRecoInWater", "TrueVtxTBRecoInWater", 100, -1000, 1000, 100, -1000, 1000, 3400, 4200, 5900);
TH3D *TrueVtxANNRecoOutWater = new TH3D ("TrueVtxANNRecoOutWater", "TrueVtxANNRecoOutWater", 100, -1000, 1000, 100, -1000, 1000, 3400, 4200, 5900);
TH3D *TrueVtxTBRecoOutWater = new TH3D ("TrueVtxTBRecoOutWater", "TrueVtxTBRecoOutWater", 100, -1000, 1000, 100, -1000, 1000, 3400, 4200, 5900);

TH1D *ANNRecoInWater = new TH1D ("ANNRecoInWater", "ANNRecoInWater", 3400, 4200, 5900);
TH1D *TBRecoInWater = new TH1D ("TBRecoInWater", "TBRecoInWater", 3400, 4200, 5900);
TH1D *ANNRecoOutWater = new TH1D ("ANNRecoOutWater", "ANNRecoOutWater", 3400, 4200, 5900);
TH1D *TBRecoOutWater = new TH1D ("TBRecoOutWater", "TBRecoOutWater", 3400, 4200, 5900);

TH1D *ANNRecoInTarget1 = new TH1D ("ANNRecoInTarget1", "ANNRecoInTarget1", 3400, 4200, 5900);
TH1D *TBRecoInTarget1 = new TH1D ("TBRecoInTarget1", "TBRecoInTarget1", 3400, 4200, 5900);
TH1D *ANNRecoOutTarget1 = new TH1D ("ANNRecoOutTarget1", "ANNRecoOutTarget1", 3400, 4200, 5900);
TH1D *TBRecoOutTarget1 = new TH1D ("TBRecoOutTarget1", "TBRecoOutTarget1", 3400, 4200, 5900);

TH1D *ANNRecoInTarget2 = new TH1D ("ANNRecoInTarget2", "ANNRecoInTarget2", 3400, 4200, 5900);
TH1D *TBRecoInTarget2 = new TH1D ("TBRecoInTarget2", "TBRecoInTarget2", 3400, 4200, 5900);
TH1D *ANNRecoOutTarget2 = new TH1D ("ANNRecoOutTarget2", "ANNRecoOutTarget2", 3400, 4200, 5900);
TH1D *TBRecoOutTarget2 = new TH1D ("TBRecoOutTarget2", "TBRecoOutTarget2", 3400, 4200, 5900);

TH1D *ANNRecoInTarget3 = new TH1D ("ANNRecoInTarget3", "ANNRecoInTarget3", 3400, 4200, 5900);
TH1D *TBRecoInTarget3 = new TH1D ("TBRecoInTarget3", "TBRecoInTarget3", 3400, 4200, 5900);
TH1D *ANNRecoOutTarget3 = new TH1D ("ANNRecoOutTarget3", "ANNRecoOutTarget3", 3400, 4200, 5900);
TH1D *TBRecoOutTarget3 = new TH1D ("TBRecoOutTarget3", "TBRecoOutTarget3", 3400, 4200, 5900);

TH1D *ANNRecoInTarget4 = new TH1D ("ANNRecoInTarget4", "ANNRecoInTarget4", 3400, 4200, 5900);
TH1D *TBRecoInTarget4 = new TH1D ("TBRecoInTarget4", "TBRecoInTarget4", 3400, 4200, 5900);
TH1D *ANNRecoOutTarget4 = new TH1D ("ANNRecoOutTarget4", "ANNRecoOutTarget4", 3400, 4200, 5900);
TH1D *TBRecoOutTarget4 = new TH1D ("TBRecoOutTarget4", "TBRecoOutTarget4", 3400, 4200, 5900);

TH1D *ANNRecoInTarget5 = new TH1D ("ANNRecoInTarget5", "ANNRecoInTarget5", 3400, 4200, 5900);
TH1D *TBRecoInTarget5 = new TH1D ("TBRecoInTarget5", "TBRecoInTarget5", 3400, 4200, 5900);
TH1D *ANNRecoOutTarget5 = new TH1D ("ANNRecoOutTarget5", "ANNRecoOutTarget5", 3400, 4200, 5900);
TH1D *TBRecoOutTarget5 = new TH1D ("TBRecoOutTarget5", "TBRecoOutTarget5", 3400, 4200, 5900);

TH1D *ANNRecoInTarget1Iron = new TH1D ("ANNRecoInTarget1Iron", "ANNRecoInTarget1Iron", 3400, 4200, 5900);
TH1D *TBRecoInTarget1Iron = new TH1D ("TBRecoInTarget1Iron", "TBRecoInTarget1Iron", 3400, 4200, 5900);
TH1D *ANNRecoOutTarget1Iron = new TH1D ("ANNRecoOutTarget1Iron", "ANNRecoOutTarget1Iron", 3400, 4200, 5900);
TH1D *TBRecoOutTarget1Iron = new TH1D ("TBRecoOutTarget1Iron", "TBRecoOutTarget1Iron", 3400, 4200, 5900);

TH1D *ANNRecoInTarget1Lead = new TH1D ("ANNRecoInTarget1Lead", "ANNRecoInTarget1Lead", 3400, 4200, 5900);
TH1D *TBRecoInTarget1Lead = new TH1D ("TBRecoInTarget1Lead", "TBRecoInTarget1Lead", 3400, 4200, 5900);
TH1D *ANNRecoOutTarget1Lead = new TH1D ("ANNRecoOutTarget1Lead", "ANNRecoOutTarget1Lead", 3400, 4200, 5900);
TH1D *TBRecoOutTarget1Lead = new TH1D ("TBRecoOutTarget1Lead", "TBRecoOutTarget1Lead", 3400, 4200, 5900);

TH1D *ANNRecoInTarget2Iron = new TH1D ("ANNRecoInTarget2Iron", "ANNRecoInTarget2Iron", 3400, 4200, 5900);
TH1D *TBRecoInTarget2Iron = new TH1D ("TBRecoInTarget2Iron", "TBRecoInTarget2Iron", 3400, 4200, 5900);
TH1D *ANNRecoOutTarget2Iron = new TH1D ("ANNRecoOutTarget2Iron", "ANNRecoOutTarget2Iron", 3400, 4200, 5900);
TH1D *TBRecoOutTarget2Iron = new TH1D ("TBRecoOutTarget2Iron", "TBRecoOutTarget2Iron", 3400, 4200, 5900);

TH1D *ANNRecoInTarget2Lead = new TH1D ("ANNRecoInTarget2Lead", "ANNRecoInTarget2Lead", 3400, 4200, 5900);
TH1D *TBRecoInTarget2Lead = new TH1D ("TBRecoInTarget2Lead", "TBRecoInTarget2Lead", 3400, 4200, 5900);
TH1D *ANNRecoOutTarget2Lead = new TH1D ("ANNRecoOutTarget2Lead", "ANNRecoOutTarget2Lead", 3400, 4200, 5900);
TH1D *TBRecoOutTarget2Lead = new TH1D ("TBRecoOutTarget2Lead", "TBRecoOutTarget2Lead", 3400, 4200, 5900);

TH1D *ANNRecoInTarget3Iron = new TH1D ("ANNRecoInTarget3Iron", "ANNRecoInTarget3Iron", 3400, 4200, 5900);
TH1D *TBRecoInTarget3Iron = new TH1D ("TBRecoInTarget3Iron", "TBRecoInTarget3Iron", 3400, 4200, 5900);
TH1D *ANNRecoOutTarget3Iron = new TH1D ("ANNRecoOutTarget3Iron", "ANNRecoOutTarget3Iron", 3400, 4200, 5900);
TH1D *TBRecoOutTarget3Iron = new TH1D ("TBRecoOutTarget3Iron", "TBRecoOutTarget3Iron", 3400, 4200, 5900);

TH1D *ANNRecoInTarget3Lead = new TH1D ("ANNRecoInTarget3Lead", "ANNRecoInTarget3Lead", 3400, 4200, 5900);
TH1D *TBRecoInTarget3Lead = new TH1D ("TBRecoInTarget3Lead", "TBRecoInTarget3Lead", 3400, 4200, 5900);
TH1D *ANNRecoOutTarget3Lead = new TH1D ("ANNRecoOutTarget3Lead", "ANNRecoOutTarget3Lead", 3400, 4200, 5900);
TH1D *TBRecoOutTarget3Lead = new TH1D ("TBRecoOutTarget3Lead", "TBRecoOutTarget3Lead", 3400, 4200, 5900);

TH1D *ANNRecoInTarget3Carbon = new TH1D ("ANNRecoInTarget3Carbon", "ANNRecoInTarget3Carbon", 3400, 4200, 5900);
TH1D *TBRecoInTarget3Carbon = new TH1D ("TBRecoInTarget3Carbon", "TBRecoInTarget3Carbon", 3400, 4200, 5900);
TH1D *ANNRecoOutTarget3Carbon = new TH1D ("ANNRecoOutTarget3Carbon", "ANNRecoOutTarget3Carbon", 3400, 4200, 5900);
TH1D *TBRecoOutTarget3Carbon = new TH1D ("TBRecoOutTarget3Carbon", "TBRecoOutTarget3Carbon", 3400, 4200, 5900);

TH1D *ANNRecoInTarget4Lead = new TH1D ("ANNRecoInTarget4Lead", "ANNRecoInTarget4Lead", 3400, 4200, 5900);
TH1D *TBRecoInTarget4Lead = new TH1D ("TBRecoInTarget4Lead", "TBRecoInTarget4Lead", 3400, 4200, 5900);
TH1D *ANNRecoOutTarget4Lead = new TH1D ("ANNRecoOutTarget4Lead", "ANNRecoOutTarget4Lead", 3400, 4200, 5900);
TH1D *TBRecoOutTarget4Lead = new TH1D ("TBRecoOutTarget4Lead", "TBRecoOutTarget4Lead", 3400, 4200, 5900);

TH1D *ANNRecoInTarget5Lead = new TH1D ("ANNRecoInTarget5Lead", "ANNRecoInTarget5Lead", 3400, 4200, 5900);
TH1D *TBRecoInTarget5Lead = new TH1D ("TBRecoInTarget5Lead", "TBRecoInTarget5Lead", 3400, 4200, 5900);
TH1D *ANNRecoOutTarget5Lead = new TH1D ("ANNRecoOutTarget5Lead", "ANNRecoOutTarget5Lead", 3400, 4200, 5900);
TH1D *TBRecoOutTarget5Lead = new TH1D ("TBRecoOutTarget5Lead", "TBRecoOutTarget5Lead", 3400, 4200, 5900);

TH1D *ANNRecoInTarget5Iron = new TH1D ("ANNRecoInTarget5Iron", "ANNRecoInTarget5Iron", 3400, 4200, 5900);
TH1D *TBRecoInTarget5Iron = new TH1D ("TBRecoInTarget5Iron", "TBRecoInTarget5Iron", 3400, 4200, 5900);
TH1D *ANNRecoOutTarget5Iron = new TH1D ("ANNRecoOutTarget5Iron", "ANNRecoOutTarget5Iron", 3400, 4200, 5900);
TH1D *TBRecoOutTarget5Iron = new TH1D ("TBRecoOutTarget5Iron", "TBRecoOutTarget5Iron", 3400, 4200, 5900);

TH1D *ANNRecoInIron = new TH1D ("ANNRecoInIron", "ANNRecoInIron", 3400, 4200, 5900);
TH1D *TBRecoInIron = new TH1D ("TBRecoInIron", "TBRecoInIron", 3400, 4200, 5900);
TH1D *ANNRecoOutIron = new TH1D ("ANNRecoOutIron", "ANNRecoOutIron", 3400, 4200, 5900);
TH1D *TBRecoOutIron = new TH1D ("TBRecoOutIron", "TBRecoOutIron", 3400, 4200, 5900);

TH1D *ANNRecoInLead = new TH1D ("ANNRecoInLead", "ANNRecoInLead", 3400, 4200, 5900);
TH1D *TBRecoInLead = new TH1D ("TBRecoInLead", "TBRecoInLead", 3400, 4200, 5900);
TH1D *ANNRecoOutLead = new TH1D ("ANNRecoOutLead", "ANNRecoOutLead", 3400, 4200, 5900);
TH1D *TBRecoOutLead = new TH1D ("TBRecoOutLead", "TBRecoOutLead", 3400, 4200, 5900);

TH1D *ANNRecoInCarbon = new TH1D ("ANNRecoInCarbon", "ANNRecoInCarbon", 3400, 4200, 5900);
TH1D *TBRecoInCarbon = new TH1D ("TBRecoInCarbon", "TBRecoInCarbon", 3400, 4200, 5900);
TH1D *ANNRecoOutCarbon = new TH1D ("ANNRecoOutCarbon", "ANNRecoOutCarbon", 3400, 4200, 5900);
TH1D *TBRecoOutCarbon = new TH1D ("TBRecoOutCarbon", "TBRecoOutCarbon", 3400, 4200, 5900);







TH1D *ANNTruthInWater = new TH1D ("ANNTruthInWater", "ANNTruthInWater", 3400, 4200, 5900);
TH1D *TBTruthInWater = new TH1D ("TBTruthInWater", "TBTruthInWater", 3400, 4200, 5900);
TH1D *ANNTruthOutWater = new TH1D ("ANNTruthOutWater", "ANNTruthOutWater", 3400, 4200, 5900);
TH1D *TBTruthOutWater = new TH1D ("TBTruthOutWater", "TBTruthOutWater", 3400, 4200, 5900);

TH1D *ANNTruthInTarget1 = new TH1D ("ANNTruthInTarget1", "ANNTruthInTarget1", 3400, 4200, 5900);
TH1D *TBTruthInTarget1 = new TH1D ("TBTruthInTarget1", "TBTruthInTarget1", 3400, 4200, 5900);
TH1D *ANNTruthOutTarget1 = new TH1D ("ANNTruthOutTarget1", "ANNTruthOutTarget1", 3400, 4200, 5900);
TH1D *TBTruthOutTarget1 = new TH1D ("TBTruthOutTarget1", "TBTruthOutTarget1", 3400, 4200, 5900);

TH1D *ANNTruthInTarget2 = new TH1D ("ANNTruthInTarget2", "ANNTruthInTarget2", 3400, 4200, 5900);
TH1D *TBTruthInTarget2 = new TH1D ("TBTruthInTarget2", "TBTruthInTarget2", 3400, 4200, 5900);
TH1D *ANNTruthOutTarget2 = new TH1D ("ANNTruthOutTarget2", "ANNTruthOutTarget2", 3400, 4200, 5900);
TH1D *TBTruthOutTarget2 = new TH1D ("TBTruthOutTarget2", "TBTruthOutTarget2", 3400, 4200, 5900);

TH1D *ANNTruthInTarget3 = new TH1D ("ANNTruthInTarget3", "ANNTruthInTarget3", 3400, 4200, 5900);
TH1D *TBTruthInTarget3 = new TH1D ("TBTruthInTarget3", "TBTruthInTarget3", 3400, 4200, 5900);
TH1D *ANNTruthOutTarget3 = new TH1D ("ANNTruthOutTarget3", "ANNTruthOutTarget3", 3400, 4200, 5900);
TH1D *TBTruthOutTarget3 = new TH1D ("TBTruthOutTarget3", "TBTruthOutTarget3", 3400, 4200, 5900);

TH1D *ANNTruthInTarget4 = new TH1D ("ANNTruthInTarget4", "ANNTruthInTarget4", 3400, 4200, 5900);
TH1D *TBTruthInTarget4 = new TH1D ("TBTruthInTarget4", "TBTruthInTarget4", 3400, 4200, 5900);
TH1D *ANNTruthOutTarget4 = new TH1D ("ANNTruthOutTarget4", "ANNTruthOutTarget4", 3400, 4200, 5900);
TH1D *TBTruthOutTarget4 = new TH1D ("TBTruthOutTarget4", "TBTruthOutTarget4", 3400, 4200, 5900);

TH1D *ANNTruthInTarget5 = new TH1D ("ANNTruthInTarget5", "ANNTruthInTarget5", 3400, 4200, 5900);
TH1D *TBTruthInTarget5 = new TH1D ("TBTruthInTarget5", "TBTruthInTarget5", 3400, 4200, 5900);
TH1D *ANNTruthOutTarget5 = new TH1D ("ANNTruthOutTarget5", "ANNTruthOutTarget5", 3400, 4200, 5900);
TH1D *TBTruthOutTarget5 = new TH1D ("TBTruthOutTarget5", "TBTruthOutTarget5", 3400, 4200, 5900);

TH1D *ANNTruthInTarget1Iron = new TH1D ("ANNTruthInTarget1Iron", "ANNTruthInTarget1Iron", 3400, 4200, 5900);
TH1D *TBTruthInTarget1Iron = new TH1D ("TBTruthInTarget1Iron", "TBTruthInTarget1Iron", 3400, 4200, 5900);
TH1D *ANNTruthOutTarget1Iron = new TH1D ("ANNTruthOutTarget1Iron", "ANNTruthOutTarget1Iron", 3400, 4200, 5900);
TH1D *TBTruthOutTarget1Iron = new TH1D ("TBTruthOutTarget1Iron", "TBTruthOutTarget1Iron", 3400, 4200, 5900);

TH1D *ANNTruthInTarget1Lead = new TH1D ("ANNTruthInTarget1Lead", "ANNTruthInTarget1Lead", 3400, 4200, 5900);
TH1D *TBTruthInTarget1Lead = new TH1D ("TBTruthInTarget1Lead", "TBTruthInTarget1Lead", 3400, 4200, 5900);
TH1D *ANNTruthOutTarget1Lead = new TH1D ("ANNTruthOutTarget1Lead", "ANNTruthOutTarget1Lead", 3400, 4200, 5900);
TH1D *TBTruthOutTarget1Lead = new TH1D ("TBTruthOutTarget1Lead", "TBTruthOutTarget1Lead", 3400, 4200, 5900);

TH1D *ANNTruthInTarget2Iron = new TH1D ("ANNTruthInTarget2Iron", "ANNTruthInTarget2Iron", 3400, 4200, 5900);
TH1D *TBTruthInTarget2Iron = new TH1D ("TBTruthInTarget2Iron", "TBTruthInTarget2Iron", 3400, 4200, 5900);
TH1D *ANNTruthOutTarget2Iron = new TH1D ("ANNTruthOutTarget2Iron", "ANNTruthOutTarget2Iron", 3400, 4200, 5900);
TH1D *TBTruthOutTarget2Iron = new TH1D ("TBTruthOutTarget2Iron", "TBTruthOutTarget2Iron", 3400, 4200, 5900);

TH1D *ANNTruthInTarget2Lead = new TH1D ("ANNTruthInTarget2Lead", "ANNTruthInTarget2Lead", 3400, 4200, 5900);
TH1D *TBTruthInTarget2Lead = new TH1D ("TBTruthInTarget2Lead", "TBTruthInTarget2Lead", 3400, 4200, 5900);
TH1D *ANNTruthOutTarget2Lead = new TH1D ("ANNTruthOutTarget2Lead", "ANNTruthOutTarget2Lead", 3400, 4200, 5900);
TH1D *TBTruthOutTarget2Lead = new TH1D ("TBTruthOutTarget2Lead", "TBTruthOutTarget2Lead", 3400, 4200, 5900);

TH1D *ANNTruthInTarget3Iron = new TH1D ("ANNTruthInTarget3Iron", "ANNTruthInTarget3Iron", 3400, 4200, 5900);
TH1D *TBTruthInTarget3Iron = new TH1D ("TBTruthInTarget3Iron", "TBTruthInTarget3Iron", 3400, 4200, 5900);
TH1D *ANNTruthOutTarget3Iron = new TH1D ("ANNTruthOutTarget3Iron", "ANNTruthOutTarget3Iron", 3400, 4200, 5900);
TH1D *TBTruthOutTarget3Iron = new TH1D ("TBTruthOutTarget3Iron", "TBTruthOutTarget3Iron", 3400, 4200, 5900);

TH1D *ANNTruthInTarget3Lead = new TH1D ("ANNTruthInTarget3Lead", "ANNTruthInTarget3Lead", 3400, 4200, 5900);
TH1D *TBTruthInTarget3Lead = new TH1D ("TBTruthInTarget3Lead", "TBTruthInTarget3Lead", 3400, 4200, 5900);
TH1D *ANNTruthOutTarget3Lead = new TH1D ("ANNTruthOutTarget3Lead", "ANNTruthOutTarget3Lead", 3400, 4200, 5900);
TH1D *TBTruthOutTarget3Lead = new TH1D ("TBTruthOutTarget3Lead", "TBTruthOutTarget3Lead", 3400, 4200, 5900);

TH1D *ANNTruthInTarget3Carbon = new TH1D ("ANNTruthInTarget3Carbon", "ANNTruthInTarget3Carbon", 3400, 4200, 5900);
TH1D *TBTruthInTarget3Carbon = new TH1D ("TBTruthInTarget3Carbon", "TBTruthInTarget3Carbon", 3400, 4200, 5900);
TH1D *ANNTruthOutTarget3Carbon = new TH1D ("ANNTruthOutTarget3Carbon", "ANNTruthOutTarget3Carbon", 3400, 4200, 5900);
TH1D *TBTruthOutTarget3Carbon = new TH1D ("TBTruthOutTarget3Carbon", "TBTruthOutTarget3Carbon", 3400, 4200, 5900);

TH1D *ANNTruthInTarget4Lead = new TH1D ("ANNTruthInTarget4Lead", "ANNTruthInTarget4Lead", 3400, 4200, 5900);
TH1D *TBTruthInTarget4Lead = new TH1D ("TBTruthInTarget4Lead", "TBTruthInTarget4Lead", 3400, 4200, 5900);
TH1D *ANNTruthOutTarget4Lead = new TH1D ("ANNTruthOutTarget4Lead", "ANNTruthOutTarget4Lead", 3400, 4200, 5900);
TH1D *TBTruthOutTarget4Lead = new TH1D ("TBTruthOutTarget4Lead", "TBTruthOutTarget4Lead", 3400, 4200, 5900);

TH1D *ANNTruthInTarget5Lead = new TH1D ("ANNTruthInTarget5Lead", "ANNTruthInTarget5Lead", 3400, 4200, 5900);
TH1D *TBTruthInTarget5Lead = new TH1D ("TBTruthInTarget5Lead", "TBTruthInTarget5Lead", 3400, 4200, 5900);
TH1D *ANNTruthOutTarget5Lead = new TH1D ("ANNTruthOutTarget5Lead", "ANNTruthOutTarget5Lead", 3400, 4200, 5900);
TH1D *TBTruthOutTarget5Lead = new TH1D ("TBTruthOutTarget5Lead", "TBTruthOutTarget5Lead", 3400, 4200, 5900);

TH1D *ANNTruthInTarget5Iron = new TH1D ("ANNTruthInTarget5Iron", "ANNTruthInTarget5Iron", 3400, 4200, 5900);
TH1D *TBTruthInTarget5Iron = new TH1D ("TBTruthInTarget5Iron", "TBTruthInTarget5Iron", 3400, 4200, 5900);
TH1D *ANNTruthOutTarget5Iron = new TH1D ("ANNTruthOutTarget5Iron", "ANNTruthOutTarget5Iron", 3400, 4200, 5900);
TH1D *TBTruthOutTarget5Iron = new TH1D ("TBTruthOutTarget5Iron", "TBTruthOutTarget5Iron", 3400, 4200, 5900);

TH1D *ANNTruthInIron = new TH1D ("ANNTruthInIron", "ANNTruthInIron", 3400, 4200, 5900);
TH1D *TBTruthInIron = new TH1D ("TBTruthInIron", "TBTruthInIron", 3400, 4200, 5900);
TH1D *ANNTruthOutIron = new TH1D ("ANNTruthOutIron", "ANNTruthOutIron", 3400, 4200, 5900);
TH1D *TBTruthOutIron = new TH1D ("TBTruthOutIron", "TBTruthOutIron", 3400, 4200, 5900);

TH1D *ANNTruthInLead = new TH1D ("ANNTruthInLead", "ANNTruthInLead", 3400, 4200, 5900);
TH1D *TBTruthInLead = new TH1D ("TBTruthInLead", "TBTruthInLead", 3400, 4200, 5900);
TH1D *ANNTruthOutLead = new TH1D ("ANNTruthOutLead", "ANNTruthOutLead", 3400, 4200, 5900);
TH1D *TBTruthOutLead = new TH1D ("TBTruthOutLead", "TBTruthOutLead", 3400, 4200, 5900);

TH1D *ANNTruthInCarbon = new TH1D ("ANNTruthInCarbon", "ANNTruthInCarbon", 3400, 4200, 5900);
TH1D *TBTruthInCarbon = new TH1D ("TBTruthInCarbon", "TBTruthInCarbon", 3400, 4200, 5900);
TH1D *ANNTruthOutCarbon = new TH1D ("ANNTruthOutCarbon", "ANNTruthOutCarbon", 3400, 4200, 5900);
TH1D *TBTruthOutCarbon = new TH1D ("TBTruthOutCarbon", "TBTruthOutCarbon", 3400, 4200, 5900);




//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFillEventSelection(
    PlotUtils::ChainWrapper* chain,
    std::map<std::string, std::vector<CVUniverse*> > error_bands,
    std::vector<Variable*> vars,
    std::vector<Variable2D*> vars2D,
    std::vector<Study*> studies,
    PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts,
    PlotUtils::Model<CVUniverse, MichelEvent>& model)
{
  assert(!error_bands["cv"].empty() && "\"cv\" error band is empty!  Can't set Model weight.");
  auto& cvUniv = error_bands["cv"].front();

  std::cout << "Starting MC reco loop...\n";
  const int nEntries = chain->GetEntries();
  for (int i=0; i<nEntries; ++i)
  {
    if(i%1000==0) std::cout << i << " / " << nEntries << "\r" <<std::endl;

    MichelEvent cvEvent;
    cvUniv->SetEntry(i);
    model.SetEntry(*cvUniv, cvEvent);
    const double cvWeight = model.GetWeight(*cvUniv, cvEvent);

    //=========================================
    // Systematics loop(s)
    //=========================================
    MichelEvent myevent; // make sure your event is inside the error band loop. 

    std::vector<double> ANNVtx = cvUniv->GetANNVertexVector();
    ROOT::Math::XYZTVector TrackBasedVtx = cvUniv->GetVertex();
    ROOT::Math::XYZTVector TrueVtx = cvUniv->GetTrueVertex();

    // This is where you would Access/create a Michel

    //weight is ignored in isMCSelected() for all but the CV Universe.
    if (!michelcuts.isMCSelected(*cvUniv, myevent, cvWeight).all()) continue; //all is another function that will later help me with sidebands
    PlotUtils::TargetUtils* m_TargetUtils=new PlotUtils::TargetUtils();
    //Performing vtx validation check Deborah suggested
    //if(cvUniv->hasMLPred() && cvUniv->GetANNProb()>0.2)
    double batchPOT = cvUniv->GetBatchPOT();
    double efficiency = 0.7953 - (0.01951*batchPOT); //Based on MINERvA-doc-21436
    
    //std::cout<<"batchPOT: " << batchPOT << std::endl;
    //std::cout<<"efficiency: " << efficiency << std::endl;
    if(ANNVtx.size()==3 && cvUniv->GetANNProb()>0.2)
    {
      ANNVerticesMC->Fill( ANNVtx[2], cvWeight);
      ANNVerticesGranularMC->Fill( ANNVtx[0], ANNVtx[1], ANNVtx[2], cvWeight);
      //Water
      if (m_TargetUtils->InWaterTargetVolMC(ANNVtx[0], ANNVtx[1], ANNVtx[2]))
      {
        TrueVtxANNRecoInWater->Fill( TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), cvWeight);
        ANNRecoInWater->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      else
      {
        TrueVtxANNRecoOutWater->Fill( TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), cvWeight);
        ANNRecoOutWater->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      //Tgt1
      if (m_TargetUtils->InTarget1VolMC(ANNVtx[0], ANNVtx[1], ANNVtx[2]))
      {
        ANNRecoInTarget1->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      else
      {
        ANNRecoOutTarget1->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      //Tgt2
      if (m_TargetUtils->InTarget2VolMC(ANNVtx[0], ANNVtx[1], ANNVtx[2]))
      {
        ANNRecoInTarget2->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      else
      {
        ANNRecoOutTarget2->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      //Tgt3
      if (m_TargetUtils->InTarget3VolMC(ANNVtx[0], ANNVtx[1], ANNVtx[2]))
      {
        ANNRecoInTarget3->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      else
      {
        ANNRecoOutTarget3->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      //Tgt4
      if (m_TargetUtils->InTarget4VolMC(ANNVtx[0], ANNVtx[1], ANNVtx[2]))
      {
        ANNRecoInTarget4->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      else
      {
        ANNRecoOutTarget4->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      //Tgt5
      if (m_TargetUtils->InTarget5VolMC(ANNVtx[0], ANNVtx[1], ANNVtx[2]))
      {
        ANNRecoInTarget5->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      else
      {
        ANNRecoOutTarget5->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      //Tgt1Iron
      if (m_TargetUtils->InIron1VolMC(ANNVtx[0], ANNVtx[1], ANNVtx[2]))
      {
        ANNRecoInTarget1Iron->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      else
      {
        ANNRecoOutTarget1Iron->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      //Tgt1Lead
      if (m_TargetUtils->InLead1VolMC(ANNVtx[0], ANNVtx[1], ANNVtx[2]))
      {
        ANNRecoInTarget1Lead->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      else
      {
        ANNRecoOutTarget1Lead->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      //Tgt2Iron
      if (m_TargetUtils->InIron2VolMC(ANNVtx[0], ANNVtx[1], ANNVtx[2]))
      {
        ANNRecoInTarget2Iron->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      else
      {
        ANNRecoOutTarget2Iron->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      //Tgt2Lead
      if (m_TargetUtils->InLead2VolMC(ANNVtx[0], ANNVtx[1], ANNVtx[2]))
      {
        ANNRecoInTarget2Lead->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      else
      {
        ANNRecoOutTarget2Lead->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      //Tgt3Iron
      if (m_TargetUtils->InIron3VolMC(ANNVtx[0], ANNVtx[1], ANNVtx[2]))
      {
        ANNRecoInTarget3Iron->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      else
      {
        ANNRecoOutTarget3Iron->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      //Tgt3Lead
      if (m_TargetUtils->InLead3VolMC(ANNVtx[0], ANNVtx[1], ANNVtx[2]))
      {
        ANNRecoInTarget3Lead->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      else
      {
        ANNRecoOutTarget3Lead->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      //Tgt3Carbon
      if (m_TargetUtils->InCarbon3VolMC(ANNVtx[0], ANNVtx[1], ANNVtx[2]))
      {
        ANNRecoInTarget3Carbon->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      else
      {
        ANNRecoOutTarget3Carbon->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      //Tgt4Lead
      if (m_TargetUtils->InLead4VolMC(ANNVtx[0], ANNVtx[1], ANNVtx[2]))
      {
        ANNRecoInTarget4Lead->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      else
      {
        ANNRecoOutTarget4Lead->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      //Tgt5Iron
      if (m_TargetUtils->InIron5VolMC(ANNVtx[0], ANNVtx[1], ANNVtx[2]))
      {
        ANNRecoInTarget5Iron->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      else
      {
        ANNRecoOutTarget5Iron->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      //Tgt5Lead
      if (m_TargetUtils->InLead5VolMC(ANNVtx[0], ANNVtx[1], ANNVtx[2]))
      {
        ANNRecoInTarget5Lead->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      else
      {
        ANNRecoOutTarget5Lead->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      //Material-Carbon
      if (m_TargetUtils->InCarbonTargetVolMC(ANNVtx[0], ANNVtx[1], ANNVtx[2]))
      {
        ANNRecoInCarbon->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      else
      {
        ANNRecoOutCarbon->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      //Material-Iron
      if (m_TargetUtils->InIronTargetVolMC(ANNVtx[0], ANNVtx[1], ANNVtx[2]))
      {
        ANNRecoInIron->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      else
      {
        ANNRecoOutIron->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      //Material-Lead
      if (m_TargetUtils->InLeadTargetVolMC(ANNVtx[0], ANNVtx[1], ANNVtx[2]))
      {
        ANNRecoInLead->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }
      else
      {
        ANNRecoOutLead->Fill( TrueVtx.Z(), cvWeight/efficiency);
      }







      //Water
      if (m_TargetUtils->InWaterTargetVolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
      {
        ANNTruthInWater->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      else
      {
        ANNTruthOutWater->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      //Tgt1
      if (m_TargetUtils->InTarget1VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
      {
        ANNTruthInTarget1->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      else
      {
        ANNTruthOutTarget1->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      //Tgt2
      if (m_TargetUtils->InTarget2VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
      {
        ANNTruthInTarget2->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      else
      {
        ANNTruthOutTarget2->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      //Tgt3
      if (m_TargetUtils->InTarget3VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
      {
        ANNTruthInTarget3->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      else
      {
        ANNTruthOutTarget3->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      //Tgt4
      if (m_TargetUtils->InTarget4VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
      {
        ANNTruthInTarget4->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      else
      {
        ANNTruthOutTarget4->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      //Tgt5
      if (m_TargetUtils->InTarget5VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
      {
        ANNTruthInTarget5->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      else
      {
        ANNTruthOutTarget5->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      //Tgt1Iron
      if (m_TargetUtils->InIron1VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
      {
        ANNTruthInTarget1Iron->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      else
      {
        ANNTruthOutTarget1Iron->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      //Tgt1Lead
      if (m_TargetUtils->InLead1VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
      {
        ANNTruthInTarget1Lead->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      else
      {
        ANNTruthOutTarget1Lead->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      //Tgt2Iron
      if (m_TargetUtils->InIron2VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
      {
        ANNTruthInTarget2Iron->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      else
      {
        ANNTruthOutTarget2Iron->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      //Tgt2Lead
      if (m_TargetUtils->InLead2VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
      {
        ANNTruthInTarget2Lead->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      else
      {
        ANNTruthOutTarget2Lead->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      //Tgt3Iron
      if (m_TargetUtils->InIron3VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
      {
        ANNTruthInTarget3Iron->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      else
      {
        ANNTruthOutTarget3Iron->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      //Tgt3Lead
      if (m_TargetUtils->InLead3VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
      {
        ANNTruthInTarget3Lead->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      else
      {
        ANNTruthOutTarget3Lead->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      //Tgt3Carbon
      if (m_TargetUtils->InCarbon3VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
      {
        ANNTruthInTarget3Carbon->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      else
      {
        ANNTruthOutTarget3Carbon->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      //Tgt4Lead
      if (m_TargetUtils->InLead4VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
      {
        ANNTruthInTarget4Lead->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      else
      {
        ANNTruthOutTarget4Lead->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      //Tgt5Iron
      if (m_TargetUtils->InIron5VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
      {
        ANNTruthInTarget5Iron->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      else
      {
        ANNTruthOutTarget5Iron->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      //Tgt5Lead
      if (m_TargetUtils->InLead5VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
      {
        ANNTruthInTarget5Lead->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      else
      {
        ANNTruthOutTarget5Lead->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      //Material-Carbon
      if (m_TargetUtils->InCarbonTargetVolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
      {
        ANNTruthInCarbon->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      else
      {
        ANNTruthOutCarbon->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      //Material-Iron
      if (m_TargetUtils->InIronTargetVolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
      {
        ANNTruthInIron->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      else
      {
        ANNTruthOutIron->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      //Material-Lead
      if (m_TargetUtils->InLeadTargetVolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
      {
        ANNTruthInLead->Fill( ANNVtx[2], cvWeight/efficiency);
      }
      else
      {
        ANNTruthOutLead->Fill( ANNVtx[2], cvWeight/efficiency);
      }






    }
    TBVerticesMC->Fill(TrackBasedVtx.Z(), cvWeight);
    TBVerticesGranularMC->Fill( TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z(), cvWeight);
    TrueVerticesGranular->Fill( TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), cvWeight);
    //Water
    if (m_TargetUtils->InWaterTargetVolMC(TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z()))
    {
      TrueVtxTBRecoInWater->Fill( TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), cvWeight);
      TBRecoInWater->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      //std::cout<<"TB: X: " << TrackBasedVtx.X() <<" Y: " << TrackBasedVtx.Y() <<" Z: "<< TrackBasedVtx.Z() <<std::endl;
      TrueVtxTBRecoOutWater->Fill( TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), cvWeight);
      TBRecoOutWater->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    //Tgt1
    if (m_TargetUtils->InTarget1VolMC(TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z()))
    {
      TBRecoInTarget1->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBRecoOutTarget1->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    //Tgt2
    if (m_TargetUtils->InTarget2VolMC(TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z()))
    {
      TBRecoInTarget2->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBRecoOutTarget2->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    //Tgt3
    if (m_TargetUtils->InTarget3VolMC(TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z()))
    {
      TBRecoInTarget3->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBRecoOutTarget3->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    //Tgt4
    if (m_TargetUtils->InTarget4VolMC(TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z()))
    {
      TBRecoInTarget4->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBRecoOutTarget4->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    //Tgt5
    if (m_TargetUtils->InTarget5VolMC(TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z()))
    {
      TBRecoInTarget5->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBRecoOutTarget5->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    //Tgt1Iron
    if (m_TargetUtils->InIron1VolMC(TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z()))
    {
      TBRecoInTarget1Iron->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBRecoOutTarget1Iron->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    //Tgt1Lead
    if (m_TargetUtils->InLead1VolMC(TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z()))
    {
      TBRecoInTarget1Lead->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBRecoOutTarget1Lead->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    //Tgt2Iron
    if (m_TargetUtils->InIron2VolMC(TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z()))
    {
      TBRecoInTarget2Iron->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBRecoOutTarget2Iron->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    //Tgt2Lead
    if (m_TargetUtils->InLead2VolMC(TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z()))
    {
      TBRecoInTarget2Lead->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBRecoOutTarget2Lead->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    //Tgt3Iron
    if (m_TargetUtils->InIron3VolMC(TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z()))
    {
      TBRecoInTarget3Iron->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBRecoOutTarget3Iron->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    //Tgt3Lead
    if (m_TargetUtils->InLead3VolMC(TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z()))
    {
      TBRecoInTarget3Lead->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBRecoOutTarget3Lead->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    //Tgt3Carbon
    if (m_TargetUtils->InCarbon3VolMC(TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z()))
    {
      TBRecoInTarget3Carbon->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBRecoOutTarget3Carbon->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    //Tgt4Lead
    if (m_TargetUtils->InLead4VolMC(TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z()))
    {
      TBRecoInTarget4Lead->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBRecoOutTarget4Lead->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    //Tgt5Iron
    if (m_TargetUtils->InIron5VolMC(TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z()))
    {
      TBRecoInTarget5Iron->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBRecoOutTarget5Iron->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    //Tgt5Lead
    if (m_TargetUtils->InLead5VolMC(TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z()))
    {
      TBRecoInTarget5Lead->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBRecoOutTarget5Lead->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    //Material-Carbon
    if (m_TargetUtils->InCarbonTargetVolMC(TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z()))
    {
      TBRecoInCarbon->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBRecoOutCarbon->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    //Material-Iron
    if (m_TargetUtils->InIronTargetVolMC(TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z()))
    {
      TBRecoInIron->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBRecoOutIron->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    //Material-Lead
    if (m_TargetUtils->InLeadTargetVolMC(TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z()))
    {
      TBRecoInLead->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBRecoOutLead->Fill( TrueVtx.Z(), cvWeight/efficiency);
    }
    //std::cout<<"Here6\n";









    //Water
    if (m_TargetUtils->InWaterTargetVolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
    {
      TBTruthInWater->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBTruthOutWater->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    //Tgt1
    if (m_TargetUtils->InTarget1VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
    {
      TBTruthInTarget1->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBTruthOutTarget1->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    //Tgt2
    if (m_TargetUtils->InTarget2VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
    {
      TBTruthInTarget2->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBTruthOutTarget2->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    //Tgt3
    if (m_TargetUtils->InTarget3VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
    {
      TBTruthInTarget3->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBTruthOutTarget3->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    //Tgt4
    if (m_TargetUtils->InTarget4VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
    {
      TBTruthInTarget4->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBTruthOutTarget4->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    //Tgt5
    if (m_TargetUtils->InTarget5VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
    {
      TBTruthInTarget5->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBTruthOutTarget5->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    //Tgt1Iron
    if (m_TargetUtils->InIron1VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
    {
      TBTruthInTarget1Iron->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBTruthOutTarget1Iron->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    //Tgt1Lead
    if (m_TargetUtils->InLead1VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
    {
      TBTruthInTarget1Lead->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBTruthOutTarget1Lead->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    //Tgt2Iron
    if (m_TargetUtils->InIron2VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
    {
      TBTruthInTarget2Iron->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBTruthOutTarget2Iron->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    //Tgt2Lead
    if (m_TargetUtils->InLead2VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
    {
      TBTruthInTarget2Lead->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBTruthOutTarget2Lead->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    //Tgt3Iron
    if (m_TargetUtils->InIron3VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
    {
      TBTruthInTarget3Iron->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBTruthOutTarget3Iron->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    //Tgt3Lead
    if (m_TargetUtils->InLead3VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
    {
      TBTruthInTarget3Lead->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBTruthOutTarget3Lead->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    //Tgt3Carbon
    if (m_TargetUtils->InCarbon3VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
    {
      TBTruthInTarget3Carbon->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBTruthOutTarget3Carbon->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    //Tgt4Lead
    if (m_TargetUtils->InLead4VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
    {
      TBTruthInTarget4Lead->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBTruthOutTarget4Lead->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    //Tgt5Iron
    if (m_TargetUtils->InIron5VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
    {
      TBTruthInTarget5Iron->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBTruthOutTarget5Iron->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    //Tgt5Lead
    if (m_TargetUtils->InLead5VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
    {
      TBTruthInTarget5Lead->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBTruthOutTarget5Lead->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    //Material-Carbon
    if (m_TargetUtils->InCarbonTargetVolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
    {
      TBTruthInCarbon->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBTruthOutCarbon->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    //Material-Iron
    if (m_TargetUtils->InIronTargetVolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
    {
      TBTruthInIron->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBTruthOutIron->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    //Material-Lead
    if (m_TargetUtils->InLeadTargetVolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()))
    {
      TBTruthInLead->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    else
    {
      TBTruthOutLead->Fill( TrackBasedVtx.Z(), cvWeight/efficiency);
    }
    //std::cout<<"Here6\n";








  } //End entries loop
  std::cout << "Finished MC reco loop.\n";
}

void LoopAndFillData( PlotUtils::ChainWrapper* data,
			        std::vector<CVUniverse*> data_band,
				std::vector<Variable*> vars,
                                std::vector<Variable2D*> vars2D,
                                std::vector<Study*> studies,
				PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts)

{
  std::cout << "Starting data loop...\n";
  const int nEntries = data->GetEntries();
  for (int i=0; i<data->GetEntries(); ++i) {
    for (auto universe : data_band) {
      universe->SetEntry(i);
      if(i%1000==0) std::cout << i << " / " << nEntries << "\r" << std::endl;
      MichelEvent myevent; 
      std::vector<double> ANNVtx = universe->GetANNVertexVector();
      ROOT::Math::XYZTVector TrackBasedVtx = universe->GetVertex();

      if (!michelcuts.isDataSelected(*universe, myevent).all()) continue;
      if(ANNVtx.size()==3)
      {
        ANNVerticesData->Fill(ANNVtx[2]);
        ANNVerticesGranularData->Fill(ANNVtx[0], ANNVtx[1], ANNVtx[2]);
      }
      TBVerticesGranularData->Fill(TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z());
      TBVerticesData->Fill(TrackBasedVtx.Z());
    }
  }
  std::cout << "Finished data loop.\n";
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

//==============================================================================
// Main
//==============================================================================
int main(const int argc, const char** argv)
{
  TH1::AddDirectory(false);

  //Validate input.
  //I expect a data playlist file name and an MC playlist file name which is exactly 2 arguments.
  const int nArgsExpected = 2;
  if(argc != nArgsExpected + 1) //argc is the size of argv.  I check for number of arguments + 1 because
                                //argv[0] is always the path to the executable.
  {
    std::cerr << "Expected " << nArgsExpected << " arguments, but got " << argc - 1 << "\n" << USAGE << "\n";
    return badCmdLine;
  }

  //One playlist must contain only MC files, and the other must contain only data files.
  //Only checking the first file in each playlist because opening each file an extra time
  //remotely (e.g. through xrootd) can get expensive.
  //TODO: Look in INSTALL_DIR if files not found?
  const std::string mc_file_list = argv[2],
                    data_file_list = argv[1];

  //Check that necessary TTrees exist in the first file of mc_file_list and data_file_list
  std::string reco_tree_name;
  if(!inferRecoTreeNameAndCheckTreeNames(mc_file_list, data_file_list, reco_tree_name))
  {
    std::cerr << "Failed to find required trees in MC playlist " << mc_file_list << " and/or data playlist " << data_file_list << ".\n" << USAGE << "\n";
    return badInputFile;
  }

  const bool doCCQENuValidation = (reco_tree_name == "CCQENu"); //Enables extra histograms and might influence which systematics I use.

  //const bool is_grid = false; //TODO: Are we going to put this back?  Gonzalo needs it iirc.
  PlotUtils::MacroUtil options(reco_tree_name, mc_file_list, data_file_list, "minervame1A", true);
  options.m_plist_string = util::GetPlaylist(*options.m_mc, true); //TODO: Put GetPlaylist into PlotUtils::MacroUtil

  // You're required to make some decisions
  PlotUtils::MinervaUniverse::SetNuEConstraint(true);
  PlotUtils::MinervaUniverse::SetPlaylist(options.m_plist_string); //TODO: Infer this from the files somehow?
  PlotUtils::MinervaUniverse::SetAnalysisNuPDG(14);
  PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
  PlotUtils::MinervaUniverse::SetZExpansionFaReweight(false);

  PlotUtils::MinervaUniverse::RPAMaterials(true); 

  //Now that we've defined what a cross section is, decide which sample and model
  //we're extracting a cross section for.
  PlotUtils::Cutter<CVUniverse, MichelEvent>::reco_t sidebands, preCuts;
  PlotUtils::Cutter<CVUniverse, MichelEvent>::truth_t signalDefinition, phaseSpace;

  //const double minZ = 5980, maxZ = 8422, apothem = 850; //All in mm
  const double apothem = 850; //All in mm
  preCuts.emplace_back(new reco::ZRange<CVUniverse, MichelEvent>("Active Tracker Z pos", PlotUtils::TargetProp::NukeRegion::Face, PlotUtils::TargetProp::Tracker::Back));
  preCuts.emplace_back(new reco::Apothem<CVUniverse, MichelEvent>(apothem));
  preCuts.emplace_back(new reco::MaxMuonAngle<CVUniverse, MichelEvent>(17.));
  preCuts.emplace_back(new reco::HasMINOSMatch<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::NoDeadtime<CVUniverse, MichelEvent>(1, "Deadtime"));
  preCuts.emplace_back(new reco::IsNeutrino<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::MuonEnergyMin<CVUniverse, MichelEvent>(2000.0, "EMu Min"));
  preCuts.emplace_back(new reco::MuonEnergyMax<CVUniverse, MichelEvent>(50000.0, "EMu Max"));
  preCuts.emplace_back(new reco::ANNConfidenceCut<CVUniverse, MichelEvent>(0.20));


                                                                                                                                                   
  signalDefinition.emplace_back(new truth::IsNeutrino<CVUniverse>());
  signalDefinition.emplace_back(new truth::IsCC<CVUniverse>());
                                                                                                                                                   
  phaseSpace.emplace_back(new truth::ZRange<CVUniverse>("Active Tracker Z pos", PlotUtils::TargetProp::NukeRegion::Face, PlotUtils::TargetProp::Tracker::Back));
  phaseSpace.emplace_back(new truth::Apothem<CVUniverse>(apothem));
  phaseSpace.emplace_back(new truth::MuonAngle<CVUniverse>(17.));
  phaseSpace.emplace_back(new truth::MuonEnergyMin<CVUniverse>(2000.0, "EMu Min"));
  phaseSpace.emplace_back(new truth::MuonEnergyMax<CVUniverse>(50000.0, "EMu Max"));
  //phaseSpace.emplace_back(new truth::PZMuMin<CVUniverse>(1500.));
                                                                                                                                                   
  PlotUtils::Cutter<CVUniverse, MichelEvent> mycuts(std::move(preCuts), std::move(sidebands) , std::move(signalDefinition),std::move(phaseSpace));

  std::vector<std::unique_ptr<PlotUtils::Reweighter<CVUniverse, MichelEvent>>> MnvTunev1;
  MnvTunev1.emplace_back(new PlotUtils::FluxAndCVReweighter<CVUniverse, MichelEvent>());
  MnvTunev1.emplace_back(new PlotUtils::GENIEReweighter<CVUniverse, MichelEvent>(true, false));
  MnvTunev1.emplace_back(new PlotUtils::LowRecoil2p2hReweighter<CVUniverse, MichelEvent>());
  MnvTunev1.emplace_back(new PlotUtils::MINOSEfficiencyReweighter<CVUniverse, MichelEvent>());
  MnvTunev1.emplace_back(new PlotUtils::RPAReweighter<CVUniverse, MichelEvent>());

  PlotUtils::Model<CVUniverse, MichelEvent> model(std::move(MnvTunev1));

  // Make a map of systematic universes
  // Leave out systematics when making validation histograms
  const bool doSystematics = (getenv("MNV101_SKIP_SYST") == nullptr);
  if(!doSystematics){
    std::cout << "Skipping systematics (except 1 flux universe) because environment variable MNV101_SKIP_SYST is set.\n";
    PlotUtils::MinervaUniverse::SetNFluxUniverses(2); //Necessary to get Flux integral later...  Doesn't work with just 1 flux universe though because _that_ triggers "spread errors".
  }

  std::map< std::string, std::vector<CVUniverse*> > error_bands;
  if(doSystematics) error_bands = GetStandardSystematics(options.m_mc);
  else{
    std::map<std::string, std::vector<CVUniverse*> > band_flux = PlotUtils::GetFluxSystematicsMap<CVUniverse>(options.m_mc, CVUniverse::GetNFluxUniverses());
    error_bands.insert(band_flux.begin(), band_flux.end()); //Necessary to get flux integral later...
  }
  error_bands["cv"] = {new CVUniverse(options.m_mc)};
  std::map< std::string, std::vector<CVUniverse*> > truth_bands;
  if(doSystematics) truth_bands = GetStandardSystematics(options.m_truth);
  truth_bands["cv"] = {new CVUniverse(options.m_truth)};

  std::vector<double> dansPTBins = {0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1, 1.25, 1.5, 2.5, 4.5},
                      dansPzBins = {1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20, 40, 60},
                      robsEmuBins = {0,1,2,3,4,5,7,9,12,15,18,22,36,50,75,100,120},
                      bjorkenXbins = {0.0, 0.1, 0.3, 0.5, 0.7, 0.9 , 1.1, 1.5},
                      robsRecoilBins;

  const double robsRecoilBinWidth = 50; //MeV
  for(int whichBin = 0; whichBin < 100 + 1; ++whichBin) robsRecoilBins.push_back(robsRecoilBinWidth * whichBin);

  std::vector<Variable*> vars;

  std::vector<Variable2D*> vars2D;

  std::vector<Study*> studies;

  CVUniverse* data_universe = new CVUniverse(options.m_data);
  std::vector<CVUniverse*> data_band = {data_universe};
  std::map<std::string, std::vector<CVUniverse*> > data_error_bands;
  data_error_bands["cv"] = data_band;
  
  std::vector<Study*> data_studies;

  for(auto& var: vars) var->InitializeMCHists(error_bands, truth_bands);
  for(auto& var: vars) var->InitializeDATAHists(data_band);

  for(auto& var: vars2D) var->InitializeMCHists(error_bands, truth_bands);
  for(auto& var: vars2D) var->InitializeDATAHists(data_band);

  // Loop entries and fill
  try
  {
    CVUniverse::SetTruth(false);
    LoopAndFillEventSelection(options.m_mc, error_bands, vars, vars2D, studies, mycuts, model);
    CVUniverse::SetTruth(true);
    //LoopAndFillEffDenom(options.m_truth, truth_bands, vars, vars2D, mycuts, model);
    options.PrintMacroConfiguration(argv[0]);
    std::cout << "MC cut summary:\n" << mycuts << "\n";
    mycuts.resetStats();

    CVUniverse::SetTruth(false);
    LoopAndFillData(options.m_data, data_band, vars, vars2D, data_studies, mycuts);
    std::cout << "Data cut summary:\n" << mycuts << "\n";
    
    //Write MC results
    TFile* mcOutDir = TFile::Open(MC_OUT_FILE_NAME, "RECREATE");
    if(!mcOutDir)
    {
      std::cerr << "Failed to open a file named " << MC_OUT_FILE_NAME << " in the current directory for writing histograms.\n";
      return badOutputFile;
    }

    for(auto& study: studies) study->SaveOrDraw(*mcOutDir);
    for(auto& var: vars) var->WriteMC(*mcOutDir);
    for(auto& var: vars2D) var->WriteMC(*mcOutDir);

    //Protons On Target
    auto mcPOT = new TParameter<double>("POTUsed", options.m_mc_pot);
    mcPOT->Write();


    ANNVerticesMC->SetDirectory(mcOutDir);
    TBVerticesMC->SetDirectory(mcOutDir);
    ANNVerticesMC->Write();
    TBVerticesMC->Write();

    ANNVerticesGranularMC->SetDirectory(mcOutDir);
    TBVerticesGranularMC->SetDirectory(mcOutDir);
    ANNVerticesGranularMC->Write();
    TBVerticesGranularMC->Write();

    TrueVerticesGranular->SetDirectory(mcOutDir);
    TrueVerticesGranular->Write();
    TrueVtxANNRecoInWater->SetDirectory(mcOutDir);
    TrueVtxANNRecoInWater->Write();
    TrueVtxTBRecoInWater->SetDirectory(mcOutDir);
    TrueVtxTBRecoInWater->Write();
    TrueVtxANNRecoOutWater->SetDirectory(mcOutDir);
    TrueVtxANNRecoOutWater->Write();
    TrueVtxTBRecoOutWater->SetDirectory(mcOutDir);
    TrueVtxTBRecoOutWater->Write();

    ANNRecoInWater->SetDirectory(mcOutDir);
    ANNRecoInWater->Write();
    TBRecoInWater->SetDirectory(mcOutDir);
    TBRecoInWater->Write();
    ANNRecoOutWater->SetDirectory(mcOutDir);
    ANNRecoOutWater->Write();
    TBRecoOutWater->SetDirectory(mcOutDir);
    TBRecoOutWater->Write();

    ANNRecoInTarget1->SetDirectory(mcOutDir);
    ANNRecoInTarget1->Write();
    TBRecoInTarget1->SetDirectory(mcOutDir);
    TBRecoInTarget1->Write();
    ANNRecoOutTarget1->SetDirectory(mcOutDir);
    ANNRecoOutTarget1->Write();
    TBRecoOutTarget1->SetDirectory(mcOutDir);
    TBRecoOutTarget1->Write();

    ANNRecoInTarget2->SetDirectory(mcOutDir);
    ANNRecoInTarget2->Write();
    TBRecoInTarget2->SetDirectory(mcOutDir);
    TBRecoInTarget2->Write();
    ANNRecoOutTarget2->SetDirectory(mcOutDir);
    ANNRecoOutTarget2->Write();
    TBRecoOutTarget2->SetDirectory(mcOutDir);
    TBRecoOutTarget2->Write();

    ANNRecoInTarget3->SetDirectory(mcOutDir);
    ANNRecoInTarget3->Write();
    TBRecoInTarget3->SetDirectory(mcOutDir);
    TBRecoInTarget3->Write();
    ANNRecoOutTarget3->SetDirectory(mcOutDir);
    ANNRecoOutTarget3->Write();
    TBRecoOutTarget3->SetDirectory(mcOutDir);
    TBRecoOutTarget3->Write();

    ANNRecoInTarget4->SetDirectory(mcOutDir);
    ANNRecoInTarget4->Write();
    TBRecoInTarget4->SetDirectory(mcOutDir);
    TBRecoInTarget4->Write();
    ANNRecoOutTarget4->SetDirectory(mcOutDir);
    ANNRecoOutTarget4->Write();
    TBRecoOutTarget4->SetDirectory(mcOutDir);
    TBRecoOutTarget4->Write();

    ANNRecoInTarget5->SetDirectory(mcOutDir);
    ANNRecoInTarget5->Write();
    TBRecoInTarget5->SetDirectory(mcOutDir);
    TBRecoInTarget5->Write();
    ANNRecoOutTarget5->SetDirectory(mcOutDir);
    ANNRecoOutTarget5->Write();
    TBRecoOutTarget5->SetDirectory(mcOutDir);
    TBRecoOutTarget5->Write();

    ANNRecoInTarget1Iron->SetDirectory(mcOutDir);
    ANNRecoInTarget1Iron->Write();
    TBRecoInTarget1Iron->SetDirectory(mcOutDir);
    TBRecoInTarget1Iron->Write();
    ANNRecoOutTarget1Iron->SetDirectory(mcOutDir);
    ANNRecoOutTarget1Iron->Write();
    TBRecoOutTarget1Iron->SetDirectory(mcOutDir);
    TBRecoOutTarget1Iron->Write();

    ANNRecoInTarget1Lead->SetDirectory(mcOutDir);
    ANNRecoInTarget1Lead->Write();
    TBRecoInTarget1Lead->SetDirectory(mcOutDir);
    TBRecoInTarget1Lead->Write();
    ANNRecoOutTarget1Lead->SetDirectory(mcOutDir);
    ANNRecoOutTarget1Lead->Write();
    TBRecoOutTarget1Lead->SetDirectory(mcOutDir);
    TBRecoOutTarget1Lead->Write();

    ANNRecoInTarget2Iron->SetDirectory(mcOutDir);
    ANNRecoInTarget2Iron->Write();
    TBRecoInTarget2Iron->SetDirectory(mcOutDir);
    TBRecoInTarget2Iron->Write();
    ANNRecoOutTarget2Iron->SetDirectory(mcOutDir);
    ANNRecoOutTarget2Iron->Write();
    TBRecoOutTarget2Iron->SetDirectory(mcOutDir);
    TBRecoOutTarget2Iron->Write();

    ANNRecoInTarget2Lead->SetDirectory(mcOutDir);
    ANNRecoInTarget2Lead->Write();
    TBRecoInTarget2Lead->SetDirectory(mcOutDir);
    TBRecoInTarget2Lead->Write();
    ANNRecoOutTarget2Lead->SetDirectory(mcOutDir);
    ANNRecoOutTarget2Lead->Write();
    TBRecoOutTarget2Lead->SetDirectory(mcOutDir);
    TBRecoOutTarget2Lead->Write();

    ANNRecoInTarget3Iron->SetDirectory(mcOutDir);
    ANNRecoInTarget3Iron->Write();
    TBRecoInTarget3Iron->SetDirectory(mcOutDir);
    TBRecoInTarget3Iron->Write();
    ANNRecoOutTarget3Iron->SetDirectory(mcOutDir);
    ANNRecoOutTarget3Iron->Write();
    TBRecoOutTarget3Iron->SetDirectory(mcOutDir);
    TBRecoOutTarget3Iron->Write();

    ANNRecoInTarget3Lead->SetDirectory(mcOutDir);
    ANNRecoInTarget3Lead->Write();
    TBRecoInTarget3Lead->SetDirectory(mcOutDir);
    TBRecoInTarget3Lead->Write();
    ANNRecoOutTarget3Lead->SetDirectory(mcOutDir);
    ANNRecoOutTarget3Lead->Write();
    TBRecoOutTarget3Lead->SetDirectory(mcOutDir);
    TBRecoOutTarget3Lead->Write();

    ANNRecoInTarget3Carbon->SetDirectory(mcOutDir);
    ANNRecoInTarget3Carbon->Write();
    TBRecoInTarget3Carbon->SetDirectory(mcOutDir);
    TBRecoInTarget3Carbon->Write();
    ANNRecoOutTarget3Carbon->SetDirectory(mcOutDir);
    ANNRecoOutTarget3Carbon->Write();
    TBRecoOutTarget3Carbon->SetDirectory(mcOutDir);
    TBRecoOutTarget3Carbon->Write();

    ANNRecoInTarget4Lead->SetDirectory(mcOutDir);
    ANNRecoInTarget4Lead->Write();
    TBRecoInTarget4Lead->SetDirectory(mcOutDir);
    TBRecoInTarget4Lead->Write();
    ANNRecoOutTarget4Lead->SetDirectory(mcOutDir);
    ANNRecoOutTarget4Lead->Write();
    TBRecoOutTarget4Lead->SetDirectory(mcOutDir);
    TBRecoOutTarget4Lead->Write();

    ANNRecoInTarget5Iron->SetDirectory(mcOutDir);
    ANNRecoInTarget5Iron->Write();
    TBRecoInTarget5Iron->SetDirectory(mcOutDir);
    TBRecoInTarget5Iron->Write();
    ANNRecoOutTarget5Iron->SetDirectory(mcOutDir);
    ANNRecoOutTarget5Iron->Write();
    TBRecoOutTarget5Iron->SetDirectory(mcOutDir);
    TBRecoOutTarget5Iron->Write();

    ANNRecoInTarget5Lead->SetDirectory(mcOutDir);
    ANNRecoInTarget5Lead->Write();
    TBRecoInTarget5Lead->SetDirectory(mcOutDir);
    TBRecoInTarget5Lead->Write();
    ANNRecoOutTarget5Lead->SetDirectory(mcOutDir);
    ANNRecoOutTarget5Lead->Write();
    TBRecoOutTarget5Lead->SetDirectory(mcOutDir);
    TBRecoOutTarget5Lead->Write();

    ANNRecoInIron->SetDirectory(mcOutDir);
    ANNRecoInIron->Write();
    TBRecoInIron->SetDirectory(mcOutDir);
    TBRecoInIron->Write();
    ANNRecoOutIron->SetDirectory(mcOutDir);
    ANNRecoOutIron->Write();
    TBRecoOutIron->SetDirectory(mcOutDir);
    TBRecoOutIron->Write();
  
    ANNRecoInLead->SetDirectory(mcOutDir);
    ANNRecoInLead->Write();
    TBRecoInLead->SetDirectory(mcOutDir);
    TBRecoInLead->Write();
    ANNRecoOutLead->SetDirectory(mcOutDir);
    ANNRecoOutLead->Write();
    TBRecoOutLead->SetDirectory(mcOutDir);
    TBRecoOutLead->Write();

    ANNRecoInCarbon->SetDirectory(mcOutDir);
    ANNRecoInCarbon->Write();
    TBRecoInCarbon->SetDirectory(mcOutDir);
    TBRecoInCarbon->Write();
    ANNRecoOutCarbon->SetDirectory(mcOutDir);
    ANNRecoOutCarbon->Write();
    TBRecoOutCarbon->SetDirectory(mcOutDir);
    TBRecoOutCarbon->Write();










    ANNTruthInWater->SetDirectory(mcOutDir);
    ANNTruthInWater->Write();
    TBTruthInWater->SetDirectory(mcOutDir);
    TBTruthInWater->Write();
    ANNTruthOutWater->SetDirectory(mcOutDir);
    ANNTruthOutWater->Write();
    TBTruthOutWater->SetDirectory(mcOutDir);
    TBTruthOutWater->Write();

    ANNTruthInTarget1->SetDirectory(mcOutDir);
    ANNTruthInTarget1->Write();
    TBTruthInTarget1->SetDirectory(mcOutDir);
    TBTruthInTarget1->Write();
    ANNTruthOutTarget1->SetDirectory(mcOutDir);
    ANNTruthOutTarget1->Write();
    TBTruthOutTarget1->SetDirectory(mcOutDir);
    TBTruthOutTarget1->Write();

    ANNTruthInTarget2->SetDirectory(mcOutDir);
    ANNTruthInTarget2->Write();
    TBTruthInTarget2->SetDirectory(mcOutDir);
    TBTruthInTarget2->Write();
    ANNTruthOutTarget2->SetDirectory(mcOutDir);
    ANNTruthOutTarget2->Write();
    TBTruthOutTarget2->SetDirectory(mcOutDir);
    TBTruthOutTarget2->Write();

    ANNTruthInTarget3->SetDirectory(mcOutDir);
    ANNTruthInTarget3->Write();
    TBTruthInTarget3->SetDirectory(mcOutDir);
    TBTruthInTarget3->Write();
    ANNTruthOutTarget3->SetDirectory(mcOutDir);
    ANNTruthOutTarget3->Write();
    TBTruthOutTarget3->SetDirectory(mcOutDir);
    TBTruthOutTarget3->Write();

    ANNTruthInTarget4->SetDirectory(mcOutDir);
    ANNTruthInTarget4->Write();
    TBTruthInTarget4->SetDirectory(mcOutDir);
    TBTruthInTarget4->Write();
    ANNTruthOutTarget4->SetDirectory(mcOutDir);
    ANNTruthOutTarget4->Write();
    TBTruthOutTarget4->SetDirectory(mcOutDir);
    TBTruthOutTarget4->Write();

    ANNTruthInTarget5->SetDirectory(mcOutDir);
    ANNTruthInTarget5->Write();
    TBTruthInTarget5->SetDirectory(mcOutDir);
    TBTruthInTarget5->Write();
    ANNTruthOutTarget5->SetDirectory(mcOutDir);
    ANNTruthOutTarget5->Write();
    TBTruthOutTarget5->SetDirectory(mcOutDir);
    TBTruthOutTarget5->Write();

    ANNTruthInTarget1Iron->SetDirectory(mcOutDir);
    ANNTruthInTarget1Iron->Write();
    TBTruthInTarget1Iron->SetDirectory(mcOutDir);
    TBTruthInTarget1Iron->Write();
    ANNTruthOutTarget1Iron->SetDirectory(mcOutDir);
    ANNTruthOutTarget1Iron->Write();
    TBTruthOutTarget1Iron->SetDirectory(mcOutDir);
    TBTruthOutTarget1Iron->Write();

    ANNTruthInTarget1Lead->SetDirectory(mcOutDir);
    ANNTruthInTarget1Lead->Write();
    TBTruthInTarget1Lead->SetDirectory(mcOutDir);
    TBTruthInTarget1Lead->Write();
    ANNTruthOutTarget1Lead->SetDirectory(mcOutDir);
    ANNTruthOutTarget1Lead->Write();
    TBTruthOutTarget1Lead->SetDirectory(mcOutDir);
    TBTruthOutTarget1Lead->Write();

    ANNTruthInTarget2Iron->SetDirectory(mcOutDir);
    ANNTruthInTarget2Iron->Write();
    TBTruthInTarget2Iron->SetDirectory(mcOutDir);
    TBTruthInTarget2Iron->Write();
    ANNTruthOutTarget2Iron->SetDirectory(mcOutDir);
    ANNTruthOutTarget2Iron->Write();
    TBTruthOutTarget2Iron->SetDirectory(mcOutDir);
    TBTruthOutTarget2Iron->Write();

    ANNTruthInTarget2Lead->SetDirectory(mcOutDir);
    ANNTruthInTarget2Lead->Write();
    TBTruthInTarget2Lead->SetDirectory(mcOutDir);
    TBTruthInTarget2Lead->Write();
    ANNTruthOutTarget2Lead->SetDirectory(mcOutDir);
    ANNTruthOutTarget2Lead->Write();
    TBTruthOutTarget2Lead->SetDirectory(mcOutDir);
    TBTruthOutTarget2Lead->Write();

    ANNTruthInTarget3Iron->SetDirectory(mcOutDir);
    ANNTruthInTarget3Iron->Write();
    TBTruthInTarget3Iron->SetDirectory(mcOutDir);
    TBTruthInTarget3Iron->Write();
    ANNTruthOutTarget3Iron->SetDirectory(mcOutDir);
    ANNTruthOutTarget3Iron->Write();
    TBTruthOutTarget3Iron->SetDirectory(mcOutDir);
    TBTruthOutTarget3Iron->Write();

    ANNTruthInTarget3Lead->SetDirectory(mcOutDir);
    ANNTruthInTarget3Lead->Write();
    TBTruthInTarget3Lead->SetDirectory(mcOutDir);
    TBTruthInTarget3Lead->Write();
    ANNTruthOutTarget3Lead->SetDirectory(mcOutDir);
    ANNTruthOutTarget3Lead->Write();
    TBTruthOutTarget3Lead->SetDirectory(mcOutDir);
    TBTruthOutTarget3Lead->Write();

    ANNTruthInTarget3Carbon->SetDirectory(mcOutDir);
    ANNTruthInTarget3Carbon->Write();
    TBTruthInTarget3Carbon->SetDirectory(mcOutDir);
    TBTruthInTarget3Carbon->Write();
    ANNTruthOutTarget3Carbon->SetDirectory(mcOutDir);
    ANNTruthOutTarget3Carbon->Write();
    TBTruthOutTarget3Carbon->SetDirectory(mcOutDir);
    TBTruthOutTarget3Carbon->Write();

    ANNTruthInTarget4Lead->SetDirectory(mcOutDir);
    ANNTruthInTarget4Lead->Write();
    TBTruthInTarget4Lead->SetDirectory(mcOutDir);
    TBTruthInTarget4Lead->Write();
    ANNTruthOutTarget4Lead->SetDirectory(mcOutDir);
    ANNTruthOutTarget4Lead->Write();
    TBTruthOutTarget4Lead->SetDirectory(mcOutDir);
    TBTruthOutTarget4Lead->Write();

    ANNTruthInTarget5Iron->SetDirectory(mcOutDir);
    ANNTruthInTarget5Iron->Write();
    TBTruthInTarget5Iron->SetDirectory(mcOutDir);
    TBTruthInTarget5Iron->Write();
    ANNTruthOutTarget5Iron->SetDirectory(mcOutDir);
    ANNTruthOutTarget5Iron->Write();
    TBTruthOutTarget5Iron->SetDirectory(mcOutDir);
    TBTruthOutTarget5Iron->Write();

    ANNTruthInTarget5Lead->SetDirectory(mcOutDir);
    ANNTruthInTarget5Lead->Write();
    TBTruthInTarget5Lead->SetDirectory(mcOutDir);
    TBTruthInTarget5Lead->Write();
    ANNTruthOutTarget5Lead->SetDirectory(mcOutDir);
    ANNTruthOutTarget5Lead->Write();
    TBTruthOutTarget5Lead->SetDirectory(mcOutDir);
    TBTruthOutTarget5Lead->Write();

    ANNTruthInIron->SetDirectory(mcOutDir);
    ANNTruthInIron->Write();
    TBTruthInIron->SetDirectory(mcOutDir);
    TBTruthInIron->Write();
    ANNTruthOutIron->SetDirectory(mcOutDir);
    ANNTruthOutIron->Write();
    TBTruthOutIron->SetDirectory(mcOutDir);
    TBTruthOutIron->Write();
  
    ANNTruthInLead->SetDirectory(mcOutDir);
    ANNTruthInLead->Write();
    TBTruthInLead->SetDirectory(mcOutDir);
    TBTruthInLead->Write();
    ANNTruthOutLead->SetDirectory(mcOutDir);
    ANNTruthOutLead->Write();
    TBTruthOutLead->SetDirectory(mcOutDir);
    TBTruthOutLead->Write();

    ANNTruthInCarbon->SetDirectory(mcOutDir);
    ANNTruthInCarbon->Write();
    TBTruthInCarbon->SetDirectory(mcOutDir);
    TBTruthInCarbon->Write();
    ANNTruthOutCarbon->SetDirectory(mcOutDir);
    ANNTruthOutCarbon->Write();
    TBTruthOutCarbon->SetDirectory(mcOutDir);
    TBTruthOutCarbon->Write();


    //Write data results
    TFile* dataOutDir = TFile::Open(DATA_OUT_FILE_NAME, "RECREATE");
    if(!dataOutDir)
    {
      std::cerr << "Failed to open a file named " << DATA_OUT_FILE_NAME << " in the current directory for writing histograms.\n";
      return badOutputFile;
    }

    for(auto& var: vars) var->WriteData(*dataOutDir);
    for(auto& var: vars2D) var->WriteData(*dataOutDir);

    //Protons On Target
    auto dataPOT = new TParameter<double>("POTUsed", options.m_data_pot);
    dataPOT->Write();
    
    ANNVerticesData->SetDirectory(dataOutDir);
    TBVerticesData->SetDirectory(dataOutDir);
    ANNVerticesData->Write();
    TBVerticesData->Write();

    ANNVerticesGranularData->SetDirectory(dataOutDir);
    TBVerticesGranularData->SetDirectory(dataOutDir);
    ANNVerticesGranularData->Write();
    TBVerticesGranularData->Write();

    std::cout << "Success" << std::endl;
  }
  catch(const ROOT::exception& e)
  {
    std::cerr << "Ending on a ROOT error message.  No histograms will be produced.\n"
              << "If the message talks about \"TNetXNGFile\", this could be a problem with dCache.  The message is:\n"
              << e.what() << "\n" << USAGE << "\n";
    return badFileRead;
  }

  return success;
}