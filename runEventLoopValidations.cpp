/* 
Akeem Hart 2025

To whoever stumbles across this in the future:

This was intially intended to be a one-file-contains-all special event loop just to
produce validation plots or, more specifically, plots that aren't directly utilised
in extracting a cross section but are important to check or investigate something.
This was initially derived from runEventLoop.cpp but has since significantly changed

It has the same (basic) cuts and signal definition as the cross section event loops but
this allows me to keep those event loops, which I run more often, less bloated so they
run faster, are easier to maintain/modify and have a lower memory footprint when running,
so I can make lower grid memory requests and get higher priority.

Currently main function of this macro is to produce z vertex plots initially intended
to determine how our vertexing fares with different cuts or models in different 
detector regions and between different vertexing approaches (ML vs Track Based). Now
also includes an ERecoil segmentation into 4 quartiles and 

Consequently, because this event loop is where I just throw *everything* that isn't
fed into the cross section pipeline, it's possible that it may in the future become
bloated and require fragmentation. Note to future readers as of 9 March 2026 -  It 
looks like that previous note has come to pass and this macro has become huge and
bloated. Feel free to fragment it.


Build with:
g++ -g runEventLoopValidations.cpp util/**.cpp -o runEventLoopValidations -I $CCQENU_ANA/ -I$CCQENU_INCLUDE_PATH -ggdb -Wall -fPIC -I$PLOTUTILSROOT/../include/ -I/exp/minerva/app/users/alhart/MAT_AL9/opt/include -I. `root-config --cflags --glibs --libs` -L$PLOTUTILSROOT -lMAT -lMAT-MINERvA -lUnfoldUtils -L.

Or run with root -b -l runEventLoopValidations.cpp++g
e.g root -b -l "runEventLoopValidations.cpp++g('/exp/minerva/app/users/alhart/MAT_AL9/MINERvA-101-Cross-Section/PlaylistFiles/me-playlists/DataP6/Test1A-Data.txt', '/exp/minerva/app/users/alhart/MAT_AL9/MINERvA-101-Cross-Section/PlaylistFiles/me-playlists/MCP6/Test1A-MC.txt')"
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

#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvH3D.h"

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

#include "utilities/SMSCalibrator.h"

#include <csignal>

//These 2 variables are used when doing broken down runs to speed up running on the grid
int nSubruns = 0;
int nProcess = 0;
bool verbose = false;

//std::vector<double> vertexBins = {4293.04, 4337.25, 4381.47, 4425.68, 4514.11, 4558.33, 4602.54, 4646.76, 4735.19, 4779.4, 4823.62, 4867.83, 5000.48, 5044.69, 5088.91, 5133.12, 5456.74, 5500.95, 5545.17, 5589.38, 5677.81, 5722.03, 5810.45, 5855.68, 5900.91, 5946.14, 5991.37, 6036.6, 6081.83, 6127.06, 6172.29, 6217.52, 6262.74, 6307.97, 6353.2, 6398.43, 6443.66, 6488.89, 6534.12, 6579.35, 6624.58, 6669.81, 6715.03, 6760.26, 6805.49, 6850.72, 6895.95, 6941.18, 6986.41, 7031.64, 7076.87, 7122.1, 7167.32, 7212.55, 7257.78, 7303.01, 7348.24, 7393.47, 7438.7, 7483.93, 7529.16, 7574.39, 7619.61, 7664.84, 7710.07, 7755.3, 7800.53, 7845.76, 7890.99, 7936.22, 7981.45, 8026.68, 8071.9, 8117.13, 8162.36, 8207.59, 8252.82, 8298.05, 8343.28, 8388.51, 8433.74, 8478.97, 8524.19, 8569.42, 8614.65};

// /std::vector<double> vertexBins = {4293.04, 4313.68, 4337.25, 4357.9, 4381.47, 4402.11, 4425.68, 4446.33, 4514.11, 4534.76, 4558.33, 4578.97, 4602.54, 4623.19, 4646.76, 4667.4, 4735.19, 4755.83, 4779.4, 4800.05, 4823.62, 4844.26, 4867.83, 4888.48, 5000.48, 5021.12, 5044.69, 5065.34, 5088.91, 5109.55, 5133.12, 5153.77, 5456.74, 5477.38, 5500.95, 5521.6, 5545.17, 5565.81, 5589.38, 5610.02, 5677.81, 5698.45, 5722.03, 5742.67, 5810.45, 5831.1, 5855.68, 5876.33, 5900.91, 5921.56, 5946.14, 5966.79, 5991.37, 6012.01, 6036.6, 6057.24, 6081.83, 6102.47, 6127.06, 6147.7, 6172.29, 6192.93, 6217.52, 6238.16, 6262.74, 6283.39, 6307.97, 6328.62, 6353.2, 6373.85, 6398.43, 6419.08, 6443.66, 6464.3, 6488.89, 6509.53, 6534.12, 6554.76, 6579.35, 6599.99, 6624.58, 6645.22, 6669.81, 6690.45, 6715.03, 6735.68, 6760.26, 6780.91, 6805.49, 6826.14, 6850.72, 6871.37, 6895.95, 6916.59, 6941.18, 6961.82, 6986.41, 7007.05, 7031.64, 7052.28, 7076.87, 7097.51, 7122.1, 7142.74, 7167.32, 7187.97, 7212.55, 7233.2, 7257.78, 7278.43, 7303.01, 7323.66, 7348.24, 7368.88, 7393.47, 7414.11, 7438.7, 7459.34, 7483.93, 7504.57, 7529.16, 7549.8, 7574.39, 7595.03, 7619.61, 7640.26, 7664.84, 7685.49, 7710.07, 7730.72, 7755.3, 7775.95, 7800.53, 7821.17, 7845.76, 7866.4, 7890.99, 7911.63, 7936.22, 7956.86, 7981.45, 8002.09, 8026.68, 8047.32, 8071.9, 8092.55, 8117.13, 8137.78, 8162.36, 8183.01, 8207.59, 8228.24, 8252.82, 8273.46, 8298.05, 8318.69, 8343.28, 8363.92, 8388.51, 8409.15, 8433.74, 8454.38, 8478.97, 8499.61, 8524.19, 8544.84, 8569.42, 8590.07, 8614.65, 8635.3, 8659.46, 8680.1, 8704.26, 8724.9, 8749.06, 8769.71, 8793.86, 8814.51, 8838.67, 8859.31, 8883.47, 8904.11, 8928.27, 8948.92, 8973.08, 8993.72, 9017.88, 9038.52, 9088.08, 9135.41, 9182.75, 9230.08, 9277.41, 9324.74, 9372.08, 9419.41, 9466.74, 9514.07, 9561.41, 9608.74, 9656.07, 9703.4, 9750.74, 9798.07, 9845.4, 9892.73, 9940.07, 9987.4};

std::vector<double> vertexBins = {4282.72, 4303.36, 4325.47, 4347.57, 4369.68, 4391.79, 4413.9, 4436.01, 4456.655, 4503.79, 4524.43, 4546.55, 4568.65, 4590.76, 4612.86, 4634.98, 4657.08, 4677.72, 4724.87, 4745.51, 4767.61, 4789.73, 4811.84, 4833.94, 4856.05, 4878.15, 4898.8, 4990.16, 5010.8, 5032.9, 5055.01, 5077.12, 5099.23, 5121.34, 5143.44, 5200, 5420, 5467.06, 5489.16, 5511.27, 5533.39, 5555.49, 5577.6, 5599.7, 5620.34, 5667.49, 5688.13, 5710.24, 5732.35, 5752.99, 5800.13, 5820.77, 5843.39, 5866.01, 5888.62, 5911.24, 5933.85, 5956.47, 5979.08, 6001.69, 6024.31, 6046.92, 6069.53, 6092.15, 6114.77, 6137.38, 6159.99, 6182.61, 6205.23, 6227.84, 6250.45, 6273.07, 6295.68, 6318.3, 6340.91, 6363.52, 6386.14, 6408.76, 6431.37, 6453.98, 6476.6, 6499.21, 6521.82, 6544.44, 6567.06, 6589.67, 6612.28, 6634.9, 6657.52, 6680.13, 6702.74, 6725.35, 6747.97, 6770.59, 6793.2, 6815.82, 6838.43, 6861.05, 6883.66, 6906.27, 6928.89, 6951.5, 6974.11, 6996.73, 7019.35, 7041.96, 7064.57, 7087.19, 7109.81, 7132.42, 7155.03, 7177.65, 7200.26, 7222.88, 7245.49, 7268.1, 7290.72, 7313.34, 7335.95, 7358.56, 7381.18, 7403.79, 7426.4, 7449.02, 7471.64, 7494.25, 7516.86, 7539.48, 7562.1, 7584.71, 7607.32, 7629.93, 7652.55, 7675.16, 7697.78, 7720.4, 7743.01, 7765.62, 7788.24, 7810.85, 7833.47, 7856.08, 7878.69, 7901.31, 7923.93, 7946.54, 7969.15, 7991.77, 8014.39, 8037.0, 8059.61, 8082.23, 8104.84, 8127.45, 8150.07, 8172.68, 8195.3, 8217.92, 8240.53, 8263.14, 8285.75, 8308.37, 8330.99, 8353.6, 8376.22, 8398.83, 8421.44, 8444.06, 8466.67, 8489.29, 8511.9, 8534.51, 8557.13, 8579.74, 8602.36, 8624.97, 8647.38, 8669.78, 8692.18, 8714.58, 8736.98, 8759.38, 8781.78, 8804.19, 8826.59, 8848.99, 8871.39, 8893.79, 8916.19, 8938.6, 8961.0, 8983.4, 9005.8, 9028.2, 9063.3, 9111.74, 9159.08, 9206.42, 9253.74, 9301.08, 9348.41, 9395.74, 9443.08, 9490.4, 9537.74, 9585.08, 9632.4, 9679.74, 9727.07, 9774.4, 9821.74, 9869.06, 9916.4, 9963.74, 10034.73};
std::vector<double> ErecoilBins = {0.0, 0.2, 0.4, 0.60, 0.8, 1.0, 1.20, 1.40, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0, 10.2, 10.4, 10.6, 10.8, 11.0, 11.2, 11.4, 11.6, 11.8, 12.0, 12.2, 12.4, 12.6, 12.8, 13.0, 13.2, 13.4, 13.6, 13.8, 14.0, 14.2, 14.4, 14.6, 14.8, 15.0, 15.2, 15.4, 15.6, 15.8, 16.0, 16.2, 16.4, 16.6, 16.8, 17.0, 17.2, 17.4, 17.6, 17.8, 18.0, 18.2, 18.4, 18.6, 18.8, 19.0, 19.2, 19.4, 19.6, 19.8, 20};
std::vector<double> ConfBins = {0.0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
std::vector<double> ConfBins2= {0.0, 0.05,0.1,0.15,0.2,0.25,0.3, 0.35,0.4,0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.};
std::vector<double> ResidualBins = {-400.0, -399.0, -398.0, -397.0, -396.0, -395.0, -394.0, -393.0, -392.0, -391.0, -390.0, -389.0, -388.0, -387.0, -386.0, -385.0, -384.0, -383.0, -382.0, -381.0, -380.0, -379.0, -378.0, -377.0, -376.0, -375.0, -374.0, -373.0, -372.0, -371.0, -370.0, -369.0, -368.0, -367.0, -366.0, -365.0, -364.0, -363.0, -362.0, -361.0, -360.0, -359.0, -358.0, -357.0, -356.0, -355.0, -354.0, -353.0, -352.0, -351.0, -350.0, -349.0, -348.0, -347.0, -346.0, -345.0, -344.0, -343.0, -342.0, -341.0, -340.0, -339.0, -338.0, -337.0, -336.0, -335.0, -334.0, -333.0, -332.0, -331.0, -330.0, -329.0, -328.0, -327.0, -326.0, -325.0, -324.0, -323.0, -322.0, -321.0, -320.0, -319.0, -318.0, -317.0, -316.0, -315.0, -314.0, -313.0, -312.0, -311.0, -310.0, -309.0, -308.0, -307.0, -306.0, -305.0, -304.0, -303.0, -302.0, -301.0, -300.0, -299.0, -298.0, -297.0, -296.0, -295.0, -294.0, -293.0, -292.0, -291.0, -290.0, -289.0, -288.0, -287.0, -286.0, -285.0, -284.0, -283.0, -282.0, -281.0, -280.0, -279.0, -278.0, -277.0, -276.0, -275.0, -274.0, -273.0, -272.0, -271.0, -270.0, -269.0, -268.0, -267.0, -266.0, -265.0, -264.0, -263.0, -262.0, -261.0, -260.0, -259.0, -258.0, -257.0, -256.0, -255.0, -254.0, -253.0, -252.0, -251.0, -250.0, -249.0, -248.0, -247.0, -246.0, -245.0, -244.0, -243.0, -242.0, -241.0, -240.0, -239.0, -238.0, -237.0, -236.0, -235.0, -234.0, -233.0, -232.0, -231.0, -230.0, -229.0, -228.0, -227.0, -226.0, -225.0, -224.0, -223.0, -222.0, -221.0, -220.0, -219.0, -218.0, -217.0, -216.0, -215.0, -214.0, -213.0, -212.0, -211.0, -210.0, -209.0, -208.0, -207.0, -206.0, -205.0, -204.0, -203.0, -202.0, -201.0, -200.0, -199.0, -198.0, -197.0, -196.0, -195.0, -194.0, -193.0, -192.0, -191.0, -190.0, -189.0, -188.0, -187.0, -186.0, -185.0, -184.0, -183.0, -182.0, -181.0, -180.0, -179.0, -178.0, -177.0, -176.0, -175.0, -174.0, -173.0, -172.0, -171.0, -170.0, -169.0, -168.0, -167.0, -166.0, -165.0, -164.0, -163.0, -162.0, -161.0, -160.0, -159.0, -158.0, -157.0, -156.0, -155.0, -154.0, -153.0, -152.0, -151.0, -150.0, -149.0, -148.0, -147.0, -146.0, -145.0, -144.0, -143.0, -142.0, -141.0, -140.0, -139.0, -138.0, -137.0, -136.0, -135.0, -134.0, -133.0, -132.0, -131.0, -130.0, -129.0, -128.0, -127.0, -126.0, -125.0, -124.0, -123.0, -122.0, -121.0, -120.0, -119.0, -118.0, -117.0, -116.0, -115.0, -114.0, -113.0, -112.0, -111.0, -110.0, -109.0, -108.0, -107.0, -106.0, -105.0, -104.0, -103.0, -102.0, -101.0, -100.0, -99.0, -98.0, -97.0, -96.0, -95.0, -94.0, -93.0, -92.0, -91.0, -90.0, -89.0, -88.0, -87.0, -86.0, -85.0, -84.0, -83.0, -82.0, -81.0, -80.0, -79.0, -78.0, -77.0, -76.0, -75.0, -74.0, -73.0, -72.0, -71.0, -70.0, -69.0, -68.0, -67.0, -66.0, -65.0, -64.0, -63.0, -62.0, -61.0, -60.0, -59.0, -58.0, -57.0, -56.0, -55.0, -54.0, -53.0, -52.0, -51.0, -50.0, -49.0, -48.0, -47.0, -46.0, -45.0, -44.0, -43.0, -42.0, -41.0, -40.0, -39.0, -38.0, -37.0, -36.0, -35.0, -34.0, -33.0, -32.0, -31.0, -30.0, -29.0, -28.0, -27.0, -26.0, -25.0, -24.0, -23.0, -22.0, -21.0, -20.0, -19.0, -18.0, -17.0, -16.0, -15.0, -14.0, -13.0, -12.0, -11.0, -10.0, -9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0, 100.0, 101.0, 102.0, 103.0, 104.0, 105.0, 106.0, 107.0, 108.0, 109.0, 110.0, 111.0, 112.0, 113.0, 114.0, 115.0, 116.0, 117.0, 118.0, 119.0, 120.0, 121.0, 122.0, 123.0, 124.0, 125.0, 126.0, 127.0, 128.0, 129.0, 130.0, 131.0, 132.0, 133.0, 134.0, 135.0, 136.0, 137.0, 138.0, 139.0, 140.0, 141.0, 142.0, 143.0, 144.0, 145.0, 146.0, 147.0, 148.0, 149.0, 150.0, 151.0, 152.0, 153.0, 154.0, 155.0, 156.0, 157.0, 158.0, 159.0, 160.0, 161.0, 162.0, 163.0, 164.0, 165.0, 166.0, 167.0, 168.0, 169.0, 170.0, 171.0, 172.0, 173.0, 174.0, 175.0, 176.0, 177.0, 178.0, 179.0, 180.0, 181.0, 182.0, 183.0, 184.0, 185.0, 186.0, 187.0, 188.0, 189.0, 190.0, 191.0, 192.0, 193.0, 194.0, 195.0, 196.0, 197.0, 198.0, 199.0, 200.0, 201.0, 202.0, 203.0, 204.0, 205.0, 206.0, 207.0, 208.0, 209.0, 210.0, 211.0, 212.0, 213.0, 214.0, 215.0, 216.0, 217.0, 218.0, 219.0, 220.0, 221.0, 222.0, 223.0, 224.0, 225.0, 226.0, 227.0, 228.0, 229.0, 230.0, 231.0, 232.0, 233.0, 234.0, 235.0, 236.0, 237.0, 238.0, 239.0, 240.0, 241.0, 242.0, 243.0, 244.0, 245.0, 246.0, 247.0, 248.0, 249.0, 250.0, 251.0, 252.0, 253.0, 254.0, 255.0, 256.0, 257.0, 258.0, 259.0, 260.0, 261.0, 262.0, 263.0, 264.0, 265.0, 266.0, 267.0, 268.0, 269.0, 270.0, 271.0, 272.0, 273.0, 274.0, 275.0, 276.0, 277.0, 278.0, 279.0, 280.0, 281.0, 282.0, 283.0, 284.0, 285.0, 286.0, 287.0, 288.0, 289.0, 290.0, 291.0, 292.0, 293.0, 294.0, 295.0, 296.0, 297.0, 298.0, 299.0, 300.0, 301.0, 302.0, 303.0, 304.0, 305.0, 306.0, 307.0, 308.0, 309.0, 310.0, 311.0, 312.0, 313.0, 314.0, 315.0, 316.0, 317.0, 318.0, 319.0, 320.0, 321.0, 322.0, 323.0, 324.0, 325.0, 326.0, 327.0, 328.0, 329.0, 330.0, 331.0, 332.0, 333.0, 334.0, 335.0, 336.0, 337.0, 338.0, 339.0, 340.0, 341.0, 342.0, 343.0, 344.0, 345.0, 346.0, 347.0, 348.0, 349.0, 350.0, 351.0, 352.0, 353.0, 354.0, 355.0, 356.0, 357.0, 358.0, 359.0, 360.0, 361.0, 362.0, 363.0, 364.0, 365.0, 366.0, 367.0, 368.0, 369.0, 370.0, 371.0, 372.0, 373.0, 374.0, 375.0, 376.0, 377.0, 378.0, 379.0, 380.0, 381.0, 382.0, 383.0, 384.0, 385.0, 386.0, 387.0, 388.0, 389.0, 390.0, 391.0, 392.0, 393.0, 394.0, 395.0, 396.0, 397.0, 398.0, 399.0, 400.0};
std::vector<double> multiplicityBins = {0,1,2,3,4,5,6,7,8,9,10};
//std::vector<double> trackHitsBins = {0,1,2,3,4,5,6,7,8,9,10};

double rebinNum = 10;
double nbins = 3400;
double xlow = 4200;
double xhigh = 5900;


//std::unique_ptr<PlotUtils::Cut<UNIVERSE, EVENT> >
auto ANNConfCut = new reco::ANNConfidenceCut<CVUniverse, MichelEvent>(0.40);


//==============================================================================
// Global - Declaring histograms
//==============================================================================
TH2D *TruthZVertex_vs_Erecoil_MC = new TH2D ("TruthZVertex_vs_Erecoil_MC", "TruthZVertex_vs_Erecoil_MC", vertexBins.size()-1, &vertexBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNZ_vs_Conf_vs_Erecoil_Data = new TH3D ("ANNZ_vs_Conf_vs_Erecoil_Data", "ANNZ_vs_Conf_vs_Erecoil_Data", vertexBins.size()-1, &vertexBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNZ_vs_Conf_vs_Erecoil_MC = new TH3D ("ANNZ_vs_Conf_vs_Erecoil_MC", "ANNZ_vs_Conf_vs_Erecoil_MC", vertexBins.size()-1, &vertexBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNWeightedZ_vs_Conf_vs_Erecoil_Data = new TH3D ("ANNWeightedZ_vs_Conf_vs_Erecoil_Data", "ANNWeightedZ_vs_Conf_vs_Erecoil_Data", vertexBins.size()-1, &vertexBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNWeightedZ_vs_Conf_vs_Erecoil_MC = new TH3D ("ANNWeightedZ_vs_Conf_vs_Erecoil_MC", "ANNWeightedZ_vs_Conf_vs_Erecoil_MC", vertexBins.size()-1, &vertexBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);

//With calibration
TH3D *ANNZCalibrated_vs_ConfCalibrated_vs_Erecoil_Data = new TH3D ("ANNZCalibrated_vs_ConfCalibrated_vs_Erecoil_Data", "ANNZCalibrated_vs_ConfCalibrated_vs_Erecoil_Data", vertexBins.size()-1, &vertexBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNZCalibrated_vs_ConfCalibrated_vs_Erecoil_MC = new TH3D ("ANNZCalibrated_vs_ConfCalibrated_vs_Erecoil_MC", "ANNZCalibrated_vs_ConfCalibrated_vs_Erecoil_MC", vertexBins.size()-1, &vertexBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNWeightedZCalibrated_vs_ConfCalibrated_vs_Erecoil_Data = new TH3D ("ANNWeightedZCalibrated_vs_ConfCalibrated_vs_Erecoil_Data", "ANNWeightedZCalibrated_vs_ConfCalibrated_vs_Erecoil_Data", vertexBins.size()-1, &vertexBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNWeightedZCalibrated_vs_ConfCalibrated_vs_Erecoil_MC = new TH3D ("ANNWeightedZCalibrated_vs_ConfCalibrated_vs_Erecoil_MC", "ANNWeightedZCalibrated_vs_ConfCalibrated_vs_Erecoil_MC", vertexBins.size()-1, &vertexBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH1D *ANNWeightedZCalibratedMinusANNZ_Data = new TH1D ("ANNWeightedZCalibratedMinusANNZ_Data", "ANNWeightedZCalibratedMinusANNZ_Data", 101, -100, 100);
TH1D *ANNWeightedZCalibratedMinusANNZ_MC = new TH1D ("ANNWeightedZCalibratedMinusANNZ_MC", "ANNWeightedZCalibratedMinusANNZ_MC", 101, -100, 100);
TH1D *ANNZCalibratedMinusANNZ_Data = new TH1D ("ANNZCalibratedMinusANNZ_Data", "ANNZCalibratedMinusANNZ_Data", 101, -100, 100);
TH1D *ANNZCalibratedMinusANNZ_MC = new TH1D ("ANNZCalibratedMinusANNZ_MC", "ANNZCalibratedMinusANNZ_MC", 101, -100, 100);

TH1D *ANNZCalibrated_Data = new TH1D ("ANNZCalibrated_Data", "ANNZCalibrated_Data", vertexBins.size()-1, &vertexBins[0]);
TH1D *ANNZCalibrated_MC = new TH1D ("ANNZCalibrated_MC", "ANNZCalibrated_MC", vertexBins.size()-1, &vertexBins[0]);
TH1D *ANNZCalibratedWeighted_Data = new TH1D ("ANNZCalibratedWeighted_Data", "ANNZCalibratedWeighted_Data", vertexBins.size()-1, &vertexBins[0]);
TH1D *ANNZCalibratedWeighted_MC = new TH1D ("ANNZCalibratedWeighted_MC", "ANNZCalibratedWeighted_MC", vertexBins.size()-1, &vertexBins[0]);
TH1D *ANNZ_Data = new TH1D ("ANNZ_Data", "ANNZ_Data", vertexBins.size()-1, &vertexBins[0]);
TH1D *ANNZ_MC = new TH1D ("ANNZ_MC", "ANNZ_MC", vertexBins.size()-1, &vertexBins[0]);



TH3D *ANNZCalibratedResidual_vs_ConfCalibrated_vs_Erecoil = new TH3D ("ANNZCalibratedResidual_vs_ConfCalibrated_vs_Erecoil", "ANNZCalibratedResidual_vs_ConfCalibrated_vs_Erecoil", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNZResidual_vs_Conf_vs_Erecoil = new TH3D ("ANNZResidual_vs_Conf_vs_Erecoil", "ANNZResidual_vs_Conf_vs_Erecoil", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNZResidual_vs_Conf_vs_Erecoil_Tracker = new TH3D ("ANNZResidual_vs_Conf_vs_Erecoil_Tracker", "ANNZResidual_vs_Conf_vs_Erecoil_Tracker", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNZResidual_vs_Conf_vs_Erecoil_Tgt2Iron = new TH3D ("ANNZResidual_vs_Conf_vs_Erecoil_Tgt2Iron", "ANNZResidual_vs_Conf_vs_Erecoil_Tgt2Iron", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNZResidual_vs_Conf_vs_Erecoil_Tgt2Lead = new TH3D ("ANNZResidual_vs_Conf_vs_Erecoil_Tgt2Lead", "ANNZResidual_vs_Conf_vs_Erecoil_Tgt2Lead", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNZResidual_vs_Conf_vs_Erecoil_Tgt3Iron = new TH3D ("ANNZResidual_vs_Conf_vs_Erecoil_Tgt3Iron", "ANNZResidual_vs_Conf_vs_Erecoil_Tgt3Iron", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNZResidual_vs_Conf_vs_Erecoil_Tgt3Lead = new TH3D ("ANNZResidual_vs_Conf_vs_Erecoil_Tgt3Lead", "ANNZResidual_vs_Conf_vs_Erecoil_Tgt3Lead", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNZResidual_vs_Conf_vs_Erecoil_Tgt3Carbon = new TH3D ("ANNZResidual_vs_Conf_vs_Erecoil_Tgt3Carbon", "ANNZResidual_vs_Conf_vs_Erecoil_Tgt3Carbon", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNZResidual_vs_Conf_vs_Erecoil_Tgt4Lead = new TH3D ("ANNZResidual_vs_Conf_vs_Erecoil_Tgt4Lead", "ANNZResidual_vs_Conf_vs_Erecoil_Tgt4Lead", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNZResidual_vs_Conf_vs_Erecoil_Tgt5Iron = new TH3D ("ANNZResidual_vs_Conf_vs_Erecoil_Tgt5Iron", "ANNZResidual_vs_Conf_vs_Erecoil_Tgt5Iron", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNZResidual_vs_Conf_vs_Erecoil_Tgt5Lead = new TH3D ("ANNZResidual_vs_Conf_vs_Erecoil_Tgt5Lead", "ANNZResidual_vs_Conf_vs_Erecoil_Tgt5Lead", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNZResidual_vs_Conf_vs_Erecoil_Water = new TH3D ("ANNZResidual_vs_Conf_vs_Erecoil_Water", "ANNZResidual_vs_Conf_vs_Erecoil_Water", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);

TH3D *ANNWeightedZCalibratedResidual_vs_ConfCalibrated_vs_Erecoil = new TH3D ("ANNWeightedZCalibratedResidual_vs_ConfCalibrated_vs_Erecoil", "ANNWeightedZCalibratedResidual_vs_ConfCalibrated_vs_Erecoil", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNWeightedZResidual_vs_Conf_vs_Erecoil = new TH3D ("ANNWeightedZResidual_vs_Conf_vs_Erecoil", "ANNWeightedZResidual_vs_Conf_vs_Erecoil", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tracker = new TH3D ("ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tracker", "ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tracker", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt2Iron = new TH3D ("ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt2Iron", "ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt2Iron", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt2Lead = new TH3D ("ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt2Lead", "ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt2Lead", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt3Iron = new TH3D ("ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt3Iron", "ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt3Iron", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt3Lead = new TH3D ("ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt3Lead", "ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt3Lead", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt3Carbon = new TH3D ("ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt3Carbon", "ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt3Carbon", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt4Lead = new TH3D ("ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt4Lead", "ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt4Lead", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt5Iron = new TH3D ("ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt5Iron", "ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt5Iron", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt5Lead = new TH3D ("ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt5Lead", "ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt5Lead", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ANNWeightedZResidual_vs_Conf_vs_Erecoil_Water = new TH3D ("ANNWeightedZResidual_vs_Conf_vs_Erecoil_Water", "ANNWeightedZResidual_vs_Conf_vs_Erecoil_Water", ResidualBins.size()-1, &ResidualBins[0], ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);

TH1D *ANNWeightedZ_minus_UnweightedZMC = new TH1D ("ANNWeightedZ_minus_UnweightedZMC", "ANNWeightedZ_minus_UnweightedZMC", ResidualBins.size()-1, &ResidualBins[0]);
TH1D *ANNWeightedZ_minus_UnweightedZTrackerMC = new TH1D ("ANNWeightedZ_minus_UnweightedZTrackerMC", "ANNWeightedZ_minus_UnweightedZTrackerMC", ResidualBins.size()-1, &ResidualBins[0]);
TH1D *ANNWeightedZ_minus_UnweightedZTgt2IronMC = new TH1D ("ANNWeightedZ_minus_UnweightedZTgt2IronMC", "ANNWeightedZ_minus_UnweightedZTgt2IronMC", ResidualBins.size()-1, &ResidualBins[0]);
TH1D *ANNWeightedZ_minus_UnweightedZTgt2LeadMC = new TH1D ("ANNWeightedZ_minus_UnweightedZTgt2LeadMC", "ANNWeightedZ_minus_UnweightedZTgt2LeadMC", ResidualBins.size()-1, &ResidualBins[0]);
TH1D *ANNWeightedZ_minus_UnweightedZTgt3IronMC = new TH1D ("ANNWeightedZ_minus_UnweightedZTgt3IronMC", "ANNWeightedZ_minus_UnweightedZTgt3IronMC", ResidualBins.size()-1, &ResidualBins[0]);
TH1D *ANNWeightedZ_minus_UnweightedZTgt3LeadMC = new TH1D ("ANNWeightedZ_minus_UnweightedZTgt3LeadMC", "ANNWeightedZ_minus_UnweightedZTgt3LeadMC", ResidualBins.size()-1, &ResidualBins[0]);
TH1D *ANNWeightedZ_minus_UnweightedZTgt3CarbonMC = new TH1D ("ANNWeightedZ_minus_UnweightedZTgt3CarbonMC", "ANNWeightedZ_minus_UnweightedZTgt3CarbonMC", ResidualBins.size()-1, &ResidualBins[0]);
TH1D *ANNWeightedZ_minus_UnweightedZTgt4LeadMC = new TH1D ("ANNWeightedZ_minus_UnweightedZTgt4LeadMC", "ANNWeightedZ_minus_UnweightedZTgt4LeadMC", ResidualBins.size()-1, &ResidualBins[0]);
TH1D *ANNWeightedZ_minus_UnweightedZTgt5IronMC = new TH1D ("ANNWeightedZ_minus_UnweightedZTgt5IronMC", "ANNWeightedZ_minus_UnweightedZTgt5IronMC", ResidualBins.size()-1, &ResidualBins[0]);
TH1D *ANNWeightedZ_minus_UnweightedZTgt5LeadMC = new TH1D ("ANNWeightedZ_minus_UnweightedZTgt5LeadMC", "ANNWeightedZ_minus_UnweightedZTgt5LeadMC", ResidualBins.size()-1, &ResidualBins[0]);
TH1D *ANNWeightedZ_minus_UnweightedZWaterMC = new TH1D ("ANNWeightedZ_minus_UnweightedZWaterMC", "ANNWeightedZ_minus_UnweightedZWaterMC", ResidualBins.size()-1, &ResidualBins[0]);

TH1D *ANNWeightedZ_minus_UnweightedZData = new TH1D ("ANNWeightedZ_minus_UnweightedZData", "ANNWeightedZ_minus_UnweightedZData", ResidualBins.size()-1, &ResidualBins[0]);
TH1D *ANNWeightedZ_minus_UnweightedZTrackerData = new TH1D ("ANNWeightedZ_minus_UnweightedZTrackerData", "ANNWeightedZ_minus_UnweightedZTrackerData", ResidualBins.size()-1, &ResidualBins[0]);
TH1D *ANNWeightedZ_minus_UnweightedZTgt2IronData = new TH1D ("ANNWeightedZ_minus_UnweightedZTgt2IronData", "ANNWeightedZ_minus_UnweightedZTgt2IronData", ResidualBins.size()-1, &ResidualBins[0]);
TH1D *ANNWeightedZ_minus_UnweightedZTgt2LeadData = new TH1D ("ANNWeightedZ_minus_UnweightedZTgt2LeadData", "ANNWeightedZ_minus_UnweightedZTgt2LeadData", ResidualBins.size()-1, &ResidualBins[0]);
TH1D *ANNWeightedZ_minus_UnweightedZTgt3IronData = new TH1D ("ANNWeightedZ_minus_UnweightedZTgt3IronData", "ANNWeightedZ_minus_UnweightedZTgt3IronData", ResidualBins.size()-1, &ResidualBins[0]);
TH1D *ANNWeightedZ_minus_UnweightedZTgt3LeadData = new TH1D ("ANNWeightedZ_minus_UnweightedZTgt3LeadData", "ANNWeightedZ_minus_UnweightedZTgt3LeadData", ResidualBins.size()-1, &ResidualBins[0]);
TH1D *ANNWeightedZ_minus_UnweightedZTgt3CarbonData = new TH1D ("ANNWeightedZ_minus_UnweightedZTgt3CarbonData", "ANNWeightedZ_minus_UnweightedZTgt3CarbonData", ResidualBins.size()-1, &ResidualBins[0]);
TH1D *ANNWeightedZ_minus_UnweightedZTgt4LeadData = new TH1D ("ANNWeightedZ_minus_UnweightedZTgt4LeadData", "ANNWeightedZ_minus_UnweightedZTgt4LeadData", ResidualBins.size()-1, &ResidualBins[0]);
TH1D *ANNWeightedZ_minus_UnweightedZTgt5IronData = new TH1D ("ANNWeightedZ_minus_UnweightedZTgt5IronData", "ANNWeightedZ_minus_UnweightedZTgt5IronData", ResidualBins.size()-1, &ResidualBins[0]);
TH1D *ANNWeightedZ_minus_UnweightedZTgt5LeadData = new TH1D ("ANNWeightedZ_minus_UnweightedZTgt5LeadData", "ANNWeightedZ_minus_UnweightedZTgt5LeadData", ResidualBins.size()-1, &ResidualBins[0]);
TH1D *ANNWeightedZ_minus_UnweightedZWaterData = new TH1D ("ANNWeightedZ_minus_UnweightedZWaterData", "ANNWeightedZ_minus_UnweightedZWaterData", ResidualBins.size()-1, &ResidualBins[0]);

//TH2D *ANNConf_vs_Ehad_Data = new TH2D ("ANNConf_vs_Ehad_Data", "ANNConf_vs_Ehad_Data", ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]); //Redundant, the same data can be extracted by flattening ANNZ_vs_Conf_vs_Erecoil_Data
TH2D *ANNConf_vs_Ehad_Tracker_Data = new TH2D ("ANNConf_vs_Ehad_Tracker_Data", "ANNConf_vs_Ehad_Tracker_Data", ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH2D *ANNConf_vs_Ehad_Tgt2Iron_Data = new TH2D ("ANNConf_vs_Ehad_Tgt2Iron_Data", "ANNConf_vs_Ehad_Tgt2Iron_Data", ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH2D *ANNConf_vs_Ehad_Tgt2Lead_Data = new TH2D ("ANNConf_vs_Ehad_Tgt2Lead_Data", "ANNConf_vs_Ehad_Tgt2Lead_Data", ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH2D *ANNConf_vs_Ehad_Tgt3Iron_Data = new TH2D ("ANNConf_vs_Ehad_Tgt3Iron_Data", "ANNConf_vs_Ehad_Tgt3Iron_Data", ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH2D *ANNConf_vs_Ehad_Tgt3Lead_Data = new TH2D ("ANNConf_vs_Ehad_Tgt3Lead_Data", "ANNConf_vs_Ehad_Tgt3Lead_Data", ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH2D *ANNConf_vs_Ehad_Tgt3Carbon_Data = new TH2D ("ANNConf_vs_Ehad_Tgt3Carbon_Data", "ANNConf_vs_Ehad_Tgt3Carbon_Data", ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH2D *ANNConf_vs_Ehad_Tgt4Lead_Data = new TH2D ("ANNConf_vs_Ehad_Tgt4Lead_Data", "ANNConf_vs_Ehad_Tgt4Lead_Data", ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH2D *ANNConf_vs_Ehad_Tgt5Iron_Data = new TH2D ("ANNConf_vs_Ehad_Tgt5Iron_Data", "ANNConf_vs_Ehad_Tgt5Iron_Data", ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH2D *ANNConf_vs_Ehad_Tgt5Lead_Data = new TH2D ("ANNConf_vs_Ehad_Tgt5Lead_Data", "ANNConf_vs_Ehad_Tgt5Lead_Data", ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH2D *ANNConf_vs_Ehad_Water_Data = new TH2D ("ANNConf_vs_Ehad_Water_Data", "ANNConf_vs_Ehad_Water_Data", ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);

//TH2D *ANNConf_vs_Ehad_MC = new TH2D ("ANNConf_vs_Ehad_MC", "ANNConf_vs_Ehad_MC", ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);  //Redundant, the same data can be extracted by flattening ANNZ_vs_Conf_vs_Erecoil_MC
TH2D *ANNConf_vs_Ehad_Tracker_MC = new TH2D ("ANNConf_vs_Ehad_Tracker_MC", "ANNConf_vs_Ehad_Tracker_MC", ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH2D *ANNConf_vs_Ehad_Tgt2Iron_MC = new TH2D ("ANNConf_vs_Ehad_Tgt2Iron_MC", "ANNConf_vs_Ehad_Tgt2Iron_MC", ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH2D *ANNConf_vs_Ehad_Tgt2Lead_MC = new TH2D ("ANNConf_vs_Ehad_Tgt2Lead_MC", "ANNConf_vs_Ehad_Tgt2Lead_MC", ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH2D *ANNConf_vs_Ehad_Tgt3Iron_MC = new TH2D ("ANNConf_vs_Ehad_Tgt3Iron_MC", "ANNConf_vs_Ehad_Tgt3Iron_MC", ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH2D *ANNConf_vs_Ehad_Tgt3Lead_MC = new TH2D ("ANNConf_vs_Ehad_Tgt3Lead_MC", "ANNConf_vs_Ehad_Tgt3Lead_MC", ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH2D *ANNConf_vs_Ehad_Tgt3Carbon_MC = new TH2D ("ANNConf_vs_Ehad_Tgt3Carbon_MC", "ANNConf_vs_Ehad_Tgt3Carbon_MC", ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH2D *ANNConf_vs_Ehad_Tgt4Lead_MC = new TH2D ("ANNConf_vs_Ehad_Tgt4Lead_MC", "ANNConf_vs_Ehad_Tgt4Lead_MC", ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH2D *ANNConf_vs_Ehad_Tgt5Iron_MC = new TH2D ("ANNConf_vs_Ehad_Tgt5Iron_MC", "ANNConf_vs_Ehad_Tgt5Iron_MC", ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH2D *ANNConf_vs_Ehad_Tgt5Lead_MC = new TH2D ("ANNConf_vs_Ehad_Tgt5Lead_MC", "ANNConf_vs_Ehad_Tgt5Lead_MC", ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH2D *ANNConf_vs_Ehad_Water_MC = new TH2D ("ANNConf_vs_Ehad_Water_MC", "ANNConf_vs_Ehad_Water_MC", ConfBins.size()-1, &ConfBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);

//Multiplicity - Explore (re)implementing multiplicty somehow - is that even worth doing?

//For events misreconstructed in each target, what's the probability distribution
TH1D *ANNConf_MisrecoInTgt2Iron = new TH1D ("ANNConf_MisrecoInTgt2Iron", "ANNConf_MisrecoInTgt2Iron", ConfBins.size()-1, &ConfBins[0]);
TH1D *ANNConf_MisrecoInTgt2Lead = new TH1D ("ANNConf_MisrecoInTgt2Lead", "ANNConf_MisrecoInTgt2Lead", ConfBins.size()-1, &ConfBins[0]);
TH1D *ANNConf_MisrecoInTgt3Iron = new TH1D ("ANNConf_MisrecoInTgt3Iron", "ANNConf_MisrecoInTgt3Iron", ConfBins.size()-1, &ConfBins[0]);
TH1D *ANNConf_MisrecoInTgt3Lead = new TH1D ("ANNConf_MisrecoInTgt3Lead", "ANNConf_MisrecoInTgt3Lead", ConfBins.size()-1, &ConfBins[0]);
TH1D *ANNConf_MisrecoInTgt3Carbon = new TH1D ("ANNConf_MisrecoInTgt3Carbon", "ANNConf_MisrecoInTgt3Carbon", ConfBins.size()-1, &ConfBins[0]);
TH1D *ANNConf_MisrecoInTgt4Lead = new TH1D ("ANNConf_MisrecoInTgt4Lead", "ANNConf_MisrecoInTgt4Lead", ConfBins.size()-1, &ConfBins[0]);
TH1D *ANNConf_MisrecoInTgt5Iron = new TH1D ("ANNConf_MisrecoInTgt5Iron", "ANNConf_MisrecoInTgt5Iron", ConfBins.size()-1, &ConfBins[0]);
TH1D *ANNConf_MisrecoInTgt5Lead = new TH1D ("ANNConf_MisrecoInTgt5Lead", "ANNConf_MisrecoInTgt5Lead", ConfBins.size()-1, &ConfBins[0]);
TH1D *ANNConf_MisrecoInWater = new TH1D ("ANNConf_MisrecoInWater", "ANNConf_MisrecoInWater", ConfBins.size()-1, &ConfBins[0]);

//Efficiency as a function of ANNConfidence and EHad bin

//Migration ANN only
TH2D *ANNVerticesConfusion_ByZPos = new TH2D ("ANNVerticesConfusion_ByZPos", "ANNVerticesConfusion_ByZPos", 4500, 4200, 8700, 4500, 4200, 8700);

//Investigating distance to division cut performance. We're plotting the 2D postions of the events based on 
//Investigating ANN reco performance by MC interaction channel
//Confusion matrices
TH2D *ConfusionMatrixQE = new TH2D ("ConfusionMatrixQE", "ConfusionMatrixQE", vertexBins.size()-1, &vertexBins[0], vertexBins.size()-1, &vertexBins[0]);
TH2D *ConfusionMatrixRES = new TH2D ("ConfusionMatrixRES", "ConfusionMatrixRES", vertexBins.size()-1, &vertexBins[0], vertexBins.size()-1, &vertexBins[0]);
TH2D *ConfusionMatrixDIS = new TH2D ("ConfusionMatrixDIS", "ConfusionMatrixDIS", vertexBins.size()-1, &vertexBins[0], vertexBins.size()-1, &vertexBins[0]);
TH2D *ConfusionMatrixCOHPI = new TH2D ("ConfusionMatrixCOHPI", "ConfusionMatrixCOHPI", vertexBins.size()-1, &vertexBins[0], vertexBins.size()-1, &vertexBins[0]);
TH2D *ConfusionMatrixAMNUGAMMA = new TH2D ("ConfusionMatrixAMNUGAMMA", "ConfusionMatrixAMNUGAMMA", vertexBins.size()-1, &vertexBins[0], vertexBins.size()-1, &vertexBins[0]);
TH2D *ConfusionMatrixIMD = new TH2D ("ConfusionMatrixIMD", "ConfusionMatrixIMD", vertexBins.size()-1, &vertexBins[0], vertexBins.size()-1, &vertexBins[0]);
TH2D *ConfusionMatrixNUEEL = new TH2D ("ConfusionMatrixNUEEL", "ConfusionMatrixNUEEL", vertexBins.size()-1, &vertexBins[0], vertexBins.size()-1, &vertexBins[0]);
TH2D *ConfusionMatrixMEC = new TH2D ("ConfusionMatrixMEC", "ConfusionMatrixMEC", vertexBins.size()-1, &vertexBins[0], vertexBins.size()-1, &vertexBins[0]);
TH2D *ConfusionMatrixUnknownInt = new TH2D ("ConfusionMatrixUnknownInt", "ConfusionMatrixUnknownInt", vertexBins.size()-1, &vertexBins[0], vertexBins.size()-1, &vertexBins[0]);
//Number of tracks
TH3D *ConfusionMatrixMultiplicity = new TH3D ("ConfusionMatrixMultiplicity", "ConfusionMatrixMultiplicity", vertexBins.size()-1, &vertexBins[0], vertexBins.size()-1, &vertexBins[0], multiplicityBins.size()-1, &multiplicityBins[0]);
//ERecoil
TH3D *ConfusionMatrixErecoil = new TH3D ("ConfusionMatrixErecoil", "ConfusionMatrixErecoil", vertexBins.size()-1, &vertexBins[0], vertexBins.size()-1, &vertexBins[0], ErecoilBins.size()-1, &ErecoilBins[0]);
TH3D *ConfusionMatrixConf = new TH3D ("ConfusionMatrixConf", "ConfusionMatrixConf", vertexBins.size()-1, &vertexBins[0], vertexBins.size()-1, &vertexBins[0], ConfBins.size()-1, &ConfBins[0]);

TH2D *ConfusionMatrix = new TH2D ("ConfusionMatrix", "ConfusionMatrix", vertexBins.size()-1, &vertexBins[0], vertexBins.size()-1, &vertexBins[0]);

//Confidence Breakdowns
// by interaction type
// MC
TH1D *ANNConfidenceQE_MC = new TH1D ("ANNConfidenceQE_MC", "ANNConfidenceQE_MC", ConfBins2.size()-1, &ConfBins2[0]);
TH1D *ANNConfidenceRES_MC = new TH1D ("ANNConfidenceRES_MC", "ANNConfidenceRES_MC", ConfBins2.size()-1, &ConfBins2[0]);
TH1D *ANNConfidenceDIS_MC = new TH1D ("ANNConfidenceDIS_MC", "ANNConfidenceDIS_MC", ConfBins2.size()-1, &ConfBins2[0]);
TH1D *ANNConfidenceCOHPI_MC = new TH1D ("ANNConfidenceCOHPI_MC", "ANNConfidenceCOHPI_MC", ConfBins2.size()-1, &ConfBins2[0]);
TH1D *ANNConfidenceAMNUGAMMA_MC = new TH1D ("ANNConfidenceAMNUGAMMA_MC", "ANNConfidenceAMNUGAMMA_MC", ConfBins2.size()-1, &ConfBins2[0]);
TH1D *ANNConfidenceIMD_MC = new TH1D ("ANNConfidenceIMD_MC", "ANNConfidenceIMD_MC", ConfBins2.size()-1, &ConfBins2[0]);
TH1D *ANNConfidenceNUEEL_MC = new TH1D ("ANNConfidenceNUEEL_MC", "ANNConfidenceNUEEL_MC", ConfBins2.size()-1, &ConfBins2[0]);
TH1D *ANNConfidenceMEC_MC = new TH1D ("ANNConfidenceMEC_MC", "ANNConfidenceMEC_MC", ConfBins2.size()-1, &ConfBins2[0]);
TH1D *ANNConfidenceUnknownInt_MC = new TH1D ("ANNConfidenceUnknownInt_MC", "ANNConfidenceUnknownInt_MC", ConfBins2.size()-1, &ConfBins2[0]);
// by multiplicity
TH2D *ANNConfidence_vs_NTracks_MC = new TH2D ("ANNConfidence_vs_NTracks_MC", "ANNConfidence_vs_NTracks_MC", ConfBins2.size()-1, &ConfBins2[0], multiplicityBins.size()-1, &multiplicityBins[0]);
TH2D *ANNConfidence_vs_NTracks_Data = new TH2D ("ANNConfidence_vs_NTracks_Data", "ANNConfidence_vs_NTracks_Data", ConfBins2.size()-1, &ConfBins2[0], multiplicityBins.size()-1, &multiplicityBins[0]);

TH1D *ANNConfidence_Data = new TH1D ("ANNConfidence_Data", "ANNConfidence_Data", 100, 0, 1);
TH1D *ANNConfidence_MC = new TH1D ("ANNConfidence_MC", "ANNConfidence_MC", 100, 0, 1);
TH1D *ANNCalibratedConfidence_Data = new TH1D ("ANNCalibratedConfidence_Data", "ANNCalibratedConfidence_Data", 100, 0, 1);
TH1D *ANNCalibratedConfidence_MC = new TH1D ("ANNCalibratedConfidence_MC", "ANNCalibratedConfidence_MC", 100, 0, 1);


//Iron
TH1D *ANNCalibratedConfidence_Iron_Selected = new TH1D ("ANNCalibratedConfidence_Iron_Selected", "ANNCalibratedConfidence_Iron_Selected", 100, 0, 1);
TH1D *ANNCalibratedConfidence_Iron_SelectedSignal = new TH1D ("ANNCalibratedConfidence_Iron_SelectedSignal", "ANNCalibratedConfidence_Iron_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidence_Iron_Signal = new TH1D ("ANNCalibratedConfidence_Iron_Signal", "ANNCalibratedConfidence_Iron_Signal", 100, 0, 1);

TH1D *ANNCalibratedConfidence1_Iron_SelectedSignal = new TH1D ("ANNCalibratedConfidence1_Iron_SelectedSignal", "ANNCalibratedConfidence1_Iron_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidence1_Iron_Signal = new TH1D ("ANNCalibratedConfidence1_Iron_Signal", "ANNCalibratedConfidence1_Iron_Signal", 100, 0, 1);
TH1D *ANNCalibratedConfidence1_Iron_Selected = new TH1D ("ANNCalibratedConfidence1_Iron_Selected", "ANNCalibratedConfidence1_Iron_Selected", 100, 0, 1);

TH1D *ANNCalibratedConfidence2_Iron_SelectedSignal = new TH1D ("ANNCalibratedConfidence2_Iron_SelectedSignal", "ANNCalibratedConfidence2_Iron_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidence2_Iron_Signal = new TH1D ("ANNCalibratedConfidence2_Iron_Signal", "ANNCalibratedConfidence2_Iron_Signal", 100, 0, 1);
TH1D *ANNCalibratedConfidence2_Iron_Selected = new TH1D ("ANNCalibratedConfidence2_Iron_Selected", "ANNCalibratedConfidence2_Iron_Selected", 100, 0, 1);

TH1D *ANNCalibratedConfidence3_Iron_SelectedSignal = new TH1D ("ANNCalibratedConfidence3_Iron_SelectedSignal", "ANNCalibratedConfidence3_Iron_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidence3_Iron_Signal = new TH1D ("ANNCalibratedConfidence3_Iron_Signal", "ANNCalibratedConfidence3_Iron_Signal", 100, 0, 1);
TH1D *ANNCalibratedConfidence3_Iron_Selected = new TH1D ("ANNCalibratedConfidence3_Iron_Selected", "ANNCalibratedConfidence3_Iron_Selected", 100, 0, 1);
//Lead
TH1D *ANNCalibratedConfidence_Lead_Selected = new TH1D ("ANNCalibratedConfidence_Lead_Selected", "ANNCalibratedConfidence_Lead_Selected", 100, 0, 1);
TH1D *ANNCalibratedConfidence_Lead_SelectedSignal = new TH1D ("ANNCalibratedConfidence_Lead_SelectedSignal", "ANNCalibratedConfidence_Lead_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidence_Lead_Signal = new TH1D ("ANNCalibratedConfidence_Lead_Signal", "ANNCalibratedConfidence_Lead_Signal", 100, 0, 1);

TH1D *ANNCalibratedConfidence1_Lead_SelectedSignal = new TH1D ("ANNCalibratedConfidence1_Lead_SelectedSignal", "ANNCalibratedConfidence1_Lead_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidence1_Lead_Signal = new TH1D ("ANNCalibratedConfidence1_Lead_Signal", "ANNCalibratedConfidence1_Lead_Signal", 100, 0, 1);
TH1D *ANNCalibratedConfidence1_Lead_Selected = new TH1D ("ANNCalibratedConfidence1_Lead_Selected", "ANNCalibratedConfidence1_Lead_Selected", 100, 0, 1);

TH1D *ANNCalibratedConfidence2_Lead_SelectedSignal = new TH1D ("ANNCalibratedConfidence2_Lead_SelectedSignal", "ANNCalibratedConfidence2_Lead_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidence2_Lead_Signal = new TH1D ("ANNCalibratedConfidence2_Lead_Signal", "ANNCalibratedConfidence2_Lead_Signal", 100, 0, 1);
TH1D *ANNCalibratedConfidence2_Lead_Selected = new TH1D ("ANNCalibratedConfidence2_Lead_Selected", "ANNCalibratedConfidence2_Lead_Selected", 100, 0, 1);

TH1D *ANNCalibratedConfidence3_Lead_SelectedSignal = new TH1D ("ANNCalibratedConfidence3_Lead_SelectedSignal", "ANNCalibratedConfidence3_Lead_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidence3_Lead_Signal = new TH1D ("ANNCalibratedConfidence3_Lead_Signal", "ANNCalibratedConfidence3_Lead_Signal", 100, 0, 1);
TH1D *ANNCalibratedConfidence3_Lead_Selected = new TH1D ("ANNCalibratedConfidence3_Lead_Selected", "ANNCalibratedConfidence3_Lead_Selected", 100, 0, 1);
//Carbon
TH1D *ANNCalibratedConfidence_Carbon_Selected = new TH1D ("ANNCalibratedConfidence_Carbon_Selected", "ANNCalibratedConfidence_Carbon_Selected", 100, 0, 1);
TH1D *ANNCalibratedConfidence_Carbon_SelectedSignal = new TH1D ("ANNCalibratedConfidence_Carbon_SelectedSignal", "ANNCalibratedConfidence_Carbon_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidence_Carbon_Signal = new TH1D ("ANNCalibratedConfidence_Carbon_Signal", "ANNCalibratedConfidence_Carbon_Signal", 100, 0, 1);

TH1D *ANNCalibratedConfidence1_Carbon_SelectedSignal = new TH1D ("ANNCalibratedConfidence1_Carbon_SelectedSignal", "ANNCalibratedConfidence1_Carbon_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidence1_Carbon_Signal = new TH1D ("ANNCalibratedConfidence1_Carbon_Signal", "ANNCalibratedConfidence1_Carbon_Signal", 100, 0, 1);
TH1D *ANNCalibratedConfidence1_Carbon_Selected = new TH1D ("ANNCalibratedConfidence1_Carbon_Selected", "ANNCalibratedConfidence1_Carbon_Selected", 100, 0, 1);

TH1D *ANNCalibratedConfidence2_Carbon_SelectedSignal = new TH1D ("ANNCalibratedConfidence2_Carbon_SelectedSignal", "ANNCalibratedConfidence2_Carbon_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidence2_Carbon_Signal = new TH1D ("ANNCalibratedConfidence2_Carbon_Signal", "ANNCalibratedConfidence2_Carbon_Signal", 100, 0, 1);
TH1D *ANNCalibratedConfidence2_Carbon_Selected = new TH1D ("ANNCalibratedConfidence2_Carbon_Selected", "ANNCalibratedConfidence2_Carbon_Selected", 100, 0, 1);

TH1D *ANNCalibratedConfidence3_Carbon_SelectedSignal = new TH1D ("ANNCalibratedConfidence3_Carbon_SelectedSignal", "ANNCalibratedConfidence3_Carbon_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidence3_Carbon_Signal = new TH1D ("ANNCalibratedConfidence3_Carbon_Signal", "ANNCalibratedConfidence3_Carbon_Signal", 100, 0, 1);
TH1D *ANNCalibratedConfidence3_Carbon_Selected = new TH1D ("ANNCalibratedConfidence3_Carbon_Selected", "ANNCalibratedConfidence3_Carbon_Selected", 100, 0, 1);
//Water
TH1D *ANNCalibratedConfidence_Water_Selected = new TH1D ("ANNCalibratedConfidence_Water_Selected", "ANNCalibratedConfidence_Water_Selected", 100, 0, 1);
TH1D *ANNCalibratedConfidence_Water_SelectedSignal = new TH1D ("ANNCalibratedConfidence_Water_SelectedSignal", "ANNCalibratedConfidence_Water_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidence_Water_Signal = new TH1D ("ANNCalibratedConfidence_Water_Signal", "ANNCalibratedConfidence_Water_Signal", 100, 0, 1);

TH1D *ANNCalibratedConfidence1_Water_SelectedSignal = new TH1D ("ANNCalibratedConfidence1_Water_SelectedSignal", "ANNCalibratedConfidence1_Water_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidence1_Water_Signal = new TH1D ("ANNCalibratedConfidence1_Water_Signal", "ANNCalibratedConfidence1_Water_Signal", 100, 0, 1);
TH1D *ANNCalibratedConfidence1_Water_Selected = new TH1D ("ANNCalibratedConfidence1_Water_Selected", "ANNCalibratedConfidence1_Water_Selected", 100, 0, 1);

TH1D *ANNCalibratedConfidence2_Water_SelectedSignal = new TH1D ("ANNCalibratedConfidence2_Water_SelectedSignal", "ANNCalibratedConfidence2_Water_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidence2_Water_Signal = new TH1D ("ANNCalibratedConfidence2_Water_Signal", "ANNCalibratedConfidence2_Water_Signal", 100, 0, 1);
TH1D *ANNCalibratedConfidence2_Water_Selected = new TH1D ("ANNCalibratedConfidence2_Water_Selected", "ANNCalibratedConfidence2_Water_Selected", 100, 0, 1);

TH1D *ANNCalibratedConfidence3_Water_SelectedSignal = new TH1D ("ANNCalibratedConfidence3_Water_SelectedSignal", "ANNCalibratedConfidence3_Water_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidence3_Water_Signal = new TH1D ("ANNCalibratedConfidence3_Water_Signal", "ANNCalibratedConfidence3_Water_Signal", 100, 0, 1);
TH1D *ANNCalibratedConfidence3_Water_Selected = new TH1D ("ANNCalibratedConfidence3_Water_Selected", "ANNCalibratedConfidence3_Water_Selected", 100, 0, 1);


//Iron

TH1D *ANNCalibratedConfidenceRatio2_Iron_SelectedSignal = new TH1D ("ANNCalibratedConfidenceRatio2_Iron_SelectedSignal", "ANNCalibratedConfidenceRatio2_Iron_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidenceRatio2_Iron_Signal = new TH1D ("ANNCalibratedConfidenceRatio2_Iron_Signal", "ANNCalibratedConfidenceRatio2_Iron_Signal", 100, 0, 1);
TH1D *ANNCalibratedConfidenceRatio2_Iron_Selected = new TH1D ("ANNCalibratedConfidenceRatio2_Iron_Selected", "ANNCalibratedConfidenceRatio2_Iron_Selected", 100, 0, 1);

TH1D *ANNCalibratedConfidenceRatio3_Iron_SelectedSignal = new TH1D ("ANNCalibratedConfidenceRatio3_Iron_SelectedSignal", "ANNCalibratedConfidenceRatio3_Iron_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidenceRatio3_Iron_Signal = new TH1D ("ANNCalibratedConfidenceRatio3_Iron_Signal", "ANNCalibratedConfidenceRatio3_Iron_Signal", 100, 0, 1);
TH1D *ANNCalibratedConfidenceRatio3_Iron_Selected = new TH1D ("ANNCalibratedConfidenceRatio3_Iron_Selected", "ANNCalibratedConfidenceRatio3_Iron_Selected", 100, 0, 1);
//Lead

TH1D *ANNCalibratedConfidenceRatio2_Lead_SelectedSignal = new TH1D ("ANNCalibratedConfidenceRatio2_Lead_SelectedSignal", "ANNCalibratedConfidenceRatio2_Lead_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidenceRatio2_Lead_Signal = new TH1D ("ANNCalibratedConfidenceRatio2_Lead_Signal", "ANNCalibratedConfidenceRatio2_Lead_Signal", 100, 0, 1);
TH1D *ANNCalibratedConfidenceRatio2_Lead_Selected = new TH1D ("ANNCalibratedConfidenceRatio2_Lead_Selected", "ANNCalibratedConfidenceRatio2_Lead_Selected", 100, 0, 1);

TH1D *ANNCalibratedConfidenceRatio3_Lead_SelectedSignal = new TH1D ("ANNCalibratedConfidenceRatio3_Lead_SelectedSignal", "ANNCalibratedConfidenceRatio3_Lead_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidenceRatio3_Lead_Signal = new TH1D ("ANNCalibratedConfidenceRatio3_Lead_Signal", "ANNCalibratedConfidenceRatio3_Lead_Signal", 100, 0, 1);
TH1D *ANNCalibratedConfidenceRatio3_Lead_Selected = new TH1D ("ANNCalibratedConfidenceRatio3_Lead_Selected", "ANNCalibratedConfidenceRatio3_Lead_Selected", 100, 0, 1);
//Carbon

TH1D *ANNCalibratedConfidenceRatio2_Carbon_SelectedSignal = new TH1D ("ANNCalibratedConfidenceRatio2_Carbon_SelectedSignal", "ANNCalibratedConfidenceRatio2_Carbon_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidenceRatio2_Carbon_Signal = new TH1D ("ANNCalibratedConfidenceRatio2_Carbon_Signal", "ANNCalibratedConfidenceRatio2_Carbon_Signal", 100, 0, 1);
TH1D *ANNCalibratedConfidenceRatio2_Carbon_Selected = new TH1D ("ANNCalibratedConfidenceRatio2_Carbon_Selected", "ANNCalibratedConfidenceRatio2_Carbon_Selected", 100, 0, 1);

TH1D *ANNCalibratedConfidenceRatio3_Carbon_SelectedSignal = new TH1D ("ANNCalibratedConfidenceRatio3_Carbon_SelectedSignal", "ANNCalibratedConfidenceRatio3_Carbon_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidenceRatio3_Carbon_Signal = new TH1D ("ANNCalibratedConfidenceRatio3_Carbon_Signal", "ANNCalibratedConfidenceRatio3_Carbon_Signal", 100, 0, 1);
TH1D *ANNCalibratedConfidenceRatio3_Carbon_Selected = new TH1D ("ANNCalibratedConfidenceRatio3_Carbon_Selected", "ANNCalibratedConfidenceRatio3_Carbon_Selected", 100, 0, 1);
//Water

TH1D *ANNCalibratedConfidenceRatio2_Water_SelectedSignal = new TH1D ("ANNCalibratedConfidenceRatio2_Water_SelectedSignal", "ANNCalibratedConfidenceRatio2_Water_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidenceRatio2_Water_Signal = new TH1D ("ANNCalibratedConfidenceRatio2_Water_Signal", "ANNCalibratedConfidenceRatio2_Water_Signal", 100, 0, 1);
TH1D *ANNCalibratedConfidenceRatio2_Water_Selected = new TH1D ("ANNCalibratedConfidenceRatio2_Water_Selected", "ANNCalibratedConfidenceRatio2_Water_Selected", 100, 0, 1);

TH1D *ANNCalibratedConfidenceRatio3_Water_SelectedSignal = new TH1D ("ANNCalibratedConfidenceRatio3_Water_SelectedSignal", "ANNCalibratedConfidenceRatio3_Water_SelectedSignal", 100, 0, 1);
TH1D *ANNCalibratedConfidenceRatio3_Water_Signal = new TH1D ("ANNCalibratedConfidenceRatio3_Water_Signal", "ANNCalibratedConfidenceRatio3_Water_Signal", 100, 0, 1);
TH1D *ANNCalibratedConfidenceRatio3_Water_Selected = new TH1D ("ANNCalibratedConfidenceRatio3_Water_Selected", "ANNCalibratedConfidenceRatio3_Water_Selected", 100, 0, 1);

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

    double scintillatorthickness = 10;

    TBox *tgt1 = new TBox(tgt1centre-tgt1thickness/2, ymin, tgt1centre+tgt1thickness/2, ymax);
    TBox *tgt2 = new TBox(tgt2centre-tgt2thickness/2, ymin, tgt2centre+tgt2thickness/2, ymax);
    TBox *tgt3 = new TBox(tgt3centre-tgt3thickness/2, ymin, tgt3centre+tgt3thickness/2, ymax);
    TBox *tgt4 = new TBox(tgt4centre-tgt4thickness/2, ymin, tgt4centre+tgt4thickness/2, ymax);
    TBox *tgt5 = new TBox(tgt5centre-tgt5thickness/2, ymin, tgt5centre+tgt5thickness/2, ymax);
    TBox *tgt6 = new TBox(5200, ymin, 5420, ymax);

    TBox *tgt1US = new TBox(4446.33-scintillatorthickness, ymin, 4446.33+scintillatorthickness, ymax);
    TBox *tgt2US = new TBox(4667.4-scintillatorthickness, ymin, 4667.4+scintillatorthickness, ymax);
    TBox *tgt3US = new TBox(4888.48-scintillatorthickness, ymin, 4888.48+scintillatorthickness, ymax);
    TBox *tgt4US = new TBox(5610.02-scintillatorthickness, ymin, 5610.02+scintillatorthickness, ymax);
    TBox *tgt5US = new TBox(5742.67-scintillatorthickness, ymin, 5742.67+scintillatorthickness, ymax);
    TBox *tgt6US = new TBox(5153.77-scintillatorthickness, ymin, 5153.77+scintillatorthickness, ymax);


    TBox *tgt1DS = new TBox(4514.11-scintillatorthickness, ymin, 4514.11+scintillatorthickness, ymax);
    TBox *tgt2DS = new TBox(4735.19-scintillatorthickness, ymin, 4735.19+scintillatorthickness, ymax);
    TBox *tgt3DS = new TBox(5000.48-scintillatorthickness, ymin, 5000.48+scintillatorthickness, ymax);
    TBox *tgt4DS = new TBox(5677.81-scintillatorthickness, ymin, 5677.81+scintillatorthickness, ymax);
    TBox *tgt5DS = new TBox(5810.45-scintillatorthickness, ymin, 5810.45+scintillatorthickness, ymax);
    TBox *tgt6DS = new TBox(5456.74-scintillatorthickness, ymin, 5456.74+scintillatorthickness, ymax);

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
    resLabel.SetTextSize(0.035); 
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



    tgt1US->SetFillColorAlpha(kRed,0.5);
    //tgt1US->SetFillStyle(3004);
    tgt1US->SetLineWidth(0);
    tgt1US->SetLineColor(1);
    tgt1US->Draw("same");
    tgt2US->SetFillColorAlpha(kRed,0.5);
    //tgt2US->SetFillStyle(3004);
    tgt2US->SetLineWidth(0);
    tgt2US->SetLineColor(1);
    tgt2US->Draw("same");
    tgt3US->SetFillColorAlpha(kRed,0.5);
    //tgt3US->SetFillStyle(3004);
    tgt3US->SetLineWidth(0);
    tgt3US->SetLineColor(1);
    tgt3US->Draw("same");
    tgt4US->SetFillColorAlpha(kRed,0.5);
    //tgt4US->SetFillStyle(3004);
    tgt4US->SetLineWidth(0);
    tgt4US->SetLineColor(1);
    tgt4US->Draw("same");
    tgt5US->SetFillColorAlpha(kRed,0.5);
    //tgt5US->SetFillStyle(3004);
    tgt5US->SetLineWidth(0);
    tgt5US->SetLineColor(1);
    tgt5US->Draw("same");
    tgt6US->SetFillColorAlpha(kRed,0.5);
    //tgt6US->SetFillStyle(3004);
    tgt6US->SetLineWidth(0);
    tgt6US->SetLineColor(1);
    tgt6US->Draw("same");


    tgt1DS->SetFillColorAlpha(kGreen,0.5);
    //tgt1DS->SetFillStyle(3004);
    tgt1DS->SetLineWidth(0);
    tgt1DS->SetLineColor(1);
    tgt1DS->Draw("same");
    tgt2DS->SetFillColorAlpha(kGreen,0.5);
    //tgt2DS->SetFillStyle(3004);
    tgt2DS->SetLineWidth(0);
    tgt2DS->SetLineColor(1);
    tgt2DS->Draw("same");
    tgt3DS->SetFillColorAlpha(kGreen,0.5);
    //tgt3DS->SetFillStyle(3004);
    tgt3DS->SetLineWidth(0);
    tgt3DS->SetLineColor(1);
    tgt3DS->Draw("same");
    tgt4DS->SetFillColorAlpha(kGreen,0.5);
    //tgt4DS->SetFillStyle(3004);
    tgt4DS->SetLineWidth(0);
    tgt4DS->SetLineColor(1);
    tgt4DS->Draw("same");
    tgt5DS->SetFillColorAlpha(kGreen,0.5);
    //tgt5DS->SetFillStyle(3004);
    tgt5DS->SetLineWidth(0);
    tgt5DS->SetLineColor(1);
    tgt5DS->Draw("same");
    tgt6DS->SetFillColorAlpha(kGreen,0.5);
    //tgt6DS->SetFillStyle(3004);
    tgt6DS->SetLineWidth(0);
    tgt6DS->SetLineColor(1);
    tgt6DS->Draw("same");

}

void DrawDataMC(TH1D* DataIn, TH1D* MCIn, std::string name, double mcpot, double datapot, bool ZoomTargets = false)
{
  double xaxismin = 4280;
  double xaxismax = ZoomTargets ? 6000 : 8400;
  TCanvas *canvas = new TCanvas("canvas", "canvas", 1500, 1000);
  TH1D *DataHist = (TH1D*)DataIn->Clone();
  TH1D *MCHist = (TH1D*)MCIn->Clone();
  DataHist->GetXaxis()->SetRangeUser(xaxismin, xaxismax);
  MCHist->GetXaxis()->SetRangeUser(xaxismin, xaxismax);
  MCHist->SetTitle("MC");
  DataHist->SetTitle("Data");
  name.erase(std::remove(name.begin(), name.end(), ' '), name.end()); //remove spaces
  name.erase(std::remove(name.begin(), name.end(), '-'), name.end()); //remove dashes
  PlotUtils::MnvPlotter plotter;
  plotter.ApplyStyle(PlotUtils::kDefaultStyle);
  //plotter.ApplyStyle(PlotUtils::kCCQENuStyle);
  plotter.DrawDataMC(DataHist, MCHist, datapot/mcpot, "TL", true);
  double ymin = gPad->GetFrame()->GetY1();
  double ymax = gPad->GetFrame()->GetY2();
  PlotRegions(ymin, ymax, mcpot, datapot);
  canvas->SaveAs((std::string("./")+name+std::string(".png")).c_str());
  canvas->Clear();
  //canvas->SetLogy(0);
  plotter.DrawDataMCRatio(DataHist, MCHist, datapot/mcpot, true, 0.0, 2.0, "Data/MC Ratio");
  ymin = gPad->GetFrame()->GetY1();
  ymax = gPad->GetFrame()->GetY2();
  PlotRegions(ymin, ymax, mcpot, datapot);
  canvas->SaveAs((std::string("./")+name+std::string("Ratio.png")).c_str());

  canvas->Clear();
  canvas->SetLogy();
  DataHist->SetMinimum(10);
  MCHist->SetMinimum(10);
  plotter.DrawDataMC(DataHist, MCHist, datapot/mcpot, "TL", true);
  ymin = gPad->GetFrame()->GetY1();
  ymax = gPad->GetFrame()->GetY2();
  PlotRegions(std::pow(10,ymin), std::pow(10,ymax), mcpot, datapot);
  TH1D* tmpMC = (TH1D*)canvas->GetPrimitive("tmpMC"); //The name of the MC histogram from DrawDataMC, which is drawn first
  canvas->SaveAs((std::string("./")+name+std::string("Log.png")).c_str());
}



void DrawFullAndEmptyDataMC(TH1D* DataInFull, TH1D* MCInFull, TH1D* DataInEmpty, TH1D* MCInEmpty, std::string name, double mcpotfull, double datapotfull, double mcpotempty, double datapotempty, bool ZoomTargets = false)
{
        TCanvas *canvas = new TCanvas("canvas", "canvas", 1500, 1000);
        //TH1D *DataHist = (TH1D*)DataIn->Clone();
        //TH1D *MCHist = (TH1D*)MCIn->Clone();
        MCInFull->SetTitle("MC");
        DataInFull->SetTitle("Data");
        MCInEmpty->SetTitle("MC");
        DataInEmpty->SetTitle("Data");
        double xaxismin = 4280;
        double xaxismax = ZoomTargets ? 6000 : 8400;
        DataInFull->GetXaxis()->SetRangeUser(xaxismin, xaxismax);
        DataInEmpty->GetXaxis()->SetRangeUser(xaxismin, xaxismax);
        MCInFull->GetXaxis()->SetRangeUser(xaxismin, xaxismax);
        MCInEmpty->GetXaxis()->SetRangeUser(xaxismin, xaxismax);
        DataInFull->SetXTitle("Z Position (mm)");
        DataInEmpty->SetXTitle("Z Position (mm)");
        MCInFull->SetXTitle("Z Position (mm)");
        MCInEmpty->SetXTitle("Z Position (mm)");
        PlotUtils::MnvPlotter plotter;
        plotter.ApplyStyle(PlotUtils::kCCQENuStyle);

        //Draw filled linear plot
        plotter.DrawDataMC(DataInFull, MCInFull, datapotfull/mcpotfull, "TL", true);
        double ymin = gPad->GetFrame()->GetY1();
        double ymax = gPad->GetFrame()->GetY2();
        PlotRegions(ymin, ymax, 0, 0);
        canvas->SaveAs((std::string("./")+name+std::string("Filled.png")).c_str());
        
        //Draw empty linear plot
        canvas->Clear();
        plotter.DrawDataMC(DataInEmpty, MCInEmpty, datapotempty/mcpotempty, "TL", true);
        ymin = gPad->GetFrame()->GetY1();
        ymax = gPad->GetFrame()->GetY2();
        PlotRegions(ymin, ymax, 0, 0);
        canvas->SaveAs((std::string("./")+name+std::string("Empty.png")).c_str());


        //Draw filled ratio plot
        canvas->Clear();
        plotter.DrawDataMCRatio(DataInFull, MCInFull, datapotfull/mcpotfull, true, 0, 2.0, "Data MC Ratio");
        ymin = gPad->GetFrame()->GetY1();
        ymax = gPad->GetFrame()->GetY2();
        PlotRegions(ymin, ymax, 0, 0);
        canvas->SaveAs((std::string("./")+name+std::string("FilledRatio.png")).c_str());


        //Draw empty ratio plot
        canvas->Clear();
        plotter.DrawDataMCRatio(DataInEmpty, MCInEmpty, datapotempty/mcpotempty, true, 0.0, 2.0, "Data MC Ratio");
        ymin = gPad->GetFrame()->GetY1();
        ymax = gPad->GetFrame()->GetY2();
        PlotRegions(ymin, ymax, 0, 0);
        canvas->SaveAs((std::string("./")+name+std::string("FilledRatio.png")).c_str());

        //Draw filled log plot
        canvas->Clear();
        canvas->SetLogy();
        DataInFull->SetMinimum(1);
        MCInFull->SetMinimum(1);
        plotter.DrawDataMC(DataInFull, MCInFull, datapotfull/mcpotfull, "TL", true);
        ymin = gPad->GetFrame()->GetY1();
        ymax = gPad->GetFrame()->GetY2();
        PlotRegions(std::pow(10,ymin), std::pow(10,ymax), 0, 0);
        canvas->SaveAs((std::string("./")+name+std::string("FilledLog.png")).c_str());

        //Draw empty log plot
        canvas->Clear();
        canvas->SetLogy();
        DataInEmpty->SetMinimum(1);
        MCInEmpty->SetMinimum(1);
        plotter.DrawDataMC(DataInEmpty, MCInEmpty, datapotempty/mcpotempty, "TL", true);
        ymin = gPad->GetFrame()->GetY1();
        ymax = gPad->GetFrame()->GetY2();
        PlotRegions(std::pow(10,ymin), std::pow(10,ymax), 0, 0);
        canvas->SaveAs((std::string("./")+name+std::string("EmptyLog.png")).c_str());


        //Drawing split plots, with both filled and empty in the same plot
        //Preparing canvas
        canvas->Clear();
        canvas->SetLogy(false); //Unsetting log axis
        TPad *p1 = new TPad("p1","top",0,0.45,1,0.9); // xlow,yloy, xup,yup
        TPad *p2 = new TPad("p2","bottom",0,0,1,0.45);
        auto tex3 = new TLatex(.5,.95,"Title2");
        tex3->SetTextAlign(22);
        tex3->SetTextSize(0.025);
        std::string title =" Z vertex distributions - Calibrated (and weighted) (top) and uncalibrated (bottom)";
        tex3->DrawLatexNDC(.5,.95, title.c_str());
        //Draw linear plots
        p1->Draw();
        p1->cd();
        //Draw filled linear plot
        plotter.DrawDataMC(DataInFull, MCInFull, datapotfull/mcpotfull, "TL", true);
        ymin = gPad->GetFrame()->GetY1();
        ymax = gPad->GetFrame()->GetY2();
        PlotRegions(ymin, ymax, 0, 0);
        canvas->cd();
        p2->Draw();
        p2->cd();
        //Draw empty linear plot
        plotter.DrawDataMC(DataInEmpty, MCInEmpty, datapotempty/mcpotempty, "TL", true);
        ymin = gPad->GetFrame()->GetY1();
        ymax = gPad->GetFrame()->GetY2();
        PlotRegions(ymin, ymax, 0, 0);
        canvas->cd();
        canvas->SaveAs((std::string("./")+name+std::string("Split.png")).c_str());


        //Draw ratio plots
        canvas->Clear();
        p1 = new TPad("p1","top",0,0.45,1,0.9); // xlow,yloy, xup,yup
        p2 = new TPad("p2","bottom",0,0,1,0.45);

        tex3->DrawLatexNDC(.5,.95, title.c_str());
        //Draw filled ratio plot
        p1->Draw();
        p1->cd();
        plotter.DrawDataMCRatio(DataInFull, MCInFull, datapotfull/mcpotfull, true, 0.0, 2.0, "Data MC Ratio");
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
        plotter.DrawDataMCRatio(DataInEmpty, MCInEmpty, datapotempty/mcpotempty, true, 0.0, 2.0, "Data MC Ratio");
        ymin = gPad->GetFrame()->GetY1();
        ymax = gPad->GetFrame()->GetY2();
        PlotRegions(ymin, ymax, 0, 0);
        auto tex2 = new TLatex(0.5,0.5,"Title1");
        tex2->SetTextAlign(22);
        tex2->SetTextSize(0.25);
        tex2->Draw("SAME");
        canvas->cd();

        canvas->SaveAs((std::string("./")+name+std::string("RatioSplit.png")).c_str());




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
        plotter.DrawDataMC(DataInFull, MCInFull, datapotfull/mcpotfull, "TL", true);
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
        plotter.DrawDataMC(DataInEmpty, MCInEmpty, datapotempty/mcpotempty, "TL", true);
        ymin = gPad->GetFrame()->GetY1();
        ymax = gPad->GetFrame()->GetY2();
        std::cout << "Empty ymin: " << ymin << " ymax: " << ymax << std::endl;
        PlotRegions(std::pow(10,ymin), std::pow(10,ymax), 0, 0);
        canvas->cd();
        canvas->SaveAs((std::string("./")+name+std::string("LogSplit.png")).c_str());
}

//==============================================================================
// End - Plotting Functions
//==============================================================================

  int getTgtCodeFromCalibrator(CVUniverse *universe, int ord, bool useExtendedTarget )
  {
    double vtx_x, vtx_y, vtx_z;
    //First step, what target does the tuple say it's in
    int ANNTgtCode = universe->GetCalibratedANNTargetCode(ord);
    //std::cout<<"ANNTgtCode: " << ANNTgtCode << " LatexName: " << universe->LatexName()<< " sigma: " << universe->GetSigma() << std::endl;
    if (ANNTgtCode > 0) return ANNTgtCode; //Maybe do a more explicit check for if it's a valid target code?

    //Water is a special case, due to (I think?) a bug in the MAT tuples. To figure out if an event happened in water we cant use the regular tuple branch, instead we check if the ANN reconstruction put in the water segment (36)
    //if (universe->GetANNSegment() == 36 ) return 6000;

    ROOT::Math::XYZVector Vtx = universe->GetCalibratedANNVertex();
    vtx_x = Vtx.X();
    vtx_y = Vtx.Y();
    vtx_z = Vtx.Z();
    PlotUtils::TargetUtils tgtUtil;
    if (tgtUtil.InWaterTargetVolMC( vtx_x, vtx_y, vtx_z, 850. )) return 6000;
    return -999;
  }

double GetANNProbFromCalibrator(CVUniverse *universe, int ord) {
  SMSCalibrator::CalibratedProbs cal = universe->GetANNSegmentsCalibrated();
  double excludedProb = 0;
  for (int i = 0; i < ord; i++)
  {
    excludedProb += cal.probs[i];
  }
  
  double prob = cal.probs[ord]/(1-excludedProb);
  //std::cout<<"prob: " << prob << " sigma: " << T::m_nsigma << std::endl;
  return prob;
}


//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFillMC(
    PlotUtils::ChainWrapper* chain,
    std::map<std::string, std::vector<CVUniverse*> > error_bands,
    PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts,
    PlotUtils::Model<CVUniverse, MichelEvent>& model,
    SMSCalibrator& calibrator
    )
{
  //calibrator.cacheEvt = -1;
  std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
  assert(!error_bands["cv"].empty() && "\"cv\" error band is empty!  Can't set Model weight.");
  auto& cvUniv = error_bands["cv"].front();

  std::cout << "Starting MC reco loop...\n";

  //const int nEntries = chain->GetEntries();
  const int nEntries = 10000;
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
      std::cout << i << " / " << nEntries << ". Time Elapsed: " << selapsed <<" Est Remaining: " << sleft <<"          \r" <<std::flush;
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
        //std::cout<<"Here1\n";
        universe->SetCalibrator(&calibrator);
        //std::cout<<"Here2\n";
        bool hasMLPred=universe->hasMLPred();
        ROOT::Math::XYZVector ANNVtx = universe->GetANNVertex();
        ROOT::Math::XYZTVector TrackBasedVtx = universe->GetVertex();
        ROOT::Math::XYZTVector TrueVtx = universe->GetTrueVertex();

        double WeightedANNZ = universe->GetANNSegmentsZPosWeighted();
        double ANNZResidual = ANNVtx.Z() - TrueVtx.Z();
        double ANNWeightedZResidual = WeightedANNZ - TrueVtx.Z();
        double ANNWeightedMinusUnweighted = WeightedANNZ - ANNVtx.Z();

        PlotUtils::TargetUtils* m_TargetUtils=new PlotUtils::TargetUtils();
        const double weight = model.GetWeight(*universe, myevent); //Only calculate the per-universe weight for events that will actually use it.
        double ANNProb = universe->GetANNProb(false);
        double erecoil = universe->GetRecoilE()/pow(10,3);
        double curvsig = 1/universe->GetMuonQPErr();
        double pmu = universe->GetPmu()/1000;
        int multiplicity = universe->GetMultiplicity();
        double ANNTopProbCalibrated = universe->GetANNProb();
        double ANNProb1Calibrated = universe->GetCalibratedANNProbAdjusted(0);
        double ANNProb2Calibrated = universe->GetCalibratedANNProbAdjusted(1);
        double ANNProb3Calibrated = universe->GetCalibratedANNProbAdjusted(2);

        double ANNProbFromCalibratorRaw1 = universe->GetCalibratedANNProb(0);
        double ANNProbFromCalibratorRaw2 = universe->GetCalibratedANNProb(1);
        double ANNProbFromCalibratorRaw3 = universe->GetCalibratedANNProb(2);
        const bool isSignal = michelcuts.isSignal(*universe, weight);

        // This is where you would Access/create a Michel

        if (isSignal && hasMLPred)
        {
          int tgtCode = util::getTgtCode(universe, true, true ); //The second bool (true) turns on the extended target

          if ((tgtCode == 2026) || (tgtCode == 3026) || (tgtCode == 5026))
          {
            ANNCalibratedConfidence_Iron_Signal->Fill(ANNTopProbCalibrated, weight);
            ANNCalibratedConfidence1_Iron_Signal->Fill(ANNProb1Calibrated, weight);
            ANNCalibratedConfidence2_Iron_Signal->Fill(ANNProb2Calibrated, weight);
            ANNCalibratedConfidence3_Iron_Signal->Fill(ANNProb3Calibrated, weight);

            ANNCalibratedConfidenceRatio2_Iron_Signal->Fill(ANNProbFromCalibratorRaw2/ANNProbFromCalibratorRaw1, weight);
            ANNCalibratedConfidenceRatio3_Iron_Signal->Fill(ANNProbFromCalibratorRaw3/ANNProbFromCalibratorRaw2, weight);

          }
          if ((tgtCode == 2082) || (tgtCode == 3082) || (tgtCode == 4082) || (tgtCode == 5082))
          {
            ANNCalibratedConfidence_Lead_Signal->Fill(ANNTopProbCalibrated, weight);
            ANNCalibratedConfidence1_Lead_Signal->Fill(ANNProb1Calibrated, weight);
            ANNCalibratedConfidence2_Lead_Signal->Fill(ANNProb2Calibrated, weight);
            ANNCalibratedConfidence3_Lead_Signal->Fill(ANNProb3Calibrated, weight);

            ANNCalibratedConfidenceRatio2_Lead_Signal->Fill(ANNProbFromCalibratorRaw2/ANNProbFromCalibratorRaw1, weight);
            ANNCalibratedConfidenceRatio3_Lead_Signal->Fill(ANNProbFromCalibratorRaw3/ANNProbFromCalibratorRaw2, weight);
          }
          if ((tgtCode == 3006))
          {
            ANNCalibratedConfidence_Carbon_Signal->Fill(ANNTopProbCalibrated, weight);
            ANNCalibratedConfidence1_Carbon_Signal->Fill(ANNProb1Calibrated, weight);
            ANNCalibratedConfidence2_Carbon_Signal->Fill(ANNProb2Calibrated, weight);
            ANNCalibratedConfidence3_Carbon_Signal->Fill(ANNProb3Calibrated, weight);

            ANNCalibratedConfidenceRatio2_Carbon_Signal->Fill(ANNProbFromCalibratorRaw2/ANNProbFromCalibratorRaw1, weight);
            ANNCalibratedConfidenceRatio3_Carbon_Signal->Fill(ANNProbFromCalibratorRaw3/ANNProbFromCalibratorRaw2, weight);
          }
          if ((tgtCode == 6000))
          {
            ANNCalibratedConfidence_Water_Signal->Fill(ANNTopProbCalibrated, weight);
            ANNCalibratedConfidence1_Water_Signal->Fill(ANNProb1Calibrated, weight);
            ANNCalibratedConfidence2_Water_Signal->Fill(ANNProb2Calibrated, weight);
            ANNCalibratedConfidence3_Water_Signal->Fill(ANNProb3Calibrated, weight);

            ANNCalibratedConfidenceRatio2_Water_Signal->Fill(ANNProbFromCalibratorRaw2/ANNProbFromCalibratorRaw1, weight);
            ANNCalibratedConfidenceRatio3_Water_Signal->Fill(ANNProbFromCalibratorRaw3/ANNProbFromCalibratorRaw2, weight);
          }
        }

        //weight is ignored in isMCSelected() for all but the CV Universe.
        if (!michelcuts.isMCSelected(*universe, myevent, weight).all()) continue;   //All cuts except ANN confidence cut

        if(hasMLPred)
        {
          ANNZ_vs_Conf_vs_Erecoil_MC->Fill( ANNVtx.Z(), ANNProb, erecoil, weight);
          ANNWeightedZ_vs_Conf_vs_Erecoil_MC->Fill( WeightedANNZ, ANNProb, erecoil, weight);

          //Calibrated
          //std::cout<<"Here3\n";
          double calibratedZ = universe->GetANNZPos();
          double calibratedWeightedZ = universe->GetANNZPosWeightedCalibrator();
          //std::cout<<"calibratedWeightedZ: " << calibratedWeightedZ <<std::endl;
          //std::cout<<"Here3a\n";
          ANNZCalibrated_vs_ConfCalibrated_vs_Erecoil_MC->Fill( calibratedZ, ANNTopProbCalibrated, erecoil, weight);
          //std::cout<<"Here3b\n";
          ANNWeightedZCalibrated_vs_ConfCalibrated_vs_Erecoil_MC->Fill( calibratedWeightedZ, ANNTopProbCalibrated, erecoil, weight);
          //std::cout<<"Here4\n";
          ANNZCalibratedResidual_vs_ConfCalibrated_vs_Erecoil->Fill( (calibratedZ - TrueVtx.Z()), ANNTopProbCalibrated, weight);
          ANNWeightedZCalibratedResidual_vs_ConfCalibrated_vs_Erecoil->Fill( (calibratedWeightedZ - TrueVtx.Z()), ANNTopProbCalibrated, weight);
          ANNWeightedZCalibratedMinusANNZ_MC->Fill( (universe->GetANNZPosWeightedCalibrator()-ANNVtx.Z()), weight);
          ANNZCalibratedMinusANNZ_MC->Fill( (calibratedZ-ANNVtx.Z()), weight);
          ANNZCalibrated_MC->Fill( calibratedZ, weight);
          ANNZCalibratedWeighted_MC->Fill( (universe->GetANNZPosWeightedCalibrator()), weight);
          ANNZ_MC->Fill( (ANNVtx.Z()), weight);


          if (m_TargetUtils->InIron2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNConf_vs_Ehad_Tgt2Iron_MC->Fill( ANNProb, erecoil, weight);
          if (m_TargetUtils->InLead2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNConf_vs_Ehad_Tgt2Lead_MC->Fill( ANNProb, erecoil, weight);
          if (m_TargetUtils->InIron3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNConf_vs_Ehad_Tgt3Iron_MC->Fill( ANNProb, erecoil, weight);
          if (m_TargetUtils->InLead3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNConf_vs_Ehad_Tgt3Lead_MC->Fill( ANNProb, erecoil, weight);
          if (m_TargetUtils->InCarbon3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNConf_vs_Ehad_Tgt3Carbon_MC->Fill( ANNProb, erecoil, weight);
          if (m_TargetUtils->InLead4VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNConf_vs_Ehad_Tgt4Lead_MC->Fill( ANNProb, erecoil, weight);
          if (m_TargetUtils->InIron5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNConf_vs_Ehad_Tgt5Iron_MC->Fill( ANNProb, erecoil, weight);
          if (m_TargetUtils->InLead5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNConf_vs_Ehad_Tgt5Lead_MC->Fill( ANNProb, erecoil, weight);
          if (m_TargetUtils->InWaterTargetVolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNConf_vs_Ehad_Water_MC->Fill( ANNProb, erecoil, weight);
          if (m_TargetUtils->InTracker(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNConf_vs_Ehad_Tracker_MC->Fill( ANNProb, erecoil, weight);
          
          //For events misreconstructed in each target, what's the probability distribution
          if (m_TargetUtils->InIron2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true) && !m_TargetUtils->InIron2VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true)) ANNConf_MisrecoInTgt2Iron->Fill( ANNProb);
          if (m_TargetUtils->InLead2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true) && !m_TargetUtils->InLead2VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true)) ANNConf_MisrecoInTgt2Lead->Fill( ANNProb);
          if (m_TargetUtils->InIron3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true) && !m_TargetUtils->InIron3VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true)) ANNConf_MisrecoInTgt3Iron->Fill( ANNProb);
          if (m_TargetUtils->InLead3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true) && !m_TargetUtils->InLead3VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true)) ANNConf_MisrecoInTgt3Lead->Fill( ANNProb);
          if (m_TargetUtils->InCarbon3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true) && !m_TargetUtils->InCarbon3VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true)) ANNConf_MisrecoInTgt3Carbon->Fill( ANNProb);
          if (m_TargetUtils->InLead4VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()) && !m_TargetUtils->InLead4VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z())) ANNConf_MisrecoInTgt4Lead->Fill( ANNProb);
          if (m_TargetUtils->InIron5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true) && !m_TargetUtils->InIron5VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true)) ANNConf_MisrecoInTgt5Iron->Fill( ANNProb);
          if (m_TargetUtils->InLead5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true) && !m_TargetUtils->InLead5VolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), 850, true)) ANNConf_MisrecoInTgt5Lead->Fill( ANNProb);
          if (m_TargetUtils->InWaterTargetVolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()) && !m_TargetUtils->InWaterTargetVolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z())) ANNConf_MisrecoInWater->Fill( ANNProb);


          ANNZResidual_vs_Conf_vs_Erecoil->Fill( ANNZResidual, ANNProb, weight);
          if (m_TargetUtils->InIron2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZResidual_vs_Conf_vs_Erecoil_Tgt2Iron->Fill( ANNZResidual, ANNProb, weight);
          if (m_TargetUtils->InLead2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZResidual_vs_Conf_vs_Erecoil_Tgt2Lead->Fill( ANNZResidual, ANNProb, weight);
          if (m_TargetUtils->InIron3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZResidual_vs_Conf_vs_Erecoil_Tgt3Iron->Fill( ANNZResidual, ANNProb, weight);
          if (m_TargetUtils->InLead3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZResidual_vs_Conf_vs_Erecoil_Tgt3Lead->Fill( ANNZResidual, ANNProb, weight);
          if (m_TargetUtils->InCarbon3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZResidual_vs_Conf_vs_Erecoil_Tgt3Carbon->Fill( ANNZResidual, ANNProb, weight);
          if (m_TargetUtils->InLead4VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNZResidual_vs_Conf_vs_Erecoil_Tgt4Lead->Fill( ANNZResidual, ANNProb, weight);
          if (m_TargetUtils->InIron5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZResidual_vs_Conf_vs_Erecoil_Tgt5Iron->Fill( ANNZResidual, ANNProb, weight);
          if (m_TargetUtils->InLead5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNZResidual_vs_Conf_vs_Erecoil_Tgt5Lead->Fill( ANNZResidual, ANNProb, weight);
          if (m_TargetUtils->InWaterTargetVolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNZResidual_vs_Conf_vs_Erecoil_Water->Fill( ANNZResidual, ANNProb, weight);
          if (m_TargetUtils->InTracker(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNZResidual_vs_Conf_vs_Erecoil_Tracker->Fill( ANNZResidual, ANNProb, weight);

          ANNWeightedZResidual_vs_Conf_vs_Erecoil->Fill( (ANNWeightedZResidual), ANNProb, weight);
          if (m_TargetUtils->InIron2VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt2Iron->Fill( (ANNWeightedZResidual), ANNProb, weight);
          if (m_TargetUtils->InLead2VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt2Lead->Fill( (ANNWeightedZResidual), ANNProb, weight);
          if (m_TargetUtils->InIron3VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt3Iron->Fill( (ANNWeightedZResidual), ANNProb, weight);
          if (m_TargetUtils->InLead3VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt3Lead->Fill( (ANNWeightedZResidual), ANNProb, weight);
          if (m_TargetUtils->InCarbon3VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt3Carbon->Fill( (ANNWeightedZResidual), ANNProb, weight);
          if (m_TargetUtils->InLead4VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ))  ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt4Lead->Fill( (ANNWeightedZResidual), ANNProb, weight);
          if (m_TargetUtils->InIron5VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt5Iron->Fill( (ANNWeightedZResidual), ANNProb, weight);
          if (m_TargetUtils->InLead5VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt5Lead->Fill( (ANNWeightedZResidual), ANNProb, weight);
          if (m_TargetUtils->InWaterTargetVolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ))  ANNWeightedZResidual_vs_Conf_vs_Erecoil_Water->Fill( (ANNWeightedZResidual), ANNProb, weight);
          if (m_TargetUtils->InTracker(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ))  ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tracker->Fill( (ANNWeightedZResidual), ANNProb, weight);

          if (ANNConfCut->passesCut(*universe, myevent, weight)) //ANN confidence cut
          {
            ANNWeightedZ_minus_UnweightedZMC->Fill( ANNWeightedMinusUnweighted, weight);
            if (m_TargetUtils->InIron2VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZ_minus_UnweightedZTgt2IronMC->Fill( ANNWeightedMinusUnweighted, weight);
            if (m_TargetUtils->InLead2VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZ_minus_UnweightedZTgt2LeadMC->Fill( ANNWeightedMinusUnweighted, weight);
            if (m_TargetUtils->InIron3VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZ_minus_UnweightedZTgt3IronMC->Fill( ANNWeightedMinusUnweighted, weight);
            if (m_TargetUtils->InLead3VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZ_minus_UnweightedZTgt3LeadMC->Fill( ANNWeightedMinusUnweighted, weight);
            if (m_TargetUtils->InCarbon3VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZ_minus_UnweightedZTgt3CarbonMC->Fill( ANNWeightedMinusUnweighted, weight);
            if (m_TargetUtils->InLead4VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ))  ANNWeightedZ_minus_UnweightedZTgt4LeadMC->Fill( ANNWeightedMinusUnweighted, weight);
            if (m_TargetUtils->InIron5VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZ_minus_UnweightedZTgt5IronMC->Fill( ANNWeightedMinusUnweighted, weight);
            if (m_TargetUtils->InLead5VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZ_minus_UnweightedZTgt5LeadMC->Fill( ANNWeightedMinusUnweighted, weight);
            if (m_TargetUtils->InWaterTargetVolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ))  ANNWeightedZ_minus_UnweightedZWaterMC->Fill( ANNWeightedMinusUnweighted, weight);
            if (m_TargetUtils->InTracker(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ))  ANNWeightedZ_minus_UnweightedZTrackerMC->Fill( ANNWeightedMinusUnweighted, weight);
          }
          ANNConfidence_MC->Fill(ANNProb, weight);
          ANNCalibratedConfidence_MC->Fill(ANNTopProbCalibrated, weight);

          int tgtCode = util::getTgtCode(universe, false, true ); //The second bool (true) turns on the extended target

          //Weighted Average
          int ANNtgtCode = util::getTgtCode(universe, false, true ); //The second bool (true) turns on the extended target
          if ((ANNtgtCode == 2026) || (ANNtgtCode == 3026) || (ANNtgtCode == 5026))
          {
            ANNCalibratedConfidence_Iron_Selected->Fill(ANNTopProbCalibrated, weight);
          }
          if ((ANNtgtCode == 2082) || (ANNtgtCode == 3082) || (ANNtgtCode == 4082) || (ANNtgtCode == 5082))
          {
            ANNCalibratedConfidence_Lead_Selected->Fill(ANNTopProbCalibrated, weight);
          }
          if ((ANNtgtCode == 3006))
          {
            ANNCalibratedConfidence_Carbon_Selected->Fill(ANNTopProbCalibrated, weight);
          }
          if ((ANNtgtCode == 6000))
          {
            ANNCalibratedConfidence_Water_Selected->Fill(ANNTopProbCalibrated, weight);
          }
          
          //Top plane
          ANNtgtCode =   getTgtCodeFromCalibrator(universe, 0, true ); //The bool (true) turns on the extended target
          //std::cout<<"Top plane target code: " << ANNtgtCode << std::endl;
          if ((ANNtgtCode == 2026) || (ANNtgtCode == 3026) || (ANNtgtCode == 5026))
          {
            ANNCalibratedConfidence1_Iron_Selected->Fill(ANNProb1Calibrated, weight);
          }
          if ((ANNtgtCode == 2082) || (ANNtgtCode == 3082) || (ANNtgtCode == 4082) || (ANNtgtCode == 5082))
          {
            ANNCalibratedConfidence1_Lead_Selected->Fill(ANNProb1Calibrated, weight);
          }
          if ((ANNtgtCode == 3006))
          {
            ANNCalibratedConfidence1_Carbon_Selected->Fill(ANNProb1Calibrated, weight);
          }
          if ((ANNtgtCode == 6000))
          {
            ANNCalibratedConfidence1_Water_Selected->Fill(ANNProb1Calibrated, weight);
          }
          //2nd plane
          ANNtgtCode =   getTgtCodeFromCalibrator(universe, 1, true ); //The bool (true) turns on the extended target
          //std::cout<<"2nd plane target code: " << ANNtgtCode << std::endl;
          if ((ANNtgtCode == 2026) || (ANNtgtCode == 3026) || (ANNtgtCode == 5026))
          {
            ANNCalibratedConfidence2_Iron_Selected->Fill(ANNProb2Calibrated, weight);
            ANNCalibratedConfidenceRatio2_Iron_Selected->Fill(ANNProbFromCalibratorRaw2/ANNProbFromCalibratorRaw1, weight);
          }
          if ((ANNtgtCode == 2082) || (ANNtgtCode == 3082) || (ANNtgtCode == 4082) || (ANNtgtCode == 5082))
          {
            ANNCalibratedConfidence2_Lead_Selected->Fill(ANNProb2Calibrated, weight);
            ANNCalibratedConfidenceRatio2_Lead_Selected->Fill(ANNProbFromCalibratorRaw2/ANNProbFromCalibratorRaw1, weight);
          }
          if ((ANNtgtCode == 3006))
          {
            ANNCalibratedConfidence2_Carbon_Selected->Fill(ANNProb2Calibrated, weight);
            ANNCalibratedConfidenceRatio2_Carbon_Selected->Fill(ANNProbFromCalibratorRaw2/ANNProbFromCalibratorRaw1, weight);
          }
          if ((ANNtgtCode == 6000))
          {
            ANNCalibratedConfidence2_Water_Selected->Fill(ANNProb2Calibrated, weight);
            ANNCalibratedConfidenceRatio2_Water_Selected->Fill(ANNProbFromCalibratorRaw2/ANNProbFromCalibratorRaw1, weight);
          }
          //3rd plane
          ANNtgtCode =   getTgtCodeFromCalibrator(universe, 2, true ); //The bool (true) turns on the extended target
          //std::cout<<"3rd plane target code: " << ANNtgtCode << std::endl;
          if ((ANNtgtCode == 2026) || (ANNtgtCode == 3026) || (ANNtgtCode == 5026))
          {
            ANNCalibratedConfidence3_Iron_Selected->Fill(ANNProb3Calibrated, weight);
            ANNCalibratedConfidenceRatio3_Iron_Selected->Fill(ANNProbFromCalibratorRaw3/ANNProbFromCalibratorRaw2, weight);
          }
          if ((ANNtgtCode == 2082) || (ANNtgtCode == 3082) || (ANNtgtCode == 4082) || (ANNtgtCode == 5082))
          {
            ANNCalibratedConfidence3_Lead_Selected->Fill(ANNProb3Calibrated, weight);
            ANNCalibratedConfidenceRatio3_Lead_Selected->Fill(ANNProbFromCalibratorRaw3/ANNProbFromCalibratorRaw2, weight);
          }
          if ((ANNtgtCode == 3006))
          {
            ANNCalibratedConfidence3_Carbon_Selected->Fill(ANNProb3Calibrated, weight);
            ANNCalibratedConfidenceRatio3_Carbon_Selected->Fill(ANNProbFromCalibratorRaw3/ANNProbFromCalibratorRaw2, weight);
          }
          if ((ANNtgtCode == 6000))
          {
            ANNCalibratedConfidence3_Water_Selected->Fill(ANNProb3Calibrated, weight);
            ANNCalibratedConfidenceRatio3_Water_Selected->Fill(ANNProbFromCalibratorRaw3/ANNProbFromCalibratorRaw2, weight);
          }
        }

        //if (!ANNConfCut->passesCut(*universe, myevent, weight);) continue;   //ANN confidence cut

        //Performing vtx validation check Deborah suggested
        //double batchPOT = cvUniv->GetBatchPOT();
        //double efficiency = 0.5563 - (0.01353*batchPOT); //Based on MINERvA-doc-21436
        //weight/=efficiency; //Not doing efficiency correction right now
        //Hadron Energy Spectrum plot
        //int multiplicity = cvUniv->GetMultiplicity();

        //-----------------------------------------------------------
        // Making a note of events that may be interesting to view in Arachne
        //-----------------------------------------------------------
        if ( std::abs(TrueVtx.Z() - ANNVtx.Z()) > 2000 && ANNProb>0.4 && verbose)
        {
          std::cout<<"Difference between true and ann reco z larger than 2000: \tANN Segment: " << universe->GetANNSegment() << "\tTrue Segment: " << universe->GetTruthSegment() << " mc_run: " << universe->GetInt("mc_run") << " mc_subrun: " << universe->GetInt("mc_subrun")  << " mc_nthEvtInSpill: " << universe->GetInt("mc_nthEvtInSpill") << " mc_nthEvtInFile: " << universe->GetInt("mc_nthEvtInFile") << " ev_global_gate: " << universe->GetInt("ev_global_gate")  << std::endl;
          std::string ArachneLink = "https://mnvevdgpvm02.fnal.gov/Arachne/?det=SIM_minerva&recoVer=v22r1p1&run="+std::to_string(universe->GetInt("mc_run"))+"&subrun="+std::to_string(universe->GetInt("mc_subrun"))+"&gate="+std::to_string(universe->GetInt("mc_nthEvtInFile")+1)+"&slice=-1";
          std::cout<<"Arachne Link: " << ArachneLink << std::endl;
        }
        ANNVerticesConfusion_ByZPos->Fill(TrueVtx.Z(), ANNVtx.Z(), weight);
        int intType = universe->GetInteractionType();
        if (intType == 1)
        {
          ConfusionMatrixQE->Fill(TrueVtx.Z(), ANNVtx.Z(), weight);
          ANNConfidenceQE_MC->Fill(ANNProb, weight);
        }
        else if (intType == 2)
        {
          ConfusionMatrixRES->Fill(TrueVtx.Z(), ANNVtx.Z(), weight);
          ANNConfidenceRES_MC->Fill(ANNProb, weight);
        }
        else if (intType == 3)
        {
          ConfusionMatrixDIS->Fill(TrueVtx.Z(), ANNVtx.Z(), weight);
          ANNConfidenceDIS_MC->Fill(ANNProb, weight);
        }
        else if (intType == 4)
        {
          ConfusionMatrixCOHPI->Fill(TrueVtx.Z(), ANNVtx.Z(), weight);
          ANNConfidenceCOHPI_MC->Fill(ANNProb, weight);
        }
        else if (intType == 5)
        {
          ConfusionMatrixAMNUGAMMA->Fill(TrueVtx.Z(), ANNVtx.Z(), weight);
          ANNConfidenceAMNUGAMMA_MC->Fill(ANNProb, weight);
        }
        else if (intType == 6)
        {
          ConfusionMatrixIMD->Fill(TrueVtx.Z(), ANNVtx.Z(), weight);
          ANNConfidenceIMD_MC->Fill(ANNProb, weight);
        }
        else if (intType == 7)
        {
          ConfusionMatrixNUEEL->Fill(TrueVtx.Z(), ANNVtx.Z(), weight);
          ANNConfidenceNUEEL_MC->Fill(ANNProb, weight);
        }
        else if (intType == 8)
        {
          ConfusionMatrixMEC->Fill(TrueVtx.Z(), ANNVtx.Z(), weight);
          ANNConfidenceMEC_MC->Fill(ANNProb, weight);
        }
        else
        {
          ConfusionMatrixUnknownInt->Fill(TrueVtx.Z(), ANNVtx.Z(), weight);
          ANNConfidenceUnknownInt_MC->Fill(ANNProb, weight);
        }
        ANNConfidence_vs_NTracks_MC->Fill(ANNProb, multiplicity, weight);

        ConfusionMatrixMultiplicity->Fill(TrueVtx.Z(), ANNVtx.Z(), multiplicity, weight);
        ConfusionMatrixErecoil->Fill(TrueVtx.Z(), ANNVtx.Z(), erecoil, weight);
        ConfusionMatrixConf->Fill(TrueVtx.Z(), ANNVtx.Z(), ANNProb, weight);
        ConfusionMatrix->Fill(TrueVtx.Z(), ANNVtx.Z(),weight);
        if (isSignal) // If it is signal
        {
          if (hasMLPred)
          {

            int tgtCode = util::getTgtCode(universe, true, true ); //The second bool (true) turns on the extended target
            //Weighted Average
            int ANNtgtCode = util::getTgtCode(universe, false, true ); //The second bool (true) turns on the extended target
            if (ANNtgtCode == tgtCode) //If signal (tgtCode) and selected (ANNtgtCode) match
            {
              if ((tgtCode == 2026) || (tgtCode == 3026) || (tgtCode == 5026))
              {
                ANNCalibratedConfidence_Iron_SelectedSignal->Fill(ANNTopProbCalibrated, weight);
              }
              if ((tgtCode == 2082) || (tgtCode == 3082) || (tgtCode == 4082) || (tgtCode == 5082))
              {
                ANNCalibratedConfidence_Lead_SelectedSignal->Fill(ANNTopProbCalibrated, weight);
              }
              if ((tgtCode == 3006))
              {
                ANNCalibratedConfidence_Carbon_SelectedSignal->Fill(ANNTopProbCalibrated, weight);
              }
              if ((tgtCode == 6000))
              {
                ANNCalibratedConfidence_Water_SelectedSignal->Fill(ANNTopProbCalibrated, weight);
              }
            }
            //Top plane
            ANNtgtCode =   getTgtCodeFromCalibrator(universe, 0, true ); //The bool (true) turns on the extended target
            if (ANNtgtCode == tgtCode) //If signal (tgtCode) and selected (ANNtgtCode) match
            {
              if ((tgtCode == 2026) || (tgtCode == 3026) || (tgtCode == 5026))
              {
                ANNCalibratedConfidence1_Iron_SelectedSignal->Fill(ANNTopProbCalibrated, weight);
              }
              if ((tgtCode == 2082) || (tgtCode == 3082) || (tgtCode == 4082) || (tgtCode == 5082))
              {
                ANNCalibratedConfidence1_Lead_SelectedSignal->Fill(ANNTopProbCalibrated, weight);
              }
              if ((tgtCode == 3006))
              {
                ANNCalibratedConfidence1_Carbon_SelectedSignal->Fill(ANNTopProbCalibrated, weight);
              }
              if ((tgtCode == 6000))
              {
                ANNCalibratedConfidence1_Water_SelectedSignal->Fill(ANNTopProbCalibrated, weight);
              }
            }
            //2nd plane
            ANNtgtCode =   getTgtCodeFromCalibrator(universe, 1, true ); //The bool (true) turns on the extended target
            if (ANNtgtCode == tgtCode) //If signal (tgtCode) and selected (ANNtgtCode) match
            {
              if ((tgtCode == 2026) || (tgtCode == 3026) || (tgtCode == 5026))
              {
                ANNCalibratedConfidence2_Iron_SelectedSignal->Fill(ANNProb2Calibrated, weight);
                ANNCalibratedConfidenceRatio2_Iron_SelectedSignal->Fill(ANNProbFromCalibratorRaw2/ANNProbFromCalibratorRaw1, weight);
              }
              if ((tgtCode == 2082) || (tgtCode == 3082) || (tgtCode == 4082) || (tgtCode == 5082))
              {
                ANNCalibratedConfidence2_Lead_SelectedSignal->Fill(ANNProb2Calibrated, weight);
                ANNCalibratedConfidenceRatio2_Lead_SelectedSignal->Fill(ANNProbFromCalibratorRaw2/ANNProbFromCalibratorRaw1, weight);
              }
              if ((tgtCode == 3006))
              {
                ANNCalibratedConfidence2_Carbon_SelectedSignal->Fill(ANNProb2Calibrated, weight);
                ANNCalibratedConfidenceRatio2_Carbon_SelectedSignal->Fill(ANNProbFromCalibratorRaw2/ANNProbFromCalibratorRaw1, weight);
              }
              if ((tgtCode == 6000))
              {
                ANNCalibratedConfidence2_Water_SelectedSignal->Fill(ANNProb2Calibrated, weight);
                ANNCalibratedConfidenceRatio2_Water_SelectedSignal->Fill(ANNProbFromCalibratorRaw2/ANNProbFromCalibratorRaw1, weight);
              }
            }
            //3rd plane
            ANNtgtCode =   getTgtCodeFromCalibrator(universe, 2, true ); //The bool (true) turns on the extended target
            if (ANNtgtCode == tgtCode) //If signal (tgtCode) and selected (ANNtgtCode) match
            {
              if ((tgtCode == 2026) || (tgtCode == 3026) || (tgtCode == 5026))
              {
                ANNCalibratedConfidence3_Iron_SelectedSignal->Fill(ANNProb3Calibrated, weight);
                ANNCalibratedConfidenceRatio3_Iron_SelectedSignal->Fill(ANNProbFromCalibratorRaw3/ANNProbFromCalibratorRaw2, weight);
              }
              if ((tgtCode == 2082) || (tgtCode == 3082) || (tgtCode == 4082) || (tgtCode == 5082))
              {
                ANNCalibratedConfidence3_Lead_SelectedSignal->Fill(ANNProb3Calibrated, weight);
                ANNCalibratedConfidenceRatio3_Lead_SelectedSignal->Fill(ANNProbFromCalibratorRaw3/ANNProbFromCalibratorRaw2, weight);
              }
              if ((tgtCode == 3006))
              {
                ANNCalibratedConfidence3_Carbon_SelectedSignal->Fill(ANNProb3Calibrated, weight);
                ANNCalibratedConfidenceRatio3_Carbon_SelectedSignal->Fill(ANNProbFromCalibratorRaw3/ANNProbFromCalibratorRaw2, weight);
              }
              if ((tgtCode == 6000))
              {
                ANNCalibratedConfidence3_Water_SelectedSignal->Fill(ANNProb3Calibrated, weight);
                ANNCalibratedConfidenceRatio3_Water_SelectedSignal->Fill(ANNProbFromCalibratorRaw3/ANNProbFromCalibratorRaw2, weight);
              }
            }

          }
        }
      }
    }
  } //End entries loop
  std::cout << "Finished MC reco loop.\n"; 
}

void LoopAndFillData( PlotUtils::ChainWrapper* data, std::vector<CVUniverse*> data_band, PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts,
    SMSCalibrator& calibrator)
{
  //calibrator.cacheEvt = -1;
  std::cout << "Starting data loop...\n";
  std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
  const int nEntries = data->GetEntries();
  //const int nEntries = 10000;
  int  startNum = 0;
  int endNum = nEntries;
  if (nSubruns != 0 )
  {
    std::raise(SIGINT);
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
  std::cout<<" endNum: "<< endNum << std::endl;
  for (int i = startNum; i < endNum; ++i)
  {
    for (auto universe : data_band) {
      //std::cout<<"i: " << i << std::endl;
      universe->SetEntry(i);
      universe->SetCalibrator(&calibrator);
      if(i%1000==0)
      { 
        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
        double mselapsed = std::chrono::duration_cast<std::chrono::milliseconds> (now - startTime).count();
        double msperevent = mselapsed/i;
        int sleft = (nEntries-i)*msperevent/1000;
        double selapsed = ((int)mselapsed)/1000;
        std::cout << i << " / " << nEntries << ". Time Elapsed: " << selapsed <<" Est Remaining: " << sleft <<"          \r" <<std::flush;
      }
      MichelEvent myevent; 

      bool hasMLPred=universe->hasMLPred();
      //std::raise(SIGINT);
      ROOT::Math::XYZVector ANNVtx = universe->GetANNVertex();
      //std::cout<<"ANNVtx.Z1() : " << ANNVtx.Z()<< std::endl;
      ROOT::Math::XYZTVector TrackBasedVtx = universe->GetVertex();

      double WeightedANNZ = universe->GetANNSegmentsZPosWeighted();
      double ANNWeightedMinusUnweighted = WeightedANNZ - ANNVtx.Z();

      double ANNProb = universe->GetANNProb(false);
      double erecoil = universe->GetRecoilE()/pow(10,3);
      double curvsig = 1/universe->GetMuonQPErr();

      PlotUtils::TargetUtils* m_TargetUtils=new PlotUtils::TargetUtils();

      if (!michelcuts.isDataSelected(*universe, myevent).all()) continue;   //All cuts except ANN confidence cut

      if(hasMLPred)
      {
        ANNZ_vs_Conf_vs_Erecoil_Data->Fill( ANNVtx.Z(), ANNProb, erecoil);
        ANNWeightedZ_vs_Conf_vs_Erecoil_Data->Fill( WeightedANNZ, ANNProb, erecoil);

        //Calibrated
        double ANNTopProbCalibrated = universe->GetANNProb();
        double calibratedZ = universe->GetANNZPos();
        ANNZCalibrated_vs_ConfCalibrated_vs_Erecoil_Data->Fill( calibratedZ, ANNTopProbCalibrated, erecoil);
        ANNWeightedZCalibrated_vs_ConfCalibrated_vs_Erecoil_Data->Fill( universe->GetANNZPosWeightedCalibrator(), ANNTopProbCalibrated, erecoil);
        ANNWeightedZCalibratedMinusANNZ_Data->Fill( (universe->GetANNZPosWeightedCalibrator()-ANNVtx.Z()));
        ANNZCalibratedMinusANNZ_Data->Fill( (calibratedZ-ANNVtx.Z()));

        ANNZCalibrated_Data->Fill( calibratedZ);
        ANNZCalibratedWeighted_Data->Fill( (universe->GetANNZPosWeightedCalibrator()));
        ANNZ_Data->Fill( (ANNVtx.Z()));
        //std::cout<<"ANNVtx.Z2() : " << ANNVtx.Z()<< std::endl;



        if (m_TargetUtils->InIron2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNConf_vs_Ehad_Tgt2Iron_Data->Fill( ANNProb, erecoil);
        if (m_TargetUtils->InLead2VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNConf_vs_Ehad_Tgt2Lead_Data->Fill( ANNProb, erecoil);
        if (m_TargetUtils->InIron3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNConf_vs_Ehad_Tgt3Iron_Data->Fill( ANNProb, erecoil);
        if (m_TargetUtils->InLead3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNConf_vs_Ehad_Tgt3Lead_Data->Fill( ANNProb, erecoil);
        if (m_TargetUtils->InCarbon3VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNConf_vs_Ehad_Tgt3Carbon_Data->Fill( ANNProb, erecoil);
        if (m_TargetUtils->InLead4VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNConf_vs_Ehad_Tgt4Lead_Data->Fill( ANNProb, erecoil);
        if (m_TargetUtils->InIron5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNConf_vs_Ehad_Tgt5Iron_Data->Fill( ANNProb, erecoil);
        if (m_TargetUtils->InLead5VolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z(), 850, true))  ANNConf_vs_Ehad_Tgt5Lead_Data->Fill( ANNProb, erecoil);
        if (m_TargetUtils->InWaterTargetVolMC(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNConf_vs_Ehad_Water_Data->Fill( ANNProb, erecoil);
        if (m_TargetUtils->InTracker(ANNVtx.X(), ANNVtx.Y(), ANNVtx.Z()))  ANNConf_vs_Ehad_Tracker_Data->Fill( ANNProb, erecoil);

        if (ANNConfCut->passesCut(*universe, myevent)){ //ANN confidence cut
          ANNWeightedZ_minus_UnweightedZData->Fill( ANNWeightedMinusUnweighted);
          if (m_TargetUtils->InIron2VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZ_minus_UnweightedZTgt2IronData->Fill( ANNWeightedMinusUnweighted);
          if (m_TargetUtils->InLead2VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZ_minus_UnweightedZTgt2LeadData->Fill( ANNWeightedMinusUnweighted);
          if (m_TargetUtils->InIron3VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZ_minus_UnweightedZTgt3IronData->Fill( ANNWeightedMinusUnweighted);
          if (m_TargetUtils->InLead3VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZ_minus_UnweightedZTgt3LeadData->Fill( ANNWeightedMinusUnweighted);
          if (m_TargetUtils->InCarbon3VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZ_minus_UnweightedZTgt3CarbonData->Fill( ANNWeightedMinusUnweighted);
          if (m_TargetUtils->InLead4VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ))  ANNWeightedZ_minus_UnweightedZTgt4LeadData->Fill( ANNWeightedMinusUnweighted);
          if (m_TargetUtils->InIron5VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZ_minus_UnweightedZTgt5IronData->Fill( ANNWeightedMinusUnweighted);
          if (m_TargetUtils->InLead5VolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ, 850, true))  ANNWeightedZ_minus_UnweightedZTgt5LeadData->Fill( ANNWeightedMinusUnweighted);
          if (m_TargetUtils->InWaterTargetVolMC(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ))  ANNWeightedZ_minus_UnweightedZWaterData->Fill( ANNWeightedMinusUnweighted);
          if (m_TargetUtils->InTracker(ANNVtx.X(), ANNVtx.Y(), WeightedANNZ))  ANNWeightedZ_minus_UnweightedZTrackerData->Fill( ANNWeightedMinusUnweighted);
        }
        int multiplicity = universe->GetMultiplicity();
        std::cout<<"i: " << i << " ANNProb: " << ANNProb << " ANNVtx.Z(): " << ANNVtx.Z() <<std::endl;
        ANNConfidence_vs_NTracks_Data->Fill(ANNProb, multiplicity);
        ANNConfidence_Data->Fill(ANNProb);
        ANNCalibratedConfidence_Data->Fill(ANNTopProbCalibrated);
      }
      //if (!ANNConfCut->passesCut(*universe, myevent)) continue;   //ANN confidence cut
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
  if(argc > nArgsExpected + 1 || argc < nArgsExpected) //argc is the size of argv.  I check for number of arguments + 1 because
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
  //std::cout<<"Here\n";
  PlotUtils::MacroUtil options(reco_tree_name, mc_file_list, data_file_list, "minervame1A", true); //minervame1A is just a placeholder, it gets overwritted immediately below
  //std::cout<<"Here1\n";
  options.m_plist_string = util::GetPlaylist(*options.m_mc, true); //TODO: Put GetPlaylist into PlotUtils::MacroUtil
  //std::cout<<"Here2\n";

  // You're required to make some decisions
  PlotUtils::MinervaUniverse::SetNuEConstraint(true);
  PlotUtils::MinervaUniverse::SetPlaylist(options.m_plist_string); //TODO: Infer this from the files somehow?
  int nuoranu = nuOrAntiNuMode(options.m_plist_string);
  int nupdg;
  if (nuoranu == 1)
  {
    nupdg = 14;
  }
  else if (nuoranu == 2)
  {
    nupdg = -14;
  }
  PlotUtils::MinervaUniverse::SetAnalysisNuPDG(nupdg);
  //PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
  PlotUtils::MinervaUniverse::SetZExpansionFaReweight(false);


  int filledOrEmpty = util::filledOrEmptyMEPlaylist(options.m_plist_string); //filled = 1, empty = 2
  SMSCalibrator calibrator = util::loadCalibrator((nupdg == 14), (filledOrEmpty==1));


  std::vector<int> testSegments = {199, 200, 198};
  std::vector<double> testProbs = {0.906584, 0.0906537, 0.00170366};
  SMSCalibrator::CalibratedProbs cal = calibrator.Top3Calibrated(testSegments, testProbs, -1);
  std::cout<<"Test Calibrated Probability " << cal.segments[0]<< ":" << cal.probs[0]<< " " << cal.segments[1]<< ":" << cal.probs[1]<< " " << cal.segments[2]<< ":" << cal.probs[2]<< "\n";

  //return 0;
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
  int tuneA = 4;
  int tuneB = 3;
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
  
  // Make a map of systematic universes
  // Leave out systematics when making validation histograms
  const bool doSystematics = (getenv("MNV101_SKIP_SYST") == nullptr);
  if (!doSystematics)
  {
    std::cout << "Skipping systematics (except 1 flux universe) because environment variable MNV101_SKIP_SYST is set.\n";
    PlotUtils::MinervaUniverse::SetNFluxUniverses(2); // Necessary to get Flux integral later...  Doesn't work with just 1 flux universe though because _that_ triggers "spread errors".
  }

  std::map<std::string, std::vector<CVUniverse *>> error_bands;
  /* if (doSystematics)
    error_bands = GetStandardSystematics(options.m_mc);
  else
  {
    std::map<std::string, std::vector<CVUniverse *>> band_flux = PlotUtils::GetFluxSystematicsMap<CVUniverse>(options.m_mc, CVUniverse::GetNFluxUniverses());
    error_bands.insert(band_flux.begin(), band_flux.end()); // Necessary to get flux integral later...
  } */
  error_bands["cv"] = {new CVUniverse(options.m_mc)};
  std::map<std::string, std::vector<CVUniverse *>> truth_bands;
  //if (doSystematics) truth_bands = GetStandardSystematics(options.m_truth);
  truth_bands["cv"] = {new CVUniverse(options.m_truth)};

  CVUniverse *data_universe = new CVUniverse(options.m_data);
  std::vector<CVUniverse *> data_band = {data_universe};
  std::map<std::string, std::vector<CVUniverse *>> data_error_bands;
  data_error_bands["cv"] = data_band;

  // Loop entries and fill
  //try
  {
    std::cout << "Staring event loops\n";
    CVUniverse::SetTruth(false);
    LoopAndFillMC(options.m_mc, error_bands, mycuts, model, calibrator);
    CVUniverse::SetTruth(true);
    //LoopAndFillEffDenom(options.m_truth, truth_bands, mycuts, model); //Stay commented out
    options.PrintMacroConfiguration(argv[0]);
    std::cout << "MC cut summary:\n" << mycuts << "\n";
    mycuts.resetStats();
    CVUniverse::SetTruth(false);
    LoopAndFillData(options.m_data, data_band, mycuts, calibrator);
    std::cout << "Data cut summary:\n" << mycuts << "\n";

    std::string outFileName = "VertexValidations.root";
    if (nSubruns != 0 ) outFileName = "VertexValidations_n"+std::to_string(nProcess)+ ".root";
    TFile* OutDir = TFile::Open(outFileName.c_str(), "RECREATE");

    TruthZVertex_vs_Erecoil_MC->Write();

    ANNZ_vs_Conf_vs_Erecoil_Data->Write();
    ANNZ_vs_Conf_vs_Erecoil_MC->Write();
    ANNWeightedZ_vs_Conf_vs_Erecoil_Data->Write();
    ANNWeightedZ_vs_Conf_vs_Erecoil_MC->Write();

    ANNZCalibrated_vs_ConfCalibrated_vs_Erecoil_Data->Write();
    ANNZCalibrated_vs_ConfCalibrated_vs_Erecoil_MC->Write();
    ANNWeightedZCalibrated_vs_ConfCalibrated_vs_Erecoil_Data->Write();
    ANNWeightedZCalibrated_vs_ConfCalibrated_vs_Erecoil_MC->Write();

    ANNWeightedZCalibratedMinusANNZ_Data->Write();
    ANNWeightedZCalibratedMinusANNZ_MC->Write();

    ANNZCalibratedMinusANNZ_Data->Write();
    ANNZCalibratedMinusANNZ_MC->Write();

    ANNZCalibrated_Data->Write();
    ANNZCalibrated_MC->Write();
    ANNZCalibratedWeighted_Data->Write();
    ANNZCalibratedWeighted_MC->Write();
    ANNZ_Data->Write();
    ANNZ_MC->Write();

    ANNZCalibratedResidual_vs_ConfCalibrated_vs_Erecoil->Write();
    ANNZResidual_vs_Conf_vs_Erecoil->Write();
    ANNZResidual_vs_Conf_vs_Erecoil_Tracker->Write();
    ANNZResidual_vs_Conf_vs_Erecoil_Tgt2Iron->Write();
    ANNZResidual_vs_Conf_vs_Erecoil_Tgt2Lead->Write();
    ANNZResidual_vs_Conf_vs_Erecoil_Tgt3Iron->Write();
    ANNZResidual_vs_Conf_vs_Erecoil_Tgt3Lead->Write();
    ANNZResidual_vs_Conf_vs_Erecoil_Tgt3Carbon->Write();
    ANNZResidual_vs_Conf_vs_Erecoil_Tgt4Lead->Write();
    ANNZResidual_vs_Conf_vs_Erecoil_Tgt5Iron->Write();
    ANNZResidual_vs_Conf_vs_Erecoil_Tgt5Lead->Write();
    ANNZResidual_vs_Conf_vs_Erecoil_Water->Write();

    ANNWeightedZCalibratedResidual_vs_ConfCalibrated_vs_Erecoil->Write();
    ANNWeightedZResidual_vs_Conf_vs_Erecoil->Write();
    ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tracker->Write();
    ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt2Iron->Write();
    ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt2Lead->Write();
    ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt3Iron->Write();
    ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt3Lead->Write();
    ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt3Carbon->Write();
    ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt4Lead->Write();
    ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt5Iron->Write();
    ANNWeightedZResidual_vs_Conf_vs_Erecoil_Tgt5Lead->Write();
    ANNWeightedZResidual_vs_Conf_vs_Erecoil_Water->Write();

    ANNWeightedZ_minus_UnweightedZMC->Write();
    ANNWeightedZ_minus_UnweightedZTrackerMC->Write();
    ANNWeightedZ_minus_UnweightedZTgt2IronMC->Write();
    ANNWeightedZ_minus_UnweightedZTgt2LeadMC->Write();
    ANNWeightedZ_minus_UnweightedZTgt3IronMC->Write();
    ANNWeightedZ_minus_UnweightedZTgt3LeadMC->Write();
    ANNWeightedZ_minus_UnweightedZTgt3CarbonMC->Write();
    ANNWeightedZ_minus_UnweightedZTgt4LeadMC->Write();
    ANNWeightedZ_minus_UnweightedZTgt5IronMC->Write();
    ANNWeightedZ_minus_UnweightedZTgt5LeadMC->Write();
    ANNWeightedZ_minus_UnweightedZWaterMC->Write();

    ANNWeightedZ_minus_UnweightedZData->Write();
    ANNWeightedZ_minus_UnweightedZTrackerData->Write();
    ANNWeightedZ_minus_UnweightedZTgt2IronData->Write();
    ANNWeightedZ_minus_UnweightedZTgt2LeadData->Write();
    ANNWeightedZ_minus_UnweightedZTgt3IronData->Write();
    ANNWeightedZ_minus_UnweightedZTgt3LeadData->Write();
    ANNWeightedZ_minus_UnweightedZTgt3CarbonData->Write();
    ANNWeightedZ_minus_UnweightedZTgt4LeadData->Write();
    ANNWeightedZ_minus_UnweightedZTgt5IronData->Write();
    ANNWeightedZ_minus_UnweightedZTgt5LeadData->Write();
    ANNWeightedZ_minus_UnweightedZWaterData->Write();

    ANNConf_vs_Ehad_Tracker_Data->Write();
    ANNConf_vs_Ehad_Tgt2Iron_Data->Write();
    ANNConf_vs_Ehad_Tgt2Lead_Data->Write();
    ANNConf_vs_Ehad_Tgt3Iron_Data->Write();
    ANNConf_vs_Ehad_Tgt3Lead_Data->Write();
    ANNConf_vs_Ehad_Tgt3Carbon_Data->Write();
    ANNConf_vs_Ehad_Tgt4Lead_Data->Write();
    ANNConf_vs_Ehad_Tgt5Iron_Data->Write();
    ANNConf_vs_Ehad_Tgt5Lead_Data->Write();
    ANNConf_vs_Ehad_Water_Data->Write();

    ANNConf_vs_Ehad_Tracker_MC->Write();
    ANNConf_vs_Ehad_Tgt2Iron_MC->Write();
    ANNConf_vs_Ehad_Tgt2Lead_MC->Write();
    ANNConf_vs_Ehad_Tgt3Iron_MC->Write();
    ANNConf_vs_Ehad_Tgt3Lead_MC->Write();
    ANNConf_vs_Ehad_Tgt3Carbon_MC->Write();
    ANNConf_vs_Ehad_Tgt4Lead_MC->Write();
    ANNConf_vs_Ehad_Tgt5Iron_MC->Write();
    ANNConf_vs_Ehad_Tgt5Lead_MC->Write();
    ANNConf_vs_Ehad_Water_MC->Write();


    ANNConf_MisrecoInTgt2Iron->Write();
    ANNConf_MisrecoInTgt2Lead->Write();
    ANNConf_MisrecoInTgt3Iron->Write();
    ANNConf_MisrecoInTgt3Lead->Write();
    ANNConf_MisrecoInTgt3Carbon->Write();
    ANNConf_MisrecoInTgt4Lead->Write();
    ANNConf_MisrecoInTgt5Iron->Write();
    ANNConf_MisrecoInTgt5Lead->Write();
    ANNConf_MisrecoInWater->Write();

    ANNVerticesConfusion_ByZPos->Write();

  //Investigating distance to division cut performance. We're plotting the 2D postions of the events based on 
  //Investigating ANN reco performance by MC interaction channel
  //Confusion matrices
  ConfusionMatrixQE->Write();
  ConfusionMatrixRES->Write();
  ConfusionMatrixDIS->Write();
  ConfusionMatrixCOHPI->Write();
  ConfusionMatrixAMNUGAMMA->Write();
  ConfusionMatrixIMD->Write();
  ConfusionMatrixNUEEL->Write();
  ConfusionMatrixMEC->Write();
  ConfusionMatrixUnknownInt->Write();
  //Number of tracks
  ConfusionMatrixMultiplicity->Write();

  ConfusionMatrixErecoil->Write();
  ConfusionMatrixConf->Write();
  ConfusionMatrix->Write();

  //Confidence
  ANNConfidenceQE_MC->Write();
  ANNConfidenceRES_MC->Write();
  ANNConfidenceDIS_MC->Write();
  ANNConfidenceCOHPI_MC->Write();
  ANNConfidenceAMNUGAMMA_MC->Write();
  ANNConfidenceIMD_MC->Write();
  ANNConfidenceNUEEL_MC->Write();
  ANNConfidenceMEC_MC->Write();
  ANNConfidenceUnknownInt_MC->Write();

  ANNConfidence_Data->Write();
  ANNConfidence_MC->Write();
  ANNCalibratedConfidence_Data->Write();

  ANNCalibratedConfidence_MC->Write();
  
  ANNCalibratedConfidence_Iron_Selected->Write();
  ANNCalibratedConfidence1_Iron_Selected->Write();
  ANNCalibratedConfidence2_Iron_Selected->Write();
  ANNCalibratedConfidence3_Iron_Selected->Write();
  ANNCalibratedConfidence_Iron_SelectedSignal->Write();
  ANNCalibratedConfidence1_Iron_SelectedSignal->Write();
  ANNCalibratedConfidence2_Iron_SelectedSignal->Write();
  ANNCalibratedConfidence3_Iron_SelectedSignal->Write();
  ANNCalibratedConfidence_Iron_Signal->Write();
  ANNCalibratedConfidence1_Iron_Signal->Write();
  ANNCalibratedConfidence2_Iron_Signal->Write();
  ANNCalibratedConfidence3_Iron_Signal->Write();

  ANNCalibratedConfidence_Lead_Selected->Write();
  ANNCalibratedConfidence1_Lead_Selected->Write();
  ANNCalibratedConfidence2_Lead_Selected->Write();
  ANNCalibratedConfidence3_Lead_Selected->Write();
  ANNCalibratedConfidence_Lead_SelectedSignal->Write();
  ANNCalibratedConfidence1_Lead_SelectedSignal->Write();
  ANNCalibratedConfidence2_Lead_SelectedSignal->Write();
  ANNCalibratedConfidence3_Lead_SelectedSignal->Write();
  ANNCalibratedConfidence_Lead_Signal->Write();
  ANNCalibratedConfidence1_Lead_Signal->Write();
  ANNCalibratedConfidence2_Lead_Signal->Write();
  ANNCalibratedConfidence3_Lead_Signal->Write();

  ANNCalibratedConfidence_Carbon_Selected->Write();
  ANNCalibratedConfidence1_Carbon_Selected->Write();
  ANNCalibratedConfidence2_Carbon_Selected->Write();
  ANNCalibratedConfidence3_Carbon_Selected->Write();
  ANNCalibratedConfidence_Carbon_SelectedSignal->Write();
  ANNCalibratedConfidence1_Carbon_SelectedSignal->Write();
  ANNCalibratedConfidence2_Carbon_SelectedSignal->Write();
  ANNCalibratedConfidence3_Carbon_SelectedSignal->Write();
  ANNCalibratedConfidence_Carbon_Signal->Write();
  ANNCalibratedConfidence1_Carbon_Signal->Write();
  ANNCalibratedConfidence2_Carbon_Signal->Write();
  ANNCalibratedConfidence3_Carbon_Signal->Write();

  ANNCalibratedConfidence_Water_Selected->Write();
  ANNCalibratedConfidence1_Water_Selected->Write();
  ANNCalibratedConfidence2_Water_Selected->Write();
  ANNCalibratedConfidence3_Water_Selected->Write();
  ANNCalibratedConfidence_Water_SelectedSignal->Write();
  ANNCalibratedConfidence1_Water_SelectedSignal->Write();
  ANNCalibratedConfidence2_Water_SelectedSignal->Write();
  ANNCalibratedConfidence3_Water_SelectedSignal->Write();
  ANNCalibratedConfidence_Water_Signal->Write();
  ANNCalibratedConfidence1_Water_Signal->Write();
  ANNCalibratedConfidence2_Water_Signal->Write();
  ANNCalibratedConfidence3_Water_Signal->Write();




  //Iron

  ANNCalibratedConfidenceRatio2_Iron_SelectedSignal->Write();
  ANNCalibratedConfidenceRatio2_Iron_Signal->Write();
  ANNCalibratedConfidenceRatio2_Iron_Selected->Write();

  ANNCalibratedConfidenceRatio3_Iron_SelectedSignal->Write();
  ANNCalibratedConfidenceRatio3_Iron_Signal->Write();
  ANNCalibratedConfidenceRatio3_Iron_Selected->Write();
  //Lead

  ANNCalibratedConfidenceRatio2_Lead_SelectedSignal->Write();
  ANNCalibratedConfidenceRatio2_Lead_Signal->Write();
  ANNCalibratedConfidenceRatio2_Lead_Selected->Write();

  ANNCalibratedConfidenceRatio3_Lead_SelectedSignal->Write();
  ANNCalibratedConfidenceRatio3_Lead_Signal->Write();
  ANNCalibratedConfidenceRatio3_Lead_Selected->Write();
  //Carbon

  ANNCalibratedConfidenceRatio2_Carbon_SelectedSignal->Write();
  ANNCalibratedConfidenceRatio2_Carbon_Signal->Write();
  ANNCalibratedConfidenceRatio2_Carbon_Selected->Write();

  ANNCalibratedConfidenceRatio3_Carbon_SelectedSignal->Write();
  ANNCalibratedConfidenceRatio3_Carbon_Signal->Write();
  ANNCalibratedConfidenceRatio3_Carbon_Selected->Write();
  //Water

  ANNCalibratedConfidenceRatio2_Water_SelectedSignal->Write();
  ANNCalibratedConfidenceRatio2_Water_Signal->Write();
  ANNCalibratedConfidenceRatio2_Water_Selected->Write();

  ANNCalibratedConfidenceRatio3_Water_SelectedSignal->Write();
  ANNCalibratedConfidenceRatio3_Water_Signal->Write();
  ANNCalibratedConfidenceRatio3_Water_Selected->Write();


  // by multiplicity
  ANNConfidence_vs_NTracks_MC->Write();
  ANNConfidence_vs_NTracks_Data->Write();

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

    {
      std::cout<<"potData: " << potData << " potMC " << potMC <<std::endl;
      TCanvas *canvas = new TCanvas("canvas", "canvas", 2000, 1000);
      ANNWeightedZCalibratedMinusANNZ_Data->SetTitle("Calibrated Z Positions (weighted average) - Uncalibrated most probable position");
      ANNWeightedZCalibratedMinusANNZ_Data->SetYTitle("Percentage of events");
      ANNWeightedZCalibratedMinusANNZ_Data->SetXTitle("Distance (mm)");
      ANNWeightedZCalibratedMinusANNZ_Data->Scale(100/ANNWeightedZCalibratedMinusANNZ_Data->GetEntries());
      ANNWeightedZCalibratedMinusANNZ_Data->SetLineWidth(3);
      ANNWeightedZCalibratedMinusANNZ_Data->SetLineColor(kBlue);
      ANNWeightedZCalibratedMinusANNZ_Data->SetStats(0);
      ANNWeightedZCalibratedMinusANNZ_Data->Draw("HIST");
      ANNWeightedZCalibratedMinusANNZ_MC->SetLineWidth(3);
      ANNWeightedZCalibratedMinusANNZ_MC->Scale(100/ANNWeightedZCalibratedMinusANNZ_MC->GetEntries());
      ANNWeightedZCalibratedMinusANNZ_MC->SetLineColor(kRed);
      ANNWeightedZCalibratedMinusANNZ_MC->Draw("HIST SAME");
      TLegend *legend = new TLegend(0.8,0.7,0.9,0.8);
      legend->SetBorderSize(0);
      legend->AddEntry(ANNWeightedZCalibratedMinusANNZ_Data,"Data","l");
      legend->AddEntry(ANNWeightedZCalibratedMinusANNZ_MC,"MC","l");
      legend->Draw();
      canvas->SaveAs("WeightedZMinusANNZ.png");
    }
    {
      TCanvas *canvas = new TCanvas("canvas", "canvas", 2000, 1000);
      ANNZCalibratedMinusANNZ_Data->SetTitle("Most probable z position Calibrated - Uncalibrated");
      ANNZCalibratedMinusANNZ_Data->SetYTitle("Percentage of events");
      ANNZCalibratedMinusANNZ_Data->SetXTitle("Distance (mm)");
      ANNZCalibratedMinusANNZ_Data->Scale(100/ANNZCalibratedMinusANNZ_Data->GetEntries());
      ANNZCalibratedMinusANNZ_Data->SetLineWidth(3);
      ANNZCalibratedMinusANNZ_Data->SetLineColor(kBlue);
      ANNZCalibratedMinusANNZ_Data->SetStats(0);
      ANNZCalibratedMinusANNZ_Data->Draw("HIST");
      ANNZCalibratedMinusANNZ_MC->SetLineWidth(3);
      ANNZCalibratedMinusANNZ_MC->Scale(100/ANNZCalibratedMinusANNZ_MC->GetEntries());
      ANNZCalibratedMinusANNZ_MC->SetLineColor(kRed);
      ANNZCalibratedMinusANNZ_MC->Draw("HIST SAME");
      TLegend *legend = new TLegend(0.8,0.7,0.9,0.8);
      legend->SetBorderSize(0);
      legend->AddEntry(ANNZCalibratedMinusANNZ_Data,"Data","l");
      legend->AddEntry(ANNZCalibratedMinusANNZ_MC,"MC","l");
      legend->Draw();
      canvas->SaveAs("ZMinusANNZ.png");
    }

      //TCanvas *canvas = new TCanvas("canvas", "canvas", 2000, 1000);
      DrawDataMC(ANNZCalibrated_Data, ANNZCalibrated_MC, "ANNCalibratedZPositions", potMC, potData);
      DrawDataMC(ANNZCalibrated_Data, ANNZCalibrated_MC, "ANNCalibratedZPositionsZoomed", potMC, potData, true);
      DrawDataMC(ANNZCalibratedWeighted_Data, ANNZCalibratedWeighted_MC, "ANNCalibratedWeightedZPositions", potMC, potData);
      DrawDataMC(ANNZCalibratedWeighted_Data, ANNZCalibratedWeighted_MC, "ANNCalibratedWeightedZPositionsZoomed", potMC, potData, true);
      DrawDataMC(ANNZ_Data, ANNZ_MC, "ANNUncalibratedZPositions", potMC, potData);
      DrawDataMC(ANNZ_Data, ANNZ_MC, "ANNUncalibratedZPositionsZoomed", potMC, potData, true);

      DrawFullAndEmptyDataMC(ANNZCalibratedWeighted_Data, ANNZCalibratedWeighted_MC,ANNZ_Data, ANNZ_MC, "ANNCalibratedAndUncalibratedSplit", potMC, potData, potMC, potData);
      DrawFullAndEmptyDataMC(ANNZCalibratedWeighted_Data, ANNZCalibratedWeighted_MC,ANNZ_Data, ANNZ_MC, "ANNCalibratedAndUncalibratedSplitZoomed", potMC, potData, potMC, potData, true);

      TH1D* selected = (TH1D*)ANNCalibratedConfidence_Iron_Selected->GetCumulative(false);
      TH1D* selected1 = (TH1D*)ANNCalibratedConfidence1_Iron_Selected->GetCumulative(false);
      TH1D* selected2 = (TH1D*)ANNCalibratedConfidence2_Iron_Selected->GetCumulative(false);
      TH1D* selected3 = (TH1D*)ANNCalibratedConfidence3_Iron_Selected->GetCumulative(false);

      TH1D* selectedRatio2 = (TH1D*)ANNCalibratedConfidenceRatio2_Iron_Selected->GetCumulative(false);
      TH1D* selectedRatio3 = (TH1D*)ANNCalibratedConfidenceRatio3_Iron_Selected->GetCumulative(false);

      TH1D* selectedSignal = (TH1D*)ANNCalibratedConfidence_Iron_SelectedSignal->GetCumulative(false);
      TH1D* selectedSignal1 = (TH1D*)ANNCalibratedConfidence1_Iron_SelectedSignal->GetCumulative(false);
      TH1D* selectedSignal2 = (TH1D*)ANNCalibratedConfidence2_Iron_SelectedSignal->GetCumulative(false);
      TH1D* selectedSignal3 = (TH1D*)ANNCalibratedConfidence3_Iron_SelectedSignal->GetCumulative(false);

      TH1D* selectedSignalRatio2 = (TH1D*)ANNCalibratedConfidenceRatio2_Iron_SelectedSignal->GetCumulative(false);
      TH1D* selectedSignalRatio3 = (TH1D*)ANNCalibratedConfidenceRatio3_Iron_SelectedSignal->GetCumulative(false);

      double signal = ANNCalibratedConfidence_Iron_Signal->GetSumOfWeights();
      //TH1D* signal2 = (TH1D*)ANNCalibratedConfidence2_Iron_Signal->GetCumulative(false);
      //TH1D* signal3 = (TH1D*)ANNCalibratedConfidence3_Iron_Signal->GetCumulative(false);

      TH1D* purity = (TH1D*)selectedSignal->Clone("Purity");
      TH1D* purity1 = (TH1D*)selectedSignal1->Clone("Purity");
      TH1D* purity2 = (TH1D*)selectedSignal2->Clone("Purity");
      TH1D* purity3 = (TH1D*)selectedSignal3->Clone("Purity");

      TH1D* purityRatio2 = (TH1D*)selectedSignalRatio2->Clone("Purity");
      TH1D* purityRatio3 = (TH1D*)selectedSignalRatio3->Clone("Purity");

      TH1D* recall = (TH1D*)selectedSignal->Clone("Recall");
      TH1D* recall1 = (TH1D*)selectedSignal1->Clone("Recall");
      TH1D* recall2 = (TH1D*)selectedSignal2->Clone("Recall");
      TH1D* recall3 = (TH1D*)selectedSignal3->Clone("Recall");

      TH1D* recallRatio2 = (TH1D*)selectedSignalRatio2->Clone("Recall");
      TH1D* recallRatio3 = (TH1D*)selectedSignalRatio3->Clone("Recall");

      TH1D* F1Score = (TH1D*)selectedSignal->Clone("F1 Score");
      TH1D* F1Score1 = (TH1D*)selectedSignal1->Clone("F1 Score");
      TH1D* F1Score2 = (TH1D*)selectedSignal2->Clone("F1 Score");
      TH1D* F1Score3 = (TH1D*)selectedSignal3->Clone("F1 Score");

      TH1D* F1ScoreRatio2 = (TH1D*)selectedSignalRatio2->Clone("F1 Score");
      TH1D* F1ScoreRatio3 = (TH1D*)selectedSignalRatio3->Clone("F1 Score");

      purity->Divide(selected);
      purity1->Divide(selected1);
      purity2->Divide(selected2);
      purity3->Divide(selected3);

      purityRatio2->Divide(selectedRatio2);
      purityRatio3->Divide(selectedRatio3);

      recall->Scale(1/signal);
      recall1->Scale(1/signal);
      recall2->Scale(1/signal);
      recall3->Scale(1/signal);

      recallRatio2->Scale(1/signal);
      recallRatio3->Scale(1/signal);

      TH1D* PrecisionRecallSum = (TH1D*)recall->Clone("Precision + Recall");
      TH1D* PrecisionRecallSum1 = (TH1D*)recall1->Clone("Precision + Recall");
      TH1D* PrecisionRecallSum2 = (TH1D*)recall2->Clone("Precision + Recall");
      TH1D* PrecisionRecallSum3 = (TH1D*)recall3->Clone("Precision + Recall");

      TH1D* PrecisionRecallSumRatio2 = (TH1D*)recall2->Clone("Precision + Recall");
      TH1D* PrecisionRecallSumRatio3 = (TH1D*)recall3->Clone("Precision + Recall");

      PrecisionRecallSum->Add(purity);
      PrecisionRecallSum1->Add(purity1);
      PrecisionRecallSum2->Add(purity2);
      PrecisionRecallSum3->Add(purity3);

      PrecisionRecallSumRatio2->Add(purityRatio2);
      PrecisionRecallSumRatio3->Add(purityRatio3);

      F1Score->Multiply(purity, recall, 2);
      F1Score1->Multiply(purity1, recall1, 2);
      F1Score2->Multiply(purity2, recall2, 2);
      F1Score3->Multiply(purity3, recall3, 2);

      F1ScoreRatio2->Multiply(purityRatio2, recallRatio2, 2);
      F1ScoreRatio3->Multiply(purityRatio3, recallRatio3, 2);

      F1Score->Divide(PrecisionRecallSum);
      F1Score1->Divide(PrecisionRecallSum1);
      F1Score2->Divide(PrecisionRecallSum2);
      F1Score3->Divide(PrecisionRecallSum3);

      F1ScoreRatio2->Divide(PrecisionRecallSumRatio2);
      F1ScoreRatio3->Divide(PrecisionRecallSumRatio3);

    {
      std::cout<<"potData: " << potData << " potMC " << potMC <<std::endl;
      TCanvas *canvas = new TCanvas("canvas", "canvas", 2000, 1000);
      selectedSignal1->SetTitle("Selected Signal");
      selectedSignal1->SetYTitle("Number of events");
      selectedSignal1->SetXTitle("Calibrated ANN Confidence");
      selectedSignal1->SetLineWidth(3);
      selectedSignal1->SetLineColor(kBlue);
      selectedSignal1->SetStats(0);
      selectedSignal1->Draw("HIST");

      selectedSignal2->SetLineWidth(3);
      selectedSignal2->SetLineColor(kRed);
      selectedSignal2->Draw("HIST SAME");

      selectedSignal3->SetLineWidth(3);
      selectedSignal3->SetLineColor(kGreen);
      selectedSignal3->Draw("HIST SAME");

      /* selectedSignalRatio2->SetLineWidth(3);
      selectedSignalRatio2->SetLineColor(kOrange);
      selectedSignalRatio2->Draw("HIST SAME");

      selectedSignalRatio3->SetLineWidth(3);
      selectedSignalRatio3->SetLineColor(kPink);
      selectedSignalRatio3->Draw("HIST SAME"); */

      TLegend *legend = new TLegend(0.8,0.7,0.9,0.8);
      legend->SetBorderSize(0);
      legend->AddEntry(selectedSignal1,"Top Plane","l");
      legend->AddEntry(selectedSignal2,"2nd Plane","l");
      legend->AddEntry(selectedSignal3,"3rd Plane","l");
      //legend->AddEntry(selectedSignalRatio2,"2nd Plane/1st Plane Ratio","l");
      //legend->AddEntry(selectedSignalRatio3,"3rd Plane/2nd Plane Ratio","l");
      legend->Draw();
      canvas->SaveAs("SelectedSignal.png");
    }

    {
      std::cout<<"potData: " << potData << " potMC " << potMC <<std::endl;
      TCanvas *canvas = new TCanvas("canvas", "canvas", 2000, 1000);
      selected1->SetTitle("Selected");
      selected1->SetYTitle("Number of events");
      selected1->SetXTitle("Calibrated ANN Confidence");
      selected1->SetLineWidth(3);
      selected1->SetLineColor(kBlue);
      selected1->SetStats(0);
      selected1->Draw("HIST");

      selected2->SetLineWidth(3);
      selected2->SetLineColor(kRed);
      selected2->Draw("HIST SAME");

      selected3->SetLineWidth(3);
      selected3->SetLineColor(kGreen);
      selected3->Draw("HIST SAME");

      /* selectedRatio2->SetLineWidth(3);
      selectedRatio2->SetLineColor(kOrange);
      selectedRatio2->Draw("HIST SAME");

      selectedRatio3->SetLineWidth(3);
      selectedRatio3->SetLineColor(kPink);
      selectedRatio3->Draw("HIST SAME"); */

      TLegend *legend = new TLegend(0.8,0.7,0.9,0.8);
      legend->SetBorderSize(0);
      legend->AddEntry(selected1,"Top Plane","l");
      legend->AddEntry(selected2,"2nd Plane","l");
      legend->AddEntry(selected3,"3rd Plane","l");
      //legend->AddEntry(selectedRatio2,"2nd Plane/1st Plane Ratio","l");
      //legend->AddEntry(selectedRatio3,"3rd Plane/2nd Plane Ratio","l");
      legend->Draw();
      canvas->SaveAs("Selected.png");
    }
    /* std::cout<<"Total number of signal events: " << signal << ":" << signal2 << ":" << signal3 << std::endl;
    {
      std::cout<<"potData: " << potData << " potMC " << potMC <<std::endl;
      TCanvas *canvas = new TCanvas("canvas", "canvas", 2000, 1000);
      signal->SetTitle("Signal");
      signal->SetYTitle("Number of events");
      signal->SetXTitle("Calibrated ANN Confidence");
      signal->SetLineWidth(3);
      signal->SetLineColor(kBlue);
      signal->SetStats(0);
      signal->Draw("HIST");

      signal2->SetLineWidth(3);
      signal2->SetLineColor(kRed);
      signal2->Draw("HIST SAME");

      signal3->SetLineWidth(3);
      signal3->SetLineColor(kGreen);
      signal3->Draw("HIST SAME");

      TLegend *legend = new TLegend(0.8,0.7,0.9,0.8);
      legend->SetBorderSize(0);
      legend->AddEntry(signal,"Top Plane","l");
      legend->AddEntry(signal2,"2nd Plane","l");
      legend->AddEntry(signal3,"3rd Plane","l");
      legend->Draw();
      canvas->SaveAs("Signal.png");
    } */

    {
      std::cout<<"potData: " << potData << " potMC " << potMC <<std::endl;
      TCanvas *canvas = new TCanvas("canvas", "canvas", 2000, 1000);
      purity1->SetTitle("Purity");
      purity1->SetYTitle("Number of events");
      purity1->SetXTitle("Calibrated ANN Confidence");
      purity1->SetLineWidth(3);
      purity1->SetLineColor(kBlue);
      purity1->SetStats(0);
      purity1->Draw("HIST");

      purity2->SetLineWidth(3);
      purity2->SetLineColor(kRed);
      purity2->Draw("HIST SAME");

      purity3->SetLineWidth(3);
      purity3->SetLineColor(kGreen);
      purity3->Draw("HIST SAME");

      /* purityRatio2->SetLineWidth(3);
      purityRatio2->SetLineColor(kOrange);
      purityRatio2->Draw("HIST SAME");

      purityRatio3->SetLineWidth(3);
      purityRatio3->SetLineColor(kPink);
      purityRatio3->Draw("HIST SAME"); */

      TLegend *legend = new TLegend(0.8,0.7,0.9,0.8);
      legend->SetBorderSize(0);
      legend->AddEntry(purity1,"Top Plane","l");
      legend->AddEntry(purity2,"2nd Plane","l");
      legend->AddEntry(purity3,"3rd Plane","l");
      //legend->AddEntry(purityRatio2,"2nd Plane/1st Plane Ratio","l");
      //legend->AddEntry(purityRatio3,"3rd Plane/2nd Plane Ratio","l");
      legend->Draw();
      canvas->SaveAs("Purity.png");
    }

    {
      std::cout<<"potData: " << potData << " potMC " << potMC <<std::endl;
      TCanvas *canvas = new TCanvas("canvas", "canvas", 2000, 1000);
      recall1->SetTitle("Recall");
      recall1->SetYTitle("Number of events");
      recall1->SetXTitle("Calibrated ANN Confidence");
      recall1->SetLineWidth(3);
      recall1->SetLineColor(kBlue);
      recall1->SetStats(0);
      recall1->Draw("HIST");

      recall2->SetLineWidth(3);
      recall2->SetLineColor(kRed);
      recall2->Draw("HIST SAME");

      recall3->SetLineWidth(3);
      recall3->SetLineColor(kGreen);
      recall3->Draw("HIST SAME");

      /* recallRatio2->SetLineWidth(3);
      recallRatio2->SetLineColor(kOrange);
      recallRatio2->Draw("HIST SAME");

      recallRatio3->SetLineWidth(3);
      recallRatio3->SetLineColor(kPink);
      recallRatio3->Draw("HIST SAME"); */

      TLegend *legend = new TLegend(0.8,0.7,0.9,0.8);
      legend->SetBorderSize(0);
      legend->AddEntry(recall1,"Top Plane","l");
      legend->AddEntry(recall2,"2nd Plane","l");
      legend->AddEntry(recall3,"3rd Plane","l");
      //legend->AddEntry(recallRatio2,"2nd Plane/1st Plane Ratio","l");
      //legend->AddEntry(recallRatio3,"3rd Plane/2nd Plane Ratio","l");
      legend->Draw();
      canvas->SaveAs("Recall.png");
    }

    {
      std::cout<<"potData: " << potData << " potMC " << potMC <<std::endl;
      TCanvas *canvas = new TCanvas("canvas", "canvas", 2000, 1000);
      F1Score1->SetTitle("F1Score");
      F1Score1->SetYTitle("Number of events");
      F1Score1->SetXTitle("Calibrated ANN Confidence");
      F1Score1->SetLineWidth(3);
      F1Score1->SetLineColor(kBlue);
      F1Score1->SetStats(0);
      F1Score1->Draw("HIST");

      F1Score2->SetLineWidth(3);
      F1Score2->SetLineColor(kRed);
      F1Score2->Draw("HIST SAME");

      F1Score3->SetLineWidth(3);
      F1Score3->SetLineColor(kGreen);
      F1Score3->Draw("HIST SAME");

      /* F1ScoreRatio2->SetLineWidth(3);
      F1ScoreRatio2->SetLineColor(kOrange);
      F1ScoreRatio2->Draw("HIST SAME");

      F1ScoreRatio3->SetLineWidth(3);
      F1ScoreRatio3->SetLineColor(kPink);
      F1ScoreRatio3->Draw("HIST SAME"); */

      TLegend *legend = new TLegend(0.8,0.7,0.9,0.8);
      legend->SetBorderSize(0);
      legend->AddEntry(F1Score1,"Top Plane","l");
      legend->AddEntry(F1Score2,"2nd Plane","l");
      legend->AddEntry(F1Score3,"3rd Plane","l");
      //legend->AddEntry(F1ScoreRatio2,"2nd Plane/1st Plane Ratio","l");
      //legend->AddEntry(F1ScoreRatio3,"3rd Plane/2nd Plane Ratio","l");
      legend->Draw();
      canvas->SaveAs("F1Score.png");
    }

    /* ANNZCalibrated_Data
    ANNZCalibrated_MC
    ANNZCalibratedWeighted_Data
    ANNZCalibratedWeighted_MC
    ANNZ_Data
    ANNZ_MC
    */
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