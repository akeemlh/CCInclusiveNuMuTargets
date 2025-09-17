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
#include "PlotUtils/CrashOnROOTMessage.h" //Sets up ROOT's debug callbacks by itself
#include "PlotUtils/Cutter.h"
#include "PlotUtils/Model.h"
#include "PlotUtils/FluxAndCVReweighter.h"
#include "PlotUtils/GENIEReweighter.h"
#include "PlotUtils/LowRecoil2p2hReweighter.h"
#include "PlotUtils/RPAReweighter.h"
#include "PlotUtils/MINOSEfficiencyReweighter.h"
#include "PlotUtils/TargetUtils.h"
#include "util/NukeUtils.h"
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


std::vector<double> vertexBins = {4293.04, 4337.25, 4381.47, 4425.68, 4514.11, 4558.33, 4602.54, 4646.76, 4735.19, 4779.4, 4823.62, 4867.83, 5000.48, 5044.69, 5088.91, 5133.12, 5456.74, 5500.95, 5545.17, 5589.38, 5677.81, 5722.03, 5810.45, 5855.68, 5900.91, 5946.14, 5991.37, 6036.6, 6081.83, 6127.06, 6172.29, 6217.52, 6262.74, 6307.97, 6353.2, 6398.43, 6443.66, 6488.89, 6534.12, 6579.35, 6624.58, 6669.81, 6715.03, 6760.26, 6805.49, 6850.72, 6895.95, 6941.18, 6986.41, 7031.64, 7076.87, 7122.1, 7167.32, 7212.55, 7257.78, 7303.01, 7348.24, 7393.47, 7438.7, 7483.93, 7529.16, 7574.39, 7619.61, 7664.84, 7710.07, 7755.3, 7800.53, 7845.76, 7890.99, 7936.22, 7981.45, 8026.68, 8071.9, 8117.13, 8162.36, 8207.59, 8252.82, 8298.05, 8343.28, 8388.51, 8433.74, 8478.97, 8524.19, 8569.42, 8614.65};

double rebinNum = 10;
double nbins = 3400;
double xlow = 4200;
double xhigh = 5900;


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
    PlotUtils::Model<CVUniverse, MichelEvent>& model,
    TH1D& ANNVerticesFilledMC_ByModule,
    TH1D& ANNVerticesFilledMC_ByZPos,
    TH1D& TBVerticesFilledMC_ByModule,
    TH1D& TBVerticesFilledMC_ByZPos,
    TH1D& TruthVerticesFilledMC_ByModule,
    TH1D& TruthVerticesFilledMC_ByZPos,
    TH2D& ANNVerticesFilledMCERecoil_ByModule,
    TH2D& ANNVerticesFilledMCMultiplicity_ByModule,
    TH2D& TBVerticesFilledMCERecoil_ByModule,
    TH2D& TBVerticesFilledMCMultiplicity_ByModule,
    TH2D& TruthVerticesFilledMCERecoil_ByModule,
    TH2D& TruthVerticesFilledMCMultiplicity_ByModule
    )
{
  std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
  assert(!error_bands["cv"].empty() && "\"cv\" error band is empty!  Can't set Model weight.");
  auto& cvUniv = error_bands["cv"].front();

  std::cout << "Starting MC reco loop...\n";

  const int nEntries = chain->GetEntries();
  //const int nEntries = 10000;
  for (int i=0; i<nEntries; ++i)
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
      //std::cout<<"Here4\n";
      std::vector<CVUniverse*> error_band_universes = band.second;
      int univCount = 0;
      for (auto universe : error_band_universes)
      {
        univCount++; // Put the iterator right at the start so it's executed even in paths that lead to a continue, don't forget to subtract by 1 when we use it
        MichelEvent myevent; // make sure your event is inside the error band loop. 
        // Tell the Event which entry in the TChain it's looking at
        universe->SetEntry(i);
        std::vector<double> ANNVtx = cvUniv->GetANNVertexVector();
        ROOT::Math::XYZTVector TrackBasedVtx = cvUniv->GetVertex();
        ROOT::Math::XYZTVector TrueVtx = cvUniv->GetTrueVertex();
        // This is where you would Access/create a Michel

        //weight is ignored in isMCSelected() for all but the CV Universe.
        if (!michelcuts.isMCSelected(*universe, myevent, cvWeight).all()) continue; //all is another function that will later help me with sidebands
        const double weight = model.GetWeight(*universe, myevent); //Only calculate the per-universe weight for events that will actually use it.
        PlotUtils::TargetUtils* m_TargetUtils=new PlotUtils::TargetUtils();

        //Performing vtx validation check Deborah suggested
        double batchPOT = cvUniv->GetBatchPOT();
        double efficiency = 0.5563 - (0.01353*batchPOT); //Based on MINERvA-doc-21436
        //cvWeight/=efficiency; //Not doing efficiency correction right now

        //Hadron Energy Spectrum plot
        double erecoil = cvUniv->GetRecoilE()/pow(10,3);


        if (m_TargetUtils->InWaterTargetVolMC(TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z()) && (m_TargetUtils->InWaterTargetVolMC(TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z())))
        {
          //These events may be visually interesting - check out in event display
          //std::cout<<"Check in arachne: \n" << " ev_run: " << cvUniv->GetInt("ev_run") << " ev_subrun: " << cvUniv->GetInt("ev_subrun") << " ev_gate: " << cvUniv->GetInt("ev_gate") << std::endl;
        }
        double ANNProb = cvUniv->GetANNProb();
        double ANNX, ANNY, ANNZ;
        if(ANNVtx.size()==3)
        {
          ANNX = ANNVtx[0];
          ANNY = ANNVtx[1];
          ANNZ = ANNVtx[2];
        }

        if(ANNVtx.size()==3)
        {
          ANNVerticesFilledMC_ByModule.Fill( ANNVtx[2], cvWeight);
          ANNVerticesFilledMC_ByZPos.Fill( ANNVtx[2], cvWeight);

          ANNVerticesFilledMCERecoil_ByModule.Fill( ANNVtx[2], erecoil, cvWeight);
          ANNVerticesFilledMCMultiplicity_ByModule.Fill( ANNVtx[2], cvUniv->GetMultiplicity(), cvWeight);
        }
        TBVerticesFilledMC_ByModule.Fill(TrackBasedVtx.Z(), cvWeight);
        TBVerticesFilledMC_ByZPos.Fill( TrackBasedVtx.Z(), cvWeight);
        TruthVerticesFilledMC_ByModule.Fill( TrueVtx.Z(), cvWeight);
        TruthVerticesFilledMC_ByZPos.Fill( TrueVtx.Z(), cvWeight);

        TBVerticesFilledMCERecoil_ByModule.Fill( TrackBasedVtx.Z(), erecoil, cvWeight);
        TBVerticesFilledMCMultiplicity_ByModule.Fill( TrackBasedVtx.Z(), cvUniv->GetMultiplicity(), cvWeight);
        TruthVerticesFilledMCERecoil_ByModule.Fill( TrueVtx.Z(), erecoil, cvWeight);
        TruthVerticesFilledMCMultiplicity_ByModule.Fill( TrueVtx.Z(), cvUniv->GetMultiplicity(), cvWeight);

      }
    }
  } //End entries loop
  std::cout << "Finished MC reco loop.\n";
}

void LoopAndFillData( PlotUtils::ChainWrapper* data,
			        std::vector<CVUniverse*> data_band,
				PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts,
        TH1D& ANNVerticesFilledData_ByModule,
        TH1D& ANNVerticesFilledData_ByZPos,
        TH1D& TBVerticesFilledData_ByZPos,
        TH1D& TBVerticesFilledData_ByModule,
        TH2D& ANNVerticesFilledDataERecoil_ByModule,
        TH2D& TBVerticesFilledDataERecoil_ByModule,
        TH2D& ANNVerticesFilledDataMultiplicity_ByModule,
        TH2D& TBVerticesFilledDataMultiplicity_ByModule
        )

{
  std::cout << "Starting data loop...\n";
  std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
  const int nEntries = data->GetEntries();
  //const int nEntries = 10000;
  for (int i=0; i<nEntries; ++i) {
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
      std::vector<double> ANNVtx = universe->GetANNVertexVector();
      ROOT::Math::XYZTVector TrackBasedVtx = universe->GetVertex();
      double batchPOT = universe->GetBatchPOT();
      //Incorporate batch POT efficiency scaling as applied above
      double erecoil = universe->GetRecoilE()/pow(10,3);

      if (!michelcuts.isDataSelected(*universe, myevent).all()) continue;

      double ANNX, ANNY, ANNZ;
      if(ANNVtx.size()==3)
      {
        ANNX = ANNVtx[0];
        ANNY = ANNVtx[1];
        ANNZ = ANNVtx[2];
      }

       if(ANNVtx.size()==3)
      {
        ANNVerticesFilledData_ByModule.Fill(ANNVtx[2]);
        ANNVerticesFilledData_ByZPos.Fill(ANNVtx[2]);

        ANNVerticesFilledDataERecoil_ByModule.Fill( ANNVtx[2], erecoil);
        ANNVerticesFilledDataMultiplicity_ByModule.Fill( ANNVtx[2], universe->GetMultiplicity());
      }
      TBVerticesFilledData_ByZPos.Fill(TrackBasedVtx.Z());
      TBVerticesFilledData_ByModule.Fill(TrackBasedVtx.Z());

      TBVerticesFilledDataERecoil_ByModule.Fill( TrackBasedVtx.Z(), erecoil);
      TBVerticesFilledDataMultiplicity_ByModule.Fill( TrackBasedVtx.Z(), universe->GetMultiplicity());
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
int main(const int argc, char** argv)
{
  TH1::AddDirectory(false);

  //Validate input.
  //I expect a data playlist file name and an MC playlist file name which is exactly 2 arguments.
  const int nArgsExpected = 5;
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
  const std::string mc_file_list_empty = argv[4],
                    data_file_list_empty = argv[3],
                    mc_file_list_filled = argv[2],
                    data_file_list_filled = argv[1];

  //Check that necessary TTrees exist in the first file of mc_file_list_filled and data_file_list_filled
  std::string reco_tree_name;
  if(!inferRecoTreeNameAndCheckTreeNames(mc_file_list_filled, data_file_list_filled, reco_tree_name))
  {
    std::cerr << "Failed to find required trees in MC playlist " << mc_file_list_filled << " and/or data playlist " << data_file_list_filled << ".\n" << USAGE << "\n";
    return badInputFile;
  }

  const bool doCCQENuValidation = (reco_tree_name == "CCQENu"); //Enables extra histograms and might influence which systematics I use.

  //const bool is_grid = false; //TODO: Are we going to put this back?  Gonzalo needs it iirc.
  PlotUtils::MacroUtil optionsFilled(reco_tree_name, mc_file_list_filled, data_file_list_filled, "minervame1A", true); //minervame1A is just a placeholder, it gets overwritted immediately below
  PlotUtils::MacroUtil optionsEmpty(reco_tree_name, mc_file_list_empty, data_file_list_empty, "minervame1A", true); //minervame1A is just a placeholder, it gets overwritted immediately below
  optionsFilled.m_plist_string = util::GetPlaylist(*optionsFilled.m_mc, true); //TODO: Put GetPlaylist into PlotUtils::MacroUtil

  // You're required to make some decisions
  PlotUtils::MinervaUniverse::SetNuEConstraint(true);
  PlotUtils::MinervaUniverse::SetPlaylist(optionsFilled.m_plist_string); //TODO: Infer this from the files somehow?
  int nuoranu = nuOrAntiNuMode(optionsFilled.m_plist_string);
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

  //const double minZ = 5980, maxZ = 8422, apothem = 850; //All in mm
  const double apothem = 850; //All in mm
  preCuts.emplace_back(new reco::ZRangeANN<CVUniverse, MichelEvent>("Active Tracker Z pos", PlotUtils::TargetProp::NukeRegion::Face, PlotUtils::TargetProp::Tracker::Back));
  preCuts.emplace_back(new reco::Apothem<CVUniverse, MichelEvent>(apothem));
  preCuts.emplace_back(new reco::MaxMuonAngle<CVUniverse, MichelEvent>(17.));
  preCuts.emplace_back(new reco::HasMINOSMatch<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::NoDeadtime<CVUniverse, MichelEvent>(1, "Deadtime"));
  preCuts.emplace_back(new reco::IsNeutrino<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::MuonEnergyMin<CVUniverse, MichelEvent>(2000.0, "EMu Min"));
  preCuts.emplace_back(new reco::MuonEnergyMax<CVUniverse, MichelEvent>(50000.0, "EMu Max"));
  preCuts.emplace_back(new reco::ANNConfidenceCut<CVUniverse, MichelEvent>(0.40));


                                                                                                                                   
  signalDefinition.emplace_back(new truth::IsNeutrino<CVUniverse>());
  signalDefinition.emplace_back(new truth::IsCC<CVUniverse>());
                                                                                                                                                   
  phaseSpace.emplace_back(new truth::ZRange<CVUniverse>("Active Tracker Z pos", PlotUtils::TargetProp::NukeRegion::Face, PlotUtils::TargetProp::Tracker::Back));
  phaseSpace.emplace_back(new truth::Apothem<CVUniverse>(apothem));
  phaseSpace.emplace_back(new truth::MuonAngle<CVUniverse>(17.));
  phaseSpace.emplace_back(new truth::MuonEnergyMin<CVUniverse>(2000.0, "EMu Min"));
  phaseSpace.emplace_back(new truth::MuonEnergyMax<CVUniverse>(50000.0, "EMu Max"));
                                                                                                                                                   
  PlotUtils::Cutter<CVUniverse, MichelEvent> mycuts(std::move(preCuts), std::move(sidebands) , std::move(signalDefinition),std::move(phaseSpace));

  std::vector<std::unique_ptr<PlotUtils::Reweighter<CVUniverse, MichelEvent>>> MnvTunev1;
  MnvTunev1.emplace_back(new PlotUtils::FluxAndCVReweighter<CVUniverse, MichelEvent>());
  MnvTunev1.emplace_back(new PlotUtils::GENIEReweighter<CVUniverse, MichelEvent>(true, false));
  MnvTunev1.emplace_back(new PlotUtils::LowRecoil2p2hReweighter<CVUniverse, MichelEvent>());
  MnvTunev1.emplace_back(new PlotUtils::MINOSEfficiencyReweighter<CVUniverse, MichelEvent>());
  MnvTunev1.emplace_back(new PlotUtils::RPAReweighter<CVUniverse, MichelEvent>());
 
  PlotUtils::Model<CVUniverse, MichelEvent> model(std::move(MnvTunev1));

  PlotUtils::MinervaUniverse::SetNFluxUniverses(2); //Necessary to get Flux integral later...  Doesn't work with just 1 flux universe though because _that_ triggers "spread errors".


  //==============================================================================
  // Global - Declaring histograms
  //==============================================================================
  TH1D *ANNVerticesFilledMC_ByModule = new TH1D ("ANNVerticesFilledMC_ByModule", "ANNVerticesFilledMC_ByModule", vertexBins.size()-1, &vertexBins[0]);
  TH1D *TBVerticesFilledMC_ByModule = new TH1D ("TBVerticesFilledMC_ByModule", "TBVerticesFilledMC_ByModule", vertexBins.size()-1, &vertexBins[0]);
  TH1D *ANNVerticesFilledData_ByModule = new TH1D ("ANNVerticesFilledData_ByModule", "ANNVerticesFilledData_ByModule", vertexBins.size()-1, &vertexBins[0]);
  TH1D *TBVerticesFilledData_ByModule = new TH1D ("TBVerticesFilledData_ByModule", "TBVerticesFilledData_ByModule", vertexBins.size()-1, &vertexBins[0]);
  TH1D *TruthVerticesFilledMC_ByModule = new TH1D ("TruthVerticesFilledMC_ByModule", "TruthVerticesFilledMC_ByModule", vertexBins.size()-1, &vertexBins[0]);
  TH1D *ANNVerticesFilledMC_ByZPos = new TH1D ("ANNVerticesFilledMC_ByZPos", "ANNVerticesFilledMC_ByZPos", 3400, 4200, 5900);
  TH1D *TBVerticesFilledMC_ByZPos = new TH1D ("TBVerticesFilledMC_ByZPos", "TBVerticesFilledMC_ByModuleTruthVerticesFilledMC_ByZPosHighRes", 3400, 4200, 5900);
  TH1D *ANNVerticesFilledData_ByZPos = new TH1D ("ANNVerticesFilledData_ByZPos", "ANNVerticesFilledData_ByZPos", 3400, 4200, 5900);
  TH1D *TBVerticesFilledData_ByZPos = new TH1D ("TBVerticesFilledData_ByZPos", "TBVerticesFilledData_ByZPos", 3400, 4200, 5900);
  TH1D *TruthVerticesFilledMC_ByZPos = new TH1D ("TruthVerticesFilledMC_ByZPos", "TruthVerticesFilledMC_ByZPos", 3400, 4200, 5900);
  //ERecoil
  TH2D *TruthVerticesFilledMCERecoil_ByModule = new TH2D ("TruthVerticesFilledMCERecoil_ByModule", "TruthVerticesFilledMCERecoil_ByModule", vertexBins.size()-1, &vertexBins[0], 100, 0, 50);
  TH2D *ANNVerticesFilledMCERecoil_ByModule = new TH2D ("ANNVerticesFilledMCERecoil_ByModule", "ANNVerticesFilledMCERecoil_ByModule", vertexBins.size()-1, &vertexBins[0], 100, 0, 50);
  TH2D *TBVerticesFilledMCERecoil_ByModule = new TH2D ("TBVerticesFilledMCERecoil_ByModule", "TBVerticesFilledMCERecoil_ByModule", vertexBins.size()-1, &vertexBins[0], 100, 0, 50);
  TH2D *ANNVerticesFilledDataERecoil_ByModule = new TH2D ("ANNVerticesFilledDataERecoil_ByModule", "ANNVerticesFilledDataERecoil_ByModule", vertexBins.size()-1, &vertexBins[0], 100, 0, 50);
  TH2D *TBVerticesFilledDataERecoil_ByModule = new TH2D ("TBVerticesFilledDataERecoil_ByModule", "TBVerticesFilledDataERecoil_ByModule", vertexBins.size()-1, &vertexBins[0], 100, 0, 50);
  //Multiplicity
  TH2D *TruthVerticesFilledMCMultiplicity_ByModule = new TH2D ("TruthVerticesFilledMCMultiplicity_ByModule", "TruthVerticesFilledMCMultiplicity_ByModule", vertexBins.size()-1, &vertexBins[0], 10, 0, 10);
  TH2D *ANNVerticesFilledMCMultiplicity_ByModule = new TH2D ("ANNVerticesFilledMCMultiplicity_ByModule", "ANNVerticesFilledMCMultiplicity_ByModule", vertexBins.size()-1, &vertexBins[0], 10, 0, 10);
  TH2D *TBVerticesFilledMCMultiplicity_ByModule = new TH2D ("TBVerticesFilledMCMultiplicity_ByModule", "TBVerticesFilledMCMultiplicity_ByModule", vertexBins.size()-1, &vertexBins[0], 10, 0, 10);
  TH2D *ANNVerticesFilledDataMultiplicity_ByModule = new TH2D ("ANNVerticesFilledDataMultiplicity_ByModule", "ANNVerticesFilledDataMultiplicity_ByModule", vertexBins.size()-1, &vertexBins[0], 10, 0, 10);
  TH2D *TBVerticesFilledDataMultiplicity_ByModule = new TH2D ("TBVerticesFilledDataMultiplicity_ByModule", "TBVerticesFilledDataMultiplicity_ByModule", vertexBins.size()-1, &vertexBins[0], 10, 0, 10);


  TH1D *ANNVerticesEmptyMC_ByModule = new TH1D ("ANNVerticesEmptyMC_ByModule", "ANNVerticesEmptyMC_ByModule", vertexBins.size()-1, &vertexBins[0]);
  TH1D *TBVerticesEmptyMC_ByModule = new TH1D ("TBVerticesEmptyMC_ByModule", "TBVerticesEmptyMC_ByModule", vertexBins.size()-1, &vertexBins[0]);
  TH1D *ANNVerticesEmptyData_ByModule = new TH1D ("ANNVerticesEmptyData_ByModule", "ANNVerticesEmptyData_ByModule", vertexBins.size()-1, &vertexBins[0]);
  TH1D *TBVerticesEmptyData_ByModule = new TH1D ("TBVerticesEmptyData_ByModule", "TBVerticesEmptyData_ByModule", vertexBins.size()-1, &vertexBins[0]);
  TH1D *TruthVerticesEmptyMC_ByModule = new TH1D ("TruthVerticesEmptyMC_ByModule", "TruthVerticesEmptyMC_ByModule", vertexBins.size()-1, &vertexBins[0]);
  TH1D *ANNVerticesEmptyMC_ByZPos = new TH1D ("ANNVerticesEmptyMC_ByZPos", "ANNVerticesEmptyMC_ByZPos", 3400, 4200, 5900);
  TH1D *TBVerticesEmptyMC_ByZPos = new TH1D ("TBVerticesEmptyMC_ByZPos", "TBVerticesEmptyMC_ByModuleTruthVerticesEmptyMC_ByZPosHighRes", 3400, 4200, 5900);
  TH1D *ANNVerticesEmptyData_ByZPos = new TH1D ("ANNVerticesEmptyData_ByZPos", "ANNVerticesEmptyData_ByZPos", 3400, 4200, 5900);
  TH1D *TBVerticesEmptyData_ByZPos = new TH1D ("TBVerticesEmptyData_ByZPos", "TBVerticesEmptyData_ByZPos", 3400, 4200, 5900);
  TH1D *TruthVerticesEmptyMC_ByZPos = new TH1D ("TruthVerticesEmptyMC_ByZPos", "TruthVerticesEmptyMC_ByZPos", 3400, 4200, 5900);
  //ERecoil
  TH2D *TruthVerticesEmptyMCERecoil_ByModule = new TH2D ("TruthVerticesEmptyMCERecoil_ByModule", "TruthVerticesEmptyMCERecoil_ByModule", vertexBins.size()-1, &vertexBins[0], 100, 0, 50);
  TH2D *ANNVerticesEmptyMCERecoil_ByModule = new TH2D ("ANNVerticesEmptyMCERecoil_ByModule", "ANNVerticesEmptyMCERecoil_ByModule", vertexBins.size()-1, &vertexBins[0], 100, 0, 50);
  TH2D *TBVerticesEmptyMCERecoil_ByModule = new TH2D ("TBVerticesEmptyMCERecoil_ByModule", "TBVerticesEmptyMCERecoil_ByModule", vertexBins.size()-1, &vertexBins[0], 100, 0, 50);
  TH2D *ANNVerticesEmptyDataERecoil_ByModule = new TH2D ("ANNVerticesEmptyDataERecoil_ByModule", "ANNVerticesEmptyDataERecoil_ByModule", vertexBins.size()-1, &vertexBins[0], 100, 0, 50);
  TH2D *TBVerticesEmptyDataERecoil_ByModule = new TH2D ("TBVerticesEmptyDataERecoil_ByModule", "TBVerticesEmptyDataERecoil_ByModule", vertexBins.size()-1, &vertexBins[0], 100, 0, 50);
  //Multiplicity
  TH2D *TruthVerticesEmptyMCMultiplicity_ByModule = new TH2D ("TruthVerticesEmptyMCMultiplicity_ByModule", "TruthVerticesEmptyMCMultiplicity_ByModule", vertexBins.size()-1, &vertexBins[0], 10, 0, 10);
  TH2D *ANNVerticesEmptyMCMultiplicity_ByModule = new TH2D ("ANNVerticesEmptyMCMultiplicity_ByModule", "ANNVerticesEmptyMCMultiplicity_ByModule", vertexBins.size()-1, &vertexBins[0], 10, 0, 10);
  TH2D *TBVerticesEmptyMCMultiplicity_ByModule = new TH2D ("TBVerticesEmptyMCMultiplicity_ByModule", "TBVerticesEmptyMCMultiplicity_ByModule", vertexBins.size()-1, &vertexBins[0], 10, 0, 10);
  TH2D *ANNVerticesEmptyDataMultiplicity_ByModule = new TH2D ("ANNVerticesEmptyDataMultiplicity_ByModule", "ANNVerticesEmptyDataMultiplicity_ByModule", vertexBins.size()-1, &vertexBins[0], 10, 0, 10);
  TH2D *TBVerticesEmptyDataMultiplicity_ByModule = new TH2D ("TBVerticesEmptyDataMultiplicity_ByModule", "TBVerticesEmptyDataMultiplicity_ByModule", vertexBins.size()-1, &vertexBins[0], 10, 0, 10);



  /* TH1D *ANNVerticesFilledMC_ByZPosQ1 = new TH1D ("ANNVerticesFilledMC_ByZPosQ1", "ANNVerticesFilledMC_ByZPosQ1", 3400, 4200, 5900);
  TH1D *TBVerticesFilledMC_ByZPosQ1 = new TH1D ("TBVerticesFilledMC_ByZPosQ1", "TBVerticesFilledMC_ByZPosQ1", 3400, 4200, 5900);
  TH1D *ANNVerticesFilledData_ByZPosQ1 = new TH1D ("ANNVerticesFilledData_ByZPosQ1", "ANNVerticesFilledData_ByZPosQ1", 3400, 4200, 5900);
  TH1D *TBVerticesFilledData_ByZPosQ1 = new TH1D ("TBVerticesFilledData_ByZPosQ1", "TBVerticesFilledData_ByZPosQ1", 3400, 4200, 5900);
  TH1D *TruthVerticesFilledMC_ByZPosQ1 = new TH1D ("TruthVerticesFilledMC_ByZPosQ1", "TruthVerticesFilledMC_ByZPosQ1", 3400, 4200, 5900);
  TH1D *ANNVerticesFilledMC_ByZPosQ2 = new TH1D ("ANNVerticesFilledMC_ByZPosQ2", "ANNVerticesFilledMC_ByZPosQ2", 3400, 4200, 5900);
  TH1D *TBVerticesFilledMC_ByZPosQ2 = new TH1D ("TBVerticesFilledMC_ByZPosQ2", "TBVerticesFilledMC_ByZPosQ2", 3400, 4200, 5900);
  TH1D *ANNVerticesFilledData_ByZPosQ2 = new TH1D ("ANNVerticesFilledData_ByZPosQ2", "ANNVerticesFilledData_ByZPosQ2", 3400, 4200, 5900);
  TH1D *TBVerticesFilledData_ByZPosQ2 = new TH1D ("TBVerticesFilledData_ByZPosQ2", "TBVerticesFilledData_ByZPosQ2", 3400, 4200, 5900);
  TH1D *TruthVerticesFilledMC_ByZPosQ2 = new TH1D ("TruthVerticesFilledMC_ByZPosQ2", "TruthVerticesFilledMC_ByZPosQ2", 3400, 4200, 5900);
  TH1D *ANNVerticesFilledMC_ByZPosQ3 = new TH1D ("ANNVerticesFilledMC_ByZPosQ3", "ANNVerticesFilledMC_ByZPosQ3", 3400, 4200, 5900);
  TH1D *TBVerticesFilledMC_ByZPosQ3 = new TH1D ("TBVerticesFilledMC_ByZPosQ3", "TBVerticesFilledMC_ByZPosQ3", 3400, 4200, 5900);
  TH1D *ANNVerticesFilledData_ByZPosQ3 = new TH1D ("ANNVerticesFilledData_ByZPosQ3", "ANNVerticesFilledData_ByZPosQ3", 3400, 4200, 5900);
  TH1D *TBVerticesFilledData_ByZPosQ3 = new TH1D ("TBVerticesFilledData_ByZPosQ3", "TBVerticesFilledData_ByZPosQ3", 3400, 4200, 5900);
  TH1D *TruthVerticesFilledMC_ByZPosQ3 = new TH1D ("TruthVerticesFilledMC_ByZPosQ3", "TruthVerticesFilledMC_ByZPosQ3", 3400, 4200, 5900);
  TH1D *ANNVerticesFilledMC_ByZPosQ4 = new TH1D ("ANNVerticesFilledMC_ByZPosQ4", "ANNVerticesFilledMC_ByZPosQ4", 3400, 4200, 5900);
  TH1D *TBVerticesFilledMC_ByZPosQ4 = new TH1D ("TBVerticesFilledMC_ByZPosQ4", "TBVerticesFilledMC_ByZPosQ4", 3400, 4200, 5900);
  TH1D *ANNVerticesFilledData_ByZPosQ4 = new TH1D ("ANNVerticesFilledData_ByZPosQ4", "ANNVerticesFilledData_ByZPosQ4", 3400, 4200, 5900);
  TH1D *TBVerticesFilledData_ByZPosQ4 = new TH1D ("TBVerticesFilledData_ByZPosQ4", "TBVerticesFilledData_ByZPosQ4", 3400, 4200, 5900);
  TH1D *TruthVerticesFilledMC_ByZPosQ4 = new TH1D ("TruthVerticesFilledMC_ByZPosQ4", "TruthVerticesFilledMC_ByZPosQ4", 3400, 4200, 5900); */
  //==============================================================================
  // End - Declaring histograms
  //==============================================================================


  // Loop entries and fill
  try
  {
    std::map< std::string, std::vector<CVUniverse*> > error_bandsFilled;
    std::map<std::string, std::vector<CVUniverse*> > band_fluxFilled = PlotUtils::GetFluxSystematicsMap<CVUniverse>(optionsFilled.m_mc, CVUniverse::GetNFluxUniverses());
    //error_bandsFilled.insert(band_fluxFilled.begin(), band_fluxFilled.end()); //Necessary to get flux integral later...
    error_bandsFilled["cv"] = {new CVUniverse(optionsFilled.m_mc)};
    CVUniverse* data_universFilled = new CVUniverse(optionsFilled.m_data);
    std::vector<CVUniverse*> data_bandFilled = {data_universFilled};
    CVUniverse::SetTruth(false);
    LoopAndFillMC(optionsFilled.m_mc, error_bandsFilled, mycuts, model, *ANNVerticesFilledMC_ByModule, *ANNVerticesFilledMC_ByZPos, *TBVerticesFilledMC_ByModule, *TBVerticesFilledMC_ByZPos, *TruthVerticesFilledMC_ByModule, *TruthVerticesFilledMC_ByZPos, *ANNVerticesFilledMCERecoil_ByModule, *ANNVerticesFilledMCMultiplicity_ByModule, *TBVerticesFilledMCERecoil_ByModule, *TBVerticesFilledMCMultiplicity_ByModule, *TruthVerticesFilledMCERecoil_ByModule, *TruthVerticesFilledMCMultiplicity_ByModule);
    mycuts.resetStats();
    CVUniverse::SetTruth(false);
    LoopAndFillData(optionsFilled.m_data, data_bandFilled, mycuts, *ANNVerticesFilledData_ByModule, *ANNVerticesFilledData_ByZPos, *TBVerticesFilledData_ByZPos, *TBVerticesFilledData_ByModule, *ANNVerticesFilledDataERecoil_ByModule, *TBVerticesFilledDataERecoil_ByModule, *ANNVerticesFilledDataMultiplicity_ByModule, *TBVerticesFilledDataMultiplicity_ByModule);
    //std::cout << "Data cut summary:\n" << mycuts << "\n";
    
    std::map< std::string, std::vector<CVUniverse*> > error_bandsEmpty;
    std::map<std::string, std::vector<CVUniverse*> > band_fluxEmpty = PlotUtils::GetFluxSystematicsMap<CVUniverse>(optionsEmpty.m_mc, CVUniverse::GetNFluxUniverses());
    //error_bandsEmpty.insert(band_fluxEmpty.begin(), band_fluxEmpty.end()); //Necessary to get flux integral later...
    error_bandsEmpty["cv"] = {new CVUniverse(optionsFilled.m_mc)};
    CVUniverse* data_universeEmpty = new CVUniverse(optionsFilled.m_data);
    std::vector<CVUniverse*> data_bandEmpty = {data_universeEmpty};
    CVUniverse::SetTruth(false);
    LoopAndFillMC(optionsEmpty.m_mc, error_bandsEmpty, mycuts, model, *ANNVerticesEmptyMC_ByModule, *ANNVerticesEmptyMC_ByZPos, *TBVerticesEmptyMC_ByModule, *TBVerticesEmptyMC_ByZPos, *TruthVerticesEmptyMC_ByModule, *TruthVerticesEmptyMC_ByZPos, *ANNVerticesEmptyMCERecoil_ByModule, *ANNVerticesEmptyMCMultiplicity_ByModule, *TBVerticesEmptyMCERecoil_ByModule, *TBVerticesEmptyMCMultiplicity_ByModule, *TruthVerticesEmptyMCERecoil_ByModule, *TruthVerticesEmptyMCMultiplicity_ByModule);
    mycuts.resetStats();
    CVUniverse::SetTruth(false);
    LoopAndFillData(optionsEmpty.m_data, data_bandEmpty, mycuts, *ANNVerticesEmptyData_ByModule, *ANNVerticesEmptyData_ByZPos, *TBVerticesEmptyData_ByZPos, *TBVerticesEmptyData_ByModule, *ANNVerticesEmptyDataERecoil_ByModule, *TBVerticesEmptyDataERecoil_ByModule, *ANNVerticesEmptyDataMultiplicity_ByModule, *TBVerticesEmptyDataMultiplicity_ByModule);


    TFile* OutDir = TFile::Open("VertexValidations.root", "RECREATE");
    //Filled
    ANNVerticesFilledMC_ByModule->Write();
    TBVerticesFilledMC_ByModule->Write();
    ANNVerticesFilledData_ByModule->Write();
    TBVerticesFilledData_ByModule->Write();
    TruthVerticesFilledMC_ByModule->Write();
    ANNVerticesFilledMC_ByZPos->Write();
    TBVerticesFilledMC_ByZPos->Write();
    ANNVerticesFilledData_ByZPos->Write();
    TBVerticesFilledData_ByZPos->Write();
    TruthVerticesFilledMC_ByZPos->Write();
    //ERecoil
    TruthVerticesFilledMCERecoil_ByModule->Write();
    ANNVerticesFilledMCERecoil_ByModule->Write();
    TBVerticesFilledMCERecoil_ByModule->Write();
    ANNVerticesFilledDataERecoil_ByModule->Write();
    TBVerticesFilledDataERecoil_ByModule->Write();
    //Multiplicity
    TruthVerticesFilledMCMultiplicity_ByModule->Write();
    ANNVerticesFilledMCMultiplicity_ByModule->Write();
    TBVerticesFilledMCMultiplicity_ByModule->Write();
    ANNVerticesFilledDataMultiplicity_ByModule->Write();
    TBVerticesFilledDataMultiplicity_ByModule->Write();

    //Empty
    ANNVerticesEmptyMC_ByModule->Write();
    TBVerticesEmptyMC_ByModule->Write();
    ANNVerticesEmptyData_ByModule->Write();
    TBVerticesEmptyData_ByModule->Write();
    TruthVerticesEmptyMC_ByModule->Write();
    ANNVerticesEmptyMC_ByZPos->Write();
    TBVerticesEmptyMC_ByZPos->Write();
    ANNVerticesEmptyData_ByZPos->Write();
    TBVerticesEmptyData_ByZPos->Write();
    TruthVerticesEmptyMC_ByZPos->Write();
    //ERecoil
    TruthVerticesEmptyMCERecoil_ByModule->Write();
    ANNVerticesEmptyMCERecoil_ByModule->Write();
    TBVerticesEmptyMCERecoil_ByModule->Write();
    ANNVerticesEmptyDataERecoil_ByModule->Write();
    TBVerticesEmptyDataERecoil_ByModule->Write();
    //Multiplicity
    TruthVerticesEmptyMCMultiplicity_ByModule->Write();
    ANNVerticesEmptyMCMultiplicity_ByModule->Write();
    TBVerticesEmptyMCMultiplicity_ByModule->Write();
    ANNVerticesEmptyDataMultiplicity_ByModule->Write();
    TBVerticesEmptyDataMultiplicity_ByModule->Write();

    auto EmptyMCPOT = new TParameter<double>("EmptyMCPOT", optionsEmpty.m_mc_pot);
    EmptyMCPOT->Write();
    auto EmptyDataPOT = new TParameter<double>("EmptyDataPOT", optionsEmpty.m_data_pot);
    EmptyDataPOT->Write();
    auto FilledMCPOT = new TParameter<double>("FilledMCPOT", optionsFilled.m_mc_pot);
    FilledMCPOT->Write();
    auto FilledDataPOT = new TParameter<double>("FilledDataPOT", optionsFilled.m_data_pot);
    FilledDataPOT->Write();
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
  catch(const ROOT::exception& e)
  {
    std::cerr << "Ending on a ROOT error message.  No histograms will be produced.\n"
              << "If the message talks about \"TNetXNGFile\", this could be a problem with dCache.  The message is:\n"
              << e.what() << "\n" << USAGE << "\n";
    return badFileRead;
  }

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