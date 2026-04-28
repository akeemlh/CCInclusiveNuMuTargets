#define MC_OUT_FILE_NAME_BASE "runEventLoopTrackerMC"
#define DATA_OUT_FILE_NAME_BASE "runEventLoopTrackerData"
#define MIGRATION_2D_OUT_FILE_NAME_BASE "runEventLoopTracker2DMigration"

#define USAGE \
"\n*** USAGE ***\n"\
"runEventLoopTracker <dataPlaylist.txt> <mcPlaylist.txt> <optional petal(0-11)>\n\n"\
"To run eventLoop without petals being consideration - i.e if you're not planning \n"\
"to perform a daisy reweigt then pass -1 as the petal\n"\
"If petal is not provided, default behaviour is to run over everything, this can take a while\n"\
"Additionally, if running on the grid, it is possibly more time efficient to only run 1 or a small number of\n"\
"petals at at time and distribute the petals you're interested in over several nodes\n\n"\
"*** Explanation ***\n"\
"Reduce MasterAnaDev AnaTuples to event selection histograms to extract a\n"\
"single-differential inclusive cross section for the 2021 MINERvA 101 tutorial.\n\n"\
"*** The Input Files ***\n"\
"Playlist files are plaintext files with 1 file name per line.  Filenames may be\n"\
"xrootd URLs or refer to the local filesystem.  The first playlist file's\n"\
"entries will be treated like data, and the second playlist's entries must\n"\
"have the \"Truth\" tree to use for calculating the efficiency denominator.\n\n"\
"*** Output ***\n"\
"Produces three files with the base name, " MC_OUT_FILE_NAME_BASE ", " DATA_OUT_FILE_NAME_BASE " and\n"\
" " MIGRATION_2D_OUT_FILE_NAME_BASE " all histograms needed for the ExtractCrossSection program also built by this\n"\
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
#include "PlotUtils/LowQ2PiReweighter.h"
#include "PlotUtils/AMUDISReweighter.h"
#include "PlotUtils/SuSAFromValencia2p2hReweighter.h"
#include "PlotUtils/FSIReweighter.h"
#include "PlotUtils/TargetUtils.h"

#include "util/COHPionReweighter.h"
#include "util/AsitsRobFineFluxReweighter.h"
#include "util/DiffractiveReweighter.h"
#include "PlotUtils/BodekRitchieReweighter.h"
#include "PlotUtils/PionReweighter.h"
#include "util/Q2TpiReweighter.h"

#include "util/NukeUtils.h"
#pragma GCC diagnostic pop

//ROOT includes
#include "TParameter.h"

#include <chrono>

#include "Math/Vector3D.h"
#include "TH3D.h"

//c++ includes
#include <iostream>
#include <cstdlib> //getenv()
#include <fstream>
#include <sstream> //reading input files
#include <deque>


#include "utilities/SMSCalibrator.h"

bool useCalibrator = true;

//These 2 variables are used when doing broken down runs to speed up running on the grid
int nSubruns = 0;
int nProcess = 0;

//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFillEventSelection(
    PlotUtils::ChainWrapper* chain,
    std::map<std::string, std::vector<CVUniverse*> > error_bands,
    std::vector<std::vector<std::unique_ptr<Variable>>> &vars_vec,
    std::vector<std::vector<std::unique_ptr<Variable2D>>> &vars2D_vec,
    std::vector<std::vector<Study *>> &studies_vec,
    std::vector<PlotUtils::Cutter<CVUniverse, MichelEvent>*> &michelcuts_vec,
    PlotUtils::Model<CVUniverse, MichelEvent> &model,
    std::vector<int> &inpetal_vec,
    SMSCalibrator& calibrator)
{
  //calibrator.cacheEvt = -1;
  const int WINDOW_SIZE = 10; // number of samples to average over
  std::deque<std::pair<int, std::chrono::steady_clock::time_point>> timeWindow;
  std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();  assert(!error_bands["cv"].empty() && "\"cv\" error band is empty!  Can't set Model weight.");
  assert(!error_bands["cv"].empty() && "\"cv\" error band is empty!  Can't set Model weight.");
  auto& cvUniv = error_bands["cv"].front();

  std::cout << "Starting MC reco loop...\n";
  const int nEntries = chain->GetEntries();
  //const int nEntries = 10000;
  //for (int i=19598; i<nEntries; ++i)
  PlotUtils::TargetUtils util;
  for (int i = 0; i < nEntries; ++i)
  {
    if (i % 1000 == 0 && i > 0)
    {
        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();

        // Add current sample to the sliding window
        timeWindow.push_back({i, now});
        if (timeWindow.size() > WINDOW_SIZE)
            timeWindow.pop_front();

        // Overall elapsed time
        double msElapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - startTime).count();
        double sElapsed   = msElapsed / 1000.0;

        // Dynamic rate: events/sec over the sliding window
        double eventsPerSec = 0.0;
        double msRemaining  = 0.0;
        if (timeWindow.size() >= 2)
        {
            auto& oldest = timeWindow.front();
            double msWindow  = std::chrono::duration_cast<std::chrono::milliseconds>(now - oldest.second).count();
            int    evWindow  = i - oldest.first;
            eventsPerSec = (msWindow > 0) ? (eventsPerSec = evWindow / (msWindow / 1000.0)) : 0.0;
            msRemaining  = (eventsPerSec > 0) ? ((nEntries - i) / eventsPerSec * 1000.0) : 0.0;
        }
        else
        {
            // Fallback to overall average until window is populated
            double msPerEvent = msElapsed / i;
            msRemaining = (nEntries - i) * msPerEvent;
        }

        int sRemaining = static_cast<int>(msRemaining / 1000.0);

        // Rate indicator: pick a sensible unit
        std::string rateStr;
        if (eventsPerSec >= 1000.0)
            rateStr = std::to_string(static_cast<int>(eventsPerSec / 1000.0)) + "k ev/s";
        else
            rateStr = std::to_string(static_cast<int>(eventsPerSec)) + " ev/s";

        std::cout << i << " / " << nEntries
                  << "  Elapsed: "   << std::fixed << std::setprecision(1) << sElapsed  << "s"
                  << "  Remaining: " << sRemaining << "s"
                  << "  Rate: "      << rateStr
                  << "          \r"  << std::flush;
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
        if (useCalibrator) universe->SetCalibrator(&calibrator);
        univCount++; // Put the iterator right at the start so it's executed even in paths that lead to a continue, don't forget to subtract by 1 when we use it
        MichelEvent myevent; // make sure your event is inside the error band loop. 
    
        // Tell the Event which entry in the TChain it's looking at
        universe->SetEntry(i);
         
        // This is where you would Access/create a Michel

        for (int pet = 0; pet < inpetal_vec.size(); pet++) 
        {
          std::vector<std::unique_ptr<Variable>>& vars=vars_vec[pet];
          std::vector<std::unique_ptr<Variable2D>>& vars2D = vars2D_vec[pet];
          std::vector<Study *>& studies = studies_vec[pet];
          PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts = *michelcuts_vec[pet];
          int inpetal = inpetal_vec[pet];


          int intType = universe->GetInteractionType();
          if (intType == 3 && !util::IsTrueDIS(*universe)) intType = 30; //SoftDIS
          //weight is ignored in isMCSelected() for all but the CV Universe.
          if (!michelcuts.isMCSelected(*universe, myevent, cvWeight).all()) continue; //all is another function that will later help me with sidebands
          const double weight = model.GetWeight(*universe, myevent); //Only calculate the per-universe weight for events that will actually use it.
          for(auto& var: vars)
          { 
            var->selectedMCReco->FillUniverse(universe, var->GetRecoValue(*universe), weight); //"Fake data" for closure
            (*var->m_intChannels)[intType].FillUniverse(universe, var->GetRecoValue(*universe), weight);
          }
          for(auto& var: vars2D)
          {
            var->selectedMCReco->FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight); //"Fake data" for closure
            (*var->m_intChannels)[intType].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
            //(*var->m_interactionChannels)[intType].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
          }
          const bool isSignal = michelcuts.isSignal(*universe, weight);
          if(isSignal)
          {
            for(auto& study: studies) study->SelectedSignal(*universe, myevent, weight);

            for(auto& var: vars)
            {
              //Cross section components
              var->efficiencyNumerator->FillUniverse(universe, var->GetTrueValue(*universe), weight);
              var->migration->FillUniverse(universe, var->GetRecoValue(*universe), var->GetTrueValue(*universe), weight);
              var->selectedSignalReco->FillUniverse(universe, var->GetRecoValue(*universe), weight); //Efficiency numerator in reco variables.  Useful for warping studies.
            }
            for(auto& var: vars2D)
            {
              var->efficiencyNumerator->FillUniverse(universe, var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
              var->migration->Fill(var->GetRecoValueX(*universe), var->GetRecoValueY(*universe),var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
              var->selectedSignalReco->FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight); //Efficiency numerator in reco variables.  Useful for warping studies.
            }
          }
          else
          {
            int bkgd_ID = -1;
            if (universe->GetCurrent() == 2) bkgd_ID = 0;
            else if (universe->GetTruthNuPDG() == -(universe->GetAnalysisNuPDG())) bkgd_ID = 1;
            for(auto& var: vars)
            {
              (*var->m_backgroundHists)[bkgd_ID].FillUniverse(universe, var->GetRecoValue(*universe), weight);
            }
            for(auto& var: vars2D)
            {
              (*var->m_backgroundHists)[bkgd_ID].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
            } 
          }
        } //End petal loop
      } // End band's universe loop
    } // End Band loop
  } //End entries loop
  std::cout << "Finished MC reco loop.\n";
}

void LoopAndFillData( PlotUtils::ChainWrapper* data,
    std::vector<CVUniverse*> data_band,
    std::vector<std::vector<std::unique_ptr<Variable>>> &vars_vec,
    std::vector<std::vector<std::unique_ptr<Variable2D>>> &vars2D_vec,
    std::vector<std::vector<Study *>> &studies_vec,
    std::vector<PlotUtils::Cutter<CVUniverse, MichelEvent>*> &michelcuts_vec,
    std::vector<int> &inpetal_vec,
    SMSCalibrator& calibrator)
{
  //calibrator.cacheEvt = -1;
  std::cout << "Starting data loop...\n";
  const int WINDOW_SIZE = 10; // number of samples to average over
  std::deque<std::pair<int, std::chrono::steady_clock::time_point>> timeWindow;
  std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
  const int nEntries = data->GetEntries();
  PlotUtils::TargetUtils util;
  //const int nEntries = 10000;
  for (int i = 0; i < nEntries; ++i)
  {
    for (auto universe : data_band) {
      universe->SetEntry(i);
      if (useCalibrator) universe->SetCalibrator(&calibrator);
      if (i % 1000 == 0 && i > 0)
      {
          std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();

          // Add current sample to the sliding window
          timeWindow.push_back({i, now});
          if (timeWindow.size() > WINDOW_SIZE)
              timeWindow.pop_front();

          // Overall elapsed time
          double msElapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - startTime).count();
          double sElapsed   = msElapsed / 1000.0;

          // Dynamic rate: events/sec over the sliding window
          double eventsPerSec = 0.0;
          double msRemaining  = 0.0;
          if (timeWindow.size() >= 2)
          {
              auto& oldest = timeWindow.front();
              double msWindow  = std::chrono::duration_cast<std::chrono::milliseconds>(now - oldest.second).count();
              int    evWindow  = i - oldest.first;
              eventsPerSec = (msWindow > 0) ? (eventsPerSec = evWindow / (msWindow / 1000.0)) : 0.0;
              msRemaining  = (eventsPerSec > 0) ? ((nEntries - i) / eventsPerSec * 1000.0) : 0.0;
          }
          else
          {
              // Fallback to overall average until window is populated
              double msPerEvent = msElapsed / i;
              msRemaining = (nEntries - i) * msPerEvent;
          }

          int sRemaining = static_cast<int>(msRemaining / 1000.0);

          // Rate indicator: pick a sensible unit
          std::string rateStr;
          if (eventsPerSec >= 1000.0)
              rateStr = std::to_string(static_cast<int>(eventsPerSec / 1000.0)) + "k ev/s";
          else
              rateStr = std::to_string(static_cast<int>(eventsPerSec)) + " ev/s";

          std::cout << i << " / " << nEntries
                    << "  Elapsed: "   << std::fixed << std::setprecision(1) << sElapsed  << "s"
                    << "  Remaining: " << sRemaining << "s"
                    << "  Rate: "      << rateStr
                    << "          \r"  << std::flush;
      }
      MichelEvent myevent; 
      for (int pet = 0; pet < inpetal_vec.size(); pet++) 
      {
        std::vector<std::unique_ptr<Variable>>& vars=vars_vec[pet];
        std::vector<std::unique_ptr<Variable2D>>& vars2D = vars2D_vec[pet];
        std::vector<Study *>& studies = studies_vec[pet];
        PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts = *michelcuts_vec[pet];
        int inpetal = inpetal_vec[pet];

        if (!michelcuts.isDataSelected(*universe, myevent).all()) continue;
        for(auto& study: studies) study->Selected(*universe, myevent, 1); 
        //std::cout<<"petal: " <<petal << "X: " << universe->GetANNVertex().X() << "\tY: "<<universe->GetANNVertex().Y()<<std::endl;
        for(auto& var: vars)
        {
          var->dataHist->FillUniverse(universe, var->GetRecoValue(*universe, myevent.m_idx), 1);
        }

        for(auto& var: vars2D)
        {
          var->dataHist->FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), 1);
        }
      } //End of petal loop
    }
  }
  std::cout << "Finished data loop.\n";
}

void LoopAndFillEffDenom( PlotUtils::ChainWrapper* truth,
    std::map<std::string, std::vector<CVUniverse*> > truth_bands,
    std::vector<std::vector<std::unique_ptr<Variable>>> &vars_vec,
    std::vector<std::vector<std::unique_ptr<Variable2D>>> &vars2D_vec,
    std::vector<PlotUtils::Cutter<CVUniverse, MichelEvent>*> &michelcuts_vec,
    PlotUtils::Model<CVUniverse, MichelEvent> &model,
    std::vector<int> &inpetal_vec)
{
  const int WINDOW_SIZE = 10; // number of samples to average over
  std::deque<std::pair<int, std::chrono::steady_clock::time_point>> timeWindow;
  std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
  assert(!truth_bands["cv"].empty() && "\"cv\" error band is empty!  Could not set Model entry.");
  auto& cvUniv = truth_bands["cv"].front();
  PlotUtils::TargetUtils util;
  std::cout << "Starting efficiency denominator loop...\n";
  const int nEntries = truth->GetEntries();
  //const int nEntries = 10000;

  for (int i = 0; i < nEntries; ++i)
  {
    if (i % 1000 == 0 && i > 0)
    {
        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();

        // Add current sample to the sliding window
        timeWindow.push_back({i, now});
        if (timeWindow.size() > WINDOW_SIZE)
            timeWindow.pop_front();

        // Overall elapsed time
        double msElapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - startTime).count();
        double sElapsed   = msElapsed / 1000.0;

        // Dynamic rate: events/sec over the sliding window
        double eventsPerSec = 0.0;
        double msRemaining  = 0.0;
        if (timeWindow.size() >= 2)
        {
            auto& oldest = timeWindow.front();
            double msWindow  = std::chrono::duration_cast<std::chrono::milliseconds>(now - oldest.second).count();
            int    evWindow  = i - oldest.first;
            eventsPerSec = (msWindow > 0) ? (eventsPerSec = evWindow / (msWindow / 1000.0)) : 0.0;
            msRemaining  = (eventsPerSec > 0) ? ((nEntries - i) / eventsPerSec * 1000.0) : 0.0;
        }
        else
        {
            // Fallback to overall average until window is populated
            double msPerEvent = msElapsed / i;
            msRemaining = (nEntries - i) * msPerEvent;
        }

        int sRemaining = static_cast<int>(msRemaining / 1000.0);

        // Rate indicator: pick a sensible unit
        std::string rateStr;
        if (eventsPerSec >= 1000.0)
            rateStr = std::to_string(static_cast<int>(eventsPerSec / 1000.0)) + "k ev/s";
        else
            rateStr = std::to_string(static_cast<int>(eventsPerSec)) + " ev/s";

        std::cout << i << " / " << nEntries
                  << "  Elapsed: "   << std::fixed << std::setprecision(1) << sElapsed  << "s"
                  << "  Remaining: " << sRemaining << "s"
                  << "  Rate: "      << rateStr
                  << "          \r"  << std::flush;
    }

    MichelEvent cvEvent;
    cvUniv->SetEntry(i);
    model.SetEntry(*cvUniv, cvEvent);
    const double cvWeight = model.GetWeight(*cvUniv, cvEvent);

    //=========================================
    // Systematics loop(s)
    //=========================================
    for (auto band : truth_bands)
    {
      std::vector<CVUniverse*> truth_band_universes = band.second;
      for (auto universe : truth_band_universes)
      {
        MichelEvent myevent; //Only used to keep the Model happy

        // Tell the Event which entry in the TChain it's looking at
        universe->SetEntry(i);

        for (int pet = 0; pet < inpetal_vec.size(); pet++) 
        {
          std::vector<std::unique_ptr<Variable>>& vars=vars_vec[pet];
          std::vector<std::unique_ptr<Variable2D>>& vars2D = vars2D_vec[pet];
          PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts = *michelcuts_vec[pet];
          int inpetal = inpetal_vec[pet];
          int intType = universe->GetInteractionType();
          if (intType == 3 && !util::IsTrueDIS(*universe)) intType = 30; //SoftDIS

          if (!michelcuts.isEfficiencyDenom(*universe, cvWeight)) continue; //Weight is ignored for isEfficiencyDenom() in all but the CV universe 
          const double weight = model.GetWeight(*universe, myevent); //Only calculate the weight for events that will use it
          for(auto& var: vars)
          {
            var->efficiencyDenominator->FillUniverse(universe, var->GetTrueValue(*universe), weight);
            (*var->m_intChannelsEffDenom)[intType].FillUniverse(universe, var->GetTrueValue(*universe), weight);
          }

          for(auto& var: vars2D)
          {
            var->efficiencyDenominator->FillUniverse(universe, var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
            (*var->m_intChannelsEffDenom)[intType].FillUniverse(universe, var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
          }
        } //End petal loop
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
  std::cout<<"Here2\n";
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


std::string breakUpInputFileList(std::string inputFile, int numTotal, int thisNum, int& numFiles) {
    std::ifstream in(inputFile);
    if (!in.is_open()) {
        std::cerr << "Error opening input file: " << inputFile << std::endl;
        return "";
    }

    // Read all lines into memory
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
      if (line.length()>0)
      {
        lines.push_back(line);
      }
    }
    in.close();

    size_t totalLines = lines.size();
    if (totalLines == 0) {
        std::cerr << "Input file is empty!" << std::endl;
        return "";
    }

    // Split into chunks
    size_t chunkSize = std::floor(totalLines / numTotal);
    if (chunkSize < 1) chunkSize = 1;
    std::string outfile = inputFile+"tmp"+std::to_string(thisNum)+".txt";
    std::ofstream out(outfile);
    if (!out.is_open()) {
        std::cerr << "Error opening output file: " << outfile << std::endl;
        return "";
    }
    int start = thisNum * chunkSize;
    int end = (thisNum +1)* chunkSize;
    if (thisNum == numTotal-1 ) end = lines.size();
    int nFilesTmp = 0;
    for (size_t j = start; (j < end && j < totalLines); ++j) {
        out << lines[j] << "\n";
        nFilesTmp++;
    }
    numFiles = nFilesTmp;
    out.close();

    return outfile;
}

/* 
// Extracts the filename from a full directory path or xrootd URL
std::string getBaseName(const std::string& path) {
    size_t pos = path.find_last_of("/\\");
    return (pos == std::string::npos) ? path : path.substr(pos + 1);
}

// Trims whitespace and hidden carriage returns (CRITICAL for shell commands)
std::string trimLine(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (std::string::npos == first) return "";
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, (last - first + 1));
}

std::string breakUpInputFileList(std::string inputFile, int numTotal, int thisNum, int& numFiles) {
    std::ifstream in(inputFile);
    if (!in.is_open()) {
        std::cerr << "Error opening input file: " << inputFile << std::endl;
        return "";
    }

    // Read and trim all lines into memory
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        std::string trimmedLine = trimLine(line);
        if (trimmedLine.length() > 0) {
            lines.push_back(trimmedLine);
        }
    }
    in.close();

    size_t totalLines = lines.size();
    if (totalLines == 0) {
        std::cerr << "Input file is empty!" << std::endl;
        return "";
    }

    // Determine the local scratch directory for grid jobs
    const char* scratchDirPtr = std::getenv("_CONDOR_SCRATCH_DIR");
    std::string scratchDir = scratchDirPtr ? std::string(scratchDirPtr) : "."; // Fallback to current dir if run locally

    // Split into chunks
    size_t chunkSize = std::floor(totalLines / numTotal);
    if (chunkSize < 1) chunkSize = 1;
    
    std::string outfile = inputFile + "tmp" + std::to_string(thisNum) + ".txt";
    std::ofstream out(outfile);
    if (!out.is_open()) {
        std::cerr << "Error opening output file: " << outfile << std::endl;
        return "";
    }

    int start = thisNum * chunkSize;
    int end = (thisNum + 1) * chunkSize;
    if (thisNum == numTotal - 1) end = lines.size();
    
    int nFilesTmp = 0;
    
    for (size_t j = start; (j < end && j < totalLines); ++j) {
        std::string remotePath = lines[j];
        std::string baseName = getBaseName(remotePath);
        std::string localPath = scratchDir + "/" + baseName;

        // Construct the ifdh cp command
        std::string cmd = "ifdh cp " + remotePath + " " + localPath;
        std::cout << "Executing: " << cmd << std::endl;

        // Execute the copy
        int retCode = std::system(cmd.c_str());

        if (retCode == 0) {
            // Copy succeeded, write the local path to the playlist
            out << localPath << "\n";
            std::cout << "Successfully staged file locally: " << localPath << std::endl;
        } else {
            // Copy failed, fallback to the remote path so the job doesn't crash entirely
            std::cerr << "Warning: ifdh cp failed (Return Code: " << retCode << ") for " << remotePath << "\n";
            std::cerr << "Falling back to streaming remote path directly..." << std::endl;
            out << remotePath << "\n";
        }
        
        nFilesTmp++;
    }
    
    numFiles = nFilesTmp;
    out.close();

    return outfile;
} */

//==============================================================================
// Main
//==============================================================================
int main(const int argc, const char** argv)
{
  TH1::AddDirectory(false);

  //Validate input.
  if(argc < 3)
  {
    std::cerr << "Expected 3 or more arguments, but got " << argc - 1 << "\n" << USAGE << "\n";
    return badCmdLine;
  } 

  std::cout<<"Here1\n";
  //One playlist must contain only MC files, and the other must contain only data files.
  //Only checking the first file in each playlist because opening each file an extra time
  //remotely (e.g. through xrootd) can get expensive.
  //TODO: Look in INSTALL_DIR if files not found?
  std::string mc_file_list = argv[2],
                    data_file_list = argv[1];
  std::vector<int> petals = {};
  std::cout<<"Here2\n";
  if (argc>3)
  {
    for (int i = 3; i<argc; i++)
    {
      int petalToAdd = std::stoi(argv[i]);
      std::cout<<"Adding petal " << petalToAdd << std::endl;
      petals.push_back(petalToAdd);
    }  
  }
  else //If no argument is given, do all petals
  {
      std::cout<<"No petal specified. Doing all" << std::endl;
      petals={-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
  }


  //Check that necessary TTrees exist in the first file of mc_file_list and data_file_list
  std::string reco_tree_name;
  if(!inferRecoTreeNameAndCheckTreeNames(mc_file_list, data_file_list, reco_tree_name))
  {
    std::cerr << "Failed to find required trees in MC playlist " << mc_file_list << " and/or data playlist " << data_file_list << ".\n" << USAGE << "\n";
    return badInputFile;
  }
  char* numSubruns = getenv("NumGridSubruns");
  char* numProcess = getenv("PROCESS");

  int nDataFiles = 0;
  int nMCFiles = 0;
  if (numSubruns != nullptr && numProcess != nullptr)
  {
    nSubruns = std::stoi(numSubruns);
    nProcess = std::stoi(numProcess);
    std::cout<<"Processing files in " << numSubruns << " job runs. Currently working on: " << numProcess << "/" << numSubruns << std::endl;
    mc_file_list= breakUpInputFileList(mc_file_list, nSubruns , nProcess, nMCFiles);
    data_file_list= breakUpInputFileList(data_file_list, nSubruns , nProcess, nDataFiles);
  }

  const bool doCCQENuValidation = (reco_tree_name == "CCQENu"); //Enables extra histograms and might influence which systematics I use.
  std::cout<<"Test0\n";
  //const bool is_grid = false; //TODO: Are we going to put this back?  Gonzalo needs it iirc.
  PlotUtils::MacroUtil options(reco_tree_name, mc_file_list, data_file_list, "minervame1A", true); //minervame1A is just a placeholder, it gets overwritted immediately below
  options.m_plist_string = util::GetPlaylist(*options.m_mc, true); //TODO: Put GetPlaylist into PlotUtils::MacroUtil
  std::cout<<"Test1\n";
  // You're required to make some decisions
  PlotUtils::MinervaUniverse::SetNuEConstraint(true);
  PlotUtils::MinervaUniverse::SetPlaylist(options.m_plist_string); //TODO: Infer this from the files somehow?
  int nuoranu = util::nuOrAntiNuMode(options.m_plist_string);
  int nupdg;
  if (nuoranu==1) nupdg = 14;
  else if (nuoranu==2) nupdg = -14;
  PlotUtils::MinervaUniverse::SetAnalysisNuPDG(nupdg);
  PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
  PlotUtils::MinervaUniverse::SetZExpansionFaReweight(false);

  PlotUtils::MinervaUniverse::RPAMaterials(true); 
  std::cout<<"Test2\n";

  int filledOrEmpty = util::filledOrEmptyMEPlaylist(options.m_plist_string); //filled = 1, empty = 2
  SMSCalibrator calibrator = util::loadCalibrator((nupdg == 14), (filledOrEmpty==1));
  // ============================ End - Setting up analysis  ============================
  
  // ============================ Tune  ============================

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
  MnvTune.emplace_back(new PlotUtils::AsitsRobFineFluxReweighter<CVUniverse, MichelEvent>());
  MnvTune.emplace_back(new PlotUtils::FluxAndCVReweighter<CVUniverse, MichelEvent>());
  //Setting the A component of mnvtune vA.B.C
  if (tuneA == 1)
  {
    MnvTune.emplace_back(new PlotUtils::GENIEReweighter<CVUniverse, MichelEvent>(true, false));
    MnvTune.emplace_back(new PlotUtils::MINOSEfficiencyReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::RPAReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::LowRecoil2p2hReweighter<CVUniverse, MichelEvent>());
  }
  if (tuneA == 2)
  {
    MnvTune.emplace_back(new PlotUtils::GENIEReweighter<CVUniverse, MichelEvent>(true, false));
    MnvTune.emplace_back(new PlotUtils::MINOSEfficiencyReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::RPAReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::LowRecoil2p2hReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::LowQ2PiReweighter<CVUniverse, MichelEvent>("JOINT")); //Is JOINT the correct option for mnvtune2?
  }
  if (tuneA == 3)
  {
    MnvTune.emplace_back(new PlotUtils::GENIEReweighter<CVUniverse, MichelEvent>(true, false));
    MnvTune.emplace_back(new PlotUtils::MINOSEfficiencyReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::RPAReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::SuSAFromValencia2p2hReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::BodekRitchieReweighter<CVUniverse, MichelEvent>(2)); //Is 2 the right mode?
  }
  if (tuneA == 4)
  {
    PlotUtils::MinervaUniverse::SetDeuteriumGeniePiTune(true);
    PlotUtils::MinervaUniverse::SetReadoutVolume("Tracker");
    PlotUtils::MinervaUniverse::SetMHRWeightNeutronCVReweight( true );
    PlotUtils::MinervaUniverse::SetMHRWeightElastics( true );
    PlotUtils::MinervaUniverse::SetNonResPiReweight(true);
    MnvTune.emplace_back(new PlotUtils::GENIEReweighter<CVUniverse, MichelEvent>(true, true));
    MnvTune.emplace_back(new PlotUtils::MINOSEfficiencyReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::RPAReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::LowRecoil2p2hReweighter<CVUniverse, MichelEvent>());
  }
  //Setting the B component of mnvtune vA.B.C
  if (tuneB == 3)
  {
    MnvTune.emplace_back(new PlotUtils::LowQ2PiReweighter<CVUniverse, MichelEvent>("MENU1PI"));
    MnvTune.emplace_back(new PlotUtils::DiffractiveReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::COHPionReweighter<CVUniverse, MichelEvent>());
  }
  if (tuneB == 6)
  {
    MnvTune.emplace_back(new PlotUtils::LowQ2PiReweighter<CVUniverse, MichelEvent>("MENU1PI"));
    MnvTune.emplace_back(new PlotUtils::DiffractiveReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::COHPionReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::PionReweighter<CVUniverse, MichelEvent>());
  }
  if (tuneB == 7)
  {
    MnvTune.emplace_back(new PlotUtils::LowQ2PiReweighter<CVUniverse, MichelEvent>("MENU1PI"));
    MnvTune.emplace_back(new PlotUtils::DiffractiveReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::COHPionReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::PionReweighter<CVUniverse, MichelEvent>());
    MnvTune.emplace_back(new PlotUtils::Q2TpiReweighter<CVUniverse, MichelEvent>());
  }
  //Setting the C component of mnvtune vA.B.C
  if (tuneC == 1)
  {
    MnvTune.emplace_back(new PlotUtils::FSIReweighter<CVUniverse, MichelEvent>(true, false));
  }

  // ============================ End - Tune  ============================

  // ============================  Warps  ============================
  const bool NO_2P2H_WARP = (getenv("NO_2P2H_WARP") != nullptr);
  if (NO_2P2H_WARP)
  {
    std::cout << "Turning off LowRecoil2p2hReweighter because environment variable NO_2P2H_WARP is set.\n";
  }
  const bool AMU_DIS_WARP = (getenv("AMU_DIS_WARP") != nullptr);
  if (AMU_DIS_WARP)
  {
    std::cout << "Turning on AMUDISReweighter because environment variable AMU_DIS_WARP is set.\n";
  }
  const bool LOW_Q2_PION_WARP = (getenv("LOW_Q2_PION_WARP") != nullptr);
  if (LOW_Q2_PION_WARP)
  {
    std::cout << "Turning on LowQ2PiReweighter because environment variable LOW_Q2_PION_WARP is set.\n";
  }
  const bool SUSA_2P2H_WARP = (getenv("SUSA_2P2H_WARP") != nullptr);
  if (SUSA_2P2H_WARP)
  {
    std::cout << "Replace LowRecoil2p2hReweighter with SuSAFromValencia2p2hReweighter because environment variable SUSA_2P2H_WARP is set.\n";
  }
  const bool NO_LOW_Q2Pi = (getenv("NO_LOW_Q2Pi") != nullptr);
  if (NO_LOW_Q2Pi)
  {
    std::cout << "Turning off LowQ2PiReweighter because environment variable NO_LOW_Q2Pi is set.\n";
  }
  const bool MARESNOM = (getenv("MA_RES_NOMINAL") != nullptr);
  if (MARESNOM)
  {
    std::cout << "Setting M_A^RES AND CCResNorm to their GENIE nominal values 1.12 and 1.0, respectively.\n";
  }
  if (NO_2P2H_WARP) //Removing LowRecoil2p2hReweighter
  {
    auto it = find_if(MnvTune.begin(), MnvTune.end(), [] (auto& w) { return w->GetName() == "LowRecoil2p2hTune"; } );
    if (it!=MnvTune.end())
    {
      MnvTune.erase(it);
      std::cout<<"Applying NO_2P2H_WARP - removing found LowRecoil2p2hReweighter\n";
    }
    else std::cout<<"Warning - NO_2P2H_WARP - Could not apply warp since there were no 2p2h reweighters found\n";
  }
  if (AMU_DIS_WARP)
  {
    std::cout<<"Applying AMU_DIS_WARP - Applying AMUDISReweighter\n";
    MnvTune.emplace_back(new PlotUtils::AMUDISReweighter<CVUniverse, MichelEvent>());
  }
  if (LOW_Q2_PION_WARP)  // Low Q2 pion suppression (mnvtunev2)
  {
    std::cout<<"Applying LOW_Q2_PION_WARP - Applying LowQ2PiReweighter\n";
    MnvTune.emplace_back(new PlotUtils::LowQ2PiReweighter<CVUniverse, MichelEvent>("JOINT"));
  }
  if (SUSA_2P2H_WARP) //Replacing LowRecoil2p2hReweighter with SuSAFromValencia2p2hReweighter
  {
    auto it = find_if(MnvTune.begin(), MnvTune.end(), [] (auto& w) { return w->GetName() == "LowRecoil2p2hTune"; } );
      if (it!=MnvTune.end())
      {
        MnvTune.erase(it);
        std::cout << "Removing existing LowRecoil2p2hTune.\n";
      }
    MnvTune.emplace_back(new PlotUtils::SuSAFromValencia2p2hReweighter<CVUniverse, MichelEvent>());
  }
  if (NO_LOW_Q2Pi) //Removing LowRecoil2p2hReweighter
  {
    auto it = find_if(MnvTune.begin(), MnvTune.end(), [] (auto& w) { return w->GetName() == "LowQ2Pi_MENU1PI"; } );
    if (it!=MnvTune.end())
    {
      MnvTune.erase(it);
      std::cout<<"Applying NO_LOW_Q2Pi - removing found LowQ2Pi_MENU1PI\n";
    }
    else std::cout<<"Warning - NO_LOW_Q2Pi - Could not apply warp since there were no LowQ2Pi_MENU1PI reweighters found\n";
  }
  if (MARESNOM)
  {
    auto it = find_if(MnvTune.begin(), MnvTune.end(), [] (auto& w) { return w->GetName() == "GENIE"; } );
      if (it!=MnvTune.end())
      {
        MnvTune.erase(it);
        std::cout << "Removing model baseline GENIE tune.\n";
      }
    MnvTune.emplace_back(new PlotUtils::GENIEReweighter<CVUniverse, MichelEvent>(true, false));
  }
  const bool NO_CALIBRATION = (getenv("NO_CALIBRATION") != nullptr);
  if (NO_CALIBRATION)
  {
    std::cout << "Turning off SMS ANN Calibration.\n";
    useCalibrator = false;
  }
  std::cout<<"Tune components applied:\n";
  for (auto&& t : MnvTune) std::cout<< "\t"<< t->GetName() <<std::endl;
  
  // ============================ End - Warps  =============================

  // ============================  Error Bands  ============================

  // Make a map of systematic universes
  // Leave out systematics when making validation histograms
  const bool doSystematics = (getenv("MNV101_SKIP_SYST") == nullptr);
  if (!doSystematics)
  {
    std::cout << "Skipping systematics (except 1 flux universe) because environment variable MNV101_SKIP_SYST is set.\n";
    PlotUtils::MinervaUniverse::SetNFluxUniverses(2); // Necessary to get Flux integral later...  Doesn't work with just 1 flux universe though because _that_ triggers "spread errors".
  }

  std::map<std::string, std::vector<CVUniverse *>> error_bands;
  if (doSystematics) 
  {
    error_bands = GetStandardSystematics(options.m_mc);
    std::map<std::string, std::vector<CVUniverse *>> modelSyst = GetModelSystematics(options.m_mc, MnvTune);
    error_bands.insert(modelSyst.begin(), modelSyst.end());
  }
  else
  {
    std::map<std::string, std::vector<CVUniverse *>> band_flux = PlotUtils::GetFluxSystematicsMap<CVUniverse>(options.m_mc, CVUniverse::GetNFluxUniverses());
    error_bands.insert(band_flux.begin(), band_flux.end()); // Necessary to get flux integral later...
  }
  error_bands["cv"] = {new CVUniverse(options.m_mc)};

  std::map<std::string, std::vector<CVUniverse *>> truth_bands;
  if (doSystematics)
  {
    truth_bands = GetStandardSystematics(options.m_truth);
    std::map<std::string, std::vector<CVUniverse *>> modelSyst = GetModelSystematics(options.m_truth, MnvTune);
    truth_bands.insert(modelSyst.begin(), modelSyst.end());
  }
  else{
    std::map<std::string, std::vector<CVUniverse*> > band_flux = PlotUtils::GetFluxSystematicsMap<CVUniverse>(options.m_truth, CVUniverse::GetNFluxUniverses());
    truth_bands.insert(band_flux.begin(), band_flux.end()); //Necessary to get cross-section later.
  }
  truth_bands["cv"] = {new CVUniverse(options.m_truth)};

  CVUniverse *data_universe = new CVUniverse(options.m_data);
  std::vector<CVUniverse *> data_band = {data_universe};
  std::map<std::string, std::vector<CVUniverse *>> data_error_bands;
  data_error_bands["cv"] = data_band;

  PlotUtils::Model<CVUniverse, MichelEvent> model(std::move(MnvTune));
  // ============================  End - Error Bands  ============================

  std::vector<std::vector<std::unique_ptr<Variable>>> varsVector;
  std::vector<std::vector<std::unique_ptr<Variable2D>>> vars2DVector;
  std::vector<PlotUtils::Cutter<CVUniverse, MichelEvent>*> cutsVector;
  std::vector<std::vector<Study*>> studiesVector;
  std::vector<std::vector<Study*>> dataStudiesVector;


  std::vector<double> RecoilBins, segmentBins, angleBins;
  //const double RecoilBinWidth = 50; //MeV
  //for(int whichBin = 0; whichBin < 100 + 1; ++whichBin) RecoilBins.push_back(RecoilBinWidth * whichBin);

  const double numsegments = 180;
  for (double whichBin = 0; whichBin < numsegments; whichBin++)
    segmentBins.push_back(whichBin - 0.5);

  const double nAngleBins = 170;
  for (double whichBin = 0; whichBin < nAngleBins+1; whichBin++)
    angleBins.push_back(0.1*whichBin);

  //std::function<double(const CVUniverse&)> ANNRecoilEGeV = [](const CVUniverse& univ) { return univ.GetANNRecoilE()/1000;};
  std::function<double(const CVUniverse&)> q0TrueGeV = [](const CVUniverse& univ) { return univ.Getq0True()/1000;};
  std::function<double(const CVUniverse &)> muonAngleDegrees = [](const CVUniverse &univ) { return (univ.GetThetamu() * 180 / M_PI); };
  std::function<double(const CVUniverse &)> muonAngleDegreesTruth = [](const CVUniverse &univ) { return (univ.GetDouble("truth_muon_theta")* 180 / M_PI); };
  std::function<double(const CVUniverse &)> ANNConfPercentage = [](const CVUniverse &univ) { return univ.GetANNProb() * 100; };
  std::function<double(const CVUniverse &)> ANNConfPercentageUncalib = [](const CVUniverse &univ) { return univ.GetANNProb(false) * 100; };
  std::function<double(const CVUniverse &)> nullFunc = [](const CVUniverse &univ) { return 0; };  
  std::function<double(const CVUniverse &)> GetANNSegment = [](const CVUniverse &univ) { return univ.GetANNSegment(); };
  std::function<double(const CVUniverse &)> GetANNProb = [](const CVUniverse &univ) { return univ.GetANNProb(); };

  for (int i = 0; i < petals.size(); i++)
  {
    int ptl = petals[i];
    std::cout << "Trying petal: " << ptl << std::endl;

    //Now that we've defined what a cross section is, decide which sample and model
    //we're extracting a cross section for.
    PlotUtils::Cutter<CVUniverse, MichelEvent>::reco_t sidebands, preCuts;
    PlotUtils::Cutter<CVUniverse, MichelEvent>::truth_t signalDefinition, phaseSpace;

    preCuts = util::GetAnalysisCuts(nupdg);
    //preCuts.emplace_back(new reco::ZRangeANN<CVUniverse, MichelEvent>("Z pos in active tracker", PlotUtils::TargetProp::Tracker::Face, PlotUtils::TargetProp::Tracker::Back));
    preCuts.emplace_back(new reco::ZRangeANN<CVUniverse, MichelEvent>("Z pos in active tracker", 5980, 8422));  //Using Anezka's tracker Z definition rather than the the TargetUtils one
    preCuts.emplace_back(new reco::IsInPetal<CVUniverse, MichelEvent>(ptl));

                                                                                                                                                      
    if (nupdg>0) signalDefinition.emplace_back(new truth::IsNeutrino<CVUniverse>());
    else if (nupdg<0) signalDefinition.emplace_back(new truth::IsAntiNeutrino<CVUniverse>());
    signalDefinition.emplace_back(new truth::IsCC<CVUniverse>());
                                                                                                                                                      
    phaseSpace = util::GetPhaseSpace();
    phaseSpace.emplace_back(new truth::ZRange<CVUniverse>("Z pos in active tracker", PlotUtils::TargetProp::Tracker::Face, PlotUtils::TargetProp::Tracker::Back));
    phaseSpace.emplace_back(new truth::IsInPetal<CVUniverse>(ptl));

    cutsVector.push_back(new PlotUtils::Cutter<CVUniverse, MichelEvent>(std::move(preCuts), std::move(sidebands), std::move(signalDefinition), std::move(phaseSpace)));

    // ----------------- Variables ----------------- 
    std::vector<std::unique_ptr<Variable>> vars;
    std::vector<std::unique_ptr<Variable2D>> vars2D;
    vars.push_back(std::make_unique<Variable>("pTmu", "p_{T, #mu} [GeV/c]", util::PTBins, &CVUniverse::GetANNMuonPTGeV, &CVUniverse::GetMuonPTTrue));
    vars.push_back(std::make_unique<Variable>("pZmu", "p_{||, #mu} [GeV/c]", util::PzBins, &CVUniverse::GetANNMuonPzGeV, &CVUniverse::GetMuonPzTrue));
    vars.push_back(std::make_unique<Variable>("Emu", "E_{#mu} [GeV]", util::EmuBins, &CVUniverse::GetANNEmuGeV, &CVUniverse::GetElepTrueGeV));
    vars.push_back(std::make_unique<Variable>("Enu", "E_{#nu} [GeV]", util::EnuBins, &CVUniverse::GetANNEnuGeV, &CVUniverse::GetEnuTrueGeV));
    vars.push_back(std::make_unique<Variable>("Erecoil", "E_{recoil} [GeV]", util::Erecoilbins, &CVUniverse::GetANNRecoilEGeV, q0TrueGeV)); // TODO: q0 is not the same as recoil energy without a spline correction
    vars.push_back(std::make_unique<Variable>("BjorkenX", "X", util::bjorkenXbins, &CVUniverse::GetBjorkenX, &CVUniverse::GetCalcBjorkenXTrue));
    //vars.push_back(std::make_unique<Variable>("BjorkenXTest", "X_Test", util::bjorkenXbinsTest, &CVUniverse::GetBjorkenX, &CVUniverse::GetCalcBjorkenXTrue));
    //vars.push_back(std::make_unique<Variable>("BjorkenY", "Y", util::bjorkenYbins, &CVUniverse::GetBjorkenY, &CVUniverse::GetCalcBjorkenYTrue));
    vars.push_back(std::make_unique<Variable>("BjorkenX_2", "X_2", util::bjorkenXbins2, &CVUniverse::GetBjorkenX, &CVUniverse::GetCalcBjorkenXTrue));
    vars.push_back(std::make_unique<Variable>("BjorkenX_4", "X_4", util::bjorkenXbins4, &CVUniverse::GetBjorkenX, &CVUniverse::GetCalcBjorkenXTrue));
    vars.push_back(std::make_unique<Variable>("segment", "segmentNum", segmentBins, GetANNSegment, &CVUniverse::GetTruthSegment)); // Just used for plotting events by detector position tbh - not for any actual physics
    //vars.push_back(std::make_unique<Variable>("beamAngle", "Angle", angleBins, muonAngleDegrees, muonAngleDegreesTruth));              // Neutrino angle 
    vars2D.push_back(std::make_unique<Variable2D>("pTmu_pZmu", *vars[0], *vars[1]));
    //vars2D.push_back(std::make_unique<Variable2D>("Emu_Erecoil", *vars[2], *vars[3]));
    //vars2D.push_back(new Variable2D("BjorkenX_BjorkenY", *vars[4], *vars[5]));
    //Probe some more variables?????

    //Testing ANN probability distribution
    //std::vector<double> probBins;
    //const double nProbBins = 20;
    //for (double whichBin = 0; whichBin < nProbBins+1; whichBin++)
    //{
    //  probBins.push_back((100/nProbBins)*whichBin);
    //}
    //vars.push_back(std::make_unique<Variable>("ANNConfidence", "ANN Confidence [%]", probBins, ANNConfPercentage, nullFunc));
    //vars.push_back(std::make_unique<Variable>("ANNConfidenceUncalib", "ANN Confidence (Uncalibrated) [%]", probBins, ANNConfPercentageUncalib, nullFunc));
    


    for(auto& var: vars) var->InitializeMCHists(error_bands, truth_bands);
    for(auto& var: vars) var->InitializeDATAHists(data_band);

    for(auto& var: vars2D) var->InitializeMCHists(error_bands, truth_bands);
    for(auto& var: vars2D) var->InitializeDATAHists(data_band);
    //std::vector<Study*> studies;
    //std::vector<Study*> data_studies;
    varsVector.push_back(std::move(vars));
    vars2DVector.push_back(std::move(vars2D));
    studiesVector.push_back(std::vector<Study*>()); //No studies right now 
    dataStudiesVector.push_back(std::vector<Study*>());
  }

  // Loop entries and fill
  //try
  {
    CVUniverse::SetTruth(false);
    std::cout<<"Daisy reweight test0\n";
    LoopAndFillEventSelection(options.m_mc, error_bands, varsVector, vars2DVector, studiesVector, cutsVector, model, petals, calibrator);
    CVUniverse::SetTruth(true);
    std::cout<<"Daisy reweight test1\n";
    LoopAndFillEffDenom(options.m_truth, truth_bands, varsVector, vars2DVector, cutsVector, model, petals);  //No calibrator needed with Eff Denom since it only uses truth values
      std::cout<<"Daisy reweight test2\n";
    options.PrintMacroConfiguration(argv[0]);
    for (int i = 0; i < petals.size(); i++)
    {
      std::cout << "Tracker Daisy MC cut summary for petal: "<< petals[i] <<"\n"
                << *cutsVector[i] << "\n";
      cutsVector[i]->resetStats();
    }

    CVUniverse::SetTruth(false);
    LoopAndFillData(options.m_data, data_band, varsVector, vars2DVector, dataStudiesVector, cutsVector, petals, calibrator);
    for (int i = 0; i < petals.size(); i++)
    {
      std::cout << "Tracker Daisy Data cut summary for petal: "<< petals[i] <<"\n"
                << *cutsVector[i] << "\n";
    }

    // ============================  Saving Out  ============================

    auto playlistStr = new TNamed("PlaylistUsed", options.m_plist_string);

    for (int i = 0; i < petals.size(); i++)
    {
      int ptl = petals[i];
      std::string mcOutFileName = std::string(MC_OUT_FILE_NAME_BASE) +"_petal_"+ std::to_string(ptl)+".root";
      if (nSubruns != 0 ) mcOutFileName = std::string(MC_OUT_FILE_NAME_BASE) +"_petal_"+ std::to_string(ptl) + "_n"+nProcess+ ".root";
      //if (ptl==-1) mcOutFileName = std::string(MC_OUT_FILE_NAME_BASE) + ".root";
      //else mcOutFileName = std::string(MC_OUT_FILE_NAME_BASE) +"_petal_"+ std::to_string(ptl)+".root";
      //Write MC results
      TFile* mcOutDir = TFile::Open(mcOutFileName.c_str(), "RECREATE");
      if(!mcOutDir)
      {
        std::cerr << "Failed to open a file named " << mcOutFileName << " in the current directory for writing histograms.\n";
        return badOutputFile;
      }

      for(auto& study: studiesVector[i]) study->SaveOrDraw(*mcOutDir);
      for(auto& var: varsVector[i]) var->WriteMC(*mcOutDir);
      for(auto& var: vars2DVector[i]) var->WriteMC(*mcOutDir);

      //Playlist name - Used for flux calculations later on
      playlistStr->Write();
      //Protons On Target
      double mcpot = options.m_mc_pot;
      auto mcPOT = new TParameter<double>("POTUsed", mcpot);
      mcPOT->Write();

      PlotUtils::TargetUtils targetInfo;
      assert(error_bands["cv"].size() == 1 && "List of error bands must contain a universe named \"cv\" for the flux integral.");

      //I suggest getting the flux integral and the fiducial nucleons in the ExtractCrossSection step and leaving this file purely for event loop things
      //That also makes debugging a lot easier in the case of issues with your flux or nnucleons
      /* for(const auto& var: varsVector[i])
      {
        //Flux integral only if systematics are being done (temporary solution)
        util::GetFluxIntegral(*error_bands["cv"].front(), var->efficiencyNumerator->hist)->Write((var->GetName() + "_reweightedflux_integrated").c_str());
        //Always use MC number of nucleons for cross section
        auto nNucleons = new TParameter<double>((var->GetName() + "_fiducial_nucleons").c_str(), targetInfo.GetTrackerNNucleons(PlotUtils::TargetProp::Tracker::Face, PlotUtils::TargetProp::Tracker::Back, true, util::apothem));
        nNucleons->Write();
      }
      for(const auto& var: vars2DVector[i])
      {
        //Flux integral only if systematics are being done (temporary solution)
        util::GetFluxIntegral(*error_bands["cv"].front(), var->efficiencyNumerator->hist)->Write((var->GetName() + "_reweightedflux_integrated").c_str());
        //Always use MC number of nucleons for cross section
        auto nNucleons = new TParameter<double>((var->GetName() + "_fiducial_nucleons").c_str(), targetInfo.GetTrackerNNucleons(PlotUtils::TargetProp::Tracker::Face, PlotUtils::TargetProp::Tracker::Back, true, util::apothem));
        nNucleons->Write();
      } */
      mcOutDir->Close();

      playlistStr = new TNamed("PlaylistUsed", options.m_plist_string);

      //Write data results

      std::string dataOutFileName = std::string(DATA_OUT_FILE_NAME_BASE) +"_petal_"+ std::to_string(ptl)+".root";
      if (nSubruns != 0 ) dataOutFileName = std::string(DATA_OUT_FILE_NAME_BASE) +"_petal_"+ std::to_string(ptl) + "_n"+nProcess+ ".root";

      //if (ptl==-1) dataOutFileName = std::string(DATA_OUT_FILE_NAME_BASE) + ".root";
      //else dataOutFileName = std::string(DATA_OUT_FILE_NAME_BASE) +"_petal_"+ std::to_string(ptl)+".root";
      TFile* dataOutDir = TFile::Open(dataOutFileName.c_str(), "RECREATE");
      if(!dataOutDir)
      {
        std::cerr << "Failed to open a file named " << dataOutFileName << " in the current directory for writing histograms.\n";
        return badOutputFile;
      }

      for(auto& var: varsVector[i]) var->WriteData(*dataOutDir);
      for(auto& var: vars2DVector[i]) var->WriteData(*dataOutDir);

      playlistStr->Write();
      //Protons On Target
      double datapot = options.m_data_pot;
      auto dataPOT = new TParameter<double>("POTUsed", datapot);
      dataPOT->Write();

      dataOutDir->Close();

      //Saving 2D migration matrices
      //Putting this right at the end in case of a crash
      std::string migrationOutDirName = std::string(MIGRATION_2D_OUT_FILE_NAME_BASE) +"_petal_"+ std::to_string(ptl)+".root";
      if (nSubruns != 0 ) migrationOutDirName = std::string(MIGRATION_2D_OUT_FILE_NAME_BASE) +"_petal_"+ std::to_string(ptl) + "_n"+nProcess+ ".root";

      //if (ptl==-1) migrationOutDirName = std::string(MIGRATION_2D_OUT_FILE_NAME_BASE) + ".root";
      //else migrationOutDirName = std::string(MIGRATION_2D_OUT_FILE_NAME_BASE) +"_petal_"+ std::to_string(ptl)+".root";
      TFile* migrationOutDir = TFile::Open(migrationOutDirName.c_str(), "RECREATE");
      if(!migrationOutDir)
      {
        std::cerr << "Failed to open a file named " << migrationOutDirName << " in the current directory for writing histograms.\n";
        return badOutputFile;
      }
      for(auto& var: vars2DVector[i])var->WriteMigration(*migrationOutDir);

      migrationOutDir->Close();

      options.m_data->tree->PrintCacheStats();
      options.m_mc->tree->PrintCacheStats();
      options.m_truth->tree->PrintCacheStats();
      
      std::cout << "Success" << std::endl;
    }
    // ============================  End Saving Out  ============================
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