#define MC_OUT_FILE_NAME "runEventLoopTrackerMC"
#define DATA_OUT_FILE_NAME "runEventLoopTrackerData"
#define MIGRATION_OUT_FILE_NAME "runEventLoopTracker2DMigration"

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

#include "util/NukeUtils.h"
#pragma GCC diagnostic pop

//ROOT includes
#include "TParameter.h"

#include "Math/Vector3D.h"
#include "TH3D.h"

//c++ includes
#include <iostream>
#include <cstdlib> //getenv()
#include <fstream>
#include <sstream> //reading input files

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
  //const int nEntries = 10000;
  //for (int i=19598; i<nEntries; ++i)
  PlotUtils::TargetUtils util;
  for (int i=0; i<nEntries; ++i)
  {
    if(i%1000==0) std::cout << i << " / " << nEntries << "\r" <<std::flush;
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
         
        // This is where you would Access/create a Michel

        //weight is ignored in isMCSelected() for all but the CV Universe.
        if (!michelcuts.isMCSelected(*universe, myevent, cvWeight).all()) continue; //all is another function that will later help me with sidebands
        const double weight = model.GetWeight(*universe, myevent); //Only calculate the per-universe weight for events that will actually use it.
        //Get daisy petal
        int petal = util.GetDaisyPetal(universe->GetANNVertex().X(), universe->GetANNVertex().Y() ); // *10 to convert from cm to mm
        //if (petal>0) std::cout<<"Petal: " <<petal << std::endl;
        for(auto& var: vars)
        { 
          var->selectedMCReco->FillUniverse(universe, var->GetRecoValue(*universe), weight); //"Fake data" for closure
          (*var->m_intChannels)[universe->GetInteractionType()].FillUniverse(universe, var->GetRecoValue(*universe), weight);
          if (petal >= 0 && petal < 12)
          {
            var->selectedMCRecoDaisy[petal]->FillUniverse(universe, var->GetRecoValue(*universe), weight); //"Fake data" for closure
            (*(var->ChannelsDaisy[petal]))[universe->GetInteractionType()].FillUniverse(universe, var->GetRecoValue(*universe), weight);
          }
        }
        for(auto& var: vars2D)
        {
          var->selectedMCReco->FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight); //"Fake data" for closure
          (*var->m_intChannels)[universe->GetInteractionType()].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
          if (petal >= 0 && petal < 12)
          {
            var->selectedMCRecoDaisy[petal]->FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight); //"Fake data" for closure
            (*(var->ChannelsDaisy[petal]))[universe->GetInteractionType()].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
          }
          //(*var->m_interactionChannels)[universe->GetInteractionType()].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
        }
        const bool isSignal = michelcuts.isSignal(*universe, weight);        if(isSignal)
        {
          for(auto& study: studies) study->SelectedSignal(*universe, myevent, weight);

          for(auto& var: vars)
          {
            //Cross section components
            var->efficiencyNumerator->FillUniverse(universe, var->GetTrueValue(*universe), weight);
            var->migration->FillUniverse(universe, var->GetRecoValue(*universe), var->GetTrueValue(*universe), weight);
            var->selectedSignalReco->FillUniverse(universe, var->GetRecoValue(*universe), weight); //Efficiency numerator in reco variables.  Useful for warping studies.
            //Daisy reweight
            if (petal >= 0 && petal < 12)
            {
              var->EffNumDaisy[petal]->FillUniverse(universe, var->GetTrueValue(*universe), weight);
              var->MigrationDaisy[petal]->FillUniverse(universe, var->GetRecoValue(*universe), var->GetTrueValue(*universe), weight);
              var->selectedSignalRecoDaisy[petal]->FillUniverse(universe, var->GetRecoValue(*universe), weight);
            }
          }
          for(auto& var: vars2D)
          {
            var->efficiencyNumerator->FillUniverse(universe, var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
            var->migration->Fill(var->GetRecoValueX(*universe), var->GetRecoValueY(*universe),var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
            var->selectedSignalReco->FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight); //Efficiency numerator in reco variables.  Useful for warping studies.
            //Daisy reweight
            if (petal >= 0 && petal < 12)
            {
              var->EffNumDaisy[petal]->FillUniverse(universe, var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
              var->MigrationDaisy[petal]->Fill(var->GetRecoValueX(*universe), var->GetRecoValueY(*universe),var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
              var->selectedSignalRecoDaisy[petal]->FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
            }
          }
        }
        else
        {
          int bkgd_ID = -1;
          if (universe->GetCurrent()==2)bkgd_ID=0;
          else bkgd_ID=1;

          for(auto& var: vars)
          {
            (*var->m_backgroundHists)[bkgd_ID].FillUniverse(universe, var->GetRecoValue(*universe), weight);
            if (petal >= 0 && petal < 12) (*(var->BackgroundsDaisy[petal]))[bkgd_ID].FillUniverse(universe, var->GetRecoValue(*universe), weight);
          }
          for(auto& var: vars2D)
          {
            (*var->m_backgroundHists)[bkgd_ID].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
            if (petal >= 0 && petal < 12) (*(var->BackgroundsDaisy[petal]))[bkgd_ID].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
          } 
        }
      } // End band's universe loop
    } // End Band loop
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
  PlotUtils::TargetUtils util;
  //const int nEntries = 10000;
  for (int i=0; i<data->GetEntries(); ++i) {
    for (auto universe : data_band) {
      universe->SetEntry(i);
      if(i%1000==0) std::cout << i << " / " << nEntries << "\r" << std::flush;
      MichelEvent myevent; 
      if (!michelcuts.isDataSelected(*universe, myevent).all()) continue;
      //Get daisy petal
      int petal = util.GetDaisyPetal(universe->GetANNVertex().X(), universe->GetANNVertex().Y() );
      for(auto& study: studies) study->Selected(*universe, myevent, 1); 
      //std::cout<<"petal: " <<petal << "X: " << universe->GetANNVertex().X() << "\tY: "<<universe->GetANNVertex().Y()<<std::endl;
      for(auto& var: vars)
      {
        var->dataHist->FillUniverse(universe, var->GetRecoValue(*universe, myevent.m_idx), 1);
        if (petal >= 0 && petal < 12)  var->dataDaisy[petal]->FillUniverse(*universe, var->GetRecoValue(*universe), 1); //"Fake data" for closure
      }

      for(auto& var: vars2D)
      {
        var->dataHist->FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), 1);
        if (petal >= 0 && petal < 12)  var->dataDaisy[petal]->FillUniverse(*universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), 1); //"Fake data" for closure
      }
    }
  }
  std::cout << "Finished data loop.\n";
}

void LoopAndFillEffDenom( PlotUtils::ChainWrapper* truth,
    				std::map<std::string, std::vector<CVUniverse*> > truth_bands,
    				std::vector<Variable*> vars,
                                std::vector<Variable2D*> vars2D,
    				PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts,
                                PlotUtils::Model<CVUniverse, MichelEvent>& model)
{
  assert(!truth_bands["cv"].empty() && "\"cv\" error band is empty!  Could not set Model entry.");
  auto& cvUniv = truth_bands["cv"].front();
  PlotUtils::TargetUtils util;
  std::cout << "Starting efficiency denominator loop...\n";
  const int nEntries = truth->GetEntries();
  //const int nEntries = 10000;
  for (int i=0; i<nEntries; ++i)
  {
    if(i%1000==0) std::cout << i << " / " << nEntries << "\r" << std::flush;

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

        if (!michelcuts.isEfficiencyDenom(*universe, cvWeight)) continue; //Weight is ignored for isEfficiencyDenom() in all but the CV universe 
        const double weight = model.GetWeight(*universe, myevent); //Only calculate the weight for events that will use it
        
        int petal = util.GetDaisyPetal(universe->GetTrueVertex().X(), universe->GetTrueVertex().Y() ); // *10 to convert from cm to mm
        //Fill efficiency denominator now: 
        for(auto var: vars)
        {
          var->efficiencyDenominator->FillUniverse(universe, var->GetTrueValue(*universe), weight);
          (*var->m_intChannelsEffDenom)[universe->GetInteractionType()].FillUniverse(universe, var->GetTrueValue(*universe), weight);
          //Daisy reweight
          if (petal >= 0 && petal < 12)
          {
            var->EffDenomDaisy[petal]->FillUniverse(universe, var->GetTrueValue(*universe), weight);
            (*(var->EffDenomDaisyIntChannels[petal]))[universe->GetInteractionType()].FillUniverse(universe, var->GetTrueValue(*universe), weight);
          }
        }

        for(auto var: vars2D)
        {
          var->efficiencyDenominator->FillUniverse(universe, var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
          (*var->m_intChannelsEffDenom)[universe->GetInteractionType()].FillUniverse(universe, var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
          if (petal >= 0 && petal < 12)
          {
            var->EffDenomDaisy[petal]->FillUniverse(universe, var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
            (*(var->EffDenomDaisyIntChannels[petal]))[universe->GetInteractionType()].FillUniverse(universe, var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
          }
        }
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
  //Now that we've defined what a cross section is, decide which sample and model
  //we're extracting a cross section for.
  PlotUtils::Cutter<CVUniverse, MichelEvent>::reco_t sidebands, preCuts;
  PlotUtils::Cutter<CVUniverse, MichelEvent>::truth_t signalDefinition, phaseSpace;

  //const double minZ = 5980, maxZ = 8422, apothem = 850; //All in mm
  const double apothem = 850; //All in mm
  preCuts.emplace_back(new reco::ZRangeANN<CVUniverse, MichelEvent>("Active Tracker Z pos", PlotUtils::TargetProp::Tracker::Face, PlotUtils::TargetProp::Tracker::Back));
  preCuts.emplace_back(new reco::Apothem<CVUniverse, MichelEvent>(apothem));
  preCuts.emplace_back(new reco::MaxMuonAngle<CVUniverse, MichelEvent>(17.));
  preCuts.emplace_back(new reco::HasMINOSMatch<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::NoDeadtime<CVUniverse, MichelEvent>(1, "Deadtime"));
  if (nupdg>0)  preCuts.emplace_back(new reco::IsNeutrino<CVUniverse, MichelEvent>()); //Used minos curvature
  else if (nupdg<0)  preCuts.emplace_back(new reco::IsAntiNeutrino<CVUniverse, MichelEvent>()); //Used minos curvature
  preCuts.emplace_back(new reco::MuonCurveSignificance<CVUniverse, MichelEvent>(5));
  preCuts.emplace_back(new reco::MuonEnergyMin<CVUniverse, MichelEvent>(2000.0, "EMu Min"));
  preCuts.emplace_back(new reco::MuonEnergyMax<CVUniverse, MichelEvent>(50000.0, "EMu Max"));
  preCuts.emplace_back(new reco::ANNConfidenceCut<CVUniverse, MichelEvent>(0.20));


                                                                                                                                                   
  if (nupdg>0) signalDefinition.emplace_back(new truth::IsNeutrino<CVUniverse>());
  else if (nupdg<0) signalDefinition.emplace_back(new truth::IsAntiNeutrino<CVUniverse>());
  signalDefinition.emplace_back(new truth::IsCC<CVUniverse>());
                                                                                                                                                   
  phaseSpace.emplace_back(new truth::ZRange<CVUniverse>("Active Tracker Z pos", PlotUtils::TargetProp::Tracker::Face, PlotUtils::TargetProp::Tracker::Back));
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
  vars.push_back(new Variable("tracker_pTmu", "p_{T, #mu} [GeV/c]", dansPTBins, &CVUniverse::GetMuonPT, &CVUniverse::GetMuonPTTrue));
  vars.push_back(new Variable("tracker_pZmu", "p_{||, #mu} [GeV/c]", dansPzBins, &CVUniverse::GetMuonPz, &CVUniverse::GetMuonPzTrue));
  vars.push_back(new Variable("tracker_Emu", "E_{#mu} [GeV]", robsEmuBins, &CVUniverse::GetEmuGeV, &CVUniverse::GetElepTrueGeV));
  vars.push_back(new Variable("tracker_Erecoil", "E_{recoil}", robsRecoilBins, &CVUniverse::GetRecoilE, &CVUniverse::Getq0True)); //TODO: q0 is not the same as recoil energy without a spline correction
  vars.push_back(new Variable("tracker_BjorkenX", "X", bjorkenXbins, &CVUniverse::GetBjorkenX, &CVUniverse::GetBjorkenXTrue));
  vars2D.push_back(new Variable2D("tracker_pTmu_pZmu", *vars[1], *vars[0]));

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
    std::cout<<"Daisy reweight test0\n";
    LoopAndFillEventSelection(options.m_mc, error_bands, vars, vars2D, studies, mycuts, model);
    CVUniverse::SetTruth(true);
    std::cout<<"Daisy reweight test1\n";
    LoopAndFillEffDenom(options.m_truth, truth_bands, vars, vars2D, mycuts, model);
      std::cout<<"Daisy reweight test2\n";
    options.PrintMacroConfiguration(argv[0]);
    std::cout << "MC cut summary:\n" << mycuts << "\n";
    mycuts.resetStats();

    CVUniverse::SetTruth(false);
    LoopAndFillData(options.m_data, data_band, vars, vars2D, data_studies, mycuts);
    std::cout << "Data cut summary:\n" << mycuts << "\n";
    std::string filepath;

    //Write MC results
    filepath = std::string(MC_OUT_FILE_NAME)+".root";
    TFile* mcOutDir = TFile::Open(filepath.c_str(), "RECREATE");
    if(!mcOutDir)
    {
      std::cerr << "Failed to open a file named " << filepath << " in the current directory for writing histograms.\n";
      return badOutputFile;
    }

    for(auto& study: studies) study->SaveOrDraw(*mcOutDir);
    for(auto& var: vars) var->WriteMC(*mcOutDir);
    for(auto& var: vars2D) var->WriteMC(*mcOutDir);

    //Protons On Target
    auto mcPOT = new TParameter<double>("POTUsed", options.m_mc_pot);
    mcPOT->Write();

    //Write MC Daisy results
    filepath = std::string(MC_OUT_FILE_NAME)+"Daisy.root";
    TFile* mcDaisyOutDir = TFile::Open(filepath.c_str(), "RECREATE");
    if(!mcDaisyOutDir)
    {
      std::cerr << "Failed to open a file named " << filepath << " in the current directory for writing histograms.\n";
      return badOutputFile;
    }

    for(auto& var: vars) var->WriteDaisyMC(*mcDaisyOutDir);
    for(auto& var: vars2D) var->WriteDaisyMC(*mcDaisyOutDir);
  
    //Protons On Target
    mcPOT->Write();

    PlotUtils::TargetUtils targetInfo;
    assert(error_bands["cv"].size() == 1 && "List of error bands must contain a universe named \"cv\" for the flux integral.");

    for(const auto& var: vars)
    {
      //Flux integral only if systematics are being done (temporary solution)
      util::GetFluxIntegral(*error_bands["cv"].front(), var->efficiencyNumerator->hist)->Write((var->GetName() + "_reweightedflux_integrated").c_str());
      //Always use MC number of nucleons for cross section
      auto nNucleons = new TParameter<double>((var->GetName() + "_fiducial_nucleons").c_str(), targetInfo.GetTrackerNNucleons(PlotUtils::TargetProp::Tracker::Face, PlotUtils::TargetProp::Tracker::Back, true, apothem));
      nNucleons->Write();
    }
    for(const auto& var: vars2D)
    {
      //Flux integral only if systematics are being done (temporary solution)
      util::GetFluxIntegral(*error_bands["cv"].front(), var->efficiencyNumerator->hist)->Write((var->GetName() + "_reweightedflux_integrated").c_str());
      //Always use MC number of nucleons for cross section
      auto nNucleons = new TParameter<double>((var->GetName() + "_fiducial_nucleons").c_str(), targetInfo.GetTrackerNNucleons(PlotUtils::TargetProp::Tracker::Face, PlotUtils::TargetProp::Tracker::Back, true, apothem));
      nNucleons->Write();
    }

    //Write data results
    filepath = std::string(DATA_OUT_FILE_NAME)+".root";
    TFile* dataOutDir = TFile::Open(filepath.c_str(), "RECREATE");
    if(!dataOutDir)
    {
      std::cerr << "Failed to open a file named " << filepath << " in the current directory for writing histograms.\n";
      return badOutputFile;
    }

    for(auto& var: vars) var->WriteData(*dataOutDir);
    for(auto& var: vars2D) var->WriteData(*dataOutDir);

    //Protons On Target
    auto dataPOT = new TParameter<double>("POTUsed", options.m_data_pot);
    dataPOT->Write();


    //Write data daisy results
    filepath = std::string(DATA_OUT_FILE_NAME)+"Daisy.root";
    TFile* daisyDataOutDir = TFile::Open(filepath.c_str(), "RECREATE");
    if(!daisyDataOutDir)
    {
      std::cerr << "Failed to open a file named " << filepath << " in the current directory for writing histograms.\n";
      return badOutputFile;
    }

    for(auto& var: vars) var->WriteDaisyData(*daisyDataOutDir);
    for(auto& var: vars2D) var->WriteDaisyData(*daisyDataOutDir);

    //Protons On Target
    dataPOT->Write();


    //Write Migration results
    filepath = std::string(MIGRATION_OUT_FILE_NAME)+".root";
    TFile* migrationOutDir = TFile::Open(filepath.c_str(), "RECREATE");
    if(!migrationOutDir)
    {
      std::cerr << "Failed to open a file named " << filepath << " in the current directory for writing histograms.\n";
      return badOutputFile;
    }
    for(auto& var: vars2D)var->WriteMigration(*migrationOutDir);

    //Write Migration results
    filepath = std::string(MIGRATION_OUT_FILE_NAME)+"Daisy.root";
    TFile* migrationDaisyOutDir = TFile::Open(filepath.c_str(), "RECREATE");
    if(!migrationDaisyOutDir)
    {
      std::cerr << "Failed to open a file named " << filepath << " in the current directory for writing histograms.\n";
      return badOutputFile;
    }
    for(auto& var: vars2D)var->WriteDaisyMigration(*migrationDaisyOutDir);

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