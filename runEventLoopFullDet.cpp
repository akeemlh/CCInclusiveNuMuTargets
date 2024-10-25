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
    if(ANNVtx.size()==3)
    {
      ANNVerticesMC->Fill( ANNVtx[2], cvWeight);
      ANNVerticesGranularMC->Fill( ANNVtx[0], ANNVtx[1], ANNVtx[2], cvWeight);
      if (m_TargetUtils->InWaterTargetVolMC(ANNVtx[0], ANNVtx[1], ANNVtx[2]))
      {
        TrueVtxANNRecoInWater->Fill( TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), cvWeight);
      }
      else
      {
        TrueVtxANNRecoOutWater->Fill( TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), cvWeight);
      }
    }
    TBVerticesMC->Fill(TrackBasedVtx.Z(), cvWeight);
    TBVerticesGranularMC->Fill( TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z(), cvWeight);
    TrueVerticesGranular->Fill( TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), cvWeight);
    if (m_TargetUtils->InWaterTargetVolMC(TrackBasedVtx.X(), TrackBasedVtx.Y(), TrackBasedVtx.Z()))
    {
      TrueVtxTBRecoInWater->Fill( TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), cvWeight);
    }
    else
    {
      std::cout<<"TB: X: " << TrackBasedVtx.X() <<" Y: " << TrackBasedVtx.Y() <<" Z: "<< TrackBasedVtx.Z() <<std::endl;
      TrueVtxTBRecoOutWater->Fill( TrueVtx.X(), TrueVtx.Y(), TrueVtx.Z(), cvWeight);
    }
    std::cout<<"Here6\n";


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