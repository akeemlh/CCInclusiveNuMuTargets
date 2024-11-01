#define MC_OUT_FILE_NAME "runEventLoopTargetsMC.root"
#define DATA_OUT_FILE_NAME "runEventLoopTargetsData.root"

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
#include "util/Variable2DNuke.h"
#include "util/Variable1DNuke.h"
#include "util/GetFluxIntegral.h"
#include "util/GetPlaylist.h"
#include "cuts/SignalDefinition.h"
#include "cuts/q3RecoCut.h"
#include "studies/Study.h"
#include "studies/PerEventVarByGENIELabel2D.h"
#include "studies/WaterTargetIntOrigin2D.h"
#include "util/NukeUtils.h"
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
#include "fstream"

#include "Math/Vector3D.h"
#include "TH3D.h"

//c++ includes
#include <iostream>
#include <cstdlib> //getenv()

bool PseudoTargetBuffer = false; //To exlclude the plane either end of a nuclear target

//Treat the plastic between planes as though they were nuclear targets themselves
int getPlasticPseudoTargetCode(int mod, int plane)
{
  if (PseudoTargetBuffer)
  {
    if (mod == -2 && plane == 2) return -1;
    if (mod == 0 && plane == 1) return -1;
    if (mod == 3 && plane == 2) return -1;
    if (mod == 5 && plane == 1) return -1;
    if (mod == 8 && plane == 2) return -1;
    if (mod == 11 && plane == 1) return -1;
    if (mod == 14 && plane == 2) return -1;
    if (mod == 15 && plane == 1) return -1;
    if (mod == 18 && plane == 2) return -1;
    if (mod == 20 && plane == 1) return -1;
    if (mod == 21 && plane == 2) return -1;
  }
  if (mod >= -5 && mod <= -2)  return 7;
  else if (mod >= 0 && mod <= 3)  return 8;
  else if (mod >= 5 && mod <= 8)  return 9;
  else if (mod >= 11 && mod <= 14)  return 10;
  else if (mod >= 15 && mod <= 18)  return 11;
  else if (mod >= 20 && mod <= 21)  return 12;
  else return -1;
}

//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFillEventSelection(
    PlotUtils::ChainWrapper* chain,
    std::map<std::string, std::vector<CVUniverse*> > error_bands,
    std::vector<Variable1DNuke*> vars,
    std::vector<Variable2DNuke*> vars2D,
    std::vector<Study*> studies,
    PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts,
    PlotUtils::Model<CVUniverse, MichelEvent>& model)
{
  assert(!error_bands["cv"].empty() && "\"cv\" error band is empty!  Can't set Model weight.");
  auto& cvUniv = error_bands["cv"].front();

  std::cout << "Starting MC reco loop...\n";
  const int nEntries = chain->GetEntries();
  //const int nEntries = 10000;
  for (int i=0; i<nEntries; ++i)
  {
    //std::cout << i << " / " << nEntries << "\n";
    if(i%1000==0) std::cout << i << " / " << nEntries << "\n";
    //if(i%1000==0) std::cout << i << " / " << nEntries << "\r" <<std::flush;

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
        
        int pseudoTarget = getPlasticPseudoTargetCode(universe->GetTruthVtxModule(), universe->GetTruthVtxPlane());
        int truthTgtCode = universe->GetTruthTargetCode();
        int effectiveTgtCode = (pseudoTarget==-1) ? truthTgtCode: pseudoTarget;
        //If this has a module num 6 it came from water target
        bool inWaterSegment = (universe->GetTruthVtxModule()==6);
        //if (inWaterSegment) std::cout << "MC says in Water, entry: " << i << std::endl;
        //if (universe->GetANNSegment()==36 || universe->GetANNVtxModule()==6) std::cout << "Reco says in Water, entry: " << i << std::endl;
        //Q: Very rarely we get annTgtCode==1000. What is that? Example in 1A MC playlist file, entry i = 68260, also 1P 947009. Answer: When material is unknown, z is left as 0 so 1000 means target 1 unknown material
        int code = inWaterSegment ? -999 : effectiveTgtCode;

        //To do: use universe->hasMLPred()
        //Nuke Target Study
        const double weight = model.GetWeight(*universe, myevent); //Only calculate the per-universe weight for events that will actually use it.
        for(auto& var: vars)
        {
          if (truthTgtCode>0 || inWaterSegment || code>0) //If this event occurs inside a nuclear target
          {
            //Plot events that occur within the nuclear targets grouped by which target they occur in
            (*var->m_HistsByTgtCodeMC)[code].FillUniverse(universe, var->GetRecoValue(*universe), weight);
            if (util::TgtCodeLabelsNuke.count(code)!=0) (*var->m_intChannelsByTgtCode[code])[universe->GetInteractionType()].FillUniverse(universe, var->GetRecoValue(*universe), weight);
          }
          if (!(truthTgtCode>0) && !inWaterSegment) //If our vertex is not in the nuclear targets
          {
            int tmpModCode = (universe->GetANNVtxModule()*10)+universe->GetANNVtxPlane();
            auto USTgtID = util::USModPlaCodeToTgtId.find(tmpModCode);
            auto DSTgtID = util::DSModPlaCodeToTgtId.find(tmpModCode);
            if (USTgtID != util::USModPlaCodeToTgtId.end()) //Is event reconstructed immediately upstream of a nuclear target
            {
              //std::cout<<"Found in US target "<< USTgtID->second<< "\n";
              //Check where it really interacted
              int tmpTruthModCode = (universe->GetTruthVtxModule()*10)+universe->GetTruthVtxPlane();
              //Checking if this interaction truthfully originated in an upstream or downstream sideband
              //Checking if this interaction truthfully originated in the neigbouring nuclear target
              //If this has a non-zero target code then it actually originated in a nuke target
              //If this has a segment num 36 it came from water target
              int truthTgtID = universe->GetTruthTargetID();
              if (truthTgtID > 0 || inWaterSegment) 
              {
                (*var->m_sidebandHistSetUSMC[USTgtID->second])[2].FillUniverse(universe, var->GetRecoValue(*universe), weight);
              }
              else if (util::isUSPlane(tmpModCode)>0) //If originated immediately upstream of nuke target
              {
                (*var->m_sidebandHistSetUSMC[USTgtID->second])[0].FillUniverse(universe, var->GetRecoValue(*universe), weight);
              }
              else if (util::isDSPlane(tmpModCode)>0) //If originated immediately downstream of nuke target
              {
                (*var->m_sidebandHistSetUSMC[USTgtID->second])[1].FillUniverse(universe, var->GetRecoValue(*universe), weight);
              }
              else //Fill "Other" histogram if this event didn't really have a vtx in the nuke target or sideband
              {
                (*var->m_sidebandHistSetUSMC[USTgtID->second])[-1].FillUniverse(universe, var->GetRecoValue(*universe), weight);
              }
            }
            else if (DSTgtID != util::DSModPlaCodeToTgtId.end()) //Or is event reconstructed immediately downstream of a nuclear target
            {
              //Check where it really interacted
              int tmpTruthModCode = (universe->GetTruthVtxModule()*10)+universe->GetTruthVtxPlane();
              //Checking if this interaction truthfully originated in an upstream or downstream sideband
              //Checking if this interaction truthfully originated in the neigbouring nuclear target
              //If this has a non-zero target code then it actually originated in a nuke target
              //If this has a segment num 36 it came from water target
              int truthTgtID = universe->GetTruthTargetID();
              if (truthTgtID > 0 || inWaterSegment) 
              {
                (*var->m_sidebandHistSetDSMC[DSTgtID->second])[2].FillUniverse(universe, var->GetRecoValue(*universe), weight);
              }
              else if (util::isUSPlane(tmpModCode)>0) //Is event reconstructed immediately upstream of a nuclear target
              {
                (*var->m_sidebandHistSetDSMC[DSTgtID->second])[0].FillUniverse(universe, var->GetRecoValue(*universe), weight);
              }
              else if (util::isDSPlane(tmpModCode)>0) //Is event reconstructed immediately upstream of a nuclear target
              {
                (*var->m_sidebandHistSetDSMC[DSTgtID->second])[1].FillUniverse(universe, var->GetRecoValue(*universe), weight);
              }
              else //Fill "Other" histogram if this event didn't really have a vtx in the nuke target or sideband
              {
                (*var->m_sidebandHistSetDSMC[DSTgtID->second])[-1].FillUniverse(universe, var->GetRecoValue(*universe), weight);
              }
            }
          }
        }
        for(auto& var: vars2D)
        {
          if (truthTgtCode>0 || inWaterSegment || code>0) //If this event occurs inside a nuclear target
          {
            //Plot events that occur within the nuclear targets grouped by which target they occur in
            (*var->m_HistsByTgtCodeMC)[code].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
            if (util::TgtCodeLabelsNuke.count(code)!=0) (*var->m_intChannelsByTgtCode[code])[universe->GetInteractionType()].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
          }
          if (!(truthTgtCode>0) && !inWaterSegment) //If our vertex is not in the nuclear targets
          {
            int tmpModCode = (universe->GetANNVtxModule()*10)+universe->GetANNVtxPlane();
            auto USTgtID = util::USModPlaCodeToTgtId.find(tmpModCode);
            auto DSTgtID = util::DSModPlaCodeToTgtId.find(tmpModCode);
            if (USTgtID != util::USModPlaCodeToTgtId.end()) //Is event reconstructed immediately upstream of a nuclear target
            {
              //std::cout<<"Found in US target "<< USTgtID->second<< "\n";
              //Check where it really interacted
              int tmpTruthModCode = (universe->GetTruthVtxModule()*10)+universe->GetTruthVtxPlane();
              //Checking if this interaction truthfully originated in an upstream or downstream sideband
              //Checking if this interaction truthfully originated in the neigbouring nuclear target
              //If this has a non-zero target code then it actually originated in a nuke target
              //If this has a segment num 36 it came from water target
              int truthTgtID = universe->GetTruthTargetID();
              if (truthTgtID > 0 || inWaterSegment) 
              {
                (*var->m_sidebandHistSetUSMC[USTgtID->second])[2].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
              }
              else if (util::isUSPlane(tmpModCode)>0) //If originated immediately upstream of nuke target
              {
                (*var->m_sidebandHistSetUSMC[USTgtID->second])[0].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
              }
              else if (util::isDSPlane(tmpModCode)>0) //If originated immediately downstream of nuke target
              {
                (*var->m_sidebandHistSetUSMC[USTgtID->second])[1].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
              }
              else //Fill "Other" histogram if this event didn't really have a vtx in the nuke target or sideband
              {
                (*var->m_sidebandHistSetUSMC[USTgtID->second])[-1].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
              }
            }
            else if (DSTgtID != util::DSModPlaCodeToTgtId.end()) //Or is event reconstructed immediately downstream of a nuclear target
            {
              //Check where it really interacted
              int tmpTruthModCode = (universe->GetTruthVtxModule()*10)+universe->GetTruthVtxPlane();
              //Checking if this interaction truthfully originated in an upstream or downstream sideband
              //Checking if this interaction truthfully originated in the neigbouring nuclear target
              //If this has a non-zero target code then it actually originated in a nuke target
              //If this has a segment num 36 it came from water target
              int truthTgtID = universe->GetTruthTargetID();
              if (truthTgtID > 0 || inWaterSegment) 
              {
                (*var->m_sidebandHistSetDSMC[DSTgtID->second])[2].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
              }
              else if (util::isUSPlane(tmpModCode)>0) //Is event reconstructed immediately upstream of a nuclear target
              {
                (*var->m_sidebandHistSetDSMC[DSTgtID->second])[0].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
              }
              else if (util::isDSPlane(tmpModCode)>0) //Is event reconstructed immediately upstream of a nuclear target
              {
                (*var->m_sidebandHistSetDSMC[DSTgtID->second])[1].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
              }
              else //Fill "Other" histogram if this event didn't really have a vtx in the nuke target or sideband
              {
                (*var->m_sidebandHistSetDSMC[DSTgtID->second])[-1].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
              }
            }
          }
        }
        const bool isSignal = michelcuts.isSignal(*universe, weight);

        if(isSignal)
        {
          for(auto& study: studies) study->SelectedSignal(*universe, myevent, weight);
          for(auto& var: vars)
          {
            //Cross section components
            if (truthTgtCode>0 || inWaterSegment || code>0) //If this event occurs inside a nuclear target
            {
              //Plot events that occur within the nuclear targets grouped by which target they occur in
              (*var->m_HistsByTgtCodeEfficiencyNumerator)[code].FillUniverse(universe, var->GetRecoValue(*universe), weight);
              (*var->m_HistsByTgtCodeMigration)[code].FillUniverse(universe, var->GetRecoValue(*universe), var->GetTrueValue(*universe), weight);
            }
          }

          for(auto& var: vars2D)
          {
            //Cross section components
            if (truthTgtCode>0 || inWaterSegment || code>0) //If this event occurs inside a nuclear target
            {
              //Plot events that occur within the nuclear targets grouped by which target they occur in
              (*var->m_HistsByTgtCodeEfficiencyNumerator)[code].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
              //var->FillMigration(universe, code, var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
              (*var->m_HistsByTgtCodeMigration)[code].Fill(var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), var->GetTrueValueX(*universe), var->GetTrueValueY(*universe),universe->ShortName(),(univCount-1),weight);
              //(var->m_migration)[code].Fill(var->GetRecoValueX(*universe), var->GetRecoValueY(*universe),var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);

            }
          }
        }
        else
        {
          int bkgd_ID = -1;
          if (universe->GetCurrent()==2)bkgd_ID=0;
          else bkgd_ID=1;

          //Cross section components
          int ANNTgtCode = universe->GetANNTargetCode();
          int TruthTgtCode = universe->GetTruthTargetCode();
          //std::cout<<"ANNTgt " <<ANNTgtCode << " TruthTgtCode " << TruthTgtCode << std::endl;
          if (ANNTgtCode!=TruthTgtCode || (universe->GetTruthVtxModule()==6 ^ (universe->GetANNSegment()==36 || universe->GetANNVtxModule()==6)))
          {
            bkgd_ID=2;
            //std::cout<<"Wrong material background, ANN target code " << ANNTgtCode << " truth target code " << TruthTgtCode << std::endl;
          }
          /* else
          {
            std::cout<<"ANN and truth agree on target code " << ANNTgtCode <<std::endl;
          } */

          if (truthTgtCode>0 || inWaterSegment || code>0) //If this event occurs inside a nuclear target
          {
            for(auto& var: vars) (*var->m_bkgsByTgtCode[code])[bkgd_ID].FillUniverse(universe, var->GetRecoValue(*universe), weight);
            for(auto& var: vars2D) (*var->m_bkgsByTgtCode[code])[bkgd_ID].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
          }
        }
      } // End band's universe loop
    } // End Band loop
  } //End entries loop
  std::cout << "Finished MC reco loop.\n";
}

void LoopAndFillData( PlotUtils::ChainWrapper* data,
			        std::vector<CVUniverse*> data_band,
				std::vector<Variable1DNuke*> vars,
                                std::vector<Variable2DNuke*> vars2D,
                                std::vector<Study*> studies,
				PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts)
{
  std::cout << "Starting data loop...\n";
  //const int nEntries = 10000;
  const int nEntries = data->GetEntries();
  for (int i=0; i<data->GetEntries(); ++i) {
    for (auto universe : data_band) {
      universe->SetEntry(i);
      if(i%1000==0) std::cout << i << " / " << nEntries << "\r" << std::flush;
      MichelEvent myevent; 
      if (!michelcuts.isDataSelected(*universe, myevent).all()) continue;

      int pseudoTarget = getPlasticPseudoTargetCode(universe->GetANNVtxModule(), universe->GetANNVtxPlane());
      int annTgtCode = universe->GetANNTargetCode();
      int effectiveTgtCode = (pseudoTarget==-1) ? annTgtCode: pseudoTarget;
      //If this has a segment num 36 it came from water target
      bool inWaterSegment = (universe->GetANNSegment()==36 || universe->GetANNVtxModule()==6);
      //Q: Very rarely we get annTgtCode==1000. What is that? Example in 1A MC playlist file, entry i = 68260, also 1P 947009. Answer: When material is unknown, z is left as 0 so 1000 means target 1 unknown material
      int code = inWaterSegment ? -999 : effectiveTgtCode;


      for(auto& study: studies) study->Selected(*universe, myevent, 1); 
      for(auto& var: vars)
      {
        if (annTgtCode>0 || inWaterSegment || code>0 ) //If this event occurs inside a nuclear target
        {
          //Plot events that occur within the nuclear targets grouped by which target they occur in
          (*var->m_HistsByTgtCodeData)[code].FillUniverse(universe, var->GetRecoValue(*universe, myevent.m_idx), 1);
        }
        if (!(annTgtCode>0) && !inWaterSegment) //If our vertex is not in the nuclear targets
        {
          int tmpModCode = (universe->GetANNVtxModule()*10)+universe->GetANNVtxPlane();
          int USModNum = util::isUSPlane(tmpModCode);
          int DSModNum = util::isDSPlane(tmpModCode);
          if (USModNum>0) //Is event reconstructed immediately upstream of a nuclear target
          {
            (*var->m_sidebandHistsUSData)[USModNum].FillUniverse(universe, var->GetRecoValue(*universe), 1);
          }
          else if (DSModNum>0) //Or is event reconstructed immediately downstream of a nuclear target
          {
            (*var->m_sidebandHistsDSData)[DSModNum].FillUniverse(universe, var->GetRecoValue(*universe), 1);
          }
        }
      }
      //Nuke Target Study
      for(auto& var: vars2D)
      {
        if (annTgtCode>0 || inWaterSegment || code>0 ) //If this event occurs inside a nuclear target
        {
          //Plot events that occur within the nuclear targets grouped by which target they occur in
          (*var->m_HistsByTgtCodeData)[code].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), 1);
        }
        if (!(annTgtCode>0) && !inWaterSegment) //If our vertex is not in the nuclear targets
        {
          int tmpModCode = (universe->GetANNVtxModule()*10)+universe->GetANNVtxPlane();
          int USModNum = util::isUSPlane(tmpModCode);
          int DSModNum = util::isDSPlane(tmpModCode);
          if (USModNum>0) //Is event reconstructed immediately upstream of a nuclear target
          {
            (*var->m_sidebandHistsUSData)[USModNum].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), 1);
          }
          else if (DSModNum>0) //Or is event reconstructed immediately downstream of a nuclear target
          {
            (*var->m_sidebandHistsDSData)[DSModNum].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), 1);
          }
        }
      }
    }
  }
  std::cout << "Finished data loop.\n";
}

void LoopAndFillEffDenom( PlotUtils::ChainWrapper* truth,
    				std::map<std::string, std::vector<CVUniverse*> > truth_bands,
    				std::vector<Variable1DNuke*> vars,
                                std::vector<Variable2DNuke*> vars2D,
    				PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts,
                                PlotUtils::Model<CVUniverse, MichelEvent>& model)
{
  assert(!truth_bands["cv"].empty() && "\"cv\" error band is empty!  Could not set Model entry.");
  auto& cvUniv = truth_bands["cv"].front();

  std::cout << "Starting efficiency denominator loop...\n";
  //const int nEntries = 10000;
  const int nEntries = truth->GetEntries();
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

        int pseudoTarget = getPlasticPseudoTargetCode(universe->GetTruthVtxModule(), universe->GetTruthVtxPlane()); ///Wrong!!!!!!!!!!!!!!!!11
        int truthTgtCode = universe->GetTruthTargetCode();
        int effectiveTgtCode = (pseudoTarget==-1) ? truthTgtCode: pseudoTarget;
        int truthTgtID = universe->GetTruthTargetID();
        //If this has a segment num 36 it came from water target
        bool inWaterSegment = (truthTgtID==6);
        int code = inWaterSegment ? -999 : effectiveTgtCode;
        //Fill efficiency denominator now: 
        for(auto var: vars)
        {
          (*var->m_HistsByTgtCodeEfficiencyDenominator)[code].FillUniverse(universe, var->GetTrueValue(*universe), weight);
        }

        for(auto var: vars2D)
        {
          (*var->m_HistsByTgtCodeEfficiencyDenominator)[code].FillUniverse(universe, var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
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
  std:: cout << reco_tree_name << std::endl;

  std::string playlistname = "minervame1A";
  size_t mc_label = mc_file_list.find("/MC/");
  if (mc_label == string::npos) std::cout << "FAILED TO DETERMINE WHICH PLAYLIST THIS IS. DEFAULTING TO ME1A FOR FLUX\n";
  else if (mc_file_list.substr(mc_label+4, 4) == "Test" ) std::cout << "IDENTIFIED TEST PLAYLIST. DEFAULTING TO ME1A FOR FLUX\n";
  else
  {
    std::string playlistkey = mc_file_list.substr(mc_label+4, 2);
    std::cout<<"Identified Playlist " << playlistkey << std::endl;
    playlistname = "minervame"+playlistkey;
  }
  
  std::cout<<"Getting playlist\n";
  std::cout << mc_file_list.substr(mc_label+4, 4)<<std::endl;
  std::cout << mc_file_list.substr(mc_label+4, 2)<<std::endl;
  std::cout << playlistname << std::endl;
  //const bool is_grid = false; //TODO: Are we going to put this back?  Gonzalo needs it iirc.
  PlotUtils::MacroUtil options(reco_tree_name, mc_file_list, data_file_list, playlistname, true);
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
  PlotUtils::Cutter<CVUniverse, MichelEvent>::reco_t nukeSidebands, nukePreCut;
  PlotUtils::Cutter<CVUniverse, MichelEvent>::truth_t nukeSignalDefinition, nukePhaseSpace;

  const double apothem = 850; //All in mm
  nukePreCut.emplace_back(new reco::ZRange<CVUniverse, MichelEvent>("Nuclear Targets Z pos", PlotUtils::TargetProp::NukeRegion::Face, PlotUtils::TargetProp::NukeRegion::Back));
  nukePreCut.emplace_back(new reco::Apothem<CVUniverse, MichelEvent>(apothem));
  nukePreCut.emplace_back(new reco::MaxMuonAngle<CVUniverse, MichelEvent>(17.));
  nukePreCut.emplace_back(new reco::HasMINOSMatch<CVUniverse, MichelEvent>());
  nukePreCut.emplace_back(new reco::NoDeadtime<CVUniverse, MichelEvent>(1, "Deadtime"));
  nukePreCut.emplace_back(new reco::IsNeutrino<CVUniverse, MichelEvent>());
  nukePreCut.emplace_back(new reco::MuonEnergyMin<CVUniverse, MichelEvent>(2000.0, "EMu Min"));
  nukePreCut.emplace_back(new reco::MuonEnergyMax<CVUniverse, MichelEvent>(50000.0, "EMu Max"));
  nukePreCut.emplace_back(new reco::ANNConfidenceCut<CVUniverse, MichelEvent>(0.20));

  //nukeSidebands.emplace_back(new reco::ZRange<CVUniverse, MichelEvent>("Test sideband z pos", 0, 1000000000000.0));
  //nukeSidebands.emplace_back(new reco::USScintillator<CVUniverse, MichelEvent>());
  //nukeSidebands.emplace_back(new reco::DSScintillator<CVUniverse, MichelEvent>());

                                                                                                                   
  nukeSignalDefinition.emplace_back(new truth::IsNeutrino<CVUniverse>());
  nukeSignalDefinition.emplace_back(new truth::IsCC<CVUniverse>());
                                                        
  nukePhaseSpace.emplace_back(new truth::ZRange<CVUniverse>("Nuclear Targets Z pos", PlotUtils::TargetProp::NukeRegion::Face, PlotUtils::TargetProp::NukeRegion::Back));
  nukePhaseSpace.emplace_back(new truth::Apothem<CVUniverse>(apothem));
  nukePhaseSpace.emplace_back(new truth::MuonAngle<CVUniverse>(17.));
  nukePhaseSpace.emplace_back(new truth::MuonEnergyMin<CVUniverse>(2000.0, "EMu Min"));
  nukePhaseSpace.emplace_back(new truth::MuonEnergyMax<CVUniverse>(50000.0, "EMu Max"));
  //nukePhaseSpace.emplace_back(new truth::PZMuMin<CVUniverse>(1500.));
                                                                                                                                                   
  PlotUtils::Cutter<CVUniverse, MichelEvent> nukeCuts(std::move(nukePreCut), std::move(nukeSidebands) , std::move(nukeSignalDefinition),std::move(nukePhaseSpace));

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

  std::vector<Variable1DNuke*> nukeVars;
  std::vector<Variable2DNuke*> nukeVars2D;

  nukeVars.push_back(new Variable1DNuke("nuke_pTmu", "p_{T, #mu} [GeV/c]", dansPTBins, &CVUniverse::GetMuonPT, &CVUniverse::GetMuonPTTrue));
  nukeVars.push_back(new Variable1DNuke("nuke_pZmu", "p_{||, #mu} [GeV/c]", dansPzBins, &CVUniverse::GetMuonPz, &CVUniverse::GetMuonPzTrue));
  nukeVars.push_back(new Variable1DNuke("nuke_Emu", "E_{#mu} [GeV]", robsEmuBins, &CVUniverse::GetEmuGeV, &CVUniverse::GetElepTrueGeV));
  nukeVars.push_back(new Variable1DNuke("nuke_Erecoil", "E_{recoil}", robsRecoilBins, &CVUniverse::GetRecoilE, &CVUniverse::Getq0True)); //TODO: q0 is not the same as recoil energy without a spline correction
  nukeVars.push_back(new Variable1DNuke("nuke_BjorkenX", "X", bjorkenXbins, &CVUniverse::GetBjorkenX, &CVUniverse::GetBjorkenXTrue));
  nukeVars2D.push_back(new Variable2DNuke("nuke_pTmu_pZmu", *nukeVars[1], *nukeVars[0]));

  std::vector<Study*> studies;
  std::function<double(const CVUniverse&, const MichelEvent&)> ptmu = [](const CVUniverse& univ, const MichelEvent& /* evt */) { return univ.GetMuonPT();};
  std::function<double(const CVUniverse&, const MichelEvent&)> pzmu = [](const CVUniverse& univ, const MichelEvent& /* evt */) { return univ.GetMuonPz();};

  studies.push_back(new PerEventVarByGENIELabel2D(pzmu, ptmu, std::string("pzmu_vs_ptmu_GENIE_labels"), std::string("GeV/c"), dansPzBins, dansPTBins, error_bands));
  studies.push_back(new WaterTargetIntOrigin2D(pzmu, ptmu, std::string("pzmu_vs_ptmu_water_breakdown"), std::string("GeV/c"), dansPzBins, dansPTBins, error_bands));

  CVUniverse* data_universe = new CVUniverse(options.m_data);
  std::vector<CVUniverse*> data_band = {data_universe};
  std::map<std::string, std::vector<CVUniverse*> > data_error_bands;
  data_error_bands["cv"] = data_band;
  
  std::vector<Study*> data_studies;
  //data_studies.push_back(new PerEventVarByGENIELabel2D(ptmu, pzmu, std::string("ptmu_vs_pzmu"), std::string("GeV/c"), dansPTBins, dansPzBins, data_error_bands));
  //Wouldn't make sense to do a PerEventVarByGENIELabel2D study for data since data wont have the GENIE simulation labels

  for(auto& var: nukeVars) var->InitializeMCHists(error_bands, truth_bands);
  for(auto& var: nukeVars) var->InitializeDATAHists(data_band);

  for(auto& var: nukeVars2D) var->InitializeMCHists(error_bands, truth_bands);
  for(auto& var: nukeVars2D) var->InitializeDATAHists(data_band);

  // Loop entries and fill
  try
  {
    CVUniverse::SetTruth(false);
    LoopAndFillEventSelection(options.m_mc, error_bands, nukeVars, nukeVars2D, studies, nukeCuts, model);
    CVUniverse::SetTruth(true);
    LoopAndFillEffDenom(options.m_truth, truth_bands, nukeVars, nukeVars2D, nukeCuts, model);
    options.PrintMacroConfiguration(argv[0]);
    std::cout << "Nuclear Target MC cut summary:\n" << nukeCuts << "\n";
    nukeCuts.resetStats();

    CVUniverse::SetTruth(false);
    LoopAndFillData(options.m_data, data_band, nukeVars, nukeVars2D, data_studies, nukeCuts);
    std::cout << "Nuclear Target Data cut summary:\n" << nukeCuts << "\n";



    //Write MC results
    TFile* mcOutDir = TFile::Open(MC_OUT_FILE_NAME, "RECREATE");
    if(!mcOutDir)
    {
      std::cerr << "Failed to open a file named " << MC_OUT_FILE_NAME << " in the current directory for writing histograms.\n";
      return badOutputFile;
    }
    std::cout<<"Saving " << studies.size()<< " studies\n";
    for(auto& study: studies) study->SaveOrDraw(*mcOutDir);
    std::cout<<"Saved studies\n";

    for(auto& var: nukeVars) var->WriteMC(*mcOutDir);
    for(auto& var: nukeVars2D) var->WriteMC(*mcOutDir);

    //Protons On Target
    auto mcPOT = new TParameter<double>("POTUsed", options.m_mc_pot);
    mcPOT->Write();

    PlotUtils::TargetUtils targetInfo;
    assert(error_bands["cv"].size() == 1 && "List of error bands must contain a universe named \"cv\" for the flux integral.");

    for(auto& var: nukeVars) 
    {
      //Flux integral only if systematics are being done (temporary solution)
      //Always use MC number of nucleons for cross section
      std::vector<int> targetCodes = {1026, 1082, 2026, 2082, 3006, 3026, 3082, 4082, 5026, 5082};
      for (int code : targetCodes)
      {
        int tgtZ = code%1000;
        int tgtID = (code-tgtZ)/1000;
        auto nNucleons = new TParameter<double>((var->GetName() +"_"+std::to_string(code)+"_fiducial_nucleons").c_str(), targetInfo.GetPassiveTargetNNucleons( tgtID, tgtZ, true));
        nNucleons->Write();
        util::GetFluxIntegral(*error_bands["cv"].front(), (*var->m_HistsByTgtCodeEfficiencyNumerator)[code].hist)->Write((var->GetName()+ "_"+std::to_string(code) + "_reweightedflux_integrated").c_str());

      }
      //For water
      auto nNucleons = new TParameter<double>((var->GetName() + "_TargetWater_fiducial_nucleons").c_str(), targetInfo.GetPassiveTargetNNucleons( 6, 1, true));
      nNucleons->Write();
      util::GetFluxIntegral(*error_bands["cv"].front(), (*var->m_HistsByTgtCodeEfficiencyNumerator)[-999].hist)->Write((var->GetName()+ "_TargetWater_reweightedflux_integrated").c_str());


      //For Pseudotargets
      if (PseudoTargetBuffer)
      {
        //Target 7
        auto nNucleonsTgt7 = new TParameter<double>((var->GetName() + "_Target7_fiducial_nucleons").c_str(), targetInfo.GetTrackerNNucleons( 7, true));
        nNucleonsTgt7->Write();
        util::GetFluxIntegral(*error_bands["cv"].front(), (*var->m_HistsByTgtCodeEfficiencyNumerator)[7].hist)->Write((var->GetName() + "_Target7_reweightedflux_integrated").c_str());
        //Target 8
        auto nNucleonsTgt8 = new TParameter<double>((var->GetName() + "_Target8_fiducial_nucleons").c_str(), targetInfo.GetTrackerNNucleons( 6, true));
        nNucleonsTgt8->Write();
        util::GetFluxIntegral(*error_bands["cv"].front(), (*var->m_HistsByTgtCodeEfficiencyNumerator)[8].hist)->Write((var->GetName() + "_Target8_reweightedflux_integrated").c_str());
        //Target 9
        auto nNucleonsTgt9 = new TParameter<double>((var->GetName() + "_Target9_fiducial_nucleons").c_str(), targetInfo.GetTrackerNNucleons( 6, true));
        nNucleonsTgt9->Write();
        util::GetFluxIntegral(*error_bands["cv"].front(), (*var->m_HistsByTgtCodeEfficiencyNumerator)[9].hist)->Write((var->GetName() + "_Target9_reweightedflux_integrated").c_str());
        //Target 10
        auto nNucleonsTgt10 = new TParameter<double>((var->GetName() + "_Target10_fiducial_nucleons").c_str(), targetInfo.GetTrackerNNucleons( 6, true));
        nNucleonsTgt10->Write();
        util::GetFluxIntegral(*error_bands["cv"].front(), (*var->m_HistsByTgtCodeEfficiencyNumerator)[10].hist)->Write((var->GetName() + "_Target10_reweightedflux_integrated").c_str());
        //Target 11
        auto nNucleonsTgt11 = new TParameter<double>((var->GetName() + "_Target11_fiducial_nucleons").c_str(), targetInfo.GetTrackerNNucleons( 6, true));
        nNucleonsTgt11->Write();
        util::GetFluxIntegral(*error_bands["cv"].front(), (*var->m_HistsByTgtCodeEfficiencyNumerator)[11].hist)->Write((var->GetName() + "_Target11_reweightedflux_integrated").c_str());
        //Target 12
        auto nNucleonsTgt12 = new TParameter<double>((var->GetName() + "_Target12_fiducial_nucleons").c_str(), targetInfo.GetTrackerNNucleons( 2, true));
        nNucleonsTgt12->Write();
        util::GetFluxIntegral(*error_bands["cv"].front(), (*var->m_HistsByTgtCodeEfficiencyNumerator)[12].hist)->Write((var->GetName() + "_Target12_reweightedflux_integrated").c_str());

      }
      else
      {
        //Target 7
        auto nNucleonsTgt7 = new TParameter<double>((var->GetName() + "_Target7_fiducial_nucleons").c_str(), targetInfo.GetTrackerNNucleons( 8, true));
        nNucleonsTgt7->Write();
        util::GetFluxIntegral(*error_bands["cv"].front(), (*var->m_HistsByTgtCodeEfficiencyNumerator)[7].hist)->Write((var->GetName() + "_Target7_reweightedflux_integrated").c_str());
        //Target 8
        auto nNucleonsTgt8 = new TParameter<double>((var->GetName() + "_Target8_fiducial_nucleons").c_str(), targetInfo.GetTrackerNNucleons( 8, true));
        nNucleonsTgt8->Write();
        util::GetFluxIntegral(*error_bands["cv"].front(), (*var->m_HistsByTgtCodeEfficiencyNumerator)[8].hist)->Write((var->GetName() + "_Target8_reweightedflux_integrated").c_str());
        //Target 9
        auto nNucleonsTgt9 = new TParameter<double>((var->GetName() + "_Target9_fiducial_nucleons").c_str(), targetInfo.GetTrackerNNucleons( 8, true));
        nNucleonsTgt9->Write();
        util::GetFluxIntegral(*error_bands["cv"].front(), (*var->m_HistsByTgtCodeEfficiencyNumerator)[9].hist)->Write((var->GetName() + "_Target9_reweightedflux_integrated").c_str());
        //Target 10
        auto nNucleonsTgt10 = new TParameter<double>((var->GetName() + "_Target10_fiducial_nucleons").c_str(), targetInfo.GetTrackerNNucleons( 8, true));
        nNucleonsTgt10->Write();
        util::GetFluxIntegral(*error_bands["cv"].front(), (*var->m_HistsByTgtCodeEfficiencyNumerator)[10].hist)->Write((var->GetName() + "_Target10_reweightedflux_integrated").c_str());
        //Target 11
        auto nNucleonsTgt11 = new TParameter<double>((var->GetName() + "_Target11_fiducial_nucleons").c_str(), targetInfo.GetTrackerNNucleons( 8, true));
        nNucleonsTgt11->Write();
        util::GetFluxIntegral(*error_bands["cv"].front(), (*var->m_HistsByTgtCodeEfficiencyNumerator)[11].hist)->Write((var->GetName() + "_Target11_reweightedflux_integrated").c_str());
        //Target 12
        auto nNucleonsTgt12 = new TParameter<double>((var->GetName() + "_Target12_fiducial_nucleons").c_str(), targetInfo.GetTrackerNNucleons( 4, true));
        nNucleonsTgt12->Write();
        util::GetFluxIntegral(*error_bands["cv"].front(), (*var->m_HistsByTgtCodeEfficiencyNumerator)[12].hist)->Write((var->GetName() + "_Target12_reweightedflux_integrated").c_str());

      }


    }

    //Write data results
    TFile* dataOutDir = TFile::Open(DATA_OUT_FILE_NAME, "RECREATE");
    if(!dataOutDir)
    {
      std::cerr << "Failed to open a file named " << DATA_OUT_FILE_NAME << " in the current directory for writing histograms.\n";
      return badOutputFile;
    }

    for(auto& var: nukeVars) var->WriteData(*dataOutDir);

    for(auto& var: nukeVars2D) var->WriteData(*dataOutDir);

    for(auto& study: data_studies) study->SaveOrDraw(*dataOutDir);

    //Protons On Target
    auto dataPOT = new TParameter<double>("POTUsed", options.m_data_pot);
    dataPOT->Write();

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

//To Do: Run this twice but also for tracker producing an separate root file