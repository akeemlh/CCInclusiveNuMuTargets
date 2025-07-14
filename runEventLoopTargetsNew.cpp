#define MC_OUT_FILE_NAME_BASE "runEventLoopTargetsMC"
#define DATA_OUT_FILE_NAME_BASE "runEventLoopTargetsData"
#define MIGRATION_2D_OUT_FILE_NAME_BASE "runEventLoopTargets2DMigration"

#define USAGE \
"\n*** USAGE ***\n"\
"runEventLoop <dataPlaylist.txt> <mcPlaylist.txt> <optional target code>\n\n"\
"for water target use 6000 for target code\n"\
"for pseudo targets use 7-12 for target code\n"\
"for other targets, provide the target number and z of material of interest in format num*1000+material z.\n"\
"for other targets, provide the target number and z of material of interest in format num*1000+material z.\n"\
"E.G. for target 2 iron: runEventLoop <dataPlaylist.txt> <mcPlaylist.txt> 2026.\n"\
"E.G. for water target: runEventLoop <dataPlaylist.txt> <mcPlaylist.txt> 6000.\n"\
"E.G. for pseudotarget 8: runEventLoop <dataPlaylist.txt> <mcPlaylist.txt> 8.\n"\
"If target code is not provided, default behaviour is to run over everything, this can be very memory intensive\n"\
"Additionally, if running on the grid, it is possibly more time efficient to only run 1 or a small number of\n"\
"targets here and distribute the targets you're interested in over several nodes\n\n"\
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
#include "util/Variable2DNukeNew.h"
#include "util/Variable1DNukeNew.h"
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
#include "PlotUtils/LowQ2PiReweighter.h"
#include "PlotUtils/AMUDISReweighter.h"
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

bool PseudoTargetBuffer = true; //To exlclude the plane immediately after either end of a nuclear target

//Treat the plastic between planes as though they were nuclear targets themselves
//
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

double GetTruthSegment(CVUniverse* universe)
{
  return util::planecode(universe->GetTruthVtxModule(), universe->GetTruthVtxPlane());
}

int getTgtCode(CVUniverse* universe, int mode /*mc = 0, ANN = 1, TB = 2*/)
{
  double vtx_x, vtx_y, vtx_z;
  int mod, plane;
  if (mode == 0) //Truth
  {
    mod = universe->GetTruthVtxModule();
    plane = universe->GetTruthVtxPlane();
    int pseudoTarget = getPlasticPseudoTargetCode(mod, plane);
    if (pseudoTarget != -1) return pseudoTarget;
    ROOT::Math::XYZTVector Vtx = universe->GetTrueVertex();
    vtx_x = Vtx.X();
    vtx_y = Vtx.Y();
    vtx_z = Vtx.Z();
  } 
  else if (mode == 1) //ANN
  {
    mod = universe->GetANNVtxModule();
    plane = universe->GetANNVtxPlane();
    int pseudoTarget = getPlasticPseudoTargetCode(mod, plane);
    if (pseudoTarget != -1) return pseudoTarget;
    ROOT::Math::XYZVector Vtx = universe->GetANNVertex();
    vtx_x = Vtx.X();
    vtx_y = Vtx.Y();
    vtx_z = Vtx.Z();
  }
  else if (mode == 2) //TB
  {
    mod = universe->GetMADVtxModule();
    plane = universe->GetMADVtxPlane();
    int pseudoTarget = getPlasticPseudoTargetCode(mod, plane);
    if (pseudoTarget != -1) return pseudoTarget;
    ROOT::Math::XYZTVector Vtx = universe->GetVertex();
    vtx_x = Vtx.X();
    vtx_y = Vtx.Y();
    vtx_z = Vtx.Z();
  }
  return util::getTargetCodeFromVtxInfo(vtx_x, vtx_y, vtx_z, mod, plane);
}

int getSidebandTgtCode(CVUniverse* universe, int mode /*mc = 0, ANN = 1, TB = 2*/, int USorDS /*Some planes are both the US planes of one target but the DS planes of another and so need to be treated twice*/)
{
  PlotUtils::TargetUtils tgtUtil;
  //tgtUtil.SetDistToDivCut( ??? ); //What was Anezka's + what is it in GenieXSecExtract/ Should this be tweaked?
  double vtx_x, vtx_y, vtx_z;
  int mod, plane;
  if (mode == 0) //Truth
  {
    mod = universe->GetTruthVtxModule();
    plane = universe->GetTruthVtxPlane();
    ROOT::Math::XYZTVector Vtx = universe->GetTrueVertex();
    vtx_x = Vtx.X();
    vtx_y = Vtx.Y();
    vtx_z = Vtx.Z();
  } 
  else if (mode == 1) //ANN
  {
    mod = universe->GetANNVtxModule();
    plane = universe->GetANNVtxPlane();
    ROOT::Math::XYZVector Vtx = universe->GetANNVertex();
    vtx_x = Vtx.X();
    vtx_y = Vtx.Y();
    vtx_z = Vtx.Z();
  }
  else if (mode == 2) //TB
  {
    mod = universe->GetMADVtxModule();
    plane = universe->GetMADVtxPlane();
    ROOT::Math::XYZTVector Vtx = universe->GetVertex();
    vtx_x = Vtx.X();
    vtx_y = Vtx.Y();
    vtx_z = Vtx.Z();
  }

  //Removing planes immediately up and downstream of targets
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

  //if (mod >= -5 && mod <= -2)  return 7; //Not using upstream of target 1 because of significant rock muon contamination
  /* if (mod == 17 && plane == 2)
 { std::cout<< "mod: " << mod << "\n";

  std::cout<< "plane: " << plane << "\n";

  std::cout<< "vtx_x: " << vtx_x << "\n";

  std::cout<< "vtx_y: " << vtx_y << "\n";

  std::cout<< "vtx_z: " << vtx_z << "\n";} */
  if (mod >= 0 && mod <= 3) //DS of target 1 (Iron and Lead) and US of target 2 (Iron and Lead)
  {
    if (USorDS) //If true we're looking for US Planes therefore we're looking at US of tgt2
    {
      if (tgtUtil.InIron2VolMC( vtx_x, vtx_y, tgtUtil.GetTarget2CenterZMC(), 850., true )) return 2026;
      else if (tgtUtil.InLead2VolMC( vtx_x, vtx_y, tgtUtil.GetTarget2CenterZMC(), 850., true )) return 2082;
    }
    else //If false we're looking for DS Planes therefore we're looking at DS of tgt1
    { 
      if (tgtUtil.InIron1VolMC( vtx_x, vtx_y, tgtUtil.GetTarget1CenterZMC(), 850., true )) return 1026;
      else if (tgtUtil.InLead1VolMC( vtx_x, vtx_y, tgtUtil.GetTarget1CenterZMC(), 850., true )) return 1082;
    }
  } 
  if (mod >= 5 && mod <= 8) //DS of target 2 (Iron and Lead) and US of target 3 (Carbon, iron and Lead)
  {
    if (USorDS) //If true we're looking for US Planes therefore we're looking at US of tgt3
    {
      if (tgtUtil.InCarbon3VolMC( vtx_x, vtx_y, tgtUtil.GetTarget3CenterZMC(), 850., true )) return 3006;
      else if (tgtUtil.InIron3VolMC( vtx_x, vtx_y, tgtUtil.GetTarget3CenterZMC(), 850., true )) return 3026;
      else if (tgtUtil.InLead3VolMC( vtx_x, vtx_y, tgtUtil.GetTarget3CenterZMC(), 850., true )) return 3082;
    }
    else //If false we're looking for DS Planes therefore we're looking at DS of tgt2
    { 
      if (tgtUtil.InIron2VolMC( vtx_x, vtx_y, tgtUtil.GetTarget2CenterZMC(), 850., true )) return 2026;
      else if (tgtUtil.InLead2VolMC( vtx_x, vtx_y, tgtUtil.GetTarget2CenterZMC(), 850., true )) return 2082;
    }
  } 
  if (mod >= 11 && mod <= 14) //DS of target 3 (Carbon, iron and Lead) and US of water target
  {
    if (USorDS) //If true we're looking for US Planes therefore we're looking at US of water
    {
      if (tgtUtil.InWaterTargetVolMC( vtx_x, vtx_y, (PlotUtils::TargetProp::WaterTarget::Face+PlotUtils::TargetProp::WaterTarget::Back)/2, 850. )) return 6000;
    }
    else //If false we're looking for DS Planes therefore we're looking at DS of tgt3
    { 
      if (tgtUtil.InCarbon3VolMC( vtx_x, vtx_y, tgtUtil.GetTarget3CenterZMC(), 850., true )) return 3006;
      else if (tgtUtil.InIron3VolMC( vtx_x, vtx_y, tgtUtil.GetTarget3CenterZMC(), 850., true )) return 3026;
      else if (tgtUtil.InLead3VolMC( vtx_x, vtx_y, tgtUtil.GetTarget3CenterZMC(), 850., true )) return 3082;
    }
  } 
  if (mod >= 15 && mod <= 18) //DS of water target and US of target 4 ( Lead)
  {
    if (USorDS) //If true we're looking for US Planes therefore we're looking at US of target 4
    {
      if (tgtUtil.InLead4VolMC( vtx_x, vtx_y, tgtUtil.GetTarget4CenterZMC(), 850. )) return 4082;
     
    }
    else //If false we're looking for DS Planes therefore we're looking at DS of water target
    { 
      if (tgtUtil.InWaterTargetVolMC( vtx_x, vtx_y, (PlotUtils::TargetProp::WaterTarget::Face+PlotUtils::TargetProp::WaterTarget::Back)/2, 850. )) return 6000;
    }
  } 
  if (mod >= 20 && mod <= 21) //DS of target 4 (Lead) and US of target 5 (Iron and Lead)
  {
    if (USorDS) //If true we're looking for US Planes therefore we're looking at US of tgt3
    {
      if (tgtUtil.InIron5VolMC( vtx_x, vtx_y, tgtUtil.GetTarget5CenterZMC(), 850., true )) return 5026;
      else if (tgtUtil.InLead5VolMC( vtx_x, vtx_y, tgtUtil.GetTarget5CenterZMC(), 850., true )) return 5082;
    }
    else //If false we're looking for DS Planes therefore we're looking at DS of tgt4
    { 
      if (tgtUtil.InLead4VolMC( vtx_x, vtx_y, tgtUtil.GetTarget4CenterZMC(), 850. )) return 4082;
    }
  }
  if (mod >= 20 && mod <= 21) //DS of target 4 (Lead) and US of target 5 (Iron and Lead)
  {
    if (USorDS) //If true we're looking for US Planes therefore we're looking at US of tgt3
    {
      if (tgtUtil.InIron5VolMC( vtx_x, vtx_y, tgtUtil.GetTarget5CenterZMC(), 850., true )) return 5026;
      else if (tgtUtil.InLead5VolMC( vtx_x, vtx_y, tgtUtil.GetTarget5CenterZMC(), 850., true )) return 5082;
    }
    else //If false we're looking for DS Planes therefore we're looking at DS of tgt4
    { 
      if (tgtUtil.InLead4VolMC( vtx_x, vtx_y, tgtUtil.GetTarget4CenterZMC(), 850. )) return 4082;
    }
  } 
  if (mod >= 23 && mod <= 26) //DS of target 5 (Iron and Lead)
  {
    if (USorDS) //If true we're looking for US Planes therefore we're looking at US of tgt3
    {
      return -1;//This is immedaitely after target 5 and before the tracker region, so these planes are not the upstream region for any targets
    }
    else //If false we're looking for DS Planes therefore we're looking at DS of tgt4
    { 
      if (tgtUtil.InIron5VolMC( vtx_x, vtx_y, tgtUtil.GetTarget5CenterZMC(), 850., true )) return 5026;
      else if (tgtUtil.InLead5VolMC( vtx_x, vtx_y, tgtUtil.GetTarget5CenterZMC(), 850., true )) return 5082;
    }
  } 
  return -1;
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
    PlotUtils::Model<CVUniverse, MichelEvent>& model, 
    int targetCode)
{
  assert(!error_bands["cv"].empty() && "\"cv\" error band is empty!  Can't set Model weight.");
  auto& cvUniv = error_bands["cv"].front();

  std::cout << "Starting MC reco loop...\n";
  const int nEntries = chain->GetEntries();
  //std::cout<<"Here1\n";
  //const int nEntries = 10000;
  for (int i=0; i<nEntries; ++i)
  {
    if(i%1000==0) std::cout << i << " / " << nEntries << "\r" <<std::flush;
    //std::cout<<"Here2\n";
    MichelEvent cvEvent;
    cvUniv->SetEntry(i);
    model.SetEntry(*cvUniv, cvEvent);
    const double cvWeight = model.GetWeight(*cvUniv, cvEvent);
    //std::cout<<"Here3\n";
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
        //std::cout<<"Here5\n";
        universe->SetEntry(i);
        // This is where you would Access/create a Michel
        //weight is ignored in isMCSelected() for all but the CV Universe.
        //std::cout<<"Here6\n";
        if (!michelcuts.isMCSelected(*universe, myevent, cvWeight).all()) continue; //all is another function that will later help me with sidebands
        //std::cout<<"Here7\n";
        int code = getTgtCode(universe, 1);
        int truthcode = getTgtCode(universe, 0);
        //std::cout<<"Here8\n";
        int DSTargetCodeTruth = getSidebandTgtCode(universe, 0, 0); //Which, if any, target is this event truly in the DS Sideband region of?
        int USTargetCodeTruth = getSidebandTgtCode(universe, 0, 1); //Which, if any, target is this event truly in the US Sideband region of?
        //std::cout<<"Here9\n";
        int DSTargetCodeReco = getSidebandTgtCode(universe, 1, 0); //Which, if any, target is this event reconstructed in the DS Sideband region of?
        int USTargetCodeReco = getSidebandTgtCode(universe, 1, 1); //Which, if any, target is this event reconstructed in the US Sideband region of?
        //std::cout<<"Here10\n";
        //To do: use universe->hasMLPred()
        //Nuke Target Study
        const double weight = model.GetWeight(*universe, myevent); //Only calculate the per-universe weight for events that will actually use it.
        //Checking if events that are reconstructed outside of our targets of interest occur in our sideband region, which we are also interested in
        //std::cout<<"Here11\n";
        if (code!=targetCode && code<1000) //If our vertex is not in the nuclear targets
        {
          if (USTargetCodeReco==targetCode) //Get true origins of the events reconstructed in the upsteam region of tgt x
          {
            if (USTargetCodeTruth==USTargetCodeReco) //If the event truly occurred in the US region of tgt x
            {
              for(auto& var: vars) (*var->m_sidebandHistSetUSMC)[0].FillUniverse(universe, var->GetRecoValue(*universe), weight);
              for(auto& var: vars2D) (*var->m_sidebandHistSetUSMC)[0].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
            }
            else if (DSTargetCodeTruth==USTargetCodeReco) //If the event truly occurred in the DS region of tgt x
            {
              for(auto& var: vars) (*var->m_sidebandHistSetUSMC)[1].FillUniverse(universe, var->GetRecoValue(*universe), weight);
              for(auto& var: vars2D) (*var->m_sidebandHistSetUSMC)[1].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
            }
            else if (truthcode==USTargetCodeReco) //If the event truly occurred in the signal region (inside) of tgt x
            {
              for(auto& var: vars) (*var->m_sidebandHistSetUSMC)[2].FillUniverse(universe, var->GetRecoValue(*universe), weight);
              for(auto& var: vars2D) (*var->m_sidebandHistSetUSMC)[2].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
            }
            else 
            {
              for(auto& var: vars) (*var->m_sidebandHistSetUSMC)[-1].FillUniverse(universe, var->GetRecoValue(*universe), weight);
              for(auto& var: vars2D) (*var->m_sidebandHistSetUSMC)[-1].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
            }
          }
          //std::cout<<"Here12\n";
          if (DSTargetCodeReco==targetCode) //Get true origins of the events reconstructed in the downsteam region of tgt x
          {
            if (USTargetCodeTruth==DSTargetCodeReco) //If the event truly occurred in the US region of tgt x
            {
              for(auto& var: vars) (*var->m_sidebandHistSetDSMC)[0].FillUniverse(universe, var->GetRecoValue(*universe), weight);
              for(auto& var: vars2D) (*var->m_sidebandHistSetDSMC)[0].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
            }
            else if (DSTargetCodeTruth==DSTargetCodeReco) //If the event truly occurred in the DS region of tgt x
            {
              for(auto& var: vars) (*var->m_sidebandHistSetDSMC)[1].FillUniverse(universe, var->GetRecoValue(*universe), weight);
              for(auto& var: vars2D) (*var->m_sidebandHistSetDSMC)[1].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
            }
            else if (truthcode==DSTargetCodeReco) //If the event truly occurred in the signal region (inside) of tgt x
            {
              for(auto& var: vars) (*var->m_sidebandHistSetDSMC)[2].FillUniverse(universe, var->GetRecoValue(*universe), weight);
              for(auto& var: vars2D) (*var->m_sidebandHistSetDSMC)[2].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
            }
            else 
            {
              for(auto& var: vars) (*var->m_sidebandHistSetDSMC)[-1].FillUniverse(universe, var->GetRecoValue(*universe), weight);
              for(auto& var: vars2D) (*var->m_sidebandHistSetDSMC)[-1].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
            }
          }
        }

        //std::cout<<"Here13\n";
        if (code!=targetCode) continue; //From here on, we only care about events that are reconstructed in our target of interest
        //std::cout<<"Here14\n";
        for(auto& var: vars)
        {
          (*var->selectedMCReco).FillUniverse(universe, var->GetRecoValue(*universe), weight);
          (*var->m_interactionTypeHists)[universe->GetInteractionType()].FillUniverse(universe, var->GetRecoValue(*universe), weight);
        }
        //std::cout<<"Here15\n";
        for(auto& var: vars2D)
        {
          //std::cout<<"Here15a\n";
          (*var->selectedMCReco).FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe) , weight);
          //std::cout<<"Here15b\n";
          (*var->m_interactionTypeHists)[universe->GetInteractionType()].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
          //std::cout<<"Here15c\n";
        }
        //std::cout<<"Here16\n";
        const bool isSignal = michelcuts.isSignal(*universe, weight);
        if(isSignal && truthcode==targetCode) //If it is signal and is in nuclear target of interest
        {
          //std::cout<<"Here17\n";
          for(auto& study: studies) study->SelectedSignal(*universe, myevent, weight);
          //std::cout<<"Here18\n";
          for(auto& var: vars)
          {
            //Cross section components
            (*var->efficiencyNumerator).FillUniverse(universe, var->GetTrueValue(*universe), weight);
            (*var->migration).FillUniverse(universe, var->GetRecoValue(*universe), var->GetTrueValue(*universe), weight);
            (*var->selectedSignalReco).FillUniverse(universe, var->GetRecoValue(*universe), weight);
          }
          //std::cout<<"Here19\n";
          for(auto& var: vars2D)
          {
            //Cross section components
            (*var->efficiencyNumerator).FillUniverse(universe, var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
            (*var->migration).Fill(var->GetRecoValueX(*universe), var->GetRecoValueY(*universe),var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
            (*var->selectedSignalReco).FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
          }
        }
        //------------------------------------------------------
        //Backgrounds
        //------------------------------------------------------
        else if( !isSignal) //i.e if it passed the cuts but isn't a signal event/what we're looking for it must be a background
        {
          int bkgd_ID = -1;
          if (universe->GetCurrent()==2)bkgd_ID=0;
          else bkgd_ID=1;
          if (code>0) //If this event is reconstructed inside a nuclear target
          {
            for(auto& var: vars) (*var->m_backgroundHists)[bkgd_ID].FillUniverse(universe, var->GetRecoValue(*universe), weight);
            for(auto& var: vars2D) (*var->m_backgroundHists)[bkgd_ID].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
          }
        }
        else if ( code!=truthcode) //If the event truly came from somewhere other than where it was reconstructed
        {
          if (code == USTargetCodeTruth) //If it was reconstructed in tgt x but actually came from the US sideband region of tgt x
          {
            for(auto& var: vars) (*var->m_backgroundHists)[2].FillUniverse(universe, var->GetRecoValue(*universe), weight);
            for(auto& var: vars2D) (*var->m_backgroundHists)[2].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
          }
          else if (code == DSTargetCodeTruth) //If it was reconstructed in tgt x but actually came from the DS sideband region of tgt x
          {
            for(auto& var: vars) (*var->m_backgroundHists)[3].FillUniverse(universe, var->GetRecoValue(*universe), weight);
            for(auto& var: vars2D) (*var->m_backgroundHists)[3].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
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
				PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts, 
        int targetCode)
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

      int code = getTgtCode(universe, 1);
      int DSTargetCodeReco = getSidebandTgtCode(universe, 1, 0); //Which, if any, target is this event reconstructed in the DS Sideband region of?
      int USTargetCodeReco = getSidebandTgtCode(universe, 1, 1); //Which, if any, target is this event reconstructed in the US Sideband region of?

      for(auto& study: studies) study->Selected(*universe, myevent, 1); 
      //Checking if events that are reconstructed outside of our target of interest occur in our sideband region, which we are also interested in
      if (code!=targetCode && code<1000) //If our vertex is not in the nuclear targets
      {
        if (USTargetCodeReco==targetCode) //Get true origins of the events reconstructed in the upsteam region of tgt x
        {
          for(auto& var: vars) (*var->m_US_Sideband_Data).FillUniverse(universe, var->GetRecoValue(*universe), 1);
          for(auto& var: vars2D) (*var->m_US_Sideband_Data).FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), 1);
        }
        else if (DSTargetCodeReco==targetCode) //Get true origins of the events reconstructed in the upsteam region of tgt x
        {
          for(auto& var: vars) (*var->m_DS_Sideband_Data).FillUniverse(universe, var->GetRecoValue(*universe), 1);
          for(auto& var: vars2D) (*var->m_DS_Sideband_Data).FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), 1);
        }
      }

      if (code!=targetCode) continue; //From here on, we only care about events that are reconstructed in our target of interest

      for(auto& var: vars) (*var->dataHist).FillUniverse(universe, var->GetRecoValue(*universe, myevent.m_idx), 1);
      for(auto& var: vars2D) (*var->dataHist).FillUniverse(universe, var->GetRecoValueX(*universe, myevent.m_idx), var->GetRecoValueY(*universe, myevent.m_idx), 1);
    }
  }
  std::cout << "Finished data loop.\n";
}

void LoopAndFillEffDenom( PlotUtils::ChainWrapper* truth,
    std::map<std::string, std::vector<CVUniverse*> > truth_bands,
    std::vector<Variable1DNuke*> vars,
    std::vector<Variable2DNuke*> vars2D,
    PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts,
    PlotUtils::Model<CVUniverse, MichelEvent>& model, 
    int targetCode)
{
  assert(!truth_bands["cv"].empty() && "\"cv\" error band is empty!  Could not set Model entry.");
  auto& cvUniv = truth_bands["cv"].front();

  std::cout << "Starting efficiency denominator loop...\n";
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

        int code = getTgtCode(universe, 0);
        if (code!=targetCode) continue; //From here on, we only care about events that are really in our target of interest
        //Fill efficiency denominator now: 
        for(auto var: vars)
        {
          (*var->efficiencyDenominator).FillUniverse(universe, var->GetTrueValue(*universe), weight);
          (*var->m_intChannelsEffDenom)[universe->GetInteractionType()].FillUniverse(universe, var->GetTrueValue(*universe), weight);
        }
        for(auto var: vars2D)
        {
          (*var->m_intChannelsEffDenom)[universe->GetInteractionType()].FillUniverse(universe, var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
          (*var->efficiencyDenominator).FillUniverse(universe, var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
        }
      }
    }
  }
  std::cout << "Finished efficiency denominator loop.\n";
}

//Returns false if recoTreeName could not be inferred
bool inferRecoTreeNameAndCheckTreeNames(const std::string& mcPlaylistName, const std::string& dataPlaylistName, std::string& recoTreeName)
{
  std::cout<<"Here123.1\n";
  const std::vector<std::string> knownTreeNames = {"Truth", "Meta"};
  bool areFilesOK = false;
  std::ifstream playlist(mcPlaylistName);
  std::string firstFile = "";
  playlist >> firstFile;
  std::cout<<"Here123.2\n";
  std::cout<<"firstFile "<<firstFile<<"\n";
  auto testFile = TFile::Open(firstFile.c_str());
  if(!testFile)
  {
    std::cout << "Failed to open the first MC file at " << firstFile << "\n";
    return false;
  }

  //Does the MC playlist have the Truth tree?  This is needed for the efficiency denominator.
  const auto truthTree = testFile->Get("Truth");
  if(truthTree == nullptr || !truthTree->IsA()->InheritsFrom(TClass::GetClass("TTree")))
  {
    std::cout << "Could not find the \"Truth\" tree in MC file named " << firstFile << "\n";
    return false;
  }
  std::cout<<"Here123.3\n";
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
  std::cout<<"Here123.4\n";
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
  std::cout<<"Here123.5\n";
  const auto recoTree = testFile->Get(recoTreeName.c_str());
  if(recoTree == nullptr || !recoTree->IsA()->InheritsFrom(TClass::GetClass("TTree")))
  {
    std::cerr << "Could not find the \"" << recoTreeName << "\" tree in data file named " << firstFile << "\n";
    return false;
  }
  std::cout<<"Here123.6\n";

  return areFilesOK;
}

//==============================================================================
// Main
//==============================================================================
int main(const int argc, const char** argv)
{
  std::cout<<"Running event loop\n";
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
  const std::string mc_file_list = argv[2],
                    data_file_list = argv[1];
  std::vector<int> targets = {};
  std::cout<<"Here2\n";
  if (argc>3)
  {
    for (int i = 3; i<argc; i++)
    {
      int tgtToAdd = std::stoi(argv[i]);
      std::cout<<"Adding target " << tgtToAdd << std::endl;
      targets.push_back(tgtToAdd);
    }  
  }
  else //If no argument is given, do all targets
  {
    for (auto code: util::TgtCodeLabelsNuke)
    {
      int tgtToAdd = code.first;
      std::cout<<"Adding target " << tgtToAdd << std::endl;
      targets.push_back(tgtToAdd);
    }
  }
  std::cout<<"Here3\n";

  //Check that necessary TTrees exist in the first file of mc_file_list and data_file_list
  std::string reco_tree_name;
  if(!inferRecoTreeNameAndCheckTreeNames(mc_file_list, data_file_list, reco_tree_name))
  {
    std::cout << "Failed to find required trees in MC playlist " << mc_file_list << " and/or data playlist " << data_file_list << ".\n" << USAGE << "\n";
    return badInputFile;
  }
  std::cout<<"Here4\n";
  const bool doCCQENuValidation = (reco_tree_name == "CCQENu"); //Enables extra histograms and might influence which systematics I use.
  std:: cout << reco_tree_name << std::endl;

  std::string playlistname = "minervame1A";
  size_t mc_label = mc_file_list.find("/MC/");
  if (mc_label == std::string::npos) std::cout << "FAILED TO DETERMINE WHICH PLAYLIST THIS IS. DEFAULTING TO ME1A FOR FLUX\n";
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
  std::cout<<"playlistname: " << playlistname << std::endl; //isn't all this unnecessary?
  PlotUtils::MacroUtil options(reco_tree_name, mc_file_list, data_file_list, playlistname, true);
  options.m_plist_string = util::GetPlaylist(*options.m_mc, true); //TODO: Put GetPlaylist into PlotUtils::MacroUtil

  // You're required to make some decisions
  PlotUtils::MinervaUniverse::SetNuEConstraint(true);
  PlotUtils::MinervaUniverse::SetPlaylist(options.m_plist_string); //TODO: Infer this from the files somehow?
  int nuoranu = util::nuOrAntiNuMode(options.m_plist_string);
  int nupdg;
  if (nuoranu==1) nupdg = 14;
  else if (nuoranu==2) nupdg = -14;
  std::cout<<" nupdg: " <<  nupdg <<std::endl; 
  PlotUtils::MinervaUniverse::SetAnalysisNuPDG(nupdg);
  PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
  PlotUtils::MinervaUniverse::SetZExpansionFaReweight(false);

  PlotUtils::MinervaUniverse::RPAMaterials(true); 

  //Now that we've defined what a cross section is, decide which sample and model
  //we're extracting a cross section for.
  PlotUtils::Cutter<CVUniverse, MichelEvent>::reco_t nukeSidebands, nukePreCut;
  PlotUtils::Cutter<CVUniverse, MichelEvent>::truth_t nukeSignalDefinition, nukePhaseSpace;

  const double apothem = 850; //All in mm
  nukePreCut.emplace_back(new reco::ZRangeANN<CVUniverse, MichelEvent>("Nuclear Targets Z pos", PlotUtils::TargetProp::NukeRegion::Face, PlotUtils::TargetProp::NukeRegion::Back));
  nukePreCut.emplace_back(new reco::Apothem<CVUniverse, MichelEvent>(apothem));
  nukePreCut.emplace_back(new reco::MaxMuonAngle<CVUniverse, MichelEvent>(17.));
  nukePreCut.emplace_back(new reco::HasMINOSMatch<CVUniverse, MichelEvent>());
  nukePreCut.emplace_back(new reco::NoDeadtime<CVUniverse, MichelEvent>(1, "Deadtime"));
  if (nupdg>0)  nukePreCut.emplace_back(new reco::IsNeutrino<CVUniverse, MichelEvent>()); //Used minos curvature
  else if (nupdg<0)  nukePreCut.emplace_back(new reco::IsAntiNeutrino<CVUniverse, MichelEvent>()); //Used minos curvature
  nukePreCut.emplace_back(new reco::MuonCurveSignificance<CVUniverse, MichelEvent>(5));
  nukePreCut.emplace_back(new reco::MuonEnergyMin<CVUniverse, MichelEvent>(2000.0, "EMu Min"));
  nukePreCut.emplace_back(new reco::MuonEnergyMax<CVUniverse, MichelEvent>(50000.0, "EMu Max"));
  nukePreCut.emplace_back(new reco::ANNConfidenceCut<CVUniverse, MichelEvent>(0.20));

  //nukeSidebands.emplace_back(new reco::ZRange<CVUniverse, MichelEvent>("Test sideband z pos", 0, 1000000000000.0));
  //nukeSidebands.emplace_back(new reco::USScintillator<CVUniverse, MichelEvent>());
  //nukeSidebands.emplace_back(new reco::DSScintillator<CVUniverse, MichelEvent>());

                                                                                                                   
  if (nupdg>0) nukeSignalDefinition.emplace_back(new truth::IsNeutrino<CVUniverse>());
  else if (nupdg<0) nukeSignalDefinition.emplace_back(new truth::IsAntiNeutrino<CVUniverse>());
  nukeSignalDefinition.emplace_back(new truth::IsCC<CVUniverse>());
                                                        
  nukePhaseSpace.emplace_back(new truth::ZRange<CVUniverse>("Nuclear Targets Z pos", PlotUtils::TargetProp::NukeRegion::Face, PlotUtils::TargetProp::NukeRegion::Back));
  nukePhaseSpace.emplace_back(new truth::Apothem<CVUniverse>(apothem));
  nukePhaseSpace.emplace_back(new truth::MuonAngle<CVUniverse>(17.));
  nukePhaseSpace.emplace_back(new truth::MuonEnergyMin<CVUniverse>(2000.0, "EMu Min"));
  nukePhaseSpace.emplace_back(new truth::MuonEnergyMax<CVUniverse>(50000.0, "EMu Max"));
  //nukePhaseSpace.emplace_back(new truth::PZMuMin<CVUniverse>(1500.));
                                                                                                                                                   
  PlotUtils::Cutter<CVUniverse, MichelEvent> nukeCuts(std::move(nukePreCut), std::move(nukeSidebands) , std::move(nukeSignalDefinition),std::move(nukePhaseSpace));

  const bool NO_2P2H_WARP = (getenv("NO_2P2H_WARP") != nullptr);
  if(NO_2P2H_WARP){
    std::cout << "Turning off LowRecoil2p2hReweighter because environment variable NO_2P2H_WARP is set.\n";
  }
  const bool AMU_DIS_WARP = (getenv("AMU_DIS_WARP") != nullptr);
  if(AMU_DIS_WARP){
  std::cout << "Turning on AMUDISReweighter because environment variable AMU_DIS_WARP is set.\n";
  }
  const bool LOW_Q2_PION_WARP = (getenv("LOW_Q2_PION_WARP") != nullptr);
  if(LOW_Q2_PION_WARP){
  std::cout << "Turning on LowQ2PiReweighter because environment variable LOW_Q2_PION_WARP is set.\n";
  }


  std::vector<std::unique_ptr<PlotUtils::Reweighter<CVUniverse, MichelEvent>>> MnvTunev1;
  MnvTunev1.emplace_back(new PlotUtils::FluxAndCVReweighter<CVUniverse, MichelEvent>());
  MnvTunev1.emplace_back(new PlotUtils::GENIEReweighter<CVUniverse, MichelEvent>(true, false));
  // standard - include this tune.  warping study, comment out 2p2h
  if ( !NO_2P2H_WARP ) MnvTunev1.emplace_back(new PlotUtils::LowRecoil2p2hReweighter<CVUniverse, MichelEvent>());
  MnvTunev1.emplace_back(new PlotUtils::MINOSEfficiencyReweighter<CVUniverse, MichelEvent>());
  MnvTunev1.emplace_back(new PlotUtils::RPAReweighter<CVUniverse, MichelEvent>());
  // for a warping study, AMU DIS reweighter
  if ( AMU_DIS_WARP ) MnvTunev1.emplace_back(new PlotUtils::AMUDISReweighter<CVUniverse, MichelEvent>());
  // for a warping study, Low Q2 pion suppression (mnvtunev2)
  if ( LOW_Q2_PION_WARP ) MnvTunev1.emplace_back(new PlotUtils::LowQ2PiReweighter<CVUniverse, MichelEvent>("JOINT"));

  PlotUtils::Model<CVUniverse, MichelEvent> model(std::move(MnvTunev1));
  // Make a map of systematic universes
  // Leave out systematics when making validation histograms
  const bool doSystematics = (getenv("MNV101_SKIP_SYST") == nullptr);
  if(!doSystematics){
    std::cout << "Skipping systematics (except 1 flux universe) because environment variable MNV101_SKIP_SYST is set.\n";
    PlotUtils::MinervaUniverse::SetNFluxUniverses(2); //Necessary to get Flux integral later...  Doesn't work with just 1 flux universe though because _that_ triggers "spread errors".
  }

  std::map< std::string, std::vector<CVUniverse*> > error_bands;
  if(doSystematics)  error_bands = GetStandardSystematics(options.m_mc);
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
                      robsRecoilBins, segmentBins;

  const double robsRecoilBinWidth = 50; //MeV
  for(int whichBin = 0; whichBin < 100 + 1; ++whichBin) robsRecoilBins.push_back(robsRecoilBinWidth * whichBin);

  const double numsegments = 180;
  for(double whichBin = 0; whichBin < numsegments; whichBin++) segmentBins.push_back(whichBin - 0.5);

  std::vector<Variable1DNuke*> nukeVars;
  std::vector<Variable2DNuke*> nukeVars2D;

  nukeVars.push_back(new Variable1DNuke("nuke_pTmu", "p_{T, #mu} [GeV/c]", dansPTBins, &CVUniverse::GetMuonPT, &CVUniverse::GetMuonPTTrue));
  nukeVars.push_back(new Variable1DNuke("nuke_pZmu", "p_{||, #mu} [GeV/c]", dansPzBins, &CVUniverse::GetMuonPz, &CVUniverse::GetMuonPzTrue));
  nukeVars.push_back(new Variable1DNuke("nuke_Emu", "E_{#mu} [GeV]", robsEmuBins, &CVUniverse::GetEmuGeV, &CVUniverse::GetElepTrueGeV));
  nukeVars.push_back(new Variable1DNuke("nuke_Erecoil", "E_{recoil}", robsRecoilBins, &CVUniverse::GetRecoilE, &CVUniverse::Getq0True)); //TODO: q0 is not the same as recoil energy without a spline correction
  nukeVars.push_back(new Variable1DNuke("nuke_BjorkenX", "X", bjorkenXbins, &CVUniverse::GetBjorkenX, &CVUniverse::GetBjorkenXTrue));
  nukeVars.push_back(new Variable1DNuke("nuke_segment", "segmentNum", segmentBins, &CVUniverse::GetANNSegment, &CVUniverse::GetTruthSegment)); //Just used for plotting events by detector position tbh - not for any actual physics
  nukeVars2D.push_back(new Variable2DNuke("nuke_pTmu_pZmu", *nukeVars[1], *nukeVars[0]));

  std::vector<Study*> studies;
  std::function<double(const CVUniverse&, const MichelEvent&)> ptmu = [](const CVUniverse& univ, const MichelEvent& /* evt */) { return univ.GetMuonPT();};
  std::function<double(const CVUniverse&, const MichelEvent&)> pzmu = [](const CVUniverse& univ, const MichelEvent& /* evt */) { return univ.GetMuonPz();};

  //studies.push_back(new PerEventVarByGENIELabel2D(pzmu, ptmu, std::string("pzmu_vs_ptmu_GENIE_labels"), std::string("GeV/c"), dansPzBins, dansPTBins, error_bands));
  //studies.push_back(new WaterTargetIntOrigin2D(pzmu, ptmu, std::string("pzmu_vs_ptmu_water_breakdown"), std::string("GeV/c"), dansPzBins, dansPTBins, error_bands));

  CVUniverse* data_universe = new CVUniverse(options.m_data);
  std::vector<CVUniverse*> data_band = {data_universe};
  std::map<std::string, std::vector<CVUniverse*> > data_error_bands;
  data_error_bands["cv"] = data_band;
  
  std::vector<Study*> data_studies;
  //data_studies.push_back(new PerEventVarByGENIELabel2D(ptmu, pzmu, std::string("ptmu_vs_pzmu"), std::string("GeV/c"), dansPTBins, dansPzBins, data_error_bands));
  //Wouldn't make sense to do a PerEventVarByGENIELabel2D study for data since data wont have the GENIE simulation labels

  for (auto tgt : targets)
  {
    for(auto& var: nukeVars) var->InitializeMCHists(error_bands, truth_bands);
    for(auto& var: nukeVars) var->InitializeDATAHists(data_band);

    for(auto& var: nukeVars2D) var->InitializeMCHists(error_bands, truth_bands);
    for(auto& var: nukeVars2D) var->InitializeDATAHists(data_band);

    // Loop entries and fill
    try
    {
      std::cout << "Staring event loops\n";
      CVUniverse::SetTruth(false);
      LoopAndFillEventSelection(options.m_mc, error_bands, nukeVars, nukeVars2D, studies, nukeCuts, model, tgt);
      CVUniverse::SetTruth(true);
      LoopAndFillEffDenom(options.m_truth, truth_bands, nukeVars, nukeVars2D, nukeCuts, model, tgt);
      options.PrintMacroConfiguration(argv[0]);
      std::cout << "Nuclear Target MC cut summary:\n" << nukeCuts << "\n";
      nukeCuts.resetStats();

      CVUniverse::SetTruth(false);
      LoopAndFillData(options.m_data, data_band, nukeVars, nukeVars2D, data_studies, nukeCuts, tgt);
      std::cout << "Nuclear Target Data cut summary:\n" << nukeCuts << "\n";


      std::string mcOutFileName = MC_OUT_FILE_NAME_BASE + std::to_string(tgt)+".root";
      //Write MC results
      TFile* mcOutDir = TFile::Open(mcOutFileName.c_str(), "RECREATE");
      if(!mcOutDir)
      {
        std::cerr << "Failed to open a file named " << mcOutFileName << " in the current directory for writing histograms.\n";
        return badOutputFile;
      }
      std::cout<<"Saving " << studies.size()<< " studies\n";
      for(auto& study: studies) study->SaveOrDraw(*mcOutDir);
      std::cout<<"Saved studies\n";

      for(auto& var: nukeVars) var->WriteMC(*mcOutDir);
      std::cout<<"Saved 1D Variables\n";
      for(auto& var: nukeVars2D) var->WriteMC(*mcOutDir);
      std::cout<<"Saved 2D Variables\n";


      //Protons On Target
      auto mcPOT = new TParameter<double>("POTUsed", options.m_mc_pot);
      mcPOT->Write();

      PlotUtils::TargetUtils targetInfo;
      assert(error_bands["cv"].size() == 1 && "List of error bands must contain a universe named \"cv\" for the flux integral.");

      for(auto& var: nukeVars) 
      {
        //Flux integral only if systematics are being done (temporary solution)
        util::GetFluxIntegral(*error_bands["cv"].front(), var->efficiencyNumerator->hist)->Write((var->GetName() + "_reweightedflux_integrated").c_str());
        //Always use MC number of nucleons for cross section
        //This may not even be necessary since we can always pull the same information in the extract cross section script as long ad we have the target information, which we do
        TParameter<double>* nNucleons;
        if (tgt==6000) nNucleons = new TParameter<double>((var->GetName() +"_fiducial_nucleons").c_str(), targetInfo.GetPassiveTargetNNucleons( 6, 1, true));
        else if (tgt>=7 && tgt<=11)
        {
          nNucleons = new TParameter<double>((var->GetName() + "_fiducial_nucleons").c_str(), targetInfo.GetTrackerNNucleons( 6, true));
        }
        else if (tgt == 12) new TParameter<double>((var->GetName() + "_fiducial_nucleons").c_str(), targetInfo.GetTrackerNNucleons( 2, true));
        else
        {
          int tgtZ = tgt%1000;
          int tgtID = (tgt-tgtZ)/1000;
          nNucleons = new TParameter<double>((var->GetName() +"_fiducial_nucleons").c_str(), targetInfo.GetPassiveTargetNNucleons( tgtID, tgtZ, true));
        }
        nNucleons->Write();
      }

      //Write data results
      std::string dataOutFileName = DATA_OUT_FILE_NAME_BASE + std::to_string(tgt)+".root";
      TFile* dataOutDir = TFile::Open(dataOutFileName.c_str(), "RECREATE");
      if(!dataOutDir)
      {
        std::cerr << "Failed to open a file named " << dataOutFileName << " in the current directory for writing histograms.\n";
        return badOutputFile;
      }

      for(auto& var: nukeVars) var->WriteData(*dataOutDir);
      for(auto& var: nukeVars2D) var->WriteData(*dataOutDir);

      for(auto& study: data_studies) study->SaveOrDraw(*dataOutDir);

      //Protons On Target
      auto dataPOT = new TParameter<double>("POTUsed", options.m_data_pot);
      dataPOT->Write();

      //Saving 2D migration matrices
      //Putting this right at the end in case of a crash
      std::string migrationOutDirName = MIGRATION_2D_OUT_FILE_NAME_BASE + std::to_string(tgt)+".root";
      TFile* migrationOutDir = TFile::Open(migrationOutDirName.c_str(), "RECREATE");
      if(!dataOutDir)
      {
        std::cerr << "Failed to open a file named " << migrationOutDirName << " in the current directory for writing histograms.\n";
        return badOutputFile;
      }
      for(auto& var: nukeVars2D) 
      {
        var->WriteMigration(*migrationOutDir); //Save 2D migration to separate files, because it's huge
      } 

      std::cout << "Success" << std::endl;
    }
    catch(const ROOT::exception& e)
    {
      std::cerr << "Ending on a ROOT error message.  No histograms will be produced.\n"
                << "If the message talks about \"TNetXNGFile\", this could be a problem with dCache.  The message is:\n"
                << e.what() << "\n" << USAGE << "\n";
      return badFileRead;
    }
  }
  return success;
}