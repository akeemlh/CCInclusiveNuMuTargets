#ifndef Systematics_h
#define Systematics_h

//==============================================================================
// Get Several standard MINERvA systematics
//==============================================================================

#include "event/CVUniverse.h"
#include "PlotUtils/FluxSystematics.h"
#include "PlotUtils/GenieSystematics.h"
#include "PlotUtils/MinosEfficiencySystematics.h"
#include "PlotUtils/MuonSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "PlotUtils/MuonResolutionSystematics.h"
#include "PlotUtils/AngleSystematics.h"
#include "PlotUtils/TargetMassSystematics.h"
#include "PlotUtils/MLVertexSystematics.h"
#include "PlotUtils/CoherentPiSystematics.h"
#include "PlotUtils/DiffractiveUncSystematics.h"
#include "universes/ResponseSystematics.h"
#include "PlotUtils/MLClassificationSystematics.h"

#include "weighters/Reweighter.h"

typedef std::map<std::string, std::vector<CVUniverse*>> UniverseMap;

UniverseMap GetStandardSystematics(PlotUtils::ChainWrapper* chain)
{
  // return map
  UniverseMap error_bands;

  // CV
  error_bands[std::string("cv")].push_back(new CVUniverse(chain));

  //========================================================================
  // FLUX
  //========================================================================
  UniverseMap bands_flux =
      PlotUtils::GetFluxSystematicsMap<CVUniverse>(chain, CVUniverse::GetNFluxUniverses());
  error_bands.insert(bands_flux.begin(), bands_flux.end());

  //========================================================================
  // GENIE
  //========================================================================
  // Standard
  UniverseMap bands_genie =
      PlotUtils::GetGenieSystematicsMap<CVUniverse>(chain); //PlotUtils::GetStandardGenieSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_genie.begin(), bands_genie.end());

  //========================================================================
  // Muons
  //========================================================================
  // Muon reco in MINERvA -- Catchall systematic for pmu reco in minerva.
  // Lateral-only. Shifts pmu.
  UniverseMap bands_muon_minerva =
      PlotUtils::GetMinervaMuonSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_muon_minerva.begin(), bands_muon_minerva.end());

  // Muons in MINOS -- Catchall systematic for wiggle solution -- correlates
  // flux universes and minos muon momentum reco.
  // Lateral AND Vertical systematic. Shifts Pmu and GetFluxAndCVUniverse.
  //
  // Expect a non-zero systematic even when no pmu involved.
  UniverseMap bands_muon_minos =
     PlotUtils::GetMinosMuonSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_muon_minos.begin(), bands_muon_minos.end());

  UniverseMap bands_muon_resolution = PlotUtils::GetMuonResolutionSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_muon_resolution.begin(), bands_muon_resolution.end());

  //Muon Angle Resolution
  UniverseMap muonAngle_resolutions = PlotUtils::GetMuonAngleResolutionSystematicsMap<CVUniverse>(chain);
  error_bands.insert(muonAngle_resolutions.begin(),muonAngle_resolutions.end());

  UniverseMap bands_geant = PlotUtils::GetGeantHadronSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_geant.begin(), bands_geant.end());

  // Beam angle
  UniverseMap bands_angle = PlotUtils::GetAngleSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_angle.begin(), bands_angle.end());

  // Target Mass Systematics
  UniverseMap targetMass_systematics = PlotUtils::GetTargetMassSystematicsMap<CVUniverse>(chain);
  error_bands.insert(targetMass_systematics.begin(), targetMass_systematics.end());

  // Hadron inelastics cross sections
  //TODO: There's some special recoil function I need to write for the response systematics to work correctly
  /*UniverseMap bands_response = PlotUtils::GetResponseSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_response.begin(), bands_response.end());*/

  //========================================================================
  // Particle Response Systematics (Detector)
  //========================================================================

  const bool use_ID = true;
  const bool use_OD = true;
  std::string name_tag = "allNonMuonClusters";
  const bool use_neutron = false;
  const bool use_new = false;
  const bool use_proton = true;
  //const bool use_nucl = false;
  //const bool use_tracker = false;
  //const bool use_ecal = false;
  //const bool use_hcal = false;
  // ID, OD, nucl, tracker, ecal, hcal
  UniverseMap response_systematics = PlotUtils::GetResponseSystematicsMap<CVUniverse>(chain, use_ID, use_OD, name_tag, use_neutron, use_new, use_proton); //, use_ID, use_OD, use_nucl, use_tracker, use_ecal, use_hcal);

  error_bands.insert(response_systematics.begin(), response_systematics.end());

  //========================================================================
  // ML Vertex Systematics
  //========================================================================
  //Efficiency
  UniverseMap ml_vertex_systematics = PlotUtils::GetMLVertexSystematicsMap<CVUniverse>(chain);
  error_bands.insert(ml_vertex_systematics.begin(), ml_vertex_systematics.end());
 
  //Classification
  //UniverseMap ml_vertex_classification_systematics = PlotUtils::GetMLVertexClassificationSystematicsMap<CVUniverse>(chain);
  //error_bands.insert(ml_vertex_classification_systematics.begin(), ml_vertex_classification_systematics.end());

  return error_bands;
}

using Tune = std::unique_ptr<PlotUtils::Reweighter<CVUniverse, MichelEvent>>;

UniverseMap GetModelSystematics(PlotUtils::ChainWrapper* chain, std::vector<Tune>& model)
{
  // return map
  UniverseMap error_bands;
  for (Tune& t : model)
  {
    if (t->GetName() == "RPA")
    {
        // RPA
        UniverseMap bands_rpa = PlotUtils::GetRPASystematicsMap<CVUniverse>(chain);
        error_bands.insert(bands_rpa.begin(), bands_rpa.end());
        std::cout<<"Applying RPA reweighter systematics\n";
    }
    else if (t->GetName() == "LowRecoil2p2hTune")
    {
        // 2P2H
        UniverseMap bands_2p2h = PlotUtils::Get2p2hSystematicsMap<CVUniverse>(chain);
        error_bands.insert(bands_2p2h.begin(), bands_2p2h.end());
        std::cout<<"Applying 2P2H reweighter systematics\n";
    }
    else if (t->GetName() == "LowQ2Pi_MENU1PI") //If using tune x.3.z which uses the A. Bercellie low Q2 pion production suppression
    {
        UniverseMap lowq2_systematics = PlotUtils::GetLowQ2PiSystematicsMap<CVUniverse>(chain);
        error_bands.insert(lowq2_systematics.begin(), lowq2_systematics.end());
        std::cout<<"Applying A. Bercellie low Q2 pion production suppression reweighter systematics\n";
    }
    else if (t->GetName() == "DiffractiveReweighter") //If using tune x.3.z which uses the A. Bercellie low Q2 pion production suppression
    {
        // Diffractive uncertainty
        UniverseMap diff_systematics = PlotUtils::GetDiffractiveUncSystematicsMap<CVUniverse>(chain);
        error_bands.insert(diff_systematics.begin(), diff_systematics.end());
        std::cout<<"Applying diffractive reweighter systematics\n";
    }
    else if (t->GetName() == "COHPionReweighter") //If using tune x.3.z which uses the A. Bercellie low Q2 pion production suppression
    {
        // Coherent uncertainty
        UniverseMap coh_systematics = PlotUtils::GetCoherentPiSystematicsMap<CVUniverse>(chain);
        error_bands.insert(coh_systematics.begin(), coh_systematics.end());
        std::cout<<"Applying coherent reweighter systematics\n";
    }
    else if (t->GetName() == "MINOSEfficiency")
    {
        // MINOS Efficiency
        UniverseMap bands_minoseff = PlotUtils::GetMinosEfficiencySystematicsMap<CVUniverse>(chain);
        error_bands.insert(bands_minoseff.begin(), bands_minoseff.end());
        std::cout<<"Applying MINOS Efficiency reweighter systematics\n";
    }
    else if (t->GetName() == "LowRecPionReweight")
    {
        // Q2Tpi Efficiency
        UniverseMap bands_UntrackedPion = PlotUtils::GetUntrackedPionSystematicsMap<CVUniverse>(chain);
        error_bands.insert(bands_UntrackedPion.begin(), bands_UntrackedPion.end());
        std::cout<<"Applying Untracked Pion reweighter systematics\n";
    }
    else if (t->GetName() == "Q2TpiReweighterReweighter")
    {
        // Q2Tpi Efficiency
        UniverseMap bands_Q2Tpi = PlotUtils::GetChargedPionTuneSystematicMap<CVUniverse>(chain);
        error_bands.insert(bands_Q2Tpi.begin(), bands_Q2Tpi.end());
        std::cout<<"Applying Q2-Tpi reweighter systematics\n";
    }
    //Other model reweighters/tunes not handled here:
    //      FSIReweighter - Used in mnvTune vX.Y.1 - Already handled by GetGenieSystematicsMap
    //      GENIEReweighter - Used in all minerva tunes to date - already handled by GetGenieSystematicsMap
    //      FluxAndCVReweighter - Used in all minerva tunes to date - already handled by GetFluxSystematicsMap
    //      BodekRitchieReweighter - Used in mnvTune vX.3.Z - I dont think is actually handled anywhere? Will ignore anyway since it has no bearing on my analysis
    //      SuSAFromValencia2p2hReweighter - Used in mnvTune vX.3.Z - I dont think is actually handled anywhere? Will ignore anyway since it has no bearing on my analysis
  }
  return error_bands;
}

#endif  // Systematics_h
