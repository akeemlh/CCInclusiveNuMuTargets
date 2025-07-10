// =============================================================================
// Base class for an un-systematically shifted (i.e. CV) universe. Implement
// "Get" functions for all the quantities that you need for your analysis.
//
// This class inherits from PU::MinervaUniverse, which in turn inherits from
// PU::BaseUniverse. PU::BU defines the interface with anatuples.
// 
// Within the class, "WeightFunctions" and "MuonFunctions" are included to gain
// access to standardized weight and muon variable getters. See:
// https://cdcvs.fnal.gov/redmine/projects/minerva-sw/wiki/MinervaUniverse_Structure_
// for a full list of standardized functions you can use. In general, if a
// standard version of a function is available, you should be using it.
// =============================================================================
#ifndef CVUNIVERSE_H
#define CVUNIVERSE_H

#include <iostream>

#include "PlotUtils/MinervaUniverse.h"
#include "utilities/PhysicsVariables.h"
#include "Math/Vector3D.h"

class CVUniverse : public PlotUtils::MinervaUniverse {

  public:
  #include "PlotUtils/MuonFunctions.h" // GetMinosEfficiencyWeight
  #include "PlotUtils/TruthFunctions.h" //Getq3True
  // ========================================================================
  // Constructor/Destructor
  // ========================================================================
  CVUniverse(PlotUtils::ChainWrapper* chw, double nsigma = 0)
      : PlotUtils::MinervaUniverse(chw, nsigma) {}

  virtual ~CVUniverse() {}

  // ========================================================================
  // Quantities defined here as constants for the sake of below. Definition
  // matched to Dan's CCQENuInclusiveME variables from:
  // `/minerva/app/users/drut1186/cmtuser/Minerva_v22r1p1_OrigCCQENuInc/Ana/CCQENu/ana_common/include/CCQENuUtils.h`
  // ========================================================================
  static constexpr double M_n = 939.56536;
  static constexpr double M_p = 938.272013;
  static constexpr double M_nucleon = (1.5*M_n+M_p)/2.5;

  static constexpr int PDG_n = 2112;
  static constexpr int PDG_p = 2212;

  // ========================================================================
  // Write a "Get" function for all quantities access by your analysis.
  // For composite quantities (e.g. Enu) use a calculator function.
  //
  // In order to properly calculate muon variables and systematics use the
  // various functions defined in MinervaUniverse.
  // E.g. GetPmu, GetEmu, etc.
  // ========================================================================

  // Quantities only needed for cuts
  // Although unlikely, in principle these quanties could be shifted by a
  // systematic. And when they are, they'll only be shifted correctly if we
  // write these accessor functions.
  
  //Muon kinematics
  double GetMuonPT() const //GeV/c
  {
    return GetPmu()/1000. * sin(GetThetamu());
  }

  double GetMuonPz() const //GeV/c
  {
    return GetPmu()/1000. * cos(GetThetamu());
  }

  double GetMuonPTTrue() const //GeV/c
  {
    return GetPlepTrue()/1000. * sin(GetThetalepTrue());
  }

  double GetMuonPzTrue() const //GeV/c
  {
    return GetPlepTrue()/1000. * cos(GetThetalepTrue());
  }

  double GetEmuGeV() const //GeV
  {
    return GetEmu()/1000.;
  }

  double GetElepTrueGeV() const //GeV
  {
    return GetElepTrue()/1000.;
  }

  int GetInteractionType() const {
    return GetInt("mc_intType");
  }

  int GetTargetNucleon() const {
    return GetInt("mc_targetNucleon");
  }
  
  double GetBjorkenXTrue() const {
    return GetDouble("mc_Bjorkenx");
  }

  double GetBjorkenYTrue() const {
    return GetDouble("mc_Bjorkeny");
  }

  virtual bool IsMinosMatchMuon() const {
    return GetInt("has_interaction_vertex") == 1;
  }
  
  ROOT::Math::XYZTVector GetVertex() const
  {
    ROOT::Math::XYZTVector result;
    result.SetCoordinates(GetVec<double>("vtx").data());
    return result;
  }

  ROOT::Math::XYZTVector GetTrueVertex() const
  {
    ROOT::Math::XYZTVector result;
    result.SetCoordinates(GetVec<double>("mc_vtx").data());
    return result;
  }

  virtual int GetTDead() const {
    return GetInt("phys_n_dead_discr_pair_upstream_prim_track_proj");
  }
  
  //TODO: If there was a spline correcting Eavail, it might not really be Eavail.
  //      Our energy correction spline, one of at least 2 I know of, corrects q0
  //      so that we get the right neutrino energy in an inclusive sample.  So,
  //      this function could be correcting for neutron energy which Eavail should
  //      not do.
  virtual double GetEavail() const {
    return GetDouble("recoilE_SplineCorrected");
  }
  
  virtual double GetQ2Reco() const{
    return GetDouble("qsquared_recoil");
  }

  //GetRecoilE is designed to match the NSF validation suite
  virtual double GetRecoilE() const {
    return GetVecElem("recoil_summed_energy", 0);
  }
  
  virtual double Getq3() const{
    double eavail = GetEavail()/pow(10,3);
    double q2 = GetQ2Reco() / pow(10,6);
    double q3mec = sqrt(eavail*eavail + q2);
    return q3mec;
  }
   
  virtual int GetCurrent() const { return GetInt("mc_current"); }

  virtual int GetTruthNuPDG() const { return GetInt("mc_incoming"); }

  virtual double GetMuonQP() const {
    return GetDouble((GetAnaToolName() + "_minos_trk_qp").c_str());
  }

  virtual double GetMuonQPErr() const { //Relative error on charge over momentum for minos track
    return GetDouble((GetAnaToolName() + "_minos_trk_eqp_qp").c_str());
  }

  //Some functions to match CCQENuInclusive treatment of DIS weighting. Name matches same Dan area as before.
  virtual double GetTrueExperimentersQ2() const {
    double Enu = GetEnuTrue(); //MeV
    double Emu = GetElepTrue(); //MeV
    double thetaMu = GetThetalepTrue();
    return 4.0*Enu*Emu*pow(sin(thetaMu/2.0),2.0);//MeV^2
  }

  virtual double CalcTrueExperimentersQ2(double Enu, double Emu, double thetaMu) const{
    return 4.0*Enu*Emu*pow(sin(thetaMu/2.0),2.0);//MeV^2
  }

  virtual double GetTrueExperimentersW() const {
    double nuclMass = M_nucleon;
    int struckNucl = GetTargetNucleon();
    if (struckNucl == PDG_n){
      nuclMass=M_n;
    }
    else if (struckNucl == PDG_p){
      nuclMass=M_p;
    }
    double Enu = GetEnuTrue();
    double Emu = GetElepTrue();
    double thetaMu = GetThetalepTrue();
    double Q2 = CalcTrueExperimentersQ2(Enu, Emu, thetaMu);
    return TMath::Sqrt(pow(nuclMass,2) + 2.0*(Enu-Emu)*nuclMass - Q2);
  }

  virtual double GetInclQ2Reco() const{ 
    double branchval = GetDouble("MasterAnaDev_Q2_Inclusive");
    double retval = (branchval==-9999) ? 0 : branchval;
    return retval;
  }

  virtual double GetANNRecoilE() const {
    double branchval =  GetDouble("MasterAnaDev_ANN_recoil_E");
    double retval = (branchval==-1) ? 0 : branchval;
    return retval;
  }

  double GetBjorkenX() const //GeV
  {
    double q2 = GetInclQ2Reco();
    double M = (M_p + M_n) /2.0; //Not exact. To do - weight it based on target material p n ratio
    double r = GetANNRecoilE();
    double xbj = xBjorken(q2, M, r);
    return xbj;
  }

  double GetANNProb() const { return GetVecElem("ANN_plane_probs", 0); }
  double GetTruthMuE() const { return GetDouble("truth_muon_E") ; }

  int GetTargetZ() const {
    return GetInt("MasterAnaDev_targetZ");
  }

  int GetANNTargetZ() const { return GetInt("MasterAnaDev_ANN_targetZ"); }
  int GetTruthTargetZ() const { return GetInt("truth_targetZ");}
  int GetMCTargetZ() const { return GetInt("mc_targetZ");}

  int GetANNTargetID() const {return GetInt("MasterAnaDev_ANN_targetID");}
  int GetTruthTargetID() const {return GetInt("truth_targetID");}

  int GetANNTargetCode() const {return GetInt("MasterAnaDev_ANN_target_code");}
  int GetTruthTargetCode() const {return GetInt("truth_target_code");}

  int GetANNVtxPlane() const {return GetVecElemInt("ANN_vtx_planes", 0);}
  int GetTruthVtxPlane() const {return GetInt("truth_vtx_plane");}
  int GetMADVtxPlane() const {return GetInt("MasterAnaDev_vtx_plane");}

  int GetANNVtxModule() const {return GetVecElemInt("ANN_vtx_modules", 0);}
  int GetTruthVtxModule() const {return GetInt("truth_vtx_module");}
  int GetMADVtxModule() const {return GetInt("MasterAnaDev_vtx_module");}
  
  ROOT::Math::XYZVector GetANNVertex() const //Dangerous since ANN_vtx doesnt always have a 3-vector in it, calling SetCoordinates without a 3 vector gives a segfault
  {
    ROOT::Math::XYZVector result;
    //std::cout<<"GetInt(\"ANN_vtx_sz\"): " << GetInt("ANN_vtx_sz") <<std::endl;
    if (GetInt("ANN_vtx_sz") != 3) return ROOT::Math::XYZVector(-1.0,-1.0,-1.0); //To prevent segfaults when accessing events for which there is no ANN_vtx 3-vector
    result.SetCoordinates(GetVec<double>("ANN_vtx").data());
    return result;
  }

  std::vector<double> GetANNVertexVector() const { return GetVecDouble("ANN_vtx");}

  double GetANNSegment() const
  {
    return GetVecElemInt("ANN_segments", 0);
  }

 //Taken from Oscar's code (/exp/minerva/app/users/omorenop/cmtuser/git-Mat/Personal/Test/InclusiveUtils.h)
  double GetTruthSegment() const //Dangerous since ANN_vtx doesnt always have a 3-vector in it, calling SetCoordinates without a 3 vector gives a segfault
  {
    int mdl = GetTruthVtxModule();
    int plane = GetTruthVtxPlane();
    if (mdl == -5) {
        if (plane == 1) {return 1;}
        if (plane == 2) {return 2;}}
    if (mdl == -4) {
        if (plane == 1) {return 3;}
        if (plane == 2) {return 4;}}
    if (mdl == -3) {
        if (plane == 1) {return 5;}
        if (plane == 2) {return 6;}}
    if (mdl == -2) {
        if (plane == 1) {return 7;}
        if (plane == 2) {return 8;}}
    if (mdl == -1) {
        if (plane == 1) {return 9;} //target 1
        if (plane == 2) {return 9;}}
    if (mdl == 0) {
        if (plane == 1) {return 10;}
        if (plane == 2) {return 11;}}
    if (mdl == 1) {
        if (plane == 1) {return 12;}
        if (plane == 2) {return 13;}}
    if (mdl == 2) {
        if (plane == 1) {return 14;}
        if (plane == 2) {return 15;}}
    if (mdl == 3) {
        if (plane == 1) {return 16;}
        if (plane == 2) {return 17;}}
    if (mdl == 4) {
        if (plane == 1) {return 18;}} //Target 2!}
    if (mdl == 5) {
        if (plane == 1) {return 19;}
        if (plane == 2) {return 20;}}
    if (mdl == 6) {
        if (plane == 1) {return 21;}
        if (plane == 2) {return 22;}
    }
    if (mdl == 7) {
        if (plane == 1) {return 23;}
        if (plane == 2) {return 24;}}
    if (mdl == 8) {
        if (plane == 1) {return 25;}
        if (plane == 2) {return 26;}}
    if (mdl == 9) {
        if (plane == 1) {return 27;} //target 3!
        if (plane == 2) {return 27;}}
    //water?
    //if (mdl == 10) {
    //    if (plane == 1) {
    //        return 1; //weird!!
    //    }
    //    if (plane == 2) {
    //        return 0; //weird!!
    //    }
    //}

    if (mdl == 11) {
        if (plane == 1) {return 28;}
        if (plane == 2) {return 29;}}
    if (mdl == 12) {
        if (plane == 1) {return 30;}
        if (plane == 2) {return 31;}}
    if (mdl == 13) {
        if (plane == 1) {return 32;}
        if (plane == 2) {return 33;}}
    if (mdl == 14) {
        if (plane == 1) {return 34;}
        if (plane == 2) {return 35;}}
    //Here goes Water!!!
    if (mdl == -999) {
        if(plane = -999) {return 36;}}
    if (mdl == 15) {
        if (plane == 1) {return 37;}
        if (plane == 2) {return 38;}}
    if (mdl == 16) {
        if (plane == 1) {return 39;}
        if (plane == 2) {return 40;}}
    if (mdl == 17) {
        if (plane == 1) {return 41;}
        if (plane == 2) {return 42;}}
    if (mdl == 18) {
        if (plane == 1) {return 43;}
        if (plane == 2) {return 44;}}
    if (mdl == 19) {
        if (plane == 1) {return 45; }//Target 4!!}
        if (plane == 2) {return 45;}} //Target 4!!}}
    if (mdl == 20) {
        if (plane == 1) {return 46;}
        if (plane == 2) {return 47;}}
    if (mdl == 21) {
        if (plane == 1) {return 48;}
        if (plane == 2) {return 49;}  }
    if (mdl == 22) {
        if (plane == 1) {return 50; }//target 5!!}
        if (plane == 2) {return 50;}} //target 5!!}  }
    if (mdl == 23) {
        if (plane == 1) {return 51;}
        if (plane == 2) {return 52;}}
    if (mdl == 24) {
        if (plane == 1) {return 53;}
        if (plane == 2) {return 54;}  }
    if (mdl == 25) {
        if (plane == 1) {return 55;}
        if (plane == 2) {return 56;}  }
    if (mdl == 26) {
        if (plane == 1) {return 57;}
        if (plane == 2) {return 58;}  }
    if (mdl == 27) {
        if (plane == 1) {return 59;}
        if (plane == 2) {return 60;}  }
    if (mdl == 28) {
        if (plane == 1) {return 61;}
        if (plane == 2) {return 62;}  }
    if (mdl == 29) {
        if (plane == 1) {return 63;}
        if (plane == 2) {return 64;}  }
    if (mdl == 30) {
        if (plane == 1) {return 65;}
        if (plane == 2) {return 66;}  }
    if (mdl == 31) {
        if (plane == 1) {return 67;}
        if (plane == 2) {return 68;}  }
    if (mdl == 32) {
        if (plane == 1) {return 69;}
        if (plane == 2) {return 70;}  }
    if (mdl == 33) {
        if (plane == 1) {return 71;}
        if (plane == 2) {return 72;}  }
    if (mdl == 34) {
        if (plane == 1) {return 73;}
        if (plane == 2) {return 74;}  }
    if (mdl == 35) {
        if (plane == 1) {return 75;}
        if (plane == 2) {return 76;}  }
    if (mdl == 36) {
        if (plane == 1) {return 77;}
        if (plane == 2) {return 78;}  }
    if (mdl == 37) {
        if (plane == 1) {return 79;}
        if (plane == 2) {return 80;}  }
    if (mdl == 38) {
        if (plane == 1) {return 81;}
        if (plane == 2) {return 82;}  }
    if (mdl == 39) {
        if (plane == 1) {return 83;}
        if (plane == 2) {return 84;}  }
    if (mdl == 40) {
        if (plane == 1) {return 85;}
        if (plane == 2) {return 86;}  }
    if (mdl == 41) {
        if (plane == 1) {return 87;}
        if (plane == 2) {return 88;}  }
    if (mdl == 42) {
        if (plane == 1) {return 89;}
        if (plane == 2) {return 90;}  }
    if (mdl == 43) {
        if (plane == 1) {return 91;}
        if (plane == 2) {return 92;}  }
    if (mdl == 44) {
        if (plane == 1) {return 93;}
        if (plane == 2) {return 94;}  }
    if (mdl == 45) {
        if (plane == 1) {return 95;}
        if (plane == 2) {return 96;}  }
    if (mdl == 46) {
        if (plane == 1) {return 97;}
        if (plane == 2) {return 98;}  }
    if (mdl == 47) {
        if (plane == 1) {return 99;}
        if (plane == 2) {return 100;}  }
    if (mdl == 48) {
        if (plane == 1) {return 101;}
        if (plane == 2) {return 102;}  }
    if (mdl == 49) {
        if (plane == 1) {return 103;}
        if (plane == 2) {return 104;}  }
    if (mdl == 50) {
        if (plane == 1) {return 105;}
        if (plane == 2) {return 106;}  }
    if (mdl == 51) {
        if (plane == 1) {return 107;}
        if (plane == 2) {return 108;}  }
    if (mdl == 52) {
        if (plane == 1) {return 109;}
        if (plane == 2) {return 110;}  }
    if (mdl == 53) {
        if (plane == 1) {return 111;}
        if (plane == 2) {return 112;}  }
    if (mdl == 54) {
        if (plane == 1) {return 113;}
        if (plane == 2) {return 114;}  }
    if (mdl == 55) {
        if (plane == 1) {return 115;}
        if (plane == 2) {return 116;}  }
    if (mdl == 56) {
        if (plane == 1) {return 117;}
        if (plane == 2) {return 118;}  }
    if (mdl == 57) {
        if (plane == 1) {return 119;}
        if (plane == 2) {return 120;}  }
    if (mdl == 58) {
        if (plane == 1) {return 121;}
        if (plane == 2) {return 122;}  }
    if (mdl == 59) {
        if (plane == 1) {return 123;}
        if (plane == 2) {return 124;}  }
    if (mdl == 60) {
        if (plane == 1) {return 125;}
        if (plane == 2) {return 126;}  }
    if (mdl == 61) {
        if (plane == 1) {return 127;}
        if (plane == 2) {return 128;}  }
    if (mdl == 62) {
        if (plane == 1) {return 129;}
        if (plane == 2) {return 130;}  }
    if (mdl == 63) {
        if (plane == 1) {return 131;}
        if (plane == 2) {return 132;}  }
    if (mdl == 64) {
        if (plane == 1) {return 133;}
        if (plane == 2) {return 134;}  }
    if (mdl == 65) {
        if (plane == 1) {return 135;}
        if (plane == 2) {return 136;}  }
    if (mdl == 66) {
        if (plane == 137) {return 1;}
        if (plane == 2) {return 138;}  }
    if (mdl == 67) {
        if (plane == 1) {return 139;}
        if (plane == 2) {return 140;}  }
    if (mdl == 68) {
        if (plane == 1) {return 141;}
        if (plane == 2) {return 142;}  }
    if (mdl == 69) {
        if (plane == 1) {return 143;}
        if (plane == 2) {return 144;}  }
    if (mdl == 70) {
        if (plane == 1) {return 145;}
        if (plane == 2) {return 146;}  }
    if (mdl == 71) {
        if (plane == 1) {return 147;}
        if (plane == 2) {return 148;}  }
    if (mdl == 72) {
        if (plane == 1) {return 149;}
        if (plane == 2) {return 150;}  }
    if (mdl == 73) {
        if (plane == 1) {return 151;}
        if (plane == 2) {return 152;}  }
    if (mdl == 74) {
        if (plane == 1) {return 153;}
        if (plane == 2) {return 154;}  }
    if (mdl == 75) {
        if (plane == 1) {return 155;}
        if (plane == 2) {return 156;}
    }
    if (mdl == 76) {
        if (plane == 1) {return 157;}
        if (plane == 2) {return 158;}  }
    if (mdl == 77) {
        if (plane == 1) {return 159;}
        if (plane == 2) {return 160;}  }
    if (mdl == 78) {
        if (plane == 1) {return 161;}
        if (plane == 2) {return 162;}  }
    if (mdl == 79) {
        if (plane == 1) {return 163;}
        if (plane == 2) {return 164;}  }
    if (mdl == 80) {
        if (plane == 1) {return 165;}
        if (plane == 2) {return 166;}  }
    if (mdl == 81) {
        if (plane == 1) {return 167;}
        if (plane == 2) {return 168;}  }
    if (mdl == 82) {
        if (plane == 1) {return 169;}
        if (plane == 2) {return 170;}  }
    if (mdl == 83) {
        if (plane == 1) {return 171;}
        if (plane == 2) {return 172;}  }
    if (mdl == 84) {
        if (plane == 1) {return 173;}
        if (plane == 2) {return 174;}  }
    if (mdl == 85) {
        if (plane == 1) {return 175;}
        if (plane == 2) {return 176;}  }
    if (mdl == 86) {
        if (plane == 1) {return 177;}
        if (plane == 2) {return 178;}  }
    if (mdl == 87) {
        if (plane == 1) {return 179;}
        if (plane == 2) {return 180;}  }
    if (mdl == 88) {
        if (plane == 1) {return 181;}
        if (plane == 2) {return 182;}  }
    if (mdl == 89) {
        if (plane == 1) {return 183;}
        if (plane == 2) {return 184;}  }
    if (mdl == 90) {
        if (plane == 1) {return 185;}
        if (plane == 2) {return 186;}  }
    if (mdl == 91) {
        if (plane == 1) {return 187;}
        if (plane == 2) {return 188;}  }
    if (mdl == 92) {
        if (plane == 1) {return 189;}
        if (plane == 2) {return 190;}  }
    if (mdl == 93) {
        if (plane == 1) {return 191;}
        if (plane == 2) {return 192;}  }
    if (mdl == 94) {
        if (plane == 1) {return 193;}
        if (plane == 2) {return 194;}  }
    if (mdl == 95) {
        if (plane == 1) {return 195;}
        if (plane == 2) {return 195;}  }
    if (mdl == 96) {
        if (plane == 1) {return 196;}
        if (plane == 2) {return 196;}  }
    if (mdl == 97) {
        if (plane == 1) {return 197;}
        if (plane == 2) {return 197;}  }
    if (mdl == 98) {
        if (plane == 1) {return 198;} //not needed}
        if (plane == 2) {return 198;}  }
    if (mdl == 99) {
        if (plane == 1) {return 199;} //not needed}
        if (plane == 2) {return 199;}  }
    if (mdl == 100) {
        if (plane == 2) {return 200;}  }
    if (mdl == 101) {
        if (plane == 2) {return 201;}  }
    if (mdl == 102) {
        if (plane == 2) {return 202;}  }
    if (mdl == 103) {
        if (plane == 2) {return 203;}  }
    if (mdl == 104) {
        if (plane == 2) {return 204;}  }
    if (mdl == 105) {
        if (plane == 2) {return 205;}  }
    if (mdl == 106) {
        if (plane == 2) {return 206;}  }
    if (mdl == 107) {
        if (plane == 2) {return 207;}  }
    if (mdl == 108) {
        if (plane == 2) {return 208;}  }
    if (mdl == 109) {
        if (plane == 2) {return 209;}  }
    if (mdl == 110) {
        if (plane == 2) {return 210;}  }
    if (mdl == 111) {
        if (plane == 2) {return 211;}  }
    if (mdl == 112) {
        if (plane == 2) {return 212;}  }
    if (mdl == 113) {
        if (plane == 2) {return 213;}  }
    if (mdl == 114) {
        if (plane == 2) {return 214;}  }
    if (mdl == 115) {
        if (plane == 2) {return 2;}  }
    if (mdl == 116) {
        if (plane == 2) {return 0;}  }
    if (mdl == 117) {
        if (plane == 2) {return 1;}  }
    if (mdl == 118) {
        if (plane == 2) {return 0;}  }
    if (mdl == 119) {
        if (plane == 2) {return 2;}  }

    return -999;
  }

  bool hasMLPred() const {return (GetInt("hasMLPrediction")==1);}

  int GetMultiplicity() const {return GetInt("multiplicity");}

  //Still needed for some systematics to compile, but shouldn't be used for reweighting anymore.
  protected:
  #include "PlotUtils/WeightFunctions.h" // Get*Weight
};

#endif
