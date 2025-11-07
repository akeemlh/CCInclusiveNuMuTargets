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
#include "PlotUtils/CaloCorrection.h"

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

  double GetANNPmu() const //MeV/c
  {
    return sqrt(pow(GetVecElem("muon_corrected_p",0),2)+pow(GetVecElem("muon_corrected_p",1),2)+pow(GetVecElem("muon_corrected_p",2),2));
  }

  double GetANNPmuGeV() const //GeV/c
  {
    return GetANNPmu()/1000;
  }

  double GetANNMuonPT() const //MeV/c
  {
    return GetANNPmu() * sin(GetThetamu());
  }

  double GetANNMuonPTGeV() const //GeV/c
  {
    return GetANNMuonPT()/1000;
  }

  double GetANNMuonPz() const //MeV/c
  {
    return GetANNPmu() * cos(GetThetamu());
  }

  double GetANNMuonPzGeV() const //GeV/c
  {
    return GetANNMuonPz()/1000;
  }

  double GetANNEmu() const //MeV
  {
    return GetVecElem("muon_corrected_p", 3);
  }

  double GetANNEmuGeV() const //GeV
  {
    return GetANNEmu()/1000;
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
    return GetDouble("MasterAnaDev_ANN_recoil_E");
  }

  virtual double GetANNRecoilEGeV() const { //GeV
    return GetANNRecoilE()/1000;
  }

  virtual double GetANNEnu() const { //GeV
    return GetANNRecoilE() + GetANNEmu();
  }

  virtual double GetANNEnuGeV() const { //GeV
    return GetANNRecoilEGeV() + GetANNEmuGeV();
  }










 
double getZPosFromSegment(int segment) const
{
    if (segment ==1) return 4293.04;
    if (segment == 2) return 4313.68;
    if (segment == 3) return 4337.25;
    if (segment == 4) return 4357.9;
    if (segment == 5) return 4381.47;
    if (segment == 6) return 4402.11;
    if (segment == 7) return 4425.68;
    if (segment == 8) return 4446.33;
    if (segment == 9) return 4481.21; //Target 1
    if (segment == 10) return 4514.11;
    if (segment == 11) return 4534.76;
    if (segment == 12) return 4558.33;
    if (segment == 13) return 4578.97;
    if (segment == 14) return 4602.54;
    if (segment == 15) return 4623.19;
    if (segment == 16) return 4646.76;
    if (segment == 17) return 4667.4;
    if (segment == 18) return 4702.29; //Target 2
    if (segment == 19) return 4735.19;
    if (segment == 20) return 4755.83;
    if (segment == 21) return 4779.4;
    if (segment == 22) return 4800.05;
    if (segment == 23) return 4823.62;
    if (segment == 24) return 4844.26;
    if (segment == 25) return 4867.83;
    if (segment == 26) return 4888.48;
    if (segment == 27) return 4923.36; //Target 3
    if (segment == 28) return 5000.48;
    if (segment == 29) return 5021.12;
    if (segment == 30) return 5044.69;
    if (segment == 31) return 5065.34;
    if (segment == 32) return 5088.91;
    if (segment == 33) return 5109.55;
    if (segment == 34) return 5133.12;
    if (segment == 35) return 5153.77;
    if (segment == 36) return 5310; //Water Target -- Find better number, this is a guesstimate
    if (segment == 37) return 5456.74;
    if (segment == 38) return 5477.38;
    if (segment == 39) return 5500.95;
    if (segment == 40) return 5521.6;
    if (segment == 41) return 5545.17;
    if (segment == 42) return 5565.81;
    if (segment == 43) return 5589.38;
    if (segment == 44) return 5610.02;
    if (segment == 45) return 5644.91; //Target 4
    if (segment == 46) return 5677.81;
    if (segment == 47) return 5698.45;
    if (segment == 48) return 5722.03;
    if (segment == 49) return 5742.67;
    if (segment == 50) return 5777.55; //Target 5
    if (segment == 51) return 5810.45;
    if (segment == 52) return 5831.1;
    if (segment == 53) return 5855.68;
    if (segment == 54) return 5876.33;
    if (segment == 55) return 5900.91;
    if (segment == 56) return 5921.56;
    if (segment == 57) return 5946.14;
    if (segment == 58) return 5966.79;
    if (segment == 59) return 5991.37;
    if (segment == 60) return 6012.01;
    if (segment == 61) return 6036.6;
    if (segment == 62) return 6057.24;
    if (segment == 63) return 6081.83;
    if (segment == 64) return 6102.47;
    if (segment == 65) return 6127.06;
    if (segment == 66) return 6147.7;
    if (segment == 67) return 6172.29;
    if (segment == 68) return 6192.93;
    if (segment == 69) return 6217.52;
    if (segment == 70) return 6238.16;
    if (segment == 71) return 6262.74;
    if (segment == 72) return 6283.39;
    if (segment == 73) return 6307.97;
    if (segment == 74) return 6328.62;
    if (segment == 75) return 6353.2;
    if (segment == 76) return 6373.85;
    if (segment == 77) return 6398.43;
    if (segment == 78) return 6419.08;
    if (segment == 79) return 6443.66;
    if (segment == 80) return 6464.3;
    if (segment == 81) return 6488.89;
    if (segment == 82) return 6509.53;
    if (segment == 83) return 6534.12;
    if (segment == 84) return 6554.76;
    if (segment == 85) return 6579.35;
    if (segment == 86) return 6599.99;
    if (segment == 87) return 6624.58;
    if (segment == 88) return 6645.22;
    if (segment == 89) return 6669.81;
    if (segment == 90) return 6690.45;
    if (segment == 91) return 6715.03;
    if (segment == 92) return 6735.68;
    if (segment == 93) return 6760.26;
    if (segment == 94) return 6780.91;
    if (segment == 95) return 6805.49;
    if (segment == 96) return 6826.14;
    if (segment == 97) return 6850.72;
    if (segment == 98) return 6871.37;
    if (segment == 99) return 6895.95;
    if (segment == 100) return 6916.59;
    if (segment == 101) return 6941.18;
    if (segment == 102) return 6961.82;
    if (segment == 103) return 6986.41;
    if (segment == 104) return 7007.05;
    if (segment == 105) return 7031.64;
    if (segment == 106) return 7052.28;
    if (segment == 107) return 7076.87;
    if (segment == 108) return 7097.51;
    if (segment == 109) return 7122.1;
    if (segment == 110) return 7142.74;
    if (segment == 111) return 7167.32;
    if (segment == 112) return 7187.97;
    if (segment == 113) return 7212.55;
    if (segment == 114) return 7233.2;
    if (segment == 115) return 7257.78;
    if (segment == 116) return 7278.43;
    if (segment == 117) return 7303.01;
    if (segment == 118) return 7323.66;
    if (segment == 119) return 7348.24;
    if (segment == 120) return 7368.88;
    if (segment == 121) return 7393.47;
    if (segment == 122) return 7414.11;
    if (segment == 123) return 7438.7;
    if (segment == 124) return 7459.34;
    if (segment == 125) return 7483.93;
    if (segment == 126) return 7504.57;
    if (segment == 127) return 7529.16;
    if (segment == 128) return 7549.8;
    if (segment == 129) return 7574.39;
    if (segment == 130) return 7595.03;
    if (segment == 131) return 7619.61;
    if (segment == 132) return 7640.26;
    if (segment == 133) return 7664.84;
    if (segment == 134) return 7685.49;
    if (segment == 135) return 7710.07;
    if (segment == 136) return 7730.72;
    if (segment == 137) return 7755.3;
    if (segment == 138) return 7775.95;
    if (segment == 139) return 7800.53;
    if (segment == 140) return 7821.17;
    if (segment == 141) return 7845.76;
    if (segment == 142) return 7866.4;
    if (segment == 143) return 7890.99;
    if (segment == 144) return 7911.63;
    if (segment == 145) return 7936.22;
    if (segment == 146) return 7956.86;
    if (segment == 147) return 7981.45;
    if (segment == 148) return 8002.09;
    if (segment == 149) return 8026.68;
    if (segment == 150) return 8047.32;
    if (segment == 151) return 8071.9;
    if (segment == 152) return 8092.55;
    if (segment == 153) return 8117.13;
    if (segment == 154) return 8137.78;
    if (segment == 155) return 8162.36;
    if (segment == 156) return 8183.01;
    if (segment == 157) return 8207.59;
    if (segment == 158) return 8228.24;
    if (segment == 159) return 8252.82;
    if (segment == 160) return 8273.46;
    if (segment == 161) return 8298.05;
    if (segment == 162) return 8318.69;
    if (segment == 163) return 8343.28;
    if (segment == 164) return 8363.92;
    if (segment == 165) return 8388.51;
    if (segment == 166) return 8409.15;
    if (segment == 167) return 8433.74;
    if (segment == 168) return 8454.38;
    if (segment == 169) return 8478.97;
    if (segment == 170) return 8499.61;
    if (segment == 171) return 8524.19;
    if (segment == 172) return 8544.84;
    if (segment == 173) return 8569.42;
    if (segment == 174) return 8590.07;
    if (segment == 175) return 8614.65;
    if (segment == 176) return 8635.3;
    if (segment == 177) return 8659.46;
    if (segment == 178) return 8680.1;
    if (segment == 179) return 8704.26;
    if (segment == 180) return 8724.9;
    if (segment == 181) return 8749.06;
    if (segment == 182) return 8769.71;
    if (segment == 183) return 8793.86;
    if (segment == 184) return 8814.51;
    if (segment == 185) return 8838.67;
    if (segment == 186) return 8859.31;
    if (segment == 187) return 8883.47;
    if (segment == 188) return 8904.11;
    if (segment == 189) return 8928.27;
    if (segment == 190) return 8948.92;
    if (segment == 191) return 8973.08;
    if (segment == 192) return 8993.72;
    if (segment == 193) return 9017.88;
    if (segment == 194) return 9038.52;
    if (segment == 195) return 9088.08;
    if (segment == 196) return 9135.41;
    if (segment == 197) return 9182.75;
    if (segment == 198) return 9230.08;
    if (segment == 199) return 9277.41;
    if (segment == 200) return 9324.74;
    if (segment == 201) return 9372.08;
    if (segment == 202) return 9419.41;
    if (segment == 203) return 9466.74;
    if (segment == 204) return 9514.07;
    if (segment == 205) return 9561.41;
    if (segment == 206) return 9608.74;
    if (segment == 207) return 9656.07;
    if (segment == 208) return 9703.4;
    if (segment == 209) return 9750.74;
    if (segment == 210) return 9798.07;
    if (segment == 211) return 9845.4;
    if (segment == 212) return 9892.73;
    if (segment == 213) return 9940.07;
    if (segment == 214) return 9987.4;
}











bool IsInHexagonTrue(double apothem = 850. ) const
{
    double x = GetVecElem("mc_vtx",0);
    double y = GetVecElem("mc_vtx",1);
    double lenOfSide = apothem*(2/sqrt(3)); 
    double slope     = (lenOfSide/2.0)/apothem;
    double xp        = fabs(x);
    double yp        = fabs(y);
    
    if( (xp*xp + yp*yp) < apothem*apothem )             return true;
    else if( xp <= apothem && yp*yp < lenOfSide/2.0 )   return true; 
    else if( xp <= apothem && yp < lenOfSide-xp*slope ) return true;

    return false;
}

// Functions for coherent pion reweight
// Anezka 24/04/2023 using Mehreen's def

double GetTrueHighEpi() const {
    int nFSpi = GetInt("mc_nFSPart");
    double pionE = -1.0;
    double pionKE = -1.0;
    for (int i = 0; i < nFSpi; i++){
        int pdg = GetVecElem("mc_FSPartPDG",i);
        if(pdg != -211) continue;
        double energy = GetVecElem("mc_FSPartE", i);
        double mass = 139.569;
        double tpi = energy - mass;
        if (tpi >= pionKE){
            pionKE = tpi;
            pionE = energy;
        }
    } 
    //std::cout << "Printing energy of pion " << pionE << std::endl;
    return pionE; // MeV
}

double thetaWRTBeam(double x, double y, double z) const{
    double pyp = -1.0 *sin( MinervaUnits::numi_beam_angle_rad)*z + cos( MinervaUnits::numi_beam_angle_rad )*y;
    double pzp = cos( MinervaUnits::numi_beam_angle_rad )*z + sin( MinervaUnits::numi_beam_angle_rad )*y;
    double denom2 = pow(x,2) + pow(pyp,2) + pow(pzp,2);
    if( 0. == denom2 ) return -9999.;
    else return acos(pzp / sqrt(denom2) );
}

double GetTrueAngleHighTpi() const {
    int nFSpi = GetInt("mc_nFSPart");
    double angle = -9999.; //WRTbeam and in degrees
    double pionKE = 0.0;
    int idk = -9999;
    for (int i = 0; i < nFSpi; i++){
        int pdg = GetVecElem("mc_FSPartPDG",i);
        if(pdg != -211) continue;
        double energy = GetVecElem("mc_FSPartE", i);
        double mass = 139.569;
        double tpi = energy - mass;
        if (tpi >= pionKE) {
            pionKE = tpi;
            TVector3 pimomentumvec(GetVecElem("mc_FSPartPx", i), GetVecElem("mc_FSPartPy", i),GetVecElem("mc_FSPartPz", i));
            double deg_wrtb = thetaWRTBeam(pimomentumvec.X(), pimomentumvec.Y(), pimomentumvec.Z()); //rad
        
            angle = deg_wrtb; //*180./M_PI;
        }
    }
    //Making sure angle is only between 0 and pi
    if (angle < 0.0) angle = -1.0*angle;
    if (angle > M_PI) angle = 2.0*M_PI - angle; 
    return angle*180./M_PI; // Degrees
}

// from Mehreen, Ben and Everardo but update with constant defaults from TargetUtils
bool IsInPlastic() const {
  if (!IsInHexagonTrue())
    return false;  // This is in the calorimeters

  double mc_vtx_z = GetVecElem("mc_vtx", 2);
  if (mc_vtx_z > 8467.0) return false;  // Ditto

  int mc_nuclei = GetInt("mc_targetZ");
  // In the carbon target?  center+offset
  if (fabs(mc_vtx_z - PlotUtils::TargetUtils::Get().GetTarget3CarbonCenterZMC()) <=
          PlotUtils::TargetProp::ThicknessMC::Tgt3::C / 2 &&
      mc_nuclei == 6)
    return false;

  // In the water target?
  if (PlotUtils::TargetUtils::Get().InWaterTargetMC(
                            GetVecElem("mc_vtx", 0), GetVecElem("mc_vtx", 1),
                            GetVecElem("mc_vtx", 2), GetInt("mc_targetZ")))
    return false;

  // Finally, do you have target material?  I'm going to say lead/iron isn't a
  // big consideration elsewhere in the detector
  if (mc_nuclei == 26 || mc_nuclei == 82) return false;

  return true;
}

double GetDiffractiveUncWeight() const {
    // inspired by Mehreen: https://github.com/MinervaExpt/LowRecoilPions/blob/902f51bd72e1dff74d26e0df7158f27750947521/event/CVUniverse.h#L787
    // diffractive weight
    // Note, this assumes you're not using the diffractive model in GENIE
    // As of 03/2021, we don't really trust our diffractive model, so
    // as a rough approximation, weight every coherent event
    // (diffractive is coherent on hydrogen) by 1.4368.
    // Coherent xsec scales by A^(1/3), and 1/(12^(1/3)) = 0.4368
    // for water 1/(16^(1/3)) = 0.3969
    
    if(GetInt("mc_intType") == 4){
        if(PlotUtils::TargetUtils::Get().InWaterTargetMC(GetVecElem("mc_vtx",0), GetVecElem("mc_vtx",1), GetVecElem("mc_vtx",0), GetInt("mc_targetZ"))){
            return 1.3969;
        }
        else if (IsInPlastic()){
            return 1.4368;
        }
        else{
            return 1.;
        }
    }
    else{
        return 1.;
    }

}



    virtual double GetCOHPionWeight() const {
        double weight = 1.0;
        if(GetInt("mc_intType") != 4) return 1.0;
        if(GetInt("mc_intType") == 4){
            //int npi = GetTrueNPionsinEvent();
            //if (npi == 0) return 1.0;
            double angle = GetTrueAngleHighTpi();//*180./M_PI; //this is now in degrees
            //Angle is already in degrees GetTrueAngleHighTpi()
            //std::cout << "Printing Angle for COH event " << angle << std::endl;
            double epi = GetTrueHighEpi()/1000.; // This is supposed to be the energy of the pion!!!! 
            if (epi < 0) return 1.0;
            else {
        weight *= GetCoherentPiWeight(angle, epi); //Inputs are in Degrees and GeV
            //std::cout << "Printing COHerent weight for COH event for angle: " << angle << " degrees and KE : " << KE << " GeV. And weight: " << weight << std::endl;
        return weight;
            }
            }
            else return 1.0;
    }


    virtual double GetDiffractiveWeight() const { return GetDiffractiveUncWeight(); } //Redundant with below - could remove/refactor away

    double GetBjorkenX() const
    {
        double Emu = GetANNEmuGeV();
        double mu_massGev = MinervaUnits::M_mu/1000;
        double qsquared =  2.0 * (GetANNRecoilEGeV()+Emu) * (Emu - GetANNPmuGeV() * cos(GetThetamu())) - (mu_massGev*mu_massGev);
        double M = (M_p + M_n) /(2.0*1000); //Not exact. To do - weight it based on target material p n ratio!!!!!!
        double r = GetANNRecoilEGeV();
        double xbj = xBjorken(qsquared, M, r);
        //double q2 = GetInclQ2Reco();
        //double r2 = GetANNRecoilE();
        //double M2 = (M_p + M_n) /2.0;
        //double xbj2 = xBjorken(q2, M2, r2);
        //std::cout<<"Bjorken X calculation: xbj: " <<xbj <<" xbj2 " << xbj2<<" qsquared: "<< qsquared << " M: " << M << " r: " << r << " Emu "<<Emu<< " GetANNRecoilEGeV() " <<GetANNRecoilEGeV() <<" GetANNPmuGeV() " <<GetANNPmuGeV() << " GetThetamu() " << GetThetamu()<< " mu_massGev " << mu_massGev<< std::endl;
        return xbj;
    }

    /* double GetBjorkenX() const
    {
        double q2 = GetInclQ2Reco();
        double M = (M_p + M_n) /2.0; //Not exact. To do - weight it based on target material p n ratio!!!!!!
        double r = GetANNRecoilE();
        double xbj = xBjorken(q2, M, r);
        return xbj;
    } */

    double GetBjorkenY() const
    {
        return GetANNRecoilE()/GetANNEnu();
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
    if (GetInt("ANN_vtx_sz") != 3) return ROOT::Math::XYZVector(-999,-999,-999); //To prevent segfaults when accessing events for which there is no ANN_vtx 3-vector
    result.SetCoordinates(GetVec<double>("ANN_vtx").data());
    return result;
  }

  std::vector<double> GetANNVertexVector() const { return GetVecDouble("ANN_vtx");}

  double GetANNSegment() const
  {
    return GetVecElemInt("ANN_segments", 0);
  }

  double GetANNSegmentsWeighted() const
  {
    return GetVecElemInt("ANN_segments", 0)*GetVecElem("ANN_plane_probs", 0) + GetVecElemInt("ANN_segments", 1)*GetVecElem("ANN_plane_probs", 1);
  }

  double GetANNSegmentsZPosWeighted(double cutoff) const
  {
    double zPosFromSegment0 = getZPosFromSegment(GetVecElemInt("ANN_segments", 0));
    double zPosFromSegment1 = getZPosFromSegment(GetVecElemInt("ANN_segments", 1));
    double conf0 = GetVecElem("ANN_plane_probs", 0);
    double conf1 = GetVecElem("ANN_plane_probs", 1);
    std::cout<<"cutoff: "<<cutoff<< " zPosFromSegment0: " << zPosFromSegment0 << " conf0 "<< conf0<< " zPosFromSegment1: " << zPosFromSegment1 << " conf1 "<< conf1 <<  std::endl;
    //std::cout<<"zPosFromSegment1: " << zPosFromSegment1 << " conf1 "<< conf1 <<  std::endl;
    //if (conf1 < 0.2 ) return 0;
    if (conf1 >= cutoff) return (zPosFromSegment0*conf0*conf0+zPosFromSegment1*conf1*conf1)/(conf0*conf0+conf1*conf1);
    else return zPosFromSegment0;
  }

   double GetANNSegmentsZPosWeighted2(int cutoff) const
  {
    if (GetInt("hasMLPrediction")!=1) return -999;
    double zPosFromSegment0 = getZPosFromSegment(GetVecElemInt("ANN_segments", 0));
    double zPosFromSegment1 = getZPosFromSegment(GetVecElemInt("ANN_segments", 1));
    double conf0 = GetVecElem("ANN_plane_probs", 0);
    double conf1 = GetVecElem("ANN_plane_probs", 1);
    //std::cout<<"cutoff: "<<cutoff<< " Segment0: " << GetVecElemInt("ANN_segments", 0) << " conf0 "<< conf0<< " Segment1: " << GetVecElemInt("ANN_segments", 1) << " zPosFromSegment0: " << zPosFromSegment0 << " conf0 "<< conf0<< " zPosFromSegment1: " << zPosFromSegment1 << " conf1 "<< conf1 <<  std::endl;
    if (std::abs(GetVecElemInt("ANN_segments", 0)-GetVecElemInt("ANN_segments", 1)) < cutoff)
    {
        return (zPosFromSegment0*conf0+zPosFromSegment1*conf1)/(conf0+conf1);
        //return (zPosFromSegment0*conf0*conf0+zPosFromSegment1*conf1*conf1)/(conf0*conf0+conf1*conf1);
    }
    else return zPosFromSegment0;
  }


    // ========================================================================
    //  ENERGY FUNCTIONS
    // ========================================================================
    /* double ApplyCaloTuning(double calRecoilE) const;
    double GetCalRecoilEnergy() const{ return GetANNRecoilE();}
    double GetNonCalRecoilEnergy() const{
        return 0;
    }
    double GetRecoilEnergy() const{
        return ApplyCaloTuning(GetCalRecoilEnergy()) + GetNonCalRecoilEnergy();  //  LOOK AT PlotUtils/SystCalcs/RecoilEnergyFunctions.h

    } */


    ////////////////////////////////////////////////////////////////////////
    // LOOK HERE
    // Need to define these functions in your personal CVUniverse

    //virtual double GetCalRecoilEnergy() const {
    //  std::cout << "GetCalRecoilEnergy() should be implemented in CVUniverse to"<< std::endl;
    //  std::cout << "     return all recoil energy found calorimetrically (MeV)"
    //  std::cout << " If all of your recoil energy is find calorimetrically (such as E available
    //  std::cout << "   or by spline correction, just put that branch/ calculation here"
    //            << std::endl;
    //      
    //  return calorimetric_recoil_energy;
    //}
    //
    //virtual double GetNonCalRecoilEnergy() const {
    //  std::cout << "GetNonCalRecoilEnergy() should be implemented in CVUniverse to"<< std::endl;
    //  std::cout << "     return all recoil energy not found calorimetrically like dEdX (MeV)"
    //            << std::endl;
    //  return non_calorimetric_recoil_energy;
    //}

 //Taken from Oscar's code (/exp/minerva/app/users/omorenop/cmtuser/git-Mat/Personal/Test/InclusiveUtils.h)
  double GetTruthSegment() const //Dangerous since ANN_vtx doesnt always have a 3-vector in it, calling SetCoordinates without a 3 vector gives a segfault
  {
    int mdl = GetTruthVtxModule();
    int plane = GetTruthVtxPlane();

    int tgtZ = GetTruthTargetZ();
    int tgtID = GetTruthTargetID();

    if (tgtID == 1 ) {return 9;} //target 1!
    if (tgtID == 2 ) {return 18;} //target 2!
    if (tgtID == 3 ) {return 27;} //target 3!
    if (tgtID == 4 ) {return 45;} //target 4!
    if (tgtID == 5 ) {return 50;} //target 5!
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
        if (plane == 1) {return 137;}
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
