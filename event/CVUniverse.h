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
#include "utilities/SMSCalibrator.h"
#include <cmath>

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

  SMSCalibrator *calibrator =nullptr;
  bool calibratorSet = false;

  void SetCalibrator(SMSCalibrator* cal)
  {
    calibrator = cal;
    calibratorSet = true;
  }

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

  virtual double GetANNEmu() const //MeV
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

  //Akeem - Note to self, check this over? Is this what we want and is it overridden in systematic universes?
  virtual double GetCalRecoilEnergy() const {
    //return GetDouble("part_response_total_recoil_passive_allNonMuonClusters_id")+ GetDouble("part_response_total_recoil_passive_allNonMuonClusters_od"); // in MeV
    return GetDouble((GetAnaToolName() + "_ANN_recoil_E").c_str());  // [MeV]
    // usual inclusive recoil energy should just be calorimetric energy
    // instead of Gonzalo's return GetDouble("Recoil_Ecalo") / 1000.0;  // [GeV]
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




ROOT::Math::XYZVector adjustDNNVertex( double zPos ) const
{
    if (zPos <= 0) return ROOT::Math::XYZVector( -999, -999, -999 );
    std::vector<double> annVtx = GetVecDouble("ANN_vtx");
    double mu_px = GetDouble("MasterAnaDev_muon_Px");
    double mu_py = GetDouble("MasterAnaDev_muon_Py");
    double mu_pz = GetDouble("MasterAnaDev_muon_Pz");
    double mu_time = GetDouble("muon_trackVertexTime");

    double projx = annVtx[0] + ( mu_px / mu_pz ) * ( zPos - annVtx[2] );
    double projy = annVtx[1] + ( mu_py / mu_pz ) * ( zPos - annVtx[2] );
    ROOT::Math::XYZVector outVtx = ROOT::Math::XYZVector( projx, projy, zPos );
    return outVtx;
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
    else
    {
        //std::cout<<"Invalid segment passed to getZPosFromSegment: " << segment << std::endl;
        return -999;
    }
}



int getNearestSegmentToZPos(double zpos) const  //A very rough way of getting the segment number from a given z vertex
{
    //The centre of the scintillator planes, from (https://cdcvs.fnal.gov/redmine/projects/minerva/wiki/Z_Positions_of_Planes_in_the_Full_MINERvA_Detector)
    std::vector<double> plane_centres = {4293.04,
    4313.68,
    4337.25,
    4357.9,
    4381.47,
    4402.11,
    4425.68,
    4446.33,
    4514.11,
    4534.76,
    4558.33,
    4578.97,
    4602.54,
    4623.19,
    4646.76,
    4667.4,
    4735.19,
    4755.83,
    4779.4,
    4800.05,
    4823.62,
    4844.26,
    4867.83,
    4888.48,
    5000.48,
    5021.12,
    5044.69,
    5065.34,
    5088.91,
    5109.55,
    5133.12,
    5153.77,
    5456.74,
    5477.38,
    5500.95,
    5521.6,
    5545.17,
    5565.81,
    5589.38,
    5610.02,
    5677.81,
    5698.45,
    5722.03,
    5742.67,
    5810.45,
    5831.1,
    5855.68,
    5876.33,
    5900.91,
    5921.56,
    5946.14,
    5966.79,
    5991.37,
    6012.01,
    6036.6,
    6057.24,
    6081.83,
    6102.47,
    6127.06,
    6147.7,
    6172.29,
    6192.93,
    6217.52,
    6238.16,
    6262.74,
    6283.39,
    6307.97,
    6328.62,
    6353.2,
    6373.85,
    6398.43,
    6419.08,
    6443.66,
    6464.3,
    6488.89,
    6509.53,
    6534.12,
    6554.76,
    6579.35,
    6599.99,
    6624.58,
    6645.22,
    6669.81,
    6690.45,
    6715.03,
    6735.68,
    6760.26,
    6780.91,
    6805.49,
    6826.14,
    6850.72,
    6871.37,
    6895.95,
    6916.59,
    6941.18,
    6961.82,
    6986.41,
    7007.05,
    7031.64,
    7052.28,
    7076.87,
    7097.51,
    7122.1,
    7142.74,
    7167.32,
    7187.97,
    7212.55,
    7233.2,
    7257.78,
    7278.43,
    7303.01,
    7323.66,
    7348.24,
    7368.88,
    7393.47,
    7414.11,
    7438.7,
    7459.34,
    7483.93,
    7504.57,
    7529.16,
    7549.8,
    7574.39,
    7595.03,
    7619.61,
    7640.26,
    7664.84,
    7685.49,
    7710.07,
    7730.72,
    7755.3,
    7775.95,
    7800.53,
    7821.17,
    7845.76,
    7866.4,
    7890.99,
    7911.63,
    7936.22,
    7956.86,
    7981.45,
    8002.09,
    8026.68,
    8047.32,
    8071.9,
    8092.55,
    8117.13,
    8137.78,
    8162.36,
    8183.01,
    8207.59,
    8228.24,
    8252.82,
    8273.46,
    8298.05,
    8318.69,
    8343.28,
    8363.92,
    8388.51,
    8409.15,
    8433.74,
    8454.38,
    8478.97,
    8499.61,
    8524.19,
    8544.84,
    8569.42,
    8590.07,
    8614.65,
    8635.3,
    8659.46,
    8680.1,
    8704.26,
    8724.9,
    8749.06,
    8769.71,
    8793.86,
    8814.51,
    8838.67,
    8859.31,
    8883.47,
    8904.11,
    8928.27,
    8948.92,
    8973.08,
    8993.72,
    9017.88,
    9038.52,
    9088.08,
    9135.41,
    9182.75,
    9230.08,
    9277.41,
    9324.74,
    9372.08,
    9419.41,
    9466.74,
    9514.07,
    9561.41,
    9608.74,
    9656.07,
    9703.4,
    9750.74,
    9798.07,
    9845.4,
    9892.73,
    9940.07,
    9987.4};

    PlotUtils::TargetUtils util;
    if (util.InTarget1ZMC (zpos, 26)) return 9;
    else if (util.InTarget2ZMC (zpos, 26)) return 18;
    else if (util.InTarget3ZMC (zpos, 6)) return 27;
    else if (util.InWaterTargetZMC (zpos)) return 36;
    else if (util.InTarget4ZMC (zpos, 26)) return 45;
    else if (util.InTarget5ZMC (zpos, 26)) return 50;

    if (zpos < 4290.11) return 0;
    if (zpos > 10034.73) return 0;
    auto it = std::lower_bound(plane_centres.begin(), plane_centres.end(), zpos);

    // Case: zpos is smaller than or equal to the first element
    if (it == plane_centres.begin()) {
        return 0;
    }

    // Case: zpos is strictly larger than the last element
    if (it == plane_centres.end()) {
        return plane_centres.size() - 1;
    }

    // Case: zpos is between two elements
    // We compare zpos to the element at *it (the upper bound)
    // and the element at *(it - 1) (the lower bound)
    double valAbove = *it;
    double valBelow = *(it - 1);

    int retVal = -999;
    if (std::abs(zpos - valBelow) < std::abs(zpos - valAbove)) {
        // The element BEFORE the iterator is closer
        int idx = std::distance(plane_centres.begin(), it) - 1;
        retVal = (idx+1); // +1 because segments and indices are shifted by 1, i.e plane_centres[0] is segment 1 not 0
    } else {
        // The element AT the iterator is closer
        int idx = std::distance(plane_centres.begin(), it);
        retVal = (idx+1); // +1 because segments and indices are shifted by 1, i.e plane_centres[0] is segment 1 not 0
    }
    //Shifting by +1 for every target this event is more downstream than to account for the fact that the nuclear targets
    //take up a segment number each but are not represented in the vector of scintillator plane centre positions above
    if (retVal > 9) retVal += 1;
    if (retVal > 18) retVal += 1;
    if (retVal > 27) retVal += 1;
    if (retVal > 36) retVal += 1;
    if (retVal > 45) retVal += 1;
    if (retVal > 50) retVal += 1;
    return retVal;
}



int getAdjustedTargetCode(double xpos, double ypos, double zpos) const
{
    PlotUtils::TargetUtils util;
    if (util.InIron1VolMC(xpos, ypos, zpos)) return 1026;
    else if (util.InLead1VolMC(xpos, ypos, zpos)) return 1082;
    else if (util.InIron2VolMC(xpos, ypos, zpos)) return 2026;
    else if (util.InLead2VolMC(xpos, ypos, zpos)) return 2082;
    else if (util.InCarbon3VolMC(xpos, ypos, zpos)) return 3006;
    else if (util.InIron3VolMC(xpos, ypos, zpos)) return 3026;
    else if (util.InLead3VolMC(xpos, ypos, zpos)) return 3082;
    else if (util.InLead4VolMC(xpos, ypos, zpos )) return 4082;
    else if (util.InIron5VolMC(xpos, ypos, zpos)) return 5026;
    else if (util.InLead5VolMC(xpos, ypos, zpos)) return 5082;
    else return -999; //If you're not in a target you don't get a target code
}








//Taken from MasterAnaDev.cpp
int GetTargetFromSegment( int segment, int& vtx_module, int& vtx_plane ) const
{
    int targetID = -1;
    vtx_module = -999;
    vtx_plane = -999;

    if( segment == 0 ){
        vtx_module = -5;
        vtx_plane = 0;
        targetID = -1;
    }
    else if( segment == 1 ){
        vtx_module = -5;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 2 ){
        vtx_module = -5;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 3 ){
        vtx_module = -4;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 4 ){
        vtx_module = -4;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 5 ){
        vtx_module = -3;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 6 ){
        vtx_module = -3;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 7 ){
        vtx_module = -2;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 8 ){
        vtx_module = -2;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 9 ){
        vtx_module = -1;
        vtx_plane = 1;
        targetID = 1;
    }
    else if( segment == 10 ){
        vtx_module = 0;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 11 ){
        vtx_module = 0;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 12 ){
        vtx_module = 1;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 13 ){
        vtx_module = 1;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 14 ){
        vtx_module = 2;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 15 ){
        vtx_module = 2;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 16 ){
        vtx_module = 3;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 17 ){
        vtx_module = 3;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 18 ){
        vtx_module = 4;
        vtx_plane = 1;
        targetID = 2;
    }
    else if( segment == 19 ){
        vtx_module = 5;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 20 ){
        vtx_module = 5;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 21 ){
        vtx_module = 6;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 22 ){
        vtx_module = 6;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 23 ){
        vtx_module = 7;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 24 ){
        vtx_module = 7;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 25 ){
        vtx_module = 8;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 26 ){
        vtx_module = 8;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 27 ){
        vtx_module = 9;
        vtx_plane = 2;
        targetID = 3;
    }
    else if( segment == 28 ){
        vtx_module = 11;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 29 ){
        vtx_module = 11;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 30 ){
        vtx_module = 12;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 31 ){
        vtx_module = 12;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 32 ){
        vtx_module = 13;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 33 ){
        vtx_module = 13;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 34 ){
        vtx_module = 14;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 35 ){
        vtx_module = 14;
        vtx_plane = 2;
        targetID = -1;
    }
    //segment 36 is water!!
    /*  else if( segment == 36 ){
        vtx_module = 15;
        vtx_plane = 1;
        targetID = 6;
    }*/
    else if( segment == 36 ){
        vtx_module = -999;
        vtx_plane = -999;
        targetID = 6;
    }
    else if( segment == 37 ){
        vtx_module = 15;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 38 ){
        vtx_module = 15;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 39 ){
        vtx_module = 16;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 40 ){
        vtx_module = 16;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 41 ){
        vtx_module = 17;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 42 ){
        vtx_module = 17;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 43 ){
        vtx_module = 18;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 44 ){
        vtx_module = 18;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 45 ){
        vtx_module = 19;
        vtx_plane = 1;
        targetID = 4;
    }
    else if( segment == 46 ){
        vtx_module = 20;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 47 ){
        vtx_module = 20;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 48 ){
        vtx_module = 21;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 49 ){
        vtx_module = 21;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 50 ){
        vtx_module = 22;
        vtx_plane = 1;
        targetID = 5;
    }
    else if( segment == 51 ){
        vtx_module = 23;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 52 ){
        vtx_module = 23;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 53 ){
        vtx_module = 24;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 54 ){
        vtx_module = 24;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 55 ){
        vtx_module = 25;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 56 ){
        vtx_module = 25;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 57 ){
        vtx_module = 26;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 58 ){
        vtx_module = 26;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 59 ){
        vtx_module = 27;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 60 ){
        vtx_module = 27;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 61 ){
        vtx_module = 28;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 62 ){
        vtx_module = 28;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 63 ){
        vtx_module = 29;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 64 ){
        vtx_module = 29;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 65 ){
        vtx_module = 30;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 66 ){
        vtx_module = 30;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 67 ){
        vtx_module = 31;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 68 ){
        vtx_module = 31;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 69 ){
        vtx_module = 32;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 70 ){
        vtx_module = 32;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 71 ){
        vtx_module = 33;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 72 ){
        vtx_module = 33;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 73 ){
        vtx_module = 34;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 74 ){
        vtx_module = 34;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 75 ){
        vtx_module = 35;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 76 ){
        vtx_module = 35;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 77 ){
        vtx_module = 36;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 78 ){
        vtx_module = 36;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 79 ){
        vtx_module = 37;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 80 ){
        vtx_module = 37;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 81 ){
        vtx_module = 38;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 82 ){
        vtx_module = 38;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 83 ){
        vtx_module = 39;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 84 ){
        vtx_module = 39;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 85 ){
        vtx_module = 40;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 86 ){
        vtx_module = 40;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 87 ){
        vtx_module = 41;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 88 ){
        vtx_module = 41;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 89 ){
        vtx_module = 42;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 90 ){
        vtx_module = 42;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 91 ){
        vtx_module = 43;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 92 ){
        vtx_module = 43;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 93 ){
        vtx_module = 44;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 94 ){
        vtx_module = 44;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 95 ){
        vtx_module = 45;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 96 ){
        vtx_module = 45;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 97 ){
        vtx_module = 46;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 98 ){
        vtx_module = 46;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 99 ){
        vtx_module = 47;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 100 ){
        vtx_module = 47;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 101 ){
        vtx_module = 48;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 102 ){
        vtx_module = 48;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 103 ){
        vtx_module = 49;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 104 ){
        vtx_module = 49;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 105 ){
        vtx_module = 50;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 106 ){
        vtx_module = 50;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 107 ){
        vtx_module = 51;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 108 ){
        vtx_module = 51;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 109 ){
        vtx_module = 52;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 110 ){
        vtx_module = 52;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 111 ){
        vtx_module = 53;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 112 ){
        vtx_module = 53;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 113 ){
        vtx_module = 54;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 114 ){
        vtx_module = 54;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 115 ){
        vtx_module = 55;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 116 ){
        vtx_module = 55;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 117 ){
        vtx_module = 56;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 118 ){
        vtx_module = 56;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 119 ){
        vtx_module = 57;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 120 ){
        vtx_module = 57;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 121 ){
        vtx_module = 58;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 122 ){
        vtx_module = 58;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 123 ){
        vtx_module = 59;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 124 ){
        vtx_module = 59;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 125 ){
        vtx_module = 60;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 126 ){
        vtx_module = 60;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 127 ){
        vtx_module = 61;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 128 ){
        vtx_module = 61;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 129 ){
        vtx_module = 62;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 130 ){
        vtx_module = 62;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 131 ){
        vtx_module = 63;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 132 ){
        vtx_module = 63;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 133 ){
        vtx_module = 64;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 134 ){
        vtx_module = 64;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 135 ){
        vtx_module = 65;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 136 ){
        vtx_module = 65;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 137 ){
        vtx_module = 66;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 138 ){
        vtx_module = 66;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 139 ){
        vtx_module = 67;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 140 ){
        vtx_module = 67;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 141 ){
        vtx_module = 68;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 142 ){
        vtx_module = 68;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 143 ){
        vtx_module = 69;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 144 ){
        vtx_module = 69;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 145 ){
        vtx_module = 70;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 146 ){
        vtx_module = 70;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 147 ){
        vtx_module = 71;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 148 ){
        vtx_module = 71;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 149 ){
        vtx_module = 72;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 150 ){
        vtx_module = 72;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 151 ){
        vtx_module = 73;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 152 ){
        vtx_module = 73;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 153 ){
        vtx_module = 74;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 154 ){
        vtx_module = 74;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 155 ){
        vtx_module = 75;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 156 ){
        vtx_module = 75;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 157 ){
        vtx_module = 76;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 158 ){
        vtx_module = 76;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 159 ){
        vtx_module = 77;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 160 ){
        vtx_module = 77;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 161 ){
        vtx_module = 78;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 162 ){
        vtx_module = 78;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 163 ){
        vtx_module = 79;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 164 ){
        vtx_module = 79;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 165 ){
        vtx_module = 80;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 166 ){
        vtx_module = 80;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 167 ){
        vtx_module = 81;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 168 ){
        vtx_module = 81;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 169 ){
        vtx_module = 82;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 170 ){
        vtx_module = 82;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 171 ){
        vtx_module = 83;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 172 ){
        vtx_module = 83;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 173 ){
        vtx_module = 84;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 174 ){
        vtx_module = 84;
        vtx_plane = 2;
        targetID = -1;
    }
    //Oscar added this part for DSCAL region
    else if( segment == 175 ){
        vtx_module = 85;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 176 ){
        vtx_module = 85;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 177 ){
        vtx_module = 86;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 178 ){
        vtx_module = 86;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 179 ){
        vtx_module = 87;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 180 ){
        vtx_module = 87;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 181 ){
        vtx_module = 88;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 182 ){
        vtx_module = 88;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 183 ){
        vtx_module = 89;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 184 ){
        vtx_module = 89;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 185 ){
        vtx_module = 90;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 186 ){
        vtx_module = 90;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 187 ){
        vtx_module = 91;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 188 ){
        vtx_module = 91;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 189 ){
        vtx_module = 92;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 190 ){
        vtx_module = 92;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 191 ){
        vtx_module = 93;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 192 ){
        vtx_module = 93;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 193 ){
        vtx_module = 94;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 194){
        vtx_module = 94;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 195 ){
        vtx_module = 95;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 196 ){
        vtx_module = 96;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 197 ){
        vtx_module = 97;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 198 ){
        vtx_module = 98;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 199 ){
        vtx_module = 99;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 200 ){
        vtx_module = 100;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 201 ){
        vtx_module = 101;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 202 ){
        vtx_module = 102;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 203 ){
        vtx_module = 103;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 204 ){
        vtx_module = 104;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 205 ){
        vtx_module = 105;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 206 ){
        vtx_module = 106;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 207 ){
        vtx_module = 107;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 208 ){
        vtx_module = 108;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 209 ){
        vtx_module = 109;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 210 ){
        vtx_module = 110;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 211 ){
        vtx_module = 111;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 212 ){
        vtx_module = 112;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 213 ){
        vtx_module = 113;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 214 ){
        vtx_module = 114;
        vtx_plane = 2;
        targetID = -1;
    }
    if( targetID == -999 ) std::cerr << "Can't find any module/plane for given segment " << segment << std::endl;
    return targetID;
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

  bool hasMLPred() const {return (GetInt("hasMLPrediction")==1);}

  double GetTruthMuE() const { return GetDouble("truth_muon_E") ; }

  int GetTargetZ() const {
    return GetInt("MasterAnaDev_targetZ");
  }

  double GetANNTargetDistToDiv() const {return GetDouble("MasterAnaDev_ANN_target_dist_to_division");}
  double GetTBTargetDistToDiv() const {return GetDouble("MasterAnaDev_target_dist_to_division");}
  double GetTruthTargetDistToDiv() const {return GetDouble("truth_target_dist_to_division");}

  int GetANNTargetZ() const { return GetInt("MasterAnaDev_ANN_targetZ"); }
  int GetTruthTargetZ() const { return GetInt("truth_targetZ");}
  int GetMCTargetZ() const { return GetInt("mc_targetZ");}

  int GetANNTargetID() const {return GetInt("MasterAnaDev_ANN_targetID");}
  int GetTruthTargetID() const {return GetInt("truth_targetID");}

  int GetTruthTargetCode() const {return GetInt("truth_target_code");}

  int GetTruthVtxPlane() const {return GetInt("truth_vtx_plane");}
  int GetMADVtxPlane() const {return GetInt("MasterAnaDev_vtx_plane");}

  int GetTruthVtxModule() const {return GetInt("truth_vtx_module");}
  int GetMADVtxModule() const {return GetInt("MasterAnaDev_vtx_module");}


  SMSCalibrator::CalibratedProbs GetANNSegmentsCalibrated() const
  {
    int currentEntry = GetEntry();
    std::vector<int> segments = {};
    std::vector<double> probs = {};
    if (!calibrator->isCurrentCachedEntry(currentEntry)) //We actually don't need to load the segments and plane_probs if this event has already been evaluated and cached, saves us a few I/O operations
    {    
        segments = GetVecInt("ANN_segments");
        probs = GetVecDouble("ANN_plane_probs");
        if (!calibrator) {
            std::cerr << "CRITICAL ERROR: Calibrator pointer is NULL in Universe " 
                    << ShortName() << std::endl;
            throw std::runtime_error("Calibrator not set");
        }
        //std::cout<<"Uncalib - " << segments[0] << ":" << probs[0] << ", "<< segments[1] << ":" << probs[1] << ", "<< segments[2] << ":" << probs[2] << "\n";
        int seg0 = segments[0];
        //std::cout<<"seg0: " << seg0 << std::endl;
    }
    SMSCalibrator::CalibratedProbs cal = calibrator->Top3Calibrated(segments, probs, currentEntry);
    //std::cout<<"Calib - " << cal.segments[0] << ":" << cal.probs[0] << ", "<< cal.segments[1] << ":" << cal.probs[1] << ", "<< cal.segments[2] << ":" << cal.probs[2] << "\n";
    return cal;
  }

  //------------------------------------------------------------------
  //Calibrator methods - Used when for some reason you want to select a specific classification (0-2)
  //I imagine the one would rarely want to use this so you don't mess with the normal shifting of classifications
  //During the MLVertexClassificationUniverse shifts, but included here for completeness and diagnostics
  //------------------------------------------------------------------

  double GetCalibratedANNProb(int i = 0)
  {
    return GetANNSegmentsCalibrated().probs[i];
  }
  double GetCalibratedANNProbAdjusted(int i = 0)
  {
    SMSCalibrator::CalibratedProbs c = GetANNSegmentsCalibrated();
    if (i == 0) return c.probs[0];
    if (i == 1) return c.probs[1]/(1-c.probs[0]);
    if (i == 2) return c.probs[2]/(1-(c.probs[0]+c.probs[1]));
  }
  int GetCalibratedANNSegment(int i = 0)
  {
    return GetANNSegmentsCalibrated().segments[i];
  }
  double GetCalibratedANNZPos(int i = 0)
  {
    return getZPosFromSegment(GetANNSegmentsCalibrated().segments[i]);
  }
  virtual ROOT::Math::XYZVector GetCalibratedANNVertex(int i = 0)
  {
    double calibratedZ = GetCalibratedANNZPos(i);
    //std::cout<<"calibratedZ: " << calibratedZ << std::endl;
    if (calibratedZ == -999) return ROOT::Math::XYZVector(-999, -999, -999);
    ROOT::Math::XYZVector adjustedVertex = adjustDNNVertex(calibratedZ);
    return adjustedVertex;
  }
  int GetCalibratedANNTargetCode(int i = 0) {
    ROOT::Math::XYZVector vtx = GetCalibratedANNVertex(i);
    //std::cout<<"vtx.Z(): " << vtx.Z() << std::endl;
    if (vtx.Z() == -999) return -999;
    int tgtCode = getAdjustedTargetCode(vtx.X(), vtx.Y(), vtx.Z());
    //std::cout<<"tgtCode: " << tgtCode << std::endl;
    return tgtCode;
  }
  //------------------------------------------------------------------
  //End - Calibrator methods - Used when for some reason you want to select a specific classification (0-2)
  //------------------------------------------------------------------

  //------------------------------------------------------------------
  //Calibrator methods - that I expect you to actually use
  //------------------------------------------------------------------

  virtual double GetANNProb(bool useCalibrator = true) const
  {
    if (calibratorSet && useCalibrator)
    {
        return GetANNSegmentsCalibrated().probs[0];
    }
    else
    {
        return GetVecElem("ANN_plane_probs", 0);
    }
  }

  double GetANNZPosWeightedCalibrator() const
  {
    SMSCalibrator::CalibratedProbs cal = GetANNSegmentsCalibrated();
    double totalProb = 0;
    double totalZPos = 0;
    if (cal.probs[0] == -999) return -999;
    for (int i = 0; i < 3; i++) //Only average over valid segments
    {
        double zpos = getZPosFromSegment(cal.segments[i]);
        if (zpos>0)
        {
            totalZPos+=zpos*cal.probs[i];
            totalProb+=cal.probs[i];
        }
    }
    return totalZPos/totalProb;
  }

  int GetANNSegmentWeighted() const
  {
    double calibratedZ = GetANNZPos();
    if (calibratedZ == -999) return -999;
    int uncalibratedSegment = GetANNSegmentsCalibrated().segments[0];
    double unweightedZ = getZPosFromSegment(uncalibratedSegment);
    //If the difference between the weighted and unweighted top plane is less than the thicnkess of
    // a scintillator plane then no need to make any changes, we can just return that plane
    if (std::abs(calibratedZ-unweightedZ) < 20) return uncalibratedSegment; //20 is not the exact width of a plane but is a good-enough conservative rounding
    else
    {
        return getNearestSegmentToZPos(calibratedZ);
    }
  }

  virtual int GetANNSegment(bool useCalibrator = true) const
  {
    if (calibratorSet && useCalibrator)
    {

        //return GetANNSegmentWeighted(); //Weighted
        return GetANNSegmentsCalibrated().segments[0]; //Unweighted
    }
    else
    {
        return GetVecElem("ANN_segments", 0);
    }
  }

  int GetANNVtxPlane(bool useCalibrator = true) const 
  {
    if (calibratorSet && useCalibrator)
    {
        int segment = GetANNSegment();
        int mod;
        int plane;
        GetTargetFromSegment( segment, mod, plane );
        return plane;
    }
    else
    {
        return GetVecElemInt("ANN_vtx_planes", 0);
    }
  }


  int GetANNVtxModule(bool useCalibrator = true) const 
  {
    if (calibratorSet && useCalibrator)
    {
        int segment = GetANNSegment();
        int mod;
        int plane;
        GetTargetFromSegment( segment, mod, plane );
        return mod;
    }
    else
    {
        return GetVecElemInt("ANN_vtx_modules", 0);
    }
  }

  virtual double GetANNZPos(bool useCalibrator = true) const
  {
    if (calibratorSet && useCalibrator)
    {
        //return GetANNZPosWeightedCalibrator(); //Weighted
        return getZPosFromSegment(GetANNSegment()); //Unweighted
    }
    else
    {
        return GetVecElem("ANN_vtx", 2);
    }
  }

  virtual ROOT::Math::XYZVector GetANNVertex(bool useCalibrator = true) const
  {
    if (calibratorSet && useCalibrator)
    {
        double calibratedZ = GetANNZPos();
        if (calibratedZ == -999) return ROOT::Math::XYZVector(-999, -999, -999);
        ROOT::Math::XYZVector adjustedVertex = adjustDNNVertex(calibratedZ);
        return adjustedVertex;
    }
    else
    {
        ROOT::Math::XYZVector result;
        //std::cout<<"GetInt(\"ANN_vtx_sz\"): " << GetInt("ANN_vtx_sz") <<std::endl;
        if (GetInt("ANN_vtx_sz") != 3) return ROOT::Math::XYZVector(-999,-999,-999); //To prevent segfaults when accessing events for which there is no ANN_vtx 3-vector
        result.SetCoordinates(GetVec<double>("ANN_vtx").data());
        return result;
    }
  }


  virtual int GetANNTargetCode(bool useCalibrator = true) const
  {
    if (calibratorSet && useCalibrator)
    {
        ROOT::Math::XYZVector vtx = GetANNVertex();
        if (vtx.Z() == -999) return -999;
        PlotUtils::TargetUtils util;
        int seg = GetANNSegment();
        std::string name = LatexName();
        //The nuclear targets
        if (seg == 9)
        {
            if (util.InIron1VolMC(vtx.X(), vtx.Y(), util.GetTarget1CenterZMC())) return 1026;
            else if (util.InLead1VolMC(vtx.X(), vtx.Y(), util.GetTarget1CenterZMC())) return 1082;
        } 
        else if (seg == 18)
        {
            if (util.InIron2VolMC(vtx.X(), vtx.Y(), util.GetTarget2CenterZMC())) return 2026;
            else if (util.InLead2VolMC(vtx.X(), vtx.Y(), util.GetTarget2CenterZMC())) return 2082;
        } 
        else if (seg == 27)
        {
            if (util.InCarbon3VolMC(vtx.X(), vtx.Y(), util.GetTarget3CenterZMC())) return 3006;
            else if (util.InIron3VolMC(vtx.X(), vtx.Y(), util.GetTarget3CenterZMC())) return 3026;
            else if (util.InLead3VolMC(vtx.X(), vtx.Y(), util.GetTarget3CenterZMC())) return 3082;
        } 
        else if (seg == 36) {return 6000;}
        else if (seg == 45) {return 4082;}
        else if (seg == 50)
        {
            if (util.InIron5VolMC(vtx.X(), vtx.Y(), util.GetTarget5CenterZMC())) return 5026;
            else if (util.InLead5VolMC(vtx.X(), vtx.Y(), util.GetTarget5CenterZMC())) return 5082;
        }
        else
        {
            return seg;
        }
        return -999;
    }
    else
    {
        return GetInt("MasterAnaDev_ANN_target_code");
    }
  }

  virtual std::vector<double> GetANNVertexVector(bool useCalibrator = true) const
  {
    if (calibratorSet && useCalibrator)
    {
        ROOT::Math::XYZVector vtx = GetANNVertex();
        std::vector outVtx = {vtx.X(), vtx.Y(), vtx.Z()};
        return outVtx;
    }
    else
    {
        return GetVecDouble("ANN_vtx");
    }
  }





  //Stuff to do with weighting

  double GetANNSegmentsZPosWeighted() const
  {
    if (!hasMLPred()) return -999;
    double zPosFromSegment0 = getZPosFromSegment(GetVecElemInt("ANN_segments", 0));
    double zPosFromSegment1 = getZPosFromSegment(GetVecElemInt("ANN_segments", 1));
    double zPosFromSegment2 = getZPosFromSegment(GetVecElemInt("ANN_segments", 2));
    double conf0 = GetVecElem("ANN_plane_probs", 0);
    double conf1 = GetVecElem("ANN_plane_probs", 1);
    double conf2 = GetVecElem("ANN_plane_probs", 2);
    return (zPosFromSegment0*conf0+zPosFromSegment1*conf1+zPosFromSegment2*conf2)/(conf0+conf1+conf2);
  }


  //------------------------------------------------------------------
  //End Calibrator methods
  //------------------------------------------------------------------

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
  /* double GetTruthSegment() const //Dangerous since ANN_vtx doesnt always have a 3-vector in it, calling SetCoordinates without a 3 vector gives a segfault
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
  } */

    double GetTruthSegment() const {
        // 1. Check Targets (Highest Priority)
        // Using a small switch or array for immediate return.
        int tgtID = GetTruthTargetID();
        switch (tgtID) {
            case 1: return 9;
            case 2: return 18;
            case 3: return 27;
            case 4: return 45;
            case 5: return 50;
        }

        // 2. Retrieve Vertex Info
        // (tgtZ was removed as it was unused in the original code)
        int mdl = GetTruthVtxModule();
        int plane = GetTruthVtxPlane();

        // 3. Handle Special/Error Cases
        // FIXED BUG: Original code had "if(plane = -999)" which is an assignment.
        // Assuming you meant to check if this is the "Water" case.
        if (mdl == -999) {
            // Original logic implied returning 36 if mdl is -999,
            // regardless of plane (due to the assignment bug evaluating to true).
            return 36; 
        }

        // 4. Bounds Check
        // Valid mdl range in your logic is -5 to 119.
        // Valid plane is 1 or 2.
        if (plane < 1 || plane > 2) return -999;
        if (mdl < -5 || mdl > 119) return -999;

        // 5. Lookup Table (O(1) access)
        // The table handles the mapping. We offset mdl by +5 so -5 becomes index 0.
        // [index][0] = plane 1, [index][1] = plane 2
        static const int segment_lut[125][2] = {
            /* mdl -5 */ {1, 2},
            /* mdl -4 */ {3, 4},
            /* mdl -3 */ {5, 6},
            /* mdl -2 */ {7, 8},
            /* mdl -1 */ {9, 9},
            /* mdl  0 */ {10, 11},
            /* mdl  1 */ {12, 13},
            /* mdl  2 */ {14, 15},
            /* mdl  3 */ {16, 17},
            /* mdl  4 */ {18, -999}, // Plane 2 was missing in original -> -999
            /* mdl  5 */ {19, 20},
            /* mdl  6 */ {21, 22},
            /* mdl  7 */ {23, 24},
            /* mdl  8 */ {25, 26},
            /* mdl  9 */ {27, 27},
            /* mdl 10 */ {-999, -999}, // Missing in original
            /* mdl 11 */ {28, 29},
            /* mdl 12 */ {30, 31},
            /* mdl 13 */ {32, 33},
            /* mdl 14 */ {34, 35},
            /* mdl 15 */ {37, 38},
            /* mdl 16 */ {39, 40},
            /* mdl 17 */ {41, 42},
            /* mdl 18 */ {43, 44},
            /* mdl 19 */ {45, 45},
            /* mdl 20 */ {46, 47},
            /* mdl 21 */ {48, 49},
            /* mdl 22 */ {50, 50},
            /* mdl 23 */ {51, 52},
            /* mdl 24 */ {53, 54},
            /* mdl 25 */ {55, 56},
            /* mdl 26 */ {57, 58},
            /* mdl 27 */ {59, 60},
            /* mdl 28 */ {61, 62},
            /* mdl 29 */ {63, 64},
            /* mdl 30 */ {65, 66},
            /* mdl 31 */ {67, 68},
            /* mdl 32 */ {69, 70},
            /* mdl 33 */ {71, 72},
            /* mdl 34 */ {73, 74},
            /* mdl 35 */ {75, 76},
            /* mdl 36 */ {77, 78},
            /* mdl 37 */ {79, 80},
            /* mdl 38 */ {81, 82},
            /* mdl 39 */ {83, 84},
            /* mdl 40 */ {85, 86},
            /* mdl 41 */ {87, 88},
            /* mdl 42 */ {89, 90},
            /* mdl 43 */ {91, 92},
            /* mdl 44 */ {93, 94},
            /* mdl 45 */ {95, 96},
            /* mdl 46 */ {97, 98},
            /* mdl 47 */ {99, 100},
            /* mdl 48 */ {101, 102},
            /* mdl 49 */ {103, 104},
            /* mdl 50 */ {105, 106},
            /* mdl 51 */ {107, 108},
            /* mdl 52 */ {109, 110},
            /* mdl 53 */ {111, 112},
            /* mdl 54 */ {113, 114},
            /* mdl 55 */ {115, 116},
            /* mdl 56 */ {117, 118},
            /* mdl 57 */ {119, 120},
            /* mdl 58 */ {121, 122},
            /* mdl 59 */ {123, 124},
            /* mdl 60 */ {125, 126},
            /* mdl 61 */ {127, 128},
            /* mdl 62 */ {129, 130},
            /* mdl 63 */ {131, 132},
            /* mdl 64 */ {133, 134},
            /* mdl 65 */ {135, 136},
            /* mdl 66 */ {137, 138},
            /* mdl 67 */ {139, 140},
            /* mdl 68 */ {141, 142},
            /* mdl 69 */ {143, 144},
            /* mdl 70 */ {145, 146},
            /* mdl 71 */ {147, 148},
            /* mdl 72 */ {149, 150},
            /* mdl 73 */ {151, 152},
            /* mdl 74 */ {153, 154},
            /* mdl 75 */ {155, 156},
            /* mdl 76 */ {157, 158},
            /* mdl 77 */ {159, 160},
            /* mdl 78 */ {161, 162},
            /* mdl 79 */ {163, 164},
            /* mdl 80 */ {165, 166},
            /* mdl 81 */ {167, 168},
            /* mdl 82 */ {169, 170},
            /* mdl 83 */ {171, 172},
            /* mdl 84 */ {173, 174},
            /* mdl 85 */ {175, 176},
            /* mdl 86 */ {177, 178},
            /* mdl 87 */ {179, 180},
            /* mdl 88 */ {181, 182},
            /* mdl 89 */ {183, 184},
            /* mdl 90 */ {185, 186},
            /* mdl 91 */ {187, 188},
            /* mdl 92 */ {189, 190},
            /* mdl 93 */ {191, 192},
            /* mdl 94 */ {193, 194},
            /* mdl 95 */ {195, 195},
            /* mdl 96 */ {196, 196},
            /* mdl 97 */ {197, 197},
            /* mdl 98 */ {198, 198},
            /* mdl 99 */ {199, 199},
            /* mdl 100*/ {-999, 200}, // Plane 1 missing in original -> -999
            /* mdl 101*/ {-999, 201},
            /* mdl 102*/ {-999, 202},
            /* mdl 103*/ {-999, 203},
            /* mdl 104*/ {-999, 204},
            /* mdl 105*/ {-999, 205},
            /* mdl 106*/ {-999, 206},
            /* mdl 107*/ {-999, 207},
            /* mdl 108*/ {-999, 208},
            /* mdl 109*/ {-999, 209},
            /* mdl 110*/ {-999, 210},
            /* mdl 111*/ {-999, 211},
            /* mdl 112*/ {-999, 212},
            /* mdl 113*/ {-999, 213},
            /* mdl 114*/ {-999, 214},
            /* mdl 115*/ {-999, 2},
            /* mdl 116*/ {-999, 0},
            /* mdl 117*/ {-999, 1},
            /* mdl 118*/ {-999, 0},
            /* mdl 119*/ {-999, 2}
        };

        // Access table: mdl + 5 to handle negative indices. 
        // Plane - 1 to make it 0-indexed.
        return segment_lut[mdl + 5][plane - 1];
    }

  int GetMultiplicity() const {return GetInt("multiplicity");}

  //Still needed for some systematics to compile, but shouldn't be used for reweighting anymore.
  protected:
  #include "PlotUtils/WeightFunctions.h" // Get*Weight
};

#endif
