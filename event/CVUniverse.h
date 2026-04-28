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
#include "util/MuonEnergyCorrection.h"

class CVUniverse : public PlotUtils::MinervaUniverse {

  public:
  #include "PlotUtils/MuonFunctions.h" // GetMinosEfficiencyWeight
  #include "PlotUtils/TruthFunctions.h" //Getq3True
  #include "PlotUtils/WeightFunctions.h" // Get*Weight
  
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
    if (cal!=nullptr)
    {
        calibrator = cal;
        calibratorSet = true;
    }

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

  //double GetANNPmu() const //MeV/c
  //{
  //  return sqrt(pow(GetVecElem("muon_corrected_p",0),2)+pow(GetVecElem("muon_corrected_p",1),2)+pow(GetVecElem("muon_corrected_p",2),2));
  //}

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

  //virtual double GetANNEmu() const //MeV
  //{
  //  return GetVecElem("muon_corrected_p", 3);
  //}

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

  double GetEnuTrueGeV() const //GeV
  {
    return GetEnuTrue()/1000.;
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

  virtual int IsMinosMatchMuon() const { return GetInt("isMinosMatchTrack") == 1; }

  virtual bool HasInteractionVertex() const {
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
  

  virtual bool InANNFiducial() const {
    return GetInt((GetAnaToolName() + "_ANN_in_fiducial_area").c_str()) == 1;
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

  virtual int GetHelicity() const { return GetInt((GetAnaToolName() + "_nuHelicity").c_str());}

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
    ROOT::Math::XYZVector ANNVtx = GetANNVertex(false);
    if (ANNVtx.X() == -999) return ANNVtx;
    double mu_px = GetDouble("MasterAnaDev_muon_Px");
    double mu_py = GetDouble("MasterAnaDev_muon_Py");
    double mu_pz = GetDouble("MasterAnaDev_muon_Pz");
    double mu_time = GetDouble("muon_trackVertexTime");

    double projx = ANNVtx.X() + ( mu_px / mu_pz ) * ( zPos - ANNVtx.Z() );
    double projy = ANNVtx.Y() + ( mu_py / mu_pz ) * ( zPos - ANNVtx.Z() );
    ROOT::Math::XYZVector outVtx = ROOT::Math::XYZVector( projx, projy, zPos );
    return outVtx;
}


double getZPosFromSegment(int segment) const
{
    static constexpr std::array<double, 215> zPositions = {
        -999.0,   // 0 (Unused, segment indices start at 1)
        4293.04,  // 1
        4313.68,  // 2
        4337.25,  // 3
        4357.9,   // 4
        4381.47,  // 5
        4402.11,  // 6
        4425.68,  // 7
        4446.33,  // 8
        4481.21,  // 9  (Target 1)
        4514.11,  // 10
        4534.76,  // 11
        4558.33,  // 12
        4578.97,  // 13
        4602.54,  // 14
        4623.19,  // 15
        4646.76,  // 16
        4667.4,   // 17
        4702.29,  // 18 (Target 2)
        4735.19,  // 19
        4755.83,  // 20
        4779.4,   // 21
        4800.05,  // 22
        4823.62,  // 23
        4844.26,  // 24
        4867.83,  // 25
        4888.48,  // 26
        4923.36,  // 27 (Target 3)
        5000.48,  // 28
        5021.12,  // 29
        5044.69,  // 30
        5065.34,  // 31
        5088.91,  // 32
        5109.55,  // 33
        5133.12,  // 34
        5153.77,  // 35
        5310.0,   // 36 (Water Target -- Find better number)
        5456.74,  // 37
        5477.38,  // 38
        5500.95,  // 39
        5521.6,   // 40
        5545.17,  // 41
        5565.81,  // 42
        5589.38,  // 43
        5610.02,  // 44
        5644.91,  // 45 (Target 4)
        5677.81,  // 46
        5698.45,  // 47
        5722.03,  // 48
        5742.67,  // 49
        5777.55,  // 50 (Target 5)
        5810.45,  // 51
        5831.1,   // 52
        5855.68,  // 53
        5876.33,  // 54
        5900.91,  // 55
        5921.56,  // 56
        5946.14,  // 57
        5966.79,  // 58
        5991.37,  // 59
        6012.01,  // 60
        6036.6,   // 61
        6057.24,  // 62
        6081.83,  // 63
        6102.47,  // 64
        6127.06,  // 65
        6147.7,   // 66
        6172.29,  // 67
        6192.93,  // 68
        6217.52,  // 69
        6238.16,  // 70
        6262.74,  // 71
        6283.39,  // 72
        6307.97,  // 73
        6328.62,  // 74
        6353.2,   // 75
        6373.85,  // 76
        6398.43,  // 77
        6419.08,  // 78
        6443.66,  // 79
        6464.3,   // 80
        6488.89,  // 81
        6509.53,  // 82
        6534.12,  // 83
        6554.76,  // 84
        6579.35,  // 85
        6599.99,  // 86
        6624.58,  // 87
        6645.22,  // 88
        6669.81,  // 89
        6690.45,  // 90
        6715.03,  // 91
        6735.68,  // 92
        6760.26,  // 93
        6780.91,  // 94
        6805.49,  // 95
        6826.14,  // 96
        6850.72,  // 97
        6871.37,  // 98
        6895.95,  // 99
        6916.59,  // 100
        6941.18,  // 101
        6961.82,  // 102
        6986.41,  // 103
        7007.05,  // 104
        7031.64,  // 105
        7052.28,  // 106
        7076.87,  // 107
        7097.51,  // 108
        7122.1,   // 109
        7142.74,  // 110
        7167.32,  // 111
        7187.97,  // 112
        7212.55,  // 113
        7233.2,   // 114
        7257.78,  // 115
        7278.43,  // 116
        7303.01,  // 117
        7323.66,  // 118
        7348.24,  // 119
        7368.88,  // 120
        7393.47,  // 121
        7414.11,  // 122
        7438.7,   // 123
        7459.34,  // 124
        7483.93,  // 125
        7504.57,  // 126
        7529.16,  // 127
        7549.8,   // 128
        7574.39,  // 129
        7595.03,  // 130
        7619.61,  // 131
        7640.26,  // 132
        7664.84,  // 133
        7685.49,  // 134
        7710.07,  // 135
        7730.72,  // 136
        7755.3,   // 137
        7775.95,  // 138
        7800.53,  // 139
        7821.17,  // 140
        7845.76,  // 141
        7866.4,   // 142
        7890.99,  // 143
        7911.63,  // 144
        7936.22,  // 145
        7956.86,  // 146
        7981.45,  // 147
        8002.09,  // 148
        8026.68,  // 149
        8047.32,  // 150
        8071.9,   // 151
        8092.55,  // 152
        8117.13,  // 153
        8137.78,  // 154
        8162.36,  // 155
        8183.01,  // 156
        8207.59,  // 157
        8228.24,  // 158
        8252.82,  // 159
        8273.46,  // 160
        8298.05,  // 161
        8318.69,  // 162
        8343.28,  // 163
        8363.92,  // 164
        8388.51,  // 165
        8409.15,  // 166
        8433.74,  // 167
        8454.38,  // 168
        8478.97,  // 169
        8499.61,  // 170
        8524.19,  // 171
        8544.84,  // 172
        8569.42,  // 173
        8590.07,  // 174
        8614.65,  // 175
        8635.3,   // 176
        8659.46,  // 177
        8680.1,   // 178
        8704.26,  // 179
        8724.9,   // 180
        8749.06,  // 181
        8769.71,  // 182
        8793.86,  // 183
        8814.51,  // 184
        8838.67,  // 185
        8859.31,  // 186
        8883.47,  // 187
        8904.11,  // 188
        8928.27,  // 189
        8948.92,  // 190
        8973.08,  // 191
        8993.72,  // 192
        9017.88,  // 193
        9038.52,  // 194
        9088.08,  // 195
        9135.41,  // 196
        9182.75,  // 197
        9230.08,  // 198
        9277.41,  // 199
        9324.74,  // 200
        9372.08,  // 201
        9419.41,  // 202
        9466.74,  // 203
        9514.07,  // 204
        9561.41,  // 205
        9608.74,  // 206
        9656.07,  // 207
        9703.4,   // 208
        9750.74,  // 209
        9798.07,  // 210
        9845.4,   // 211
        9892.73,  // 212
        9940.07,  // 213
        9987.4    // 214
    };

    if (segment >= 1 && segment < zPositions.size()) {
        return zPositions[segment];
    }

    std::cerr << "Warning: getZPosFromSegment queried with out-of-bounds segment: " << segment << "\n";
    return -999.0; 
}



int getNearestSegmentToZPos(double zpos) const  //A very rough way of getting the segment number from a given z vertex
{
    //The centre of the scintillator planes, from (https://cdcvs.fnal.gov/redmine/projects/minerva/wiki/Z_Positions_of_Planes_in_the_Full_MINERvA_Detector)
    static const std::vector<double> plane_centres = {4293.04,
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
struct SegmentInfo {
    int vtx_module;
    int vtx_plane;
    int targetID;
};
int GetTargetFromSegment( int segment, int& vtx_module, int& vtx_plane ) const
{

    static const std::array<SegmentInfo, 215> segmentMap = {{
        /* 00 - 09 */ {-5, 0, -1}, {-5, 1, -1}, {-5, 2, -1}, {-4, 1, -1}, {-4, 2, -1}, {-3, 1, -1}, {-3, 2, -1}, {-2, 1, -1}, {-2, 2, -1}, {-1, 1,  1},
        /* 10 - 19 */ { 0, 1, -1}, { 0, 2, -1}, { 1, 1, -1}, { 1, 2, -1}, { 2, 1, -1}, { 2, 2, -1}, { 3, 1, -1}, { 3, 2, -1}, { 4, 1,  2}, { 5, 1, -1},
        /* 20 - 29 */ { 5, 2, -1}, { 6, 1, -1}, { 6, 2, -1}, { 7, 1, -1}, { 7, 2, -1}, { 8, 1, -1}, { 8, 2, -1}, { 9, 2,  3}, {11, 1, -1}, {11, 2, -1},
        /* 30 - 39 */ {12, 1, -1}, {12, 2, -1}, {13, 1, -1}, {13, 2, -1}, {14, 1, -1}, {14, 2, -1}, {-999, -999, 6}, {15, 1, -1}, {15, 2, -1}, {16, 1, -1},
        /* 40 - 49 */ {16, 2, -1}, {17, 1, -1}, {17, 2, -1}, {18, 1, -1}, {18, 2, -1}, {19, 1,  4}, {20, 1, -1}, {20, 2, -1}, {21, 1, -1}, {21, 2, -1},
        /* 50 - 59 */ {22, 1,  5}, {23, 1, -1}, {23, 2, -1}, {24, 1, -1}, {24, 2, -1}, {25, 1, -1}, {25, 2, -1}, {26, 1, -1}, {26, 2, -1}, {27, 1, -1},
        /* 60 - 69 */ {27, 2, -1}, {28, 1, -1}, {28, 2, -1}, {29, 1, -1}, {29, 2, -1}, {30, 1, -1}, {30, 2, -1}, {31, 1, -1}, {31, 2, -1}, {32, 1, -1},
        /* 70 - 79 */ {32, 2, -1}, {33, 1, -1}, {33, 2, -1}, {34, 1, -1}, {34, 2, -1}, {35, 1, -1}, {35, 2, -1}, {36, 1, -1}, {36, 2, -1}, {37, 1, -1},
        /* 80 - 89 */ {37, 2, -1}, {38, 1, -1}, {38, 2, -1}, {39, 1, -1}, {39, 2, -1}, {40, 1, -1}, {40, 2, -1}, {41, 1, -1}, {41, 2, -1}, {42, 1, -1},
        /* 90 - 99 */ {42, 2, -1}, {43, 1, -1}, {43, 2, -1}, {44, 1, -1}, {44, 2, -1}, {45, 1, -1}, {45, 2, -1}, {46, 1, -1}, {46, 2, -1}, {47, 1, -1},
        /* 100 - 109 */ {47, 2, -1}, {48, 1, -1}, {48, 2, -1}, {49, 1, -1}, {49, 2, -1}, {50, 1, -1}, {50, 2, -1}, {51, 1, -1}, {51, 2, -1}, {52, 1, -1},
        /* 110 - 119 */ {52, 2, -1}, {53, 1, -1}, {53, 2, -1}, {54, 1, -1}, {54, 2, -1}, {55, 1, -1}, {55, 2, -1}, {56, 1, -1}, {56, 2, -1}, {57, 1, -1},
        /* 120 - 129 */ {57, 2, -1}, {58, 1, -1}, {58, 2, -1}, {59, 1, -1}, {59, 2, -1}, {60, 1, -1}, {60, 2, -1}, {61, 1, -1}, {61, 2, -1}, {62, 1, -1},
        /* 130 - 139 */ {62, 2, -1}, {63, 1, -1}, {63, 2, -1}, {64, 1, -1}, {64, 2, -1}, {65, 1, -1}, {65, 2, -1}, {66, 1, -1}, {66, 2, -1}, {67, 1, -1},
        /* 140 - 149 */ {67, 2, -1}, {68, 1, -1}, {68, 2, -1}, {69, 1, -1}, {69, 2, -1}, {70, 1, -1}, {70, 2, -1}, {71, 1, -1}, {71, 2, -1}, {72, 1, -1},
        /* 150 - 159 */ {72, 2, -1}, {73, 1, -1}, {73, 2, -1}, {74, 1, -1}, {74, 2, -1}, {75, 1, -1}, {75, 2, -1}, {76, 1, -1}, {76, 2, -1}, {77, 1, -1},
        /* 160 - 169 */ {77, 2, -1}, {78, 1, -1}, {78, 2, -1}, {79, 1, -1}, {79, 2, -1}, {80, 1, -1}, {80, 2, -1}, {81, 1, -1}, {81, 2, -1}, {82, 1, -1},
        /* 170 - 179 */ {82, 2, -1}, {83, 1, -1}, {83, 2, -1}, {84, 1, -1}, {84, 2, -1}, {85, 1, -1}, {85, 2, -1}, {86, 1, -1}, {86, 2, -1}, {87, 1, -1},
        /* 180 - 189 */ {87, 2, -1}, {88, 1, -1}, {88, 2, -1}, {89, 1, -1}, {89, 2, -1}, {90, 1, -1}, {90, 2, -1}, {91, 1, -1}, {91, 2, -1}, {92, 1, -1},
        /* 190 - 199 */ {92, 2, -1}, {93, 1, -1}, {93, 2, -1}, {94, 1, -1}, {94, 2, -1}, {95, 2, -1}, {96, 2, -1}, {97, 2, -1}, {98, 2, -1}, {99, 2, -1},
        /* 200 - 209 */ {100, 2, -1}, {101, 2, -1}, {102, 2, -1}, {103, 2, -1}, {104, 2, -1}, {105, 2, -1}, {106, 2, -1}, {107, 2, -1}, {108, 2, -1}, {109, 2, -1},
        /* 210 - 214 */ {110, 2, -1}, {111, 2, -1}, {112, 2, -1}, {113, 2, -1}, {114, 2, -1}
    }};

    if (segment >= 0 && segment < segmentMap.size()) {
        vtx_module = segmentMap[segment].vtx_module;
        vtx_plane  = segmentMap[segment].vtx_plane;
        return segmentMap[segment].targetID;
    }

    std::cerr << "Can't find a module/plane for given segment " << segment << std::endl;
    vtx_module = -999;
    vtx_plane  = -999;
    return -1;
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

// Following methods taken from Everardo's code:


// Get TRUE Hadron Quantities (always MeV, radians)
// NOTE: These 'truth_pi_*' containers all have 20 elements
// Their non-default elements correspond to all of and only the CHARGED
// (N.B. + & -) pions from the event's TG4Trajectory primary
// trajectories.
// Pion indices here do not correspond to reco hadron indices.
// Be careful with these functions. Make sure you understand the input
// and output.
// input: truth index. output: may be pip or pim.
double GetTpiTrue(int idx) const {
  double t_pi_E = GetVecElem("truth_pi_E", idx);
  if (t_pi_E == -1.) {
    std::cerr << "CVU::GetTpiTrue: Default energy.\n"
                 "Tried to access truth pion energy for a nonexistent "
                 "truth pion trajectory.\n";
    throw t_pi_E;
  }
  double chargedPiMass = 139.569;
  return t_pi_E - chargedPiMass;
}


int GetNChargedPionsTrue() const {
  return GetInt("truth_N_pip") + GetInt("truth_N_pim");
}

// input: truth index. output: may be pip or pim.
int GetPiChargeTrue(int idx) const {
  int t_pi_charge = GetVecElem("truth_pi_charge", idx);
  if (t_pi_charge == 0) {
    std::cerr << "CVU::GetPiChargeTrue: Default charge.\n"
                 "Tried to access truth pion charge for a nonexistent "
                 "truth pion trajectory.\n";
    throw t_pi_charge;
  }
  return t_pi_charge;
}

// !!!!!! Vector filled with pip and pim !!!!!!!
std::vector<double> GetTpiTrueVec() const {
  std::vector<double> ret;
  const int n_true_pions = GetNChargedPionsTrue();
  for (int idx = 0; idx < n_true_pions; ++idx) {
    ret.push_back(GetTpiTrue(idx));
  }
  return ret;
}
//==============================================================================
// Get and set pion candidates
//==============================================================================
int GetHighestEnergyTruePionIndex() const {
  std::vector<double> tpi_vec = GetTpiTrueVec();  // pip and pim
  const int n_true_pions = GetNChargedPionsTrue();
  // if (n_true_pions != tpi_vec.size()) {
  //  std::cerr << "GetHighestEnergyTruePionIndex Error: True pion size
  //  mismatch\n"; std::exit(1);
  //}

  int reigning_idx = -1;
  double reigning_tpi = 0;
  for (int idx = 0; idx < n_true_pions; ++idx) {
    if (tpi_vec[idx] > reigning_tpi && GetPiChargeTrue(idx) > 0.) {
      reigning_idx = idx;
      reigning_tpi = tpi_vec[idx];
    }
  }

  // This solution below doesn't account for pion charge.
  // see cache/AlgTests for some better attempts at std::alg/lambda solutions
  // TruePionIdx reigning_idx = distance(tpi_vec, max_element(tpi_vec.begin(),
  // tpi_vec.end()));

  return reigning_idx;
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

// from MINERvA Document 33762-v3
double GetAsitsRobFineFluxWeight() const
{
    // Implements updated parameterization of Rob Fine's Thesis-era flux correction
    // to the nu+e flux. This does not alter the uncertainties.
    int inNuPdg = GetTruthNuPDG();
    
    if (inNuPdg != 14 && inNuPdg != -14) return 1.0; // Apply only to muon neutrinos and muon antineutrinos

    double inEnu = GetEnuTrueGeV();
    double ratio = 1.0;

    if (inEnu < 4.5)
    {
        ratio = 1.095;
    }
    else if (inEnu < 12.78)
    {
        // Polynomial shape from 4.5 to 12.78 GeV
        double enupoly[6] = {
            -4.11608956e+00, 3.66971952e+00, -9.73158862e-01,
            1.21054305e-01, -7.13034149e-03, 1.60919755e-04
        };
        double tempratio = enupoly[0];
        double powenu = inEnu;
        for (int i = 1; i < 6; ++i)
        {
            tempratio += enupoly[i] * powenu;
            powenu *= inEnu;
        }
        ratio = tempratio;
    }
    else if (inEnu < 26.0)
    {
        ratio = 1.17;
    }
    else if (inEnu < 29.0)
    {
        ratio = 1.17 + (1.09 - 1.17) * (inEnu - 26.0) / (29.0 - 26.0);
    }
    else if (inEnu < 48.0)
    {
        ratio = 1.09;
    }
    else if (inEnu < 80.0)
    {
        ratio = 1.09 + (2.5 - 1.09) * (inEnu - 48.0) / (80.0 - 48.0);
    }
    else if (inEnu < 100.0)
    {
        ratio = 2.5 + (1.5 - 2.5) * (inEnu - 80.0) / (100.0 - 80.0);
    }
    else
    {
        ratio = 1.5;
    }

    return ratio;
}


    virtual double GetCOHPionWeight() const {
        double weight = 1.0;
        if(GetInt("mc_intType") == 4)
        {
            //int npi = GetTrueNPionsinEvent();
            //if (npi == 0) return 1.0;
            double angle = GetTrueAngleHighTpi();//*180./M_PI; //this is now in degrees
            //Angle is already in degrees GetTrueAngleHighTpi()
            //std::cout << "Printing Angle for COH event " << angle << std::endl;
            double epi = GetTrueHighEpi()/1000.; // This is supposed to be the energy of the pion!!!! 
            if (epi < 0 | angle < 0) return 1.0;
            else {
                weight *= GetCoherentPiWeight(angle, epi); //Inputs are in Degrees and GeV
                    //std::cout << "Printing COHerent weight for COH event for angle: " << angle << " degrees and KE : " << KE << " GeV. And weight: " << weight << std::endl;
                return weight;
            }
        }
        return 1.0;
    }

    double GetEhadTrue() const { return GetEnuTrue() - GetElepTrue(); }


    double CalcWexp(const double Q2, const double Ehad) const {
        double mass_hadron = (MinervaUnits::M_p + MinervaUnits::M_n) / 2; //Original mass for Everardos weighter

        //double za = GetANNTargetEffectiveZARatio(); //implementing an effective nucleon weight from actual target p/n ratio
        //if (za <= 0.0) za = 0.5; //Fallback
        //double mass_hadron = (za * MinervaUnits::M_p + (1.0 - za) * MinervaUnits::M_n);
        double W = pow(mass_hadron, 2.0) - Q2 + 2.0 * (mass_hadron)*Ehad;
        W = W > 0 ? sqrt(W) : 0.0;
        return W;
    }

    virtual double GetWexpTrue() const {
        return CalcWexp(GetQ2True(), GetEhadTrue());
    }

    virtual double GetCCOnePiWeight() const {
        double weight = 1.0;
        double wExpTrue = GetWexpTrue();
        int current = GetCurrent();
        int nuPDG = GetTruthNuPDG();
        int nPosPi = GetInt("truth_N_pip");
        int nPos0 = GetInt("truth_N_pi0");
        //int nNegPi = GetInt("truth_N_pin");
        int nK0 = GetInt("truth_N_K0");
        int nKCharged = GetInt("truth_N_chargedK");

        if (wExpTrue  > 1.4 && current == 1 && nuPDG == 14 && nPosPi == 1 && nPos0 == 0 && nK0 == 0 && nKCharged == 0)
        {
            int idx = (int)GetHighestEnergyTruePionIndex();
            weight *= GetChargedPionTuneWeight(GetTpiTrue(idx), GetQ2True());
        }
        return weight;
    }


    virtual double GetDiffractiveWeight() const { return GetDiffractiveUncWeight(); } //Redundant with below - could remove/refactor away

    double GetERecoilTrue() const
    {
        double Enu = GetEnuTrue();
        double Emu = GetElepTrue();
        double Erecoil = Enu-Emu;
        return Erecoil;
    }

    double GetBjorkenX() const
    {
        double Emu = GetANNEmuGeV();
        double mu_massGev = MinervaUnits::M_mu/1000;
        double qsquared =  2.0 * (GetANNRecoilEGeV()+Emu) * (Emu - GetANNPmuGeV() * cos(GetThetamu())) - (mu_massGev*mu_massGev);

        double za = GetANNTargetEffectiveZARatio();

        if (za <= 0.0) za = 0.5;

        //double M = (M_p + M_n) /(2.0*1000); //Not exact. To do - weight it based on target material p n ratio!!!!!!
        double M = (za * MinervaUnits::M_p + (1.0 - za) * MinervaUnits::M_n) / 1000.0;
        double r = GetANNRecoilEGeV();
        double xbj = xBjorken(qsquared, M, r);
        //double q2 = GetInclQ2Reco();
        //double r2 = GetANNRecoilE();
        //double M2 = (M_p + M_n) /2.0;
        //double xbj2 = xBjorken(q2, M2, r2);
        //std::cout<<"Bjorken X calculation: xbj: " <<xbj <<" xbj2 " << xbj2<<" qsquared: "<< qsquared << " M: " << M << " r: " << r << " Emu "<<Emu<< " GetANNRecoilEGeV() " <<GetANNRecoilEGeV() <<" GetANNPmuGeV() " <<GetANNPmuGeV() << " GetThetamu() " << GetThetamu()<< " mu_massGev " << mu_massGev<< std::endl;
        return xbj;
    }

    double GetCalcBjorkenXTrue() const {//Get calculated true Bjorken X
        double nuclMass = M_nucleon;
        int struckNucl = GetTargetNucleon();
        if (struckNucl == PDG_n){
            nuclMass=MinervaUnits::M_n;
        }
        else if (struckNucl == PDG_p){
            nuclMass=MinervaUnits::M_p;
        }
        double Enu = GetEnuTrue();
        double Emu = GetElepTrue();
        double Erecoil = Enu-Emu;
        double thetaMu = GetThetalepTrue();
        double Q2 = GetQ2True();
        
        double xbj = xBjorken(Q2, nuclMass, Erecoil);
        return xbj;
    }

    double GetBjorkenY() const
    {
        return GetANNRecoilE()/GetANNEnu();
    }

    double GetCalcBjorkenYTrue() const
    {
        double Enu = GetEnuTrue();
        double Emu = GetElepTrue();
        double Erecoil = Enu-Emu;
        return Erecoil/Enu;
    }

    virtual double GetTrueLowestTpiEvent() const {
        double tpi = 999999.;
        int pdgsize = GetInt("mc_nFSPart");
        for (int i = 0; i< pdgsize; i++)
            {
                int pdg = GetVecElem("mc_FSPartPDG", i);
                if (pdg != 211) continue;
            double energy = GetVecElem("mc_FSPartE", i);
            double momentumx = GetVecElem("mc_FSPartPx", i);
            double momentumy = GetVecElem("mc_FSPartPy", i);
            double momentumz = GetVecElem("mc_FSPartPz", i);
                double pionmomentum = TMath::Sqrt(pow(momentumx, 2) + pow(momentumy,2)+pow(momentumz,2));
            double pionmass = TMath::Sqrt(pow(energy, 2) - pow(pionmomentum, 2));  
            double KE = energy - pionmass;
            if (tpi > KE) tpi = KE;     
            }
        ///if(tpi < 99999.) std::cout << "Lowest Energy Primary Pion KE is " << tpi << std::endl;
        return tpi;
    } 

    virtual int GetTrueNPionsinEvent() const {
        int npion = 0;
        int pdgsize = GetInt("mc_nFSPart");
        for (int i = 0; i< pdgsize; i++)
        {
            int pdg = GetVecElem("mc_FSPartPDG", i);
            if (pdg == 211) npion++;
        }
        return npion;

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

  int GetTruthTargetCode() const {
    return GetInt("truth_target_code");
    //int ID = GetTruthTargetID();
    //int TgtZ = GetMCTargetZ();
    //if (ID != 0 && TgtZ != 0) return (ID*1000 + TgtZ);
  }

  int GetTruthVtxPlane() const {return GetInt("truth_vtx_plane");}
  int GetMADVtxPlane() const {return GetInt("MasterAnaDev_vtx_plane");}

  int GetTruthVtxModule() const {return GetInt("truth_vtx_module");}
  int GetMADVtxModule() const {return GetInt("MasterAnaDev_vtx_module");}



    //Correction to muon energy from shifting the vertex position by the ANN calibration
    double getMuonEnergyCorrection() const
    {
        double Emu = GetANNEmu(false); //Uncorrected EMu as input
        //std::cout<<"Original GetANNEmu(): " << GetANNEmu() <<std::endl;
        if ((Emu < 10.) || (Emu > 1e8)) return 0; //Outside of range that we can correct.
        int uncalibratedSegment = GetANNSegment(false);
        int calibratedSegment = GetANNSegment();
        int nSegmentsTransitted = uncalibratedSegment - calibratedSegment;
        if (nSegmentsTransitted == 0) return 0.0;
        int nTrackerSegmentsTransitted = nSegmentsTransitted; //Will need to adjust this when adding nuclear targets

        //Nuclear target transits
        //Very approximate because I'm in a rush and don't have time to do the math.
        //To be robust we need to caluclate intersections. It's not impossible that in the multimaterial targets an angled track
        //could partially transit multiple materials or that due to the track angle the material it transits doesn't match the vertex shadow, like we're doing here
        ROOT::Math::XYZVector vtx = GetANNVertex();
        PlotUtils::TargetUtils targetInfo;
        double distIronTransitted = 0;
        double distCarbonTransitted = 0;
        double distLeadTransitted = 0;
        double distWaterTransitted = 0; //This is very sketchy since the water isn't constant thickness all over :(
        if (((calibratedSegment >=9) && (uncalibratedSegment <= 9)) || ((calibratedSegment <= 9) && (uncalibratedSegment >=9))) //If this transits target 1
        {
            double scale = ((calibratedSegment == 9) || (uncalibratedSegment == 9)) ? 0.5 : 1.0; //If it shifts to inside this target, since we're taking the z pos
            //as being in the centre it only actually passes through half the material 
            if (targetInfo.InIron1VolMC(vtx.X(), vtx.Y(), targetInfo.GetTarget1CenterZMC()) ) distIronTransitted+=scale*PlotUtils::TargetProp::Thickness::Tgt1::Fe;
            else if (targetInfo.InLead1VolMC(vtx.X(), vtx.Y(), targetInfo.GetTarget1CenterZMC()) ) distLeadTransitted+=scale*PlotUtils::TargetProp::Thickness::Tgt1::Pb;
            nTrackerSegmentsTransitted--; //Remove this segment since it's a target not tracker segment
        }
        if (((calibratedSegment >=18) && (uncalibratedSegment <= 18)) || ((calibratedSegment <= 18) && (uncalibratedSegment >=18))) //If this transits target 1
        {
            double scale = ((calibratedSegment == 9) || (uncalibratedSegment == 9)) ? 0.5 : 1.0; //If it shifts to inside this target, since we're taking the z pos
            //as being in the centre it only actually passes through half the material 
            if (targetInfo.InIron2VolMC(vtx.X(), vtx.Y(), targetInfo.GetTarget2CenterZMC()) ) distIronTransitted+=scale*PlotUtils::TargetProp::Thickness::Tgt2::Fe;
            else if (targetInfo.InLead2VolMC(vtx.X(), vtx.Y(), targetInfo.GetTarget2CenterZMC()) ) distLeadTransitted+=scale*PlotUtils::TargetProp::Thickness::Tgt2::Pb;
            nTrackerSegmentsTransitted--; //Remove this segment since it's a target not tracker segment
        }
        if (((calibratedSegment >=27) && (uncalibratedSegment <= 27)) || ((calibratedSegment <= 27) && (uncalibratedSegment >=27))) //If this transits target 1
        {
            double scale = ((calibratedSegment == 9) || (uncalibratedSegment == 9)) ? 0.5 : 1.0; //If it shifts to inside this target, since we're taking the z pos
            //as being in the centre it only actually passes through half the material 
            if (targetInfo.InIron3VolMC(vtx.X(), vtx.Y(), targetInfo.GetTarget3CenterZMC()) ) distIronTransitted+=scale*PlotUtils::TargetProp::Thickness::Tgt3::Fe;
            else if (targetInfo.InLead3VolMC(vtx.X(), vtx.Y(), targetInfo.GetTarget3CenterZMC()) ) distLeadTransitted+=scale*PlotUtils::TargetProp::Thickness::Tgt3::Pb;
            nTrackerSegmentsTransitted--; //Remove this segment since it's a target not tracker segment
        }
        if (((calibratedSegment >=27) && (uncalibratedSegment <= 27)) || ((calibratedSegment <= 27) && (uncalibratedSegment >=27))) //If this transits target 1
        {
            double scale = ((calibratedSegment == 9) || (uncalibratedSegment == 9)) ? 0.5 : 1.0; //If it shifts to inside this target, since we're taking the z pos
            //as being in the centre it only actually passes through half the material 
            if (targetInfo.InWaterTargetVolMC(vtx.X(), vtx.Y(), targetInfo.GetTarget3CenterZMC()) ) distWaterTransitted+=scale*(PlotUtils::TargetProp::WaterTarget::Back - PlotUtils::TargetProp::WaterTarget::Face);
            nTrackerSegmentsTransitted--; //Remove this segment since it's a target not tracker segment
        }
        if (((calibratedSegment >=45) && (uncalibratedSegment <= 45)) || ((calibratedSegment <= 45) && (uncalibratedSegment >=45))) //If this transits target 1
        {
            double scale = ((calibratedSegment == 9) || (uncalibratedSegment == 9)) ? 0.5 : 1.0; //If it shifts to inside this target, since we're taking the z pos
            //as being in the centre it only actually passes through half the material 
            distLeadTransitted+=scale*PlotUtils::TargetProp::Thickness::Tgt4::Pb;
            nTrackerSegmentsTransitted--; //Remove this segment since it's a target not tracker segment
        }
        if (((calibratedSegment >=50) && (uncalibratedSegment <= 50)) || ((calibratedSegment <= 50) && (uncalibratedSegment >=50))) //If this transits target 1
        {
            double scale = ((calibratedSegment == 9) || (uncalibratedSegment == 9)) ? 0.5 : 1.0; //If it shifts to inside this target, since we're taking the z pos
            //as being in the centre it only actually passes through half the material 
            if (targetInfo.InIron5VolMC(vtx.X(), vtx.Y(), targetInfo.GetTarget5CenterZMC()) ) distIronTransitted+=scale*PlotUtils::TargetProp::Thickness::Tgt5::Fe;
            else if (targetInfo.InLead5VolMC(vtx.X(), vtx.Y(), targetInfo.GetTarget5CenterZMC()) ) distLeadTransitted+=scale*PlotUtils::TargetProp::Thickness::Tgt5::Pb;
            nTrackerSegmentsTransitted--; //Remove this segment since it's a target not tracker segment
        }
        double trackerSegmentThickness = 1.7; //cm
        double distanceTravelledInScint = trackerSegmentThickness * nTrackerSegmentsTransitted * cos(GetThetamu());
        double distanceTravelledInIron = distIronTransitted * cos(GetThetamu());
        double distanceTravelledInLead = distLeadTransitted * cos(GetThetamu());
        double distanceTravelledInCarbon = distCarbonTransitted * cos(GetThetamu());
        double distanceTravelledInWater = distWaterTransitted * cos(GetThetamu());

        double energyCorrection = MuonECorrection::muonEnergyCorrection(MuonECorrection::Material::Polystyrene, distanceTravelledInScint, Emu);
        energyCorrection += MuonECorrection::muonEnergyCorrection(MuonECorrection::Material::Iron, distanceTravelledInIron, Emu);
        energyCorrection += MuonECorrection::muonEnergyCorrection(MuonECorrection::Material::Lead, distanceTravelledInLead, Emu);
        energyCorrection += MuonECorrection::muonEnergyCorrection(MuonECorrection::Material::GraphiteC, distanceTravelledInCarbon, Emu);
        energyCorrection += MuonECorrection::muonEnergyCorrection(MuonECorrection::Material::Water, distanceTravelledInWater, Emu);
        //energyCorrection += muonEnergyCorrection(Material::TitaniumDioxide, distanceTravelledInWater, Emu);
        return energyCorrection;
    }



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
        for (int c = 0; c < 3; c++)
        {
            int seg = segments[c];
            int prob = probs[c];
            //The ANN model occasionally puts events in segment 215 (doesn't actually exist) as an overflow. Our calibrator doesn't know how to handle overflows so just putting it in plane 214, which has no real use in my analysis anyway. Not the best solution but luckily the ANN does this very very rarely (only a few times in the entire ME RHC and FHC playlists
            if ((seg>=215 || seg<0))
            { //Known edge case, what if there already is a result for segment 214? This doesn't happen even once in my sample but it's not inconceivable. In future may want to do something more robust.
                segments[c] = 214;
                probs[c] = prob;
            }
        }
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


  double GetANNTargetEffectiveZARatio(bool useCalibrator = true) const 
  {
    int tgt = GetANNTargetCode(useCalibrator);
    double z;
    double a;
    if (tgt == 6000)
    {
        z = 10.0; // 8 O + 2 H
        a  = 15.999 + (2 * 1.0078); // ~ 16 O + 2 H
    }
    else if (tgt > 1000)
    {
        z = tgt%1000;
        if (z == 6) a = 12.0107;
        else if (z == 26) a = 55.8450;
        else if (z == 82) a = 207.20;
        else return 0;
    }
    else //Assume plastic scintillator (CH) otherwise
    {
        z = 7.0;  // 6 C + 1 H
        a = 12.011 + 1.0078;  // ~ 12 C + 1 H
    }
    return (z/a);
  }



    double GetANNEmu(bool useCalibrator = true) const {
        double ANNEmu = GetVecElem("muon_corrected_p", 3);
        if (calibratorSet && useCalibrator)
        {
            double shift_val = getMuonEnergyCorrection();
            //std::cout<<"getMuonEnergyCorrection(): " << shift_val << std::endl;
            return shift_val+ANNEmu;
        }
        else
        {
            return ANNEmu;
        }
    }

    double GetANNPmu(bool useCalibrator = true) const {
        double pmu = sqrt(pow(GetVecElem("muon_corrected_p",0),2)+pow(GetVecElem("muon_corrected_p",1),2)+pow(GetVecElem("muon_corrected_p",2),2));
        if (calibratorSet && useCalibrator)
        {
            double p_sum_squared = std::pow(GetANNEmu(),2);
            double m_mu = MinervaUnits::M_mu;
            double p_squared = p_sum_squared - (m_mu * m_mu);
            double p = sqrt(p_squared);
            //std::cout<<"Corrected Pmu: " << p <<" original Pmu: " << pmu << std::endl;
            return sqrt(p_squared);
        }
        else
        {
            return pmu;
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

  //double GetRecoilEnergy() const { return ApplyCaloTuning(GetCalRecoilEnergy()) + GetNonCalRecoilEnergy();  //  LOOK AT PlotUtils/SystCalcs/RecoilEnergyFunctions.h
  //double GetEnu()           const {return GetEmu()+ GetRecoilEnergy();} 
  //double GetEnuGeV()           const {return GetEnu()/1000;} 
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
};

#endif
