#ifndef MUON_ENERGY_CORRECTION
#define MUON_ENERGY_CORRECTION

// =============================================================================
//
// Muon energy-loss correction using tabulated stopping powers from:
//   Groom, Mokhov & Striganov, "Muon Stopping Power and Range Tables
//   10 MeV - 100 TeV", Atomic Data and Nuclear Data Tables, Vol. 76,
//   No. 2, July 2001.
//
// Usage:
//   double dE = muonEnergyCorrection(Material::Polystyrene, delta_x_cm, T_MeV);
//
// Arguments:
//   material    : one of the Material enum values
//   delta_x_cm  : signed path-length correction [cm]
//                 +ve => muon travelled MORE (add energy to reco vertex)
//                 -ve => muon travelled LESS (subtract energy from reco vertex)
//   T_muon_MeV  : muon kinetic energy at the vertex [MeV]
//
// Returns:
//   delta_E_MeV : energy to ADD to reconstructed interaction energy [MeV]
//                 Apply as: E_corrected = E_reco + delta_E
//
// Notes:
//   - "Total" column used (ionisation + brems + pair prod + photonuclear).
//   - Interpolation is linear in log10(T), matching the table's log spacing.
//   - Valid range: 10 MeV <= T_muon_MeV <= 100 TeV for all materials.
//   - Track-angle correction (1/cos θ) must be applied by the caller if
//     delta_x_cm is a projected (beam-axis) displacement rather than true
//     path length.
// =============================================================================

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

namespace MuonECorrection {

// ---------------------------------------------------------------------------
// Supported materials
// ---------------------------------------------------------------------------
enum class Material {
    Polystyrene,    // (C6H5CHCH2)n,  rho = 1.060 g/cm^3
    GraphiteC,      // Carbon graphite, rho = 2.210 g/cm^3
    Water,          // Liquid H2O,    rho = 1.000 g/cm^3
    Iron,           // Fe,            rho = 7.874 g/cm^3
    Lead,           // Pb,            rho = 11.350 g/cm^3
    TitaniumDioxide // TiO2,          rho = 4.260 g/cm^3
};

// Human-readable name (useful for logging / error messages)
inline std::string materialName(Material m) {
    switch (m) {
        case Material::Polystyrene:     return "Polystyrene [(C6H5CHCH2)n]";
        case Material::GraphiteC:       return "Carbon graphite (C)";
        case Material::Water:           return "Liquid water (H2O)";
        case Material::Iron:            return "Iron (Fe)";
        case Material::Lead:            return "Lead (Pb)";
        case Material::TitaniumDioxide: return "Titanium dioxide (TiO2)";
    }
    return "Unknown";
}

// ---------------------------------------------------------------------------
// Internal table structure
// ---------------------------------------------------------------------------
struct MaterialData {
    double      rho;          // density [g/cm^3]
    int         N;            // number of table entries
    const double* T_MeV;     // kinetic energies [MeV]
    const double* dEdx_mass; // total mass stopping power [MeV cm^2/g]
};

// ---------------------------------------------------------------------------
// Table data — all values from Groom et al. 2001, "Total" column
// ---------------------------------------------------------------------------

// ---- Polystyrene (C6H5CHCH2)n, rho = 1.060 g/cm^3 -----------------------
static constexpr int N_PS = 45;
static constexpr double T_PS[N_PS] = {
    10.,14.,20.,30.,40.,80.,100.,140.,200.,300.,
    318.,400.,800.,1.e3,1.4e3,2.e3,3.e3,4.e3,8.e3,1.e4,
    1.4e4,2.e4,3.e4,4.e4,8.e4,1.e5,1.4e5,2.e5,3.e5,4.e5,
    8.e5,1.e6,1.18e6,1.4e6,2.e6,3.e6,4.e6,8.e6,1.e7,1.4e7,
    2.e7,3.e7,4.e7,8.e7,1.e8
};
static constexpr double dEdx_PS[N_PS] = {
    7.803,6.084,4.749,3.683,3.144,
    2.359,2.211,2.058,1.971,1.937,
    1.936,1.943,2.015,2.049,2.103,
    2.163,2.231,2.279,2.391,2.427,
    2.480,2.538,2.607,2.662,2.825,
    2.895,3.022,3.203,3.489,3.770,
    4.885,5.445,5.953,6.558,8.246,
    11.059,13.904,25.400,31.208,42.844,
    60.458,89.984,119.748,240.029,300.688
};

// ---- Carbon graphite (C), rho = 2.210 g/cm^3 ----------------------------
static constexpr int N_C = 45;
static constexpr double T_C[N_C] = {
    10.,14.,20.,30.,40.,80.,100.,140.,200.,300.,
    317.,400.,800.,1.e3,1.4e3,2.e3,3.e3,4.e3,8.e3,1.e4,
    1.4e4,2.e4,3.e4,4.e4,8.e4,1.e5,1.4e5,2.e5,3.e5,4.e5,
    8.e5,1.e6,1.057e6,1.4e6,2.e6,3.e6,4.e6,8.e6,1.e7,1.4e7,
    2.e7,3.e7,4.e7,8.e7,1.e8
};
static constexpr double dEdx_C[N_C] = {
    7.114,5.548,4.329,3.353,2.858,
    2.124,1.989,1.851,1.772,1.743,
    1.742,1.749,1.817,1.848,1.899,
    1.955,2.020,2.065,2.171,2.205,
    2.255,2.310,2.376,2.429,2.589,
    2.657,2.785,2.966,3.255,3.540,
    4.673,5.244,5.406,6.379,8.102,
    10.972,13.875,25.606,31.532,43.402,
    61.366,91.472,121.812,244.371,306.157
};

// ---- Liquid water (H2O), rho = 1.000 g/cm^3 -----------------------------
static constexpr int N_W = 45;
static constexpr double T_W[N_W] = {
    10.,14.,20.,30.,40.,80.,100.,140.,200.,300.,
    317.,400.,800.,1.e3,1.4e3,2.e3,3.e3,4.e3,8.e3,1.e4,
    1.4e4,2.e4,3.e4,4.e4,8.e4,1.e5,1.4e5,2.e5,3.e5,4.e5,
    8.e5,1.e6,1.029e6,1.4e6,2.e6,3.e6,4.e6,8.e6,1.e7,1.4e7,
    2.e7,3.e7,4.e7,8.e7,1.e8
};
static constexpr double dEdx_W[N_W] = {
    7.902,6.166,4.817,3.738,3.192,
    2.398,2.256,2.104,2.015,1.981,
    1.981,1.988,2.064,2.099,2.156,
    2.218,2.290,2.341,2.459,2.497,
    2.553,2.615,2.690,2.749,2.931,
    3.009,3.156,3.365,3.698,4.028,
    5.341,6.003,6.097,7.318,9.315,
    12.640,16.003,29.584,36.442,50.168,
    70.932,105.697,140.712,282.051,353.261
};

// ---- Iron (Fe), rho = 7.874 g/cm^3 --------------------------------------
static constexpr int N_FE = 45;
static constexpr double T_FE[N_FE] = {
    10.,14.,20.,30.,40.,80.,100.,140.,200.,274.,
    300.,400.,800.,1.e3,1.4e3,2.e3,3.e3,4.e3,8.e3,1.e4,
    1.4e4,2.e4,3.e4,4.e4,8.e4,1.e5,1.4e5,2.e5,3.e5,3.47e5,
    4.e5,8.e5,1.e6,1.4e6,2.e6,3.e6,4.e6,8.e6,1.e7,1.4e7,
    2.e7,3.e7,4.e7,8.e7,1.e8
};
static constexpr double dEdx_FE[N_FE] = {
    5.494,4.321,3.399,2.654,2.274,
    1.717,1.616,1.516,1.463,1.451,
    1.453,1.467,1.548,1.582,1.637,
    1.697,1.766,1.816,1.935,1.975,
    2.038,2.112,2.213,2.302,2.620,
    2.773,3.076,3.536,4.293,4.652,
    5.063,8.190,9.781,12.958,17.799,
    25.857,34.006,66.835,83.371,116.400,
    166.242,249.364,332.844,668.005,836.153
};

// ---- Lead (Pb), rho = 11.350 g/cm^3 -------------------------------------
static constexpr int N_PB = 45;
static constexpr double T_PB[N_PB] = {
    10.,14.,20.,30.,40.,80.,100.,140.,200.,226.,
    300.,400.,800.,1.e3,1.4e3,2.e3,3.e3,4.e3,8.e3,1.e4,
    1.4e4,2.e4,3.e4,4.e4,8.e4,1.e5,1.4e5,1.41e5,2.e5,3.e5,
    4.e5,8.e5,1.e6,1.4e6,2.e6,3.e6,4.e6,8.e6,1.e7,1.4e7,
    2.e7,3.e7,4.e7,8.e7,1.e8
};
static constexpr double dEdx_PB[N_PB] = {
    3.823,3.054,2.436,1.928,1.666,
    1.283,1.215,1.152,1.124,1.123,
    1.131,1.152,1.238,1.272,1.325,
    1.385,1.454,1.504,1.635,1.683,
    1.765,1.871,2.033,2.190,2.810,
    3.126,3.763,3.782,4.750,6.396,
    8.087,14.997,18.526,25.569,36.310,
    54.186,72.258,144.979,181.567,254.601,
    364.707,548.073,732.025,1469.266,1838.603
};

// ---- Titanium dioxide (TiO2), rho = 4.260 g/cm^3 ------------------------
static constexpr int N_TIO2 = 45;
static constexpr double T_TIO2[N_TIO2] = {
    10.,14.,20.,30.,40.,80.,100.,140.,200.,288.,
    300.,400.,800.,1.e3,1.4e3,2.e3,3.e3,4.e3,8.e3,1.e4,
    1.4e4,2.e4,3.e4,4.e4,8.e4,1.e5,1.4e5,2.e5,3.e5,4.e5,
    5.154e5,8.e5,1.e6,1.4e6,2.e6,3.e6,4.e6,8.e6,1.e7,1.4e7,
    2.e7,3.e7,4.e7,8.e7,1.e8
};
static constexpr double dEdx_TIO2[N_TIO2] = {
    6.055,4.746,3.723,2.902,2.486,
    1.867,1.755,1.641,1.580,1.562,
    1.562,1.574,1.650,1.683,1.735,
    1.794,1.861,1.908,2.020,2.057,
    2.114,2.179,2.262,2.333,2.574,
    2.686,2.904,3.228,3.758,4.294,
    4.917,6.457,7.555,9.745,13.080,
    18.632,24.246,46.885,58.296,81.109,
    115.562,173.067,230.866,463.300,580.060
};

// ---------------------------------------------------------------------------
// Dispatch table
// ---------------------------------------------------------------------------
inline MaterialData getMaterialData(Material m) {
    switch (m) {
        case Material::Polystyrene:
            return {1.060, N_PS,   T_PS,   dEdx_PS};
        case Material::GraphiteC:
            return {2.210, N_C,    T_C,    dEdx_C};
        case Material::Water:
            return {1.000, N_W,    T_W,    dEdx_W};
        case Material::Iron:
            return {7.874, N_FE,   T_FE,   dEdx_FE};
        case Material::Lead:
            return {11.350, N_PB,  T_PB,   dEdx_PB};
        case Material::TitaniumDioxide:
            return {4.260, N_TIO2, T_TIO2, dEdx_TIO2};
    }
    throw std::invalid_argument("muonEnergyCorrection: unknown material");
}

// ---------------------------------------------------------------------------
// Core interpolation — log-linear in T, matching the table's log spacing
// ---------------------------------------------------------------------------
inline double interpolateDeDx(const MaterialData& mat, double T_MeV) {
    const int    N     = mat.N;
    const double* T    = mat.T_MeV;
    const double* dEdx = mat.dEdx_mass;

    if (T_MeV <= T[0])     return dEdx[0];
    if (T_MeV >= T[N - 1]) return dEdx[N - 1];

    // Binary search for bracketing interval
    int lo = 0, hi = N - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        (T_MeV >= T[mid] ? lo : hi) = mid;
    }

    const double logT    = std::log10(T_MeV);
    const double logT_lo = std::log10(T[lo]);
    const double logT_hi = std::log10(T[hi]);
    const double span    = logT_hi - logT_lo;

    if (span < 1e-12) return dEdx[lo];

    const double frac = (logT - logT_lo) / span;
    return dEdx[lo] + frac * (dEdx[hi] - dEdx[lo]);
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Returns the mass stopping power [MeV cm^2/g] for a given material and
/// muon kinetic energy, interpolated from the Groom et al. 2001 tables.
inline double massStoppingPower(Material m, double T_muon_MeV) {
    if (T_muon_MeV < 10.0 || T_muon_MeV > 1.0e8)
        throw std::out_of_range(
            "massStoppingPower: T_muon_MeV of "+ std::to_string(T_muon_MeV)+ " is out of table range [10 MeV, 100 TeV]"
            " for material: " + materialName(m));
    return interpolateDeDx(getMaterialData(m), T_muon_MeV);
}

/// Returns the linear stopping power [MeV/cm] for a given material and
/// muon kinetic energy.
inline double linearStoppingPower(Material m, double T_muon_MeV) {
    const MaterialData mat = getMaterialData(m);
    if (T_muon_MeV < 10.0 || T_muon_MeV > 1.0e8)
        throw std::out_of_range(
            "linearStoppingPower: T_muon_MeV of "+ std::to_string(T_muon_MeV)+ " is out of table range [10 MeV, 100 TeV]"
            " for material: " + materialName(m));
    return interpolateDeDx(mat, T_muon_MeV) * mat.rho;
}

/// Compute the muon vertex energy correction.
///
/// @param m           Target material (Material enum)
/// @param delta_x_cm  Signed path-length correction [cm]
///                    +ve: muon went further → add energy to E_reco
///                    -ve: muon went less far → subtract energy from E_reco
/// @param T_muon_MeV  Muon kinetic energy at the vertex [MeV]
/// @return            Energy correction delta_E [MeV] to add to E_reco
///
/// NOTE: If delta_x_cm is a projected distance along the beam axis rather
///       than the true muon path length, the caller must divide by cos(theta)
///       before passing it here, where theta is the muon track angle w.r.t.
///       the beam axis.
inline double muonEnergyCorrection(Material m,
                                   double delta_x_cm,
                                   double T_muon_MeV)
{
    return linearStoppingPower(m, T_muon_MeV) * delta_x_cm;
}

}
#endif //MUON_ENERGY_CORRECTION