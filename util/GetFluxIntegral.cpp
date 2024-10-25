//File: GetFluxIntegral.cpp
//Brief: Provides the integrated flux needed at the end of a differential cross section extraction.
//       Matches the binning of an input histogram.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//Includes from this package
#include "util/GetFluxIntegral.h"
#include "event/CVUniverse.h"

//PlotUtils includes
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/FluxReweighter.h"

class CVUniverse;

namespace PlotUtils
{
  class MnvH1D;
  class MnvH2D;
}

namespace util
{
  PlotUtils::MnvH1D* GetFluxIntegral(const CVUniverse& univ, PlotUtils::MnvH1D* templateHist, const double Emin /*GeV*/, const double Emax /*GeV*/)
  {
    const bool useMuonCorrelations = true;
    assert(!(useMuonCorrelations && (univ.GetAnalysisNuPDG() < 0)) && "Muon momentum correlations are not yet ready for ME antineutrino analyses!");

    auto& frw = PlotUtils::flux_reweighter(univ.GetPlaylist(), univ.GetAnalysisNuPDG(), univ.UseNuEConstraint(), univ.GetNFluxUniverses());
    return frw.GetIntegratedFluxReweighted(univ.GetAnalysisNuPDG(), templateHist, Emin, Emax, useMuonCorrelations);
  }
  PlotUtils::MnvH2D* GetFluxIntegral(const CVUniverse& univ, PlotUtils::MnvH2D* templateHist, const double Emin /*GeV*/, const double Emax /*GeV*/)
  {
    const bool useMuonCorrelations = true;
    assert(!(useMuonCorrelations && (univ.GetAnalysisNuPDG() < 0)) && "Muon momentum correlations are not yet ready for ME antineutrino analyses!");

    auto& frw = PlotUtils::flux_reweighter(univ.GetPlaylist(), univ.GetAnalysisNuPDG(), univ.UseNuEConstraint(), univ.GetNFluxUniverses());
    return frw.GetIntegratedFluxReweighted(univ.GetAnalysisNuPDG(), templateHist, Emin, Emax, useMuonCorrelations);
  }
}
