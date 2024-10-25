//File: GetFluxIntegral.h
//Brief: Provides the integrated flux needed at the end of a differential cross section extraction.
//       Matches the binning of an input histogram.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef UTIL_GETFLUXINTEGRAL_H
#define UTIL_GETFLUXINTEGRAL_H

class CVUniverse;

namespace PlotUtils
{
  class MnvH1D;
  class MnvH2D;
}

namespace util
{
  PlotUtils::MnvH1D* GetFluxIntegral(const CVUniverse& univ, PlotUtils::MnvH1D* templateHist, const double Emin = 0 /*GeV*/, const double Emax = 100 /*GeV*/);
  PlotUtils::MnvH2D* GetFluxIntegral(const CVUniverse& univ, PlotUtils::MnvH2D* templateHist, const double Emin = 0 /*GeV*/, const double Emax = 100 /*GeV*/);
}

#endif //UTIL_GETFLUXINTEGRAL_H
