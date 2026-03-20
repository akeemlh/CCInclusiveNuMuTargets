#ifndef STUDY_H
#define STUDY_H

//ROOT includes
#include "TDirectory.h"

class MichelEvent;
class CVUniverse;

class Study
{
  public:
    Study() {} //TODO: Any base constructor needed?

    void Selected(const CVUniverse& univ, const MichelEvent& evt, const double weight)
    {
      fillSelected(univ, evt, weight);
    }
    
    
    void SelectedSignal(const CVUniverse& univ, const MichelEvent& evt, const double weight)
    {
      fillSelectedSignal(univ, evt, weight);
    }

    void TruthSignal(const CVUniverse& univ, const MichelEvent& evt, const double weight)
    {
      fillTruthSignal(univ, evt, weight);
    }

    //Find Andrew if you need to know how to overload functions for drawing.
    //Only need this when you write a new Study.
    virtual void SaveOrDraw(TDirectory& outDir) = 0;

  private:
    using Hist = PlotUtils::HistWrapper<CVUniverse>;

    virtual void fillSelected(const CVUniverse& univ, const MichelEvent& evt, const double weight) = 0;
    virtual void fillSelectedSignal(const CVUniverse& univ, const MichelEvent& evt, const double weight) = 0;
    virtual void fillTruthSignal(const CVUniverse& univ, const MichelEvent& evt, const double weight) = 0;
};

#endif //STUDY_H
