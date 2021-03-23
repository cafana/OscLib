// n.b. Stan sets up some type traits that need to be loaded before Eigen is.
// Since Eigen gets dragged in via IOscCalc.h we have to get Stan set up before
// that is included.
#ifdef OSCLIB_STAN
#include "OscLib/Stan.h"
#endif

#include "OscLib/OscCalcPMNS.h"

namespace osc
{
  // --------------------------------------------------------------------------
  template <typename T>
  _OscCalcPMNS<T>::_OscCalcPMNS()
      : fMixDirty(true), fDmDirty(true), fPropDirty(true), fPrevAnti(0)
  {
  }

  // --------------------------------------------------------------------------
  template<typename T>
  _OscCalcPMNS<T>::~_OscCalcPMNS()
  {
  }

  // --------------------------------------------------------------------------
  template <typename T>
  _IOscCalcAdjustable<T>* _OscCalcPMNS<T>::Copy() const
  {
    return new _OscCalcPMNS<T>(*this);
  }

  // --------------------------------------------------------------------------
  template <typename T>
  T _OscCalcPMNS<T>::P(int flavBefore, int flavAfter, double E)
  {
    const int anti = (flavBefore > 0) ? +1 : -1;
    assert(flavAfter/anti > 0);
    if(anti != fPrevAnti) fPropDirty = true;

    int i = -1, j = -1;
    if(abs(flavBefore) == 12) i = 0;
    if(abs(flavBefore) == 14) i = 1;
    if(abs(flavBefore) == 16) i = 2;
    if(abs(flavAfter) == 12) j = 0;
    if(abs(flavAfter) == 14) j = 1;
    if(abs(flavAfter) == 16) j = 2;
    assert(i >= 0 && j >= 0);

    if(fMixDirty){
      fPMNS.SetMix(this->fTh12, this->fTh23, this->fTh13, this->fdCP);
      fMixDirty = false;
    }
    if(fDmDirty){
      fPMNS.SetDeltaMsqrs(this->fDmsq21, this->fDmsq32);
      fDmDirty = false;
    }

    if(fPropDirty || E != fPrevE){
      fPMNS.Reset();
      // Assume Z/A=0.5
      const double Ne = this->fRho/2;
      fPMNS.PropMatter(this->fL, E, Ne, anti);

      fPropDirty = false;
      fPrevE = E;
      fPrevAnti = anti;
    }

    return fPMNS.P(i, j);
  }
}

//---------------------------------------------------------------------------
template class osc::_OscCalcPMNS<double>;

#ifdef OSCLIB_STAN
template class osc::_OscCalcPMNS<stan::math::var>;
#endif
