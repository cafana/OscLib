#include "OscLib/func/OscCalculatorPMNS.h"

namespace osc
{
  // --------------------------------------------------------------------------
  template <typename T>
  _OscCalculatorPMNS<T>::_OscCalculatorPMNS()
      : fMixDirty(true), fDmDirty(true), fPropDirty(true), fPrevAnti(0)
  {
  }

  // --------------------------------------------------------------------------
  template<typename T>
  _OscCalculatorPMNS<T>::~_OscCalculatorPMNS()
  {
  }

  // --------------------------------------------------------------------------
  template <typename T>
  _IOscCalculatorAdjustable<T>* _OscCalculatorPMNS<T>::Copy() const
  {
    return new _OscCalculatorPMNS<T>(*this);
  }

  // --------------------------------------------------------------------------
  template <typename T>
  T _OscCalculatorPMNS<T>::P(int flavBefore, int flavAfter, double E)
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
template class osc::_OscCalculatorPMNS<double>;

#ifndef DARWINBUILD
#include "Utilities/func/StanVar.h"
  template class osc::_OscCalculatorPMNS<stan::math::var>;
#endif
