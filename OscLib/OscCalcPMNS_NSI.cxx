// n.b. Stan sets up some type traits that need to be loaded before Eigen is.
// Since Eigen gets dragged in via IOscCalc.h we have to get Stan set up before
// that is included.
#ifdef OSCLIB_STAN
#include "OscLib/Stan.h"
#endif

#include "OscLib/OscCalcPMNS_NSI.h"

#include "OscLib/Constants.h"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>

namespace
{
  template <typename U> double ValAsDouble(const U& x) { return x; }

#ifdef OSCLIB_STAN
  template <> double ValAsDouble<stan::math::var>(const stan::math::var& x) { return x.val(); }
#endif
}

namespace osc
{
  template <typename T>
  _OscCalcPMNS_NSI<T>::_OscCalcPMNS_NSI()
    : fMixDirty(true), fDmDirty(true), fPropDirty(true), fEpsDirty(true), fPrevAnti(0)
  {
  }
  
  //---------------------------------------------------------------------------
  // 2016-10-24 - Getting all oscillation parameters (Standard and NSI)
  template <typename T>
  std::vector<double> _OscCalcPMNS_NSI<T>::GetState() const
  {
    std::vector<double> state;
    state.push_back(this->fL);
    state.push_back(this->fRho);
    state.push_back(ValAsDouble(this->fDmsq21));
    state.push_back(ValAsDouble(this->fDmsq32));
    state.push_back(ValAsDouble(this->fTh12));
    state.push_back(ValAsDouble(this->fTh13));
    state.push_back(ValAsDouble(this->fTh23));
    state.push_back(ValAsDouble(this->fdCP));
    state.push_back(ValAsDouble(fEps_ee));
    state.push_back(ValAsDouble(fEps_emu));
    state.push_back(ValAsDouble(fEps_etau));
    state.push_back(ValAsDouble(fEps_mumu));
    state.push_back(ValAsDouble(fEps_mutau));
    state.push_back(ValAsDouble(fEps_tautau));
    state.push_back(ValAsDouble(fDelta_emu));
    state.push_back(ValAsDouble(fDelta_etau));
    state.push_back(ValAsDouble(fDelta_mutau));
    
    return state;
  }
  //---------------------------------------------------------------------------
  
  //---------------------------------------------------------------------------
  // 2016-10-24 - Setting all oscillation parameters (Standard and NSI) - or the State
  template <typename T>
  void _OscCalcPMNS_NSI<T>::SetState(std::vector<double> state)
  {
    int iState(0);
    fMixDirty = true ;
    SetL(state[iState++]);
    SetRho(state[iState++]);
    SetDmsq21(state[iState++]);
    SetDmsq32(state[iState++]);
    SetTh12(state[iState++]);
    SetTh13(state[iState++]);
    SetTh23(state[iState++]);
    SetdCP(state[iState++]);
    SetEps_ee(state[iState++]);
    SetEps_emu(state[iState++]);
    SetEps_etau(state[iState++]);
    SetEps_mumu(state[iState++]);
    SetEps_mutau(state[iState++]);
    SetEps_tautau(state[iState++]);
    SetDelta_emu(state[iState++]);
    SetDelta_etau(state[iState++]);
    SetDelta_mutau(state[iState++]);
  }
  //---------------------------------------------------------------------------

  template <typename T>
  _OscCalcPMNS_NSI<T>::~_OscCalcPMNS_NSI()
  {
  }

  template <typename T>
  T _OscCalcPMNS_NSI<T>::P(int flavBefore, int flavAfter, double E)
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
      fPMNS_NSI.SetMix(this->fTh12, this->fTh23, this->fTh13, this->fdCP);
      fMixDirty = false;
    }
    if(fDmDirty){
      fPMNS_NSI.SetDeltaMsqrs(this->fDmsq21, this->fDmsq32);
      fDmDirty = false;
    }
    if(fEpsDirty){
      fPMNS_NSI.SetNSI(fEps_ee,    fEps_emu,    fEps_etau,
                       fEps_mumu,  fEps_mutau,  fEps_tautau,
                       fDelta_emu, fDelta_etau, fDelta_mutau);
      fEpsDirty = false;
    }


    fPMNS_NSI.ResetToFlavour(i);
    const double Ne = this->fRho * constants::kZPerA;
    fPMNS_NSI.PropMatter(this->fL, E, Ne, anti);
    return fPMNS_NSI.P(j);
  }

  //---------------------------------------------------------------------------
  template <typename T>
  _IOscCalcAdjustable<T>* _OscCalcPMNS_NSI<T>::Copy() const
  {
    return new _OscCalcPMNS_NSI<T>(*this);
  }
  
  //---------------------------------------------------------------------------
  const OscCalcPMNS_NSI* DowncastToNSI(const IOscCalc* calc)
  {
    const OscCalcPMNS_NSI* calc_nsi = dynamic_cast<const OscCalcPMNS_NSI*>(calc);
    if(calc_nsi) return calc_nsi;
    else             std::cout << "Input calculator was not of type OscCalcPMNS_NSI." << std::endl;
    return nullptr; // If the cast failed, calc_nsi should be nullptr anyway (?)                    
  }

  //---------------------------------------------------------------------------                     
  OscCalcPMNS_NSI* DowncastToNSI(IOscCalc* calc)
  {
    OscCalcPMNS_NSI* calc_nsi = dynamic_cast<OscCalcPMNS_NSI*>(calc);
    if(calc_nsi) return calc_nsi;
    else             std::cout << "Input calculator was not of type OscCalcPMNS_NSI." << std::endl;
    return nullptr; // If the cast failed, calc_nsi should be nullptr anyway (?)                    
  }
  
  template class _OscCalcPMNS_NSI<double>;

#ifdef OSCLIB_STAN
  template class _OscCalcPMNS_NSI<stan::math::var>;
#endif

} // namespace
