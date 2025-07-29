// n.b. Stan sets up some type traits that need to be loaded before Eigen is.
// Since Eigen gets dragged in via IOscCalc.h we have to get Stan set up before
// that is included.
#ifdef OSCLIB_STAN
#include "OscLib/Stan.h"
#endif

#include "OscLib/OscCalcPMNSOpt.h"

#include "OscLib/Constants.h"

namespace osc
{
  //---------------------------------------------------------------------------
  template <typename T>
  _OscCalcPMNSOpt<T>::_OscCalcPMNSOpt()
      : fMixIdx(0), fDmIdx(0), fLRIdx(0)
  {
  }

  //---------------------------------------------------------------------------
  template <typename T>
  _OscCalcPMNSOpt<T>::_OscCalcPMNSOpt(const osc::_OscCalcPMNSOpt<T>& calc)
    : _IOscCalcAdjustable<T>(calc), fMixIdx(0), fDmIdx(0), fLRIdx(0)
  {
    // NB we don't copy any of the cache entries over
  }

  //---------------------------------------------------------------------------
  template <typename T>
  _OscCalcPMNSOpt<T>::~_OscCalcPMNSOpt()
  {
  }

  //---------------------------------------------------------------------------
  template <typename T>
  _IOscCalcAdjustable<T>* _OscCalcPMNSOpt<T>::Copy() const
  {
    return new _OscCalcPMNSOpt<T>(*this);
  }

  //---------------------------------------------------------------------------
  template <typename T>
  T _OscCalcPMNSOpt<T>::P(int flavBefore, int flavAfter, double E)
  {
    // Normal usage of a calculator in a fit is to configure one set of
    // oscillation parameters and then calculate probabilities for all flavours
    // and for a range of energies. The underlying PMNSOpt object has to redo
    // its expensive diagonalization every time the energy changes, but not
    // when simply switching to a different flavour. Unfortunately callers in
    // CAFAna wind up making the calls in the wrong order to take advantage of
    // this. So, we maintain a map of calculators for different energies (also
    // for neutrinos and antineutrinos) and select the right one if available,
    // on the hope that it's already configured to do what we're asked, since
    // it might already have done a different flavour. This works because
    // callers often pass the same set of energies (bin centers) over and over
    // again.

    // If the caches get too large clear them. The caller is probably giving
    // the precise energy of each event?
    for(int i = 0; i < 2; ++i){
      if(fPMNSOpt[i].size() > 10000){
        fPMNSOpt[i].clear();
      }
    }

    const int anti = (flavBefore > 0) ? +1 : -1;
    assert(flavAfter/anti > 0);

    Val_t& calc = fPMNSOpt[(1+anti)/2][E];
    if(!calc.pmns) calc.pmns.reset(new _PMNSOpt<T>);

    int i = -1, j = -1;
    if(abs(flavBefore) == 12) i = 0;
    if(abs(flavBefore) == 14) i = 1;
    if(abs(flavBefore) == 16) i = 2;
    if(abs(flavAfter) == 12) j = 0;
    if(abs(flavAfter) == 14) j = 1;
    if(abs(flavAfter) == 16) j = 2;
    assert(i >= 0 && j >= 0);

    const bool dirty = (fMixIdx > calc.mixIdx ||
                        fDmIdx > calc.dmIdx ||
                        fLRIdx > calc.lrIdx);

    // If the calculator has out of date parameters update them
    if(fMixIdx > calc.mixIdx){
      calc.pmns->SetMix(this->fTh12, this->fTh23, this->fTh13, this->fdCP);
      calc.mixIdx = fMixIdx;
    }
    if(fDmIdx > calc.dmIdx){
      calc.pmns->SetDeltaMsqrs(this->fDmsq21, this->fDmsq32);
      calc.dmIdx = fDmIdx;
    }

    if(dirty){
      // Cache results for all nine flavour combinations
      for(int ii = 0; ii < 3; ++ii){
        calc.pmns->ResetToFlavour(ii);
        const double Ne = this->fRho * constants::kZPerA;
        calc.pmns->PropMatter(this->fL, E, Ne, anti);
        for(int jj = 0; jj < 3; ++jj){
          calc.P[ii][jj] = calc.pmns->P(jj);
        }
      }

      calc.lrIdx = fLRIdx;
    }

    return calc.P[i][j];
  }
}

//---------------------------------------------------------------------------
template class osc::_OscCalcPMNSOpt<double>;

#ifdef OSCLIB_STAN
template class osc::_OscCalcPMNSOpt<stan::math::var>;
#endif
