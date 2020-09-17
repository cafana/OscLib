#ifndef OSC_OSCCALCULATORPMNSOPT_H
#define OSC_OSCCALCULATORPMNSOPT_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// \file OscCalcPMNSOpt.h                                               //
//                                                                      //
// Adapt the PMNSOpt calculator to standard interface                   //
// <c.backhouse@ucl.ac.uk>						//
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "OscLib/IOscCalc.h"
#include "OscLib/PMNSOpt.h"

#include <cassert>
#include <unordered_map>

namespace osc
{
  /// \brief Optimized version of \ref OscCalcPMNS
  ///
  /// Adapt the \ref PMNSOpt calculator to standard interface
  template <typename T>
  class _OscCalcPMNSOpt: public _IOscCalcAdjustable<T>
  {
    public:
    using _IOscCalc<T>::P;
      _OscCalcPMNSOpt();
      virtual ~_OscCalcPMNSOpt();

      _IOscCalcAdjustable<T>* Copy() const override;

      T P(int flavBefore, int flavAfter, double E) override;

      void SetL     (double   L     ) override {++fLRIdx;  this->fL      = L;}
      void SetRho   (double   rho   ) override {++fLRIdx;  this->fRho    = rho;}
      void SetDmsq21(const T& dmsq21) override {++fDmIdx;  this->fDmsq21 = dmsq21;}
      void SetDmsq32(const T& dmsq32) override {++fDmIdx;  this->fDmsq32 = dmsq32;}
      void SetTh13  (const T& th13  ) override {++fMixIdx; this->fTh13   = th13;}
      void SetTh12  (const T& th12  ) override {++fMixIdx; this->fTh12   = th12;}
      void SetTh23  (const T& th23  ) override {++fMixIdx; this->fTh23   = th23;}
      void SetdCP   (const T& dCP   ) override {++fMixIdx; this->fdCP    = dCP;}

      TMD5* GetParamsHash() const override
      {
        return _IOscCalcAdjustable<T>::GetParamsHashDefault("PMNSOpt");
      }

    protected:
      // How many times the mixing parameters and splittings have been set
      long fMixIdx;
      long fDmIdx;
      long fLRIdx;

      struct Val_t
      {
        Val_t() : mixIdx(-1), dmIdx(-1), lrIdx(-1), pmns(0) {}

        // How many times the mixing parameters and splittings had been set when
        // 'pmns' was last updated. If too small then 'pmns' must be updated
        // before use.
        long mixIdx;
        long dmIdx;
        long lrIdx;
        T    P[3][3]; ///< Cache of oscillation probabilities
        _PMNSOpt<T>* pmns;  ///< The calculator itself
      };

    std::unordered_map<double, Val_t> fPMNSOpt[2]; // [anti][E]
  };
  typedef _OscCalcPMNSOpt<double> OscCalcPMNSOpt;

} // namespace

#endif
