#ifndef OSC_OSCCALCULATORPMNS_NSI_H
#define OSC_OSCCALCULATORPMNS_NSI_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// \file OscCalcPMNS_NSI.h                                              //
//                                                                      //
// Adapt the PMNS_NSI calculator to standard interface                  //
// <c.backhouse@ucl.ac.uk>						//
// Modifications by MAAceroO <marioacero@mail.uniatlantico.edu.co>      //
//////////////////////////////////////////////////////////////////////////

#include "OscLib/IOscCalc.h"
#include "OscLib/PMNS_NSI.h"

namespace osc
{
  /// \brief Optimized version of \ref OscCalcPMNS
  ///
  /// Adapt the \ref PMNS_NSI calculator to standard interface
  template <typename T>
  class _OscCalcPMNS_NSI: public _IOscCalcAdjustable<T>
  {
  public:
    using _IOscCalc<T>::P;
    _OscCalcPMNS_NSI();
    virtual ~_OscCalcPMNS_NSI();

    _IOscCalcAdjustable<T>* Copy() const override;

    virtual TMD5* GetParamsHash() const override;

    virtual T P(int flavBefore, int flavAfter, double E) override;

    // Standard oscillation parameters
    virtual void SetL     (double   L     ) override {fPropDirty = true; this->fL      = L;}
    virtual void SetRho   (double   rho   ) override {fPropDirty = true; this->fRho    = rho;}
    virtual void SetDmsq21(const T& dmsq21) override {fDmDirty   = true; this->fDmsq21 = dmsq21;}
    virtual void SetDmsq32(const T& dmsq32) override {fDmDirty   = true; this->fDmsq32 = dmsq32;}
    virtual void SetTh12  (const T& th12  ) override {fMixDirty  = true; this->fTh12   = th12;}
    virtual void SetTh13  (const T& th13  ) override {fMixDirty  = true; this->fTh13   = th13;}
    virtual void SetTh23  (const T& th23  ) override {fMixDirty  = true; this->fTh23   = th23;}
    virtual void SetdCP   (const T& dCP   ) override {fMixDirty  = true; this->fdCP    = dCP;}

    // Non-Standard Interactions parameters (three real -diagonal- and thee complex -off-diagonal: norm and phase-)
    virtual void SetEps_ee      (const T& eps_ee     ){fEpsDirty = true; fEps_ee      = eps_ee;}
    virtual void SetEps_emu     (const T& eps_emu    ){fEpsDirty = true; fEps_emu     = eps_emu;}
    virtual void SetEps_etau    (const T& eps_etau   ){fEpsDirty = true; fEps_etau    = eps_etau;}
    virtual void SetEps_mumu    (const T& eps_mumu   ){fEpsDirty = true; fEps_mumu    = eps_mumu;}
    virtual void SetEps_mutau   (const T& eps_mutau  ){fEpsDirty = true; fEps_mutau   = eps_mutau;}
    virtual void SetEps_tautau  (const T& eps_tautau ){fEpsDirty = true; fEps_tautau  = eps_tautau;}
    virtual void SetDelta_emu   (const T& Delta_emu  ){fEpsDirty = true; fDelta_emu   = Delta_emu;}
    virtual void SetDelta_etau  (const T& Delta_etau ){fEpsDirty = true; fDelta_etau  = Delta_etau;}
    virtual void SetDelta_mutau (const T& Delta_mutau){fEpsDirty = true; fDelta_mutau = Delta_mutau;}

     // Setting the state (oscillation parameters)
    void SetState(std::vector<double> state);

    // Getters
    T GetEps_ee()      const { return fEps_ee; }
    T GetEps_emu()     const { return fEps_emu; }
    T GetEps_etau()    const { return fEps_etau; }
    T GetEps_mumu()    const { return fEps_mumu; }
    T GetEps_mutau()   const { return fEps_mutau; }
    T GetEps_tautau()  const { return fEps_tautau; }
    T GetDelta_emu()   const { return fDelta_emu; }
    T GetDelta_etau()  const { return fDelta_etau; }
    T GetDelta_mutau() const { return fDelta_mutau; }

     // Getting the state (oscillation parameters)
    std::vector<double> GetState() const;

  protected:
    _PMNS_NSI<T> fPMNS_NSI;

    bool fMixDirty;
    bool fDmDirty;
    bool fPropDirty;
    bool fEpsDirty;
    double fPrevE;
    int fPrevAnti;

    T fEps_ee;
    T fEps_mumu;
    T fEps_tautau;
    T fEps_emu;
    T fEps_etau;
    T fEps_mutau;
    T fDelta_emu;
    T fDelta_etau;
    T fDelta_mutau;

  };

  typedef _OscCalcPMNS_NSI<double> OscCalcPMNS_NSI;

  const OscCalcPMNS_NSI* DowncastToNSI(const IOscCalc* calc);
  OscCalcPMNS_NSI* DowncastToNSI(IOscCalc* calc);

  template <typename T>
  const _OscCalcPMNS_NSI<T>* DowncastToNSI(const _IOscCalc<T>* calc)
  {
    return dynamic_cast<const _OscCalcPMNS_NSI<T>*>(calc);
  }

  template <typename T>
  _OscCalcPMNS_NSI<T>* DowncastToNSI(_IOscCalc<T>* calc)
  {
    return dynamic_cast<_OscCalcPMNS_NSI<T>*>(calc);
  }

} // namespace

#endif
