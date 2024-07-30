
#pragma once

#include <string>

#include "OscLib/IOscCalc.h"

#include "TMD5.h"

class TMD5;

namespace osc {
  
  /// This class is an OscLib adapter for the NuFast algorithm (arXiv:2405.02400v1).
  template<typename T>
  class _OscCalcNuFastDumb : public _IOscCalcAdjustable<T> {
  public:
    _OscCalcNuFastDumb(void);
    virtual ~_OscCalcNuFastDumb(void);
    
    virtual _OscCalcNuFastDumb<T>* Copy(void) const override;
    
    // Define setter functions for oscillation parameters.
    virtual void SetL     (double L       ) override { this->fL = L; this->fDirty = true; }
    virtual void SetRho   (double rho     ) override { this->fRho = rho; this->fDirty = true; }
    virtual void SetDmsq21(const T& dmsq21) override { this->fDmsq21 = dmsq21; this->fDirty = true; }
    virtual void SetDmsq32(const T& dmsq32) override { this->fDmsq32 = dmsq32; this->fDirty = true; }
    virtual void SetTh12  (const T& th12  ) override { this->fTh12 = th12; this->fDirty = true; }
    virtual void SetTh13  (const T& th13  ) override { this->fTh13 = th13; this->fDirty = true; }
    virtual void SetTh23  (const T& th23  ) override { this->fTh23 = th23; this->fDirty = true; }
    virtual void SetdCP   (const T& dCP   ) override { this->fdCP = dCP; this->fDirty = true; }
    /// Electron density in matter used by NuFast (default = 0.5).
    virtual void SetYe    (double Ye      ) { this->fYe = Ye; this->fDirty = true; }
    /// Number of Newton steps to improve eigen{value}{vector} convergence (default = 0 seems fine, use 1 or 2 if you need super accuracy).
    virtual void SetNNewton(double nNewton) { this->fNNewton = nNewton; this->fDirty = true; }
    
    // We don't need new getter functions except for Ye and number of Newton corrections (specific to NuFast).
    virtual double GetYe(void) const { return this->fYe; }
    virtual double GetNNewton(void) const { return this->fNNewton; }
    
    /// Implement the NuFast algorithm.
    virtual T P(int flavBefore, int flavAfter, double E) override;
    
    /// Specialized printer function.
    virtual void Print(const std::string& prefix = "") const override;
    
    /// Implementation of params hash.
    virtual TMD5* GetParamsHash(void) const override {
      return _IOscCalcAdjustable<T>::GetParamsHashDefault("NuFastDumb");
    }
  protected:
    /// Copy constructor.
    _OscCalcNuFastDumb(const _OscCalcNuFastDumb<T>& other);
    /// Probability computation helper function.
    inline void RecomputeProbabilityMatrix(double E, bool isAnti);
    
    /// Electron density in matter. NuFast recommends 0.5, and this choice matches other calcs.
    double fYe;
    /// Number of Newton steps to use to improve eigen{value}{vector} slns. Increase me for more accuracy and worse runtime.
    double fNNewton;
    /// Current probabilities. Ex fCurrProbs[1][0] is P(mu->e). Updated only when parameters are reset.
    // TODO Am I also updated when a different energy is used?
    mutable T fCurrProbs[3][3];
    /// Has a parameter been (re)set since we last computed probabilites?
    mutable bool fDirty;
    /// What was the energy that we just used to calculate the probability matrix?
    mutable double fCurrEnergy;
    
    // Constants for unit conversions, using Fermi constant for matter effect.
    double const eVsqkm_to_GeV_over4 = 1e-9 / 1.97327e-7 * 1e3 / 4;
    double const YerhoE2a = 1.52e-4;
  };
  
  typedef _OscCalcNuFastDumb<double> OscCalcNuFastDumb;
}

