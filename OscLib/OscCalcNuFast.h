
#pragma once

#include <string>

#include "OscLib/IOscCalc.h"

#include "TMD5.h"

class TMD5;

namespace osc {
  
  /// This class is an OscLib adapter for the NuFast algorithm (arXiv:2405.02400v1).
  template<typename T>
  class _OscCalcNuFast : public _IOscCalcAdjustable<T> {
  public:
    _OscCalcNuFast(void);
    virtual ~_OscCalcNuFast(void);
    
    virtual _OscCalcNuFast<T>* Copy(void) const override;
    
    // Define setter functions for oscillation parameters.
    virtual void SetL     (double L       ) override { this->fDirty = true; this->fL = L; }
    virtual void SetRho   (double rho     ) override { this->fDirty = true; this->fRho = rho; }
    virtual void SetDmsq21(const T& dmsq21) override { this->fDirty = true; this->fDmsq21 = dmsq21; }
    virtual void SetDmsq32(const T& dmsq32) override { this->fDirty = true; this->fDmsq32 = dmsq32; }
    virtual void SetTh12  (const T& th12  ) override;
    virtual void SetTh13  (const T& th13  ) override;
    virtual void SetTh23  (const T& th23  ) override;
    virtual void SetdCP   (const T& dCP   ) override;
    /// Electron density in matter used by NuFast (default = 0.5).
    virtual void SetYe    (double Ye      ) { this->fDirty = true; this->fYe = Ye; }
    /// Number of Newton steps to improve eigen{value}{vector} convergence (default = 0 seems fine, use 1 or 2 if you need super accuracy).
    virtual void SetNNewton(double nNewton) { this->fDirty = true; this->fNNewton = nNewton; }
    
    // We don't need new getter functions except for Ye and number of Newton corrections (specific to NuFast).
    virtual double GetYe(void) const { return this->fYe; }
    virtual double GetNNewton(void) const { return this->fNNewton; }
    
    /// Implement the NuFast algorithm.
    virtual T P(int flavBefore, int flavAfter, double E) override;
    
    /// Specialized printer function.
    virtual void Print(const std::string& prefix = "") const override;
    
    /// Implementation of params hash.
    virtual TMD5* GetParamsHash(void) const override {
      return _IOscCalcAdjustable<T>::GetParamsHashDefault("NuFast");
    }
  protected:
    /// Copy constructor.
    _OscCalcNuFast(const _OscCalcNuFast<T>& other);
    /// Probability computation helper function.
    inline void RecomputeEnergyIndependentTerms(void);
    inline void RecomputeProbabilityMatrix(double E, bool isAnti);
    
    /// Electron density in matter. NuFast recommends 0.5, and this choice matches other calcs.
    double fYe;
    /// Number of Newton steps to use to improve eigen{value}{vector} slns. Increase me for more accuracy and worse runtime.
    double fNNewton;
    /// Current probabilities. Ex fCurrProbs[1][0] is P(mu->e). Updated only when parameters are reset.
    // TODO Am I also updated when a different energy is used?
    mutable T fCurrProbs[3][3];
    /// Store some precomputations that don't depend on energy.
    mutable struct {
      T s12sq, s13sq, s23sq; // Sine squared of mixing angles. Immediately updated on reset of relevant angle.
      T sind, cosd; // Sine and cosine of delta CP. Immediately updated on reset of dCP.
      T Jrr; // Used for Jarlskog invariant. Depends on all mixing angles.
      T Jmatter; // Jarlskog invariant. Depends on all mixing angles and delta.
      // What follows are some energy-independent terms of the PMNS matrix.
    } fPre;
    /// Do we need to recompute the nontrivial energy-independent terms?
    mutable bool fDirty;
    /// What was the energy that we just used to calculate the probability matrix?
    mutable double fCurrEnergy;
    
    // Constants for unit conversions, using Fermi constant for matter effect.
    double const eVsqkm_to_GeV_over4 = 1e-9 / 1.97327e-7 * 1e3 / 4;
    double const YerhoE2a = 1.52e-4;
  };
  
  typedef _OscCalcNuFast<double> OscCalcNuFast;
}

