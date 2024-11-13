
#ifndef OSCLIBCALCNUFAST_H
#define OSCLIBCALCNUFAST_H

// \file OscCalcNuFast.h
// \brief Implementation of the NuFast algorithm by P. Denton (BNL) and S. Parke (FNAL).
//        The algorithm makes a few approximations specific to LBL matter effects to make probability computations faster.
//        See arXiv:2405.02400.
// \author <cullenms@fnal.gov>

#include <string>

#include "OscLib/Cache.h"
#include "OscLib/IOscCalc.h"

#include "TMD5.h"

namespace Eigen {
	template<class T> using ArrayX = Eigen::Array<T, Eigen::Dynamic, 1>;
}

namespace osc {
  /// This class is an OscLib adapter for the NuFast algorithm (arXiv:2405.02400v1).
  template<typename T>
  class _OscCalcNuFast : public _IOscCalcAdjustable<T>,
                         protected analytic::ProbCache<double, T>,
                         protected analytic::ProbCache<Eigen::ArrayXd, Eigen::ArrayX<T>> {
  // Warning: The cache for Stan-templated calculators is dangerous to use because Stan uses its own cache internally,
  // and the Stan cache can be cleared without warning. To avoid the issue, CAFAna's StanFitter forcibly clears
  // this cache for calculatos between proposed steps. We leave it to the user to use the cache for Stan types safely.
  
  public:
    _OscCalcNuFast(void);
    virtual ~_OscCalcNuFast(void) = default;
    
    virtual _OscCalcNuFast<T>* Copy(void) const override;
    
    // Define setter functions for oscillation parameters.
    virtual void SetL     (double L       ) override { this->fIsDirty = true; this->fL = L; }
    virtual void SetRho   (double rho     ) override { this->fIsDirty = true; this->fRho = rho; }
    virtual void SetDmsq21(const T& dmsq21) override { this->fIsDirty = true; this->fDmsq21 = dmsq21; }
    // Like other calculators, indicate normal [inverted] mass ordering with positive [negative] dmsq32.
    virtual void SetDmsq32(const T& dmsq32) override { this->fIsDirty = true; this->fDmsq32 = dmsq32; }
    virtual void SetTh12  (const T& th12  ) override;
    virtual void SetTh13  (const T& th13  ) override;
    virtual void SetTh23  (const T& th23  ) override;
    virtual void SetdCP   (const T& dCP   ) override;
    /// Electron density in matter used by NuFast (default = 0.5).
    virtual void SetYe    (double Ye      ) { this->fIsDirty = true; this->fYe = Ye; }
    /// Number of Newton steps to improve eigen{value}{vector} convergence (default = 0 seems fine, use 1 or 2 if you need super accuracy).
    virtual void SetNNewton(double nNewton) { this->fIsDirty = true; this->fNNewton = nNewton; }
    
    // We don't need new getter functions except for Ye and number of Newton corrections (specific to NuFast).
    virtual double GetYe(void) const { return this->fYe; }
    virtual double GetNNewton(void) const { return this->fNNewton; }
    
    /// Implement the NuFast algorithm.
    virtual T P(int from, int to, double E) override;
    virtual Eigen::ArrayX<T> P(int from, int to, const std::vector<double>& E) override;
    virtual Eigen::ArrayX<T> P(int from, int to, const Eigen::ArrayXd& E) override;
    
    /// Specialized printer function.
    virtual void Print(const std::string& prefix = "") const override;
    
    /// Implementation of params hash.
    virtual TMD5* GetParamsHash(void) const override {
      return _IOscCalcAdjustable<T>::GetParamsHashDefault("NuFast");
    }
  protected:
    void ClearProbCaches(void) {
      analytic::ProbCache<double, T>::clear();
      analytic::ProbCache<Eigen::ArrayXd, Eigen::ArrayX<T>>::clear();
    }
    /// Copy constructor.
    _OscCalcNuFast(const _OscCalcNuFast<T>& other);
    /// from and to aren't necessary for computing the matrix but are used for the return (so the cache doesn't need to be searched for the result).
    template<class VT, class KVT> VT RecomputeProbabilityMatrix(int from, int to, KVT E);
    // Master probability function called with some combination of primitives, stan::math::vars, and Eigen objects.
    template<class VT, class KVT> VT _P(int from, int to, KVT E);
    
    /// Electron density in matter. NuFast recommends 0.5, and this choice matches other calcs.
    double fYe;
    /// Number of Newton steps to use to improve eigen{value}{vector} slns. Increase me for more accuracy and worse runtime.
    unsigned int fNNewton;
    /// Store some precomputations that don't depend on energy.
    struct {
      T s12, s13, s23; // Sines of mixing angles.
      T s12sq, s13sq, s23sq; // Sine squared of mixing angles. Immediately updated on reset of relevant angle.
      T c12sq, c13sq, c23sq; // Cosine squared of mixing angles. Immediately updated on reset of relevant angle.
      T sind, cosd; // Sine and cosine of delta CP. Immediately updated on reset of dCP.
    } fPre;
    
    /// Do we need to recompute the nontrivial energy-independent terms?
    mutable bool fIsDirty;
    
    // Constants for unit conversions, using Fermi constant for matter effect.
    double const eVsqkm_to_GeV_over4 = 1e-9 / 1.97327e-7 * 1e3 / 4;
    double const YerhoE2a = 1.52e-4;
  };
  
  typedef _OscCalcNuFast<double> OscCalcNuFast;
}

#endif

