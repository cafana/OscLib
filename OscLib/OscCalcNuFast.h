
#ifndef OSCLIBCALCNUFAST_H
#define OSCLIBCALCNUFAST_H

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
  public:
    _OscCalcNuFast(void);
    virtual ~_OscCalcNuFast(void);
    
    virtual _OscCalcNuFast<T>* Copy(void) const override;
    
    // Define setter functions for oscillation parameters.
    virtual void SetL     (double L       ) override { this->fIsDirty = true; this->fL = L; }
    virtual void SetRho   (double rho     ) override { this->fIsDirty = true; this->fRho = rho; }
    virtual void SetDmsq21(const T& dmsq21) override { this->fIsDirty = true; this->fDmsq21 = dmsq21; }
    virtual void SetDmsq32(const T& dmsq32) override { this->fIsDirty = true; this->fDmsq32 = dmsq32; }
    virtual void SetTh12  (const T& th12  ) override;
    virtual void SetTh13  (const T& th13  ) override;
    virtual void SetTh23  (const T& th23  ) override;
    virtual void SetdCP   (const T& dCP   ) override;
    /// Electron density in matter used by NuFast (default = 0.5).
    virtual void SetYe    (double Ye      ) { this->fIsDirty = true; this->fYe = Ye; }
    /// Number of Newton steps to improve eigen{value}{vector} convergence (default = 0 seems fine, use 1 or 2 if you need super accuracy).
    virtual void SetNNewton(double nNewton) { this->fIsDirty = true; this->fNNewton = nNewton; }
    virtual void SetAggressive(bool isAggressive) { this->fIsAggressive = isAggressive; }
    
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
    /// Probability computation helper function.
    void inline __attribute__((always_inline)) RecomputeEnergyIndependentTerms(void);
    /// from and to aren't necessary for computing the matrix but are used for the return (so the cache doesn't need to be searched for the result).
    template<class VT, class KVT> inline __attribute__((always_inline)) VT RecomputeProbabilityMatrix(int from, int to, KVT E);
    // Master probability function called with some combination of primitives, stan::math::vars, and Eigen objects.
    template<class VT, class KVT> inline __attribute__((always_inline)) VT _P(int from, int to, KVT E);
    
    /// Electron density in matter. NuFast recommends 0.5, and this choice matches other calcs.
    double fYe;
    /// Number of Newton steps to use to improve eigen{value}{vector} slns. Increase me for more accuracy and worse runtime.
    double fNNewton;
    /// Should we aggressively precompute terms not dependent on energy?
    bool fIsAggressive;
    /// Store some precomputations that don't depend on energy.
    struct {
      T Dmsq31; // Squared mass splitting between m3 and m1.
      T Dmsqee;
      T s12, s13, s23; // Sines of mixing angles.
      T s12sq, s13sq, s23sq; // Sine squared of mixing angles. Immediately updated on reset of relevant angle.
      T c12sq, c13sq, c23sq; // Cosine squared of mixing angles. Immediately updated on reset of relevant angle.
      T sind, cosd; // Sine and cosine of delta CP. Immediately updated on reset of dCP.
      T c13sqxs12sq, c13sqxs23sq, c12sqxc23sq; // This is cos(th13)^2 * sin(th12)^2, etc.
		  T s13sqxs12sqxs23sq;
      T Jrr; // Used for Jarlskog invariant.
      T Jmatter_first; // Jarlskog times masses before an energy-dependent term is used.
      T Um2sq_first;
      T See, Tee;
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

