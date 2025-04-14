
#ifdef OSCLIB_STAN
#include "OscLib/Stan.h"
#endif

#include "OscLib/OscCalcNuFast.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace osc {
  
  template<typename T>
  _OscCalcNuFast<T>::_OscCalcNuFast(void) :
    fYe(0.5), fNNewton(0), fIsDirty(true)
    {}
  
  template<typename T>
  _OscCalcNuFast<T>::_OscCalcNuFast(const _OscCalcNuFast<T>& other) :
    _IOscCalcAdjustable<T>(other), fYe(other.GetYe()), fNNewton(other.GetNNewton()), fIsDirty(true)
    {}
  
  template<typename T>
  _OscCalcNuFast<T>* _OscCalcNuFast<T>::Copy(void) const { return new _OscCalcNuFast<T>(*this); }
  
  template<typename T>
  void _OscCalcNuFast<T>::SetTh12(const T& th12) {
    this->fTh12 = th12;
    // Precomute sine squared and cosine squared of th12.
    this->fPre.s12 = sin(th12);
    this->fPre.s12sq = this->fPre.s12*this->fPre.s12;
    this->fPre.c12sq = 1-this->fPre.s12sq;
    this->fIsDirty = true;
  }
  
  template<typename T>
  void _OscCalcNuFast<T>::SetTh13(const T& th13) {
    this->fTh13 = th13;
    this->fPre.s13 = sin(th13);
    this->fPre.s13sq = this->fPre.s13*this->fPre.s13;
    this->fPre.c13sq = 1-this->fPre.s13sq;
    this->fIsDirty = true;
  }
  
  template<typename T>
  void _OscCalcNuFast<T>::SetTh23(const T& th23) {
    this->fTh23 = th23;
    this->fPre.s23 = sin(th23);
    this->fPre.s23sq = this->fPre.s23*this->fPre.s23;
    this->fPre.c23sq = 1-this->fPre.s23sq;
    this->fIsDirty = true;
  }
  
  template<typename T>
  void _OscCalcNuFast<T>::SetdCP(const T& dCP) {
    this->fdCP = dCP;
    // Precompute the sine and cosine of delta CP.
    this->fPre.sind = sin(dCP);
    this->fPre.cosd = cos(dCP);
    this->fIsDirty = true;
  }
  
  // The following code is adapted from the NuFast oscillation calculator (Probability_Matter_LBL).
  template<typename T> template<class VT, class KVT>
  VT _OscCalcNuFast<T>::_P(int from, int to, KVT E) {
    // First, make sure the PDG codes have the same sign.
    assert((from > 0 && to > 0) || (from < 0 && to < 0) );
    // Nice trick from analytic: If we need antineutrinos, just flip flavor signs and send E -> -E. The NuFast algorithm already uses -E means antineutrinos.
    if ( from < 0 ) { return P(-from,-to,-E); }
    
    // The probability cache was invalidated if the user updated any parameters.
    if ( this->fIsDirty ) {
      ClearProbCaches();
      this->fIsDirty = false;
    }
    
    // Try to find the requested energy in the probabilities cache. E is possibly negative (for antineutrinos).
    typename analytic::ProbCache<KVT,VT>::iterator it = analytic::ProbCache<KVT,VT>::find(E);
    if ( it != analytic::ProbCache<KVT,VT>::end() ) {
      return it->second.P(from,to);
    }
    
    // Need to do energy-dependent computation.
    return this->RecomputeProbabilityMatrix<VT,KVT>(from,to,E);
  }
  
  template<typename T> T _OscCalcNuFast<T>::P(int from, int to, double E) {
    //std::cout << "NuFast: Calling double overload." << std::endl;
    return _P<T>(from,to,E);
  }
  
  template<typename T> Eigen::ArrayX<T> _OscCalcNuFast<T>::P(int from, int to, const Eigen::ArrayXd& E)
  {
    //std::cout << "NuFast: Calling Eigen::ArrayXd overload." << std::endl;
    return _P<Eigen::ArrayX<T>>(from,to,E);
  }
  
  template<typename T> Eigen::ArrayX<T> _OscCalcNuFast<T>::P(int from, int to, const std::vector<double>& E)
  {
    //std::cout << "NuFast: Calling std::vector<double> overload (size = " << E.size() << ")." << std::endl;
    return P(from,to,Eigen::Map<const Eigen::ArrayXd>(E.data(),E.size()));
  }
  
  template<typename T> template<class VT, class KVT>
  VT _OscCalcNuFast<T>::RecomputeProbabilityMatrix(int from, int to, KVT E) {
    // Now we need to interface internally used NuFast variables with the OscCalc members.
    #define s12sq this->fPre.s12sq
    #define s13sq this->fPre.s13sq
    #define s23sq this->fPre.s23sq
    #define sind this->fPre.sind
    #define cosd this->fPre.cosd
    #define c12sq this->fPre.c12sq
    #define c13sq this->fPre.c13sq
    #define c23sq this->fPre.c23sq
    
    #define delta this->fdCP
    #define Dmsq21 this->fDmsq21
    #define Dmsq32 this->fDmsq32
    #define L this->fL
    #define rho this->fRho
    #define Ye this->fYe
    #define N_Newton this->fNNewton
    #define probs_returned this->fCurrProbs
    
    // --------------------------------------------------------------------- //
    // First calculate useful simple functions of the oscillation parameters //
    // --------------------------------------------------------------------- //
    
    // NOvA doesn't use dmsq31 by convention, but we can compute it if we know dmsq32 and dmsq21.
    // This works for both mass orderings because we indicate NO [IO] with positive [negative] dmsq32.
    // Therefore dmsq32 and dmsq21 will have the same [opposite] sign for NO [IO].
    const T
    Dmsq31 = Dmsq32 + Dmsq21, // dmsq31 = msq3 - msq1 = (msq3-msq2) + (msq2-msq1) = dmsq32 + dmsq21.
    Dmsqee = Dmsq31 - s12sq * Dmsq21,
    
    c13sqxs12sq = c13sq * s12sq,
    c13sqxs23sq = c13sq * s23sq,
    c12sqxc23sq = c12sq * c23sq,
    s13sqxs12sqxs23sq = s13sq * s12sq * s23sq,
    
    Jrr = sqrt(c12sqxc23sq * s13sqxs12sqxs23sq),
    Jmatter_first = 8 * Jrr * c13sq * sind 
                      * Dmsq21 * Dmsq31 * (Dmsq31 - Dmsq21),
    
    Um2sq_first = c12sqxc23sq + s13sqxs12sqxs23sq - 2 * Jrr * cosd,
    See = Dmsq21+Dmsq31 - Dmsq21 * c13sqxs12sq - Dmsq31 * s13sq,
    Tee = Dmsq21 * Dmsq31 * (c13sq - c13sqxs12sq);
    
    const VT
    Amatter = (Ye*rho*YerhoE2a) * E,
    
    // calculate A, B, C, See, Tee, and part of Tmm
    C = Amatter * Tee,
    A = Dmsq21+Dmsq31 + Amatter,
    
    // ---------------------------------- //
    // Get lambda3 from lambda+ of MP/DMP //
    // ---------------------------------- //
    xmat = Amatter / Dmsqee;
    VT lambda3 = Dmsq31 + (0.5 * Dmsqee) * (xmat - 1 + sqrt((1-xmat) * (1-xmat) + (4*s13sq) * xmat));
    
    // ---------------------------------------------------------------------------- //
    // Newton iterations to improve lambda3 arbitrarily, if needed, (B needed here) //
    // ---------------------------------------------------------------------------- //
    if ( N_Newton > 0 ) {
      const VT B = Dmsq21*Dmsq31 + Amatter * See; // B is only needed for N_Newton >= 1
      for ( unsigned int i = 0 ; i < N_Newton ; i++ ) {
        lambda3 = (lambda3 * lambda3 * (lambda3 + lambda3 - A) + C) / (lambda3 * (2 * (lambda3 - A) + lambda3) + B); // this strange form prefers additions to multiplications
      }
    }
    
    // ------------------- //
    // Get  Delta lambda's //
    // ------------------- //
    const VT
    Dlambda21 = sqrt((A - lambda3) * (A - lambda3) - 4 * C / lambda3),
    lambda2 = 0.5 * (A - lambda3 + Dlambda21),
    Dlambda32 = lambda3 - lambda2,
    Dlambda31 = Dlambda32 + Dlambda21,
    
    // ----------------------- //
    // Use Rosetta for Veisq's //
    // ----------------------- //
    // denominators	  
    PiDlambdaInv = 1 / (Dlambda31 * Dlambda32 * Dlambda21),
    Xp3 = PiDlambdaInv * Dlambda21,
    Xp2 = -PiDlambdaInv * Dlambda31,
    
    // numerators
    Ue3sq = (lambda3 * (lambda3 - See) + Tee) * Xp3,
    Ue2sq = (lambda2 * (lambda2 - See) + Tee) * Xp2,
    
    Smm = A - Dmsq21 * Um2sq_first - Dmsq31 * c13sqxs23sq,
    Tmm = Dmsq21*Dmsq31 * (1 - c13sqxs23sq - Um2sq_first) + Amatter * (See + Smm - A),
    
    Um3sq = (lambda3 * (lambda3 - Smm) + Tmm) * Xp3,
    Um2sq = (lambda2 * (lambda2 - Smm) + Tmm) * Xp2,
    
    // ------------- //
    // Use NHS for J //
    // ------------- //
    Jmatter = Jmatter_first * PiDlambdaInv,
    
    // ----------------------- //
    // Get all elements of Usq //
    // ----------------------- //
    Ue1sq = 1 - Ue3sq - Ue2sq,
    Um1sq = 1 - Um3sq - Um2sq,
    
    Ut3sq = 1 - Um3sq - Ue3sq,
    Ut2sq = 1 - Um2sq - Ue2sq,
    Ut1sq = 1 - Um1sq - Ue1sq,
    
    // ----------------------- //
    // Get the kinematic terms //
    // ----------------------- //
    Lover4E = (eVsqkm_to_GeV_over4 * L) / E,
    
    D21 = Dlambda21 * Lover4E,
    D32 = Dlambda32 * Lover4E,
	  
    sinD21 = sin(D21),
    sinD31 = sin(D32 + D21),
    sinD32 = sin(D32),
    
    triple_sin = sinD21 * sinD31 * sinD32,
    
    sinsqD21_2 = 2 * sinD21 * sinD21,
    sinsqD31_2 = 2 * sinD31 * sinD31,
    sinsqD32_2 = 2 * sinD32 * sinD32,
    
    // ------------------------------------------------------------------- //
    // Calculate the three necessary probabilities, separating CPC and CPV //
    // ------------------------------------------------------------------- //
    Pme_CPC = (Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2
            + (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2
	          + (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2,
    Pme_CPV = -Jmatter * triple_sin,
    
    Pmm = 1 - 2 * (Um2sq * Um1sq * sinsqD21_2
        + Um3sq * Um1sq * sinsqD31_2
        + Um3sq * Um2sq * sinsqD32_2),
    
    Pee = 1 - 2 * (Ue2sq * Ue1sq * sinsqD21_2
        + Ue3sq * Ue1sq * sinsqD31_2
        + Ue3sq * Ue2sq * sinsqD32_2);
    
    // Store the probabilities in the cache.
    const analytic::Probs<VT> ps(Pee,Pme_CPC - Pme_CPV,Pme_CPC + Pme_CPV,Pmm);
    analytic::ProbCache<KVT, VT>::emplace(E,ps);
    // Return the probability.
    return ps.P(from,to);
    
    // Undefine all the things.
    #undef s12sq
    #undef s13sq
    #undef s23sq
    #undef sind
    #undef cosd
    #undef c12sq
    #undef c13sq
    #undef c23sq
    
    #undef delta
    #undef Dmsq21
    #undef Dmsq32
    #undef L
    #undef rho
    #undef Ye
    #undef N_Newton
    #undef probs_returned
  }
  
  template<typename T>
  void _OscCalcNuFast<T>::Print(const std::string& prefix) const {
    std::cout << prefix << "dmsq21 = " << this->fDmsq21 << " eV^2\n"
              << prefix << "dmsq32 = " << this->fDmsq32 << " eV^2\n"
              << prefix << "th12 = " << this->fTh12 << "\n"
              << prefix << "th13 = " << this->fTh13 << "\n"
              << prefix << "th23 = " << this->fTh23 << "\n"
              << prefix << "dCP = " << this->fdCP << "\n"
              << prefix << "L = " << this->fL << " km\n"
              << prefix << "rho = " << this->fRho << " g/cm^3\n"
              << prefix << "Ye = " << this->fYe << "\n"
              << prefix << "N Newton = " << this->fNNewton << "\n" << std::endl;
  }
}

template class osc::_OscCalcNuFast<double>;
#ifdef OSCLIB_STAN
  template class osc::_OscCalcNuFast<stan::math::var>;
#endif

