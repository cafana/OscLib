
#ifdef OSCLIB_STAN
#include "OscLib/Stan.h"
#endif

#include "OscLib/OscCalcNuFast.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <tuple>
#include <vector>

namespace osc {
  template<typename T>
  _OscCalcNuFast<T>::_OscCalcNuFast(void) :
    fYe(0.5), fNNewton(1), fIsAggressive(true), fIsDirty(true)
    {}
  
  template<typename T>
  _OscCalcNuFast<T>::_OscCalcNuFast(const _OscCalcNuFast<T>& other) :
    _IOscCalcAdjustable<T>(other), fYe(other.GetYe()), fNNewton(other.GetNNewton()), fIsAggressive(true), fIsDirty(true)
    {}
  
  template<typename T>
  _OscCalcNuFast<T>::~_OscCalcNuFast(void) {}
  
  template<typename T>
  _OscCalcNuFast<T>* _OscCalcNuFast<T>::Copy(void) const { return new _OscCalcNuFast<T>(*this); }
  
  template<typename T>
  void _OscCalcNuFast<T>::SetTh12(const T& th12) {
    this->fIsDirty = true;
    this->fTh12 = th12;
    // Precomute sine squared of th12.
    this->fPre.s12 = sin(th12);
    this->fPre.s12sq = this->fPre.s12*this->fPre.s12;
  }
  
  template<typename T>
  void _OscCalcNuFast<T>::SetTh13(const T& th13) {
    this->fIsDirty = true;
    this->fTh13 = th13;
    this->fPre.s13 = sin(th13);
    this->fPre.s13sq = this->fPre.s13*this->fPre.s13;
  }
  
  template<typename T>
  void _OscCalcNuFast<T>::SetTh23(const T& th23) {
    this->fIsDirty = true;
    this->fTh23 = th23;
    this->fPre.s23 = sin(th23);
    this->fPre.s23sq = this->fPre.s23*this->fPre.s23;
  }
  
  template<typename T>
  void _OscCalcNuFast<T>::SetdCP(const T& dCP) {
    this->fIsDirty = true;
    this->fdCP = dCP;
    // Precompute the sine and cosine of delta CP.
    this->fPre.sind = sin(dCP);
    this->fPre.cosd = cos(dCP);
  }
  
  // The following code is a copy of the NuFast oscillation calculator (Probability_Matter_LBL) except pointing to class members.
  
  // Probability_Matter_LBL calculates all nine oscillation probabilities including
  // the matter effect in an optimized, fast, and efficient way. The precision can
  // be controlled with N_Newton. For many applications N_Newton=0 may be enough,
  // but many years of DUNE or HK-LBL may require N_Newton=1. This code may be
  // suitable for atmospheric neutrinos. The code is standalone.
  //
  // Inputs:
  //   mixing angles (usual parameterization)
  //   phase (usual parameterization) make Dmsq31 positive/negative for the NO/IO
  //   Delta msq's (eV^2)
  //   L (km)
  //   E (GeV) positive for neutrinos, negative for antineutrinos
  //   rho (g/cc)
  //   Ye: electron fraction, typically around 0.5
  //   N_Newton: number of Newton's method iterations to do. should be zero, one, two (or higher)
  // Outputs:
  //   probs_returned is all nine oscillation probabilities: e.g. probs_returned[1][0] is mu->e
  
  // The following is the function signature out-of-the-box from NuFast.
  //void Probability_Matter_LBL(double s12sq, double s13sq, double s23sq,
  //                            double delta, double Dmsq21, double Dmsq31,
  //                            double L, double E, double rho,
  //                            double Ye, int N_Newton, double (*probs_returned)[3][3])
  //{
  template<typename T> template<class VT, class KVT>
  VT _OscCalcNuFast<T>::_P(int from, int to, KVT E) {
    // First, make sure the PDG codes have the same sign.
    assert((from > 0 && to > 0) || (from < 0 && to < 0) );
    // Nice trick from analytic: If we need antineutrinos, just flip flavor signs and send E -> -E. The NuFast algorithm already uses -E means antineutrinos.
    if ( from < 0 ) { return P(-from,-to,-E); }
    
    // Recompute energy independent PMNS terms if we've updated any parameters.
    if ( !this->fIsAggressive || this->fIsDirty ) {
      //std::cout << "NuFast: Recomputing energy independent terms." << std::endl;
      this->RecomputeEnergyIndependentTerms();
    }
    
    // Try to find the requested energy in the probabilities cache. E is possibly negative (for antineutrinos).
    typename analytic::ProbCache<KVT,VT>::iterator it = analytic::ProbCache<KVT,VT>::find(E);
    if ( it != analytic::ProbCache<KVT,VT>::end() ) {
      //std::cout << "NuFast: Just using probability cache." << std::endl;
      return it->second.P(from,to);
    }
    
    // Need to do energy-dependent computation.
    //std::cout << "NuFast: Recomputing probability matrix and returning result." << std::endl;
    return this->RecomputeProbabilityMatrix<VT,KVT>(from,to,E);
    
    // Recompute the probability matrix if necessary (either a parameter has been updated or we used a different energy before).
    //std::cout << "E = " << E << ", currE = " << this->fCurrEnergy << ", diff = " << this->fCurrEnergy - (isAnti ? -E : E) << std::endl;
    /*
    if ( std::abs(this->fCurrEnergy - (isAnti ? -E : E)) > std::numeric_limits<double>::epsilon() ) {
    //  std::cout << "Nufast info: Recomputing probability matrix..." << std::endl;
      this->RecomputeProbabilityMatrix(E,isAnti);
    }
    */
    
    // We have good PDG codes. Now return the corresponding element of the probability matrix.
    /*
    unsigned int indexBefore = -1, indexAfter = -1;
    // I'd use std::unordered_map for pdg->index, but this should be faster.
    if ( flavBefore == +12 || flavBefore == -12 ) { indexBefore = 0; }
    if ( flavBefore == +14 || flavBefore == -14 ) { indexBefore = 1; }
    if ( flavBefore == +16 || flavBefore == -16 ) { indexBefore = 2; }
    if ( flavAfter  == +12 || flavAfter  == -12 ) { indexAfter  = 0; }
    if ( flavAfter  == +14 || flavAfter  == -14 ) { indexAfter  = 1; }
    if ( flavAfter  == +16 || flavAfter  == -16 ) { indexAfter  = 2; }
    assert( indexBefore != -1 && indexAfter != -1 ); // This checks that we got PDG codes that are actually neutrinos.
    */
    // Probability matrix:
    //std::cout << "E = " << E << ", isAnti = " << isAnti << std::endl;
    //std::cout << this->fCurrProbs[0][0] << " " << this->fCurrProbs[0][1] << " " << this->fCurrProbs[0][2] << std::endl
    //          << this->fCurrProbs[1][0] << " " << this->fCurrProbs[1][1] << " " << this->fCurrProbs[1][2] << std::endl
    //          << this->fCurrProbs[2][0] << " " << this->fCurrProbs[2][1] << " " << this->fCurrProbs[2][2] << std::endl;
    //std::cout << "Returning prob " << indexBefore << " -> " << indexAfter << ": " << this->fCurrProbs[indexBefore][indexAfter] << std::endl;
    
    //return this->fCurrProbs[indexBefore][indexAfter];
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
  
  template<typename T>
  void _OscCalcNuFast<T>::RecomputeEnergyIndependentTerms(void) {
    // NOvA doesn't use dmsq31 by convention, but we can compute it if we know dmsq32 and dmsq21.
    this->fPre.Dmsq31 = this->fDmsq32 + this->fDmsq21; // dmsq31 = msq3 - msq1 = (msq3-msq2) + (msq2-msq1) = dmsq32 + dmsq21.
    this->fPre.Dmsqee = this->fPre.Dmsq31 - this->fPre.s12sq * this->fDmsq21;
    this->fPre.c13sqxs12sq = (1-this->fPre.s13sq) * this->fPre.s12sq;
    this->fPre.c13sqxs23sq = (1-this->fPre.s13sq) * this->fPre.s23sq;
    this->fPre.c12sqxc23sq = (1-this->fPre.s12sq) * (1-this->fPre.s23sq);
    this->fPre.s13sqxs12sqxs23sq = this->fPre.s13sq * this->fPre.s12sq * this->fPre.s23sq;
    this->fPre.Jrr = sqrt(this->fPre.c12sqxc23sq * this->fPre.s13sqxs12sqxs23sq);
    this->fPre.Jmatter_first = 8 * this->fPre.Jrr * (1-this->fPre.s13sq) * this->fPre.sind 
                                 * this->fDmsq21 * this->fPre.Dmsq31 * (this->fPre.Dmsq31 - this->fDmsq21);
    this->fPre.Um2sq_first = this->fPre.c12sqxc23sq + this->fPre.s13sqxs12sqxs23sq - 2 * this->fPre.Jrr * this->fPre.cosd;
    this->fPre.See = this->fDmsq21+this->fPre.Dmsq31 - this->fDmsq21 * this->fPre.c13sqxs12sq - this->fPre.Dmsq31 * this->fPre.s13sq;
    this->fPre.Tee = this->fDmsq21 * this->fPre.Dmsq31 * (1 - this->fPre.s13sq - this->fPre.c13sqxs12sq);
    ClearProbCaches(); // The caches were invalidated if this step was necessary.
    this->fIsDirty = false;
  }
  
  template<typename T> template<class VT, class KVT>
  VT _OscCalcNuFast<T>::RecomputeProbabilityMatrix(int from, int to, KVT E) {
    // Now we need to interface internally used NuFast variables with the OscCalc members.
    #define s12sq this->fPre.s12sq
    #define s13sq this->fPre.s13sq
    #define s23sq this->fPre.s23sq
    #define sind this->fPre.sind
    #define cosd this->fPre.cosd
    #define c12sq (1-s12sq)
    #define c13sq (1-s13sq)
    #define c23sq (1-s23sq)
    
    #define delta this->fdCP
    #define Dmsq21 this->fDmsq21
    #define L this->fL
    #define rho this->fRho
    #define Ye this->fYe
    #define N_Newton this->fNNewton
    #define probs_returned this->fCurrProbs
    
    // Aggressive precomputation. For terms that don't depend on energy, refer to the cached precomputations.
    #define Dmsq31 this->fPre.Dmsq31
    #define Dmsqee this->fPre.Dmsqee
    #define c13sqxs12sq this->fPre.c13sqxs12sq
    #define c13sqxs23sq this->fPre.c13sqxs23sq
    #define c12sqxc23sq this->fPre.c12sqxc23sq
    #define s13sqxs12sqxs23sq this->fPre.s13sqxs12sqxs23sq
    #define Jrr this->fPre.Jrr
    #define Jmatter_first this->fPre.Jmatter_first
    #define Um2sq_first this->fPre.Um2sq_first
    #define See this->fPre.See
    #define Tee this->fPre.Tee
    
//    VT /*c13sq, sind, cosd, Jrr,*/ Jmatter, /*Dmsqee,*/ Amatter;
//    VT Ue1sq, Ue2sq, Ue3sq, Um1sq, Um2sq, Um3sq, Ut1sq, Ut2sq, Ut3sq;
//    VT A, B, C;
//    VT /*See, Tee,*/ Smm, Tmm;
//   VT xmat, lambda2, lambda3, Dlambda21, Dlambda31, Dlambda32;
//    VT Xp2, Xp3, PiDlambdaInv/*, tmp*/;
//    VT Lover4E, D21, D32;
//    VT sinD21, sinD31, sinD32;
//    VT sinsqD21_2, sinsqD31_2, sinsqD32_2, triple_sin;
//    VT Pme_CPC, Pme_CPV, Pmm, Pee;
    
    // --------------------------------------------------------------------- //
    // First calculate useful simple functions of the oscillation parameters //
    // --------------------------------------------------------------------- //
    
    const VT Amatter = (Ye*rho*YerhoE2a) * E;
    
    // calculate A, B, C, See, Tee, and part of Tmm
    const VT C = Amatter * Tee;
    const VT A = Dmsq21+Dmsq31 + Amatter;
    
    // ---------------------------------- //
    // Get lambda3 from lambda+ of MP/DMP //
    // ---------------------------------- //
    const VT xmat = Amatter / Dmsqee;
    VT lambda3 = Dmsq31 + (0.5 * Dmsqee) * (xmat - 1 + sqrt((1-xmat) * (1-xmat) + (4*s13sq) * xmat));
    
    // ---------------------------------------------------------------------------- //
    // Newton iterations to improve lambda3 arbitrarily, if needed, (B needed here) //
    // ---------------------------------------------------------------------------- //
    const VT B = Dmsq21*Dmsq31 + Amatter * See; // B is only needed for N_Newton >= 1
    for (int i = 0; i < N_Newton; i++) {
      lambda3 = (lambda3 * lambda3 * (lambda3 + lambda3 - A) + C) / (lambda3 * (2 * (lambda3 - A) + lambda3) + B); // this strange form prefers additions to multiplications
    }
    
    // ------------------- //
    // Get  Delta lambda's //
    // ------------------- //
    const VT Dlambda21 = sqrt((A - lambda3) * (A - lambda3) - 4 * C / lambda3);
    const VT lambda2 = 0.5 * (A - lambda3 + Dlambda21);
    const VT Dlambda32 = lambda3 - lambda2;
    const VT Dlambda31 = Dlambda32 + Dlambda21;
    
    // ----------------------- //
    // Use Rosetta for Veisq's //
    // ----------------------- //
    // denominators	  
    const VT PiDlambdaInv = 1 / (Dlambda31 * Dlambda32 * Dlambda21);
    const VT Xp3 = PiDlambdaInv * Dlambda21;
    const VT Xp2 = -PiDlambdaInv * Dlambda31;
    
    // numerators
    const VT Ue3sq = (lambda3 * (lambda3 - See) + Tee) * Xp3;
    const VT Ue2sq = (lambda2 * (lambda2 - See) + Tee) * Xp2;
    
    const VT Smm = A - Dmsq21 * Um2sq_first - Dmsq31 * c13sqxs23sq;
    const VT Tmm = Dmsq21*Dmsq31 * (1 - c13sqxs23sq - Um2sq_first) + Amatter * (See + Smm - A);
    
    const VT Um3sq = (lambda3 * (lambda3 - Smm) + Tmm) * Xp3;
    const VT Um2sq = (lambda2 * (lambda2 - Smm) + Tmm) * Xp2;
    
    // ------------- //
    // Use NHS for J //
    // ------------- //
    const VT Jmatter = Jmatter_first * PiDlambdaInv;
    
    // ----------------------- //
    // Get all elements of Usq //
    // ----------------------- //
    const VT Ue1sq = 1 - Ue3sq - Ue2sq;
    const VT Um1sq = 1 - Um3sq - Um2sq;
    
    const VT Ut3sq = 1 - Um3sq - Ue3sq;
    const VT Ut2sq = 1 - Um2sq - Ue2sq;
    const VT Ut1sq = 1 - Um1sq - Ue1sq;
    
    // ----------------------- //
    // Get the kinematic terms //
    // ----------------------- //
    const VT Lover4E = (eVsqkm_to_GeV_over4 * L) / E;
    
    const VT D21 = Dlambda21 * Lover4E;
    const VT D32 = Dlambda32 * Lover4E;
	  
    const VT sinD21 = sin(D21);
    const VT sinD31 = sin(D32 + D21);
    const VT sinD32 = sin(D32);
    
    const VT triple_sin = sinD21 * sinD31 * sinD32;
    
    const VT sinsqD21_2 = 2 * sinD21 * sinD21;
    const VT sinsqD31_2 = 2 * sinD31 * sinD31;
    const VT sinsqD32_2 = 2 * sinD32 * sinD32;
    
    // ------------------------------------------------------------------- //
    // Calculate the three necessary probabilities, separating CPC and CPV //
    // ------------------------------------------------------------------- //
    const VT Pme_CPC = (Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2
                     + (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2
	                   + (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2;
    const VT Pme_CPV = -Jmatter * triple_sin;
    
    const VT Pmm = 1 - 2 * (Um2sq * Um1sq * sinsqD21_2
                 + Um3sq * Um1sq * sinsqD31_2
                 + Um3sq * Um2sq * sinsqD32_2);
    
    const VT Pee = 1 - 2 * (Ue2sq * Ue1sq * sinsqD21_2
                 + Ue3sq * Ue1sq * sinsqD31_2
                 + Ue3sq * Ue2sq * sinsqD32_2);
    
    // ---------------------------- //
    // Assign all the probabilities //
    // ---------------------------- //
    /*
    probs_returned[0][0] = Pee;                                             // Pee
    probs_returned[0][1] = Pme_CPC - Pme_CPV;                               // Pem
    probs_returned[0][2] = 1 - Pee - probs_returned[0][1];                  // Pet
    
    probs_returned[1][0] = Pme_CPC + Pme_CPV;                               // Pme
    probs_returned[1][1] = Pmm;                                             // Pmm
    probs_returned[1][2] = 1 - probs_returned[1][0] - Pmm;                  // Pmt
    
    probs_returned[2][0] = 1 - Pee - probs_returned[1][0];                  // Pte
    probs_returned[2][1] = 1 - probs_returned[0][1] - Pmm;                  // Ptm
    probs_returned[2][2] = 1 - probs_returned[0][2] - probs_returned[1][2]; // Ptt
    */
    
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
    #undef L
    #undef rho
    #undef Ye
    #undef N_Newton
    #undef probs_returned
		
    #undef c13sqxs12sq
		#undef c13sqxs23sq
		#undef c12sqxc23sq
		#undef s13sqxs12sqxs23sq
    #undef Jrr
    #undef Jmatter
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

