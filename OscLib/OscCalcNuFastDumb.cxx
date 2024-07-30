
#ifdef OSCLIB_STAN
#include "OscLib/Stan.h"
#endif

#include "OscLib/OscCalcNuFastDumb.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

namespace osc {
  template<typename T>
  _OscCalcNuFastDumb<T>::_OscCalcNuFastDumb(void) :
    fYe(0.5), fNNewton(1), fDirty(true), fCurrEnergy(0)
    {}
  
  template<typename T>
  _OscCalcNuFastDumb<T>::_OscCalcNuFastDumb(const _OscCalcNuFastDumb<T>& other) :
    _IOscCalcAdjustable<T>(other), fYe(other.GetYe()), fNNewton(other.GetNNewton()), fDirty(true), fCurrEnergy(0)
    {}
  
  template<typename T>
  _OscCalcNuFastDumb<T>::~_OscCalcNuFastDumb(void) {}
  
  template<typename T>
  _OscCalcNuFastDumb<T>* _OscCalcNuFastDumb<T>::Copy(void) const { return new _OscCalcNuFastDumb<T>(*this); }
  
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
  // We instead use a function signature consistent with OscLib conventions.
  template<typename T>
  T _OscCalcNuFastDumb<T>::P(int flavBefore, int flavAfter, double E) {
    // First, make sure the PDG codes have the same sign.
    static bool hasWarnedAboutSignChange = false;
    if ( (flavBefore > 0 && flavAfter < 0) || (flavBefore < 0 && flavAfter > 0) ) {
      // They're asking to change neutrino sign? Warn and return no probability.
      if ( !hasWarnedAboutSignChange ) {
        std::cout << "NuFast warning: Requesting probability with change in sign (" << flavBefore << " -> " << flavAfter << ")" << std::endl;
        hasWarnedAboutSignChange = true;
      }
      return 0;
    }
    
    bool isAnti = flavBefore < 0;
    // Recompute the probability matrix if necessary (either a parameter has been updated or we used a different energy before).
    //std::cout << "E = " << E << ", currE = " << this->fCurrEnergy << ", diff = " << this->fCurrEnergy - (isAnti ? -E : E) << std::endl;
    if ( this->fDirty || std::abs(this->fCurrEnergy - (isAnti ? -E : E)) > std::numeric_limits<double>::epsilon() ) {
    //  std::cout << "Nufast info: Recomputing probability matrix..." << std::endl;
      this->RecomputeProbabilityMatrix(E,isAnti);
    }
    
    // We have good PDG codes. Now return the corresponding element of the probability matrix.
    unsigned int indexBefore = -1, indexAfter = -1;
    // I'd use std::unordered_map for pdg->index, but this should be faster.
    if ( flavBefore == +12 || flavBefore == -12 ) { indexBefore = 0; }
    if ( flavBefore == +14 || flavBefore == -14 ) { indexBefore = 1; }
    if ( flavBefore == +16 || flavBefore == -16 ) { indexBefore = 2; }
    if ( flavAfter  == +12 || flavAfter  == -12 ) { indexAfter  = 0; }
    if ( flavAfter  == +14 || flavAfter  == -14 ) { indexAfter  = 1; }
    if ( flavAfter  == +16 || flavAfter  == -16 ) { indexAfter  = 2; }
    assert( indexBefore != -1 && indexAfter != -1 ); // This checks that we got PDG codes that are actually neutrinos.
    
    // Probability matrix:
    //std::cout << "E = " << E << ", isAnti = " << isAnti << std::endl;
    //std::cout << this->fCurrProbs[0][0] << " " << this->fCurrProbs[0][1] << " " << this->fCurrProbs[0][2] << std::endl
    //          << this->fCurrProbs[1][0] << " " << this->fCurrProbs[1][1] << " " << this->fCurrProbs[1][2] << std::endl
    //          << this->fCurrProbs[2][0] << " " << this->fCurrProbs[2][1] << " " << this->fCurrProbs[2][2] << std::endl;
    //std::cout << "Returning prob " << indexBefore << " -> " << indexAfter << ": " << this->fCurrProbs[indexBefore][indexAfter] << std::endl;
    
    return this->fCurrProbs[indexBefore][indexAfter];
  }
  
  template<typename T>
  void _OscCalcNuFastDumb<T>::RecomputeProbabilityMatrix(double E, bool isAnti) {
    // Internally in NuFast, positive energy is for neutrinos and negative energy is for antineutrinos.
    if ( E < 0 ) {
      std::cerr << "NuFast error: RecomputeProbabilityMatrix should be passed E >= 0. Pass isAnti = true for antineutrinos." << std::endl;
      abort();
    }
    if ( isAnti ) { E = -E; }
    this->fCurrEnergy = E; // Save the energy to test if these results can be reused in the future.
    this->fDirty = false; // Mark us as up-to-date.
    
    // Now we need to interface internally used NuFast variables with the OscCalc members.
    // The following variables and definitions are designed to efficiently interface with NuFast without rewriting its guts.
    T s12sq = sin(this->fTh12); s12sq = s12sq*s12sq; // Fast square that doesn't compute sine twice.
    T s13sq = sin(this->fTh13); s13sq = s13sq*s13sq;
    T s23sq = sin(this->fTh23); s23sq = s23sq*s23sq;
    
    #define delta this->fdCP
    #define Dmsq21 this->fDmsq21
    // NOvA doesn't use dmsq31 by convention, but we can compute it if we know dmsq32 and dmsq21.
    T Dmsq31 = this->fDmsq32 + this->fDmsq21; // dmsq31 = msq3 - msq1 = (msq3-msq2) + (msq2-msq1) = dmsq32 + dmsq21.
    #define L this->fL
    #define rho this->fRho
    #define Ye this->fYe
    #define N_Newton this->fNNewton
    #define probs_returned this->fCurrProbs
    
    // Now here begins the (almost exclusively) untouched NuFast algorithm.
    T c13sq, sind, cosd, Jrr, Jmatter, Dmsqee, Amatter;
    T Ue1sq, Ue2sq, Ue3sq, Um1sq, Um2sq, Um3sq, Ut1sq, Ut2sq, Ut3sq;
    T A, B, C;
    T See, Tee, Smm, Tmm;
    T xmat, lambda2, lambda3, Dlambda21, Dlambda31, Dlambda32;
    T Xp2, Xp3, PiDlambdaInv, tmp;
    T Lover4E, D21, D32;
    T sinD21, sinD31, sinD32;
    T sinsqD21_2, sinsqD31_2, sinsqD32_2, triple_sin;
    T Pme_CPC, Pme_CPV, Pmm, Pee;
    
    // --------------------------------------------------------------------- //
    // First calculate useful simple functions of the oscillation parameters //
    // --------------------------------------------------------------------- //
    c13sq = 1 - s13sq;
    
    // Ueisq's
    Ue2sq = c13sq * s12sq;
    Ue3sq = s13sq;
    
    // Umisq's, Utisq's and Jvac	 
    Um3sq = c13sq * s23sq;
    // Um2sq and Ut2sq are used here as temporary variables, will be properly defined later	 
    Ut2sq = s13sq * s12sq * s23sq;
    Um2sq = (1 - s12sq) * (1 - s23sq);
    
    Jrr = sqrt(Um2sq * Ut2sq);
    sind = sin(delta);
    cosd = cos(delta);
    
    Um2sq = Um2sq + Ut2sq - 2 * Jrr * cosd;
    Jmatter = 8 * Jrr * c13sq * sind;
    Amatter = Ye * rho * E * YerhoE2a;
    Dmsqee = Dmsq31 - s12sq * Dmsq21;
    
    // calculate A, B, C, See, Tee, and part of Tmm
    A = Dmsq21 + Dmsq31; // temporary variable
    See = A - Dmsq21 * Ue2sq - Dmsq31 * Ue3sq;
    Tmm = Dmsq21 * Dmsq31; // using Tmm as a temporary variable	  
    Tee = Tmm * (1 - Ue3sq - Ue2sq);
    C = Amatter * Tee;
    A = A + Amatter;
    
    // ---------------------------------- //
    // Get lambda3 from lambda+ of MP/DMP //
    // ---------------------------------- //
    xmat = Amatter / Dmsqee;
    tmp = 1 - xmat;
    lambda3 = Dmsq31 + 0.5 * Dmsqee * (xmat - 1 + sqrt(tmp * tmp + 4 * s13sq * xmat));
    
    // ---------------------------------------------------------------------------- //
    // Newton iterations to improve lambda3 arbitrarily, if needed, (B needed here) //
    // ---------------------------------------------------------------------------- //
    B = Tmm + Amatter * See; // B is only needed for N_Newton >= 1
    for (int i = 0; i < N_Newton; i++) {
      lambda3 = (lambda3 * lambda3 * (lambda3 + lambda3 - A) + C) / (lambda3 * (2 * (lambda3 - A) + lambda3) + B); // this strange form prefers additions to multiplications
    }
    
    // ------------------- //
    // Get  Delta lambda's //
    // ------------------- //
    tmp = A - lambda3;
    Dlambda21 = sqrt(tmp * tmp - 4 * C / lambda3);
    lambda2 = 0.5 * (A - lambda3 + Dlambda21);
    Dlambda32 = lambda3 - lambda2;
    Dlambda31 = Dlambda32 + Dlambda21;
    
    // ----------------------- //
    // Use Rosetta for Veisq's //
    // ----------------------- //
    // denominators	  
    PiDlambdaInv = 1 / (Dlambda31 * Dlambda32 * Dlambda21);
    Xp3 = PiDlambdaInv * Dlambda21;
    Xp2 = -PiDlambdaInv * Dlambda31;
    
    // numerators
    Ue3sq = (lambda3 * (lambda3 - See) + Tee) * Xp3;
    Ue2sq = (lambda2 * (lambda2 - See) + Tee) * Xp2;
    
    Smm = A - Dmsq21 * Um2sq - Dmsq31 * Um3sq;
    Tmm = Tmm * (1 - Um3sq - Um2sq) + Amatter * (See + Smm - A);
    
    Um3sq = (lambda3 * (lambda3 - Smm) + Tmm) * Xp3;
    Um2sq = (lambda2 * (lambda2 - Smm) + Tmm) * Xp2;
    
    // ------------- //
    // Use NHS for J //
    // ------------- //
    Jmatter = Jmatter * Dmsq21 * Dmsq31 * (Dmsq31 - Dmsq21) * PiDlambdaInv;
    
    // ----------------------- //
    // Get all elements of Usq //
    // ----------------------- //
    Ue1sq = 1 - Ue3sq - Ue2sq;
    Um1sq = 1 - Um3sq - Um2sq;
    
    Ut3sq = 1 - Um3sq - Ue3sq;
    Ut2sq = 1 - Um2sq - Ue2sq;
    Ut1sq = 1 - Um1sq - Ue1sq;
    
    // ----------------------- //
    // Get the kinematic terms //
    // ----------------------- //
    Lover4E = eVsqkm_to_GeV_over4 * L / E;
    
    D21 = Dlambda21 * Lover4E;
    D32 = Dlambda32 * Lover4E;
	  
    sinD21 = sin(D21);
    sinD31 = sin(D32 + D21);
    sinD32 = sin(D32);
    
    triple_sin = sinD21 * sinD31 * sinD32;
    
    sinsqD21_2 = 2 * sinD21 * sinD21;
    sinsqD31_2 = 2 * sinD31 * sinD31;
    sinsqD32_2 = 2 * sinD32 * sinD32;
    
    // ------------------------------------------------------------------- //
    // Calculate the three necessary probabilities, separating CPC and CPV //
    // ------------------------------------------------------------------- //
    Pme_CPC = (Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2
            + (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2
	          + (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2;
    Pme_CPV = -Jmatter * triple_sin;
    
    Pmm = 1 - 2 * (Um2sq * Um1sq * sinsqD21_2
        + Um3sq * Um1sq * sinsqD31_2
        + Um3sq * Um2sq * sinsqD32_2);
    
    Pee = 1 - 2 * (Ue2sq * Ue1sq * sinsqD21_2
        + Ue3sq * Ue1sq * sinsqD31_2
        + Ue3sq * Ue2sq * sinsqD32_2);
    
    // ---------------------------- //
    // Assign all the probabilities //
    // ---------------------------- //
    probs_returned[0][0] = Pee;                                             // Pee
    probs_returned[0][1] = Pme_CPC - Pme_CPV;                               // Pem
    probs_returned[0][2] = 1 - Pee - probs_returned[0][1];                  // Pet
    
    probs_returned[1][0] = Pme_CPC + Pme_CPV;                               // Pme
    probs_returned[1][1] = Pmm;                                             // Pmm
    probs_returned[1][2] = 1 - probs_returned[1][0] - Pmm;                  // Pmt
    
    probs_returned[2][0] = 1 - Pee - probs_returned[1][0];                  // Pte
    probs_returned[2][1] = 1 - probs_returned[0][1] - Pmm;                  // Ptm
    probs_returned[2][2] = 1 - probs_returned[0][2] - probs_returned[1][2]; // Ptt
  }
  
  template<typename T>
  void _OscCalcNuFastDumb<T>::Print(const std::string& prefix) const {
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

template class osc::_OscCalcNuFastDumb<double>;
#ifdef OSCLIB_STAN
  template class osc::_OscCalcNuFastDumb<stan::math::var>;
#endif

