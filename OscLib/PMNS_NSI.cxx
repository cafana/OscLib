////////////////////////////////////////////////////////////////////////
// $Id: PMNS_NSI.cxx,v 1.2 2013/04/03 19:59:31 jcoelho Exp $
//
// Implementation of oscillations of neutrinos in matter in a
// three-neutrino framework with NSI.

//......................................................................
//
// Throughout I have taken:
//   - L to be the neutrino flight distance in km
//   - E to be the neutrino energy in GeV
//   - Dmij to be the differences between the mass-squares in eV^2
//   - Ne to be the electron number density in mole/cm^3
//   - theta12,theta23,theta13,deltaCP to be in radians
//   - NSI parameters are dimensionless
//
// This  class inherits from the PMNSOpt class
//
// joao.coelho@tufts.edu
////////////////////////////////////////////////////////////////////////

// n.b. Stan sets up some type traits that need to be loaded before Eigen is.
#ifdef OSCLIB_STAN
#include "OscLib/Stan.h"
#endif

#include "OscLib/PMNS_NSI.h"

// Just pull in all the cxx files from MatrixDecomp. This way we don't have
// another library to worry about building and linking everywhere.
//#include "OscLib/MatrixDecomp/zhetrd3.cxx"
//#include "OscLib/MatrixDecomp/zheevc3.cxx"
#include "OscLib/MatrixDecomp/zheevh3.h"
//#include "OscLib/MatrixDecomp/zheevq3.cxx"

#include <cstdlib>
#include <cassert>

using namespace osc;

//......................................................................

template <typename T>
_PMNS_NSI<T>::_PMNS_NSI()
{
  this->SetMix(0.,0.,0.,0.);
  this->SetDeltaMsqrs(0.,0.);
  this->SetNSI(0.,0.,0.,0.,0.,0.,0.,0.,0.);
  this->ResetToFlavour(1);
  this->fCachedNe = 0.0;
  this->fCachedE =  1.0;
  this->fCachedAnti = 1;
}

template <typename T>
_PMNS_NSI<T>::~_PMNS_NSI(){
}

//......................................................................

template <typename T>
void _PMNS_NSI<T>::SetNSI(const T& eps_ee,    const T& eps_emu,    const T& eps_etau,
                          const T& eps_mumu,  const T& eps_mutau,  const T& eps_tautau,
                          const T& delta_emu, const T& delta_etau, const T& delta_mutau)
{

  fEps_ee     = eps_ee;
  fEps_mumu   = eps_mumu;
  fEps_tautau = eps_tautau;
  fEps_emu    = eps_emu   * complex(cos(delta_emu) ,   sin(delta_emu));
  fEps_etau   = eps_etau  * complex(cos(delta_etau) ,  sin(delta_etau));
  fEps_mutau  = eps_mutau * complex(cos(delta_mutau) , sin(delta_mutau));

  fResetNSI = true;

}

//......................................................................
///
/// Solve the full Hamiltonian with Non-Standard Interactions for 
/// eigenvectors and eigenvalues.
///
template <typename T>
void _PMNS_NSI<T>::SolveHam(double E, double Ne, int anti)
{

  // Check if anything has changed before recalculating
  if(Ne!=this->fCachedNe || E!=this->fCachedE || anti!=this->fCachedAnti || !this->fBuiltHlv || fResetNSI){
    this->fCachedNe = Ne;
    this->fCachedE = E;
    this->fCachedAnti = anti;
    fResetNSI = false;
    this->BuildHlv();
  }
  else return;

  T lv = 2 * constants::kGeVToeV * E / this->fDm31;  // Osc. length in eV^-1
  double kr2GNe = constants::kMatterDensityToEffect * Ne; // Matter potential in eV

  // Finish build Hamiltonian in matter with dimension of eV
  complex A[3][3];
  for(int i=0;i<3;i++){
    A[i][i] = this->fHlv[i][i]/lv;
    for(int j=i+1;j<3;j++){
      if(anti>0) A[i][j] = this->fHlv[i][j]/lv;
      else       A[i][j] = conj(this->fHlv[i][j])/lv;
    }
  }
  if(anti>0){
    A[0][0] += kr2GNe * (1 + fEps_ee);
    A[0][1] += kr2GNe * fEps_emu;
    A[0][2] += kr2GNe * fEps_etau;
    A[1][1] += kr2GNe * fEps_mumu;
    A[1][2] += kr2GNe * fEps_mutau;
    A[2][2] += kr2GNe * fEps_tautau;
  }
  else{
    A[0][0] -= kr2GNe * (1 + fEps_ee);
    A[0][1] -= kr2GNe * conj(fEps_emu);
    A[0][2] -= kr2GNe * conj(fEps_etau);
    A[1][1] -= kr2GNe * fEps_mumu;
    A[1][2] -= kr2GNe * conj(fEps_mutau);
    A[2][2] -= kr2GNe * fEps_tautau;
  }

  // Solve Hamiltonian for eigensystem using the GLoBES method
  zheevh3(A,this->fEvec,this->fEval);

}

////////////////////////////////////////////////////////////////////////

template class osc::_PMNS_NSI<double>;

#ifdef OSCLIB_STAN
template class osc::_PMNS_NSI<stan::math::var>;
#endif
