////////////////////////////////////////////////////////////////////////
// $Id: PMNSOpt.cxx,v 1.1 2013/01/19 16:09:57 jcoelho Exp $
//
// Implementation of oscillations of neutrinos in matter in a
// three-neutrino framework. Two optimizations are relevant:
// The construction of the Hamiltonian follows DocDB-XXXX (to be posted)
// The eigensystem determination is based on the following reference:
//
//......................................................................
//
// Int. J. Mod. Phys. C       VOLUME 19, NUMBER 03            MARCH 2008
//
//     Efficient numerical diagonalization of hermitian 3x3 matrices
//
//                            Joachim Kopp
//                  Max–Planck–Institut für Kernphysik
//             Postfach 10 39 80, 69029 Heidelberg, Germany
//                    (Received 19 October 2007)
//
//                                523
//......................................................................
//
// Throughout I have taken:
//   - L to be the neutrino flight distance in km
//   - E to be the neutrino energy in GeV
//   - Dmij to be the differences between the mass-squares in eV^2
//   - Ne to be the electron number density in mole/cm^3
//   - theta12,theta23,theta13,deltaCP to be in radians
//
// The code structure follows the implementation written by M. Messier
// in the PMNS class.
//
// joao.coelho@tufts.edu
////////////////////////////////////////////////////////////////////////

// Note from templating over stan::math::var
// -----------------------------------------
// BEWARE.  stan::math::var default constructor has a nullptr in its content,
// so if you create a std::complex<stan::math::var>, both real & imag components
// are initially nullptr.
// If you write
//    std::complex<stan::math::var> zero = 0.0;
// you will now have a proper value in the real part, but the imaginary part
// will still be null!
// That can cause segfaults when we try to perform std::complex overloaded operations
// (e.g. arithmetic) that assume that the imaginary part is defined.
// Herein we attempt to explicitly call the std::complex constructor with a pair of zeros
//    std::complex<stan::math::var> zero(0.0, 0.0);
// etc. to avoid this happening.
//
// There is one unavoidable edge case, however.
// The overload of std::abs() for std::complex<> internally uses the following test:
//      if (__s == _Tp())
//        return __s;
// where __s is either the real or imaginary part of the number,
// and _Tp is the type std::complex is templated over.
// Because the way stan::math::var's operator== is currently written,
// it tries to dereference the nullptr inside the anonymous stan::math::var
// returned by calling the constructor to compare to __s's value...
// causing a segfault.
// The workaround is to simply write your own magnitude function and use that.
// i.e.
//    template<typename T>
//    double Mag(const std::complex<T>& z) { return z.real() * z.real() + z.imag() * z.imag(); }
// This was done in the zhe***3 libraries included below.

#include "OscLib/func/PMNSOpt.h"

#include "OscLib/func/MatrixDecomp/zhetrd3.h"
#include "OscLib/func/MatrixDecomp/zheevc3.h"
#include "OscLib/func/MatrixDecomp/zheevh3.h"
#include "OscLib/func/MatrixDecomp/zheevq3.h"

#include <cstdlib>
#include <cassert>

using namespace osc;

//......................................................................

template <typename T>
_PMNSOpt<T>::_PMNSOpt()
{
  this->SetMix(0.,0.,0.,0.);
  this->SetDeltaMsqrs(0.,0.);
  this->ResetToFlavour(1);
  fCachedNe = 0.0;
  fCachedE =  1.0;
  fCachedAnti = 1;
  fBuiltHlv = false;
}


//......................................................................

template <typename T>
void _PMNSOpt<T>::SetMix(const T &th12, const T &th23, const T &th13, const T &deltacp)
{

  fTheta12 = th12;
  fTheta23 = th23;
  fTheta13 = th13;
  fDeltaCP = deltacp;

  fBuiltHlv = false;

}

//......................................................................
///
/// Set the mass-splittings. These are m_2^2-m_1^2, m_3^2-m_2^2
/// and m_3^2-m_1^2 in eV^2
///
template <typename T>
void _PMNSOpt<T>::SetDeltaMsqrs(const T &dm21, const T &dm32)
{

  fDm21 = dm21;
  fDm31 = dm32 + dm21;

  if(fDm31==0) fDm31 = 1.0e-12;

  fBuiltHlv = false;

}

//......................................................................
///
/// Build H*lv, where H is the Hamiltonian in vacuum on flavour basis
/// and lv is the oscillation length
///
/// This is a dimentionless hermitian matrix, so only the
/// upper triangular part needs to be filled
///
/// The construction of the Hamiltonian avoids computing terms that
/// are simply zero. This has a big impact in the computation time.
/// This construction is described in DocDB-XXXX (to be posted)
///
template <typename T>
void _PMNSOpt<T>::BuildHlv()
{

  // Check if anything changed
  if(fBuiltHlv) return;

  // Create temp variables
  T sij, cij, h00, h11, h01;
  complex expCP(0,0), h02(0,0), h12(0,0);

  // Hamiltonian in mass base. Only one entry is variable.
  h11 = fDm21 / fDm31;

  // Rotate over theta12
  sij = sin(fTheta12);
  cij = cos(fTheta12);

  // There are 3 non-zero entries after rephasing so that h22 = 0
  h00 = h11 * sij * sij - 1;
  h01 = h11 * sij * cij;
  h11 = h11 * cij * cij - 1;

  // Rotate over theta13 with deltaCP
  sij = sin(fTheta13);
  cij = cos(fTheta13);
  expCP = complex(cos(fDeltaCP), -sin(fDeltaCP));

  // There are 5 non-zero entries after rephasing so that h22 = 0
  h02 = (-h00 * sij * cij) * expCP;
  h12 = (-h01 * sij) * expCP;
  h11 -= h00 * sij * sij;
  h00 *= cij * cij  -  sij * sij;
  h01 *= cij;

  // Finally, rotate over theta23
  sij = sin(fTheta23);
  cij = cos(fTheta23);

  // Fill the Hamiltonian rephased so that h22 = -h11
  // explicit construction of complex vars to avoid problem noted at top of file
  fHlv[0][0] = complex(h00 - 0.5 * h11, 0);
  fHlv[1][1] = complex(0.5 * h11 * (cij * cij - sij * sij)  +  2 * real(h12) * cij * sij, 0);
  fHlv[2][2] = -fHlv[1][1];

  // these are all constructed from complex variables (h02 and h12) so they are ok
  fHlv[0][1] = h02 * sij  +  h01 * cij;
  fHlv[0][2] = h02 * cij  -  h01 * sij;
  fHlv[1][2] = h12 - (h11 * cij + 2 * real(h12) * sij) * sij;

  // Tag as built
  fBuiltHlv = true;

}

//......................................................................
///
/// Solve the full Hamiltonian for eigenvectors and eigenvalues
/// This is using a method from the GLoBES software available at
/// http://www.mpi-hd.mpg.de/personalhomes/globes/3x3/
/// We should cite them accordingly
///
template <typename T>
void _PMNSOpt<T>::SolveHam(double E, double Ne, int anti)
{

  // Check if anything has changed before recalculating
  if(Ne!=fCachedNe || E!=fCachedE || anti!=fCachedAnti || !fBuiltHlv ){
    fCachedNe = Ne;
    fCachedE = E;
    fCachedAnti = anti;
    this->BuildHlv();
  }
  else return;

  auto lv = 2 * kGeV2eV*E / fDm31;  // Osc. length in eV^-1
  auto kr2GNe = kK2*M_SQRT2*kGf*Ne; // Matter potential in eV

  // Finish building Hamiltonian in matter with dimension of eV
  complex A[3][3];
  for(int i=0;i<3;i++){
    A[i][i] = fHlv[i][i]/lv;
    for(int j=i+1;j<3;j++){
      if(anti>0) A[i][j] = fHlv[i][j]/lv;
      else       A[i][j] = conj(fHlv[i][j])/lv;
    }
  }
  if(anti>0) A[0][0] += kr2GNe;
  else       A[0][0] -= kr2GNe;

  // Solve Hamiltonian for eigensystem using the GLoBES method
  zheevh3(A,fEvec,fEval);

}

///.....................................................................
///
/// Propagate the current neutrino state over a distance L in km
/// with an energy E in GeV through constant matter of density
/// Ne in mole/cm^3.
/// @param anti - +1 = neutrino case, -1 = anti-neutrino case
///
template <typename T>
void _PMNSOpt<T>::PropMatter(double L, double E, double Ne, int anti)
{

  // Solve Hamiltonian
  this->SolveHam(E, Ne, anti);

  // Store coefficients of propagation eigenstates
  complex nuComp[3];

  for(int i=0;i<3;i++){
    nuComp[i] = complex(0,0);
    for(int j=0;j<3;j++){
      nuComp[i] += fNuState[j] * conj(fEvec[j][i]);
    }
  }

  for(int i=0;i<3;i++)fNuState[i] = complex(0,0);

  // Propagate neutrino state
  for(int j=0;j<3;j++){
    auto s = sin(-fEval[j] * kKm2eV * L);
    auto c = cos(-fEval[j] * kKm2eV * L);

    complex jPart = complex(c, s) * nuComp[j];
    for(int i=0;i<3;i++){
      fNuState[i] += jPart * fEvec[i][j];
    }
  }

}

//......................................................................
///
/// Do several layers in a row. L and Ne must have the same length
///
template <typename T>
void _PMNSOpt<T>::PropMatter(const std::list<double>& L,
                           double                   E,
                           const std::list<double>& Ne,
                           int anti)
{
  if (L.size()!=Ne.size()) abort();
  auto Li  = L.begin();
  auto Lend = L.end();
  auto Ni = Ne.begin();
  for (; Li!=Lend; ++Li, ++Ni) {
    // For very low densities, use vacumm
    static const double kRhoCutoff = 1.0E-6; // Ne in moles/cm^3
    if (*Ni<kRhoCutoff) this->PropVacuum(*Li, E, anti);
    else                this->PropMatter(*Li, E, *Ni, anti);
  }
}


///.....................................................................
///
/// We know the vacuum eigensystem, so just write it explicitly
/// The eigenvalues depend on energy, so E needs to be provided in GeV
/// @param anti - +1 = neutrino case, -1 = anti-neutrino case
///
template <typename T>
void _PMNSOpt<T>::SetVacuumEigensystem(double E, int anti)
{

  T       s12, s23, s13, c12, c23, c13;
  complex expidelta(cos(fDeltaCP), anti * sin(fDeltaCP));

  s12 = sin(fTheta12);  s23 = sin(fTheta23);  s13 = sin(fTheta13);
  c12 = cos(fTheta12);  c23 = cos(fTheta23);  c13 = cos(fTheta13);

  fEvec[0][0] =  complex(c12*c13, 0);
  fEvec[0][1] =  complex(s12*c13, 0);
  fEvec[0][2] =  s13*conj(expidelta);

  fEvec[1][0] = -s12*c23-c12*s23*s13*expidelta;
  fEvec[1][1] =  c12*c23-s12*s23*s13*expidelta;
  fEvec[1][2] =  complex(s23*c13, 0);

  fEvec[2][0] =  s12*s23-c12*c23*s13*expidelta;
  fEvec[2][1] = -c12*s23-s12*c23*s13*expidelta;
  fEvec[2][2] =  complex(c23*c13, 0);

  fEval[0] = 0;
  fEval[1] = fDm21 / (2 * kGeV2eV*E);
  fEval[2] = fDm31 / (2 * kGeV2eV*E);

}

///.....................................................................
///
/// Propagate the current neutrino state over a distance L in km
/// with an energy E in GeV through vacuum
/// @param anti - +1 = neutrino case, -1 = anti-neutrino case
///
template <typename T>
void _PMNSOpt<T>::PropVacuum(double L, double E, int anti)
{

  this->SetVacuumEigensystem(E, anti);

  complex nuComp[3];

  for(int i=0;i<3;i++){
    nuComp[i] = 0;
    for(int j=0;j<3;j++){
      nuComp[i] += fNuState[j] * conj(fEvec[j][i]);
    }
  }

  const T km2EvL = kKm2eV*L;  // needed for the templated multiplication below to work
  for(int i=0;i<3;i++){
    fNuState[i] = 0;
    for(int j=0;j<3;j++){
      complex iEval(0.0,fEval[j]);
      fNuState[i] +=  exp(-iEval * km2EvL) * nuComp[j] * fEvec[i][j];
    }
  }

}

//......................................................................
///
/// Reset the neutrino state back to a pure flavour where
/// it starts
///
template <typename T>
void _PMNSOpt<T>::ResetToFlavour(int flv)
{
  int i;
  for (i=0; i<3; ++i){
    if (i==flv) fNuState[i] = complex(1, 0);
    else        fNuState[i] = complex(0, 0);
  }
}

//......................................................................
///
/// Compute oscillation probability of flavour flv
///
/// 0 = nue, 1 = numu, 2 = nutau
///
template <typename T>
T _PMNSOpt<T>::P(int flv) const
{
  assert(flv>=0 && flv<3);
  return norm(fNuState[flv]);
}

////////////////////////////////////////////////////////////////////////
// manually instantiate templates for those cases we know about.

template class osc::_PMNSOpt<double>;

#ifndef DARWINBUILD
#include "Utilities/func/Stan.h"
  template class osc::_PMNSOpt<stan::math::var>;
#endif

