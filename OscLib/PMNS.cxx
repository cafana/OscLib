////////////////////////////////////////////////////////////////////////
// $Id: PMNS.cxx,v 1.4 2012-09-20 21:42:52 greenc Exp $
//
// Implementation of oscillations of neutrinos in matter in an
// three-neutrino framework based on the following reference:
//
//......................................................................
//
// PHYSICAL REVIEW D       VOLUME 22, NUMBER 11          1 DECEMBER 1980
//
//             Matter effects on three-neutrino oscillation
//
//                      V. Barger and K. Whisnant
//    Physics Department, U. of Wisconsin, Madison, Wisconsin 53706
//
//                            S. Pakvasa
//   Physics Department, U. of Hawaii at Manoa, Honolulu, Hawaii 96822
//
//                         R.J.N. Phillips
//        Rutherford Laboratory, Chilton, Didcot, Oxon, England
//                    (Received 4 August 1980)
//
//                            22 2718
//                            --
//......................................................................
//
// Throughout I have taken:
//   - L to be the neutrino flight distance in km
//   - E to be the neutrino energy in GeV
//   - dmsqr to be the differences between the mass-squares in eV^2
//   - Ne to be the electron number density in mole/cm^3
//   - theta12,theta23,theta13,deltaCP to be in radians
//
// messier@indiana.edu
////////////////////////////////////////////////////////////////////////

// n.b. Stan sets up some type traits that need to be loaded before Eigen is.
// Since Eigen gets dragged in via IOscCalc.h we have to get Stan set up before
// that is included.
#ifdef OSCLIB_STAN
#include "OscLib/Stan.h"
#endif

#include "OscLib/Constants.h"
#include "OscLib/PMNS.h"

#include <cstdlib>
#include <iostream>
#include <cassert>

using namespace osc;


// Some useful complex numbers
template<typename T>
auto zero() { return std::complex<T>(0.0,0.0); }
template<typename T>
auto one()  { return std::complex<T>(1.0,0.0); }

// Unit conversion constants
static const double kK1 = constants::kkmTom / constants::kGeVToeV / constants::kInversemToeV / 2;
static const double kK2 = constants::kAvogadroConstant * pow(constants::kInversemToeV*100,3) / pow(constants::kGeVToeV,2);

//......................................................................

template <typename T>
_PMNS<T>::_PMNS()
{
  this->SetMix(0.,0.,0.,0.);
  this->SetDeltaMsqrs(0.,0.);
  this->Reset();
}

//......................................................................

template <typename T>
_PMNS<T>::_PMNS(const T& th12, const T& th23, const T& th13, const T& deltacp,
                const T& dms21, const T& dms32)
{
  this->SetMix(th12, th23, th13, deltacp);
  this->SetDeltaMsqrs(dms21, dms32);
  this->Reset();
}

//......................................................................

template <typename T>
void _PMNS<T>::DumpMatrix(const complex M[][3]) const
{
  std::cout 
    <<"| "<<M[0][0]<<"\t"<<M[0][1]<<"\t"<<M[0][2]<<" |\n"
    <<"| "<<M[1][0]<<"\t"<<M[1][1]<<"\t"<<M[1][2]<<" |\n"
    <<"| "<<M[2][0]<<"\t"<<M[2][1]<<"\t"<<M[2][2]<<" |\n"
    <<std::endl;
}

//......................................................................

template <typename T>
void _PMNS<T>::PrintMix() const { this->DumpMatrix(fU); }

//......................................................................

template <typename T>
void _PMNS<T>::SetMix(const T& th12, const T& th23, const T& th13, const T& deltacp)
{
  int i, j;
  T  s12, s23, s13, c12, c23, c13;
  complex idelta(0.0,deltacp);
  
  s12 = sin(th12);  s23 = sin(th23);  s13 = sin(th13);
  c12 = cos(th12);  c23 = cos(th23);  c13 = cos(th13);

  // explicitly initializing complex part
  // to circumvent problem described at top of this file
  fU[0][0] =  complex(c12*c13, 0);
  fU[0][1] =  complex(s12*c13, 0);
  fU[0][2] =  s13*exp(-idelta);
  
  fU[1][0] = -s12*c23-c12*s23*s13*exp(idelta);
  fU[1][1] =  c12*c23-s12*s23*s13*exp(idelta);
  fU[1][2] =  complex(s23*c13, 0);
  
  fU[2][0] =  s12*s23-c12*c23*s13*exp(idelta);
  fU[2][1] = -c12*s23-s12*c23*s13*exp(idelta);
  fU[2][2] =  complex(c23*c13, 0);
  
  // Compute derived forms of the mixing matrix
  for (i=0; i<3; ++i) {
    for (j=0; j<3; ++j) {
      fUtran[i][j] = fU[j][i];
      fUstar[i][j] = conj(fU[i][j]);
      fUdagg[i][j] = conj(fU[j][i]);
    }
  }
}

///.....................................................................
///
/// Initialize the mixing matrix using the older form referenced by
/// the Barger et al paper
///
/// \warning This should not be used except for testing. Use SetMix
/// above.
///
template <typename T>
void _PMNS<T>::SetMixBWCP(double th1, double th2, double th3, double d)
{
  int i, j;
  double s1, s2, s3, c1, c2, c3;
  std::complex<double> id(0.0,d);
  s1 = sin(th1);  s2 = sin(th2);  s3 = sin(th3);
  c1 = cos(th1);  c2 = cos(th2);  c3 = cos(th3);
  
  fU[0][0] =  c1;
  fU[0][1] =  s1*c3;
  fU[0][2] =  s1*s3;

  fU[1][0] = -s1*c2;
  fU[1][1] =  c1*c2*c3+s2*s3*exp(id);
  fU[1][2] =  c1*c2*s3-s2*c3*exp(id);

  fU[2][0] = -s1*s2;
  fU[2][1] =  c1*s2*c3-c2*s3*exp(id);
  fU[2][2] =  c1*s2*s3+c2*c3*exp(id);

  // Compute derived forms of the mixing matrix
  for (i=0; i<3; ++i) {
    for (j=0; j<3; ++j) {
      fUtran[i][j] = fU[j][i];
      fUstar[i][j] = conj(fU[i][j]);
      fUdagg[i][j] = conj(fU[j][i]);
    }
  }
}

//......................................................................

template <typename T>
void _PMNS<T>::PrintDeltaMsqrs() const
{
  std::cout 
    <<"|"<<fdmsqr[0][0]<<"\t"<<fdmsqr[0][1]<<"\t"<<fdmsqr[0][2]<<"|\n" 
    <<"|"<<fdmsqr[1][0]<<"\t"<<fdmsqr[1][1]<<"\t"<<fdmsqr[1][2]<<"|\n" 
    <<"|"<<fdmsqr[2][0]<<"\t"<<fdmsqr[2][1]<<"\t"<<fdmsqr[2][2]<<"|"
    << std::endl;
}

//......................................................................
///
/// Set the mass-splittings. These are m_2^2-m_1^2 and m_3^2-m_2^2 in
/// eV^2
///
template <typename T>
void _PMNS<T>::SetDeltaMsqrs(const T& dm21, const T& dm32)
{
  double eps = 5.0E-9;
  T msqr[3];
  
  msqr[0] = 0.0;
  msqr[1] = dm21;
  msqr[2] = dm21+dm32;
  
  // Degeneracies cause problems with diagonalization, so break them
  // ever so slightly
  if (dm21==0.0) {msqr[0] -= 0.5*eps; msqr[1] += 0.5*eps; }
  if (dm32==0.0) {msqr[1] -= 0.5*eps; msqr[2] += 0.5*eps; }

  // Assign the mass splittings fdmsqr_ij = msqr_i - msqr_j by
  // convention
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      // A somewhat subtle point: Barger et al. refer to the sign of
      // m1-m2 being positive which implies dm^2_12>0 and
      // dm^2_21<0. The labeling in more common use is to assume m3 is
      // heaviest such that dm_12<0 and dm_21>0. Rather than reverse
      // all the indices in all the equations, flip the signs here.
      fdmsqr[i][j] = -(msqr[i] - msqr[j]);
    }
  }
}

//......................................................................
///
/// Compute matrix multiplication A = B*C
///
template <typename T>
void _PMNS<T>::Multi(complex A[][3], const complex B[][3], const complex C[][3])
{
  int i, j, k;
  for (i=0; i<3; ++i) {
    for (j=0; j<3; ++j) {
      A[i][j] = zero<T>();
      for (k=0; k<3; ++k) {
	A[i][j] += B[i][k]*C[k][j];
      }
    }
  }
}

//......................................................................
///
/// Compute Equation 2, transition matrix for propagation through
/// vaccum, from Barger et al:
///
/// A(nua->nub) = sum_i U_ai exp(-1/2 i m_i^2 L/E) Udagg_ib
///
template <typename T>
void _PMNS<T>::EvalEqn2(complex       A[][3],
		   const complex U[][3],
		   const complex Udagg[][3],
		   const T  dmsqr[][3],
		   double L,
		   double E) 
{
  int a, b, i;
  for (a=0; a<3; ++a) {
    for (b=0; b<3; ++b) {
      A[a][b] = zero<T>();
      for (i=0; i<3; ++i) {
	complex phase(0.0,-kK1*dmsqr[i][0]*L/E);
	A[a][b] += U[a][i]*exp(phase)*Udagg[i][b];
      }
    }
  }
}

//......................................................................
///
/// Compute Eqn. 5, the matter 2*E*Hamiltonian in the mass basis
///
template <typename T>
void _PMNS<T>::EvalEqn5(complex       twoEH[][3],
		   const complex U[][3],
		   const complex Udagg[][3],
		   const T       dmsqr[][3],
		   double        E,
		   double        Gf,
		   double        Ne)
{
  int j, k;
  T         k2r2GNeE = kK2*2.0*M_SQRT2*Gf*Ne*(constants::kGeVToeV*E);
  for (k=0; k<3; ++k) {
    for (j=0; j<3; ++j) {
      twoEH[k][j] = zero<T>();
      if (k==j) twoEH[k][j] = complex(dmsqr[j][0], 0);
      twoEH[k][j] -= k2r2GNeE*U[0][j]*Udagg[k][0];
    }
  }
}

//......................................................................
///
/// Compute Equation 10, the transition matrix for propagation across
/// matter slab of constant density, from Barger et al:
///
/// A = U * X * U^(dagger) for neutrinos
///
template <typename T>
void _PMNS<T>::EvalEqn10(complex       A[][3],
		    const complex U[][3],
		    const complex X[][3],
		    const complex Udagg[][3])
{
  complex tmp[3][3];
  this->Multi(tmp, X, Udagg);
  this->Multi(A,   U, tmp);
}

//......................................................................
///
/// Evaluate Eqn. 11, the Lagragne formula for the matrix e(-iHt),
/// from Barger et al.
///
template <typename T>
void _PMNS<T>::EvalEqn11(complex X[][3],
		    double L,
		    double E, 
		    const complex twoEH[][3],
		    const T       Msqr[],
		    const T       dMsqr[][3])
{
  // The identity matrix
  static const double One[3][3] = {{1.,0.,0.},
				   {0.,1.,0.},
				   {0.,0.,1.}
  };
  int a, b, k;
  complex phase;
  complex EHM0[3][3];
  complex EHM1[3][3];
  complex EHM2[3][3];
  complex EHM21[3][3];
  complex EHM20[3][3];
  complex EHM10[3][3];

  // There are three matrices which can apper inside the product on
  // j!=k. Calculate them before hand
  for (a=0; a<3; ++a) {
    for (b=0; b<3; ++b) {
      EHM0[a][b] = twoEH[a][b]-Msqr[0]*One[a][b];
      EHM1[a][b] = twoEH[a][b]-Msqr[1]*One[a][b];
      EHM2[a][b] = twoEH[a][b]-Msqr[2]*One[a][b];
    }
  }
  this->Multi(EHM21,EHM2,EHM1);
  this->Multi(EHM20,EHM2,EHM0);
  this->Multi(EHM10,EHM1,EHM0);

  // Refer to Msqr_j as dMsqr[j][0] since only mass differences matter
  for (a=0; a<3; ++a) {
    for (b=0; b<3; ++b) {
      X[a][b] = zero<T>();
      for (k=0; k<3; ++k) {
	phase = exp(complex(0.0,-kK1*dMsqr[k][0]*L/E));
	if (k==0) {
	  X[a][b] += (EHM21[a][b]/(dMsqr[k][2]*dMsqr[k][1]))*phase;
	}
	else if (k==1) { 
	  X[a][b] += (EHM20[a][b]/(dMsqr[k][2]*dMsqr[k][0]))*phase;
	}
	else if (k==2) {
	  X[a][b] += (EHM10[a][b]/(dMsqr[k][1]*dMsqr[k][0]))*phase;
	} // Switch for product on j!=k
      } // Sum on k
    } // Loop on b
  } // Loop on a
}

//......................................................................
// 
// Find the matter eigenvalues Msqr given the variables found in
// Eqn.22. This is Eqn.21 from Barger et. al.
//
template <typename T>
void _PMNS<T>::EvalEqn21(T Msqr[],
                         const T& alpha,
                         const T& beta,
                         const T& gamma)
{
  static const auto k2PiOver3 = 2.0*M_PI/3.0;
  auto alpha2           = alpha*alpha;
  auto alpha3           = alpha*alpha2;
  auto alpha2Minus3beta = alpha2-3.0*beta;

  // Argument of the acos()
  auto arg =
    (2.0*alpha3 - 9.0*alpha*beta + 27.0*gamma)/
    (2.0*pow(alpha2Minus3beta,1.5));

  // Occasionally round off errors mean that arg wanders outside of
  // its allowed range. If its way off (1 part per billion), stop the
  // program. Otherwise, set it to its real value.
  if (fabs(arg)>1.0) {
    if (fabs(arg-1.0)>1.E-9) abort();
    if (arg<-1.0) arg = -1.00;
    if (arg>+1.0) arg = +1.00;
  }
  
  // The three roots, find the first by computing the acos() the
  // others are nearby
  auto theta0 = acos(arg)/3.0;
  auto theta1 = theta0-k2PiOver3;
  auto theta2 = theta0+k2PiOver3;
  
  // The multiplier and offset
  auto fac      = -2.0*sqrt(alpha2Minus3beta)/3.0;   // Factor in front of cos() terms
  auto constant = -alpha/3.0; // The constant offset m1^2 is irrelevant
  
  // And the eigenvalues themselves
  Msqr[0] = fac*cos(theta0) + constant;
  Msqr[1] = fac*cos(theta1) + constant;
  Msqr[2] = fac*cos(theta2) + constant;
}

///.....................................................................
///
/// Compute the values of the simplifying expressions alpha, beta, and
/// gamma. This is Eqn22 from the Barger et al paper
///
template <typename T>
void _PMNS<T>::EvalEqn22(T& alpha,
                         T& beta,
                         T& gamma,
                         double  E,
                         double  Gf,
                         double  Ne,
                         const T  dmsqr[][3],
                         const complex U[][3])
{
  // 2*sqrt(2)*Gf*Ne*E in units of eV^2
  auto k2r2EGNe = kK2*2.0*M_SQRT2*Gf*Ne*(constants::kGeVToeV*E);
  
  alpha = k2r2EGNe + dmsqr[0][1] + dmsqr[0][2];
  
  beta  =
    dmsqr[0][1]*dmsqr[0][2] + 
    k2r2EGNe*(dmsqr[0][1]*(1.0-norm(U[0][1])) + 
	      dmsqr[0][2]*(1.0-norm(U[0][2])));

  gamma = k2r2EGNe*dmsqr[0][1]*dmsqr[0][2]*norm(U[0][0]);
}

//......................................................................
///
/// Sort out the eigenvalues
///
template <typename T>
void _PMNS<T>::SortEigenvalues(T       dMsqr[][3],
                               const T dmsqr[][3],
                               const T MsqrVac[],
                               T       Msqr[])
{
  int i, j, k;
  T best;
  T MsqrTmp[3] = {-99.9E9,-99.9E9,-99.9E9};
  int flg[3] = {0,0,0};

  // Attempt to figure out which of the eigenvalues match between
  // dmsqr and MsqrVac
  for (i=0; i<3; ++i) {
    best =  1.E30;
    k    = -1;
    for (j=0; j<3; ++j) {
      using namespace std;
      auto delta = abs(MsqrVac[i] - dmsqr[j][0]);
      if (delta<best) { best = delta; k = j; }
    }
    if (best>1.E-9) abort();
    if (k<0 || k>2) abort();
    flg[k] = 1;
    MsqrTmp[i] = Msqr[k];
  }
  // Check that each eigenvalue is used
  for (i=0; i<3; ++i) if (flg[i]!=1) abort();
  
  for (i=0; i<3; ++i) {
    Msqr[i] = MsqrTmp[i];
    for (j=0; j<3; ++j) {
      dMsqr[i][j] = (MsqrTmp[i] - MsqrTmp[j]);
    }
  }
}

///.....................................................................
///
/// Update the transition matrix for a step across the vacuum
///
template <typename T>
void _PMNS<T>::PropVacuum(double L, double E, int anti)
{
  int i, j;
  complex A[3][3];
  complex Aout[3][3];

  if      (anti>0) this->EvalEqn2(A, fU,     fUdagg, fdmsqr, L, E);
  else if (anti<0) this->EvalEqn2(A, fUstar, fUtran, fdmsqr, L, E);
  else abort();
  this->Multi(Aout, A, fA);
  for (i=0; i<3; ++i) {
    for (j=0; j<3; ++j) {
      fA[i][j] = Aout[i][j];
    }
  }
}

///.....................................................................
///
/// Update the transition matrix for a step across a slab of matter
///
template <typename T>
void _PMNS<T>::PropMatter(double L, double E, double Ne, int anti)
{
//  static const double  Gf = 1.166371E-5; // G_F in units of GeV^-2
static const double  Gf = constants::kFermiConstant;
  int i, j;
  complex twoEH[3][3];
  complex X[3][3];
  T       Msqr[3];
  T       MsqrV[3];
  T       dMsqr[3][3];
  T       alpha,  beta,  gamma;
  T       alpha0, beta0, gamma0;
  complex A[3][3];
  complex Aout[3][3];
  
  // Find the transition matrix. The series of steps are to:
  if (anti>0) {
    // [1] Find the matter Hamiltonian in the mass basis...
    this->EvalEqn5(twoEH, fU, fUdagg, fdmsqr, E, Gf, Ne);

    // [2] Find the eigenvalues and sort them.
    this->EvalEqn22(alpha, beta, gamma, E, Gf, Ne, fdmsqr, fU);
    this->EvalEqn21(Msqr,  alpha, beta, gamma);

    // Repeat the above, but for vacuum (Ne=0.0) to sort out the order
    // of the eigenvalues
    this->EvalEqn22(alpha0, beta0, gamma0, E, 0.0, 0.0, fdmsqr, fU);
    this->EvalEqn21(MsqrV, alpha0, beta0, gamma0);
    this->SortEigenvalues(dMsqr, fdmsqr, MsqrV, Msqr);

    // [3] Evaluate the transition matrix
    this->EvalEqn11(X, L, E, twoEH, Msqr, dMsqr);
    this->EvalEqn10(A, fU, X, fUdagg);
  }
  else if (anti<0) {
    // As above, but make required substitutions for anti-neutrinos:
    // Gf=>-Gf, U=>Ustar, U^dagger=>U^dagger^*=U^T
    this->EvalEqn5(twoEH, fUstar, fUtran, fdmsqr, E, -Gf, Ne);
    this->EvalEqn22(alpha, beta, gamma, E, -Gf, Ne, fdmsqr, fUstar);
    this->EvalEqn21(Msqr,  alpha, beta, gamma);
    this->EvalEqn22(alpha0, beta0, gamma0, E, 0.0, 0.0, fdmsqr, fUstar);
    this->EvalEqn21(MsqrV, alpha0, beta0, gamma0);
    this->SortEigenvalues(dMsqr, fdmsqr, MsqrV, Msqr);
    this->EvalEqn11(X, L, E, twoEH, Msqr, dMsqr);
    this->EvalEqn10(A, fUstar, X, fUtran);
  }
  else abort();
  
  // [4] Apply the transition matrix to the matrix we started with...
  this->Multi(Aout, A, fA);
  for (i=0; i<3; ++i) {
    for (j=0; j<3; ++j) {
      fA[i][j] = Aout[i][j];
    }
  }
}

//......................................................................
///
/// Do several layers in a row. L and Ne must have the same length
///
template <typename T>
void _PMNS<T>::PropMatter(const std::list<double>& L,
                          double                   E,
                          const std::list<double>& Ne,
                          int anti)
{
  if (L.size()!=Ne.size()) abort();
  auto Li  (L.begin());
  auto Lend(L.end());
  auto Ni  (Ne.begin());
  for (; Li!=Lend; ++Li, ++Ni) {
    // For very low densities, use vacumm
    static const double kRhoCutoff = 1.0E-6; // Ne in moles/cm^3
    if (*Ni<kRhoCutoff) this->PropVacuum(*Li, E, anti);
    else                this->PropMatter(*Li, E, *Ni, anti);
  }
}

//......................................................................
///
/// Reset the transition matrix back to the identity matrix from where
/// it starts
///
template <typename T>
void _PMNS<T>::Reset()
{
  int i, j;
  for (i=0; i<3; ++i) {
    for (j=0; j<3; ++j) {
      if (i==j) fA[i][j] = one<T>();
      else      fA[i][j] = zero<T>();
    }
  }
}

//......................................................................
///
/// Compute oscillation probability from flavor i to flavor j
///
/// 0 = nue, 1 = numu, 2 = nutau
///
template <typename T>
T _PMNS<T>::P(int i, int j) const
{
  assert(i>=0 && i<3);
  assert(j>=0 && j<3);
  return norm(fA[i][j]);
}


////////////////////////////////////////////////////////////////////////

// explicit template instantiation for known types

template class osc::_PMNS<double>;

#ifdef OSCLIB_STAN
template class osc::_PMNS<stan::math::var>;
#endif
