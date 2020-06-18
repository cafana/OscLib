// ----------------------------------------------------------------------------
// Numerical diagonalization of 3x3 matrcies
// Copyright (C) 2006  Joachim Kopp
// ----------------------------------------------------------------------------
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
// ----------------------------------------------------------------------------
#ifndef __ZHEEVC3_H
#define __ZHEEVC3_H

#include <complex>

#include <stdio.h>
#include <math.h>
#include <complex>

// Constants
namespace
{
  const double M_SQRT3 = 1.73205080756887729352744634151;

#ifndef ZHE_SQR
    #define ZHE_SQR
  template <typename T>
  T SQR(const T& x) {return x * x; }

  template <typename T>
  T SQR_ABS(const std::complex<T> & z) { return SQR(z.real()) + SQR(z.imag()); }
#endif
}

// ----------------------------------------------------------------------------
template <typename T>
int zheevc3(std::complex<T> A[3][3], T w[3])
// ----------------------------------------------------------------------------
// Calculates the eigenvalues of a hermitian 3x3 matrix A using Cardano's
// analytical algorithm.
// Only the diagonal and upper triangular parts of A are accessed. The access
// is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The hermitian input matrix
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error
// ----------------------------------------------------------------------------
{
  T m, c1, c0;

  // Determine coefficients of characteristic poynomial. We write
  //       | a   d   f  |
  //  A =  | d*  b   e  |
  //       | f*  e*  c  |
  std::complex<T> de = A[0][1] * A[1][2];                            // d * e
  T dd = SQR_ABS(A[0][1]);                                  // d * conj(d)
  T ee = SQR_ABS(A[1][2]);                                  // e * conj(e)
  T ff = SQR_ABS(A[0][2]);                                  // f * conj(f)
  m  = real(A[0][0]) + real(A[1][1]) + real(A[2][2]);
  c1 = (real(A[0][0])*real(A[1][1])  // a*b + a*c + b*c - d*conj(d) - e*conj(e) - f*conj(f)
        + real(A[0][0])*real(A[2][2])
        + real(A[1][1])*real(A[2][2]))
       - (dd + ee + ff);
  c0 = real(A[2][2])*dd + real(A[0][0])*ee + real(A[1][1])*ff
       - real(A[0][0])*real(A[1][1])*real(A[2][2])
       - 2.0 * (real(A[0][2])*real(de) + imag(A[0][2])*imag(de));
  // c*d*conj(d) + a*e*conj(e) + b*f*conj(f) - a*b*c - 2*Re(conj(f)*d*e)

  T p, sqrt_p, q, c, s, phi;
  p = SQR(m) - 3.0*c1;
  q = m*(p - (3.0/2.0)*c1) - (27.0/2.0)*c0;
  sqrt_p = sqrt(fabs(p));

  phi = 27.0 * ( 0.25*SQR(c1)*(p - c1) + c0*(q + 27.0/4.0*c0));
  phi = (1.0/3.0) * atan2(sqrt(fabs(phi)), q);

  c = sqrt_p*cos(phi);
  s = (1.0/M_SQRT3)*sqrt_p*sin(phi);

  w[1]  = (1.0/3.0)*(m - c);
  w[2]  = w[1] + s;
  w[0]  = w[1] + c;
  w[1] -= s;

  return 0;
}


#endif
