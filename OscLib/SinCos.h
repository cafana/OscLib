#pragma once

// This file provides implementations of sincos() for all types and platforms

#include <cmath>


#ifdef OSCLIB_STAN
#include "OscLib/Stan.h"

// Stan doesn't provide sincos()
inline void sincos(const stan::math::var& x, stan::math::var* sx, stan::math::var* cx)
{
  *sx = sin(x);
  *cx = cos(x);
}
#endif


// Apple prefix it with underscores for some reason
#ifdef __APPLE_CC__
inline void sincos(double x, double* sx, double* cx)
{
  __sincos(x, sx, cx);
}
#endif


#include <Eigen/Eigen>

// Eigen also doesn't provide this
template<class T, class U> inline void
sincos(T& x,
       Eigen::Array<U, Eigen::Dynamic, 1>* sx,
       Eigen::Array<U, Eigen::Dynamic, 1>* cx)
{
  // Presumably this is faster than the commented out version below
  sx->resize(x.size());
  cx->resize(x.size());
  for(int i = 0; i < x.size(); ++i) sincos(x[i], &(*sx)[i], &(*cx)[i]);

  //  *sx = sin(x);
  //  *cx = cos(x);
}

