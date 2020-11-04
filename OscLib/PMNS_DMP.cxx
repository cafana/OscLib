// std_isnan needs to precede the header
// because it injects an important Stan overload of std::isnan
// into the std namespace that needs to be there
// BEFORE the Eigen headers are seen (due to the '#pragma once').
// PMNS_DMP.h  #includes Eigen/Eigen, so /shrug
#ifdef OSCLIB_STAN
#include "stan/math/rev/core/std_isnan.hpp"
#endif

#include "OscLib/PMNS_DMP.h"

namespace osc 
{


  template <typename T>
  void _PMNS_DMP<T>::SetParameters(_OscParameters<T> const & params) {

    _s212 = sin(2 * params.th12);
    _s213 = sin(2 * params.th13);
    _c212 = cos(2 * params.th12);
    _c213 = cos(2 * params.th13);

    _Dmsq21 = params.dmsq21;
    _Dmsq31 = params.Dmsq31();
    _cd = cos(params.deltacp);
    _sd = sin(params.deltacp);

    T s12 = sin(params.th12);
    _Dmsqee = _Dmsq21 + params.dmsq32 - s12*s12 * _Dmsq21;

    _s13 = sin(params.th13);
    _s23 = sin(params.th23);

    _s13sq = _s13*_s13;
    _s23sq = _s23*_s23;

    _c13sq = 1 - _s13sq;
    _c23sq = 1 - _s23sq;

    _c13 = sqrt(_c13sq);
    _c23 = sqrt(_c23sq);

  }

  template <typename T>
  void _PMNS_DMP<T>::SetParameters(ParamArray const & params) {
    _OscParameters<T> p = {params[0], params[1], params[2], params[3], params[4], params[5]};
    this->SetParameters(p);
  }

  template <typename T>
  typename _PMNS_DMP<T>::ProbArray _PMNS_DMP<T>::P(ParamArray const & params) {
    _OscParameters<T> p = {params[0], params[1], params[2], params[3], params[4], params[5]};
    return this->P(p);
  }

  template <typename T>
  typename _PMNS_DMP<T>::ProbArray _PMNS_DMP<T>::P(T dmsq21, T dmsq32, T th12, T th13, T th23, T delta) {
    _OscParameters<T> p = { dmsq21,  dmsq32,  th12,  th13,  th23,  delta};
    return this->P(p);
  }

  template <typename T>
  typename _PMNS_DMP<T>::ProbArray _PMNS_DMP<T>::P(_OscParameters<T> const & params) {

    this->SetParameters(params);

    const ParamArray DMSQEEA = _Dmsqee * sqrt(square(_c213*_ONE - 1./_Dmsqee * _AA) + _s213*_s213*_ONE);

    const ParamArray C2PHI   = (_Dmsqee * _c213 * _ONE - _AA) / DMSQEEA;
    const ParamArray A12     = 0.5 * (_AA + _Dmsqee*_ONE - DMSQEEA);

    const ParamArray CPHISQ   = 0.5 * (_ONE + C2PHI);
    const ParamArray SPHISQ   = 0.5 * (_ONE - C2PHI);
    const ParamArray S2PHI    = _s213 * _Dmsqee * _ONE / DMSQEEA;
    const ParamArray CPHI13SQ = CPHISQ * _c13sq + SPHISQ * _s13sq + S2PHI * _c13 * _s13;

    const ParamArray DL21 = _Dmsq21 * sqrt(square(_c212*_ONE - A12 / _Dmsq21) + CPHI13SQ * _s212 * _s212);
    const ParamArray C2PSI = (_Dmsq21 * _c212 * _ONE - A12) / DL21;
    const ParamArray CPSISQ = 0.5 * (_ONE + C2PSI);
    const ParamArray SPSISQ = 0.5 * (_ONE - C2PSI);

    const ParamArray DL31 = _Dmsq31*_ONE + 0.25 * _AA + 0.5 * (DL21 - _Dmsq21*_ONE) + 0.75 * (DMSQEEA - _Dmsqee*_ONE);

    const ParamArray JRRM = _s23 * _c23 * sqrt(SPHISQ * CPSISQ * SPSISQ);
    const ParamArray JRM = JRRM * CPHISQ;
    const ParamArray  D = -JRM * _sd;

    ProbArray DD(_ENERGIES.rows(), 3);
    DD.col(0) = DL21 * _L4E;
    DD.col(1) = DL31 * _L4E;
    DD.col(2) = DD.col(1)-DD.col(0);
    const ProbArray SDD = DD.sin();

    ProbArray PP(_ENERGIES.rows(), 9);

    // e -> e
    ParamArray   C31 = -1 * CPHISQ * SPHISQ * CPSISQ;
    ParamArray   C32 = -1 * CPHISQ * SPHISQ * SPSISQ;
    ParamArray   C21 = -1 * CPHISQ.square() * SPSISQ * CPSISQ;
    PP.col(0) =  _ONE + 4 * (C31 * SDD.col(1).square() + C32 * SDD.col(2).square() + C21 * SDD.col(0).square());

    // mu -> mu
    C31 = -1 * CPHISQ * _s23sq * (_c23sq * SPSISQ + _s23sq * SPHISQ * CPSISQ) \
         - 2 * _s23sq * JRM * _cd;
    C32 = -1 * CPHISQ * _s23sq * (_c23sq * CPSISQ + _s23sq * SPHISQ * SPSISQ) \
         + 2 * _s23sq * JRM * _cd;
    C21 = -1 * (_c23sq * CPSISQ + _s23sq * SPHISQ * SPSISQ) * (_c23sq * SPSISQ + _s23sq * SPHISQ * CPSISQ) \
         - 2 * (_c23sq - SPHISQ * _s23sq) * C2PSI * JRRM * _cd + square(2 * JRRM * _cd);
    PP.col(4) =  _ONE + 4 * (C31 * SDD.col(1).square() + C32 * SDD.col(2).square() + C21 * SDD.col(0).square());

    C31 = _s23sq          * SPHISQ * CPHISQ * CPSISQ + JRM * _cd;
    C32 = _s23sq          * SPHISQ * CPHISQ * SPSISQ - JRM * _cd;
    C21 = CPHISQ * SPSISQ * CPSISQ * (_c23sq* _ONE - SPHISQ * _s23sq) + JRM * _cd * C2PSI;
    // Pemu
    PP.col(3) = 4 * (C31 * SDD.col(1).square() + C32 * SDD.col(2).square() + C21 * SDD.col(0).square()) - 8 * D * SDD.col(0) * SDD.col(1) * SDD.col(2);
    // Pmue
    PP.col(1) = 4 * (C31 * SDD.col(1).square() + C32 * SDD.col(2).square() + C21 * SDD.col(0).square()) + 8 * D * SDD.col(0) * SDD.col(1) * SDD.col(2);

    // e->tau = 1 - ee - emu
    PP.col(6) = _ONE - PP.col(0) - PP.col(3);
    // tau->e = 1 - ee - mue
    PP.col(2) = _ONE - PP.col(0) - PP.col(1);
    // mu->tau = 1 - mue - mumu
    PP.col(7) = _ONE - PP.col(1) - PP.col(4);
    // tau->mu = 1 - emu - mumu
    PP.col(5) = _ONE - PP.col(3) - PP.col(4);
    // tau->tau = 1 - etau - mutau
    PP.col(8) = _ONE - PP.col(6) - PP.col(7);

    return PP.isNaN().select(0,PP);
  }
}


////////////////////////////////////////////////////////////////////////
// manually instantiate templates for those cases we know about.

template class osc::_PMNS_DMP<double>;

#ifdef OSCLIB_STAN
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "stan/math/rev/scal.hpp"
template class osc::_PMNS_DMP<stan::math::var>;
#endif
