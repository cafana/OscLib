// std_isnan needs to precede IOscCalc
// because it injects an important Stan overload of std::isnan
// into the std namespace that need sto be there
// BEFORE the Eigen headers are seen (due to the '#pragma once').
// IOscCalc.h  #includes Eigen/Eigen, so /shrug
#ifdef OSCLIB_STAN
#ifndef DARWINBUILD
#include "stan/math/rev/core/std_isnan.hpp"
#endif
#endif

#include "OscLib/IOscCalc.h"

#ifdef OSCLIB_STAN
#include "Utilities/StanUtils.h"
#else
// Otherwise turn this stan piece into a no-op
namespace{
  namespace util{
    template<class T> T GetValAs(T x){return x;}
  }
}
#endif

#include <Eigen/Dense>
namespace osc
{
  template<class T>
  Eigen::Matrix<T,Eigen::Dynamic,1>
  _IOscCalc<T>::P(int flavBefore, int flavAfter, const std::vector<double> &E)
  {
    Eigen::Matrix<T,Eigen::Dynamic,1> ret(E.size());
    for(auto i = 0u; i < E.size(); i++) {
      ret(i) = util::GetValAs<T>(this->P(flavBefore, flavAfter, E[i]));
    }
    return ret.array().isNaN().select(0, ret);
  }

  //---------------------------------------------------------------------------
  template<class T> _IOscCalcAdjustable<T>::~_IOscCalcAdjustable()
  {
  }

  //---------------------------------------------------------------------------
  template <typename T>
  TMD5* _IOscCalcAdjustable<T>::GetParamsHashDefault(const std::string& txt) const
  {
    TMD5* ret = new TMD5;
    ret->Update((unsigned char*)txt.c_str(), txt.size());
    const int kNumParams = 8;
    T buf[kNumParams] = {fRho, fL, fDmsq21, fDmsq32,
                         fTh12, fTh13, fTh23, fdCP};
    ret->Update((unsigned char*)buf, sizeof(T)*kNumParams);
    ret->Final();
    return ret;
  }

  template <typename T, typename U>
  void CopyParams(const osc::_IOscCalcAdjustable<T> * inCalc,
                  osc::_IOscCalcAdjustable<U> * outCalc)
  {
    assert (inCalc);
    assert (outCalc);

    outCalc->SetL(inCalc->GetL());
    outCalc->SetRho(inCalc->GetRho());

    outCalc->SetdCP(util::GetValAs<U>(inCalc->GetdCP()));
    outCalc->SetDmsq21(util::GetValAs<U>(inCalc->GetDmsq21()));
    outCalc->SetDmsq32(util::GetValAs<U>(inCalc->GetDmsq32()));
    outCalc->SetTh12(util::GetValAs<U>(inCalc->GetTh12()));
    outCalc->SetTh13(util::GetValAs<U>(inCalc->GetTh13()));
    outCalc->SetTh23(util::GetValAs<U>(inCalc->GetTh23()));
  }

  //---------------------------------------------------------------------------
  template class _IOscCalcAdjustable<double>;
  template class _IOscCalc<double>;
#ifndef DARWINBUILD
#ifdef OSCLIB_STAN
  template class _IOscCalcAdjustable<stan::math::var>;
  template class _IOscCalc<stan::math::var>;
  template void CopyParams(const osc::_IOscCalcAdjustable<double> * inCalc,
                           osc::_IOscCalcAdjustable<stan::math::var> * outCalc);
  template void CopyParams(const osc::_IOscCalcAdjustable<stan::math::var> * inCalc,
                           osc::_IOscCalcAdjustable<double> * outCalc);
  template void CopyParams(const osc::_IOscCalcAdjustable<stan::math::var> * inCalc,
                           osc::_IOscCalcAdjustable<stan::math::var> * outCalc);
#endif
#endif

  template void CopyParams(const osc::_IOscCalcAdjustable<double> * inCalc,
                           osc::_IOscCalcAdjustable<double> * outCalc);
}
