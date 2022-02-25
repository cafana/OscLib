// n.b. Stan sets up some type traits that need to be loaded before Eigen is.
// Since Eigen gets dragged in via IOscCalc.h we have to get Stan set up before
// that is included.
#ifdef OSCLIB_STAN
#include "OscLib/Stan.h"
#endif

#include "OscLib/IOscCalc.h"

namespace{
  // Most directions are implicitly defined
  template<class T, class U> T GetValAs(const U& x){return x;}

#ifdef OSCLIB_STAN
  // But not this one (since you don't want to do it by accident)
  template<> double GetValAs<double>(const stan::math::var& x){return x.val();}
#endif
}

#include <Eigen/Dense>

#include "TMD5.h"

#include <iostream>

namespace osc
{
  //---------------------------------------------------------------------------
  template<class T> Eigen::Array<T, Eigen::Dynamic, 1> _IOscCalc<T>::
  P(int flavBefore, int flavAfter, const std::vector<double>& E)
  {
    Eigen::Matrix<T,Eigen::Dynamic,1> ret(E.size());
    for(auto i = 0u; i < E.size(); i++) {
      ret(i) = this->P(flavBefore, flavAfter, E[i]);
    }
    return ret.array().isNaN().select(0, ret);
  }

  //---------------------------------------------------------------------------
  template<class T> Eigen::Array<T, Eigen::Dynamic, 1> _IOscCalc<T>::
  P(int flavBefore, int flavAfter, const Eigen::ArrayXd& E)
  {
    return P(flavBefore, flavAfter, std::vector<double>(&E[0], &E[0]+E.size()));
  }

  //---------------------------------------------------------------------------
  template<class T> _IOscCalc<T>* _NoOscillations<T>::Copy() const
  {
    return new _NoOscillations<T>;
  }
    
  //---------------------------------------------------------------------------
  template<class T> T _NoOscillations<T>::P(int from, int to, double /*E*/)
  {
    if(from == to || to == 0) return 1;
    return 0;
  }

  //---------------------------------------------------------------------------
  template<class T> TMD5* _NoOscillations<T>::GetParamsHash() const
  {
    // Always compare equal to self
    TMD5* ret = new TMD5;
    const char* txt = "NoOscillations";
    ret->Update((unsigned char*)txt, strlen(txt));
    ret->Final();
    return ret;
  }

  //---------------------------------------------------------------------------
  template<class T> void _NoOscillations<T>::Print(const std::string& prefix) const
  {
    std::cout << prefix << "No oscillations" << std::endl;
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

  //---------------------------------------------------------------------------
  template<class T> void _IOscCalcAdjustable<T>::Print(const std::string& prefix) const
  {
    std::cout << prefix << "dmsq21 = " << fDmsq21 << " eV^2\n"
              << prefix << "dmsq32 = " << fDmsq32 << " eV^2\n"
              << prefix << "th12 = " << fTh12 << "\n"
              << prefix << "th13 = " << fTh13 << "\n"
              << prefix << "th23 = " << fTh23 << "\n"
              << prefix << "dCP = " << fdCP << "\n"
              << prefix << "L = " << fL << " km\n"
              << prefix << "rho = " << fRho << " g/cm^3" << std::endl;
  }

  //---------------------------------------------------------------------------
  template <typename T, typename U>
  void CopyParams(const osc::_IOscCalcAdjustable<T> * inCalc,
                  osc::_IOscCalcAdjustable<U> * outCalc)
  {
    assert (inCalc);
    assert (outCalc);

    outCalc->SetL(inCalc->GetL());
    outCalc->SetRho(inCalc->GetRho());

    outCalc->SetdCP(GetValAs<U>(inCalc->GetdCP()));
    outCalc->SetDmsq21(GetValAs<U>(inCalc->GetDmsq21()));
    outCalc->SetDmsq32(GetValAs<U>(inCalc->GetDmsq32()));
    outCalc->SetTh12(GetValAs<U>(inCalc->GetTh12()));
    outCalc->SetTh13(GetValAs<U>(inCalc->GetTh13()));
    outCalc->SetTh23(GetValAs<U>(inCalc->GetTh23()));
  }

  //---------------------------------------------------------------------------
  template class _IOscCalcAdjustable<double>;
  template class _IOscCalc<double>;

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

  template void CopyParams(const osc::_IOscCalcAdjustable<double> * inCalc,
                           osc::_IOscCalcAdjustable<double> * outCalc);
}
