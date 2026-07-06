#ifndef OSCCALCANALYTIC_H
#define OSCCALCANALYTIC_H

#include "OscLib/Cache.h"
#include "OscLib/IOscCalc.h"

#include <functional>

namespace Eigen
{
  // Seems like an oversight to me...
  template<class T> using ArrayX = Eigen::Array<T, Eigen::Dynamic, 1>;
}

namespace osc::analytic
{
  /// std::complex takes a lot of care with inf/nan which we don't want
  template<class T, class U = T> struct cmplx
  {
    cmplx(const T& r, const U& i) : re(r), im(i) {}

    template<class A, class B> cmplx(const cmplx<A, B>& x) : re(x.re), im(x.im) {}

    template<class A, class B> cmplx<T>& operator=(const cmplx<A, B>& x)
    {
      re = x.re;
      im = x.im;
      return *this;
    }

    inline __attribute__((always_inline)) auto norm() const {return re*re + im*im;}
    inline __attribute__((always_inline)) cmplx<T, U> conj() const {return cmplx(re, -im);}

    T re; U im;
  };

  template<class T> struct Eigenvalues
  {
    cmplx<T> sume, sumxe, sumxxe;
  };

  template<class T> struct Hermitian
  {
    Hermitian() : em({}, {}), et({}, {}), mt({}, {}) {}

    inline __attribute__((always_inline))
    Eigenvalues<T> GetEigenvalues();

    T          ee;   cmplx<T> em;   cmplx<T>  et;
    /*cmplx<T> me;*/ T        mm;   cmplx<T>  mt;
    /*cmplx<T> te;   cmplx<T> tm;*/ T         tt;
  };

  template<class T> class _OscCalc: public _IOscCalcAdjustable<T>,
                                    protected ProbCache<double, T>,
                                    protected ProbCache<Eigen::ArrayXd, Eigen::ArrayX<T>>
  {
  public:
    _OscCalc();
    virtual ~_OscCalc();
    using _IOscCalc<T>::P;

    virtual _IOscCalcAdjustable<T>* Copy() const override;

    virtual void SetL(double L) override;
    virtual void SetRho(double rho) override;
    virtual void SetDmsq21(const T& dmsq21) override;
    virtual void SetDmsq32(const T& dmsq32) override;
    virtual void SetTh12(const T& th12) override;
    virtual void SetTh13(const T& th13) override;
    virtual void SetTh23(const T& th23) override;
    virtual void SetdCP(const T& dCP) override;

    virtual T P(int from, int to, double E) override;
    virtual Eigen::ArrayX<T> P(int from, int to, const std::vector<double>& E) override;
    virtual Eigen::ArrayX<T> P(int from, int to, const Eigen::ArrayXd& E) override;

    virtual TMD5* GetParamsHash() const override;

    // NSI parameters (dimensionless, relative to A_CC matter potential).
    // Default = 0 so existing 3F code is unaffected.
    // Off-diagonal are complex: set real and imaginary parts separately.
    // For antineutrinos, conjugation of off-diagonals is handled internally.
    void SetEps_ee    (double v)            { fEps_ee     = v;  ClearProbCaches(); }
    void SetEps_mumu  (double v)            { fEps_mumu   = v;  ClearProbCaches(); }
    void SetEps_tautau(double v)            { fEps_tautau = v;  ClearProbCaches(); }
    void SetEps_emu   (double re, double im){ fEps_emu_re   = re; fEps_emu_im   = im; ClearProbCaches(); }
    void SetEps_etau  (double re, double im){ fEps_etau_re  = re; fEps_etau_im  = im; ClearProbCaches(); }
    void SetEps_mutau (double re, double im){ fEps_mutau_re = re; fEps_mutau_im = im; ClearProbCaches(); }

  protected:
    void ClearProbCaches()
    {
      ProbCache<double, T>::clear();
      ProbCache<Eigen::ArrayXd, Eigen::ArrayX<T>>::clear();
    }

    template<class VT, class KVT> VT _P(int from, int to, const KVT& E);

    bool fDirty12, fDirty13, fDirty23, fDirtyCP, fDirtyMasses;

    T s12, c12, s13, c13, s23, c23, sCP, cCP;

    /*T        Ue1;*/ T        Ue2; cmplx<T>  Ue3;
    /*cmplx<T> Um1;*/ cmplx<T> Um2; T         Um3;
    /*cmplx<T> Ut1;*/ cmplx<T> Ut2; T         Ut3;

    inline __attribute__((always_inline)) void UpdatePMNS();

    T Hee;            cmplx<T>  Hem; cmplx<T>  Het;
    /*cmplx<T> Hme;*/ T         Hmm; cmplx<T>  Hmt;
    /*cmplx<T> Hte; cmplx<T>  Htm;*/ T         Htt;

    inline __attribute__((always_inline)) void UpdateHamiltonian();
    inline __attribute__((always_inline)) double Hmat();

    // NSI epsilon parameters (all double; zero by default)
    double fEps_ee = 0, fEps_mumu = 0, fEps_tautau = 0;
    double fEps_emu_re = 0,   fEps_emu_im = 0;
    double fEps_etau_re = 0,  fEps_etau_im = 0;
    double fEps_mutau_re = 0, fEps_mutau_im = 0;

  private:
    _OscCalc(const _OscCalc&) = default;
    _OscCalc& operator=(const _OscCalc&) = default;
  };
} // end namespaces

namespace osc
{
  template<class T> using _OscCalcAnalytic = osc::analytic::_OscCalc<T>;
  using OscCalcAnalytic = _OscCalcAnalytic<double>;
}

#endif
