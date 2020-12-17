#ifndef OSCCALCANALYTIC_H
#define OSCCALCANALYTIC_H

#include "OscLib/IOscCalc.h"

#include <functional>
#include <unordered_map>

namespace Eigen
{
  // Seems like an oversight to me...
  template<class T> using ArrayX = Eigen::Array<T, Eigen::Dynamic, 1>;
}

// We want to put ArrayXd into an unordered_map, so define hash and equality
namespace std
{
  template<> struct hash<Eigen::ArrayXd>
  {
    size_t operator()(const Eigen::ArrayXd& x) const
    {
      // Adapted from the boost `hash_combine` function
      // http://www.boost.org/doc/libs/1_55_0/doc/html/hash/reference.html#boost.hash_combine

      size_t seed = 0;
      for(int i = 0; i < x.size(); ++i) {
        seed ^= std::hash<double>()(x[i] + 0x9e3779b9 + (seed << 6) + (seed >> 2));
      }
      return seed;
    }
  };

  template<> struct equal_to<Eigen::ArrayXd>
  {
    bool operator()(const Eigen::ArrayXd& a, const Eigen::ArrayXd& b) const
    {
      return (a == b).all();
    }
  };
}

namespace osc::analytic
{
  /// std::complex takes a lot of care with inf/nan which we don't want
  template<class T, class U = T> struct cmplx
  {
    cmplx() {} // TODO TODO had to add this for Eigen::Array. Are we happy?

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
    //      std::array<double, 3> xs; ///< eigenvalues
    //      std::array<cmplx,  3> expixs; ///< exp(i*x)/(3*x^2+b*x+c) for each x
    // Turns out we only need these summary expressions of the eigenvalues
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

    // TODO gives some fancy Eigen type, which we might want, but leads to a lot more auto
    //    template<class U> Hermitian<decltype(T{}*U{})> operator*(U k) const {return {ee*k, em*k, et*k, mm*k, mt*k, tt*k};}
    template<class U> Hermitian<U> operator*(const U& k) const {return Hermitian<U>(ee*k, em*k, et*k, mm*k, mt*k, tt*k);}

    // TODO why can't my operator* see this? because it's a different template type?
    //  protected:
    Hermitian(const T& _ee, const cmplx<T>& _em, const cmplx<T>& _et, const T& _mm, const cmplx<T>& _mt, const T& _tt)
      : ee(_ee), em(_em), et(_et), mm(_mm), mt(_mt), tt(_tt)
    {
    }
  };

  template<class T> struct PMNS
  {
    PMNS() : e3({}, {}), m2({}, {}), t2({}, {}) {}

    /*T        e1;*/ T        e2; cmplx<T> e3;
    /*cmplx<T> m1;*/ cmplx<T> m2; T        m3;
    /*cmplx<T> t1;*/ cmplx<T> t2; T        t3;
  };

  template<class KT, class VT> class AmpCache : public std::unordered_map<KT, Eigen::Array<cmplx<VT>, 3, 3>> {};

  template<class T> class _OscCalc: public _IOscCalcAdjustable<T>,
                                    protected AmpCache<double, T>,
                                    protected AmpCache<Eigen::ArrayXd, Eigen::ArrayX<T>>
  {
  public:
    _OscCalc();
    virtual ~_OscCalc();
    using _IOscCalc<T>::P;

    virtual _IOscCalcAdjustable<T>* Copy() const override;

    // Baseline in km
    virtual void SetL(double L) override;
    // Density in g/cm^3
    virtual void SetRho(double rho) override;
    // in eV^2
    virtual void SetDmsq21(const T& dmsq21) override;
    // This is a signed quantity, use a negative value for inverted hierarchy
    virtual void SetDmsq32(const T& dmsq32) override;
    // In radians
    virtual void SetTh12(const T& th12) override;
    virtual void SetTh13(const T& th13) override;
    virtual void SetTh23(const T& th23) override;
    virtual void SetdCP(const T& dCP) override;

    virtual T P(int from, int to, double E) override;
    virtual Eigen::ArrayX<T> P(int from, int to, const std::vector<double>& E) override;
    virtual Eigen::ArrayX<T> P(int from, int to, const Eigen::ArrayXd& E) override;

    virtual TMD5* GetParamsHash() const override;

  protected:
    void ClearAmpCaches()
    {
      AmpCache<double, T>::clear();
      AmpCache<Eigen::ArrayXd, Eigen::ArrayX<T>>::clear();
    }

    /// Actual implementation of P(). VT is potentially a vector type, if a
    /// vector of energies is passed in. KVT != VT in the case T is a stan
    /// type.
    template<class VT, class KVT> VT _P(int from, int to, const KVT& E);

    template<class VT, class KVT> cmplx<VT> _Amplitude(int from, int to, const KVT& E);

    bool fDirty12, fDirty13, fDirty23, fDirtyCP, fDirtyMasses;

    T s12, c12, s13, c13, s23, c23, sCP, cCP;

    PMNS<T> U;

    inline __attribute__((always_inline)) void UpdatePMNS();

    // This is Hvac without the division by E
    Hermitian<T> H;

    inline __attribute__((always_inline)) void UpdateHamiltonian();

    inline __attribute__((always_inline)) double Hmat();

  private:
    _OscCalc(const _OscCalc&) = default;
    _OscCalc& operator=(const _OscCalc&) = default;
  };
} // end namespaces

// Public names
namespace osc
{
  template<class T> using _OscCalcAnalytic = osc::analytic::_OscCalc<T>;
  using OscCalcAnalytic = _OscCalcAnalytic<double>;
}

#endif
