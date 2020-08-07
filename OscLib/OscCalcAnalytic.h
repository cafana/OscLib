#ifndef OSCCALCANALYTIC_H
#define OSCCALCANALYTIC_H

#include "OscLib/IOscCalc.h"

#include <unordered_map>

namespace osc
{
  /// std::complex takes a lot of care with inf/nan which we don't want
  template<class T> struct cmplx
  {
    cmplx(T r, T i) : re(r), im(i) {}

    inline __attribute__((always_inline)) T norm() const {return re*re + im*im;}
    inline __attribute__((always_inline)) cmplx conj() const {return cmplx(re, -im);}

    T re, im;
  };

  template<class T> class _OscCalcAnalytic: public _IOscCalcAdjustable<T>
  {
  public:
    _OscCalcAnalytic();
    virtual ~_OscCalcAnalytic();
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

    virtual TMD5* GetParamsHash() const override;

  protected:
    bool fDirty12, fDirty13, fDirty23, fDirtyCP, fDirtyMasses;

    T s12, c12, s13, c13, s23, c23, sCP, cCP;

    /*T        Ue1;*/ T        Ue2; cmplx<T>  Ue3;
    /*cmplx<T> Um1;*/ cmplx<T> Um2; T         Um3;
    /*cmplx<T> Ut1;*/ cmplx<T> Ut2; T         Ut3;

    inline __attribute__((always_inline)) void UpdatePMNS();

    // This is Hvac without the division by E
    T Hee;            cmplx<T>  Hem; cmplx<T>  Het;
    /*cmplx<T> Hme;*/ T         Hmm; cmplx<T>  Hmt;
    /*cmplx<T> Hte; cmplx<T>  Htm;*/ T         Htt;

    inline __attribute__((always_inline)) void UpdateHamiltonian();

    inline __attribute__((always_inline)) double Hmat();

    T Mee;     cmplx<T>  Mem;   cmplx<T>  Met;
    /*cmplx  Mme;*/ T Mmm;   cmplx<T>  Mmt;
    /*cmplx  Mte;   cmplx  Mtm;*/ T Mtt;

    struct Eigenvalues
    {
      //      std::array<double, 3> xs; ///< eigenvalues
      //      std::array<cmplx,  3> expixs; ///< exp(i*x)/(3*x^2+b*x+c) for each x
      // Turns out we only need these summary expressions of the eigenvalues
      cmplx<T> sume, sumxe, sumxxe;
    };

    inline __attribute__((always_inline))
    Eigenvalues GetEigenvalues(double E);

    class Probs
    {
    public:
      Probs(T ee, T me, T em, T mm)
        : Pee(ee), Pme(me), Pem(em), Pmm(mm)
      {
      }

      inline __attribute__((always_inline)) T P(int from, int to) const;

    protected:
      T Pee, Pme, Pem, Pmm;
    };

    std::unordered_map<double, Probs> fProbCache;

  private:
    _OscCalcAnalytic(const _OscCalcAnalytic&) = default;
    _OscCalcAnalytic& operator=(const _OscCalcAnalytic&) = default;
  };

  typedef _OscCalcAnalytic<double> OscCalcAnalytic;
} // end namespace

#endif
