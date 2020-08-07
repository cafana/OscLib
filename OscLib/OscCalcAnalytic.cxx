#include "OscLib/OscCalcAnalytic.h"

#include <cassert>

// Stan doesn't provide sincos(). Make sure not to hide the fast one for double
template<class T> std::enable_if_t<!std::is_arithmetic_v<T>, void> sincos(const T& x, T* sx, T* cx)
{
  *sx = sin(x);
  *cx = cos(x);
}

namespace osc
{
  template<class T> inline __attribute__((always_inline)) cmplx<T>
  operator*(const cmplx<T>& x, const cmplx<T>& y)
  {
    return cmplx<T>(x.re*y.re - x.im*y.im, x.re*y.im + x.im*y.re);
  }

  template<class T, class U> inline __attribute__((always_inline)) cmplx<T>
  operator*(const U& x, const cmplx<T>& y)
  {
    return cmplx<T>(x*y.re, x*y.im);
  }

  template<class T, class U> inline __attribute__((always_inline)) cmplx<T>
  operator*(const cmplx<T>& x, const U& y)
  {
    return cmplx<T>(x.re*y, x.im*y);
  }

  template<class T, class U> inline __attribute__((always_inline)) cmplx<T>
  operator/(const cmplx<T>& x, const U& y)
  {
    return cmplx<T>(x.re/y, x.im/y);
  }

  template<class T> inline __attribute__((always_inline)) cmplx<T>
  operator+(const cmplx<T>& x, const cmplx<T>& y)
  {
    return cmplx<T>(x.re + y.re, x.im + y.im);
  }

  template<class T, class U> inline __attribute__((always_inline)) cmplx<T>
  operator+(const cmplx<T>& x, const U& y)
  {
    return cmplx<T>(x.re + y, x.im);
  }

  template<class T, class U> inline __attribute__((always_inline)) cmplx<T>
  operator+(const U& x, const cmplx<T>& y)
  {
    return cmplx<T>(x + y.re, y.im);
  }

  template<class T> inline __attribute__((always_inline)) cmplx<T>
  operator-(const cmplx<T>& x, const cmplx<T>& y)
  {
    return cmplx<T>(x.re - y.re, x.im - y.im);
  }

  template<class T, class U> inline __attribute__((always_inline)) cmplx<T>
  operator-(const cmplx<T>& x, const U& y)
  {
    return cmplx<T>(x.re - y, x.im);
  }

  template<class T, class U> inline __attribute__((always_inline)) cmplx<T>
  operator-(const U& x, const cmplx<T>& y)
  {
    return cmplx<T>(x - y.re, -y.im);
  }

  template<class T> inline __attribute__((always_inline)) cmplx<T>
  operator-(const cmplx<T>& x)
  {
    return cmplx<T>(-x.re, -x.im);
  }

  //---------------------------------------------------------------------------
  template<class T> _OscCalcAnalytic<T>::_OscCalcAnalytic()
    : fDirty12(true), fDirty13(true), fDirty23(true), fDirtyCP(true), fDirtyMasses(true),
      Ue3(0, 0), Um2(0, 0), Ut2(0, 0),
      Hem(0, 0), Het(0, 0), Hmt(0, 0),
      Mem(0, 0), Met(0, 0), Mmt(0, 0)
  {
  }

  //---------------------------------------------------------------------------
  template<class T> _OscCalcAnalytic<T>::~_OscCalcAnalytic()
  {
  }

  //---------------------------------------------------------------------------
  template<class T> _IOscCalcAdjustable<T>* _OscCalcAnalytic<T>::
  Copy() const
  {
    return new _OscCalcAnalytic<T>(*this);
  }

  //---------------------------------------------------------------------------
  template<class T> void _OscCalcAnalytic<T>::SetL(double L)
  {
    if(L == this->fL) return;

    this->fL = L;
    fProbCache.clear();
  }

  //---------------------------------------------------------------------------
  template<class T> void _OscCalcAnalytic<T>::SetRho(double rho)
  {
    if(rho == this->fRho) return;

    this->fRho = rho;
    fProbCache.clear();
  }

  //---------------------------------------------------------------------------
  template<class T> void _OscCalcAnalytic<T>::SetDmsq21(const T& dmsq21)
  {
    if constexpr(std::is_arithmetic_v<T>) if(dmsq21 == this->fDmsq21) return;

    this->fDmsq21 = dmsq21;
    fDirtyMasses = true;
  }

  //---------------------------------------------------------------------------
  template<class T> void _OscCalcAnalytic<T>::SetDmsq32(const T& dmsq32)
  {
    if constexpr(std::is_arithmetic_v<T>) if(dmsq32 == this->fDmsq32) return;

    this->fDmsq32 = dmsq32;
    fDirtyMasses = true;
  }

  //---------------------------------------------------------------------------
  template<class T> void _OscCalcAnalytic<T>::SetTh23(const T& th23)
  {
    if constexpr(std::is_arithmetic_v<T>) if(th23 == this->fTh23) return;

    this->fTh23 = th23;
    fDirty23 = true;
  }

  //---------------------------------------------------------------------------
  template<class T> void _OscCalcAnalytic<T>::SetTh13(const T& th13)
  {
    if constexpr(std::is_arithmetic_v<T>) if(th13 == this->fTh13) return;

    this->fTh13 = th13;
    fDirty13 = true;
  }

  //---------------------------------------------------------------------------
  template<class T> void _OscCalcAnalytic<T>::SetTh12(const T& th12)
  {
    if constexpr(std::is_arithmetic_v<T>) if(th12 == this->fTh12) return;

    this->fTh12 = th12;
    fDirty12 = true;
  }

  //---------------------------------------------------------------------------
  template<class T> void _OscCalcAnalytic<T>::SetdCP(const T& delta)
  {
    if constexpr(std::is_arithmetic_v<T>) if(delta == this->fdCP) return;

    this->fdCP = delta;
    fDirtyCP = true;
  }

  //---------------------------------------------------------------------------
  template<class T> TMD5* _OscCalcAnalytic<T>::GetParamsHash() const
  {
    return _IOscCalcAdjustable<T>::GetParamsHashDefault("Analytic");
  }

  //---------------------------------------------------------------------------
  template<class T> double _OscCalcAnalytic<T>::Hmat()
  {
    // Need to convert avogadro's constant so that the total term comes out in
    // units of inverse distance. Note that Ne will be specified in g/cm^-3
    // I put the following into Wolfram Alpha:
    // (fermi coupling constant)*((avogadro number)/cm^3)*(reduced planck constant)^2*(speed of light)^2
    // And then multiplied by 1000 because we specify L in km not m.
    const double GF = 1.368e-4;
    const double Ne = this->fRho;
    return sqrt(2)*GF*Ne;
  }

  //---------------------------------------------------------------------------
  template<class T> T sqr(T x){return x*x;}
  template<class T> T cube(T x){return x*x*x;}

  const double sqrt3 = sqrt(3);

  //---------------------------------------------------------------------------
  /// Solve x^3 + b*x^2 + c*x + d = 0
  template<class T> std::array<T, 3> SolveCubic(T b, T c, T d)
  {
    b /= 3;
    // Now solving x^3 + 3*b*x^2 + c*x + d = 0

    // https://en.wikipedia.org/wiki/Cubic_equation#Depressed_cubic
    const T p = c/3 - sqr(b);
    const T q = 2*cube(b) - b*c + d;
    // Now solving t^3 + 3*p*t + q = 0

    // https://en.wikipedia.org/wiki/Cubic_equation#Trigonometric_solution_for_three_real_roots
    //
    // Cardano's formula looks simpler, but you still have to involve trig
    // operations anyway for the cube root, and it's harder to keep everything
    // explicitly real.

    const T s = sqrt(-p);

    const T r = acos(q/(2*p*s)) / 3;

    // t_k = 2*s * cos(r/3 + 2pi*k/3)

    // Use cos(a+b) = cosa*cosb - sina*sinb to save one trig operation
    T sinr, cosr;
    sincos(r, &sinr, &cosr);

    const T t0 = 2*s*cosr;
    const T t1 = s*(sqrt3*sinr - cosr);
    const T t2 = -t0-t1;

    return {t0 - b, t1 - b, t2 - b};
  }

  //---------------------------------------------------------------------------
  template<class T> typename _OscCalcAnalytic<T>::Eigenvalues _OscCalcAnalytic<T>::GetEigenvalues(double E)
  {
    // const T a = 1; // cube term

    const T b = -Mee-Mmm-Mtt; // square term

    const T Lem = Mem.norm();
    const T Let = Met.norm();
    const T Lmt = Mmt.norm();

    // Matrix is Hermitian
    const T c =  Mee*Mmm + Mee*Mtt + Mmm*Mtt - Lem - Let - Lmt; // linear term

    const T d = (Mee*Lmt + Mmm*Let + Mtt*Lem
                      - Mee*Mmm*Mtt
                      -2*(Mem.re * Mmt.re * Met.re + Mem.re * Mmt.im * Met.im + Mem.im * Mmt.re * Met.im - Mem.im * Mmt.im * Met.re)); // const term

    const std::array<T, 3> xs = SolveCubic(b, c, d);

    // Overall phase doesn't matter, which allows us to save one exponentiation
    T c10, s10, c20, s20;
    sincos(xs[1]-xs[0], &s10, &c10);
    sincos(xs[2]-xs[0], &s20, &c20);

    const T ei0 = 1/(3*sqr(xs[0]) + 2*b*xs[0] + c);
    const cmplx<T> ei1 = cmplx(c10, s10)/(3*sqr(xs[1]) + 2*b*xs[1] + c);
    const cmplx<T> ei2 = cmplx(c20, s20)/(3*sqr(xs[2]) + 2*b*xs[2] + c);

    // TODO, is there a way to more cheaply calculate these combinations?
    return {
      ei0 + ei1 + ei2,
      xs[0]*ei0 + xs[1]*ei1 + xs[2]*ei2,
      sqr(xs[0])*ei0 + sqr(xs[1])*ei1 + sqr(xs[2])*ei2};
  }

  //---------------------------------------------------------------------------
  template<class T> void _OscCalcAnalytic<T>::UpdatePMNS()
  {
    Ue2 = s12*c13;
    Um3 = s23*c13;
    Ut3 = c23*c13;

    const cmplx<T> phase(cCP, sCP);
    Ue3 = s13*phase.conj();
    Um2 =  c12*c23-(s12*s23*s13)*phase;
    Ut2 = -c12*s23-(s12*c23*s13)*phase;
  }

  //---------------------------------------------------------------------------
  template<class T> void _OscCalcAnalytic<T>::UpdateHamiltonian()
  {
    // d1 would be zero, so all those terms drop out
    const T d2 = this->fDmsq21;
    const T d3 = this->fDmsq21 + this->fDmsq32;

    Hee = d2 * sqr(Ue2)   + d3 * Ue3.norm();
    Hmm = d2 * Um2.norm() + d3 * sqr(Um3);
    Htt = d2 * Ut2.norm() + d3 * sqr(Ut3);

    Hem = d2 * Ue2 * Um2.conj() + d3 * Ue3 * Um3;
    Het = d2 * Ue2 * Ut2.conj() + d3 * Ue3 * Ut3;
    Hmt = d2 * Um2 * Ut2.conj() + d3 * Um3 * Ut3;

    fProbCache.clear();
  }

  //---------------------------------------------------------------------------
  template<class T> T _OscCalcAnalytic<T>::Probs::P(int from, int to) const
  {
    // convert flavours to indices into matrix
    const int i0 = (from-12)/2;
    const int i1 = (to-12)/2;

    // Exploit unitarity
    switch(i0*3+i1){
    case 0: return Pee;
    case 1: return Pme;
    case 2: return 1-Pee-Pme; // Pte
    case 3: return Pem;
    case 4: return Pmm;
    case 5: return 1-Pem-Pmm; // Ptm
    case 6: return 1-Pee-Pem; // Pet
    case 7: return 1-Pme-Pmm; // Pmt
    case 8: return Pee+Pem+Pme+Pmm-1; // Ptt
    default: abort();
    }
  }

  //---------------------------------------------------------------------------
  template<class T> T _OscCalcAnalytic<T>::P(int from, int to, double E)
  {
    // -E effectively flips rho and conjugates H
    if(from < 0) return P(-from, -to, -E);

    assert(from > 0 && to > 0);

    assert(from == 12 || from == 14 || from == 16);
    assert(to == 12 || to == 14 || to == 16);

    const bool dirtyAngles = fDirty12 || fDirty13 || fDirty23 || fDirtyCP;

    if(dirtyAngles){
      if(fDirty12) sincos(this->fTh12, &s12, &c12);
      if(fDirty13) sincos(this->fTh13, &s13, &c13);
      if(fDirty23) sincos(this->fTh23, &s23, &c23);
      if(fDirtyCP) sincos(this->fdCP,  &sCP, &cCP);
      UpdatePMNS();
      UpdateHamiltonian();
    }
    else{
      if(fDirtyMasses){
        UpdateHamiltonian();
      }
      else{
        auto it = fProbCache.find(E);
        if(it != fProbCache.end()) return it->second.P(from, to);
      }
    }

    fDirty12 = fDirty13 = fDirty23 = fDirtyCP = fDirtyMasses = false;

    const double k = -this->fL * 2*1.267 / E;
    Mee = Hee * k  - this->fL * Hmat();
    Mem = Hem * k;
    Mmm = Hmm * k;
    Met = Het * k;
    Mmt = Hmt * k;
    Mtt = Htt * k;

    // Matrix exponent is based on https://www.wolframalpha.com/input/?i=matrixExp+%5B%5Br%2Cs%2Ct%5D%2C%5Bu%2Cv%2Cw%5D%2C%5Bx%2Cy%2Cz%5D%5D

    const Eigenvalues es = GetEigenvalues(E);

    const T Aee = Mmm*Mtt - Mmt.norm();
    const T Amm = Mee*Mtt - Met.norm();
    const cmplx<T> Aem = Met*Mmt.conj() - Mem*Mtt;

    const Probs ps((Aee       *es.sume - (Mmm+Mtt)  *es.sumxe + es.sumxxe).norm(),
                   (Aem.conj()*es.sume +  Mem.conj()*es.sumxe            ).norm(),
                   (Aem       *es.sume +  Mem       *es.sumxe            ).norm(),
                   (Amm       *es.sume - (Mee+Mtt)  *es.sumxe + es.sumxxe).norm());

    fProbCache.emplace(E, ps);

    return ps.P(from, to);
  }

} // namespace


// Instantiate
template class osc::_OscCalcAnalytic<double>;

#ifdef OSCLIB_STAN
#include "stan/math/rev/scal.hpp"
template class osc::_OscCalcAnalytic<stan::math::var>;
#endif
