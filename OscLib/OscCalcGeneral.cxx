/////////////////////////////////////////////////////////////////////////////
// \file    OscCalcGeneral.cxx
// \brief   Source file for general oscillation calculation
// \version $Id: OscCalcGeneral.cxx,v 1.5 2012-09-25 04:51:35 gsdavies Exp $
// \author
/////////////////////////////////////////////////////////////////////////////
#include "OscLib/OscCalcGeneral.h"

#include <boost/serialization/array_wrapper.hpp>

#pragma GCC diagnostic push
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#include "features.h"
#if __GNUC_PREREQ(9,0)
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#endif
#endif
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#pragma GCC diagnostic pop

#include <vector>

#include <cassert>
#include <complex>

#include "OscLib/Constants.h"

namespace osc
{
  const unsigned int kNumFlavours = 3;

  namespace ublas = boost::numeric::ublas;
  typedef std::complex<long double> val_t;
  typedef ublas::bounded_array<val_t, kNumFlavours> alloc_t;
  // Should be a c_matrix, but get ambiguity in triple multiplication
  typedef ublas::bounded_matrix<val_t, kNumFlavours, kNumFlavours> ComplexMat;
  typedef ublas::c_vector<val_t, kNumFlavours> ComplexVec;
  typedef ublas::unit_vector<val_t, alloc_t> UnitVec;
  const ublas::zero_matrix<val_t, alloc_t> kZeroMat(kNumFlavours, kNumFlavours);
  const ublas::identity_matrix<val_t, alloc_t> kIdentity(kNumFlavours);

  // Private data. Doing it this way avoids having to expose complex numbers,
  // matrices etc. to the user.
  struct OscCalcGeneral::Priv
  {
    Priv() :
      emutau(0),
      atmos(kIdentity),
      react(kIdentity),
      solar(kIdentity),
      pmns(kZeroMat),
      dirty(true)
    {
    }

    long double c13, s13;
    long double emutau;
    std::complex<long double> phase;

    ComplexMat atmos, react, solar;
    ComplexMat pmns;

    bool dirty;
  };


  OscCalcGeneral::OscCalcGeneral()
    : d(new OscCalcGeneral::Priv)
  {
  }

  OscCalcGeneral::~OscCalcGeneral()
  {
    delete d;
  }

  IOscCalcAdjustable* OscCalcGeneral::Copy() const
  {
    OscCalcGeneral* ret = new OscCalcGeneral(*this);
    ret->d = new Priv;
    *ret->d = *d;
    return ret;
  }

  void OscCalcGeneral::SetTh23(const double& th23)
  {
    fTh23 = th23;
    d->atmos(2, 2) = d->atmos(1, 1) = cos(th23);
    d->atmos(1, 2) = sin(th23);
    d->atmos(2, 1) = -d->atmos(1, 2);
    d->dirty = true;
  }

  void OscCalcGeneral::SetTh13(const double& th13)
  {
    fTh13 = th13;
    d->c13 = cos(th13);
    d->s13 = sin(th13);

    d->react(2, 2) = d->react(0, 0) = d->c13;
    d->react(0, 2) = d->s13*d->phase;
    d->react(2, 0) = -std::conj(d->react(0, 2));
    d->dirty = true;
  }

  void OscCalcGeneral::SetTh12(const double& th12)
  {
    fTh12 = th12;
    d->solar(1, 1) = d->solar(0, 0) = cos(th12);
    d->solar(0, 1) = sin(th12);
    d->solar(1, 0) = -d->solar(0, 1);
    d->dirty = true;
  }

  void OscCalcGeneral::SetdCP(const double& delta)
  {
    fdCP = delta;

    d->phase = std::polar((long double)1, -(long double)delta);

    d->react(2, 2) = d->react(0, 0) = d->c13;
    d->react(0, 2) = d->s13*d->phase;
    d->react(2, 0) = -std::conj(d->react(0, 2));
    d->dirty = true;
  }

  void OscCalcGeneral::SetNSIEpsilonMuTau(double emutau)
  {
    d->emutau = emutau;
  }

  double OscCalcGeneral::GetNSIEpsilonMuTau() const
  {
    return d->emutau;
  }

  TMD5* OscCalcGeneral::GetParamsHash() const
  {
    // Default isn't good enough if we need to consider NSI
    if(d->emutau) return 0;
    return IOscCalcAdjustable::GetParamsHashDefault("General");
  }


  ComplexMat GetPMNS(OscCalcGeneral::Priv* d)
  {
    if(d->dirty){
      ublas::noalias(d->pmns) = ublas::prod(d->atmos, ublas::prod<ComplexMat>(d->react, d->solar));
      d->dirty = false;
    }
    return d->pmns;
  }


  // U is the PMNS matrix
  // mSq are the squared masses. Pick an arbitrary zero and use the dmsqs
  // E is the neutrino energy
  ComplexMat VacuumHamiltonian(const ComplexMat& U,
                               const std::vector<long double>& mSq,
                               long double E)
  {
    E /= (constants::kkmTom / (constants::kInversemToeV * constants::kGeVToeV));

    assert(mSq.size() == kNumFlavours);

    ComplexMat H = kZeroMat;

    // Loop over rows and columns of the output
    for(unsigned int a = 0; a < kNumFlavours; ++a){
      for(unsigned int b = 0; b < kNumFlavours; ++b){
	// Form sum over mass states
	for(unsigned int i = 0; i < kNumFlavours; ++i){
	  H(a, b) += U(a, i)*mSq[i]/(2*E)*std::conj(U(b, i));
	}
      }
    }

    return H;
  }


  // This is calculated assuming matter neutrinos. So far, this can simply be
  // converted to antineutrinos by taking a minus sign.
  ComplexMat MatterHamiltonianComponent(long double Ne, long double emutau)
  {
    ComplexMat H = kZeroMat;

    const double k = constants::kMatterDensityToEffect / constants::kInversemToeV * constants::kkmTom * constants::kZPerA;

    H(0, 0) = k * Ne;

    // Ignoring conjugates here because we assume e_mutau is real
    H(1, 2) = H(2, 1) = emutau * k * Ne;

    return H;
  }


  /* TODO: maybe look at
     http://en.wikipedia.org/wiki/Baker%E2%80%93Campbell%E2%80%93Hausdorff_formula
     and
     http://en.wikipedia.org/wiki/Pad%C3%A9_approximant
     http://en.wikipedia.org/wiki/Pad%C3%A9_table
  */
  ComplexMat MatrixExp(const ComplexMat& m2)
  {
    // We use the identity e^(a*b) = (e^a)^b
    // First: ensure the exponent is really small, 65536 = 2^16
    ComplexMat m = m2/65536;

    // To first order e^x = 1+x
    ComplexMat ret = kIdentity+m;

    // Raise the result to the power 65536, by squaring it 16 times
    for(int n = 0; n < 16; ++n) ret = ublas::prod(ret, ret);

    return ret;
  }


  ComplexVec EvolveState(ComplexVec A, const ComplexMat& H, long double L)
  {
    const std::complex<long double> i = std::complex<long double>(0, 1);

    return ublas::prod(MatrixExp(-H*L*i), A);
  }

  void conjugate_elements(ComplexMat& m)
  {
    for(unsigned int i = 0; i < kNumFlavours; ++i)
      for(unsigned int j = 0; j < kNumFlavours; ++j)
	m(i, j) = std::conj(m(i, j));
  }

  double OscCalcGeneral::P(int from, int to, double E)
  {
    const int af = abs(from);
    const int at = abs(to);
    assert(af == 12 || af == 14 || af == 16);
    assert(at == 12 || at == 14 || at == 16);

    // No matter<->antimatter transitions
    if(af*at < 0) return 0;

    ComplexMat U = GetPMNS(d);
    // P(a->b|U) = P(abar->bbar|U*)
    if(from < 0) conjugate_elements(U);

    ComplexVec initState = UnitVec(kNumFlavours, 1); // Display accelerator bias
    if(af == 12) initState = UnitVec(kNumFlavours, 0);
    if(af == 16) initState = UnitVec(kNumFlavours, 2);

    std::vector<long double> mSq;
    mSq.push_back(0);
    mSq.push_back(fDmsq21);
    mSq.push_back(fDmsq21 + fDmsq32);

    ComplexMat H = VacuumHamiltonian(U, mSq, E);
    ComplexMat Hm = MatterHamiltonianComponent(fRho, d->emutau);
    // So far, contribution to the antineutrino Hamiltonian is just the negative
    // If there were to be any complex stuff here, would have to think about
    // how it transformed with antineutrinos.
    if(from < 0) H += -1*Hm; else H += Hm;

    ComplexVec finalState = EvolveState(initState, H, fL);

    if(at == 12) return std::norm(finalState(0));
    if(at == 14) return std::norm(finalState(1));
    if(at == 16) return std::norm(finalState(2));

    assert(0 && "Not reached");

    // should never get to this point, but it quiets a compiler error
    // a probability of 0 will be so obviously wrong that it should
    // alert people - ie every probability will be 0
    return 0.;
  }

} // namespace
