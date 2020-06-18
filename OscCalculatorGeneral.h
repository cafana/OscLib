#ifndef OSCCALCULATORGENERAL_H
#define OSCCALCULATORGENERAL_H

#include "OscLib/func/IOscCalculator.h"

namespace osc
{

  /// \brief More generic (but probably slower) oscillation calculations
  ///
  /// Calculates the oscillation probabilities from first principles
  /// (constructs the Hamiltonian and exponentiates it). The usual oscillation
  /// calculator (a power expansion of the exponentiation in alpha and
  /// sin(th13)) is probably faster, but this implementation has the advantages
  /// of simplicity and flexibility. A simple non-standard-interactions model
  /// is included.  Addition of a general four-flavour sterile neutrino model,
  /// or inclusion of varying matter densities would be relatively
  /// trivial. Also, it's always nice to have cross-checks.
  class OscCalculatorGeneral: public IOscCalculatorAdjustable
  {
  public:
    OscCalculatorGeneral();
    virtual ~OscCalculatorGeneral();
    using IOscCalculator::P;

    virtual IOscCalculatorAdjustable* Copy() const override;

    // Baseline in km
    virtual void SetL(double L) override {fL = L;}
    // Density in g/cm^3
    virtual void SetRho(double rho) override {fRho = rho;}
    // in eV^2
    virtual void SetDmsq21(const double& dmsq21) override {fDmsq21 = dmsq21;}
    // This is a signed quantity, use a negative value for inverted hierarchy
    virtual void SetDmsq32(const double& dmsq32) override {fDmsq32 = dmsq32;}
    // In radians
    virtual void SetTh12(const double& th12) override;
    virtual void SetTh13(const double& th13) override;
    virtual void SetTh23(const double& th23) override;
    virtual void SetdCP(const double& dCP) override;

    void SetNSIEpsilonMuTau(double emutau);
    double GetNSIEpsilonMuTau() const;

    virtual double P(int from, int to, double E) override;

    virtual TMD5* GetParamsHash() const override;

    struct Priv;
  protected:
    Priv* d;

  private:
    OscCalculatorGeneral(const OscCalculatorGeneral&) = default;
    OscCalculatorGeneral& operator=(const OscCalculatorGeneral&) = default;
  };

} // end namespace

#endif
