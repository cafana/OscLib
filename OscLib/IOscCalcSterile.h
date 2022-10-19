#ifndef OSC_IOSCCALCSTERILE_H
#define OSC_IOSCCALCSTERILE_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// \file IOscCalcSterile.h                                              //
//                                                                      //
// Base class for sterile oscillation calculator                        //
// <jhewes15@fnal.gov>                                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "OscLib/IOscCalc.h"

namespace osc
{
  /// \brief base class for sterile oscillation calculators
  /// In the context of a sterile oscillation calculator, a PDG code
  /// of zero corresponds to the survival probability for an active
  /// neutrino, ie. the sum of the oscillation probabilities for the
  /// three active neutrino states.
  class IOscCalcSterile: public IOscCalcAdjustable
  {
  public:
    virtual ~IOscCalcSterile() {}

    virtual void SetL  (double L  ) override {fDirty = true; fL   = L;}
    virtual void SetRho(double rho) override {fDirty = true; fRho = rho;}

    virtual void SetAngle(int i, int j, double th)    = 0;
    virtual void SetDelta(int i, int j, double delta) = 0;
    virtual void SetDm(int i, double dm)              = 0;

    double GetL()                 const override { return fL; }
    double GetRho()               const override { return fRho; }
    virtual double GetDm(int i)           const = 0;
    virtual double GetAngle(int i, int j) const = 0;
    virtual double GetDelta(int i, int j) const = 0;

    // Provide convenience methods for parameters that overlap with the
    // traditional 3 flavour parameters.
    virtual double GetDmsq21() const override { return GetDm(2)-GetDm(1); }
    virtual double GetDmsq32() const override { return GetDm(3)-GetDm(2); }
    virtual double GetTh12  () const override { return GetAngle(1, 2); }
    virtual double GetTh13  () const override { return GetAngle(1, 3); }
    virtual double GetTh23  () const override { return GetAngle(2, 3); }
    virtual double GetdCP   () const override { return GetDelta(1, 3); }

    // Delete the default Print() implementation from IOscCalcAdjustable,
    // because it won't print any of the extra sterile parameters.
    virtual void Print(const std::string& prefix = "") const override = 0;

  protected:
    // Default implementation of Print() for derived classes to use. Can't
    // override underlying Print() because this interface doesn't know the
    // number of generations.
    void PrintImpl(int nNus, const std::string& prefix = "") const;

    virtual void SetDmsq21(const double& dmsq21) override;
    virtual void SetDmsq32(const double& dmsq32) override;
    virtual void SetTh12  (const double& th12  ) override;
    virtual void SetTh13  (const double& th13  ) override;
    virtual void SetTh23  (const double& th23  ) override;
    virtual void SetdCP   (const double& dCP   ) override;

    double fRho; 
    bool   fDirty;

  };

  /// \brief version of OscCalcSterile that always returns probability of 1
  class OscCalcSterileTrivial: public IOscCalcSterile
  {
  public:
    using IOscCalcAdjustable::P;
    OscCalcSterileTrivial() {};
    virtual ~OscCalcSterileTrivial() {};
    virtual double P(int, int, double) override;

    virtual void Print(const std::string& prefix = "") const override;

  private:
    virtual IOscCalcAdjustable* Copy() const override;

    virtual void SetAngle(int, int, double) override {};
    virtual void SetDelta(int, int, double) override {};
    virtual void SetDm(int, double)         override {};

    virtual double GetDm(int)         const override { return 0; };
    virtual double GetAngle(int, int) const override { return 0; };
    virtual double GetDelta(int, int) const override { return 0; };
  };

  const IOscCalcSterile* DowncastToSterile(const IOscCalc* calc, bool quiet=false);
  IOscCalcSterile* DowncastToSterile(IOscCalc* calc, bool quiet=false);

} // namespace

#endif
