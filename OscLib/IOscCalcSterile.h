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
  class IOscCalcSterile: public IOscCalcAdjustable
  {
  public:
    virtual ~IOscCalcSterile() {};

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

  protected:

    virtual void SetDmsq21(const double& dmsq21) override;
    virtual void SetDmsq32(const double& dmsq32) override;
    virtual void SetTh12  (const double& th12  ) override;
    virtual void SetTh13  (const double& th13  ) override;
    virtual void SetTh23  (const double& th23  ) override;
    virtual void SetdCP   (const double& dCP   ) override;

    double fRho; 
    bool   fDirty;

  };

  class OscCalcSterileTrivial: public IOscCalcSterile
  {
  public:
    using IOscCalcAdjustable::P;
    OscCalcSterileTrivial() {};
    virtual ~OscCalcSterileTrivial() {};
    virtual double P(int, int, double) override;

  private:
    virtual IOscCalcAdjustable* Copy() const override;

    virtual void SetAngle(int, int, double) override {};
    virtual void SetDelta(int, int, double) override {};
    virtual void SetDm(int, double)         override {};

    virtual double GetDm(int)         const override { return 0; };
    virtual double GetAngle(int, int) const override { return 0; };
    virtual double GetDelta(int, int) const override { return 0; };
  };

  /// \brief version of OscCalcSterile that always returns probability of 1

  const IOscCalcSterile* DowncastToSterile(const IOscCalc* calc);
  IOscCalcSterile* DowncastToSterile(IOscCalc* calc);

} // namespace

#endif
