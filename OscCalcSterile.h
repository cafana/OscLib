#ifndef OSC_OSCCALCULATORSTERILE_H
#define OSC_OSCCALCULATORSTERILE_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// \file OscCalcSterile.h                                         //
//                                                                      //
// Adapt the PMNS_Sterile calculator to standard interface              //
// <aurisaam@ucmail.uc.edu>						//
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "IOscCalc.h"
#include "PMNS_Sterile.h"
#include <vector>

namespace osc
{
  /// \brief Adapt the PMNS_Sterile calculator to standard interface
  ///
  /// Adapt the \ref PMNS_Sterile calculator (3+N with matter effects) to standard interface
  class OscCalcSterile: public IOscCalcAdjustable
  {
  public:
    using IOscCalcAdjustable::P;
    OscCalcSterile();
    OscCalcSterile(const OscCalcSterile& calc);
    virtual ~OscCalcSterile();

    void SetNFlavors(int nflavors);

    virtual IOscCalcAdjustable* Copy() const override;
    
    // if flavAfter == 0, give the active fraction
    virtual double P(int flavBefore, int flavAfter, double E) override;

    virtual void SetL  (double L  ) override {fDirty = true; fL   = L;}
    virtual void SetRho(double rho) override {fDirty = true; fRho = rho;}

    void SetAngle(int i, int j, double th);
    void SetDelta(int i, int j, double delta);
    void SetDm(int i, double dm);

    void SetState(std::vector<double> state);

    //Getters
    int    GetNFlavors()          const { return fPMNS_Sterile->GetNFlavors(); }
    double GetL()                 const override { return fL; }
    double GetRho()               const override { return fRho; }
    double GetDm(int i)           const { return fPMNS_Sterile->GetDm(i); }
    double GetAngle(int i, int j) const { return fPMNS_Sterile->GetAngle(i, j); }
    double GetDelta(int i, int j) const { return fPMNS_Sterile->GetDelta(i, j); }
    std::vector<double> GetState() const;
    virtual TMD5* GetParamsHash() const override;

  protected:
    PMNS_Sterile* fPMNS_Sterile;

    virtual void SetDmsq21(const double& dmsq21) override;
    virtual void SetDmsq32(const double& dmsq32) override;
    virtual void SetTh12  (const double& th12  ) override;
    virtual void SetTh13  (const double& th13  ) override;
    virtual void SetTh23  (const double& th23  ) override;
    virtual void SetdCP   (const double& dCP   ) override;

    int    fNFlavors;
    double fRho; 
    bool   fDirty;
    double fPrevE;
    int    fPrevAnti;
    int    fPrevFlavBefore;
  };

  /// \brief Version of OscCalcSterile that always returns probability of 1
  class OscCalcSterileTrivial: public OscCalcSterile
  {
  public:
    using IOscCalc::P;
    OscCalcSterileTrivial();
    OscCalcSterileTrivial(const OscCalcSterile& calc);
    OscCalcSterileTrivial(const OscCalcSterileTrivial& calc);
    virtual ~OscCalcSterileTrivial() {};

    virtual IOscCalcAdjustable* Copy() const override;
    
    // Always return 1
    virtual double P(int flavBefore, int flavAfter, double E) override;
  };

  const OscCalcSterile* DowncastToSterile(const IOscCalc* calc);
  OscCalcSterile* DowncastToSterile(IOscCalc* calc);

} // namespace

#endif
