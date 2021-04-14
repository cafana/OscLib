#ifndef OSC_OSCCALCULATORSTERILE_H
#define OSC_OSCCALCULATORSTERILE_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// \file OscCalcSterile.h                                               //
//                                                                      //
// Adapt the PMNS_Sterile calculator to standard interface              //
// <aurisaam@ucmail.uc.edu>                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "OscLib/IOscCalcSterile.h"
#include "OscLib/PMNS_Sterile.h"

#include <vector>

namespace osc
{
  /// \brief Adapt the PMNS_Sterile calculator to standard interface
  ///
  /// Adapt the \ref PMNS_Sterile calculator (3+N with matter effects) to standard interface
  class OscCalcSterile: public IOscCalcSterile
  {
  public:
    using IOscCalcAdjustable::P;
    OscCalcSterile();
    OscCalcSterile(const OscCalcSterile& calc);
    virtual ~OscCalcSterile();

    void SetNFlavors(int nflavors);

    virtual IOscCalcAdjustable* Copy() const override;

    virtual void SetAngle(int i, int j, double th);
    virtual void SetDelta(int i, int j, double delta);
    virtual void SetDm(int i, double dm);

    virtual double GetDm(int i) const override
      { return fPMNS_Sterile->GetDm(i); }
    virtual double GetAngle(int i, int j) const override
      { return fPMNS_Sterile->GetAngle(i, j); }
    virtual double GetDelta(int i, int j) const override
      { return fPMNS_Sterile->GetDelta(i, j); }
    
    // if flavAfter == 0, give the active fraction
    virtual double P(int flavBefore, int flavAfter, double E) override;

    void SetState(std::vector<double> state);

    //Getters
    int GetNFlavors() const { return fPMNS_Sterile->GetNFlavors(); }
    std::vector<double> GetState() const;
    virtual TMD5* GetParamsHash() const override;

  protected:
    PMNS_Sterile* fPMNS_Sterile;

    int    fNFlavors;
    double fPrevE;
    int    fPrevAnti;
    int    fPrevFlavBefore;
  };

} // namespace

#endif
