#ifndef OSC_OSCCALCULATORSTERILEBEAM_H
#define OSC_OSCCALCULATORSTERILEBEAM_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// \file OscCalcSterileBeam.h                                     //
//                                                                      //
// Adapt the PMNS_Sterile calculator to standard interface              //
// <aurisaam@ucmail.uc.edu>                                             //
// <kasettisivaprasad@gmail.com>				     	//
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "IOscCalc.h"
#include "OscCalcSterile.h"
//#include "INoOscillations.h"
#include "PMNS_Sterile.h"

namespace osc
{
  /// \brief Adapt the PMNS_Sterile calculator to standard interface
  ///
  /// Adapt the \ref PMNS_Sterile calculator (3+N with matter effects) to standard interface
  class OscCalcSterileBeam: public OscCalcSterile //, public IOscCalcAdjustable
  {
  public:
    using IOscCalc::P;
    OscCalcSterileBeam();
    virtual ~OscCalcSterileBeam();

    OscCalcSterileBeam(const OscCalcSterileBeam& calc);

    virtual IOscCalcAdjustable* Copy() const override;

    std::string kBeamMode;

    virtual void SetKaonScale(double scale);
    virtual void SetPionScale(double scale);
    virtual void SetMuonScale(double scale);

    double GetKaonScale() const;
    double GetPionScale() const;
    double GetMuonScale() const;

    virtual TMD5* GetParamsHash() const override;

  protected:

    double fKaonscale;
    double fPionscale;
    double fMuonscale;

/*  PMNS_Sterile* fPMNS_Sterile;
    int    fNFlavors;
    double fRho;
    bool   fDirty;
    double fPrevE;
    int    fPrevAnti;
    int    fPrevFlavBefore;
*/

  };

  const OscCalcSterileBeam* DowncastToSterileBeam(const IOscCalc* calc);
  OscCalcSterileBeam* DowncastToSterileBeam(IOscCalc* calc);

} // namespace

#endif
