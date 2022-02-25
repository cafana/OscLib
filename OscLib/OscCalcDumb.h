//////////////////////////////////////////////////////////
/// \brief Simple oscillation probability calculator that
///        has no solar term or mass hierarchy  or delta
///        so it's some kind of average of all of those
//////////////////////////////////////////////////////////

#include "OscLib/IOscCalc.h"

namespace osc
{
  /// \brief Simple oscillation probability calculator that has no solar term
  ///        or mass hierarchy or delta so it's some kind of average of all of
  ///        those.
  class OscCalcDumb: public osc::IOscCalc
  {
  public:
    virtual IOscCalc* Copy() const override {return new OscCalcDumb;}

    using IOscCalc::P;
    virtual double P(int flavBefore, int flavAfter, double E) override;

    virtual void Print(const std::string& prefix = "") const override;
  };
}
