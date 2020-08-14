//////////////////////////////////////////////////////////
/// \brief Simple oscillation probability calculator that
///        has no solar term or mass hierarchy  or delta
///        so it's some kind of average of all of those
//////////////////////////////////////////////////////////

#include "OscLib/IOscCalc.h"

#include <cmath>

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
    virtual double P(int flavBefore, int flavAfter, double E) override
    {
      const double L      = 810.;    // km
      const double ldm    = 2.40e-3; // large delta m^2, m23
      const double ss2t13 = 0.1;     // sin^2(2theta_13)
      const double ss2t23 = 1.;      // maximal sin^2(2theta_23)
      const double ssth23 = 0.5;     // sin^2(theta_23) corresponding to above

      // Signal
      if(abs(flavAfter) == 12 && abs(flavBefore) == 14){
	// Nue appearance
	return ss2t13*ssth23*sqr(sin(1.267*ldm*L/E));
      }
      if(abs(flavAfter) == 14 && abs(flavBefore) == 14){
	// CC mu disappearance
	return 1-ss2t23*sqr(sin(1.267*ldm*L/E));
      }

      // Background
      if(abs(flavAfter) == 12 && abs(flavBefore) == 12){
	// Beam nue
	return 1-ss2t13*sqr(sin(1.267*ldm*L/E));
      }
      if(abs(flavAfter) == 14 && abs(flavBefore) == 12){
	// CC mu appearance
	return ss2t13*ssth23*sqr(sin(1.267*ldm*L/E));
      }
      if(abs(flavAfter) == 16 && abs(flavBefore) == 14){
        //numu to nutau CC appearance
        return (1-ss2t13)*ss2t23*sqr(sin(1.267*ldm*L/E));
      }
      if(abs(flavAfter) == 16 && abs(flavBefore) == 12){
        //nue to nutau CC appearance
        return ss2t13*(1-ssth23)*sqr(sin(1.267*ldm*L/E));
      }
      if(abs(flavAfter) == 16 && abs(flavBefore) == 16){
        //nutau to nutau CC disappearance
        return 1-(ss2t23-ss2t23*ss2t13+ss2t13-ss2t13*ssth23)*sqr(sin(1.267*ldm*L/E));
      }

      // Don't know what this is
      return 0;
    }

  protected:
    template<class T> T sqr(T x) const{return x*x;}
  };
}
