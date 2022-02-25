#include "OscLib/IOscCalcSterile.h"

#include <iostream>

namespace osc
{
  //---------------------------------------------------------------------------
  void IOscCalcSterile::SetDmsq21(const double&)
  {
    std::cerr << "Must use SetDm!" << std::endl;
    assert(false);
  }

  //---------------------------------------------------------------------------
  void IOscCalcSterile::SetDmsq32(const double&)
  {
    std::cerr << "Must use SetDm!" << std::endl;
    assert(false);
  }

  //---------------------------------------------------------------------------
  void IOscCalcSterile::SetTh12(const double&)
  {
    std::cerr << "Must use SetAngle!" << std::endl;
    assert(false);
  }

  //---------------------------------------------------------------------------
  void IOscCalcSterile::SetTh13(const double&)
  {
    std::cerr << "Must use SetAngle!" << std::endl;
    assert(false);
  }

  //---------------------------------------------------------------------------
  void IOscCalcSterile::SetTh23(const double&)
  {
    std::cerr << "Must use SetAngle!" << std::endl;
    assert(false);
  }

  //---------------------------------------------------------------------------
  void IOscCalcSterile::SetdCP(const double&)
  {
    std::cerr << "Must use SetDelta!" << std::endl;
    assert(false);
  }

  //---------------------------------------------------------------------------
  void IOscCalcSterile::PrintImpl(int nNus, const std::string& prefix) const
  {
    for(int i = 2; i <= nNus; ++i){
      std::cout << prefix << "dmsq" << i << "1 = " << GetDm(i) << " eV^2\n";
    }

    for(int i = 1; i < nNus; ++i){
      for(int j = i+1; j <= nNus; ++j){
        std::cout << prefix << "theta" << i << j << " = " << GetAngle(i, j) << "\n";
      }
    }

    for(int i = 1; i < nNus; ++i){
      for(int j = i+1; j <= nNus; ++j){
        std::cout << prefix << "delta" << i << j << " = " << GetDelta(i, j) << "\n";
      }
    }

    std::cout << prefix << "L = " << GetL() << " km\n"
              << prefix << "rho = " << GetRho() << " g/cm^3" << std::endl;
  }

  //---------------------------------------------------------------------------
  double OscCalcSterileTrivial::P(int, int, double)
  {
    return 1;
  }

  //---------------------------------------------------------------------------
  IOscCalcAdjustable* OscCalcSterileTrivial::Copy() const
  {
    return new OscCalcSterileTrivial();
  }

  //---------------------------------------------------------------------------
  void OscCalcSterileTrivial::Print(const std::string& prefix) const
  {
    std::cout << prefix << "Trivial sterile osc" << std::endl;
  }

  //---------------------------------------------------------------------------
  const IOscCalcSterile* DowncastToSterile(const IOscCalc* calc, bool quiet)
  {
    const IOscCalcSterile* calc_sterile
            = dynamic_cast<const IOscCalcSterile*>(calc);
    if(calc_sterile) return calc_sterile;
    if (!quiet)
      std::cout << "Input calculator was not of type IOscCalcSterile." << std::endl;
    return nullptr; // If the cast failed, calc_sterile should be nullptr anyway
  }

  //---------------------------------------------------------------------------
  IOscCalcSterile* DowncastToSterile(IOscCalc* calc, bool quiet)
  {
    IOscCalcSterile* calc_sterile
            = dynamic_cast<IOscCalcSterile*>(calc);
    if (calc_sterile) return calc_sterile;
    if (!quiet)
      std::cout << "Input calculator was not of type OscCalcSterile." << std::endl;
    return nullptr; // If the cast failed, calc_sterile should be nullptr anyway
  }

}

