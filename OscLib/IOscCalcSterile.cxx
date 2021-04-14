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
  const IOscCalcSterile* DowncastToSterile(const IOscCalc* calc)
  {
    const IOscCalcSterile* calc_sterile
            = dynamic_cast<const IOscCalcSterile*>(calc);
    if(calc_sterile) return calc_sterile;
    std::cout << "Input calculator was not of type IOscCalcSterile." << std::endl;
    return nullptr; // If the cast failed, calc_sterile should be nullptr anyway
  }

  //---------------------------------------------------------------------------
  IOscCalcSterile* DowncastToSterile(IOscCalc* calc)
  {
    IOscCalcSterile* calc_sterile
            = dynamic_cast<IOscCalcSterile*>(calc);
    if (calc_sterile) return calc_sterile;
    std::cout << "Input calculator was not of type OscCalcSterile." << std::endl;
    return nullptr; // If the cast failed, calc_sterile should be nullptr anyway
  }

}

