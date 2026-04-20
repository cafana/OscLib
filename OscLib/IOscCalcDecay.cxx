#include "OscLib/IOscCalcDecay.h"

#include <iostream>

namespace osc
{
  //---------------------------------------------------------------------------                                                                                        
  void IOscCalcDecay::SetDmsq21(const double&)
  {
    std::cerr << "Must use SetDm!" << std::endl;
    abort();
  }

  //---------------------------------------------------------------------------                                                                                        
  void IOscCalcDecay::SetDmsq32(const double&)
  {
    std::cerr << "Must use SetDm!" << std::endl;
    abort();
  }
  
  //---------------------------------------------------------------------------                                                                                        
  void IOscCalcDecay::SetTh12(const double&)
  {
    std::cerr << "Must use SetAngle!" << std::endl;
    abort();
  }
  
  //---------------------------------------------------------------------------                                                                                        
  void IOscCalcDecay::SetTh13(const double&)
  {
    std::cerr << "Must use SetAngle!" << std::endl;
    abort();
  }
  
  //---------------------------------------------------------------------------           
  void IOscCalcDecay::SetTh23(const double&)
  {
    std::cerr << "Must use SetAngle!" << std::endl;
    abort();
  }
  
  //---------------------------------------------------------------------------                                                                                        
  void IOscCalcDecay::SetdCP(const double&)
  {
    std::cerr << "Must use SetDelta!" << std::endl;
    abort();
  }
  //---------------------------------------------------------------------------                                                                                        
  void IOscCalcDecay::PrintImpl(int nNus, const std::string& prefix) const
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
  //---------------------------------------------------------------------------                                                                                        
  double OscCalcDecayTrivial::P(int, int, double)
  {
    return 1;
  }
  //---------------------------------------------------------------------------                                                                                        
  IOscCalcAdjustable* OscCalcDecayTrivial::Copy() const
  {
    return new OscCalcDecayTrivial();
  }
  
  //---------------------------------------------------------------------------                                                                                        
  void OscCalcDecayTrivial::Print(const std::string& prefix) const
  {
    std::cout << prefix << "Trivial decay  osc" << std::endl;
  }
  
  //---------------------------------------------------------------------------    
  const IOscCalcDecay* DowncastToDecay(const IOscCalc* calc, bool quiet)
  {
    const IOscCalcDecay* calc_Decay = dynamic_cast<const IOscCalcDecay*>(calc);
    if(calc_Decay) return calc_Decay;
    if (!quiet)
      std::cout << "Input calculator was not of type IOscCalcDecay." << std::endl;
    return nullptr; // If the cast failed, calc_decay should be nullptr anyway
  }

  IOscCalcDecay* DowncastToDecay(IOscCalc* calc, bool quiet)
  {
    IOscCalcDecay* calc_Decay = dynamic_cast<IOscCalcDecay*>(calc);
    if (calc_Decay) return calc_Decay;
    if (!quiet)
      std::cout << "Input calculator was not of type OscCalcDecay." << std::endl;
    return nullptr; // If the cast failed, calc_decay should be nullptr anyway                                                                                      
  }
  
}


