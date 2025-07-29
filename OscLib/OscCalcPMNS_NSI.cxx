#include "OscLib/OscCalcPMNS_NSI.h"

#include "OscLib/Constants.h"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>

namespace osc
{
  OscCalcPMNS_NSI::OscCalcPMNS_NSI()
    : fMixDirty(true), fDmDirty(true), fPropDirty(true), fEpsDirty(true), fPrevAnti(0)
  {
  }
  
  //---------------------------------------------------------------------------
  // 2016-10-24 - Getting all oscillation parameters (Standard and NSI)
  std::vector<double> OscCalcPMNS_NSI::GetState() const
  {
    std::vector<double> state;
    state.push_back(fL);
    state.push_back(fRho);
    state.push_back(fDmsq21);
    state.push_back(fDmsq32);
    state.push_back(fTh12);
    state.push_back(fTh13);
    state.push_back(fTh23);
    state.push_back(fdCP);
    state.push_back(fEps_ee);
    state.push_back(fEps_emu);
    state.push_back(fEps_etau);
    state.push_back(fEps_mumu);
    state.push_back(fEps_mutau);
    state.push_back(fEps_tautau);
    state.push_back(fDelta_emu);
    state.push_back(fDelta_etau);
    state.push_back(fDelta_mutau);
    
    return state;
  }
  //---------------------------------------------------------------------------
  
  //---------------------------------------------------------------------------
  // 2016-10-24 - Setting all oscillation parameters (Standard and NSI) - or the State
  void OscCalcPMNS_NSI::SetState(std::vector<double> state)
  {
    int iState(0);
    fMixDirty = true ;
    SetL(state[iState++]);
    SetRho(state[iState++]);
    SetDmsq21(state[iState++]);
    SetDmsq32(state[iState++]);
    SetTh12(state[iState++]);
    SetTh13(state[iState++]);
    SetTh23(state[iState++]);
    SetdCP(state[iState++]);
    SetEps_ee(state[iState++]);
    SetEps_emu(state[iState++]);
    SetEps_etau(state[iState++]);
    SetEps_mumu(state[iState++]);
    SetEps_mutau(state[iState++]);
    SetEps_tautau(state[iState++]);
    SetDelta_emu(state[iState++]);
    SetDelta_etau(state[iState++]);
    SetDelta_mutau(state[iState++]);
  }
  //---------------------------------------------------------------------------

  OscCalcPMNS_NSI::~OscCalcPMNS_NSI()
  {
  }

  double OscCalcPMNS_NSI::P(int flavBefore, int flavAfter, double E)
  {
    const int anti = (flavBefore > 0) ? +1 : -1;
    assert(flavAfter/anti > 0);
    if(anti != fPrevAnti) fPropDirty = true;

    int i = -1, j = -1;
    if(abs(flavBefore) == 12) i = 0;
    if(abs(flavBefore) == 14) i = 1;
    if(abs(flavBefore) == 16) i = 2;
    if(abs(flavAfter) == 12) j = 0;
    if(abs(flavAfter) == 14) j = 1;
    if(abs(flavAfter) == 16) j = 2;
    assert(i >= 0 && j >= 0);

    if(fMixDirty){
      fPMNS_NSI.SetMix(fTh12, fTh23, fTh13, fdCP);
      fMixDirty = false;
    }
    if(fDmDirty){
      fPMNS_NSI.SetDeltaMsqrs(fDmsq21, fDmsq32);
      fDmDirty = false;
    }
    if(fEpsDirty){
      fPMNS_NSI.SetNSI(fEps_ee,    fEps_emu,    fEps_etau,
                       fEps_mumu,  fEps_mutau,  fEps_tautau,
                       fDelta_emu, fDelta_etau, fDelta_mutau);
      fEpsDirty = false;
    }


    fPMNS_NSI.ResetToFlavour(i);
    const double Ne = fRho * constants::kZPerA;
    fPMNS_NSI.PropMatter(fL, E, Ne, anti);
    return fPMNS_NSI.P(j);
  }

  //---------------------------------------------------------------------------
  IOscCalcAdjustable* OscCalcPMNS_NSI::Copy() const
  {
    return new OscCalcPMNS_NSI(*this);
  }
  
  //---------------------------------------------------------------------------
  const OscCalcPMNS_NSI* DowncastToNSI(const IOscCalc* calc)
  {
    const OscCalcPMNS_NSI* calc_nsi = dynamic_cast<const OscCalcPMNS_NSI*>(calc);
    if(calc_nsi) return calc_nsi;
    else             std::cout << "Input calculator was not of type OscCalcPMNS_NSI." << std::endl;
    return nullptr; // If the cast failed, calc_nsi should be nullptr anyway (?)                    
  }

  //---------------------------------------------------------------------------                     
  OscCalcPMNS_NSI* DowncastToNSI(IOscCalc* calc)
  {
    OscCalcPMNS_NSI* calc_nsi = dynamic_cast<OscCalcPMNS_NSI*>(calc);
    if(calc_nsi) return calc_nsi;
    else             std::cout << "Input calculator was not of type OscCalcPMNS_NSI." << std::endl;
    return nullptr; // If the cast failed, calc_nsi should be nullptr anyway (?)                    
  }
  
} // namespace
