#include "OscLib/OscCalcCPT.h"

#include <cassert>
#include <cstdlib>

namespace osc
{

  OscCalcCPT::OscCalcCPT()
  {
    fCalc    = new OscCalcPMNSOpt;
    fBarCalc = new OscCalcPMNSOpt;
    fSigDel  = {};
  }

  OscCalcCPT::OscCalcCPT(IOscCalcAdjustable* calc,
                                     IOscCalcAdjustable* barcalc,
                                     SDMap sigdel)
  {
    fCalc    = calc;
    fBarCalc = barcalc;
    fSigDel  = sigdel;
  }

  OscCalcCPT::~OscCalcCPT()
  {
    delete fCalc;
    delete fBarCalc;
  }

  IOscCalcAdjustable* OscCalcCPT::Copy() const
  {
    return new OscCalcCPT(fCalc->Copy(),
                                fBarCalc->Copy(),
                                fSigDel);
  }

  double OscCalcCPT::P(int flavBefore, int flavAfter, double E)
  {
    assert(flavBefore*flavAfter > 0); // check for matter<->anti-matter
    return ( (flavBefore > 0) ?    fCalc->P(flavBefore, flavAfter, E)
                              : fBarCalc->P(flavBefore, flavAfter, E) );
  }

  // asymmetric setters
  void OscCalcCPT::SetL(double L, ENuSign sign)
  {
    (sign==ENuSign::kNu) ?    fCalc->SetL(L) :
                           fBarCalc->SetL(L) ;
  }

  void OscCalcCPT::SetRho(double rho, ENuSign sign)
  {
    (sign==ENuSign::kNu) ?    fCalc->SetRho(rho) :
                           fBarCalc->SetRho(rho) ;
  }

  void OscCalcCPT::SetDmsq21(double dmsq21, ENuSign sign)
  {
    (sign==ENuSign::kNu) ?    fCalc->SetDmsq21(dmsq21) :
                           fBarCalc->SetDmsq21(dmsq21) ;
  }

  void OscCalcCPT::SetDmsq32(double dmsq32, ENuSign sign)
  {
    (sign==ENuSign::kNu) ?    fCalc->SetDmsq32(dmsq32) :
                           fBarCalc->SetDmsq32(dmsq32) ;
  }

  void OscCalcCPT::SetTh12(double th12, ENuSign sign)
  {
    (sign==ENuSign::kNu) ?    fCalc->SetTh12(th12) :
                           fBarCalc->SetTh12(th12) ;
  }

  void OscCalcCPT::SetTh13(double th13, ENuSign sign)
  {
    (sign==ENuSign::kNu) ?    fCalc->SetTh13(th13) :
                           fBarCalc->SetTh13(th13) ;
  }

  void OscCalcCPT::SetTh23(double th23, ENuSign sign)
  {
    (sign==ENuSign::kNu) ?    fCalc->SetTh23(th23) :
                           fBarCalc->SetTh23(th23) ;
  }

  void OscCalcCPT::SetdCP(double dCP, ENuSign sign)
  {
    (sign==ENuSign::kNu) ?    fCalc->SetdCP(dCP) :
                           fBarCalc->SetdCP(dCP) ;
  }

  // symmetric getters
  double OscCalcCPT::GetL() const
  {
    double L = GetL(ENuSign::kNu);
    assert( L == GetL(ENuSign::kNuBar) );
    return L;
  }

  double OscCalcCPT::GetRho() const
  {
    double rho = GetRho(ENuSign::kNu);
    assert( rho == GetRho(ENuSign::kNuBar) );
    return rho;
  }

  double OscCalcCPT::GetDmsq21() const
  {
    double dmsq21 = GetDmsq21(ENuSign::kNu);
    assert( dmsq21 == GetDmsq21(ENuSign::kNuBar) );
    return dmsq21;
  }

  double OscCalcCPT::GetDmsq32() const
  {
    double dmsq32 = GetDmsq32(ENuSign::kNu);
    assert( dmsq32 == GetDmsq32(ENuSign::kNuBar) );
    return dmsq32;
  }

  double OscCalcCPT::GetTh12() const
  {
    double th12 = GetTh12(ENuSign::kNu);
    assert( th12 == GetTh12(ENuSign::kNuBar) );
    return th12;
  }

  double OscCalcCPT::GetTh13() const
  {
    double th13 = GetTh13(ENuSign::kNu);
    assert( th13 == GetTh13(ENuSign::kNuBar) );
    return th13;
  }

  double OscCalcCPT::GetTh23() const
  {
    double th23 = GetTh23(ENuSign::kNu);
    assert( th23 == GetTh23(ENuSign::kNuBar) );
    return th23;
  }

  double OscCalcCPT::GetdCP() const
  {
    double dCP = GetdCP(ENuSign::kNu);
    assert( dCP == GetdCP(ENuSign::kNuBar) );
    return dCP;
  }

  // asymmetric getters
  double OscCalcCPT::GetL(ENuSign sign) const
  {
    return ( (sign==ENuSign::kNu) ?    fCalc->GetL()
                                  : fBarCalc->GetL() );
  }

  double OscCalcCPT::GetRho(ENuSign sign) const
  {
    return ( (sign==ENuSign::kNu) ?    fCalc->GetRho()
                                  : fBarCalc->GetRho() );
  }

  double OscCalcCPT::GetDmsq21(ENuSign sign) const
  {
    return ( (sign==ENuSign::kNu) ?    fCalc->GetDmsq21()
                                  : fBarCalc->GetDmsq21() );
  }

  double OscCalcCPT::GetDmsq32(ENuSign sign) const
  {
    return ( (sign==ENuSign::kNu) ?    fCalc->GetDmsq32()
                                  : fBarCalc->GetDmsq32() );
  }

  double OscCalcCPT::GetTh12(ENuSign sign) const
  {
    return ( (sign==ENuSign::kNu) ?    fCalc->GetTh12()
                                  : fBarCalc->GetTh12() );
  }

  double OscCalcCPT::GetTh13(ENuSign sign) const
  {
    return ( (sign==ENuSign::kNu) ?    fCalc->GetTh13()
                                  : fBarCalc->GetTh13() );
  }

  double OscCalcCPT::GetTh23(ENuSign sign) const
  {
    return ( (sign==ENuSign::kNu) ?    fCalc->GetTh23()
                                  : fBarCalc->GetTh23() );
  }

  double OscCalcCPT::GetdCP(ENuSign sign) const
  {
    return ( (sign==ENuSign::kNu) ?    fCalc->GetdCP()
                                  : fBarCalc->GetdCP() );
  }


  TMD5* OscCalcCPT::GetParamsHash() const
  {
    TMD5* hash    = fCalc->GetParamsHash();
    TMD5* barhash = fBarCalc->GetParamsHash();

    // don't provide hash unless both sub-cals do
    if ( !(hash && barhash) )
    {
      delete hash;
      delete barhash;
      return 0;
    }

    TMD5* ret = new TMD5;

    //hash together sub-hashes
    ret->Update( (unsigned char*)hash->AsString(), 16);
    ret->Update( (unsigned char*)barhash->AsString(), 16);

    delete hash;
    delete barhash;

    //also hash in class name in case another CPT class has same idea
    ret->Update( (unsigned char*)"OscCalcCPT", 16);

    ret->Final();
    return ret;
  } 

  //---------------------------------------------------------------------------

  const OscCalcCPT* DowncastToCPT(const IOscCalcAdjustable* osc)
  {
    const OscCalcCPT* cpt = dynamic_cast<const OscCalcCPT*>(osc);
    assert( cpt && "Must use OscCalcCPT with CPT FitVars." );
    return cpt;
  }

  OscCalcCPT* DowncastToCPT(IOscCalcAdjustable* osc)
  {
    OscCalcCPT* cpt = dynamic_cast<OscCalcCPT*>(osc);
    assert( cpt && "Must use OscCalcCPT with CPT FitVars." );
    return cpt;
  }

} // namespace
