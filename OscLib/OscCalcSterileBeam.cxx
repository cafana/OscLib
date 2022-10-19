#include "OscLib/OscCalcSterile.h"
#include "OscLib/OscCalcSterileBeam.h"

#include "TMD5.h"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <typeinfo>

namespace osc
{
  OscCalcSterileBeam::OscCalcSterileBeam()
    : OscCalcSterile(), 
    fKaonscale(0), fPionscale(0), fMuonscale(0)
  {
    //fPMNS_Sterile = new PMNS_Sterile(fNFlavors);
  }

  //---------------------------------------------------------------------------
  OscCalcSterileBeam::~OscCalcSterileBeam()
  {
	  std::cout << " ** OscCalcSterileBeam Destructor is called ** \n";
  }

  //---------------------------------------------------------------------------
  OscCalcSterileBeam::OscCalcSterileBeam(const OscCalcSterileBeam& calc)
    : OscCalcSterile(calc)
  {
    std::cout << "copy constructor\n";
    this->fKaonscale = calc.GetKaonScale();
    this->fMuonscale = calc.GetMuonScale();
    this->fPionscale = calc.GetPionScale();
  }

  //---------------------------------------------------------------------------
  void OscCalcSterileBeam::SetKaonScale(double scale)
  {
    fKaonscale=scale;
    //assert(fKaonscale!=0.0);
  }

  //---------------------------------------------------------------------------
  void OscCalcSterileBeam::SetPionScale(double scale)
  {
    fPionscale=scale;
    //assert(fPionscale!=0.0);
  }

  //---------------------------------------------------------------------------
  void OscCalcSterileBeam::SetMuonScale(double scale)
  {
    fMuonscale=scale;
    //assert(fMuonscale!=0.0);
  }

  //---------------------------------------------------------------------------
  double OscCalcSterileBeam::GetKaonScale() const
  {
    return fKaonscale;
  }

  //---------------------------------------------------------------------------
  double OscCalcSterileBeam::GetPionScale() const
  {
    return fPionscale;
  }

  //---------------------------------------------------------------------------
  double OscCalcSterileBeam::GetMuonScale() const
  {
    return fMuonscale;
  }

  //---------------------------------------------------------------------------
  IOscCalcAdjustable* OscCalcSterileBeam::Copy() const
  {
    std::cout << " *** OscCalSterileBeam Copy() is called *** " << std::endl;
    return new OscCalcSterileBeam(*this);
  }

  //---------------------------------------------------------------------------
  TMD5* OscCalcSterileBeam::GetParamsHash() const
  {
    TMD5* ret = new TMD5;                                                            
    std::string txt = "SterileBeam";                                                 
    ret->Update((unsigned char*)txt.c_str(), txt.size());                            
    std::vector<double> buf;                                                         
    buf.push_back(fL);                                                               
    buf.push_back(fRho);                                                             
    buf.push_back(fKaonscale);                                                       
    buf.push_back(fPionscale);                                                       
    buf.push_back(fMuonscale);                                                       
    for(int i = 2; i <= fNFlavors; ++i) buf.push_back(GetDm(i));                     
    for(int j = 2; j <= fNFlavors; ++j) {                                            
      for(int i = 1; i < j; ++i) {                                                   
	buf.push_back(GetAngle(i, j));                                             
	if(i+1 != j) buf.push_back(GetDelta(i, j));                                
      }
    }
    ret->Update((unsigned char*)&buf[0], sizeof(double)*buf.size());                 
    ret->Final();                                                                    
    return ret;                                                                      
  }
  //---------------------------------------------------------------------------
  /* virtual void OscCalcSterileBeam::SetBeamMode(kBeamMode, double scaleK, double scaleP, double scaleM)
  {
    if(kBeamMode=="Kaon") {
      OscCalcSterileBeam::SetKaonScale(scaleK) ; 
      OscCalcSterileBeam::SetPionScale(1.0) ; 
      OscCalcSterileBeam::SetMuonScale(1.0) ;  }
    else if(kBeamMode=="Pion") {
      OscCalcSterileBeam::SetKaonScale(1.0) ; 
      OscCalcSterileBeam::SetPionScale(scaleP); 
      OscCalcSterileBeam::SetMuonScale(1.0) ;  }
    else if(kBeamMode=="Muon") {
      OscCalcSterileBeam::SetKaonScale(1.0) ; 
      OscCalcSterileBeam::SetPionScale(1.0) ; 
      OscCalcSterileBeam::SetMuonScale(scaleM); }
    else {
      OscCalcSterileBeam::SetKaonScale(1.0) ;  
      OscCalcSterileBeam::SetKaonScale(1.0) ;
      OscCalcSterileBeam::SetKaonScale(1.0) ;  }
  }
  
  //---------------------------------------------------------------------------
  virtual double OscCalcSterileBeam::P(int flavBefore, int flavAfter, double E)
  {
    double origProb = OscCalcSterile::P(int flavBefore, int flavAfter, double E);
    double scale = 0.0;
    double prob=0.0;
    if(kBeamMode=="Kaon" || kBeamMode=="Pion" || kBeamMode=="Muon") {
      scale = OscCalcSterileBeam::GetKaonScale() * 
	OscCalcSterileBeam::GetPionScale() * 
	OscCalcSterileBeam::GetMuonScale() ;
      prob=origProb*scale;    
    }
    return prob;
  }
  */
  //---------------------------------------------------------------------------
  const OscCalcSterileBeam* DowncastToSterileBeam(const IOscCalc* calc)
  {
    //const IOscCalc* base = new osc::OscCalcSterileBeam();
    const OscCalcSterileBeam* calc_sterile = dynamic_cast<const OscCalcSterileBeam*>(calc);
    if(calc_sterile) return calc_sterile;
    else             std::cout << "Input calculator was not of type OscCalcSterileBeam." << std::endl;
    return nullptr; // If the cast failed, calc_sterile should be nullptr anyway
  }

  //---------------------------------------------------------------------------
  OscCalcSterileBeam* DowncastToSterileBeam(IOscCalc* calc)
  {
    //IOscCalc* base = new osc::OscCalcSterileBeam();
    //std::cout << "calc's type is: " << typeid(calc).name() << std::endl;
    //if (typeid(calc) == typeid(osc::NoOscillations))
    //  {
    ///	std::cout << "calc is a osc::NoOscillations calculator!" << std::endl;
    //}
    NoOscillations* calc_noosc = dynamic_cast<NoOscillations*>(calc);
    if (calc_noosc) 
      {
	std::cout << "I was successfully cast to NoOscillations" << std::endl;
      }

    OscCalcSterileBeam* calc_sterile = dynamic_cast<OscCalcSterileBeam*>(calc);
    if(calc_sterile) return calc_sterile;
    else             std::cout << "Input calculator was not of type OscCalcSterileBeam." << std::endl;
    return nullptr; // If the cast failed, calc_sterile should be nullptr anyway
  }
} // namespace
