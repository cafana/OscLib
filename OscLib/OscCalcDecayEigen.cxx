//////////////////////////////////////////////////////////////////////////////
///
/// Implementation of oscillations of neutrinos in matter in a
/// three-neutrino framework with decay.
///
/// \author Andrea Barros - acbarros@mail.uniatlantico.edu.co
/// \author Mario Acero - marioacero@mail.uniatlantico.edu.co
//////////////////////////////////////////////////////////////////////////////

#include "OscLib/Constants.h"
#include <iostream>
#include <cassert>
#include <stdlib.h>
#include <complex>
#include "TObjString.h"
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include "TMath.h"
#include "OscLib/OscCalcDecayEigen.h"
#include "OscLib/OscCalcPMNSOptEigen.h"

using std::complex;
using std::cout;
using std::endl;
using std::list;

using namespace std;
namespace osc
{
  // Some usefule complex numbers
  const std::complex<double> OscCalcDecayEigen::zero(0, 0);
  const std::complex<double> OscCalcDecayEigen::one(1, 0);

  //---------------------------------------------------------------------------                                                                                        

  OscCalcDecayEigen::OscCalcDecayEigen() : fNumNus(3)
  {
    this->SetStdPars();
    this->ResetToFlavour(1);
    fCachedNe = 0.0;
    fCachedE =  1.0;
  }
  
  //---------------------------------------------------------------------------                                                                                        
  OscCalcDecayEigen::OscCalcDecayEigen(const OscCalcDecayEigen& calc) : fNumNus(3)
 
  {
    fRho        = calc.fRho;
    fL          = calc.fL;
    fCachedNe   = calc.fCachedNe;
    fCachedE    = calc.fCachedE;
    fDm         = calc.fDm;
    fTheta      = calc.fTheta;
    fDelta      = calc.fDelta;
    fNuState    = calc.fNuState;
    fHms        = calc.fHms;
    fHam        = calc.fHam;
    fAlpha      = calc.fAlpha;
    fHd         = calc.fHd;
    fBuffer     = calc.fBuffer;
    fIsNuBar    = calc.fIsNuBar;
  }
  
  //---------------------------------------------------------------------------                                                                                        
  OscCalcDecayEigen::~OscCalcDecayEigen()
  {}
  //---------------------------------------------------------------------------                                                                                        
  IOscCalcAdjustable* OscCalcDecayEigen::Copy() const
  {
    return new OscCalcDecayEigen(*this);
  }
  
  //---------------------------------------------------------------------------                                                                                        
  void OscCalcDecayEigen::InitializeVectors()
  {
    fDm.setZero();
    fTheta.setZero();
    fDelta.setZero();
    
    fNuState.setZero();
    fHms.setZero();
    fHam.setZero();
    fHd.setZero();
    fBuffer.setZero();
    fAlpha.setZero();
    fEval.setZero();                                                                                                                                                
    
  }

  void OscCalcDecayEigen::SaveTo(TDirectory* dir, const std::string& name) const
  {
    dir->cd();
    TObjString("OscCalcDecayEigen").Write(name.c_str());
  }
  
  /// Setting Alpha3 parameter, it must be possitive, and units are eV^2.                                                                                              
  ///                                                                                                                                                                  
  /// @param alpha3 = m3/tau3, mass and lifetime of the 3rd state in the restframe                                                                                     
  ///                                                                                                                                                                  
  void OscCalcDecayEigen::SetAlpha3(double alpha3)
  {
    if (alpha3 < 0) {
      cerr << "WARNING: Alpha3 must be positive. Doing nothing." << endl;
      return;
    }
    
    std::cout << "OscCalcDecay Alpha3 " << alpha3 << std::endl;
    fBuiltHms = fBuiltHms && (fAlpha(2) == alpha3);
    fAlpha(2) = alpha3;
    std::cout << "OscCalcDecay fAlpha(2) = " << fAlpha(2) << std::endl;

  }
  
  //.............................................................................
  ///
  /// Setting Alpha2 parameter, it must be possitive, and units are eV^2.
  ///
  /// @param alpha2 = m2/tau2, mass and lifetime of the 2nd state in the restframe
  ///

  void OscCalcDecayEigen::SetAlpha2(double alpha2)
  {
    if (alpha2 < 0) {
      cerr << "WARNING: Alpha2 must be positive. Doing nothing." << endl;
      return;
    }
    
    std::cout << "OscCalcDecay Alpha2 " << alpha2 << std::endl;
    fBuiltHms = fBuiltHms && (fAlpha(1) == alpha2);
    fAlpha(1) = alpha2;
  }
  
  //---------------------------------------------------------------------------          
  /// Reimplement SetIsNuBar to also rebuild hamiltonian.                                                                                                              
  ///                                                                                                                                                                  
  /// @param isNuBar - true if antineutrinos                                                                                                                           
  
  void OscCalcDecayEigen::SetIsNuBar(bool isNuBar)
  {
    fGotES = fGotES && (fIsNuBar == isNuBar);
    fBuiltHms = fBuiltHms && (fIsNuBar == isNuBar);
    fIsNuBar = isNuBar;
    }
  ////////////////////////////////////////////////////////////////////////////                                                                                                                                                                   
  void OscCalcDecayEigen::SetStdPars()
  {
    this->InitializeVectors();
    
    if(fNumNus>2) {
      this->SetAngle(1,2,0.587252);
      this->SetAngle(1,3,0.147969);
      this->SetAngle(2,3,0.83105);
      this->SetDm(2,7.53e-5);
      this->SetDm(3,2.4333e-3);
      this->SetL(810);
      this->SetRho(2.74);
    }
    else if(fNumNus==2){
      this->SetAngle(1,2,0.7);
      this->SetDm(2,2.4e-3);
    }
    
  }
  
  //---------------------------------------------------------------------------                                                                                        
  //---------------------------------------------------------------------------   
  void OscCalcDecayEigen::SetAngle(int i, int j, double th)
  {
    
    if (i > j) {
      cout << "Fatal Error!" << endl;
      cout << "First argument should be smaller than second argument" << endl;
      cout << "You must set reverse order (Theta" << j << i << "). " << endl;
      abort();
    }
    if (i < 1 || i > fNumNus - 1 || j < 2 || j > fNumNus) {
      cout << "Fatal Error!" << endl;
      cout << "Theta" << i << j << " not valid for " << fNumNus;
      cout << " neutrinos. Aborting" << endl;
      abort();
    }
    // Check if value is actually changing
    fBuiltHms = fBuiltHms && (fTheta(i - 1,j - 1) == th);
    
    fTheta(i-1,j-1) = th;
  }

  //---------------------------------------------------------------------------
  void OscCalcDecayEigen::SetDelta(int i, int j, double delta)
  {

    if(i>j){
      cout << "Fatal Error!" << endl;
      cout << "First argument should be smaller than second argument" << endl;
      cout << "You must set reverse order (Delta" << j << i << "). " << endl;
      abort();
    }
    if(i<1 || i>fNumNus-1 || j<2 || j>fNumNus){
      cout << "Fatal Error!" << endl;
      cout << "Delta" << i << j << " not valid for " << fNumNus;
      cout << " neutrinos. Aborting." << endl;
      abort();
    }
    if(i+1==j){
      cout << "Rotation " << i << j << " is real. Doing nothing." << endl;
      return;
    }
    // Check if value is actually changing
    fBuiltHms = fBuiltHms && (fDelta(i - 1,j - 1) == delta);
    
    fDelta(i-1,j-1) = delta;
  }
  
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  void OscCalcDecayEigen::SetDm(int i, double dm)
  {
    if (i < 2 || i > fNumNus) {
      cout << "Fatal Error!" << endl;
      cout << "Dm" << i << "1 not valid for " << fNumNus;
      cout << " neutrinos. Aborting." << endl;
      abort();
    }
    // Check if value is actually changing
    fBuiltHms = fBuiltHms && (fDm(i - 1) == dm);

    fDm(i-1) = dm;
  }
  
  //---------------------------------------------------------------------------
  //...........................................................................
  ///
  /// Get alpha3 parameter
  /// @return alpha3
  ///
  double OscCalcDecayEigen::GetAlpha3() const { return fAlpha(2); }
  ///
  /// Get alpha2 parameter
  ///
  /// @return alpha2
  ///
  double OscCalcDecayEigen::GetAlpha2() const { return fAlpha(1); } 
  //////////////////////////////////////////////////////////////////

  void OscCalcDecayEigen::RotateH(int i, int j, Eigen::Matrix3cd& Ham)
  {
    // Do nothing if angle is zero                                                                                                                                     
    if (fTheta(i,j) == 0) return;

    double fSinBuffer = sin(fTheta(i,j));
    double fCosBuffer = cos(fTheta(i,j));

    double   fHmsBufferD;
    complex  fHmsBufferC;
    
    // With Delta                                                                                                                                                      
    if (i + 1 < j) {
      complex fExpBuffer = complex(cos(fDelta(i,j)), -sin(fDelta(i,j)));

      // General case                                                                                                                                                  
      if (i > 0) {
        // Top columns                                                                                                                                                 
        for (int k = 0; k < i; k++) {
          fHmsBufferC = Ham(k,i);

          Ham(k,i) *= fCosBuffer;
          Ham(k,i) += Ham(k,j) * fSinBuffer * conj(fExpBuffer);
	  
          Ham(k,j) *= fCosBuffer;
          Ham(k,j) -= fHmsBufferC * fSinBuffer * fExpBuffer;
        }
        
	// Middle row and column                                                                                                                                       
        for (int k = i + 1; k < j; k++) {
          fHmsBufferC = Ham(k,j);
	  
          Ham(k,j) *= fCosBuffer;
          Ham(k,j) -= conj(Ham(i,k)) * fSinBuffer * fExpBuffer;

          Ham(i,k) *= fCosBuffer;
          Ham(i,k) += fSinBuffer * fExpBuffer * conj(fHmsBufferC);
        }
	
        // Nodes ij                                                                                                                                                    
	fHmsBufferC = Ham(i,i);
        fHmsBufferD = real(Ham(j,j));
	
	Ham(i,i) *= fCosBuffer * fCosBuffer;
        Ham(i,i) +=
          2 * fSinBuffer * fCosBuffer * real(Ham(i,j) * conj(fExpBuffer));
        Ham(i,i) += fSinBuffer * Ham(j,j) * fSinBuffer;
	
        Ham(j,j) *= fCosBuffer * fCosBuffer;
        Ham(j,j) += fSinBuffer * fHmsBufferC * fSinBuffer;
        Ham(j,j) -=
          2 * fSinBuffer * fCosBuffer * real(Ham(i,j) * conj(fExpBuffer));
        Ham(i,j) -= 2 * fSinBuffer * real(Ham(i,j) * conj(fExpBuffer)) *
          fSinBuffer * fExpBuffer;
        Ham(i,j) +=
          fSinBuffer * fCosBuffer * (fHmsBufferD - fHmsBufferC) * fExpBuffer;
      }
      // First rotation on j (No top columns)                                                                                                                          
      else {
        // Middle rows and columns                                                                                                                                     
        for (int k = i + 1; k < j; k++) {
          Ham(k,j) = -conj(Ham(i,k)) * fSinBuffer * fExpBuffer;

          Ham(i,k) *= fCosBuffer;
        }

        // Nodes ij                                                                                                                                                    
        fHmsBufferD = real(Ham(i,i));
	
        Ham(i,j) =
          fSinBuffer * fCosBuffer * (Ham(j,j) - fHmsBufferD) * fExpBuffer;

        Ham(i,i) *= fCosBuffer * fCosBuffer;
        Ham(i,i) += fSinBuffer * Ham(j,j) * fSinBuffer;
        
	Ham(j,j) *= fCosBuffer * fCosBuffer;
        Ham(j,j) += fSinBuffer * fHmsBufferD * fSinBuffer;
      }
    }
    // Without Delta (No middle rows or columns: j = i+1)                                                                                                              
    else {
      // General case                                                                                                                                                  
      if (i > 0) {
        // Top columns                                                                                                                                                 
        for (int k = 0; k < i; k++) {
          fHmsBufferC = Ham(k,i);
	  
          Ham(k,i) *= fCosBuffer;
          Ham(k,i) += Ham(k,j) * fSinBuffer;
	  
          Ham(k,j) *= fCosBuffer;
          Ham(k,j) -= fHmsBufferC * fSinBuffer;
        }
	// Nodes ij                                                                                                                                                    
        fHmsBufferC = Ham(i,i);
        fHmsBufferD = real(Ham(j,j));
	
        Ham(i,i) *= fCosBuffer * fCosBuffer;
        Ham(i,i) += 2 * fSinBuffer * fCosBuffer * real(Ham(i,j));
        Ham(i,i) += fSinBuffer * Ham(j,j) * fSinBuffer;
	
        Ham(j,j) *= fCosBuffer * fCosBuffer;
        Ham(j,j) += fSinBuffer * fHmsBufferC * fSinBuffer;
        Ham(j,j) -= 2 * fSinBuffer * fCosBuffer * real(Ham(i,j));

        Ham(i,j) -= 2 * fSinBuffer * real(Ham(i,j)) * fSinBuffer;
        Ham(i,j) += fSinBuffer * fCosBuffer * (fHmsBufferD - fHmsBufferC);
      }
      // First rotation (theta12)                                                                                                                                      
      else {
        Ham(i,j) = fSinBuffer * fCosBuffer * Ham(j,j);

        Ham(i,i) = fSinBuffer * Ham(j,j) * fSinBuffer;
	
        Ham(j,j) *= fCosBuffer * fCosBuffer;
      }
    }
  }
  
  
  //---------------------------------------------------------------------------                                                                                        
  void OscCalcDecayEigen::BuildHms()
  {
    
    // Check if anything changed                                                                                                                                       
    if(fBuiltHms) return;

    // Tag to recompute eigensystem                                                                                                                                    
    fGotES = false;
    for(int j=0; j<fNumNus; j++){
      // Set mass splitting                                                                                                                                            
      fHms(j,j) = fDm(j);
      fHd(j,j) = fAlpha(j);
      // Reset off-diagonal elements                                                                                                                                   
      for(int i=0; i<j; i++){
        fHms(i,j) = 0;
        fHd(i,j) = 0;
      }
      // Rotate j neutrinos                                                                                                                                            
      for(int i=0; i<j; i++){
        RotateH(i,j, fHms);
        RotateH(i,j, fHd);
      }
    }

    // Taking care of antineutrinos delta-> -delta and filling the upper triangle                                                                                      
    // because final hamiltonian will be non-hermitian.    

    for (int i = 0; i < fNumNus; i++) {
      for (int j = i + 1; j < fNumNus; j++) {
        if (fIsNuBar) {
          fHms(i,j) = conj(fHms(i,j));
          fHd(i,j)  = conj(fHd(i,j));
        }
        fHms(j,i) = conj(fHms(i,j));
        fHd(j,i)  = conj(fHd(i,j));
      }
    }
    
    const complex numi(0.0, 1.0);
    /// Construct the total Hms+Hd
    for (int j = 0; j < fNumNus; j++) {
      for (int i = 0; i < fNumNus; i++) {
        fHms(i,j) = fHms(i,j) - numi * fHd(i,j);
      }
    }
    
    // Tag as built                                                                                                                                                   
    fBuiltHms = true;
    
  }
  //---------------------------------------------------------------------------                                                                                        
  ////// Wrapper to solve non-hermitian matrix eigenvalues.                                                                                                           
  ///                                                                                                                                                                 
  /// @param A    - Input matrix                                                                                                                                       
  /// @param w    - Output eigenvalues                                                                                                                                
  ///                                                                                                                                                                  
  void OscCalcDecayEigen::complexsolver(const Eigen::Matrix3cd& A, Eigen::Vector3d& w)
  {
    Eigen::ComplexEigenSolver<Eigen::Matrix3cd> eigensolver;
    eigensolver.compute(A);
    for (int t = 0; t < w.size(); t++) {
      w(t) = eigensolver.eigenvalues()(t).real();
    }
  }
  //                                                                                                                                                                   
  //.............................................................................                                                                                        

  //  void OscCalcDecayEigen::SolveHam(double E, double Ne, int anti)
  void OscCalcDecayEigen::SolveHam(double E, double Ne)
  {
    // Check if anything has changed before recalculating
    if(Ne!=fCachedNe || E!=fCachedE || !fBuiltHms ){
      fCachedNe = Ne;
      fCachedE = E;
      this->BuildHms();
    }
    else return;
    
    double lv = (2 * constants::kGeVToeV * E); // 2E in eV	
    double kr2GNe = constants::kMatterDensityToEffect*fRho/2; // Matter potential in eV	
    
    // Finish building Hamiltonian in matter with dimension of eV	
    for (int i = 0; i < fNumNus; i++) {
      for (int j = 0; j < fNumNus; j++) { fHam(i, j) = fHms(i,j) / lv; }
    }
    
    if (!fIsNuBar){
      fHam(0, 0) += kr2GNe;}
    else{
      fHam(0, 0) -= kr2GNe;}
    
    // Solve Hamiltonian for eigenvalues using the Eigen library method                                                                                                                 
    complexsolver(fHam, fEval);

    fGotES = true;
  }


  //////////////////////////////////////////////////////////////////////////////// 
  //.............................................................................                                                                                      
  ///                                                                                                                                                                  
  /// Propagate the neutrino through a constant density path                                                                                                           
  ///                                                                                                                                                                  
  /// @param p - The neutrino path                                                                                                                                     
  ///                                                                                                                                                                  
  /////////////////////////////////////////////////////////////////////////////                                                                                     

  void OscCalcDecayEigen::PropagatePath(double L, double E, double Ne)
  {                                                                                                                                                 
    // Build Hamiltonian                                                                  
    BuildHms();
    SolveHam(E, Ne);

    double LengthIneV = constants::kkmTom/constants::kInversemToeV * L;

    // Compute evolution operator exp(-I*H*L)                                                                                                                          
    Eigen::Matrix3cd H = fHam;
    H *= complex(0, -LengthIneV);                                                                                                                                   
    H = H.exp();
    
    // Propagate using evolution operator                                                                                                                              
    for (int i = 0; i < fNumNus; i++) {
      fBuffer(i) = 0;
      for (int j = 0; j < fNumNus; j++) { fBuffer(i) += H(i,j) * fNuState(j); }
    }
    
    // Update neutrino state                                                                                                                                           
    for (int i = 0; i < fNumNus; i++) { fNuState(i) = fBuffer(i); }
  }
  //---------------------------------------------------------------------------                                                                                         
  void OscCalcDecayEigen::ResetToFlavour(int flv)
  {
    assert(flv >= 0 && flv < fNumNus);
    for (int i = 0; i < fNumNus; ++i) {
      if (i==flv) fNuState(i) = one;
      else        fNuState(i) = zero;
    }
  }
   
  double OscCalcDecayEigen::GetP(int flv) const
  {
    assert(flv >= 0 && flv < fNumNus);
    return norm(fNuState(flv));
  } 
  //--------------------------------------------------------------------------- 
  
  ///////////////////////////////////////////////////////////////////////////
  double OscCalcDecayEigen::P(int flavBefore, int flavAfter, double E)

  {
    int i = -1, j = -1;
    if(abs(flavBefore) == 12) i = 0;
    if(abs(flavBefore) == 14) i = 1;
    if(abs(flavBefore) == 16) i = 2;
    if(abs(flavAfter) == 12) j = 0;
    if(abs(flavAfter) == 14) j = 1;
    if(abs(flavAfter) == 16) j = 2;
    if(abs(flavAfter) == 0)  j = 3;
    assert(i >= 0 && j >= 0);
    ResetToFlavour(i);
    PropagatePath(fL, E, fRho);
    if (j == 3) return GetP(0) + GetP(1) + GetP(2);
    return GetP(j);
  }
  
} // namespace                                                                                                                                                         
///////////////////////////////////////////////////////////////////////////////                                                                                        








 


    
