#include "OscLib/OscCalcSterileEigen.h"
#include "OscLib/OscCalcPMNSOptEigen.h"
#include "OscLib/Constants.h"

#include <iostream>
#include <cassert>
#include <stdlib.h>

using std::complex;
using std::cout;
using std::endl;
using std::list;

namespace osc
{

  //---------------------------------------------------------------------------
  OscCalcSterileEigen::OscCalcSterileEigen() : fNumNus(4)
  {
    this->SetStdPars();
    this->ResetToFlavour(1);
    fCachedNe = 0.0;
    fCachedE =  1.0;
    fCachedAnti = 1;
  }

  //---------------------------------------------------------------------------
  OscCalcSterileEigen::OscCalcSterileEigen(const OscCalcSterileEigen& calc)
    : fNumNus(4)
  {
    fRho        = calc.fRho;
    fL          = calc.fL;
    fCachedNe   = calc.fCachedNe;
    fCachedE    = calc.fCachedE;
    fCachedAnti = calc.fCachedAnti;
    fDm         = calc.fDm;
    fTheta      = calc.fTheta;
    fDelta      = calc.fDelta;
    fNuState    = calc.fNuState;
    fHms        = calc.fHms;
    fHmsMat     = calc.fHmsMat;
  }

  //---------------------------------------------------------------------------
  OscCalcSterileEigen::~OscCalcSterileEigen()
  {}

  //---------------------------------------------------------------------------
  IOscCalcAdjustable* OscCalcSterileEigen::Copy() const
  {
    return new OscCalcSterileEigen(*this);
  }

  //---------------------------------------------------------------------------
  void OscCalcSterileEigen::InitializeVectors()
  {
    fDm.setZero();
    fTheta.setZero();
    fDelta.setZero();
    
    fNuState.setZero();
    fHms.setZero();
    fHmsMat.setZero();
  }

  //---------------------------------------------------------------------------
  void OscCalcSterileEigen::SetStdPars()
  {

    this->InitializeVectors();

    if(fNumNus>2) {
      this->SetAngle(1,2,0.6);
      this->SetAngle(1,3,0.16);
      this->SetAngle(2,3,0.7);
      this->SetDm(2,7.6e-5);
      this->SetDm(3,2.4e-3);
    }
    else if(fNumNus==2){
      this->SetAngle(1,2,0.7);
      this->SetDm(2,2.4e-3);
    }
   
  }

  //--------------------------------------------------------------------------- 
  void OscCalcSterileEigen::SetAngle(int i, int j, double th) 
  {

    if (i > j) {
      cout << "First argument should be smaller than second argument" << endl;
      cout << "Setting reverse order (Theta" << j << i << "). " << endl;
      int temp = i;
      i = j;
      j = temp;
    }
    if (i < 1 || i > fNumNus - 1 || j < 2 || j > fNumNus) {
      cout << "Theta" << i << j << " not valid for " << fNumNus;
      cout << " neutrinos. Doing nothing." << endl;
      return;
    }

    fTheta(i-1,j-1) = th;

    fBuiltHms = false;
  }

  //---------------------------------------------------------------------------
  void OscCalcSterileEigen::SetDelta(int i, int j, double delta) 
  {

    if(i>j){
      cout << "First argument should be smaller than second argument" << endl;
      cout << "Setting reverse order (Delta" << j << i << "). " << endl;
      int temp = i;
      i = j;
      j = temp;
    }
    if(i<1 || i>fNumNus-1 || j<2 || j>fNumNus){
      cout << "Delta" << i << j << " not valid for " << fNumNus;
      cout << " neutrinos. Doing nothing." << endl;
      return;
    }
    if(i+1==j){
      cout << "Rotation " << i << j << " is real. Doing nothing." << endl;
      return;
    }

    fDelta(i-1,j-1) = delta;

    fBuiltHms = false;
  }

  //---------------------------------------------------------------------------
  void OscCalcSterileEigen::SetDm(int i, double dm) 
  {
    if (i < 2 || i > fNumNus) {
      cout << "Dm" << i << "1 not valid for " << fNumNus;
      cout << " neutrinos. Doing nothing." << endl;
      return;
    }

    fDm(i-1) = dm;

    fBuiltHms = false;
  }

  //---------------------------------------------------------------------------
  void OscCalcSterileEigen::RotateH(int i,int j){

    // Do nothing if angle is zero
    if(fTheta(i,j)==0) return;

    double fSinBuffer(0), fCosBuffer(0);
    _sincos(fTheta(i,j), fSinBuffer, fCosBuffer);
    
    double  fHmsBufferD;
    complex fHmsBufferC;

    // With Delta
    if(i+1<j){
      double fSinDelta(0), fCosDelta(0);
      _sincos(fDelta(i,j),fSinDelta, fCosDelta);
      complex fExpBuffer = complex(fCosDelta, -fSinDelta);

      // General case
      if(i>0){
        // Top columns
        for(int k=0; k<i; k++){
          fHmsBufferC = fHms(k,i);

          fHms(k,i) *= fCosBuffer;
          fHms(k,i) += fHms(k,j) * fSinBuffer * conj(fExpBuffer);

          fHms(k,j) *= fCosBuffer;
          fHms(k,j) -= fHmsBufferC * fSinBuffer * fExpBuffer;
        }

        // Middle row and column
        for(int k=i+1; k<j; k++){
          fHmsBufferC = fHms(k,j);
      
          fHms(k,j) *= fCosBuffer;
          fHms(k,j) -= conj(fHms(i,k)) * fSinBuffer * fExpBuffer;

          fHms(i,k) *= fCosBuffer;
          fHms(i,k) += fSinBuffer * fExpBuffer * conj(fHmsBufferC);
        }

        // Nodes ij
        fHmsBufferC = fHms(i,i);
        fHmsBufferD = real(fHms(j,j));

        fHms(i,i) *= fCosBuffer * fCosBuffer;
        fHms(i,i) += 2 * fSinBuffer * fCosBuffer * real(fHms(i,j) * conj(fExpBuffer));
        fHms(i,i) += fSinBuffer * fHms(j,j) * fSinBuffer;

        fHms(j,j) *= fCosBuffer * fCosBuffer;
        fHms(j,j) += fSinBuffer * fHmsBufferC * fSinBuffer;
        fHms(j,j) -= 2 * fSinBuffer * fCosBuffer * real(fHms(i,j) * conj(fExpBuffer));
    
        fHms(i,j) -= 2 * fSinBuffer * real(fHms(i,j) * conj(fExpBuffer)) * fSinBuffer * fExpBuffer;
        fHms(i,j) += fSinBuffer * fCosBuffer * (fHmsBufferD - fHmsBufferC) * fExpBuffer;

      }
      // First rotation on j (No top columns)
      else{
        // Middle rows and columns
        for(int k=i+1; k<j; k++){
          fHms(k,j) = -conj(fHms(i,k)) * fSinBuffer * fExpBuffer;

          fHms(i,k) *= fCosBuffer;
        }

        // Nodes ij
        fHmsBufferD = real(fHms(i,i));

        fHms(i,j) = fSinBuffer * fCosBuffer * (fHms(j,j) - fHmsBufferD) * fExpBuffer;

        fHms(i,i) *= fCosBuffer * fCosBuffer;
        fHms(i,i) += fSinBuffer * fHms(j,j) * fSinBuffer;

        fHms(j,j) *= fCosBuffer * fCosBuffer;
        fHms(j,j) += fSinBuffer * fHmsBufferD * fSinBuffer;
      }
    
    }
    // Without Delta (No middle rows or columns: j = i+1)
    else{
      // General case
      if(i>0){
        // Top columns
        for(int k=0; k<i; k++){
          fHmsBufferC = fHms(k,i);

          fHms(k,i) *= fCosBuffer;
          fHms(k,i) += fHms(k,j) * fSinBuffer;

          fHms(k,j) *= fCosBuffer;
          fHms(k,j) -= fHmsBufferC * fSinBuffer;
        }

        // Nodes ij
        fHmsBufferC = fHms(i,i);
        fHmsBufferD = real(fHms(j,j));

        fHms(i,i) *= fCosBuffer * fCosBuffer;
        fHms(i,i) += 2 * fSinBuffer * fCosBuffer * real(fHms(i,j));
        fHms(i,i) += fSinBuffer * fHms(j,j) * fSinBuffer;

        fHms(j,j) *= fCosBuffer * fCosBuffer;
        fHms(j,j) += fSinBuffer * fHmsBufferC * fSinBuffer;
        fHms(j,j) -= 2 * fSinBuffer * fCosBuffer * real(fHms(i,j));
    
        fHms(i,j) -= 2 * fSinBuffer * real(fHms(i,j)) * fSinBuffer;
        fHms(i,j) += fSinBuffer * fCosBuffer * (fHmsBufferD - fHmsBufferC);

      }
      // First rotation (theta12)
      else{

        fHms(i,j) = fSinBuffer * fCosBuffer * fHms(j,j);

        fHms(i,i) = fSinBuffer * fHms(j,j) * fSinBuffer;

        fHms(j,j) *= fCosBuffer * fCosBuffer;
    
      }
    }

  }

  //---------------------------------------------------------------------------
  void OscCalcSterileEigen::BuildHms() 
  {

    // Check if anything changed
    if(fBuiltHms) return;

    //fHms.setZero(fNumNus, fNumNus);
    
    for(int j=0; j<fNumNus; j++){
      // Set mass splitting
      fHms(j,j) = fDm(j);
      // Reset off-diagonal elements
      for(int i=0; i<j; i++){
        fHms(i,j) = 0;
      }      
      // Rotate j neutrinos
      for(int i=0; i<j; i++){
        this->RotateH(i,j);
      }
    }

    // Tag as built
    fBuiltHms = true;

  }

  //---------------------------------------------------------------------------
  void OscCalcSterileEigen::SolveHam(double E, double Ne, int anti)
  {

    // Check if anything has changed before recalculating
    if(Ne!=fCachedNe || E!=fCachedE || anti!=fCachedAnti || !fBuiltHms ){
      fCachedNe = Ne;
      fCachedE = E;
      fCachedAnti = anti;
      this->BuildHms();
    }
    else return;

    double lv = 2 * constants::kGeVToeV * E;
    double kr2GNe = constants::kMatterDensityToEffect * Ne;

    fHmsMat = fHms;
    
    fHmsMat.triangularView<Eigen::Upper>() /= lv;
    if (anti > 0) {
      fHmsMat.adjointInPlace();
      fHmsMat(0,0) += kr2GNe;
      fHmsMat(3,3) += kr2GNe/2;
    }
    else {
      fHmsMat.transposeInPlace();
      fHmsMat(0,0) += -kr2GNe;
      fHmsMat(3,3) += -kr2GNe/2;
    }
    
    fEig.compute(fHmsMat);
  }

  //---------------------------------------------------------------------------
  void OscCalcSterileEigen::PropMatter(double L, double E, double Ne, int anti) 
  {
    // Solve Hamiltonian
    this->SolveHam(E, Ne, anti);

    fNuState = fEig.eigenvectors()*(
  				  fEig.eigenvalues().unaryExpr([L] (double x) { double sinx(0), cosx(0); _sincos(-constants::kkmTom/constants::kInversemToeV*L*x,sinx,cosx); return complex(cosx, sinx);}
  				  ).asDiagonal())*fEig.eigenvectors().adjoint()*fNuState;
  }

  //---------------------------------------------------------------------------
  void OscCalcSterileEigen::PropMatter(const list<double>& L,
                                       double E,
                                       const list<double>& Ne,
                                       int anti)
  {
    if (L.size()!=Ne.size()) abort();
    list<double>::const_iterator Li  (L.begin());
    list<double>::const_iterator Lend(L.end());
    list<double>::const_iterator Ni  (Ne.begin());
    for (; Li!=Lend; ++Li, ++Ni) {
      this->PropMatter(*Li, E, *Ni, anti);
    }
  }

  //---------------------------------------------------------------------------
  double OscCalcSterileEigen::GetP(int flv) const
  {
    assert(flv >= 0 && flv < fNumNus);
    return norm(fNuState.coeff(flv));
  }

  //---------------------------------------------------------------------------
  void OscCalcSterileEigen::ResetToFlavour(int flv) 
  {
    for (int i = 0; i < fNumNus; ++i) {
      if (i==flv) fNuState(i) = one;
      else        fNuState(i) = zero;
    }
  }

  //---------------------------------------------------------------------------
  double OscCalcSterileEigen::P(int flavBefore, int flavAfter, double E)
  {
    const int anti = (flavBefore > 0) ? +1 : -1;
    //anti must be +/- 1 but flavAfter can be zero
    assert(flavAfter/anti >= 0);
    
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
    PropMatter(fL, E, fRho/2, anti);
    if (j == 3) return GetP(0) + GetP(1) + GetP(2);
    else return GetP(j);
  }

} // namespace
