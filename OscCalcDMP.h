#pragma once

#include <complex>
#include <vector>
#include <iostream>

#include <Eigen/Eigen>

#include "OscLib/OscParameters.h"
#include "OscLib/IOscCalc.h"

using namespace Eigen;

  // Unit conversion constants
  //static const double kKm2eV  = 5.06773103202e+09; ///< km to eV^-1
  //static const double kK2     = 4.62711492217e-09; ///< mole/GeV^2/cm^3 to eV
  //static const double kGeV2eV = 1.0e+09;           ///< GeV to eV
  //static const double kGf     = 1.166371e-5; //G_F in units of GeV^-2
  //static const double eVsqkm_to_GeV = 1e-9 / 1.973269681602260e-7 * 1e3; // HS this is more like value in OscLib
  //static const double YerhoE2a = 1.52e-4;



namespace osc
{
  /// \brief Helper struct for the cache. Might not need this

  /// \brief Struct of the cache


  /// \brief A DMP based implementation of \ref OscCalcPMNS
  ///
  /// Uses DMP
  template <typename T>
  class _OscCalcDMP : public _IOscCalcAdjustable<T>
  {
  public:
    _OscCalcDMP() {}
    _OscCalcDMP(std::vector<double> energies) {
      this->SetCachedEnergies(energies);
    }
    ~_OscCalcDMP() = default;

    _IOscCalcAdjustable<T> * Copy() const override;
    
    T P(int flavBefore, int flavAfter, double E) override;
    T P(int flavBefore, int flavAfter, double E, bool fast_and_loose);

    Array<T,Dynamic,Dynamic> P() {return this->fCache.probabilities;}
    Matrix<T,Dynamic,1>      P(int flavBefore, int flavAfter, const std::vector<double> &E) override;

    void SetL     (double L       ) override {SaveLastParams(); this->fL      = L;}
    void SetRho   (double rho     ) override {SaveLastParams(); this->fRho    = rho;}
    void SetDmsq21(const T& dmsq21) override {SaveLastParams(); this->fDmsq21 = dmsq21;}
    void SetDmsq32(const T& dmsq32) override {SaveLastParams(); this->fDmsq32 = dmsq32;}
    void SetTh12  (const T& th12  ) override {SaveLastParams(); this->fTh12   = th12;}
    void SetTh13  (const T& th13  ) override {SaveLastParams(); this->fTh13   = th13;}
    void SetTh23  (const T& th23  ) override {SaveLastParams(); this->fTh23   = th23;}
    void SetdCP   (const T& dCP   ) override {SaveLastParams(); this->fdCP    = dCP;}

    void InvalidateCache() override { this->fCache.clear(); }

    std::string Name() const{ return  name;}
    void FillCache(std::vector<double> const & energies); // move back to private
    _OscCache<T> fCache; // move back to private

  private:
    int ChannelCacheIdx(int flavBefore, int flavAfter) const;
    
    std::string name =  "OscCalcDMP";
    // Fill the cache at the current parameter values 
    virtual void FillCache();

    // update fLastParams with the current parameters before changing them
    void SaveLastParams();
    void SetCachedEnergies(std::vector<double> const & energies);
    bool ParamsAreCached();
    _OscParameters<T> fLastParams;
  };

  typedef _OscCalcDMP<double> OscCalcDMP;
}
