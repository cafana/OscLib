#pragma once

#include <complex>
#include <vector>
#include <iostream>

#include <Eigen/Eigen>

#include "OscLib/OscParameters.h"
#include "OscLib/IOscCalc.h"

using namespace Eigen;

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
    Array<T,Dynamic,1>      P(int flavBefore, int flavAfter, const std::vector<double> &E) override;
    using _IOscCalcAdjustable<T>::P;

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

  private:
    void FillCache(std::vector<double> const & energies); // move back to private
    _OscCache<T> fCache; // move back to private

    TMD5* GetParamsHash() const override;

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
