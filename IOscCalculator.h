#ifndef IOSCCALCULATOR_H
#define IOSCCALCULATOR_H

#include "TMD5.h"
#include <Eigen/Eigen>
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// \file    IOscCalculator.h                                            //
//                                                                      //
// \brief   General interface to oscillation calculators                //
// \author  Christopher Backhouse - bckhouse@caltech.edu                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

/// Oscillation probability calculators
namespace osc
{
  /// General interface to oscillation calculators
  template <typename T>
  class _IOscCalculator
  {
  public:
    virtual ~_IOscCalculator() {}

    virtual _IOscCalculator* Copy() const = 0;

    /// E in GeV; flavors as PDG codes (so, neg==>antinu)
    virtual T P(int flavBefore, int flavAfter, double E) = 0;
    virtual Eigen::Matrix<T,Eigen::Dynamic,1> P(int flavBefore, int flavAfter, const std::vector<double> &E);
      
    /// \brief Use to check two calculators are in the same state
    ///
    /// \return Null means not implemented for this calculator
    virtual TMD5* GetParamsHash() const {return 0;}
  };
  typedef _IOscCalculator<double> IOscCalculator;

  /// Pass neutrinos through unchanged
  template <typename T>
  class _NoOscillations: public _IOscCalculator<T>
  {
  public:
    using _IOscCalculator<T>::P;
    virtual _IOscCalculator<T>* Copy() const override {return new _NoOscillations<T>;}
    
    virtual T P(int from, int to, double /*E*/) override
    {
      if(from == to || to == 0) return 1;
      return 0;
    }

    /// Always compare equal to self
    virtual TMD5* GetParamsHash() const override
    {
      TMD5* ret = new TMD5;
      const char* txt = "NoOscillations";
      ret->Update((unsigned char*)txt, strlen(txt));
      ret->Final();
      return ret;
    }
  };
  typedef _NoOscillations<double> NoOscillations;

  /// General interface to any calculator that lets you set the parameters
  template <typename T>
  class _IOscCalculatorAdjustable : public _IOscCalculator<T>
  {
    public:
      virtual ~_IOscCalculatorAdjustable();
      virtual _IOscCalculatorAdjustable<T>* Copy() const = 0;

      // These setters are left unimplemented here, since calculators may want
      // to compute additional values when these are set.
      virtual void SetL     (double L       ) = 0;
      virtual void SetRho   (double rho     ) = 0;
      virtual void SetDmsq21(const T& dmsq21) = 0;
      virtual void SetDmsq32(const T& dmsq32) = 0;
      virtual void SetTh12  (const T& th12  ) = 0;
      virtual void SetTh13  (const T& th13  ) = 0;
      virtual void SetTh23  (const T& th23  ) = 0;
      virtual void SetdCP   (const T& dCP   ) = 0;

      virtual double GetL     () const { return fL      ; }
      virtual double GetRho   () const { return fRho    ; }
      virtual T      GetDmsq21() const { return fDmsq21 ; }
      virtual T      GetDmsq32() const { return fDmsq32 ; }
      virtual T      GetTh12  () const { return fTh12   ; }
      virtual T      GetTh13  () const { return fTh13   ; }
      virtual T      GetTh23  () const { return fTh23   ; }
      virtual T      GetdCP   () const { return fdCP    ; }

      /// \brief Invalidate any caching used internally by the calculator.
      ///
      /// Some calculators use a cache that can become stale in ways
      /// that the calculator may not know about (e.g., Stan var clearing).
      /// Default implementation does nothing.
      virtual void InvalidateCache()  {};

    protected:
      /// \brief This is only a safe implementation if your calculator only
      /// depends on these eight parameters
      ///
      /// \param txt A string to uniquely identify your calculator class
      TMD5* GetParamsHashDefault(const std::string& txt) const;

      // Set by the user. Generally useful to derived classes
      double fRho; // density (g/cm^3); always double since never fitted
      double fL; // baseline (km);  ditto
      T      fDmsq21;
      T      fDmsq32;
      T      fTh12;
      T      fTh13;
      T      fTh23;
      T      fdCP;
  };
  typedef _IOscCalculatorAdjustable<double> IOscCalculatorAdjustable;

  //----------------------------------------------------------------------
  /// Copy parameters from one calculator to another, irrespective of their type
  template <typename T, typename U>
  void CopyParams(const osc::_IOscCalculatorAdjustable<T> * inCalc,
                  osc::_IOscCalculatorAdjustable<U> * outCalc);

} // namespace

#endif
