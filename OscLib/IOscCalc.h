#ifndef IOSCCALCULATOR_H
#define IOSCCALCULATOR_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// \file    IOscCalc.h                                                  //
//                                                                      //
// \brief   General interface to oscillation calculators                //
// \author  Christopher Backhouse - c.backhouse@ucl.ac.uk               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <Eigen/Eigen>

class TMD5;

/// Oscillation probability calculators
namespace osc
{
  /// General interface to oscillation calculators
  template <typename T>
  class _IOscCalc
  {
  public:
    virtual ~_IOscCalc() {}

    virtual _IOscCalc* Copy() const = 0;

    virtual void Print(const std::string& prefix = "") const = 0;

    /// E in GeV; flavors as PDG codes (so, neg==>antinu)
    virtual T P(int flavBefore, int flavAfter, double E) = 0;
    /// Default implementation forwards to non-vector version using a simple
    /// loop. Override if your calculator has a more efficient implementation.
    virtual Eigen::Array<T, Eigen::Dynamic, 1> P(int flavBefore, int flavAfter, const std::vector<double>& E);
    /// Default implementation forawrds to vector<double> version. Override if
    /// your calculator has a more efficient implementation.
    virtual Eigen::Array<T, Eigen::Dynamic, 1> P(int flavBefore, int flavAfter, const Eigen::ArrayXd& E);
      
    /// \brief Use to check two calculators are in the same state
    ///
    /// \return Null means not implemented for this calculator
    virtual TMD5* GetParamsHash() const {return 0;}
  };
  typedef _IOscCalc<double> IOscCalc;

  /// Pass neutrinos through unchanged
  template <typename T>
  class _NoOscillations: public _IOscCalc<T>
  {
  public:
    virtual _IOscCalc<T>* Copy() const override;

    virtual void Print(const std::string& prefix = "") const override;

    using _IOscCalc<T>::P;    
    virtual T P(int from, int to, double /*E*/) override;

    virtual TMD5* GetParamsHash() const override;
  };
  typedef _NoOscillations<double> NoOscillations;

  /// General interface to any calculator that lets you set the parameters
  template <typename T>
  class _IOscCalcAdjustable : public _IOscCalc<T>
  {
    public:
      virtual ~_IOscCalcAdjustable();
      virtual _IOscCalcAdjustable<T>* Copy() const = 0;

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

      /// Default implementation, prints all the parameters included here
      virtual void Print(const std::string& prefix = "") const override;

      /// \brief Invalidate any caching used internally by the calculator.
      ///
      /// Some calculators use a cache that can become stale in ways
      /// that the calculator may not know about (e.g., Stan var clearing).
      /// Default implementation does nothing.
      virtual void InvalidateCache() {}

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
  typedef _IOscCalcAdjustable<double> IOscCalcAdjustable;

  //----------------------------------------------------------------------
  /// Copy parameters from one calculator to another, irrespective of their type
  template <typename T, typename U>
  void CopyParams(const osc::_IOscCalcAdjustable<T> * inCalc,
                  osc::_IOscCalcAdjustable<U> * outCalc);

} // namespace

#endif
