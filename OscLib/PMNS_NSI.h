////////////////////////////////////////////////////////////////////////
/// \class PMNS_NSI
///
/// \brief Implementation of oscillations of neutrinos in matter in a
///        three-neutrino framework with NSI. 
///
/// This class inherits from the PMNSOpt class
///
/// \version $Id: PMNS_NSI.h,v 1.2 2013/04/03 19:59:31 jcoelho Exp $
///
/// @author joao.coelho@tufts.edu
////////////////////////////////////////////////////////////////////////
#ifndef PMNS_NSI_H
#define PMNS_NSI_H
#include <list>
#include <complex>

#include "OscLib/Constants.h"
#include "OscLib/PMNSOpt.h"

namespace osc {
  template <typename T>
  class _PMNS_NSI : public _PMNSOpt<T> {
  public:
    _PMNS_NSI();
    virtual ~_PMNS_NSI();
    
    void SetNSI(const T& eps_ee,     const T& eps_emu,     const T& eps_etau,
                const T& eps_mumu,   const T& eps_mutau,   const T& eps_tautau,
                const T& delta_emu=0, const T& delta_etau=0, const T& delta_mutau=0);

  protected:
    using complex = typename _PMNSOpt<T>::complex;

    /// Solve the full Hamiltonian for eigenvectors and eigenvalues
    /// @param E - neutrino energy in GeV
    /// @param Ne - electron number density of matter in mole/cm^3
    /// @param anti - +1 = neutrino case, -1 = anti-neutrino case
    virtual void SolveHam(double E, double Ne, int anti);
    
    T       fEps_ee;        ///< NSI parameter ee
    T       fEps_mumu;      ///< NSI parameter mumu
    T       fEps_tautau;    ///< NSI parameter tautau
    complex fEps_emu;       ///< NSI parameter emu
    complex fEps_etau;      ///< NSI parameter etau
    complex fEps_mutau;     ///< NSI parameter mutau
    bool    fResetNSI;      ///< True when NSI parameters are changed
  };

  // preserve older behavior
  typedef _PMNS_NSI<double> PMNS_NSI;
}
#endif
////////////////////////////////////////////////////////////////////////
