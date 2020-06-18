////////////////////////////////////////////////////////////////////////
/// \class _PMNSOpt
///
/// \brief Implementation of oscillations of neutrinos in matter in a
///        three-neutrino framework.
///
/// Two optimizations are relevant:
/// The construction of the Hamiltonian follows DocDB-XXXX (to be posted)
/// The eigensystem determination is based on the following reference:
///
///......................................................................
///
/// Int. J. Mod. Phys. C       VOLUME 19, NUMBER 03            MARCH 2008
///
///     Efficient numerical diagonalization of hermitian 3x3 matrices
///
///                            Joachim Kopp
///                  Max–Planck–Institut für Kernphysik
///             Postfach 10 39 80, 69029 Heidelberg, Germany
///                    (Received 19 October 2007)
///
///                                523
///......................................................................
///
/// The code structure follows the implementation written by M. Messier
/// in the PMNS class.
///
/// \version $Id: PMNSOpt.h,v 1.1 2013/01/19 16:09:57 jcoelho Exp $
///
/// @author joao.coelho@tufts.edu
////////////////////////////////////////////////////////////////////////
#ifndef PMNSOPT_H
#define PMNSOPT_H
#include <list>
#include <complex>

namespace osc
{

  // Unit conversion constants
  static const double kKm2eV  = 5.06773103202e+09; ///< km to eV^-1
  static const double kK2     = 4.62711492217e-09; ///< mole/GeV^2/cm^3 to eV
  static const double kGeV2eV = 1.0e+09;           ///< GeV to eV

  //G_F in units of GeV^-2
  static const double kGf     = 1.166371e-5;

  /// Optimized version of \ref PMNS
  template <typename T>
  class _PMNSOpt
  {
    public:
      _PMNSOpt();
      virtual ~_PMNSOpt() = default;

      /// Set the parameters of the PMNS matrix
      /// @param th12    - The angle theta_12 in radians
      /// @param th23    - The angle theta_23 in radians
      /// @param th13    - The angle theta_13 in radians
      /// @param deltacp - The CPV phase delta_CP in radians
      virtual void SetMix(const T &th12, const T &th23, const T &th13, const T &deltacp);

      /// Set the mass-splittings
      /// @param dm21 - m2^2-m1^2 in eV^2
      /// @param dm32 - m3^2-m2^2 in eV^2
      virtual void SetDeltaMsqrs(const T &dm21, const T &dm32);

      /// Propagate a neutrino through a slab of matter
      /// @param L - length of slab (km)
      /// @param E - neutrino energy in GeV
      /// @param Ne - electron number density of matter in mole/cm^3
      /// @param anti - +1 = neutrino case, -1 = anti-neutrino case
      virtual void PropMatter(double L, double E, double Ne, int anti=1);
      virtual void PropMatter(const std::list<double>& L,
                              double                   E,
                              const std::list<double>& Ne,
                              int anti);

      /// Propagate a neutrino through vacuum
      /// @param L - length of slab (km)
      /// @param E - neutrino energy in GeV
      /// @param anti - +1 = neutrino case, -1 = anti-neutrino case
      virtual void PropVacuum(double L, double E, int anti=1);

      /// Return the probability of final state in flavour flv
      /// @param flv - final flavor (0,1,2) = (nue,numu,nutau)
      virtual T P(int flv) const;

      /// Erase memory of neutrino propagate and reset neutrino
      /// to pure flavour flv. Preserves values of mixing and mass-splittings
      /// @param flv - final flavor (0,1,2) = (nue,numu,nutau)
      virtual void ResetToFlavour(int flv=1);

    protected:
      // A shorthand...
      typedef std::complex<T> complex;

      /// Build H*lv, where H is the Hamiltonian in vacuum on flavour basis
      /// and lv is the oscillation length
      virtual void BuildHlv();

      /// Solve the full Hamiltonian for eigenvectors and eigenvalues
      /// @param E - neutrino energy in GeV
      /// @param Ne - electron number density of matter in mole/cm^3
      /// @param anti - +1 = neutrino case, -1 = anti-neutrino case
      virtual void SolveHam(double E, double Ne, int anti);

      /// Set the eigensystem to the analytic solution of the vacuum Hamiltonian
      /// @param E - neutrino energy in GeV
      /// @param anti - +1 = neutrino case, -1 = anti-neutrino case
      virtual void SetVacuumEigensystem(double E, int anti);

      T       fDm21;          ///< m^2_2 - m^2_1 in vacuum
      T       fDm31;          ///< m^2_3 - m^2_1 in vacuum
      T       fTheta12;       ///< theta12 mixing angle
      T       fTheta23;       ///< theta23 mixing angle
      T       fTheta13;       ///< theta13 mixing angle
      T       fDeltaCP;       ///< CP violating phase
      complex fHlv[3][3];     ///< dimensionless matrix H*lv
      complex fEvec[3][3];    ///< Eigenvectors of the Hamiltonian
      T       fEval[3];       ///< Eigenvalues of the Hamiltonian
      complex fNuState[3];    ///< The neutrino current state
      double  fCachedNe;      ///< Cached electron density
      double  fCachedE;       ///< Cached neutrino energy
      int     fCachedAnti;    ///< Cached nu/nubar selector
      bool    fBuiltHlv;      ///< Tag to avoid rebuilding Hlv
  };

  // preserve older behavior
  typedef _PMNSOpt<double> PMNSOpt;
}
#endif
////////////////////////////////////////////////////////////////////////
