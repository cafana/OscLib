
#ifndef OSCLIB_CONSTANTS
#define OSCLIB_CONSTANTS

#include <cmath>

// \file Constants.h
// \brief Keep important physical constants and unit conversions in one place.
//        All values use the natural units system such that hbar = c = 1.
//        All values are taken from PDG, Physical Constants 1.1, 2024.
// \author <cullenms@fnal.gov>

namespace osc {
namespace constants {

/// G_F, the Fermi Constant in GeV^-2.
const double kFermiConstant = 1.1663788e-5;

/// N_A, Avogadro's number in mol^-1 (or if you prefer, electrons/mol when thinking about the matter effect).
const double kAvogadroConstant = 6.02214076e23;

/// c, the speed of light (unitless). Multiply by c in natural units to convert s to m.
const double kSpeedOfLight = 299792458;

/// hbar, reduced Planck's constant (unitless). Multiply by hbar in natural units to convert invese seconds to joules.
const double kReducedPlanckConstant = 1.054571817e-34;

/// e, the elementary charge (J/eV).
const double kElementaryCharge = 1.602176364e-19;

/// Unit conversion: m^-1 to eV.
const double kInversemToeV =
    kSpeedOfLight          // 1 m^-1 = SpeedOfLight s^-1.
  * kReducedPlanckConstant // 1 s^-1 = ReducedPlanckConstant J.
  / kElementaryCharge;     // ElementaryCharge J^-1 = eV.

/// Given a matter density, rho in mol/cm^3, multiply by 2*sqrt(2)*FermiConstant and convert units to eV.
const double kMatterDensityToEffect =
    sqrt(2) * kFermiConstant * 1e-18 // The 1e-18 converts the Fermi constant from GeV^-2 to eV^-2.
  * kAvogadroConstant                    // 1 mol = 1 AvogadroConstant.
  * pow(kInversemToeV*100,3);

/// Unit conversion: GeV to eV.
const double kGeVToeV = 1e9;

/// Unit conversion: km to m.
const double kkmTom = 1e3;

}
}

#endif // OSCLIB_CONSTANTS
