#pragma once  

#include <vector>
#include <iostream>

#include <Eigen/Eigen>
using namespace Eigen;

namespace osc {
  template <typename T>
  struct _OscParameters {
    T Dmsq31() const { return dmsq32 + dmsq21; }
    bool operator==(const _OscParameters<T> & rhs) const {
      return dmsq21 == rhs.dmsq21 &&
             dmsq32  == rhs.dmsq32  &&
             th12    == rhs.th12    &&
             th13    == rhs.th13    &&
             th23    == rhs.th23    &&
             deltacp == rhs.deltacp &&
             L       == rhs.L       &&
             rho     == rhs.rho;    }
    // abuse (?) placement new to overwrite using default constructor
    void reset() { new (this) _OscParameters<T>; }
    T dmsq21 = 0;
    T dmsq32 = 0;
    T th12 = 0;
    T th13 = 0;
    T th23 = 0;
    T deltacp = 0;
    double L = 0;
    double rho = 0;
  };
  typedef _OscParameters<double> OscParameters;

  template <typename T>
  inline std::ostream& operator << (std::ostream& os, _OscParameters<T> const & p)
  {
    os << "\tdelta m^2_{21}: " << p.dmsq21
       << "\tdelta m^2_{32}: " << p.dmsq32
       << "\tdelta m^2_{31}: " << p.Dmsq31()
       << "\ttheta_{12}: " << p.th12
       << "\ttheta_{13}: " << p.th13
       << "\ttheta_{23}: " << p.th23
       << "\tdelta CP: " << p.deltacp
       << "\n";
    return os;
  }

  template <typename T>
  struct _OscCache {
    bool operator==(const _OscCache & rhs) const {
      if(!(parameters == rhs.parameters)) return false;
      for(unsigned int i = 0; i < energies.size(); i++)
        if(energies[i] != rhs.energies[i]) return false;
      return true;
    }

    void clear() {
      ENERGIES.resize(0);
      probabilities.resize(0,0);
      parameters.reset();
    }
    std::vector<double> energies;
    ArrayXd ENERGIES;
    Array<T,Dynamic,Dynamic> probabilities;
    _OscParameters<T> parameters;
    unsigned int iter = 0;
  };
}
