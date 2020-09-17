#ifndef PMNSDMP_H
#define PMNSDMP_H

#include "OscLib/OscParameters.h"
#include <Eigen/Eigen>
using namespace Eigen;


namespace osc 
{

    //static double eVsqkm_to_GeV = 1e-9 / 1.97327e-7 * 1e3;
    static double eVsqkm_to_GeV = 1e-9 / 1.973269681602260e-7 * 1e3; // HS this is more like value in OscLib
    static double YerhoE2a = 1.52e-4;

    //struct OscParameters;


    //    static OscParameters p0 = {7.5e-5, 0.002425, 0.6012642166791283, 0.14887328003763659, 0.8354818739782283, 1.5 * M_PI};

    template <typename T>
    class _PMNS_DMP
    {
      public:
        using ParamArray = Array<T,Dynamic,1>;
        using ProbArray = Array<T,Dynamic,Dynamic>;

        _PMNS_DMP(ArrayXd const & energies,
                  double const rho,
                  double const L) :
        _ENERGIES(energies),
        _AA(0.5 * rho * YerhoE2a * energies),
        _L4E(eVsqkm_to_GeV * L / (4*energies)),
        _ONE(ArrayXd::Ones(energies.rows())),
        _Dmsqee(0), _Dmsq21(0), _Dmsq31(0), _c13(0), _s13(0), _c23(0), _s23(0),
        _c213(0), _s213(0), _c212(0), _s212(0),
        _c13sq(0), _s13sq(0), _c23sq(0), _s23sq(0),
        _sd(0), _cd(0)
        {}

        //double PMNS_DMP::P (int flv);
        //double PMNS_DMP::P (int a, int b);

        void SetParameters(_OscParameters<T> const & params);
        void SetParameters(ParamArray const & params);
        ProbArray P(_OscParameters<T> const & params);
        ProbArray P(ParamArray const & params);
        ProbArray P(T dmsq21, T dmsq32, T th12, T th13, T th23, T delta);

        void Print(std::ostream& os) const
        {
          os << "\tDmsqee: " << _Dmsqee
             << "\n";
        }

        T Dmsqee(){return _Dmsqee;}
        T Dmsq21(){return _Dmsq21;}
        T Dmsq31(){return _Dmsq31;}
        T c13(){   return    _c13;}
        T s13(){   return    _s13;}
        T c23(){   return    _c23;}
        T s23(){   return    _s23;}
        T c213(){  return   _c213;}
        T s213(){  return   _s213;}
        T c212(){  return   _c212;}
        T s212(){  return   _s212;}
        T c13sq(){ return  _c13sq;}
        T s13sq(){ return  _s13sq;}
        T c23sq(){ return  _c23sq;}
        T s23sq(){ return  _s23sq;}
        T sd(){    return     _sd;}
        T cd(){    return     _cd;}
        ArrayXd ENERGIES(){    return     _ENERGIES;}
        ArrayXd AA(){    return     _AA;}

      private:

        ArrayXd const _ENERGIES;
        ArrayXd const _AA;
        ArrayXd const _L4E;
        ArrayXd const _ONE;
        T _Dmsqee;
        T _Dmsq21;
        T _Dmsq31;
        T _c13;
        T _s13;
        T _c23;
        T _s23;
        T _c213;
        T _s213;
        T _c212;
        T _s212;
        T _c13sq;
        T _s13sq;
        T _c23sq;
        T _s23sq;
        T _sd;
        T _cd;
    };

    template <typename T>
    inline std::ostream& operator << (std::ostream& os, _PMNS_DMP<T> & p)
    {
        os << "\tDmsqee: " << p.Dmsqee()
           //<< "\tdelta m^2_{32}: " << p.dmsq32
           //<< "\tdelta m^2_{31}: " << p.Dmsq31()
           //<< "\ttheta_{12}: " << p.th12
           //<< "\ttheta_{13}: " << p.th13
           //<< "\ttheta_{23}: " << p.th23
           //<< "\tdelta CP: " << p.deltacp
           << "\n";
      return os;
    }
}
#endif
