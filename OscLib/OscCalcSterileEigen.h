#ifndef OSC_OSCCALCSTERILEEIGEN_H
#define OSC_OSCCALCSTERILEEIGEN_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// \file OscCalcSterile.h                                               //
//                                                                      //
// Adapt the PMNS_Sterile calculator to 3+1 oscillations with Eigen     //
// <jhewes15@fnal.gov>                                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "OscLib/IOscCalcSterile.h"
#include <Eigen/Dense>
#include <list>

namespace osc
{
  /// \brief Eigen-based 3+1 sterile oscillation calculator
  ///
  /// Adapt the \ref PMNS_Sterile calculator to Eigen, hardcoded to 3+1
  class OscCalcSterileEigen: public IOscCalcSterile
  {  
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using IOscCalcAdjustable::P;
    OscCalcSterileEigen();
    OscCalcSterileEigen(const OscCalcSterileEigen& calc);
    virtual ~OscCalcSterileEigen();

    virtual IOscCalcAdjustable* Copy() const override;

    virtual void SetAngle(int i, int j, double th) override;
    virtual void SetDelta(int i, int j, double delta) override;
    virtual void SetDm(int i, double dm) override;

    virtual double GetDm(int i) const override
      { return fDm(i-1); }
    virtual double GetAngle(int i, int j) const override
      { return fTheta(i-1, j-1); }
    virtual double GetDelta(int i, int j) const override
      { return fDelta(i-1, j-1); }

    virtual double P(int flavBefore, int flavAfter, double E) override;

    virtual void Print(const std::string& prefix = "") override
    {
      PrintImpl(4, prefix);
    }

    
    /// Set standard 3-flavor parameters
    virtual void SetStdPars();

  protected:
    
    // A shorthand...
    typedef std::complex<double> complex;
    
    virtual void InitializeVectors();
    
    /// Rotate the Hamiltonian by theta_ij and delta_ij
    virtual void RotateH(int i,int j);
    
    /// Build Hms = H*2E, where H is the Hamiltonian in vacuum on flavour basis   
    /// and E is the neutrino energy. This is effectively the matrix of masses squared.
    virtual void BuildHms();
    
    /// Solve the full Hamiltonian for eigenvectors and eigenvalues
    /// @param E - neutrino energy in GeV
    /// @param Ne - electron number density of matter in mole/cm^3
    /// @param anti - +1 = neutrino case, -1 = anti-neutrino case
    virtual void SolveHam(double E, double Ne, int anti);

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
    
    /// Return the probability of final state in flavour flv
    /// @param flv - final flavor (0,1,2) = (nue,numu,nutau)
    virtual double GetP(int flv) const;
    
    /// Erase memory of neutrino propagate and reset neutrino
    /// to pure flavour flv. Preserves values of mixing and mass-splittings
    /// @param flv - final flavor (0,1,2) = (nue,numu,nutau)
    virtual void ResetToFlavour(int flv=1);
    
    int fNumNus;

    double  fCachedNe;      ///< Cached electron density
    double  fCachedE;       ///< Cached neutrino energy
    int     fCachedAnti;    ///< Cached nu/nubar selector
    bool    fBuiltHms;      ///< Tag to avoid rebuilding Hms

    Eigen::Vector4d  fDm;      ///< m^2_i - m^2_1 in vacuum
    Eigen::Matrix4d  fTheta;   ///< theta[i][j] mixing angle
    Eigen::Matrix4d  fDelta;   ///< delta[i][j] CP violating phase
    Eigen::Vector4cd  fNuState; ///< The neutrino current state
    Eigen::Matrix4cd  fHms;     ///< matrix H*2E in eV^2
    Eigen::Matrix4cd  fHmsMat;  ///< matrix H*2E in eV^2, with matter
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix4cd> fEig;  ///< eigen solver

  };

} // namespace

#endif
