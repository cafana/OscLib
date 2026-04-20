#ifndef OSC_OSCCALCDECAYEIGEN_H
#define OSC_OSCCALCDECAYEIGEN_H

///////////////////////////////////////////////////////////////////////////////	                                                                                  
/// \class OscLib::OscCalcDecayEigen    
//brief Implementation of neutrino decay in a three-neutrino framework.                                                                                  
///                                                                                                                                                               
/// This class expands the PMNS_Fast class including the decay of the                                                                                               
/// second and third mass state of the neutrino through a decay constant                                                                                             
/// \f$\alpha_i=m_i/\tau_i\ (eV^2)\f$, where \f$m_i\f$ is the mass in the                                                                                           
/// restframe and \f$tau_i\f$ is the lifetime in the restframe.                                                                                                      
///                                                                                                                                                                  
/// Reference: https://doi.org/10.1007/JHEP04(2023)090                                                                                                                
///                                                                                                                                                                  
/// Propagation is computed by exponentiating the Hamiltonian directly                                                                                               
/// instead of solving the eigensystem. This is done with the                                                                                                         
/// [Eigen](https://eigen.tuxfamily.org/) library.                                                                                                                    
///                                                                                                                                                                   
/// \author Victor Carretero - vcarretero\@km3net.de   
/// \author Joao Coelho - jcoelho\@apc.in2p3.fr                                                                                                                       
///                                                                                                                                                                   
/// Adapt the Joao's PMNS_Decay to Novasoft, using OscCalcsterileEigen.h                                                                                               
/// as a reference                                                                                                                                                   
/// \author Andrea Barros - acbarros@mail.uniatlantico.edu.co                                                                                                        
/// \author Mario Acero - marioacero@mail.uniatlantico.edu.co                                                                                                        
///                                                                                                                                                                  
///////////////////////////////////////////////////////////////////////////////   
#include "TDirectory.h"
#include <vector>
#include <complex>
#include <list>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "OscLib/IOscCalcDecay.h"
#include "OscLib/OscCalcDecayEigen.h"
#include "OscLib/OscCalcPMNSOptEigen.h"

namespace osc
{
  class OscCalcDecayEigen: public IOscCalcDecay
  {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      using IOscCalcAdjustable::P;
    OscCalcDecayEigen();
    OscCalcDecayEigen(const OscCalcDecayEigen& calc);
    virtual ~OscCalcDecayEigen();
    
    virtual IOscCalcAdjustable* Copy() const override;
    
    void SetNFlavors(int nflavors);
    virtual void SetAngle(int i, int j, double th) override;
    virtual void SetDelta(int i, int j, double delta) override;
    virtual void SetDm(int i, double dm) override;

    virtual void   SetAlpha3(double alpha3) override;
    virtual void   SetAlpha2(double alpha2) override;
    virtual double GetAlpha3() const override;
    virtual double GetAlpha2() const override;
    
    virtual double GetDm(int i) const override
    { return fDm(i-1); }
    virtual double GetAngle(int i, int j) const override
    { return fTheta(i-1, j-1); }
    virtual double GetDelta(int i, int j) const override
    { return fDelta(i-1, j-1); }
    
    virtual double P(int flavBefore, int flavAfter, double E);
    virtual void Print(const std::string& prefix = "") const override
    {
      PrintImpl(3, prefix);
    }

    virtual void SaveTo(TDirectory* dir, const std::string& name) const;
    
    /// Set standard 3-flavor parameters                                            
    virtual void SetStdPars();
    virtual void SetIsNuBar(bool IsNuBar);

  protected:
    
    // A shorthand...                                                               
    typedef std::complex<double> complex;
    
    // Some useful complex numbers                                                 
    
    static const std::complex<double> zero; ///< zero in complex   
    
    static const std::complex<double> one;  ///< one in complex                     
    
    virtual void InitializeVectors();
    
    ///Rotate the Hamiltonian by theta_ij and delta_ij                              
    virtual void RotateH(int i, int j, Eigen::Matrix3cd& Ham);
    
    /// Wrapper to solve non-hermitian matrix eigenvalues.                          
    void complexsolver(const Eigen::Matrix3cd& A, Eigen::Vector3d& w);
    
    /// Build Hms = H*2E, where H is the Hamiltonian in vacuum on flavour basis     
    /// and E is the neutrino energy. This is effectively the matrix of masses squared.                                                                                
    virtual void BuildHms();
    
    /// Solve the full Hamiltonian for eigenvectors and eigenvalues                 
    /// @param E - neutrino energy in GeV                                           
    /// @param Ne - electron number density of matter in mole/cm^3                  
    
    //virtual void SolveHam(double E, double Ne, int anti);
    virtual void SolveHam(double E, double Ne);
    
    /// Propagation with Decay                                                      
    virtual void PropagatePath(double L, double E, double Ne);        
        
    /// Return the probability of final state in flavour flv                        
    /// @param flv - final flavor (0,1,2) = (nue,numu,nutau)                        
    virtual double GetP(int flv) const;
    
    /// Erase memory of neutrino propagate and reset neutrino                       
    /// to pure flavour flv. Preserves values of mixing and mass-splittings         
    /// @param flv - final flavor (0,1,2) = (nue,numu,nutau)                        
    virtual void ResetToFlavour(int flv=1);
    
    int fNumNus;
    double  fCachedNe;      /// Cached electron density                            
    double  fCachedE;       /// Cached neutrino energy                             
    bool    fBuiltHms;      /// Tag to avoid rebuilding Hms                        
    bool fIsNuBar;          /// anti-neutrino flag                                           
    bool fGotES;            /// Tag to avoid recalculating eigensystem

    Eigen::Vector3cd fBuffer; ///Buffer for neutrino state tranformations         
    Eigen::Vector3d fEval;    ///Eigenvalues of the Hamiltonian                      
    Eigen::Matrix3cd fHd;     ///Decay hamiltonian                                   
    Eigen::Vector3d fDm;      ///m^2_i - m^2_1 in vacuum                            
    Eigen::Matrix3d fTheta;   ///theta[i][j] mixing angle                           
    Eigen::Matrix3d fDelta;   ///delta[i][j] CP violating phase                     
    Eigen::Vector3cd fNuState; ///The neutrino current state                      
    Eigen::Matrix3cd fHms;    ///matrix H*2E in eV^2                             
    Eigen::Vector3d fAlpha;   ///alpha parameters                                   
    Eigen::Matrix3cd fHam;    ///Final hamiltonian                                   
    
  };
  
} // namespace                                                                      

#endif



                                 
                                                                            
 
