#ifndef OSC_OSCCALCULATORPMNSOPTEIGEN_H
#define OSC_OSCCALCULATORPMNSOPTEIGEN_H

#include "OscLib/IOscCalc.h"
#include "OscLib/OscParameters.h"
#include <complex>
#include <vector>


#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

#include <iostream>


  // Some useful complex numbers
  static std::complex<double> zero(0.0,0.0);
  static std::complex<double> one (1.0,0.0);

  // Unit conversion constants
  static const double kKm2eV  = 5.06773103202e+09; ///< km to eV^-1
  static const double kK2     = 4.62711492217e-09; ///< mole/GeV^2/cm^3 to eV
  static const double kGeV2eV = 1.0e+09;           ///< GeV to eV

  //G_F in units of GeV^-2
  static const double kGf     = 1.166371e-5;

namespace osc
{
  /// \brief Helper struct for the cache. Might not need this
//  struct OscParameters {
//    double Dmsq31() const { return dmsq32 + dmsq21; }
//    bool operator==(const OscParameters & rhs) const {
//      return dmsq21 == rhs.dmsq21 &&
//	dmsq32  == rhs.dmsq32  &&
//	th12    == rhs.th12    &&
//	th13    == rhs.th13    &&
//	th23    == rhs.th23    &&
//	deltacp == rhs.deltacp &&
//	L       == rhs.L       &&
//	rho     == rhs.rho;    }
//    double dmsq21 = 0;
//    double dmsq32 = 0;
//    double th12 = 0;
//    double th13 = 0;
//    double th23 = 0;
//    double deltacp = 0;
//    double L = 0;
//    double rho = 0;
//  };
//
//  /// \brief Struct of the cache
//  struct OscCache {    
//    bool operator==(const OscCache & rhs) const {
//      if(!(parameters == rhs.parameters)) return false;
//      for(unsigned int i = 0; i < energies.size(); i++) 
//	if(energies[i] != rhs.energies[i]) return false;
//      return true;
//    }    
//        
//    std::vector<double> energies;
//    Eigen::ArrayXXd probabilities;
//    OscParameters parameters;
//    unsigned int iter = 0;
//  };

  /// \brief Wrapper of a system of eigenvectors and eigenvalues
  struct EigenSystem {
    const Eigen::Matrix3cd vectors;
    const Eigen::Vector3d values;
  };

  /// \brief A re-optimized version of \ref OscCalcPMNSOpt
  ///
  /// Uses a faster caching scheme than OscCalcPMNSOpt
  class OscCalcPMNSOptEigen : public _IOscCalcAdjustable<double>
  {
  public:
    OscCalcPMNSOptEigen() {}
    OscCalcPMNSOptEigen(std::vector<double> energies) {
      this->SetCachedEnergies(energies);
    }
    ~OscCalcPMNSOptEigen() = default;

    IOscCalcAdjustable * Copy() const override;
    
    double P(int flavBefore, int flavAfter, double E) override;
    double P(int flavBefore, int flavAfter, double E, bool fast_and_loose);
    Eigen::VectorXd P(int flavBefore, int flavAfter, 
		     const std::vector<double>& E) override;
    

    void SetL     (double L            ) override {SaveLastParams(); fL      = L;}
    void SetRho   (double rho          ) override {SaveLastParams(); fRho    = rho;}
    void SetDmsq21(const double& dmsq21) override {SaveLastParams(); fDmsq21 = dmsq21;}
    void SetDmsq32(const double& dmsq32) override {SaveLastParams(); fDmsq32 = dmsq32;}
    void SetTh12  (const double& th12  ) override {SaveLastParams(); fTh12   = th12;}
    void SetTh13  (const double& th13  ) override {SaveLastParams(); fTh13   = th13;}
    void SetTh23  (const double& th23  ) override {SaveLastParams(); fTh23   = th23;}
    void SetdCP   (const double& dCP   ) override {SaveLastParams(); fdCP    = dCP;}

    TMD5* GetParamsHash() const override;
    std::string Name() const { return  name; }

  private:
    void FillCache(const std::vector<double> & energies);// move back to private
    _OscCache<double> fCache; // move back to private

    std::string name =  "OscCalcPMNSOptEigen";
    int ChannelCacheIdx(int flavBefore, int flavAfter) const;
    
    // Fill the cache at the current parameter values 
    void FillCache();


    // update fLastParams with the current parameters before changing them
    void SaveLastParams();

    void SetCachedEnergies(std::vector<double> const & energies);
    
    // Build a flavor basis hamiltonian at input parameters. No matter effects
    Eigen::Matrix3cd BuildHam(const OscParameters & params) const;

    // Add matter effects onto a flavor basis vacuum hamiltonian
    Eigen::Matrix3cd AddMatterEffects(const Eigen::Matrix3cd & ham,
				      const double & E,
				      const int & anti, 
				      const OscParameters & params) const;
    Eigen::Matrix3cd AddMatterEffects(const Eigen::Matrix3cd & ham,
				      const double & E,
				      const OscParameters & params) const;
    Eigen::Matrix3cd AddMatterEffectsAnti(const Eigen::Matrix3cd & ham,
					  const double & E,
					  const OscParameters & params) const;

    // Build a flavor basis hamiltonian at input parameters. With matter effects
    Eigen::Matrix3cd MatterHamiltonian(const double & E, 
				       const int & anti,
				       const OscParameters & params) const;
    
    // Solve the eigenvalue problem for the input hamiltonian
    EigenSystem Solve(Eigen::Matrix3cd const & ham) const;

    // rotate eigenvectors through matter
    Eigen::Matrix3cd PropMatter(Eigen::Matrix3cd const & evec,
				Eigen::Vector3d const & evals,
				const OscParameters & params) const;

    // Calculate a matrix of probabilities corresponding to all 
    // permutations of to and from states.
    // Note you probably want "to" to be propagated
    //Eigen::Array33d Probabilities(const Eigen::Matrix3cd & to,
				  //const Eigen::Matrix3cd & from) const;
//    Eigen::Array33d Probabilities(Eigen::Matrix3cd const & evec,
//				  Eigen::Vector3d const & evals,
//				  const OscParameters & params);


    void PrintArrayAddresses(const double * data, int size) const;
    void PrintMatrixAddresses(const double * data, 
			      int nrows,
			      int ncols) const;
    void PrintMatrixAddresses(const Eigen::ArrayXXd & mat) const;
    bool ParamsAreCached();
    
    OscParameters fLastParams;
    
  };
  
  

}


#endif
