#include "OscCalcPMNSOptEigen.h"
#include <cmath>

//#include "PMNSOpt.cxx"

//#include "zhetrd3.cxx"
//#include "zheevc3.cxx"
//#include "zheevh3.cxx"
//#include "zheevq3.cxx"

namespace osc {

  TMD5* 
  OscCalcPMNSOptEigen::GetParamsHash() const
  {
    std::string txt = "PMNSOptEigen";
    TMD5* ret = new TMD5;
    ret->Update((unsigned char*)txt.c_str(), txt.size());
    const int kNumParams = 8;
    double buf[kNumParams] = {fRho, fL, fDmsq21, fDmsq32,
			      fTh12, fTh13, fTh23, fdCP};
    ret->Update((unsigned char*)buf, sizeof(double)*kNumParams);
    ret->Final();
    return ret;
  }

  inline Eigen::Array33d
  Probabilities(const Eigen::Matrix3cd & to,
					const Eigen::Matrix3cd & from)
  {
    return (to * from.adjoint()).array().abs2();
  }

  void
  OscCalcPMNSOptEigen::SaveLastParams()
  {
    fLastParams.L = fL;
    fLastParams.rho = fRho;
    fLastParams.dmsq21 = fDmsq21;
    fLastParams.dmsq32 = fDmsq32;
    fLastParams.th12 = fTh12;
    fLastParams.th13 = fTh13;
    fLastParams.th23 = fTh23;
    fLastParams.deltacp = fdCP;
  }

  void _sincos(double theta, double & s, double &c)
  {
#ifdef __APPLE_CC__
    __sincos(theta, &s, &c);
#else
    sincos(theta, &s, &c);
#endif
  }

  IOscCalcAdjustable * 
  OscCalcPMNSOptEigen::Copy() const
  {
    OscCalcPMNSOptEigen * ret = new OscCalcPMNSOptEigen(*this); 
    ret->fCache.clear();
    return ret;
  }

  bool
  OscCalcPMNSOptEigen::ParamsAreCached()
  {
    return 
      fDmsq21  == this->fCache.parameters.dmsq21  &&
      fDmsq32  == this->fCache.parameters.dmsq32  &&
      fTh12    == this->fCache.parameters.th12    &&
      fTh13    == this->fCache.parameters.th13    &&
      fTh23    == this->fCache.parameters.th23    &&
      fdCP     == this->fCache.parameters.deltacp &&
      fL       == this->fCache.parameters.L       &&
      fRho     == this->fCache.parameters.rho;   
  }

  inline int
  OscCalcPMNSOptEigen::ChannelCacheIdx(int flavBefore, int flavAfter) const
  {
    // rows in the cache are arranged in the following order
    // 11 21 31 12 22 32 13 23 33 -11 -21 -31 -12 -22 -32 -13 -23 -33
    // where nue   = 1
    //       numu  = 2
    //       nutau = 3
    // and negative is for antineutrino
    int anti = flavBefore / abs(flavBefore);
    int i = (abs(flavBefore) - 12) / 2;
    int j = (abs(flavAfter)  - 12) / 2;
    int idx = (1-anti)/2*9 + (3 * j + i); 
    return idx;
  }


  Eigen::VectorXd
  OscCalcPMNSOptEigen::P(int flavBefore, int flavAfter, const std::vector<double>& E)
  {
    // Are there probabilities cached and can we use them?
    if(fCache.energies.size() != (size_t) fCache.probabilities.cols() &&
       fCache.energies.size() != 0) { // does a cache exist
      if(ParamsAreCached()) {
	if(this->fCache.energies == E)
	  return this->fCache.probabilities.col(ChannelCacheIdx(flavBefore, flavAfter));
      }      
    }
    FillCache(E);
    return this->fCache.probabilities.col(ChannelCacheIdx(flavBefore, flavAfter));
  }

  double
  OscCalcPMNSOptEigen::P(int flavBefore, int flavAfter, double E,
			       bool fast_and_loose)
  {
    if(fast_and_loose) {
      auto e_it = std::find(fCache.energies.begin(),
			    fCache.energies.end(),
			    E);
      return fCache.probabilities(e_it - fCache.energies.begin(),
				  ChannelCacheIdx(flavBefore, flavAfter));
    } 
    else {
      return P(flavBefore, flavAfter, E);
    }
  }

  double 
  OscCalcPMNSOptEigen::P(int flavBefore, int flavAfter, double E)
  {
    // Are there probabilities cached and can we use them?
    if(fCache.energies.size() != (size_t) fCache.probabilities.cols() &&
       fCache.energies.size() != 0) { // does a cache exist
      if(ParamsAreCached()) { // do current params match those cached
	auto e_it = std::find(fCache.energies.begin(),
			      fCache.energies.end(),
			      E);
	if(e_it != fCache.energies.end()) { // is the given energy cached?
	  
	  return this->fCache.probabilities(e_it - fCache.energies.begin(),
					    ChannelCacheIdx(flavBefore, flavAfter));
						   
	}
      }
    }

    
    // If we make it through this logic, just calculate it from scratch	
    const OscParameters params{fDmsq21, fDmsq32,
	fTh12,    fTh13,   fTh23,
	fdCP,     fL,      fRho};

    auto anti = 1;
    if(flavBefore < 0) anti = -1;
    EigenSystem flavor_states{Solve(MatterHamiltonian(E, anti, params))};


    // returns Eigen Matrix of probabilities for each transition
    auto probabilities = Probabilities(PropMatter(flavor_states.vectors, 
						  flavor_states.values,
						  params),
				       flavor_states.vectors);
    return probabilities((abs(flavBefore)-12)/2, (abs(flavAfter)-12)/2);

    /*
    auto row = Eigen::Map<const Eigen::RowVectorXd>(probabilities.data(),
						    probabilities.size());    
    return row(ChannelCacheIdx(flavBefore, flavAfter));	
    */
  }



  void
  OscCalcPMNSOptEigen::FillCache()
  {
    //    std::cout << "Using OscCalcPMNSOptEigen\n";
    //    std::cout << "Filling cache" << std::endl;
    const OscParameters params{fDmsq21, fDmsq32,
	fTh12,    fTh13,   fTh23,
	fdCP,     fL,      fRho};
    /// probabilities returned are in columns of energy and rows of
    /// flavors: 0-2 [nue, numu, nutau]
    ///        : 3-5 [antinue, antinumu, antinutau]
    /// Dimension two is the "from" flavor
    /// Dimension three is the "to" flavor
    /// ie. P(numu->nue) is            probabilities[energy_idx][0][1]
    ///     P(antinumu-> antinutau) is probabilities[energy_idx][4][5]
    Eigen::ArrayXXd cache = Eigen::ArrayXXd::Zero(this->fCache.energies.size(), 18);

    Eigen::Matrix3cd const ham = BuildHam(params);
    // matter hamiltonian is energy and flavor dependent
    for(int anti : {-1, 1}) {

      //#pragma omp parallel for
      for(unsigned int i = 0; i < this->fCache.energies.size(); i++) {
	if(fCache.energies[i] <= 0) continue; 

	// solve for flavor eigenstates. This part is expensive
	EigenSystem const flavor_states =
	  Solve(AddMatterEffects(ham, 
				 this->fCache.energies[i], 
				 anti,
				 params));

	// returns Eigen Matrix of probabilities for each transition
	auto const probabilities = Probabilities(PropMatter(flavor_states.vectors,
						      flavor_states.values,
						      params),
					   flavor_states.vectors);

	// reshape to a row vector
        auto const row = Eigen::Map<const Eigen::Matrix<double, 1, 9> >(probabilities.data(), probabilities.size());

	// insert the matrix of probabilties into cache as a row
	cache.block(i, (1-anti)/2 * 9, 1, 9) = row;// Eigen::Map<const Eigen::Matrix<double, 1, 9> >(probabilities.data(), probabilities.size());// row;
      }
    }

    this->fCache.probabilities = cache;
    this->fCache.parameters = params;
    this->fCache.iter++;

//    PrintMatrixAddresses(fCache.probabilities.data(),
//			 fCache.probabilities.rows(),
//			 fCache.probabilities.cols());
//    PrintMatrixAddresses(this->fCache.probabilities);
  }

  void 
  OscCalcPMNSOptEigen::PrintMatrixAddresses(const Eigen::ArrayXXd & mat) const
  {

    for(int j = 0; j < mat.cols(); j++) {
      for(int i = 0; i < mat.rows(); i++) {
	std::cout << &mat(i, j) << "  ";
      }
      std::cout << std::endl;
    }
  }

  void 
  OscCalcPMNSOptEigen::PrintArrayAddresses(const double * data,
					      int size) const
  {
    for(int i = 0; i < size; i++)
      std::cout << data + i << ", (" << data[i] << ")  ";
    std::cout << std::endl;
  }

  void 
  OscCalcPMNSOptEigen::PrintMatrixAddresses(const double * data,
					       int nrows, int ncols) const
  {
    for(int i = 0; i < nrows; i++) {
      for(int j = 0; j < ncols; j++) {
	std::cout << data + i * ncols + j << " ";
      }
      std::cout << std::endl;
    }
  }

  Eigen::Matrix3cd 
  OscCalcPMNSOptEigen::PropMatter(Eigen::Matrix3cd const & evec,
				     Eigen::Vector3d const & evals,
				     const OscParameters & params) const
  {
    const Eigen::Array3d temp = evals.array()*kKm2eV*(-1)*params.L;
    //Eigen::Array3d temp(evals.array());

    //// propagation phase
    //temp *= (-1) * kKm2eV * params.L; // NOTE
    const Eigen::Array3d SIN = temp.sin();
    const Eigen::Array3d COS = temp.cos();
    Eigen::Vector3cd PPT;
    PPT << 
      //std::complex<double>(std::cos(temp[0]), std::sin(temp[0])),  
      //std::complex<double>(std::cos(temp[1]), std::sin(temp[1])),  
      //std::complex<double>(std::cos(temp[2]), std::sin(temp[2]));  
      std::complex<double>(COS[0], SIN[0]),  
      std::complex<double>(COS[1], SIN[1]), 
      std::complex<double>(COS[2], SIN[2]);

    // propagated 
    return evec * PPT.asDiagonal();
  }


  void
  OscCalcPMNSOptEigen::FillCache(std::vector<double> const &energies)
  {
    this->SetCachedEnergies(energies);
    //    std::cout << "Set Energies " << std::endl;
    this->FillCache();
  }

  void 
  OscCalcPMNSOptEigen::SetCachedEnergies(std::vector<double> const & energies)
  {
    this->fCache.energies = energies;
  }
  

  //  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3cd> 
  /*
  EigenSystem
  OscCalcPMNSOptEigen::SolveZheevh3(const Eigen::Matrix3cd & ham,
				       bool use_eigen)
  {
    std::complex<double> A[3][3];
    std::complex<double> fEvec[3][3];
    std::complex<double> fEval[3];
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
	A[i][j] = ham(j, i);
      }
    }
    zheevh3(A, fEvec, fEval, use_eigen);

    Matrix3cd vectors(3,3);
    Array3d values(3);
    vectors << 
      fEvec[0][0], fEvec[1][0], fEvec[2][0],
      fEvec[0][1], fEvec[1][1], fEvec[2][1],
      fEvec[0][2], fEvec[1][2], fEvec[2][2];
    values << fEval[0], fEval[1], fEval[2];
    
    
    return EigenSystem{vectors, values};
    
  }
  */
  EigenSystem
  OscCalcPMNSOptEigen::Solve(Eigen::Matrix3cd const & ham) const
  {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3cd> eig;
    eig.computeDirect(ham);
    return EigenSystem{eig.eigenvectors(), eig.eigenvalues()};
  }
  
  Eigen::Matrix3cd
  OscCalcPMNSOptEigen::AddMatterEffects(const Eigen::Matrix3cd & ham,
					   const double & E,
					   const int & anti,
					   const OscParameters & params) const
  {
    double lv = 2 * kGeV2eV*E / params.Dmsq31();  // Osc. length in eV^-1
    double kr2GNe = kK2*M_SQRT2*kGf * params.rho/2; // Matter potential in eV
    
    Eigen::Matrix3cd green_eggs(ham);
    green_eggs.triangularView<Eigen::Upper>()/=lv;

    if (anti>0) {
      // FIXME we only transposehere to be consistent with the original implementation 
      // --- resulting in all that mangling of eigenvectors later
      green_eggs.transposeInPlace(); 
      green_eggs(0,0) += kr2GNe;
    }
    else {
      green_eggs.adjointInPlace();
      green_eggs(0,0) -= kr2GNe;
    }
    return green_eggs;
  } 

  Eigen::Matrix3cd
  OscCalcPMNSOptEigen::AddMatterEffectsAnti(const Eigen::Matrix3cd & ham,
					       const double & E,
					       const OscParameters & params) const
  {
    double lv = 2 * kGeV2eV*E / params.Dmsq31();  // Osc. length in eV^-1
    double kr2GNe = kK2*M_SQRT2*kGf * params.rho/2; // Matter potential in eV
    
    Eigen::Matrix3cd green_eggs(ham);
    green_eggs.triangularView<Eigen::Upper>()/=lv;

    green_eggs.adjointInPlace();
    green_eggs(0,0) -= kr2GNe;

    return green_eggs;
  } 


  Eigen::Matrix3cd
  OscCalcPMNSOptEigen::AddMatterEffects(const Eigen::Matrix3cd & ham,
					   const double & E,
					   const OscParameters & params) const
  {
    double lv = 2 * kGeV2eV*E / params.Dmsq31();  // Osc. length in eV^-1
    double kr2GNe = kK2*M_SQRT2*kGf * params.rho/2; // Matter potential in eV
    
    Eigen::Matrix3cd green_eggs(ham);
    green_eggs.triangularView<Eigen::Upper>()/=lv;

    // FIXME we only transposehere to be consistent with the original implementation 
    // --- resulting in all that mangling of eigenvectors later
    green_eggs.transposeInPlace(); 
    green_eggs(0,0) += kr2GNe;

    return green_eggs;
  } 


        


  Eigen::Matrix3cd 
  OscCalcPMNSOptEigen::MatterHamiltonian(const double & E,
					    const int & anti,
					    const OscParameters & params) const
  {
    double lv = 2 * kGeV2eV*E / params.Dmsq31();  // Osc. length in eV^-1
    double kr2GNe = kK2*M_SQRT2*kGf * params.rho/2; // Matter potential in eV

    auto HLV = BuildHam(params);
    HLV.triangularView<Eigen::Upper>()/=lv;

    if (anti>0) {
      HLV.transposeInPlace(); // FIXME we only transposehere to be consistent with the original implementation --- resulting in all that mangling of eigenvectors later
      HLV(0,0) += kr2GNe;
    }
    else {
      HLV.adjointInPlace();
      HLV(0,0) -= kr2GNe;
    }
    return HLV;
  }


  Eigen::Matrix3cd
  OscCalcPMNSOptEigen::BuildHam(const OscParameters & params) const
  {
    // Create temp variables
    double sij, cij, h00, h11, h01;
    std::complex<double> expCP, h02, h12;
  
    // Hamiltonian in mass base. Only one entry is variable.
    h11 = params.dmsq21 / params.Dmsq31();
  
    // Rotate over theta12
    _sincos(params.th12, sij, cij);

    // There are 3 non-zero entries after rephasing so that h22 = 0
    h00 = h11 * sij * sij - 1;
    h01 = h11 * sij * cij;
    h11 = h11 * cij * cij - 1;

    // Rotate over theta13 with deltaCP
    _sincos(params.deltacp, sij, cij);
    expCP = std::complex<double>(cij, -sij);

    _sincos(params.th13, sij, cij);


    // There are 5 non-zero entries after rephasing so that h22 = 0
    h02 = (-h00 * sij * cij) * expCP;
    h12 = (-h01 * sij) * expCP;
    h11 -= h00 * sij * sij;                                         
    h00 *= cij * cij  -  sij * sij;
    h01 *= cij;

    // Finally, rotate over theta23
    _sincos(params.th23, sij, cij);


    // Fill the Hamiltonian rephased so that h22 = -h11
    Eigen::Matrix3cd mHlv;
    mHlv(0,0) = h00 - 0.5 * h11;
    mHlv(1,1) = 0.5 * h11 * (cij * cij - sij * sij)  +  2 * real(h12) * cij * sij;
    mHlv(2,2) = -mHlv(1,1);
           
    mHlv(0,1) = h02 * sij  +  h01 * cij;
    mHlv(0,2) = h02 * cij  -  h01 * sij;
    mHlv(1,2) = h12 - (h11 * cij + 2 * real(h12) * sij) * sij;
  
    return mHlv;
  }
}
