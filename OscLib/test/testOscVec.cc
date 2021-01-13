#include <cmath>
#include <vector>

#include "OscLib/OscCalc.h"
#include "OscLib/OscCalcGeneral.h"
#include "OscLib/OscCalcPMNS.h"
#include "OscLib/OscCalcPMNSOpt.h"
#include "OscLib/OscCalcDMP.h"
#include "OscLib/OscCalcAnalytic.h"

#include "TMath.h"

#include <iostream>
#include <fenv.h>

int main()
{
  feenableexcept(FE_INVALID); // Spot any infs or nans early

  osc::OscCalc osc1;
  osc::OscCalcGeneral osc2;
  osc::OscCalcPMNS osc3;
  osc::OscCalcPMNSOpt osc4;
  osc::OscCalcAnalytic osc5;
  osc::OscCalcDMP osc6;

  const std::vector<osc::IOscCalcAdjustable*> oscs {&osc1, &osc2, &osc3, &osc4, &osc5, &osc6};

  const std::vector<std::string> names {"Approx", "General", "PMNS", "PMNSOpt", "Analytic", "DMP"};

  const double L = 800;
  const double rho = 3;
  const double dmsq21 = 7.6e-5;
  const double dmsq32 = 2.35e-3;
  const double th12 = 1;
  const double th13 = .15;
  const double th23 = TMath::Pi()/4-.05;
  const double delta = 1.6*TMath::Pi();

  std::vector<double> Es;
  for(double E = .1; E < 5; E *= 1.01) Es.push_back(E);
  Eigen::ArrayXd EsEig;
  EsEig.resize(Es.size());
  for(int i = 0; i < EsEig.size(); ++i) EsEig[i] = Es[i];

  for(int hie = -1; hie <= +1; hie += 2){
    for(std::size_t n = 0; n < oscs.size(); ++n){
      oscs[n]->SetL(L);
      oscs[n]->SetRho(rho);
      oscs[n]->SetDmsq21(dmsq21);
      oscs[n]->SetDmsq32(hie*dmsq32);
      oscs[n]->SetTh12(th12);
      oscs[n]->SetTh13(th13);
      oscs[n]->SetTh23(th23);
      oscs[n]->SetdCP(delta);
    }

    for(int anti = -1; anti <= +1; anti += 2){
      for(int from = 12; from <= 16; from += 2){
        for(int to = 12; to <= 16; to += 2){
          std::cout << "Calculating for " << anti*from << " to " << anti*to << " hierachy " << hie << std::endl;

          for(unsigned int oscIdx = 0; oscIdx < oscs.size(); ++oscIdx){
            std::cout << "  " << names[oscIdx] << std::endl;
            osc::IOscCalcAdjustable* osc = oscs[oscIdx];

            // Calculate one-by-one
            std::vector<double> Ps;
            for(double E: Es) Ps.push_back(osc->P(anti*from, anti*to, E));

            // Two methods of vector evaluation
            const Eigen::ArrayXd PsVec = osc->P(anti*from, anti*to, Es);
            const Eigen::ArrayXd PsEig = osc->P(anti*from, anti*to, EsEig);

            Eigen::ArrayXd PsEigLayers;
            if(osc::OscCalcAnalytic* oca = dynamic_cast<osc::OscCalcAnalytic*>(osc)){
              std::vector<osc::Layer> layers;
              layers.emplace_back(L*.6, rho);
              layers.emplace_back(L*.4, rho);
              PsEigLayers = oca->P(anti*from, anti*to, EsEig, layers);

              // Layers call disturbs these
              oca->SetL(L);
              oca->SetRho(rho);
            }

            for(unsigned int i = 0; i < Ps.size(); ++i){
              if(Ps[i] != PsVec[i]){
                std::cout << "Standard and vector probabilities differ at E = " << Es[i] << ": " << Ps[i] << " vs " << PsVec[i] << std::endl;
              }

              if(PsVec[i] != PsEig[i]){
                std::cout << "Vector and eigen probabilities somehow manage to differ at E = " << Es[i] << ": " << PsVec[i] << " vs " << PsEig[i] << std::endl;
              }

              if(PsEigLayers.size() != 0 &&
                 fabs(Ps[i]-PsEigLayers[i]) > .01){
                std::cout << "Standard and eigen layers probabilities differ at E = " << Es[i] << ": " << Ps[i] << " vs " << PsEigLayers[i] << std::endl;
              }
            } // end for i
          } // end for oscIdx
        } // end for to
      } // end for from
    } // end for anti
  } // end for hie

  return 0;
}
