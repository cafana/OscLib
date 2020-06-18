#include <cmath>
#include <vector>

#include "OscLib/OscCalc.h"
#include "OscLib/OscCalcGeneral.h"
#include "OscLib/OscCalcPMNS.h"
#include "OscLib/OscCalcPMNSOpt.h"
#include "OscLib/OscCalcDMP.h"
#include "OscLib/OscCalcPMNS_CPT.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLegend.h"
#include "TMath.h"

#include <iostream>
#include <fenv.h>

#ifndef DARWINBUILD
#include "Utilities/StanVar.h"
#endif
#include "Utilities/StanUtils.h"

int main()
{

#ifndef DARWINBUILD
  feenableexcept(FE_INVALID); // Spot any infs or nans early
#else
  std::cerr << "WARNING: OscLib/test/testOsc.cc was built on OS X where feenableexcept is unavailable" << std::endl;
#endif

  TCanvas* canvs[4];
  canvs[0] = new TCanvas("canv_norm");
  canvs[1] = new TCanvas("canv_norm_anti");
  canvs[2] = new TCanvas("canv_inv");
  canvs[3] = new TCanvas("canv_inv_anti");

  for(int n = 0; n < 4; ++n) canvs[n]->Divide(3, 3);

  TLegend* leg = new TLegend(.1, .1, .9, .9);

  osc::OscCalc osc1;
  osc::OscCalcGeneral osc2;
  osc::OscCalcPMNS osc3;
  osc::OscCalcPMNSOpt osc4;
  osc::OscCalcPMNS_CPT osc5;
  osc::_OscCalcPMNS<stan::math::var> osc6;
  osc::_OscCalcPMNSOpt<stan::math::var> osc7;
  osc::_OscCalcDMP<stan::math::var> osc8;

  std::vector<osc::IOscCalcAdjustable*> oscs {&osc1, &osc2, &osc3, &osc4, &osc5};
  std::vector<osc::_IOscCalcAdjustable<stan::math::var>*> oscStan {&osc6, &osc7, &osc8};
  std::vector<TString> names {"Approx", "General", "PMNS", "PMNSOpt", "PMNS_CPT", "PMNS_Stan", "PMNSOpt_Stan", "PMNSOpt_DMP_Stan"};
  std::vector<int> colors {kBlue, kRed, kGreen+2, kMagenta, kBlack, kCyan, kOrange, kViolet+1};

  const double L = 800;
  const double rho = 3;
  const double dmsq21 = 7.6e-5;
  const double dmsq32 = 2.35e-3;
  const double th12 = 1;
  const double th13 = .15;
  const double th23 = TMath::Pi()/4-.05;
  const double delta = 1.5*TMath::Pi();

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
      if (names[n] == "PMNS_CPT"){
          osc5.SetDmsq21Bar(dmsq21);
          osc5.SetDmsq32Bar(hie*dmsq32);
          osc5.SetTh12Bar(th12);
          osc5.SetTh13Bar(th13);
          osc5.SetTh23Bar(th23);
          osc5.SetdCPBar(delta);
      }
    }
    for(auto osc : oscStan){
      osc->SetL(L);
      osc->SetRho(rho);
      osc->SetDmsq21(dmsq21);
      osc->SetDmsq32(hie*dmsq32);
      osc->SetTh12(th12);
      osc->SetTh13(th13);
      osc->SetTh23(th23);
      osc->SetdCP(delta);

    }

    auto numCalcs = oscs.size() + oscStan.size();

    TString nus[3] = {"e", "#mu", "#tau"};
    TString antinus[3] = {"#bar{e}", "#bar{#mu}", "#bar{#tau}"};

    for(int anti = -1; anti <= +1; anti += 2){
      for(int from = 12; from <= 16; from += 2){
        for(int to = 12; to <= 16; to += 2){
          std::cout << "Calculating for " << anti*from << " to " << anti*to << " hierachy " << hie << std::endl;
          TString title = TString::Format("%s to %s",
                                          (anti < 0 ? antinus[from/2-6] : nus[from/2-6]).Data(),
                                          (anti < 0 ? antinus[to  /2-6] : nus[to  /2-6]).Data());

          TPad* canv = canvs[(1-anti)/2+(1-hie)];
          const int padIdx = (from/2-6)+3*(to/2-6)+1;
          canv->cd(padIdx);
          std::vector<TGraph*> gs;
          for(std::size_t n = 0; n < numCalcs; ++n) gs.push_back(new TGraph);
          for(double E = .1; E < 5; E *= 1.01){
            std::vector<double> Ps;
            std::size_t n = 0;
            for(; n < oscs.size(); ++n){
              Ps.push_back(oscs[n]->P(anti*from, anti*to, E));
              gs[n]->SetPoint(gs[n]->GetN(), E, Ps[n]);
            }
            for(auto osc : oscStan){
              Ps.push_back(util::GetValAs<double>(osc->P(anti * from, anti * to, E)));
              gs[n]->SetPoint(gs[n]->GetN(), E, Ps[n]);
              ++n;
            }
            for(std::size_t i = 0; i < Ps.size(); ++i){
              for(std::size_t j = i+1; j < Ps.size(); ++j){
                if(fabs(Ps[i]-Ps[j]) > .01 && E > 0.25){
                  std::cerr << "!!! Probabilities for " << title
                            << " differ at " << E << " GeV. "
                            << Ps[i] << " vs " << Ps[j]
                            << " between calculators " << i << " and " << j << std::endl;
                }
              } // end for j
            } // end for i
          } // end for E

          gs[0]->SetTitle(title);
          gs[0]->GetXaxis()->SetTitle("E (GeV)");

          for(std::size_t n = 0; n < numCalcs; ++n){
            gs[n]->SetLineColor(colors[n]);
            if(n) gs[n]->SetLineStyle(7);
            gs[n]->Draw(n ? "l same" : "al");

            if(hie == -1 && anti == -1 && from == 12 && to == 12){
              leg->AddEntry(gs[n], names[n], "l");
            }
          }

          gPad->SetLogx();
        } // end for to
      } // end for from
    } // end for anti
  } // end for hie

  TFile* f = new TFile("testOsc.root", "RECREATE");
  for(int n = 0; n < 4; ++n) canvs[n]->Write();
  TCanvas* cleg = new TCanvas;
  leg->Draw();
  cleg->Write("leg");
  f->Close();

  std::cout << "See output file testOsc.root for oscillation curves" << std::endl;

  return 0;
}
