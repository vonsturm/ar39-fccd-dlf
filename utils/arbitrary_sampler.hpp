
#include <vector>
#include <algorithm>
#include <execution>

#include "TH1D.h"
#include "TRandom2.h"

// sample arbitrary histogram
TH1D * sample_histo(TH1D * hist, int stat) {

  // random number generator
  TRandom2 r; r.SetSeed(0);

  // get cumulative distribution function
  TH1D * hcumul = (TH1D*) hist->GetCumulative();

  double norm = hist->Integral();
  hcumul->Scale(1./norm);

  TH1D * exp = (TH1D*) hist->Clone("exp"); exp->Reset();

  // sample using inverse of cumulative distribution function
  // to map uniform distribution to histogram
  double a_r[stat];
  r.RndmArray(stat, &a_r[0]);
  std::vector<double> v_r(&a_r[0], &a_r[0] + stat);
  std::vector<double> v_binx(v_r.size());

  auto tstat = [&](double const& r) {
    int idx = &r - &v_r[0];
    int binx = 0;
    hcumul->GetBinWithContent(r,binx,0,0,1);
    v_binx.at(idx) = binx;
  };

  std::for_each(std::execution::par, v_r.begin(), v_r.end(), tstat);

  for (auto b : v_binx) exp->Fill(b);

  return exp;
}
