/*
 *
 *
 *
 *
 *
 */

// C/C++
#include <iostream>
#include <vector>
#include <numeric>

// CERN root
#include "TH1.h"

// gsg
#include "gram_savitzky_golay/gram_savitzky_golay.h"

namespace smoothing {

  TH1D *  apply_mwa(TH1D * h, uint16_t nsamples, double start) {
    TH1D * h_filtered = dynamic_cast<TH1D*>( h->Clone(Form("%s_filtered_mwa",h->GetName())) );
    int begin = h->FindBin(start);
    int nb    = h->GetNbinsX();
    uint16_t & ws  = nsamples;
    std::vector<double> d;
    for (int j = begin; j < begin+ws; j++) d.push_back(h->GetBinContent(j));

    for (int i = begin; i <= nb-ws; i++) {
      double sum = accumulate(d.begin(),d.end(),0.);
      h_filtered->SetBinContent(i+nsamples,sum/ws);
      d.erase(d.begin());
      d.push_back(h->GetBinContent(i+ws-1));
    }
    return h_filtered;
  }

  TH1D *  apply_gsg(TH1D * h, uint16_t nsamples, uint16_t npoly, double start) {
    gram_sg::SavitzkyGolayFilter filter(nsamples, -nsamples, npoly, 0);
    TH1D * h_filtered = dynamic_cast<TH1D*>( h->Clone(Form("%s_filtered_gsg",h->GetName())) );
    int nb = h->GetNbinsX();
    int ws = 2*nsamples+1;
    int begin = h->FindBin(start);
    std::vector<double> d;
    for (int j = begin; j < begin+ws; j++) d.push_back(h->GetBinContent(j));

    for (int i = begin; i <= nb-ws; i++) {
      h_filtered->SetBinContent(i,filter.filter(d));
      d.erase(d.begin());
      d.push_back(h->GetBinContent(i+ws-1));
    }
    return h_filtered;
  }

  TH1D * apply_filter(TH1D * h, uint16_t filter, uint16_t nsamples, uint16_t npoly, double start) {
    switch (filter) {
      case 0 : return h;
      case 1 : return apply_mwa(h, nsamples, start);
      case 2 : return apply_gsg(h, nsamples, npoly, start);
      default : std::cout << "WARNING: Unknown filter algorithm, NO filter applied\n"; return h;
    }
  }
}