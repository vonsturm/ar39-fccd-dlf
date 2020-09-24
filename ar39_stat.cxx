// ---------------------------------------
//
// Author : K.v.Sturm
// Date   : 24.09.2020
//
// Description : Toy-MC sample Ar39 distribution and do a statistical test against possible models
//
// Usage : ./ar39_stat <options>
//
// ---------------------------------------

// C/c++
#include <cstdio>
#include <iostream>
#include <execution> // for parallel policies
#include <algorithm>
#include <type_traits>

// root cern
#include "TH1D.h"
#include "TFile.h"
#include "TRandom2.h"

// json
#include "json.hpp"

using namespace std;
using namespace nlohmann;

// dead layer model
struct dlm { int fccd; double dlf; };
// energy range
struct range { double emin; double emax; };

// fetch argument option
template<typename T>
void fetch_option(const vector<string> & args, string opt, T & var);

void usage() {
  cout << "\nChi2 test statistic sampling with Ar39 GERDA simulation\n";
  cout << "USAGE : ./ar39-stat <options>\n\n";
  cout << "OPTIONS :\n\n";
  cout << "\t\t --json       : master config file to change the models to test against\n";
  cout << "\t\t --channel -c : channel\n";
  cout << "\t\t --fccd       : fccd value in um\n";
  cout << "\t\t --dlf        : dlf value as fraction\n";
  cout << "\t\t --stat -ns   : statistics to be sampled in each toy experiment\n";
  cout << "\t\t --help -h    : print this help text\n";
  cout << "\t\t -v           : more output\n";
}

// sample from arbirary histogram distribution
// and check the chi2 distribution
TH1D * sample(TH1D * hMC, int N);

string get_filename(dlm model);

int main(int argc, char* argv[]) {

  // default input parameters
  bool verbose = false;
  bool help = false;
  
  string conf = "";
  
  int channel = 0;
  int fccd = 1000;
  double dlf = 0.5;
  int stat = 10000;
  
  int rebin = 1;
  double Emin = 45, Emax = 150;

  vector<dlm> models;
  vector<range> ranges;

  // fetch arguments
  vector<string> args(argc);
  for (int i = 0; i < argc; ++i) args.at(i) = argv[i];

  fetch_option(args, "--help", help);
  fetch_option(args, "-h",     help);

  if (help) {usage(); return 1;}

  fetch_option(args, "-v",     verbose);
  fetch_option(args, "--json", conf   );

  if (conf) {
    ifstream f_conf(conf);
    json j_conf; f_conf >> j_conf; 

    verbose = j_conf["verbose"]        .get<bool>();
    channel = j_conf["data"]["channel"].get<int>();
    fccd    = j_conf["data"]["fccd"]   .get<int>();
    dlf     = j_conf["data"]["dlf"]    .get<double>();
    stat    = j_conf["data"]["stat"]   .get<int>();
    
    for (auto f : j_conf["models"]["fccd"]) {
      for (auto d : j_conf["models"]["dlf"]) {
        models.emplace_back(f.get<int>(), d.get<double>());
      }
    }
    
    int emin_min   = j_conf["emin"]["min"]  .get<int>();
    int emin_max   = j_conf["emin"]["max"]  .get<int>();
    int emin_delta = j_conf["emin"]["delta"].get<int>();
    int emax_min   = j_conf["emax"]["min"]  .get<int>();
    int emax_max   = j_conf["emax"]["max"]  .get<int>();
    int emax_delta = j_conf["emax"]["delta"].get<int>();
    
    for (int i = emin_min; i <= emin_max; i+=emin_delta) {
      for (int j = emax_min; i <= emax_max; i+=emax_delta) {
        ranges.emplace_back((double)i, (double)j);
      }
    }
  }

  fetch_option(args, "--channel", channel);
  fetch_option(args, "-c",        channel);
  fetch_option(args, "--fccd",    fccd   );
  fetch_option(args, "--dlf",     dlf    );
  fetch_option(args, "--stat",    stat   );
  fetch_option(args, "-ns",       stat   );

  // now do stuff



  // get MC distribution
  TFile * fin = new TFile("ph2p-ar39/nplus-fccd1350um-dlf07/lar/sur_array_1/Ar39/pdf-lar-sur_array_1-Ar39.root");
  TH1D * M1_ch1 = (TH1D*) fin->Get("raw/M1_ch1");
  double mmin = M1_ch1->GetMinimum();
  double mmax = M1_ch1->GetMaximum();

  // scale to full spectrum
  TFile fdata = new TFile();

  int data_stat_ch1 = 36100; // in 45-200keV
  int N = (int) (data_stat_ch1 * M1_ch1->Integral() /
                 M1_ch1->Integral(M1_ch1->FindBin(Emin),M1_ch1->FindBin(Emax)));

  TH1D * h_dist = new TH1D("dist","dist",500,0,500);

  for (int i = 0; i < 1000; i++) {
    TH1D * s_M1_ch1 = sample(M1_ch1, N);
    M1_ch1->GetXaxis()->SetRangeUser(Emin,Emax);
    s_M1_ch1->GetXaxis()->SetRangeUser(Emin,Emax);
    double chi2 = s_M1_ch1->Chi2Test(M1_ch1, "UW CHI2");
    if (i%100 == 0) cout << i << " " << chi2 << endl;
    h_dist->Fill(chi2);
    delete s_M1_ch1; s_M1_ch1 = NULL;
    M1_ch1->GetXaxis()->SetRangeUser(mmin,mmax);
  }

  h_dist->Scale(1./h_dist->Integral());
  h_dist->SetMinimum(0);
  h_dist->SetMaximum(1);
  h_dist->Draw("hist");
  TH1D * h_cumul = (TH1D*) h_dist->GetCumulative();
  h_cumul->Draw("same hist");

  int bin68, bin90, bin95;
  h_cumul->GetBinWithContent(0.68, bin68, 0, 0, 1);
  h_cumul->GetBinWithContent(0.90, bin90, 0, 0, 1);
  h_cumul->GetBinWithContent(0.95, bin95, 0, 0, 1);

  cout << "p(0.68) : " << bin68 << endl;
  cout << "p(0.90) : "  << bin90 << endl;
  cout << "p(0.95) : " << bin95 << endl;

  fin->Close();
  fdata->Close();

  return 0;
}

// sample histogram from Emin to Emax N times
TH1D * sample(TH1D * hMC, int N) {

  // random number generator
  TRandom2 r;
  r.SetSeed(0);

  // get cumulative distribution function
  TH1D * hcumul = (TH1D*) hMC->GetCumulative();

  double norm = hMC->Integral();
  hcumul->Scale(1./norm);

  TH1D * exp = (TH1D*) hMC->Clone("exp"); exp->Reset();

  // sample using inverse of cumulative distribution function
  // to map uniform distribution to histogram
  double a_r[N];
  r.RndmArray(N, &a_r[0]);
  vector<double> v_r(&a_r[0], &a_r[0] + N);
  vector<double> v_binx(v_r.size());

  auto tstat = [&](double const& r) {
    int idx = &r - &v_r[0];
    int binx = 0;
    hcumul->GetBinWithContent(r,binx,0,0,1);
    v_binx.at(idx) = binx;
  };

  for_each(execution::par, v_r.begin(), v_r.end(), tstat);

  for (auto b : v_binx) exp->Fill(b);

  return exp;
}

template<typename T>
void fetch_option(const vector<string> & args, string opt, T & var) {
  auto result = find(args.begin(), args.end(), opt);
  if (result != args.end()) {
         if (is_same<T, int>   ::value) var = stoi(*(result+1));
    else if (is_same<T, float> ::value) var = stof(*(result+1));
    else if (is_same<T, double>::value) var = stod(*(result+1));
    else if (is_same<T, string>::value) var = *(result+1);
    else if (is_same<T, bool>  ::value) var = true;
  }
}

string get_filename(dlm model) {
  
  string name = "ph2p-ar39/nplus-fccd";
  name += to_string(model.fccd);
  name += "um-dlf";
  name += Form("%03d", model.dlf()*100); 
  name += "/lar/sur_array_1/Ar39/pdf-lar-sur_array_1-Ar39.root"
  
  return name;
}
