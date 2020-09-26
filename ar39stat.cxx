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
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <execution> // for parallel policies
#include <algorithm>
#include <type_traits>

// root cern
#include "TH1D.h"
#include "TFile.h"
#include "TRandom2.h"

// json
#include "json.hpp"

//own
#include "args_reader.cxx"
#include "arbitrary_sampler.cxx"

using namespace nlohmann;

// dead layer model
struct dlm_t {
  int fccd;
  double dlf;
  TH1D * hist;
  std::vector<std::vector<double>> chi2;
  std::vector<double> pV68;
  std::vector<double> pV90;
  std::vector<double> pV95;
  
  dlm_t(int fccd, double dlf, std::vector<std::vector<double>> chi2,
    std::vector<double> pV68, std::vector<double> pV90, std::vector<double> pV95) 
    : fccd(fccd), dlf(dlf), hist(nullptr), chi2(chi2),
      pV68(pV68), pV90(pV90), pV95(pV95) {}
};

// energy range
struct range_t {
  int ID; double emin; double emax;
  range_t(int ID, double emin, double emax) 
    : ID(ID), emin(emin), emax(emax) {}
};

void usage() {
  std::cout << "\nChi2 test statistic sampling with Ar39 GERDA simulation\n";
  std::cout << "USAGE : ./ar39-stat <options>\n\n";
  std::cout << "OPTIONS :\n\n";
  std::cout << "\t\t --json       : master config file to change the models to test against\n";
  std::cout << "\t\t -c --channel : channel\n";
  std::cout << "\t\t --fccd       : fccd value in um\n";
  std::cout << "\t\t --dlf        : dlf value as fraction\n";
  std::cout << "\t\t -s --stat    : statistics to be sampled in each toy experiment\n";
  std::cout << "\t\t -h --help    : print this help text\n";
  std::cout << "\t\t --emin       : minimum energy for chi2 test\n";
  std::cout << "\t\t --emax       : maximum energy for chi2 test\n";
  std::cout << "\t\t -b           : rebin\n";
  std::cout << "\t\t -v           : more output\n";
}

std::string get_filename(dlm_t model);
std::string get_filename(int fccd, double dlf);
std::string get_ofilename(int channel, int fccd, double dlf, range_t r);
std::string get_treename(int channel, dlm_t m);

int main(int argc, char* argv[]) {

  TH1::AddDirectory(false);

  // -------------------------------------------------------------------
  // default input parameters
  // -------------------------------------------------------------------
  bool verbose = false;
  bool help = false;

  std::string conf = "";
  
  int channel = 0;
  int fccd = 1000;
  double dlf = 0.5;
  int stat = 30000;

  int rebin = 1;
  double Emin = 45, Emax = 150;
  
  int toys = 10000;

  std::vector<dlm_t> models;
  std::vector<range_t> ranges;

  // -------------------------------------------------------------------
  // fetch arguments
  // -------------------------------------------------------------------
  std::vector<std::string> args(argc);
  for (int i = 0; i < argc; ++i) args.at(i) = argv[i];

  fetch_arg(args, "--help", help);
  fetch_arg(args, "-h",     help);

  if (help) {usage(); return 1;}

  fetch_arg(args, "-v",     verbose);
  fetch_arg(args, "--json", conf   );

  if (conf!="") {
    std::ifstream f_conf(conf);
    json j_conf; f_conf >> j_conf; 

    verbose = j_conf["verbose"]        .get<bool>();
    channel = j_conf["data"]["channel"].get<int>();
    fccd    = j_conf["data"]["fccd"]   .get<int>();
    dlf     = j_conf["data"]["dlf"]    .get<double>();
    stat    = j_conf["data"]["stat"]   .get<int>();
    rebin   = j_conf["rebin"]          .get<int>();
    toys    = j_conf["number-of-toys"] .get<int>();

    int emin_min   = j_conf["emin"]["min"]  .get<int>();
    int emin_max   = j_conf["emin"]["max"]  .get<int>();
    int emin_delta = j_conf["emin"]["delta"].get<int>();
    int emax_min   = j_conf["emax"]["min"]  .get<int>();
    int emax_max   = j_conf["emax"]["max"]  .get<int>();
    int emax_delta = j_conf["emax"]["delta"].get<int>();
    int ID = 0;

    for (int i = emin_min; i <= emin_max; i+=emin_delta) {
      for (int j = emax_min; i <= emax_max; i+=emax_delta) {
        ranges.emplace_back(ID++, (double)i, (double)j);
      }
    }

    for (auto f : j_conf["models"]["fccd"]) {
      for (auto d : j_conf["models"]["dlf"]) {
        models.emplace_back(fccd+f.get<int>(), dlf+d.get<double>(),
          std::vector<std::vector<double>>(ranges.size(), std::vector<double>(toys)), // chi2
          std::vector<double>(ranges.size()),  // pV68
          std::vector<double>(ranges.size()),  // pV90
          std::vector<double>(ranges.size())); // pV95
      }
    }
  }

  fetch_arg(args, "--channel", channel);
  fetch_arg(args, "-c",        channel);
  fetch_arg(args, "--fccd",    fccd   );
  fetch_arg(args, "--dlf",     dlf    );
  fetch_arg(args, "--stat",    stat   );
  fetch_arg(args, "--toys",    toys   );
  fetch_arg(args, "-t",        toys   );
  fetch_arg(args, "-s",        stat   );
  fetch_arg(args, "-b",        rebin  );
  fetch_arg(args, "--emin",    Emin   );
  fetch_arg(args, "--emax",    Emax   );

  if (ranges.size()<=0) ranges.emplace_back(0,Emin,Emax);

  // -------------------------------------------------------------------
  // now do the thing
  // -------------------------------------------------------------------
  TFile * fin = new TFile(get_filename(fccd, dlf).c_str());
  TH1D * M1_data = (TH1D*) fin->Get(Form("raw/M1_ch%i", channel));

  // load model histograms
  for (auto && m : models) {
    TFile fm(get_filename(m).c_str());
    m.hist = (TH1D*) fm.Get(Form("raw/M1_ch%i", channel));
    m.hist->SetName(Form("model_fccd%d_dlf%03d",m.fccd,(int)(dlf*100)));
    fm.Close();
  }

  for (int i = 0; i < toys; i++) {
    TH1D * M1_toy = sample_histo(M1_data, stat);

    for (auto && m : models) {
      if (verbose and i%(toys/10) == 0) 
        std::cout << i << " fccd " << m.fccd << " dlf " << m.dlf << std::endl;
      for (auto r : ranges) {
        if (verbose and i%(toys/10) == 0) 
          std::cout << i << " [" << r.emin << ":" << r.emax << "]" << std::endl;
        m.hist->GetXaxis()->SetRangeUser(r.emin,r.emax);
        M1_toy->GetXaxis()->SetRangeUser(r.emin,r.emax);
        m.chi2.at(r.ID).at(i) = M1_toy->Chi2Test(m.hist, "UW CHI2");
      }
    }  

    delete M1_toy; M1_toy = nullptr;
  }
  
  // model histograms are not needed anymore at this point
  for (auto && m : models) {
    delete m.hist; m.hist = nullptr;
  }

  
  // fill all chi2 distributions in histograms
  TH1D chi2_dist("chi2_dist","chi2_dist",500,0,500);

  for (auto && m : models) {
    for (auto r : ranges) {
      chi2_dist.Reset();
      for (int i = 0; i < toys; i++) {
        chi2_dist.Fill(m.chi2.at(r.ID).at(i));
      }
      chi2_dist.Scale(1./chi2_dist.Integral());
      TH1D * chi2_cumul = (TH1D*) chi2_dist.GetCumulative();

      int bin68, bin90, bin95;
      chi2_cumul->GetBinWithContent(0.68, bin68, 0, 0, 1);
      chi2_cumul->GetBinWithContent(0.90, bin90, 0, 0, 1);
      chi2_cumul->GetBinWithContent(0.95, bin95, 0, 0, 1);

      m.pV68.at(r.ID) = bin68;
      m.pV90.at(r.ID) = bin90;
      m.pV95.at(r.ID) = bin95;
      
      delete chi2_cumul; chi2_cumul = nullptr;
    }
  }
  
/*
  // dump everything to file
  for (auto r : ranges) {
    TFile of(get_ofilename(channel,fccd,dlf,r).c_str(),"RECREATE");

    for (auto m : models){
      TTree * tree = new TTree(get_treename(channel,m).c_str());
      
      tree->Branch();
      //loop
        tree->Fill();
      tree->Write();
    }
    of.Close();
  }
*/

  return 0;
}

std::string get_filename(int fccd, double dlf) {
  std::string name = "ph2p-ar39/nplus-fccd";
  name += std::to_string(fccd);
  name += "um-dlf";
  name += Form("%03d", (int)(dlf*100)); 
  name += "/lar/sur_array_1/Ar39/pdf-lar-sur_array_1-Ar39.root";

  return name;
}

std::string get_filename(dlm_t model) {
  return get_filename(model.fccd, model.dlf);
}

std::string get_ofilename(int channel, int fccd, double dlf, range_t r) {
  std::string ofname = "ar39stat_ch"+std::to_string(channel)+"_fccd"+std::to_string(fccd)+"um_dlf";
  ofname += Form("%03d_range",(int)(dlf*100));
  ofname += std::to_string((int)r.emin)+"-"+std::to_string((int)r.emax)+".root";
  
  return ofname;
}

std::string get_treename(int channel, dlm_t m) {
  std::string treename = "tree_ch"+std::to_string(channel)+"_fccd"+std::to_string(m.fccd)+"um_dlf";
  treename += Form("%03d_range",(int)(m.dlf*100));
  
  return treename;
}
