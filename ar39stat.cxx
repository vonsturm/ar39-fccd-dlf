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
#include "TTree.h"
#include "TRandom2.h"

// json
#include "json.hpp"

// progress
#include "progressbar/progressbar.hpp"

//own
#include "args_reader.cxx"
#include "arbitrary_sampler.cxx"

using namespace nlohmann;

// dead layer model
struct dlm_t {
  int ID;
  int fccd;
  double dlf;
  TH1D * hist;
  std::vector<std::vector<double>> chi2;

  dlm_t(int ID, int fccd, double dlf, std::vector<std::vector<double>> chi2) 
    : ID(ID), fccd(fccd), dlf(dlf), hist(nullptr), chi2(chi2) {}
};

// energy range
struct range_t {
  int ID; double emin; double emax;
  range_t(int ID, double emin, double emax) 
    : ID(ID), emin(emin), emax(emax) {}
};

void usage() {
  std::cout << "\nChi2 test statistic sampling with Ar39 GERDA simulation\n\n";
  std::cout << "USAGE : ./ar39-stat <options>\n\n";
  std::cout << "OPTIONS :\n";
  std::cout << "\t--json       : master config file to change the models to test against\n";
  std::cout << "\t-c --channel : channel\n";
  std::cout << "\t--fccd       : fccd value in um\n";
  std::cout << "\t--dlf        : dlf value as fraction\n";
  std::cout << "\t-s --stat    : statistics to be sampled in each toy experiment\n";
  std::cout << "\t-h --help    : print this help text\n";
  std::cout << "\t--emin       : minimum energy for chi2 test\n";
  std::cout << "\t--emax       : maximum energy for chi2 test\n";
  std::cout << "\t-b           : rebin\n";
  std::cout << "\t-v           : more output\n\n";
}

std::string get_filename(dlm_t model);
std::string get_filename(int fccd, double dlf);
std::string get_ofilename(int channel, int fccd, double dlf, range_t r, std::string dir = "");
std::string get_treename(int channel, range_t r);

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
  int stat = 10000;

  int rebin = 1;
  double Emin = 45, Emax = 150;
  
  int toys = 100;

  std::vector<dlm_t> models;
  std::vector<range_t> ranges;
  
  std::string dir = "";

  // -------------------------------------------------------------------
  // fetch arguments
  // -------------------------------------------------------------------
  std::vector<std::string> args(argc);
  for (int i = 0; i < argc; ++i) args.at(i) = argv[i];

  fetch_arg(args, "--help", help);
  fetch_arg(args, "-h",     help);

  if (help) {usage(); return 1;}

  fetch_arg(args, "-v",     verbose);
  bool found = fetch_arg(args, "--json", conf);

  if (found) {
    std::ifstream f_conf(conf);
    json j_conf; f_conf >> j_conf;

    verbose = j_conf.value("verbose",verbose);
    if (verbose) std::cout << "---Found master conf---\n";

    rebin   = j_conf.value("rebin",rebin);
    toys    = j_conf.value("number-of-toys",toys);
    dir     = j_conf.value("output-dir",dir);
    if (j_conf.contains("data")) {
      channel = j_conf["data"].value("channel",channel);
      fccd    = j_conf["data"].value("fccd",fccd);
      dlf     = j_conf["data"].value("dlf",dlf);
      stat    = j_conf["data"].value("stat",stat);
    }
    if (j_conf.contains("ranges")) {
      int emin_min   = j_conf["ranges"]["emin"]["min"]  .get<int>();
      int emin_max   = j_conf["ranges"]["emin"]["max"]  .get<int>();
      int emin_delta = j_conf["ranges"]["emin"]["delta"].get<int>();
      int emax_min   = j_conf["ranges"]["emax"]["min"]  .get<int>();
      int emax_max   = j_conf["ranges"]["emax"]["max"]  .get<int>();
      int emax_delta = j_conf["ranges"]["emax"]["delta"].get<int>();
      int ID = 0;

      for (int i = emin_min; i <= emin_max; i+=emin_delta) {
        for (int j = emax_min; j <= emax_max; j+=emax_delta) {
          ranges.emplace_back(ID++, (double)i, (double)j);
        }
      }
    }
    else ranges.emplace_back(0,Emin,Emax);

    if (j_conf.contains("models")) {
      int ID = 0;
      for (auto f : j_conf["models"]["fccd"]) {
        for (auto d : j_conf["models"]["dlf"]) {
          models.emplace_back(ID++, fccd+f.get<int>(), dlf+d.get<double>(),
            std::vector<std::vector<double>>(ranges.size(), std::vector<double>(toys))); // chi
        }
      }
    }
  }

  fetch_arg(args, "--channel", channel);
  fetch_arg(args, "-c",        channel);
  fetch_arg(args, "--fccd",    fccd   );
  fetch_arg(args, "--dlf",     dlf    );
  fetch_arg(args, "-s",        stat   );
  fetch_arg(args, "--stat",    stat   );
  fetch_arg(args, "--toys",    toys   );
  fetch_arg(args, "-t",        toys   );
  fetch_arg(args, "-b",        rebin  );
  fetch_arg(args, "-o",        dir    );
  found = fetch_arg(args, "--emin", Emin) or
          fetch_arg(args, "--emax", Emax);

  if (ranges.size()<=0 or found) {
    ranges.clear();
    ranges.emplace_back(0,Emin,Emax);
  }
  else {
    Emin = ranges.front().emin;
    Emax = ranges.back().emax;
  }

  // -------------------------------------------------------------------
  // check input parameters
  // -------------------------------------------------------------------

  if (channel<0 or channel>40) {std::cout << "Invalid channel number "           << channel << std::endl; return 1;}
  if (fccd<450 or fccd>3000)   {std::cout << "FCCD allowed range 450-3000: "     << fccd    << std::endl; return 1;}
  if (dlf<0. or dlf>1.)        {std::cout << "DLF allowed range 0.0-1.0 "        << dlf     << std::endl; return 1;}
  if (stat<100)                {std::cout << "Low statistics : "                 << stat    << std::endl; return 1;}
  if (toys<10)                 {std::cout << "Low number of toy experiments : "  << toys    << std::endl; return 1;}
  if (rebin<1 or rebin>10)     {std::cout << "Rebin allowed range 1-10 keV : "   << rebin << " keV" << std::endl; return 1;}
  if (Emin<40 or Emin>100)     {std::cout << "Emin allowed range 40-100 keV : "  << Emin  << " keV" << std::endl; return 1;}
  if (Emax<100 or Emax>200)    {std::cout << "Emax allowed range 100-200 keV : " << Emax  << " keV" << std::endl; return 1;}

  // -------------------------------------------------------------------
  // print input parameters
  // -------------------------------------------------------------------
  
  if (verbose) {
    std::cout << "Channel : " << channel << std::endl;
    std::cout << "FCCD    : " << fccd    << std::endl;
    std::cout << "DLF     : " << dlf     << std::endl;
    std::cout << "stat    : " << stat    << std::endl;
    std::cout << "toys    : " << toys    << std::endl;
    std::cout << "binning : " << rebin << " keV" << std::endl;
    std::cout << "emin    : " << Emin  << " keV" << std::endl;
    std::cout << "emax    : " << Emax  << " keV" << std::endl;
    std::cout << "output dir : " << dir << "/" << std::endl;
  }

  // -------------------------------------------------------------------
  // now do the thing
  // -------------------------------------------------------------------
  TFile fin(get_filename(fccd, dlf).c_str());
  if (fin.IsZombie() or !fin.IsOpen()) {
    std::cout << "Could not open file: " << get_filename(fccd, dlf).c_str() << std::endl;
    return 1;
  }
  TH1D * M1_data = (TH1D*) fin.Get(Form("raw/M1_ch%i", channel));
  fin.Close();
  if (!M1_data) {
    std::cout << "Histogram not found: " << Form("raw/M1_ch%i", channel) << std::endl;
    return 1;
  }
  M1_data->Rebin(rebin);

  // load model histograms
  for (auto && m : models) {
    TFile fm(get_filename(m).c_str());
    m.hist = (TH1D*) fm.Get(Form("raw/M1_ch%i", channel));
    m.hist->SetName(Form("model_fccd%d_dlf%03d",m.fccd,(int)(dlf*100)));
    fm.Close();
  }

  progressbar bar(toys);
  bar.set_done_char("-");

  for (int i = 0; i < toys; i++) {
    bar.update();
    TH1D * M1_toy = sample_histo(M1_data, stat);

    for (auto && m : models) {
      for (auto r : ranges) {
        m.hist->GetXaxis()->SetRangeUser(r.emin,r.emax);
        M1_toy->GetXaxis()->SetRangeUser(r.emin,r.emax);
        m.chi2.at(r.ID).at(i) = M1_toy->Chi2Test(m.hist, "UW CHI2");
      }
    }

    delete M1_toy; M1_toy = nullptr;
  }
  std::cout << "\n";

  // model histograms are not needed anymore at this point
  for (auto && m : models) {
    delete m.hist; m.hist = nullptr;
  }
  delete M1_data; M1_data = nullptr;

  // dump everything to file
  system(Form("mkdir -p %s",dir.c_str()));
  for (auto r : ranges) {
    TFile of(get_ofilename(channel,fccd,dlf,r,dir).c_str(),"RECREATE");
    //TTree tree(get_treename(channel,r).c_str(), get_treename(channel,r).c_str());
    TTree tree("tree", get_treename(channel,r).c_str());
    std::vector<double> v_chi2(models.size());
    for (auto m : models)
      tree.Branch(Form("chi2_%i_%03d",m.fccd,(int)(m.dlf*100)), &v_chi2.at(m.ID));

    for (int i=0; i<toys; i++) {
      for (auto m : models){
        v_chi2.at(m.ID) = m.chi2.at(r.ID).at(i);
      }
      tree.Fill();
    }

    tree.Write();
    of.Close();
  }

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

std::string get_ofilename(int channel, int fccd, double dlf, range_t r, std::string dir) {
  std::string ofname = "";
  if (dir!="") ofname = dir + "/";
  ofname += "ar39stat_ch"+std::to_string(channel)+"_fccd"+std::to_string(fccd)+"um_dlf";
  ofname += Form("%03d_range",(int)(dlf*100));
  ofname += std::to_string((int)r.emin)+"-"+std::to_string((int)r.emax)+".root";

  return ofname;
}

std::string get_treename(int channel, range_t r) {
  std::string treename = "tree_ch"+std::to_string(channel)+
    "_range"+std::to_string((int)r.emin)+"_"+std::to_string((int)r.emax);

  return treename;
}
