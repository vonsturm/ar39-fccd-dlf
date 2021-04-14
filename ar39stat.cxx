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
#include <iomanip>
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
#include "args-reader/args_reader.hpp"
#include "arbitrary_sampler.hpp"

using namespace nlohmann;

// dead layer model
struct dlm_t {
  int ID;
  int fccd;
  double dlf;
  TH1D * hist;
  std::vector<double> chi2;

  dlm_t(int ID, int fccd, double dlf, std::vector<double> chi2)
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
  std::cout << "USAGE : ./ar39-stat --json config.json <options>\n\n";
  std::cout << "OPTIONS :\n";
  std::cout << "\t-h --help          : print this help text\n";
  std::cout << "\t--json <opt>       : master config file [conf.json]\n";
  std::cout << "\t--emin <opt>       : minimum energy for chi2 test [45-100]\n";
  std::cout << "\t--emax <opt>       : maximum energy for chi2 test [100-200]\n";
  std::cout << "\t-r <opt>           : rebin\n";
  std::cout << "\t-o <opt>           : output directory\n";
  std::cout << "\t--test <opt>       : test statistics (0 = delta Chi2 => default, 1 = Chi2Test, 2 = KolmogorovTest)\n";
  std::cout << "\t-v                 : more output\n\n";
}

std::string get_filename(dlm_t model);
std::string get_filename(int fccd, double dlf);
std::string get_ofilename(bool toys, int channel, int rebin, range_t r, std::string dir);
std::string get_treename(int channel, range_t r);

int main(int argc, char* argv[]) {

  TH1::AddDirectory(false);

  // -------------------------------------------------------------------
  // default input parameters
  // -------------------------------------------------------------------
  bool verbose = false;
  bool help    = false;
  bool toys    = false;

  std::string conf   = "";
  std::string infile = "";
  std::string outdir = "";
  uint16_t channel = 0;

  uint16_t rebin = 1;
  uint16_t teststat = 0;
  range_t fit_range(10000, 45, 160);

  std::vector<dlm_t> models;

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

    j_conf.value("verbose",verbose);
    if (verbose) std::cout << "---Found master conf---\n";

    rebin  = j_conf.value("rebin",rebin);
    outdir = j_conf.value("output-dir",outdir);
    if (j_conf.contains("data")) {
      infile = j_conf["data"].value("file",infile);
      toys   = j_conf["data"].value("toys",toys);
    }
    if (j_conf.contains("fit-range")) {
      fit_range = range_t( 
        0,
        j_conf["fit-range"].value("emin",45),
        j_conf["fit-range"].value("emax",160)
      );
    }

    if (j_conf.contains("models")) {
      int ID = 0;
      int fccd_start = j_conf["models"]["fccd"]["start"].get<int>(); 
      int fccd_stop  = j_conf["models"]["fccd"]["stop"].get<int>(); 
      int fccd_step  = j_conf["models"]["fccd"]["step"].get<int>(); 
      double dlf_start = j_conf["models"]["dlf"]["start"].get<double>(); 
      double dlf_stop  = j_conf["models"]["dlf"]["stop"].get<double>(); 
      double dlf_step  = j_conf["models"]["dlf"]["step"].get<double>(); 
      for (int f = fccd_start; f <= fccd_stop; f += fccd_step ) {
        for (double d = dlf_start; d <= dlf_stop; d += dlf_step ) {
          models.emplace_back(ID++, f, d, std::vector<double>()); // chi
        }
      }
    }
  }

  fetch_arg(args, "--test", teststat);

  fetch_arg(args, "--emin", fit_range.emin);
  fetch_arg(args, "--emax", fit_range.emax);

  fetch_arg(args, "-r", rebin  );
  fetch_arg(args, "-o", outdir );

  // -------------------------------------------------------------------
  // check input parameters
  // -------------------------------------------------------------------
  if (fit_range.emin > fit_range.emax) { std::cout << "Error: Fit range minimum greater than maxium. Aborting.\n";        exit(EXIT_FAILURE); }
  if (fit_range.emin < 45)             { std::cout << "Fit range minimum too small " << fit_range.emin << " Aborting.\n"; exit(EXIT_FAILURE); }
  if (fit_range.emax > 2000)           { std::cout << "Fit range maximum too large " << fit_range.emax << " Aborting.\n"; exit(EXIT_FAILURE); }
  if (rebin<1 or rebin>10)             { std::cout << "Rebin allowed range 1-10 keV : " << rebin << " keV. Aborting.\n";  exit(EXIT_FAILURE); }

  // -------------------------------------------------------------------
  // print input parameters
  // -------------------------------------------------------------------

  if (verbose) {
    if (toys) std::cout << "Processing: gerda-factory toy MCs - " << infile << std::endl;
    else      std::cout << "Processing: Real Data - " << infile << std::endl;
    std::cout << "binning : "    << rebin          << " keV\n";
    std::cout << "range   : "    << fit_range.emin << " - " << fit_range.emax  << " keV\n";
    std::cout << "test statistics: "; 
    switch (teststat) {
      case 0  : std::cout << "delta Chi2\n";     break;
      case 1  : std::cout << "Chi2Test\n";       break;
      case 2  : std::cout << "KolmogorovTest\n"; break;
      default : std::cout << "Test statistics not implemented using default: delta Chi2\n"; break;
    }
  }

  // -------------------------------------------------------------------
  // now do the thing
  // -------------------------------------------------------------------
  TFile fin(infile.c_str());
  if (fin.IsZombie() or !fin.IsOpen()) {
    std::cout << "Could not open file: " << infile << std::endl;
    exit(EXIT_FAILURE);
  }

  // fill data/toys vector
  std::vector<TH1D*> v_data(100000);
  if (!toys) {
    v_data[0] = dynamic_cast<TH1D*>( fin.Get(Form("raw/M1_%i",channel)) );
    v_data.resize(1);
  }
  else {
    int hct = 0;
    for (auto && keyAsObj : *fin.GetListOfKeys()){
      auto key = (TKey*) keyAsObj;
      //std::cout << key->GetName() << " " << key->GetClassName() << std::endl;
      if (std::string(key->GetClassName()) == "TH1D") {
        v_data[hct++] = dynamic_cast<TH1D*>( fin.Get(key->GetName()) );
      }
    }
    v_data.resize(hct);
  }
  fin.Close();

  // load model histograms 
  for (auto && m : models) {
    TFile fm(get_filename(m).c_str());
    m.hist = (TH1D*) fm.Get(Form("raw/M1_ch%i", channel));
    m.hist->Rebin(rebin);
    //add here the ratio between the integrals
    m.hist->SetName(Form("model_fccd%d_dlf%03d",m.fccd,(int)round(m.dlf*100)));
    m.chi2.resize(v_data.size());
    fm.Close();
  }

  progressbar bar(v_data.size());
  bar.set_done_char("-");

  if (toys) std::cout << "Number of toys: " << v_data.size() << std::endl;

  int i = 0;
  for (auto && data : v_data) {
    // rebin
    data->Rebin(rebin);
    bar.update();

    double min_chi2 = 10000.;

    for (auto && m : models) {
      m.hist->GetXaxis()->SetRangeUser(fit_range.emin,fit_range.emax);
      data  ->GetXaxis()->SetRangeUser(fit_range.emin,fit_range.emax);

      if (data->Integral(fit_range.emin,fit_range.emax) < 10) {
        std::cout << "\nEmpty Hist: " << data->GetName() << " " << data->Integral() << std::endl;
        exit(EXIT_FAILURE);
      }
      if (m.hist->Integral(fit_range.emin,fit_range.emax) < 10) {
        std::cout << "\nEmpty Hist: " << m.hist->GetName() << " " << m.hist->Integral() << std::endl;
        exit(EXIT_FAILURE);
      }

      switch (teststat) {
        case 1  : m.chi2.at(i) = data->Chi2Test(m.hist, "UW CHI2"); break;
        case 2  : m.chi2.at(i) = data->KolmogorovTest(m.hist);      break;
        default : m.chi2.at(i) = data->Chi2Test(m.hist, "UW CHI2");
                  min_chi2 = m.chi2.at(i) < min_chi2 ? m.chi2.at(i) : min_chi2;
                  break;
      }
    }

    // implement delta chi2 here
    if (teststat == 0 or teststat > 2) {
      for (auto && m : models) m.chi2.at(i) -= min_chi2;
    }

    i++;
  }
  std::cout << "\n";

  // model histograms are not needed anymore at this point
  for (auto && m : models) { delete m.hist; m.hist = nullptr; }
  for (auto && d : v_data) { delete d; d = nullptr; }

  // dump everything to file
  if (outdir!="") system(Form("mkdir -p %s",outdir.c_str()));

  TFile of(get_ofilename(toys,channel,rebin,fit_range,outdir).c_str(),"RECREATE");
  TTree tree("statTree", "statTree");
  std::vector<double> v_chi2(models.size());
  int best_fccd = 0; double best_dlf = 0;
  for (auto m : models)
    tree.Branch(Form("chi2_%i_%03d",m.fccd,(int)(round(m.dlf*100))), &v_chi2.at(m.ID));
  tree.Branch("best_fccd", &best_fccd);
  tree.Branch("best_dlf",  &best_dlf);

  for (size_t i=0; i<v_data.size(); i++) {
    for (auto m : models) {
      v_chi2.at(m.ID) = m.chi2.at(i);
      if (m.chi2.at(i) == 0) { best_fccd = m.fccd; best_dlf = m.dlf; } 
    }
    tree.Fill();
  }

  tree.Write();
  of.Close();

  return 0;
}

std::string get_filename(int fccd, double dlf) {
  std::string name = "ph2p-ar39/nplus-fccd";
  name += std::to_string(fccd);
  name += "um-dlf";
  name += Form("%03d", (int)(round(dlf*100))); 
  name += "/lar/sur_array_1/Ar39/pdf-lar-sur_array_1-Ar39.root";

  return name;
}

std::string get_filename(dlm_t model) {
  return get_filename(model.fccd, model.dlf);
}

std::string get_ofilename(bool toys, int channel, int rebin, range_t r, std::string dir) {
  std::string ofname = "";
  if (dir!="") ofname = dir + "/";
  if (toys)    ofname += "toys_";
  ofname += "ar39stat_ch" + std::to_string(channel) + "_rebin" + std::to_string(rebin) + "_";
  ofname += std::to_string((int)r.emin)+"-"+std::to_string((int)r.emax)+".root";

  return ofname;
}

std::string get_treename(int channel, range_t r) {
  std::string treename = "tree_ch"+std::to_string(channel)+
    "_range"+std::to_string((int)r.emin)+"_"+std::to_string((int)r.emax);

  return treename;
}
