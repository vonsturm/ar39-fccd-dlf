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

// savitzky golay filter
#include "gram_savitzky_golay/gram_savitzky_golay.h"

//own
#include "args-reader/args_reader.hpp"
#include "arbitrary_sampler.hpp"
#include "filter.hpp"

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
  std::cout << "USAGE : ./ar39-stat --json config.json <options>\n\n";
  std::cout << "OPTIONS :\n";
  std::cout << "\t-h --help           : print this help text\n";
  std::cout << "\t--json <opt>        : master config file [conf.json]\n";
  std::cout << "\t-c --channel <opt>  : channel [0-40]\n";
  std::cout << "\t--fccd <opt>        : fccd value in um [450-3000:50] must be available\n";
  std::cout << "\t--dlf <opt>         : dlf value as fraction [0.00-1.00:0.05] must be available\n";
  std::cout << "\t-s --stat <opt>     : statistics to be sampled in each toy experiment\n";
  std::cout << "\t--range <opt>       : minimum and maximum energy for chi2 test [45-100] and [100-200]\n";
  std::cout << "\t-t --toys <opt>     : number of toy experiments\n";
  std::cout << "\t-b <opt>            : rebin\n";
  std::cout << "\t-o <opt>            : output directory\n";
  std::cout << "\t--datastat <opt>    : data statistics json file [datastat.json]\n";
  std::cout << "\t--test <opt>        : test statistics (0 = delta Chi2 => default, 1 = Chi2Test, 2 = KolmogorovTest)\n";
  std::cout << "\t--filter-gsg <2opt> : apply smoothing (savitzky-golay) to data histogram from which toys are sampled\n";
  std::cout << "                      : parameters are <number of points> and <order of polynomial> e.g. 22 3\n";
  std::cout << "\t--filter-mwa <opt>  : apply smoothing (moving window average) to data histogram from which toys are sampled\n";
  std::cout << "                      : parameters are <window width>\n";
  std::cout << "\t--filter-start <opt>: apply filter only starting from sample with energy [keV]\n";
  std::cout << "\t-v                  : more output\n\n";

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
  std::string datastat = "";

  int channel = 0;
  int fccd = 1000;
  double dlf = 0.5;
  int stat = 20000;
  int data_stat = 0;
  range_t data_range(10000, 50, 130);

  int rebin = 1;
  double Emin = 50, Emax = 130;

  int toys = 100;

  uint16_t test = 0; // 0 = delta Chi2 => default, 1 = Chi2Test, 2 = KolmogorovTest

  std::vector<dlm_t> models;
  std::vector<range_t> ranges;

  std::string dir = "";

  uint16_t filter = 0;         // 0 none, 1 moving window, 2 savitzky golay
  uint16_t nsamples = 0, npoly = 0; // in case 1 only first case 2 both
  double start = 0;               // start energy of filter algorithm

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
      if (j_conf["data"].contains("filter")) {
        std::string s_filter = j_conf["data"]["filter"].value("type","gsg");
        filter = s_filter == "mwa" ? 1 : 2;
        nsamples  = j_conf["data"]["filter"].value("width",3);
        npoly     = j_conf["data"]["filter"].value("pol-order",1);
        start     = j_conf["data"]["filter"].value("start",0);
      }
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
            std::vector<std::vector<double>>(ranges.size(), std::vector<double>())); // chi
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
  fetch_arg(args, "--test",    test   );

  // filter settings
  std::vector<uint16_t> v_filter_conf(2);
  fetch_arg(args, "--filter-start", start);
  // mwa
  found  = fetch_arg(args, "--filter-mwa", nsamples);
  filter = found ? 1 : filter;
  // gsg
  found  = fetch_arg(args, "--filter-gsg", nsamples, npoly);
  filter = found ? 2 : filter;
  
  found = fetch_arg(args, "--range", Emin, Emax);

  if (ranges.size()<=0 or found) {
    ranges.clear();
    ranges.emplace_back(0,Emin,Emax);
  }
  else {
    Emin = ranges.front().emin;
    Emax = ranges.back().emax;
  }

  // resize toys vectors
  for (auto & m : models) {
    for (auto & v : m.chi2) {
      v.resize(toys);
    }
  }

  bool found_dstat = fetch_arg(args, "--datastat", datastat);
  if (found_dstat) {
    std::ifstream f_datastat(datastat);
    json j_datastat; f_datastat >> j_datastat;
    if (j_datastat.contains("data_stat")) {
      std::string key = Form("M1_ch%i",channel);
      if (!j_datastat["data_stat"].contains(key) or 
          !j_datastat["data_stat"].contains("range")) {
            std::cerr << "--datastat option given but json file is missing necessary keys" << std::endl;
            exit(EXIT_FAILURE);
      }
      if (!j_datastat["data_stat"]["range"].contains("emin") or
          !j_datastat["data_stat"]["range"].contains("emax")) {
            std::cerr << "--datastat option given but json file is missing necessary keys" << std::endl;
            exit(EXIT_FAILURE);
      }
      data_stat = j_datastat["data_stat"][key];
      data_range.emin = j_datastat["range"]["emin"];
      data_range.emax = j_datastat["range"]["emax"];
    }
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
  if (filter == 2 and nsamples <= npoly) {
    std::cout << "WARNING: savitzky-golay filter number of samples (" << nsamples
              << ") too low for polynomial order (" << npoly << ")" << std::endl;
    return 1;
  }

  // -------------------------------------------------------------------
  // print input parameters
  // -------------------------------------------------------------------

  if (verbose) {
    std::cout << "Channel : "        << channel << std::endl;
    std::cout << "FCCD    : "        << fccd    << std::endl;
    std::cout << "DLF     : "        << dlf     << std::endl;
    std::cout << "toys    : "        << toys    << std::endl;
    std::cout << "binning : "        << rebin   << " keV" << std::endl;
    std::cout << "emin    : "        << Emin    << " keV" << std::endl;
    std::cout << "emax    : "        << Emax    << " keV" << std::endl;
    std::cout << "output dir : "     << dir     << "/" << std::endl;
    std::cout << "test statistics: "; 
    switch (test) {
      case 0  : std::cout << "delta Chi2\n";     break;
      case 1  : std::cout << "Chi2Test\n";       break;
      case 2  : std::cout << "KolmogorovTest\n"; break;
      default : std::cout << "Test statistics not implemented using default: delta Chi2\n"; break;
    }
    if (!found_dstat) std::cout << "stat    : " << stat << std::endl;
    std::cout << "filter  : "; 
    switch (filter) {
      case 0  : std::cout << "none\n"; break;
      case 1  : std::cout << "moving window average [width = " << nsamples << "]\n"; break;
      case 2  : std::cout << "savitzky-golay [width = " << nsamples << ", poly-order = " << npoly << "]\n"; break;
      default : std::cout << "unknown\n"; break;
    }
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

  if (found_dstat) {
    stat = (double)data_stat * M1_data->Integral() /
                               M1_data->Integral(M1_data->FindFixBin(data_range.emin),
                                                 M1_data->FindFixBin(data_range.emax));
    if (verbose)
      std::cout << "stat    : " << stat << " (matches data) "<< std::endl;
  }
  // rebin
  M1_data->Rebin(rebin);

  // apply filter algorithm
  TH1D * M1_smooth = smoothing::apply_filter(M1_data, filter, nsamples, npoly, start);

  // load model histograms
  for (auto && m : models) {
    TFile fm(get_filename(m).c_str());
    m.hist = (TH1D*) fm.Get(Form("raw/M1_ch%i", channel));
    m.hist->Rebin(rebin);
    //add here the ratio between the integrals
    m.hist->SetName(Form("model_fccd%d_dlf%03d",m.fccd,(int)(dlf*100)));
    fm.Close();
  }

  progressbar bar(toys);
  bar.set_done_char("-");

  for (int i = 0; i < toys; i++) {
    bar.update();
    TH1D * M1_toy = sample_histo(M1_smooth, stat);

    std::vector<double> min_chi2(ranges.size(),10000.);

    for (auto && m : models) {
      for (auto r : ranges) {
        m.hist->GetXaxis()->SetRangeUser(r.emin,r.emax);
        M1_toy->GetXaxis()->SetRangeUser(r.emin,r.emax);
        switch (test) {
          case 1  : m.chi2.at(r.ID).at(i) = M1_toy->Chi2Test(m.hist, "UW CHI2"); break;
          case 2  : m.chi2.at(r.ID).at(i) = M1_toy->KolmogorovTest(m.hist);      break;
          default : m.chi2.at(r.ID).at(i) = M1_toy->Chi2Test(m.hist, "UW CHI2");
                    min_chi2.at(r.ID) = m.chi2.at(r.ID).at(i) < min_chi2.at(r.ID) ? m.chi2.at(r.ID).at(i) : min_chi2.at(r.ID);
                    break;
        }
      }
    }

    // apply delta chi2 offset
    if (test == 0 or test > 2) {
      for (auto && m : models) {
        for (auto r : ranges) {
          m.chi2.at(r.ID).at(i) -= min_chi2.at(r.ID);
        }
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
  if (dir!="") system(Form("mkdir -p %s",dir.c_str()));
  for (auto r : ranges) {
    TFile of(get_ofilename(channel,fccd,dlf,r,dir).c_str(),"RECREATE");
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
