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
#include <limits>
#include <cmath>

// root cern
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TRandom2.h"
#include "TParameter.h"
#include "TStyle.h"
#include "TMath.h"

// json
#include "json.hpp"

// progress
#include "progressbar/progressbar.hpp"

// gerda-ar39-pdf
#include "gerda-ar39-pdf/include/gerda_ar39_pdf.hpp"

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
  range_t() : ID(0), emin(0.), emax(0.) {}
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
  std::cout << "\t-c --channel <opt> : channel\n";
  std::cout << "\t--activity <opt>   : scale to activity <opt>Bq/kg\n";
  std::cout << "\t-r <opt>           : rebin\n";
  std::cout << "\t-o <opt>           : output directory\n";
  std::cout << "\t--test <opt>       : test statistics even number plain, odd number delta ts\n";
  std::cout << "\t                   : (0|1 = Chi2Test, 2|3 = Chi2Test/NDF, 4|5 = Kolmogorov, 6/7 Chi2 by-hand)\n";
  std::cout << "\t-i --interpolate   : use gerda-ar39-pdf to interpolate between discrete pdfs in fccd and dlf\n";
  std::cout << "\t-v                 : more output\n\n";
}

std::string get_filename(dlm_t model);
std::string get_filename(int fccd, double dlf);
std::string get_ofilename(bool toys, bool interpolate, int channel, int rebin, range_t r);
std::string get_treename(int channel, range_t r);
void scale_TH1D_to_integral(TH1D * h, range_t range);
double GetChi2(TH1D * h_data, TH1D * h_model, range_t range, bool use_fixed_activity=false, double activity=1.);
double GetChi2Opt(TH1D * h_data, TH1D * h_model, range_t range);
double CalcTeststatistic(uint16_t teststat, TH1D * h_data, TH1D * h_model,
                         range_t fit_range=range_t(), bool use_fixed_activity=false, double activity=1.);
double CalcPValue(uint16_t teststat, TH1D * h_data, TH1D * h_model,
                  range_t fit_range=range_t(), bool use_fixed_activity=false, double activity=1.);

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
  std::string prefix = "";
  uint16_t channel = 0;

  double activity = 1; // in Bq/kg
  uint16_t rebin = 1;
  uint16_t teststat = 7;
  bool interpolate = false;
  range_t fit_range(10000, 45, 160);

  std::vector<dlm_t> models;

  // -------------------------------------------------------------------
  // fetch arguments
  // -------------------------------------------------------------------
  std::vector<std::string> args(argc);
  for (int i = 0; i < argc; ++i) args.at(i) = argv[i];

  args_reader::fetch_arg(args, "--help", help);
  args_reader::fetch_arg(args, "-h",     help);

  if (help) {usage(); return 0;}

  args_reader::fetch_arg(args, "-v",     verbose);
  bool found = args_reader::fetch_arg(args, "--json", conf);

  if (!found) {
    std::cout << "Error: --json option is mandatory\n\n";
    usage();
    return 1;
  }

  // read master config file
  std::ifstream f_conf(conf);
  if (!f_conf.is_open()) {
    std::cout << "Error: Could not open config file: " << conf << std::endl;
    exit(EXIT_FAILURE);
  }

  json j_conf; f_conf >> j_conf;

  j_conf.value("verbose",verbose);
  if (verbose) std::cout << "---Found master conf---\n";

  rebin  = j_conf.value("rebin",rebin);
  outdir = j_conf.value("output-dir",outdir);
  prefix = j_conf.value("output-file-prefix",prefix);

  if (j_conf.contains("data")) {
    infile  = j_conf["data"].value("file",infile);
    channel = j_conf["data"].value("channel",channel);
    toys    = j_conf["data"].value("toys",toys);
  }

  if (j_conf.contains("fit-range")) {
    fit_range = range_t(
      0,
      j_conf["fit-range"].value("emin",45),
      j_conf["fit-range"].value("emax",160)
    );
  }

  int fccd_start = 650, fccd_stop = 2400, fccd_step = 100;
  double dlf_start = 0., dlf_stop = 1., dlf_step = 0.1;
  if (j_conf.contains("models")) {
    int ID = 0;
    fccd_start = j_conf["models"]["fccd"].value("start", 650);
    fccd_stop  = j_conf["models"]["fccd"].value("stop", 2400);
    fccd_step  = j_conf["models"]["fccd"].value("step",  100);
    dlf_start  = j_conf["models"]["dlf"] .value("start", 0.);
    dlf_stop   = j_conf["models"]["dlf"] .value("stop",  1.);
    dlf_step   = j_conf["models"]["dlf"] .value("step", 0.1);
    for (int f = fccd_start; f <= fccd_stop; f += fccd_step ) {
      for (double d = dlf_start; d <= dlf_stop; d += dlf_step ) {
        models.emplace_back(ID++, f, d, std::vector<double>()); // chi
      }
    }
  }

  args_reader::fetch_arg(args, "--test", teststat);

  args_reader::fetch_arg(args, "-i", interpolate);
  args_reader::fetch_arg(args, "--interpolate", interpolate);

  args_reader::fetch_arg(args, "--emin", fit_range.emin);
  args_reader::fetch_arg(args, "--emax", fit_range.emax);
  args_reader::fetch_arg(args, "-c",     channel);
  args_reader::fetch_arg(args, "--channel",     channel);

  bool use_fixed_activity =
    args_reader::fetch_arg(args, "--activity", activity);
  args_reader::fetch_arg(args, "-r", rebin );
  args_reader::fetch_arg(args, "-o", outdir);

  // -------------------------------------------------------------------
  // check input parameters
  // -------------------------------------------------------------------
  if (fit_range.emin > fit_range.emax) { std::cout << "Error: Fit range minimum greater than maxium. Aborting.\n";        exit(EXIT_FAILURE); }
  if (fit_range.emin < 45)             { std::cout << "Fit range minimum too small " << fit_range.emin << " Aborting.\n"; exit(EXIT_FAILURE); }
  if (fit_range.emax > 2000)           { std::cout << "Fit range maximum too large " << fit_range.emax << " Aborting.\n"; exit(EXIT_FAILURE); }
  if (rebin<1 or rebin>100)            { std::cout << "Rebin allowed range 1-100 : " << rebin << " Aborting.\n";  exit(EXIT_FAILURE); }

  // -------------------------------------------------------------------
  // print input parameters
  // -------------------------------------------------------------------

  if (verbose) {
    if (toys) std::cout << "Processing: gerda-factory toy MCs - " << infile << std::endl;
    else      std::cout << "Processing: Real Data - " << infile << std::endl;
    std::cout << "channel: " << channel << "\n";
    if (use_fixed_activity) std::cout << "fixed activity: " << activity << "\n";
    std::cout << "rebin: " << rebin << "\n";
    std::cout << "range: " << fit_range.emin << " - " << fit_range.emax  << " keV\n";
    std::cout << "test statistics: ";
    switch (teststat) {
      case 0  : std::cout << "Chi2Test\n";             break;
      case 1  : std::cout << "delta Chi2Test\n";       break;
      case 2  : std::cout << "Chi2Test/NDF\n";         break;
      case 3  : std::cout << "delta Chi2Test/NDF\n";   break;
      case 4  : std::cout << "KolmogorovTest\n";       break;
      case 5  : std::cout << "delta KolmogorovTest\n"; break;
      case 6  : std::cout << "Chi2 by-hand\n";         break;
      case 7  : std::cout << "delta Chi2 by-hand\n";   break;
      default : std::cout << "Test statistics not implemented using default (7): delta Chi2 by-handt\n";
                break;
    }
    if (!interpolate)  std::cout << "DO NOT ";
    std::cout << "USE gerda-ar39-pdf to interpolate between discrete pdfs\n";
  }

  // -------------------------------------------------------------------
  // now do the thing
  // -------------------------------------------------------------------
  TFile fin(infile.c_str());
  if (fin.IsZombie() or !fin.IsOpen()) {
    std::cout << "Error: Could not open data file: " << infile << std::endl;
    exit(EXIT_FAILURE);
  }

  // fill data/toys vector
  std::vector<TH1D*> v_data(100000);
  if (!toys) {
    v_data[0] = dynamic_cast<TH1D*>( fin.Get(Form("raw/M1_ch%i",channel)) );
    v_data.resize(1);
  }
  else {
    int hct = 0;
    for (auto && keyAsObj : *fin.GetListOfKeys()){
      auto key = (TKey*) keyAsObj;
      if (std::string(key->GetClassName()) == "TH1D") {
        v_data[hct] = dynamic_cast<TH1D*>( fin.Get(key->GetName()) );
        if (teststat < 6) { // make sure histogram is unweighted
          for (int b = 1; b <= v_data[hct]->GetNbinsX(); b++) {
            double cont = v_data[hct]->GetBinContent(b);
            v_data[hct]->SetBinContent(b,0);
            for (int i = 0; i < cont; i++)
              v_data[hct]->Fill(v_data[hct]->GetBinCenter(b));
          }
        }
        hct++;
      }
    }
    v_data.resize(hct);
  }
  fin.Close();

  // Checkpoint
  if (verbose) std::cout << "Checkpoint alpha" << std::endl;

  // load model histograms
  for (auto && m : models) {
    std::string mname = Form("model_fccd%d_dlf%03d",m.fccd,(int)round(m.dlf*100));
    std::string mtitle = mname + Form(";energy[keV];cts / %.1fkeV",v_data[0]->GetBinWidth(1)*rebin);
    if (!interpolate) {
      TFile fm(get_filename(m).c_str());
      if (!fm.IsOpen()) {
        m.hist = nullptr;
        continue;
      }
      m.hist = dynamic_cast<TH1D*>( fm.Get(Form("raw/M1_ch%i", channel)) );
      m.hist->SetName(mname.c_str());
      m.hist->SetTitle(mtitle.c_str());
      fm.Close();
    }
    else {
      m.hist = new TH1D(mname.c_str(), mtitle.c_str(),
        v_data[0]->GetNbinsX(),
        v_data[0]->GetBinLowEdge(1),
        v_data[0]->GetBinLowEdge(v_data[0]->GetNbinsX()+1));
      // now fill the histogram
      int nbins = m.hist->GetNbinsX();
      for (int b = 1; b <= nbins; b++) {
        double bin_center = m.hist->GetBinCenter(b);
        if (30. <= bin_center && bin_center <= 300.) {
          double cont = gerda::ar39_pdf(channel, bin_center, m.fccd/1000., m.dlf);
          m.hist->SetBinContent(b,cont);
        }
        m.hist->SetBinError(b,0.);
      }
      if (teststat < 6) m.hist->Scale(4.e10);
    }
    m.hist->Rebin(rebin);
    m.chi2.resize(v_data.size());
  }

  // Checkpoint
  if (verbose) std::cout << "Checkpoint bravo" << std::endl;

  // Write model histograms
  if (interpolate) {
    TFile ofh((outdir+"/models_"+prefix+get_ofilename(toys,interpolate,channel,rebin,fit_range)).c_str(),"RECREATE");
    TCanvas * cm_fccd = new TCanvas("cm_fccd","cm_fccd",1000,700);
    TCanvas * cm_dlf  = new TCanvas("cm_dlf","cm_dlf",1000,700);
    int i = 0, j = 0;
    for (auto && m : models) {
      if (round(m.dlf*100) == 60) {
        cm_fccd->cd();
        m.hist->SetLineColor( i%2==0 ? kRed : kBlue );
        m.hist->DrawCopy(i == 0 ? "hist" : "hist same");
        i++;
      }
      if (m.fccd == 1800) {
        cm_dlf->cd();
        m.hist->SetLineColor( j%2==0 ? kRed : kBlue );
        m.hist->DrawCopy(j == 0 ? "hist" : "hist same");
        j++;
      }
      // write grid point histograms for dlf 0.6 to file
      if (m.fccd%50 == 0 && round(m.dlf*100) == 60) { ofh.cd(); m.hist->Write(); }
    }
    cm_fccd->Write();
    cm_dlf->Write();
    ofh.Close();
  }

  // Checkpoint
  if (verbose) std::cout << "Checkpoint charlie" << std::endl;

  progressbar bar(toys ? v_data.size() : models.size());
  bar.set_done_char("-");

  if (toys) std::cout << "Number of toys: " << v_data.size() << std::endl;

  int i = 0;
  for (auto && data : v_data) {
    // rebin
    data->Rebin(rebin);

    // check binning
    double bw_data = data->GetBinWidth(1);
    int nxbins_data = data->GetNbinsX();
    if (verbose && i==0) std::cout << "binning : " << bw_data << "keV\n";

    if (toys) bar.update();

    for (auto && m : models) {

      if (!toys) bar.update();
      if (m.hist == nullptr) continue;
      // check binning
      double bw_m = m.hist->GetBinWidth(1);
      int nxbins_m = m.hist->GetNbinsX();

      if (bw_data != bw_m) {
        std::cout << "\nError: Bin width data[" << i << "] " << bw_data << "keV"
                  <<  "/ model [fccd " << m.fccd << ", dlf" << m.dlf << "] " << bw_m << "keV are different\n";
        exit(EXIT_FAILURE);
      }
      if (nxbins_data != nxbins_m) {
        std::cout << "\nError: N Bins X data[" << i << "] " << nxbins_data
                  <<  "/ model [fccd " << m.fccd << ", dlf" << m.dlf << "] " << nxbins_m << " are different\n";
        exit(EXIT_FAILURE);
      }

      m.hist->GetXaxis()->SetRangeUser(fit_range.emin,fit_range.emax);
      data  ->GetXaxis()->SetRangeUser(fit_range.emin,fit_range.emax);

      if (teststat < 6) {
        if (data->Integral(data->FindBin(fit_range.emin),data->FindBin(fit_range.emax)) < 10) {
          std::cout << "\nEmpty Hist: " << data->GetName() << " " << data->Integral() << std::endl;
          exit(EXIT_FAILURE);
        }
        if (m.hist->Integral(m.hist->FindBin(fit_range.emin),m.hist->FindBin(fit_range.emax)) < 10) {
          std::cout << "\nEmpty Hist: " << m.hist->GetName() << " " << m.hist->Integral() << std::endl;
          exit(EXIT_FAILURE);
        }
      }

      m.chi2.at(i) = CalcTeststatistic(teststat,data,m.hist,fit_range,use_fixed_activity,activity);
    }
    i++;
  }
  std::cout << "\n";

  // Checkpoint
  if (verbose) std::cout << "Checkpoint delta" << std::endl;

  system(Form("mkdir -p %s",outdir.c_str()));
  TFile of((outdir+"/"+prefix+get_ofilename(toys,interpolate,channel,rebin,fit_range)).c_str(),"RECREATE");
  TTree tree("statTree", "statTree");
  std::vector<double> v_chi2(models.size());
  int best_fccd = 0; double best_dlf = 0., best_activity = 0., gof = -1.;
  for (auto && m : models)
    tree.Branch(Form("chi2_%i_%03d",m.fccd,(int)(round(m.dlf*100))), &v_chi2.at(m.ID));
  tree.Branch("best_fccd",     &best_fccd     );
  tree.Branch("best_dlf",      &best_dlf      );
  tree.Branch("best_activity", &best_activity );
  tree.Branch("gof",           &gof           );

  auto min_llh = std::begin(v_chi2);

  for (size_t i=0; i<v_data.size(); i++) {
    // fill model vector
    for (auto && m : models) {
      if (m.hist != nullptr) v_chi2.at(m.ID) = m.chi2.at(i);
      else                   v_chi2.at(m.ID) = std::numeric_limits<double>::max();
    }
    // find minimum
    min_llh = std::min_element(std::begin(v_chi2),std::end(v_chi2));
    // delta teststat for odd IDs
    if (teststat%2 == 1) {
      std::transform(std::begin(v_chi2), std::end(v_chi2), std::begin(v_chi2), [c=*min_llh](double x) { return x -= c; } );
    }
    // minimum corrisponds to best fit
    auto best_model = models.at(min_llh-std::begin(v_chi2));
    best_fccd = best_model.fccd;
    best_dlf  = best_model.dlf;
    // calculate activity for best fit model
    // primaries / (LT[s] * Volume[cm3] * Density[g/cm3] / 1000[g/kg])
    double sim_act = 1.e12 / (37032.160 * 188400 * 1.39);
    double model_int = best_model.hist->Integral(fit_range.emin, fit_range.emax);
    double data_int  = v_data.at(i)   ->Integral(fit_range.emin, fit_range.emax);
    best_activity = sim_act * data_int/model_int;
    //gof = v_data.at(i)->Chi2Test(best_model.hist, "UW CHI2/NDF");
    //gof = v_data.at(i)->Chi2Test(best_model.hist, "UW");
    gof = CalcPValue(teststat,v_data.at(i),best_model.hist,fit_range,use_fixed_activity,activity);
    tree.Fill();
  }

  tree.Write();
  of.Close();

  // Checkpoint
  if (verbose) std::cout << "Checkpoint echo" << std::endl;

  // some useful output
  if (!toys) {
    TFile ofh((outdir+"/bestfit_"+prefix+get_ofilename(toys,interpolate,channel,rebin,fit_range)).c_str(),"RECREATE");
/*
    // Canvas with data, best model and range in DLF
    TCanvas * c1 = new TCanvas("c_dlf","c_dlf",1000,500);
    scale_TH1D_to_integral(v_data.at(0), fit_range);
    v_data.at(0)->SetLineColor(kAzure);
    v_data.at(0)->Draw("hist");
    std::vector<TH1D*> v_hm_dlf;
    for (int i = -2; i < 3; i++) {
      v_hm_dlf.push_back( models.at(min_llh-std::begin(v_chi2)+i).hist );
      if (i==0) v_hm_dlf.back()->SetLineColor(kRed);
      else      v_hm_dlf.back()->SetLineColor(kGreen+1);
    }
    for (auto & hm : v_hm_dlf) {
      scale_TH1D_to_integral(hm, fit_range);
      hm->Draw("hist same");
    }
    c1->Write("best_fit_dlf");

    // Canvas with data, best model and range in FCCD
    TCanvas * c2 = new TCanvas("c_fccd","c_fccd",1000,500);
    v_data[0]->Draw("hist");
    std::vector<TH1D*> v_hm_fccd;
    int delta_fccd = (fccd_stop-fccd_start)/fccd_step+1;
    for (int i = -2; i < 3; i++) {
      v_hm_fccd.push_back( models.at(min_llh-std::begin(v_chi2)+i*delta_fccd).hist );
      if (i==0) v_hm_fccd.back()->SetLineColor(kRed);
      else      v_hm_fccd.back()->SetLineColor(kGreen+1);
    }
    for (auto & hm : v_hm_fccd) {
      scale_TH1D_to_integral(hm, fit_range);
      hm->Draw("hist same");
    }
    c2->Write("best_fit_fccd");
*/
    TParameter<int>    p_best_fccd("best_fccd", best_fccd);
    TParameter<double> p_best_dlf ("best_dlf",  best_dlf );
    TParameter<double> p_gof      ("gof",       gof      );
    p_best_fccd.Write();
    p_best_dlf.Write();
    p_gof.Write();

    // Canvas with TS profile
    TH2D * ts_profile = new TH2D("ts_profile","ts_profile",
      (fccd_stop-fccd_start)/fccd_step+1,fccd_start-fccd_step/2,fccd_stop+fccd_step/2,
      round((dlf_stop-dlf_start)/dlf_step)+1,dlf_start-dlf_step/2,dlf_stop+dlf_step/2);
    int N = 0;
    for (auto chi2 : v_chi2) {
      ts_profile->Fill(models.at(N).fccd, models.at(N).dlf, chi2);
      N++;
    }
    ts_profile->SetTitle("ts_profile; fccd [#mum]; dlf");
    ts_profile->Write();
    ofh.Close();
  }

  // Checkpoint
  if (verbose) std::cout << "Checkpoint foxtrot" << std::endl;

  // model and data histograms are not needed anymore at this point
  for (auto && m : models) { delete m.hist; m.hist = nullptr; }
  for (auto && d : v_data) { delete d; d = nullptr; }

  // Checkpoint
  if (verbose) std::cout << "Checkpoint golf" << std::endl;

  return 0;
}

std::string get_filename(int fccd, double dlf) {
  std::string name = "ph2p-ar39/nplus-fccd";
  name += std::to_string(fccd);
  name += "um-dlf";
  name += Form("%04d", (int)(round(dlf*1000)));
  name += "/lar/sur_array_4/Ar39/pdf-lar-sur_array_4-Ar39.root";

  return name;
}

std::string get_filename(dlm_t model) {
  return get_filename(model.fccd, model.dlf);
}

std::string get_ofilename(bool toys, bool interpolate, int channel, int rebin, range_t r) {
  std::string ofname = "ar39stat_";
  if (toys)   ofname += "toys_";
  if (interpolate) ofname += "mipol_";
  ofname += "ch" + std::to_string(channel) + "_rebin" + std::to_string(rebin) + "_";
  ofname += std::to_string((int)r.emin)+"-"+std::to_string((int)r.emax)+".root";

  return ofname;
}

std::string get_treename(int channel, range_t r) {
  std::string treename = "tree_ch"+std::to_string(channel)+
    "_range"+std::to_string((int)r.emin)+"_"+std::to_string((int)r.emax);

  return treename;
}

void scale_TH1D_to_integral(TH1D * h, range_t range) {
  h->Scale(1./h->Integral(h->FindBin(range.emin), h->FindBin(range.emax)));
}

double GetChi2(TH1D * h_data, TH1D * h_model, range_t range, bool use_fixed_activity, double activity) {

  int bmin = h_data->FindBin(range.emin);
  int bmax = h_data->FindBin(range.emax);

  double int_data  = h_data->Integral(bmin,bmax);
  double int_model = h_model->Integral(bmin,bmax);

  // do not actually scale but use a correction factor
  double sim_act = 1.e12 / (37032.160 * 188400 * 1.39);
  double cf = use_fixed_activity ? activity/sim_act : int_data/int_model;

  double chi2 = 0.;

  for (int b = bmin; b <= bmax; b++) {
    double c_data  = h_data->GetBinContent(b);
    double c_model = cf * h_model->GetBinContent(b);
    chi2 += (c_model - c_data) * (c_model - c_data) / c_model;
  }

  return chi2;
}

double GetChi2Opt(TH1D * h_data, TH1D * h_model, range_t range) {

  int bmin = h_data->FindBin(range.emin);
  int bmax = h_data->FindBin(range.emax);

  double int_data  = h_data->Integral(bmin,bmax);

  // try different scalings +-1%
  double chi2opt = -1.;

  for (int i = -10; i < 11; i++) {
    double chi2 = 0.;
    double int_model = h_model->Integral(bmin,bmax);
    h_model->Scale((1.+0.001*i)*int_data/int_model);

    for (int b = bmin; b <= bmax; b++) {
      double c_data  = h_data->GetBinContent(b);
      double c_model = h_model->GetBinContent(b);
      chi2 += (c_model - c_data) * (c_model - c_data) / c_model;
    }
    chi2opt = chi2opt < 0 ? chi2 : std::min(chi2,chi2opt);
  }

  return chi2opt;
}

double CalcTeststatistic(uint16_t teststat, TH1D * h_data, TH1D * h_model,
                         range_t fit_range, bool use_fixed_activity, double activity)
{
  double ts = DBL_MAX;

  switch (teststat) {
    case 0  :                                                      // Chi2Test
    case 1  : ts = h_data->Chi2Test(h_model,"UW CHI2");     break; // Chi2Test delta
    case 2  :                                                      // Chi2Test/NDF
    case 3  : ts = h_data->Chi2Test(h_model,"UW CHI2/NDF"); break; // Chi2Test/NDF delta
    case 4  :                                                      // KolmogorovTest
    case 5  : ts = h_data->KolmogorovTest(h_model);         break; // KolmogorovTest delta
    case 6  :                                                      // Chi2 by-hand
    case 7  :                                                      // Chi2 by-hand delta
    default : ts = GetChi2(h_data,h_model,fit_range,use_fixed_activity,activity); break; // default 7
  }

  return ts;
}

double CalcPValue(uint16_t teststat, TH1D * h_data, TH1D * h_model,
                  range_t fit_range, bool use_fixed_activity, double activity)
{
  double pValue = DBL_MAX;

  switch (teststat) {
    case 0  :                                                      // Chi2Test
    case 1  : pValue = h_data->Chi2Test(h_model,"UW");      break; // Chi2Test delta
    case 2  :                                                      // Chi2Test/NDF
    case 3  : pValue = h_data->Chi2Test(h_model,"UW");      break; // Chi2Test/NDF delta
    case 4  :                                                      // KolmogorovTest
    case 5  : pValue = h_data->KolmogorovTest(h_model,"X"); break; // KolmogorovTest delta
    case 6  :                                                      // Chi2 by-hand
    case 7  :                                                      // Chi2 by-hand delta
    default : pValue =
                TMath::Prob(
                  GetChi2(h_data,h_model,fit_range,use_fixed_activity,activity),
                  h_data->FindBin(fit_range.emax) - h_data->FindBin(fit_range.emin) - 2
                );                                          break; // default 7
  }

  return pValue;
}
