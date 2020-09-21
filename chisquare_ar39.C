

#include "./json.hpp"
using json = nlohmann::json;

int rebin = 1;

double emin = 45;
double emax = 150;

//bool pvalue = true;
bool pvalue = false;
bool draw_contour = true;

pair<double,double> get_chi2(std::string histname, std::string cycle, TH1D * data);
void write_to_json(vector<string> detnames, vector<double> best_DL, vector<double> best_fDL, string ofname);
double GetMean(TGraph * g);
double GetRMS(TGraph * g);

void chisquare_ar39() {

  gErrorIgnoreLevel = kBreak;

  int nchan = 41;

  // open file and fetch histo
//  TFile fin("nplus-dl1500um-tl06/lar/sur_array_1/Ar39/pdf-lar-sur_array_1-Ar39.root");
  TFile fin("gerda-data-bkgmodel-phaseIIplus-v07.00.root");
  vector<TH1D*> * v_data = new vector<TH1D*>(nchan);
  vector<string> v_detnames(nchan);
  vector<double> v_minpos_DL(nchan);
  vector<double> v_minpos_fDL(nchan);
  vector<double> v_min_chi2(nchan, (pvalue ? 0. : 10000.));
  vector<double> v_min_integral(nchan, 0.);
  vector<double> v_chan(nchan);
  for (int c = 0; c < nchan; c++) {
    v_chan.at(c) = c;
    v_data->at(c) = (TH1D*)fin.Get(Form("raw/M1_ch%i",c));
    if (!v_data->at(c)) continue;
    v_data->at(c)->Rebin(rebin);
    v_data->at(c)->GetXaxis()->SetRangeUser(emin,emax);
  }

  std::cout << " -> energy window " << emin << " : " << emax << endl;

  // dump all the shit to this file
  TFile * outfile = new TFile("chi2test_out.root","RECREATE");

  for (int c = 0; c < nchan; c++) {

    if (c==6 || c==8 || c==36) continue;

    TCanvas * can = new TCanvas(Form("can_ch%i",c), Form("can_ch%i",c), 1000, 700);
    TGraph2D * gr = new TGraph2D();
    gr->SetName(Form("gr_ch%i",c));
    if (pvalue) gr->SetTitle(Form("%s ch%i; FCCD [#mum]; DLF; p-value",v_data->at(c)->GetTitle(),c));
    else        gr->SetTitle(Form("%s ch%i; FCCD [#mum]; DLF; #Delta#chi^{2}",v_data->at(c)->GetTitle(),c));

    string dn = v_data->at(c)->GetTitle();
    int start = dn.find("(")+1;
    int end   = dn.find(")")-start;
    dn = dn.substr(start,end);
    v_detnames.at(c) = dn;

    int N = 0;

    for (int tl = 0; tl <= 9; tl++) {
      for (int dl = 500; dl <= 2500; dl+=100) {
        string cycle = Form("nplus-dl%ium-tl%02d", dl, tl);
        pair<double,double> chi2_int = get_chi2(Form("raw/M1_ch%i",c), cycle, v_data->at(c));
        double chi2     = chi2_int.first;
        double integral = chi2_int.second;
        if (chi2>0) {
          gr->SetPoint(N++, dl, tl*0.1, chi2);
          if ( ((chi2 < v_min_chi2.at(c)) && !pvalue) || // chi2 min
               ((chi2 > v_min_chi2.at(c)) &&  pvalue) ) {  // pvalue max
            v_min_chi2    .at(c) = chi2;
            v_minpos_DL   .at(c) = dl;
            v_minpos_fDL  .at(c) = tl*0.1;
            v_min_integral.at(c) = integral;
          }
        }
      }
    }

    if (!pvalue) {
      double offset = gr->GetZmin();

      for (int p = 0; p < gr->GetN(); p++) {
        double x, y, z;
        gr->GetPoint(p, x, y, z);
        gr->SetPoint(p, x, y, z-offset);
	  }
    }
    else {
      gr->SetMinimum(0);
      gr->SetMaximum(1);
    }

    gStyle->SetPalette(kPastel);
//    gStyle->SetPalette(kInvertedDarkBodyRadiator);
//    gStyle->SetPalette(kCool);
    gr->Draw("surf3");
    gr->Draw("same P0");
    if (!pvalue) gPad->SetLogz();
    else gr->GetZaxis()->SetRangeUser(0,1);

    TList* lc = gr->GetContourList(pvalue ? 0.1 : 4.6);

	if (lc) {
	if (!(lc->First() == NULL)) {
	  TGraph * gr_lc = (TGraph*)lc->First();
      double * cont_x = gr_lc->GetX();
      double * cont_y = gr_lc->GetY();

      TGraph2D * cont = new TGraph2D(gr_lc->GetN());
      double zz = pvalue ? 0 : 1;
      cont->SetLineColor(kRed);

      for (int i = 0; i < gr_lc->GetN(); i++) {
        cont->SetPoint(i, cont_x[i], cont_y[i], zz);
      }

      if (draw_contour) cont->Draw("same line");
    }}

    outfile->cd(); can->Write();

    if (pvalue) std::cout << "ch " << c << " -> max p-value " << gr->GetZmax() << endl;
    else        std::cout << "ch " << c << " -> min Chi2 " << gr->GetZmin() << endl;
    std::cout << "\t -> best DL " << v_minpos_DL.at(c);
    std::cout << " -> best fDL " << v_minpos_fDL.at(c) << endl;
    std::cout << "\t -> Ar39 activity " << v_min_integral.at(c) / 3.e11 * 37032160. * 1082.65 << " " << v_min_integral.at(c) << " Bq/kg" << endl;
  }

  cout << "DONE channels loop" << endl;

  TCanvas * can2 = new TCanvas("best","best",1800,700); can2->Divide(2,1);

  TGraph * gr_DL  = new TGraph(nchan, &v_chan[0], &v_minpos_DL[0]);
  gr_DL->SetTitle("best FCCD values; channel; FCCD [#mum]");
  gr_DL->SetMarkerStyle(23); gr_DL->SetMarkerSize(1); gr_DL->SetMarkerColor(kRed+1);
  can2->cd(1); gr_DL->Draw("ap");
  TGraph * gr_fDL = new TGraph(nchan, &v_chan[0], &v_minpos_fDL[0]);
  gr_fDL->SetTitle("best DL fraction values; channel; f_{DL} [#mum]");
  gr_fDL->SetMarkerStyle(22); gr_fDL->SetMarkerSize(1); gr_fDL->SetMarkerColor(kAzure+1);
  can2->cd(2); gr_fDL->Draw("ap");

  cout << "2. canvas plotted" << endl;

  TCanvas * can3 = new TCanvas("Ar39","Ar39",1800,700); can3->Divide(2,1);

  TGraph * gr_act = new TGraph(nchan, &v_chan[0], &v_min_integral[0]);
  gr_act->SetTitle("Ar39 activity at best fit; channel; activity [Bq/kg]");
  gr_act->SetMarkerStyle(22); gr_act->SetMarkerSize(1); gr_act->SetMarkerColor(kMagenta+1);
  gr_act->GetYaxis()->SetRangeUser(1,1.6);
  can3->cd(1); gr_act->Draw("ap");

  double mean_act = GetMean(gr_act);
  double rms_act = GetRMS(gr_act);

  TLine * l    = new TLine(0,mean_act,nchan,mean_act);
  TLine * l_up = new TLine(0,mean_act+rms_act,nchan,mean_act+rms_act);
  TLine * l_lo = new TLine(0,mean_act-rms_act,nchan,mean_act-rms_act);

  l->SetLineColor(kMagenta);
  l_up->SetLineColor(kMagenta); l_up->SetLineStyle(10);
  l_lo->SetLineColor(kMagenta); l_lo->SetLineStyle(10);

  l->Draw(); l_up->Draw(); l_lo->Draw();

  can3->cd(2);
  TH1D * h_act = new TH1D ("act_projY","act_projY",100,1.1,1.6);
  double * y = gr_act->GetY();
  for (int i = 0; i < gr_act->GetN(); i++) h_act->Fill(y[i]);
  h_act->SetTitle("^{39}Ar activity; activity [Bq/kg]; #");
  h_act->Draw("hist");

  can3->Modified();
  can3->Update();

  cout << "3. canvas plotted" << endl;

  // write values to json file
  write_to_json(v_detnames, v_minpos_DL, v_minpos_fDL, "chi2test_best_values.json");

  cout << "dumped everything to json" << endl;

  outfile->cd(); can2->Write(); can3->Write();
//  outfile->Close();

  cout << "done" << endl;
}

pair<double,double> get_chi2(std::string histname, std::string cycle, TH1D * data) {

  double chi2 = -1., integral = -1.;

  // open file and fetch histo
  TFile fin(("ph2p-ar39/"+cycle+"/lar/sur_array_1/Ar39/pdf-lar-sur_array_1-Ar39.root").c_str());
  TH1D * th = (TH1D*)fin.Get(histname.c_str());
  if (!th) return pair<double,double>(chi2,integral);

  // rebin
  th->Rebin(rebin);
  th->GetXaxis()->SetRangeUser(emin,emax);

  // calculate chi2
  // https://root.cern.ch/doc/master/classTH1.html#a6c281eebc0c0a848e7a0d620425090a5
  if (pvalue) chi2 = data->Chi2Test(th, "UW");
  else        chi2 = data->Chi2Test(th, "UW CHI2");
  integral = data->Integral(th->FindBin(emin), th->FindBin(emax))
             / th->Integral(th->FindBin(emin), th->FindBin(emax));

  // scale to Ar39 activity in Bq/kg
  // I(data)/I(MC) * P / LT / M_LAr
  integral *= 3.e11 / 37032160. / 1082.65;

  fin.Close();

  // return maximum
  return pair<double,double>(chi2,integral);
}

void write_to_json(vector<string> detnames, vector<double> best_DL, vector<double> best_fDL, string ofname) {

  auto j = json::object();
  j["chi2test"] = json::object();

  for (int c = 0; c < best_DL.size(); c++) {
    string s = detnames.at(c);
    j["chi2test"][s] = json::object();
    j["chi2test"][s]["fccd-nm"] = best_DL.at(c);
    j["chi2test"][s]["dlf"] = best_fDL.at(c);
  }

  // detectors without values
  j["chi2test"]["GD02D"] = json::object();
  j["chi2test"]["GD02D"]["fccd-nm"] = -1;
  j["chi2test"]["GD02D"]["dlf"] = -1;
  j["chi2test"]["ANG5"] = json::object();
  j["chi2test"]["ANG5"]["fccd-nm"] = -1;
  j["chi2test"]["ANG5"]["dlf"] = -1;
  j["chi2test"]["IC48B"] = json::object();
  j["chi2test"]["IC48B"]["fccd-nm"] = -1;
  j["chi2test"]["IC48B"]["dlf"] = -1;


  std::ofstream fout(ofname);
  fout << j.dump(4);
}

double GetMean(TGraph * g) {

  double mean = 0.;

  int n = g->GetN();
  int c = 0;
  double * y = g->GetY();

  for (int i=0; i<n; i++) {
    if (y[i] > 0.) {
      mean += y[i];
      c++;
    }
  }

  return mean/c;
}

double GetRMS(TGraph * g) {

  double mean = GetMean(g);
  double RMS = 0.;

  int n = g->GetN();
  int c = 0;
  double * y = g->GetY();

  for (int i=0; i<n; i++) {
    if (y[i] > 0) {
      RMS += y[i]*y[i];
      c++;
    }
  }

  RMS = TMath::Abs( RMS/c - mean*mean );

  return sqrt(RMS);
}
