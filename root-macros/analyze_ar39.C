

int rebin = 4;
int freec = 10000;
int freec_dl = 10100;
int freec_red = 10200;
int freec_red_dl = 10300;
int freec_green = 10400;
int freec_green_dl = 10500;
int freec_yel = 10600;
int freec_yel_dl = 10700;

bool ana_tl = true;
bool ana_dl = false;

void init_colors_blue(int nc = 20, int init = freec) {
  float f = 150/nc;
  for (int i = 0; i <= nc; ++i) new TColor(init+i, (0.+i*f)/255., (51+i*f)/255., (102+i*f)/255., "blue");
}

void init_colors_red(int nc = 20, int init = freec_red) {
  float f = 150/nc;
  for (int i = 0; i <= nc; ++i) new TColor(init+i, (102+i*f)/255., (10+i*f)/255., (0+i*f)/255., "red");
}

void init_colors_green(int nc = 20, int init = freec_green) {
  float f = 150/nc;
  for (int i = 0; i <= nc; ++i) new TColor(init+i, (0+i*f)/255., (102+i*f)/255., (102+i*f)/255., "green");
}

void init_colors_yellow(int nc = 20, int init = freec_yel) {
  float f = 150/nc;
  for (int i = 0; i <= nc; ++i) new TColor(init+i, (102+i*f)/255., (70+i*f)/255., (0+i*f)/255., "green");
}

std::vector<double> get_observables(std::string histname, std::string cycle) {

  std::vector<double> v_obs(4, -1);
  std::string simdir = "../ph2p-ar39/";

  // open file and fetch histo
  TFile fin((simdir+cycle+"/lar/sur_array_1/Ar39/pdf-lar-sur_array_1-Ar39.root").c_str());
  TH1D * th = (TH1D*)fin.Get(histname.c_str());
  if (!th) return v_obs;

  // rebin
  th->Rebin(rebin);

  // fit
  auto func = new TF1((cycle + "_fit").c_str(), "pol3", 70, 120);
  th->Fit(func, "LRQ");
  auto func_min = new TF1((cycle + "_fit_min").c_str(), "pol3", 40, 80);
  th->Fit(func_min, "LRQ");

  v_obs.at(0) = func->GetMaximumX();     auto bmax = th->FindBin(v_obs.at(0));
  v_obs.at(1) = func_min->GetMinimumX(); auto bmin = th->FindBin(v_obs.at(1));

  v_obs.at(2) = th->Integral(bmax,bmax) / th->GetBinContent(bmin);
  v_obs.at(3) = th->GetBinContent(bmax) / th->Integral(bmax-10,bmax-1);

  std::cout << cycle << " -> max at " << v_obs.at(0) << " -> min at " << v_obs.at(1);
  std::cout << " -> peak to valley ratio " << v_obs.at(2) << " -> peak to shoulder ratio " << v_obs.at(2) << std::endl;

  fin.Close();

  // return maximum
  return v_obs;
}

void analyze_ar39_tl() {

  init_colors_blue();
  init_colors_red();
  init_colors_green();
  init_colors_yellow();
  bool init = false;

  TCanvas * can = new TCanvas("can", "can", 1600, 1400); can->Divide(2,2); can->cd(1);
  std::vector<TGraph*> * v_gr = new std::vector<TGraph*>();
  std::vector<TGraph*> * v_gr_min = new std::vector<TGraph*>();
  std::vector<TGraph*> * v_gr_ratio = new std::vector<TGraph*>();
  std::vector<TGraph*> * v_gr_ps = new std::vector<TGraph*>();
  TLegend * l       = new TLegend(0.1,0.7,0.3,0.9); l->SetNColumns(3);
  TLegend * l_min   = new TLegend(0.1,0.7,0.3,0.9); l_min->SetNColumns(3);
  TLegend * l_ratio = new TLegend(0.1,0.7,0.3,0.9); l_ratio->SetNColumns(3);
  TLegend * l_ps    = new TLegend(0.1,0.7,0.3,0.9); l_ps->SetNColumns(3);

  // 2D histogram DL vs TL vs peak max
  for (int dl = 500; dl <= 2500; dl+=100) {

    // max pos
    v_gr->push_back(new TGraph(10));
    v_gr->back()->SetLineColor(freec);
    v_gr->back()->SetMarkerColor(freec++);
    v_gr->back()->SetMarkerStyle(4);
    l->AddEntry(v_gr->back(), Form("%ium",dl), "l");
    // min pos
    v_gr_min->push_back(new TGraph(10));
    v_gr_min->back()->SetLineColor(freec_green);
    v_gr_min->back()->SetMarkerColor(freec_green++);
    v_gr_min->back()->SetMarkerStyle(4);
    l_min->AddEntry(v_gr_min->back(), Form("%ium",dl), "l");
    // peak/valley
    v_gr_ratio->push_back(new TGraph(10));
    v_gr_ratio->back()->SetLineColor(freec_red);
    v_gr_ratio->back()->SetMarkerColor(freec_red++);
    v_gr_ratio->back()->SetMarkerStyle(4);
    l_ratio->AddEntry(v_gr_ratio->back(), Form("%ium",dl), "l");
    // peak/shoulder
    v_gr_ps->push_back(new TGraph(10));
    v_gr_ps->back()->SetLineColor(freec_yel);
    v_gr_ps->back()->SetMarkerColor(freec_yel++);
    v_gr_ps->back()->SetMarkerStyle(4);
    l_ps->AddEntry(v_gr_ps->back(), Form("%ium",dl), "l");

    int N = 0;

    cout << "dl" << dl << " looping tls" << endl;

    for (int tl = 0; tl <= 9; tl++) {
      string cycle = Form("nplus-dl%ium-tl%02d", dl, tl);
      std::vector<double> pm_r = get_observables("raw/M1_ch22", cycle);
      v_gr->back()->SetPoint(N, tl/10., pm_r.at(0));
      v_gr_min->back()->SetPoint(N, tl/10., pm_r.at(1));
      v_gr_ratio->back()->SetPoint(N, tl/10., pm_r.at(2));
      v_gr_ps->back()->SetPoint(N++, tl/10., pm_r.at(3));
    }
  }

  v_gr->front()->SetTitle("Peak Maximum");
  v_gr->front()->GetXaxis()->SetTitle("transition layer fraction");
  v_gr->front()->GetYaxis()->SetTitle("peak maximum [keV]");
  v_gr->front()->GetYaxis()->SetRangeUser(75,110);
  v_gr->front()->Draw("APL");
  for (auto gr : *v_gr) gr->Draw("same PL");
  l->Draw();

  can->cd(2);
  v_gr_min->front()->SetTitle("Valley Minimum");
  v_gr_min->front()->GetXaxis()->SetTitle("transition layer fraction");
  v_gr_min->front()->GetYaxis()->SetTitle("valley minimum [keV]");
  //v_gr_min->front()->GetYaxis()->SetRangeUser(75,110);
  v_gr_min->front()->Draw("APL");
  for (auto gr : *v_gr_min) gr->Draw("same PL");
  l_min->Draw();

  can->cd(3);
  v_gr_ratio->front()->SetTitle("Peak to Valley Ratio");
  v_gr_ratio->front()->GetXaxis()->SetTitle("transition layer fraction");
  v_gr_ratio->front()->GetYaxis()->SetTitle("peak/valley");
  //v_gr_ratio->front()->GetYaxis()->SetRangeUser(75,110);
  v_gr_ratio->front()->Draw("APL");
  for (auto gr : *v_gr_ratio) gr->Draw("same PL");
  l_ratio->Draw();

  can->cd(4);
  v_gr_ps->front()->SetTitle("Peak to Shoulder Ratio");
  v_gr_ps->front()->GetXaxis()->SetTitle("transition layer fraction");
  v_gr_ps->front()->GetYaxis()->SetTitle("peak/shoulder");
  //v_gr_ps->front()->GetYaxis()->SetRangeUser(75,110);
  v_gr_ps->front()->Draw("APL");
  for (auto gr : *v_gr_ps) gr->Draw("same PL");
  l_ps->Draw();

  cout << "done" << endl;
}

void analyze_ar39_dl() {

  init_colors_blue(10, freec_dl);
  init_colors_red(10, freec_red_dl);
  init_colors_green(10, freec_green_dl);
  init_colors_yellow(10, freec_yel_dl);
  bool init = false;

  TCanvas * can = new TCanvas("can_dl", "can_dl", 1600, 1400); can->Divide(2,2); can->cd(1);
  std::vector<TGraph*> * v_gr = new std::vector<TGraph*>();
  std::vector<TGraph*> * v_gr_min = new std::vector<TGraph*>();
  std::vector<TGraph*> * v_gr_ratio = new std::vector<TGraph*>();
  std::vector<TGraph*> * v_gr_ps = new std::vector<TGraph*>();
  TLegend * l       = new TLegend(0.1,0.7,0.3,0.9); l->SetNColumns(2);
  TLegend * l_min   = new TLegend(0.1,0.7,0.3,0.9); l_min->SetNColumns(2);
  TLegend * l_ratio = new TLegend(0.7,0.7,0.9,0.9); l_ratio->SetNColumns(2);
  TLegend * l_ps    = new TLegend(0.7,0.7,0.9,0.9); l_ps->SetNColumns(2);

  for (int tl = 0; tl <= 9; tl++) {

    v_gr->push_back(new TGraph());
    v_gr->back()->SetLineColor(freec_dl);
    v_gr->back()->SetMarkerColor(freec_dl++);
    v_gr->back()->SetMarkerStyle(4);
    l->AddEntry(v_gr->back(), Form("%.1f",tl/10.), "l");

    v_gr_min->push_back(new TGraph());
    v_gr_min->back()->SetLineColor(freec_green_dl);
    v_gr_min->back()->SetMarkerColor(freec_green_dl++);
    v_gr_min->back()->SetMarkerStyle(4);
    l_min->AddEntry(v_gr_min->back(), Form("%.1f",tl/10.), "l");

    v_gr_ratio->push_back(new TGraph());
    v_gr_ratio->back()->SetLineColor(freec_red_dl);
    v_gr_ratio->back()->SetMarkerColor(freec_red_dl++);
    v_gr_ratio->back()->SetMarkerStyle(4);
    l_ratio->AddEntry(v_gr_ratio->back(), Form("%.1f",tl/10.), "l");

    v_gr_ps->push_back(new TGraph());
    v_gr_ps->back()->SetLineColor(freec_yel_dl);
    v_gr_ps->back()->SetMarkerColor(freec_yel_dl++);
    v_gr_ps->back()->SetMarkerStyle(4);
    l_ps->AddEntry(v_gr_ps->back(), Form("%.1f",tl/10.), "l");

    int N = 0;

    cout << "tl" << tl << " looping dls" << endl;

    for (int dl = 500; dl <= 2500; dl+=100) {
      string cycle = Form("nplus-dl%ium-tl%02d", dl, tl);
      std::vector<double> pm_r = get_observables("raw/M1_ch22", cycle);
      v_gr->back()->SetPoint(N, dl, pm_r.at(0));
      v_gr_min->back()->SetPoint(N, dl, pm_r.at(1));
      v_gr_ratio->back()->SetPoint(N, dl, pm_r.at(2));
      v_gr_ps->back()->SetPoint(N++, dl, pm_r.at(3));
    }
  }

  v_gr->front()->SetTitle("Peak Maximum");
  v_gr->front()->GetXaxis()->SetTitle("FCCD [um]");
  v_gr->front()->GetYaxis()->SetTitle("peak maximum [keV]");
  v_gr->front()->GetYaxis()->SetRangeUser(75,110);
  v_gr->front()->Draw("APL");
  for (auto gr : *v_gr) gr->Draw("same PL");
  l->Draw();

  can->cd(2);
  v_gr_min->front()->SetTitle("Valley Minimum");
  v_gr_min->front()->GetXaxis()->SetTitle("FCCD [um]");
  v_gr_min->front()->GetYaxis()->SetTitle("valley minimum [keV]");
  v_gr_min->front()->GetYaxis()->SetRangeUser(30,70);
  v_gr_min->front()->Draw("APL");
  for (auto gr : *v_gr_min) gr->Draw("same PL");
  l_min->Draw();

  can->cd(3);
  v_gr_ratio->front()->SetTitle("Peak to Valley Ratio");
  v_gr_ratio->front()->GetXaxis()->SetTitle("FCCD [um]");
  v_gr_ratio->front()->GetYaxis()->SetTitle("peak/valley");
  v_gr_ratio->front()->GetYaxis()->SetRangeUser(0,7);
  v_gr_ratio->front()->Draw("APL");
  for (auto gr : *v_gr_ratio) gr->Draw("same PL");
  l_ratio->Draw();

  can->cd(4);
  v_gr_ps->front()->SetTitle("Peak to Shoulder Ratio");
  v_gr_ps->front()->GetXaxis()->SetTitle("FCCD [um]");
  v_gr_ps->front()->GetYaxis()->SetTitle("peak/shoulder");
//  v_gr_ps->front()->GetYaxis()->SetRangeUser(0,7);
  v_gr_ps->front()->Draw("APL");
  for (auto gr : *v_gr_ps) gr->Draw("same PL");
  l_ps->Draw();

  cout << "done" << endl;
}

void analyze_ar39() {

  if (ana_tl) analyze_ar39_tl();
  if (ana_dl) analyze_ar39_dl();
}
