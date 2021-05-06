

void proj(int ch) {

  gStyle->SetPalette(kSunset);
  gStyle->SetLabelFont(43,"XYZ");
  gStyle->SetLabelSize(13,"XYZ");
  gStyle->SetNdivisions(504, "XYZ");

  TFile * f = new TFile(Form("../output/data/bestfit_ar39stat_mipol_ch%i_rebin10_45-160.root",ch));
  TParameter<int>    * best_fccd  = (TParameter<int>*)    f->Get("best_fccd");
  TParameter<double> * best_dlf   = (TParameter<double>*) f->Get("best_dlf");

  TH2D * ts = dynamic_cast<TH2D*>( f->Get("ts_profile") );
  ts->SetLabelFont(43,"XYZ");
  ts->SetLabelSize(13,"XYZ");

  ts->SetMinimum(0); ts->SetMaximum(10);

  int bx,by,bz;
  ts->GetBinXYZ( ts->GetMinimumBin(), bx, by, bz );

  TH1D * ts_px = ts->ProjectionX("_px",by,by);
  TH1D * ts_py = ts->ProjectionY("_py",bx,bx);

  TCanvas * c = new TCanvas(Form("c_ch%i",ch),Form("c_ch%i",ch),1000,800);
  TPad * p1 = new TPad("p1","p1",0,0,0.8,0.8);
  TPad * p2 = new TPad("p2","p2",0.8,0,1,0.8);
  TPad * p3 = new TPad("p3","p3",0,0.8,0.8,1);

  p1->SetTopMargin(0); p1->SetRightMargin(0);
  p2->SetTopMargin(0); p2->SetLeftMargin(0);
  p3->SetRightMargin(0.001); p3->SetBottomMargin(0);

  ts_px->GetYaxis()->SetRangeUser(0,10);
  ts_py->GetYaxis()->SetRangeUser(0,10);

  int col = kGray;// kAzure-9;
  ts_px->SetFillStyle(1001);
  ts_py->SetFillStyle(1001);
  ts_px->SetFillColor(col);
  ts_py->SetFillColor(col);
  ts_px->SetLineColor(col);
  ts_py->SetLineColor(col);

  ts_px->GetYaxis()->SetRangeUser(0.001,10);
  ts_py->GetYaxis()->SetRangeUser(0.001,10);
  ts_px->GetYaxis()->SetTitle("#Delta#chi^{2}");
  ts_py->GetYaxis()->SetTitle("#Delta#chi^{2}");
  ts_px->GetYaxis()->SetTitleFont(43);
  ts_px->GetYaxis()->SetTitleSize(17);
  ts_px->GetYaxis()->SetTitleOffset(1.5);
  ts_py->GetYaxis()->SetTitleFont(43);
  ts_py->GetYaxis()->SetTitleSize(17);
  ts_py->GetYaxis()->SetTitleOffset(1.5);

  ts->GetXaxis()->SetRangeUser(650,2390);
  ts->GetYaxis()->SetRangeUser(0,0.99);

  TLegend * l = new TLegend(0.7,0.8,0.9,0.9);
  l->AddEntry(ts,Form("FCCD %i",best_fccd->GetVal()),"");
  l->AddEntry(ts,Form("DLF %.2f",best_dlf->GetVal()),"");

  c->cd(); p1->Draw(); p2->Draw(); p3->Draw();
  p1->cd(); ts->Draw("CONT4"); l->Draw();
  p2->cd(); ts_py->Draw("hbar");
  p3->cd(); ts_px->Draw("hist");


/*
  c->SaveAs(Form("projection_ch%i.png",ch));
  delete c;
  delete f;
*/
}
