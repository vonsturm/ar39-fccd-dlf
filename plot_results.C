#include "json.hpp"
using namespace nlohmann;

void plot_results() {

    std::ifstream fold("ged-transition-linear-bjorn.json");
    std::ifstream fnew("ged-transition-ar39.json");
    std::ifstream fchi2("chi2test_best_values.json");

    json jold, jnew, jchi2;
    fold >> jold;
    fnew >> jnew;
    fchi2 >> jchi2;

    jold = jold["nplus-transition"];
    jnew = jnew["nplus-transition"];
    jchi2 = jchi2["chi2test"];

    TH1I *fccd_old = new TH1I("fccd_old", "FCCD thickness;;thickness (#mum)", 41, 0, 41);
    TH1I *fccd_new = new TH1I("fccd_new", "fccd_new", 41, 0, 41);
    TH1I *fccd_chi2 = new TH1I("fccd_chi2", "fccd_chi2", 41, 0, 41);
    TH1D *tl_old = new TH1D("tl_old", "dead layer fraction;;DLF", 41, 0, 41);
    TH1D *tl_new = new TH1D("tl_new", "tl_new", 41, 0, 41);
    TH1D *tl_chi2 = new TH1D("tl_chi2", "tl_chi2", 41, 0, 41);
    int i = 1;

    for (auto o : jold.items()) {
        fccd_old->GetXaxis()->SetBinLabel(i, o.key().c_str());
        fccd_old->SetBinContent(i, jold[o.key()]["fccd-nm"].get<double>()/1E3);
        fccd_new->SetBinContent(i, jnew[o.key()]["fccd-nm"].get<double>()/1E3);
        fccd_chi2->SetBinContent(i, jchi2[o.key()]["fccd-nm"].get<double>());
        tl_old->GetXaxis()->SetBinLabel(i, o.key().c_str());
        tl_old->SetBinContent(i, jold[o.key()]["dlf"].get<double>());
        tl_new->SetBinContent(i, jnew[o.key()]["dlf"].get<double>());
        tl_chi2->SetBinContent(i, jchi2[o.key()]["dlf"].get<double>());
        i++;
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    c1->SetLeftMargin(0.12);
    c1->cd();

    fccd_old->GetXaxis()->LabelsOption("v");
    fccd_old->GetXaxis()->SetLabelSize(0.04);
    fccd_old->SetMarkerStyle(26);
    fccd_new->SetMarkerStyle(22);
    fccd_new->SetMarkerColor(kRed);
    fccd_chi2->SetMarkerStyle(4);
    fccd_chi2->SetMarkerColor(kBlue+1);

    fccd_old->Draw("p0");
    fccd_new->Draw("p0 same");
    fccd_chi2->Draw("p0 same");
    gPad->SetGridx();

    auto leg = new TLegend();
    leg->AddEntry(fccd_old, "official values", "p");
    leg->AddEntry(fccd_new, "from ^{39}Ar", "p");
    leg->AddEntry(fccd_chi2, "from ^{39}Ar chi2test", "p");
    leg->Draw();

    c1->SaveAs("ar39-results-fccd.pdf");

    TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
    c2->cd();

    tl_old->GetXaxis()->LabelsOption("v");
    tl_old->GetXaxis()->SetLabelSize(0.035);
    tl_old->GetYaxis()->SetRangeUser(0,1.1);
    tl_old->SetMarkerStyle(26);
    tl_new->SetMarkerStyle(22);
    tl_new->SetMarkerColor(kRed);
    tl_chi2->SetMarkerStyle(4);
    tl_chi2->SetMarkerColor(kBlue+1);

    tl_old->Draw("p0");
    tl_new->Draw("p0 same");
    tl_chi2->Draw("p0 same");
    gPad->SetGridx();

    c2->SaveAs("ar39-results-dlf.pdf");
}
