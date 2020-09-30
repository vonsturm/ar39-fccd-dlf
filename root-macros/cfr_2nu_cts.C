#include "json.hpp"
using namespace nlohmann;

void cfr_2nu_cts() {

    std::ifstream fph2("/home/gipert/code/gerda/gerda-pdfs/src/settings/phaseII/ged-mapping.json");
    std::ifstream fph2p("/home/gipert/code/gerda/gerda-pdfs/src/settings/phaseIIplus/ged-mapping.json");

    json jph2, jph2p;
    fph2  >> jph2;
    fph2p >> jph2p;

    jph2 = jph2["mapping"];
    jph2p = jph2p["mapping"];

    auto tfo1 = TFile::Open("ph2p-v1.0-rc/gedet/intrinsic_bege/2nbb/pdf-gedet-intrinsic_bege-2nbb.root");
    auto tfn1 = TFile::Open("ph2p-v1.0-rc-ar39av/gedet/intrinsic_bege/2nbb/pdf-gedet-intrinsic_bege-2nbb.root");
    auto tfo2 = TFile::Open("ph2p-v1.0-rc/gedet/intrinsic_semicoax/2nbb/pdf-gedet-intrinsic_semicoax-2nbb.root");
    auto tfn2 = TFile::Open("ph2p-v1.0-rc-ar39av/gedet/intrinsic_semicoax/2nbb/pdf-gedet-intrinsic_semicoax-2nbb.root");
    auto tfo3 = TFile::Open("ph2p-v1.0-rc/gedet/intrinsic_invcoax/2nbb/pdf-gedet-intrinsic_invcoax-2nbb.root");
    auto tfn3 = TFile::Open("ph2p-v1.0-rc-ar39av/gedet/intrinsic_invcoax/2nbb/pdf-gedet-intrinsic_invcoax-2nbb.root");

    std::map<std::string, TH1D*> h_old, h_new;

    for (auto j : jph2p.items()) {
        if (j.key()[0] == 'G') {
            h_old.emplace(j.key(), (TH1D*)tfo1->Get(("raw/M1_ch" + std::to_string(jph2p[j.key()]["channel"].get<int>())).c_str()));
            h_new.emplace(j.key(), (TH1D*)tfn1->Get(("raw/M1_ch" + std::to_string(jph2p[j.key()]["channel"].get<int>())).c_str()));
        }
        else if (j.key()[0] == 'A' or j.key()[0] == 'R') {
            h_old.emplace(j.key(), (TH1D*)tfo2->Get(("raw/M1_ch" + std::to_string(jph2p[j.key()]["channel"].get<int>())).c_str()));
            h_new.emplace(j.key(), (TH1D*)tfn2->Get(("raw/M1_ch" + std::to_string(jph2p[j.key()]["channel"].get<int>())).c_str()));
        }
        else if (j.key()[0] == 'I') {
            h_old.emplace(j.key(), (TH1D*)tfo3->Get(("raw/M1_ch" + std::to_string(jph2p[j.key()]["channel"].get<int>())).c_str()));
            h_new.emplace(j.key(), (TH1D*)tfn3->Get(("raw/M1_ch" + std::to_string(jph2p[j.key()]["channel"].get<int>())).c_str()));
        }
    }

    TH1D *av_old = new TH1D("av_old", "2#nu#beta#beta MC counts in [100, 2000] keV;;Counts", 41, 0, 41);
    TH1D *av_new = new TH1D("av_new", "av_new", 41, 0, 41);
    TH1D *av_ratio = new TH1D("av_ratio", "^{39}Ar AV / official AV", 41, 0, 41);

    int i = 1;
    for (auto o : h_old) {
        av_old->GetXaxis()->SetBinLabel(i, o.first.c_str());
        av_ratio->GetXaxis()->SetBinLabel(i, o.first.c_str());
        av_old->SetBinContent(i, o.second->Integral(100, 2000));
        av_new->SetBinContent(i, h_new[o.first]->Integral(50, 2000));
        if (av_old->GetBinContent(i) != 0) av_ratio->SetBinContent(i, av_new->GetBinContent(i)/av_old->GetBinContent(i));
        i++;
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    c1->cd();

    av_old->GetXaxis()->LabelsOption("v");
    av_old->SetMarkerStyle(26);
    av_new->SetMarkerStyle(22);
    av_new->SetMarkerColor(kRed);

    av_old->Draw("p");
    av_new->Draw("p same");
    gPad->SetGridx();

    auto leg = new TLegend();
    leg->AddEntry(av_old, "official values", "p");
    leg->AddEntry(av_new, "from ^{39}Ar", "p");
    leg->Draw();

    c1->SaveAs("ar39-results-2nbb.pdf");

    TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
    c2->cd();

    av_ratio->GetXaxis()->LabelsOption("v");
    av_ratio->SetMarkerStyle(20);
    av_ratio->GetYaxis()->SetRangeUser(0.9, 1.1);

    av_ratio->Draw("p");
    TLine *line = new TLine();
    line->DrawLine(0,1,41,1);

    gPad->SetGridx();

    c2->SaveAs("ar39-results-2nbb-ratio.pdf");
}
