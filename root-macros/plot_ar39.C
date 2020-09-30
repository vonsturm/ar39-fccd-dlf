#include "../utils/poisson_bands.hpp"

bool bands = false;
bool residuals = false;
int freec = 10000;
int rebin = 4;
double ar39_scale = 1.39; // 1.35
double emin = 45, emax = 150;

std::map<std::string, double> activity = {
    {"/lar/sur_array_1/Ar39/pdf-lar-sur_array_1-Ar39.root", 4.02E10*ar39_scale},
    {"/gedet/intrinsic_bege/2nbb/pdf-gedet-intrinsic_bege-2nbb.root", 54469},
    {"/gedet/intrinsic_semicoax/2nbb/pdf-gedet-intrinsic_semicoax-2nbb.root", 41334},
    {"/gedet/intrinsic_invcoax/2nbb/pdf-gedet-intrinsic_invcoax-2nbb.root", 30236},
    {"/lar/outside_ms/K42/pdf-lar-outside_ms-K42.root", 9.8677E07},
    {"/minishroud/ms_all/K40/pdf-minishroud-ms_all-K40.root", 7.4378E+05},
    {"/cables/cables_all/Pb214/pdf-cables-cables_all-Pb214.root", 27281},
    {"/cables/cables_all/Bi214/pdf-cables-cables_all-Bi214.root", 27281},
};

void init_colors() {
    float f = 10;
    for (int i = 0; i < 12; ++i) new TColor(freec+i, (0+i*f)/255., (51+i*f)/255., (102+i*f)/255., "blue");
}

void stack_hist(std::string histname, std::vector<std::string> cycles) {

    auto c = new TCanvas(("c_" + histname).c_str(), histname.c_str(), 600, 600);
    c->SetLeftMargin(0.13);
    auto leg = new TLegend(0.5, 0.75, 0.9 , 0.9);

    std::string file = "../data/gerda-data-bkgmodel-phaseIIplus-v07.00.root";
    auto tfd = TFile::Open(file.c_str());
    TH1* data = nullptr;
    if (tfd) {
        data = (TH1*)tfd->Get(histname.c_str());
        data->SetLineColor(kBlack);
        data->SetMarkerStyle(7);
        data->Rebin(rebin);
        data->GetXaxis()->SetRangeUser(emin, emax);
        data->SetTitle(Form("%s (%s)", data->GetTitle(), data->GetName()));
        data->SetXTitle("energy [keV]");
        data->SetYTitle(("counts / " + std::to_string(rebin) + " keV").c_str());
        data->SetTitleOffset(1.8, "Y");
        leg->AddEntry(data, "data", "p");
    }

    bool init = false;
    for (auto c : cycles) {

        // get and scale pdfs
        TH1* th = nullptr;
        for (auto bkg : activity) {
            file = c + bkg.first;
            auto tf = TFile::Open(file.c_str());
            if (!tf) continue;

            auto prim = ((TParameter<Long64_t>*)tf->Get("NumberOfPrimariesEdep"))->GetVal();
            auto _th = (TH1*)tf->Get(histname.c_str());
            _th->Scale(bkg.second/prim);

            if (!th) th = _th;
            else th->Add(_th);
        }
        if (!th) continue;

        th->Rebin(rebin);
        th->SetLineColor(freec++);
        //th->SetLineColor(kBlack);

        TH1* data_clone = nullptr;
        if (data) {
            data_clone = (TH1*)data->Clone();
            if (residuals) data_clone->Divide(th);
            if (!init) {
                data_clone->Draw("hist p0");
                init = true;
            }
            else data_clone->Draw("hist p0 same");
        }

        if (!residuals) {
            if (!init) {
                th->SetTitle(Form("%s (%s)", data_clone->GetTitle(), data_clone->GetName()));
                th->SetXTitle("energy [keV]");
                th->SetYTitle(("counts / " + std::to_string(rebin) + " keV").c_str());
                th->SetTitleOffset(1.8, "Y");
                th->Draw("hist");
                // gPad->SetLogy();
                th->GetXaxis()->SetRangeUser(emin, emax);
                init = true;
            }
            else th->Draw("hist same");
        }

        if (bands) {
            for (int i = th->FindBin(emin); i <= th->FindBin(emax); ++i) {
                poiband::draw_poisson_bands(
                    th->GetBinContent(i),
                    th->GetBinLowEdge(i),
                    th->GetBinWidth(i),
                    residuals
                );
            }
            if (!residuals) th->Draw("same hist");
        }

        if (data) data_clone->Draw("same hist p0");

        leg->AddEntry(th, c.c_str(), "l");

        auto func = new TF1((c + "_fit").c_str(), "pol3", 70, 120);
        th->Fit(func, "LRQ");
//        func->Draw("same");
        std::cout << c << " -> max at " << func->GetMaximumX() << std::endl;

        auto func_min = new TF1((c + "_fit_min").c_str(), "pol3", 40, 80);
        th->Fit(func_min, "LRQ");
//        func_min->Draw("same");
        std::cout << c << " -> min at " << func_min->GetMinimumX() << std::endl;

        th->Scale(1./th->GetBinContent(th->FindBin(func->GetMaximumX())));

    }

    leg->Draw();
}

void rainbow_plot_ar39() {

    init_colors();
    bands = false;
    residuals = false;

    std::vector<std::string> dls, tls;

    dls = {
        /* "100", "200", */ "300", "400",
        "500", "600", "700", "800", "900", "1000",
        "1100", "1200", "1300", "1400", "1500"
    };
    tls = {
        "00", "01", "02", "03", "04",
        "05", "06", "07", "08", "09",
        "10"
    };

    dls = {"1500"};
    // tls = {"10"};
    // tls = {"00"};
    // tls = {"10"};

    std::vector<std::string> cycles;
    for (auto dl : dls) {
        for (auto tl : tls) {
            cycles.push_back("nplus-dl" + dl + "um-tl" + tl);
        }
    }

    // cycles = {
    //     "pplus-1um",
    //     "pplus-10um",
    //     "pplus-20um",
    //     "pplus-50um",
    //     "pplus-100um",
    //     "pplus-200um",
    //     "pplus-500um"
    // };

    stack_hist("raw/M1_enrBEGe", cycles);
}

void plot_ar39() {

    init_colors();
    bands = true;
    residuals = false;

// chi2 analysis
    stack_hist("raw/M1_ch0",  {"nplus-dl1800um-tl06"});
    stack_hist("raw/M1_ch1",  {"nplus-dl1300um-tl07"});
    stack_hist("raw/M1_ch2",  {"nplus-dl1800um-tl06"});
    stack_hist("raw/M1_ch3",  {"nplus-dl1600um-tl06"});
    stack_hist("raw/M1_ch4",  {"nplus-dl1600um-tl06"});
    stack_hist("raw/M1_ch5",  {"nplus-dl1500um-tl07"});
    stack_hist("raw/M1_ch6",  {}); // GD02D
    stack_hist("raw/M1_ch7",  {"nplus-dl1400um-tl07"});
    stack_hist("raw/M1_ch8",  {}); // ANG5
    stack_hist("raw/M1_ch9",  {"nplus-dl2500um-tl06"}); // TL linear model not good
    stack_hist("raw/M1_ch10", {"nplus-dl1800um-tl08"});
    stack_hist("raw/M1_ch11", {"nplus-dl1800um-tl06"}); //ok
    stack_hist("raw/M1_ch12", {"nplus-dl1700um-tl06"});
    stack_hist("raw/M1_ch13", {"nplus-dl2000um-tl06"});
    stack_hist("raw/M1_ch14", {"nplus-dl1400um-tl06"});
    stack_hist("raw/M1_ch15", {"nplus-dl1400um-tl07"});
    stack_hist("raw/M1_ch16", {"nplus-dl1100um-tl06"});
    stack_hist("raw/M1_ch17", {"nplus-dl1300um-tl06"});
    stack_hist("raw/M1_ch18", {"nplus-dl1200um-tl07"});
    stack_hist("raw/M1_ch19", {"nplus-dl1200um-tl07"});
    stack_hist("raw/M1_ch20", {"nplus-dl1300um-tl07"});
    stack_hist("raw/M1_ch21", {"nplus-dl1200um-tl06"});
    stack_hist("raw/M1_ch22", {"nplus-dl1100um-tl07"});
    stack_hist("raw/M1_ch23", {"nplus-dl1200um-tl07"});
    stack_hist("raw/M1_ch24", {"nplus-dl1600um-tl06"});
    stack_hist("raw/M1_ch25", {"nplus-dl1400um-tl05"});
    stack_hist("raw/M1_ch26", {"nplus-dl1100um-tl06"});
    stack_hist("raw/M1_ch27", {"nplus-dl1800um-tl08"});
    stack_hist("raw/M1_ch28", {"nplus-dl1500um-tl08"});
    stack_hist("raw/M1_ch29", {"nplus-dl1700um-tl08"});
    stack_hist("raw/M1_ch30", {"nplus-dl1700um-tl06"});
    stack_hist("raw/M1_ch31", {"nplus-dl1500um-tl06"});
    stack_hist("raw/M1_ch32", {"nplus-dl1300um-tl07"});
    stack_hist("raw/M1_ch33", {"nplus-dl1300um-tl06"});
    stack_hist("raw/M1_ch34", {"nplus-dl1000um-tl06"});
    stack_hist("raw/M1_ch35", {"nplus-dl1700um-tl06"});
    stack_hist("raw/M1_ch36", {}); // IC48B
    stack_hist("raw/M1_ch37", {"nplus-dl800um-tl04"});
    stack_hist("raw/M1_ch38", {"nplus-dl800um-tl02"});
    stack_hist("raw/M1_ch39", {"nplus-dl1000um-tl03"});
    stack_hist("raw/M1_ch40", {"nplus-dl1200um-tl04"});

/* by-eye values Luigi
    stack_hist("raw/M1_ch0",  {"nplus-dl1600um-tl06"});
    stack_hist("raw/M1_ch1",  {"nplus-dl1300um-tl07"});
    stack_hist("raw/M1_ch2",  {"nplus-dl1600um-tl06"});
    stack_hist("raw/M1_ch3",  {"nplus-dl1600um-tl06"});
    stack_hist("raw/M1_ch4",  {"nplus-dl1500um-tl06"});
    stack_hist("raw/M1_ch5",  {"nplus-dl1600um-tl06"});
    stack_hist("raw/M1_ch6",  {}); // GD02D
    stack_hist("raw/M1_ch7",  {"nplus-dl1400um-tl06"});
    stack_hist("raw/M1_ch8",  {}); // ANG5
    stack_hist("raw/M1_ch9",  {"nplus-dl2200um-tl07"}); // TL linear model not good
    stack_hist("raw/M1_ch10", {"nplus-dl1900um-tl07"});
    stack_hist("raw/M1_ch11", {"nplus-dl1500um-tl06"});
    stack_hist("raw/M1_ch12", {"nplus-dl1500um-tl06"});
    stack_hist("raw/M1_ch13", {"nplus-dl1700um-tl06"});
    stack_hist("raw/M1_ch14", {"nplus-dl1300um-tl06"});
    stack_hist("raw/M1_ch15", {"nplus-dl1500um-tl06"});
    stack_hist("raw/M1_ch16", {"nplus-dl1400um-tl05"});
    stack_hist("raw/M1_ch17", {"nplus-dl1300um-tl06"});
    stack_hist("raw/M1_ch18", {"nplus-dl1400um-tl06"});
    stack_hist("raw/M1_ch19", {"nplus-dl1100um-tl07"});
    stack_hist("raw/M1_ch20", {"nplus-dl1300um-tl07"});
    stack_hist("raw/M1_ch21", {"nplus-dl1400um-tl06"});
    stack_hist("raw/M1_ch22", {"nplus-dl1200um-tl07"});
    stack_hist("raw/M1_ch23", {"nplus-dl1300um-tl06"});
    stack_hist("raw/M1_ch24", {"nplus-dl1600um-tl06"});
    stack_hist("raw/M1_ch25", {"nplus-dl1500um-tl05"});
    stack_hist("raw/M1_ch26", {"nplus-dl1400um-tl06"});
    stack_hist("raw/M1_ch27", {"nplus-dl2100um-tl08"});
    stack_hist("raw/M1_ch28", {"nplus-dl1800um-tl07"});
    stack_hist("raw/M1_ch29", {"nplus-dl2000um-tl08"});
    stack_hist("raw/M1_ch30", {"nplus-dl1500um-tl06"});
    stack_hist("raw/M1_ch31", {"nplus-dl1400um-tl05"});
    stack_hist("raw/M1_ch32", {"nplus-dl1300um-tl07"});
    stack_hist("raw/M1_ch33", {"nplus-dl1300um-tl06"});
    stack_hist("raw/M1_ch34", {"nplus-dl1300um-tl06"});
    stack_hist("raw/M1_ch35", {"nplus-dl1700um-tl06"});
    stack_hist("raw/M1_ch36", {}); // IC48B
    stack_hist("raw/M1_ch37", {"nplus-dl800um-tl04"});
    stack_hist("raw/M1_ch38", {"nplus-dl800um-tl02"});
    stack_hist("raw/M1_ch39", {"nplus-dl1000um-tl03"});
    stack_hist("raw/M1_ch40", {"nplus-dl1100um-tl04"});
*/
}
