// ============================================================================
// ScaleSmear_RDF.C
// Full corrected version of Matteo Bonanomi's python analysis using ROOT RDataFrame
// ============================================================================

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <iostream>
#include "CMS_lumi.h"
#include "CMS_lumi.C"

using namespace ROOT;
using namespace ROOT::VecOps;

// ------------------------------------------------------------
// Select Z flavor (mu/e/incl)
bool selectChannel(int zflav, const std::string &CHANNEL) {
    int zf = std::abs(zflav);

    if (CHANNEL == "mu") return (zf == 169);
    if (CHANNEL == "e") return (zf == 121);
    if (CHANNEL == "incl") return (zf == 169 || zf == 121);

    throw std::runtime_error("Invalid channel");
}

// ------------------------------------------------------------
// Compute genEventSumw from ROOT::TFile
double getGenEventSumw(const std::string &fname) {
    TFile f(fname.c_str());
    if (f.IsZombie()) {
        std::cerr << "Cannot open: " << fname << std::endl;
        return 0;
    }

    TTree *Runs = (TTree*)f.Get("Runs");
    if (!Runs) {
        std::cerr << "No Runs tree in: " << fname << std::endl;
        return 0;
    }

    double genEventSumw = 0;
    Long64_t genEventCountSum = 0;

    double sumw;
    Long64_t cnt;

    Runs->SetBranchAddress("genEventSumw", &sumw);
    Runs->SetBranchAddress("genEventCount", &cnt);

    for (Long64_t i=0; i<Runs->GetEntries(); i++) {
        Runs->GetEntry(i);
        genEventSumw     += sumw;
        genEventCountSum += cnt;
    }

    std::cout << "File: " << fname 
              << "  genEventSumw=" << genEventSumw
              << "  genEventCount=" << genEventCountSum << std::endl;

    return genEventSumw;  // we only use sumw downstream
}

// ------------------------------------------------------------
// Lambda to extract Z1 mass and flavor
auto pickZ1 = [](short bestZIdx, const RVec<float>& Zmass, const RVec<int>& Zflav) {
    if (bestZIdx < 0 || bestZIdx >= (int)Zmass.size()) return -1.f;
    return Zmass[bestZIdx];
};

auto pickFlav = [](short bestZIdx, const RVec<int>& Zflav) {
    if (bestZIdx < 0 || bestZIdx >= (int)Zflav.size()) return -999;
    return Zflav[bestZIdx];
};

// ============================================================================
// MAIN
// ============================================================================
void ScaleSmear_RDF() {

    // --------------------------
    // Config
    // --------------------------
    std::string YEAR    = "2022preEE";
    std::string CHANNEL = "mu";
    double lumi         = 7.98;  // fb^-1 for 2022 preEE

    // --------------------------
    // File paths
    // --------------------------
    std::string base = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2022_MC";
    std::string baseData = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2022_Data";

    std::string fname_dy = base + "/DYJetsToLL/ZZ4lAnalysis.root";
    std::string fname_tt = base + "/TTto2L2Nu/ZZ4lAnalysis.root";
    std::string fname_wz = base + "/WZto3LNu/ZZ4lAnalysis.root";
    std::string fname_dt = baseData + "/Data_eraCD_preEE.root";

    // --------------------------
    // Read genEventSumw
    // --------------------------
    double cnt_dy = getGenEventSumw(fname_dy);
    double cnt_tt = getGenEventSumw(fname_tt);
    double cnt_wz = getGenEventSumw(fname_wz);

    if (cnt_dy == 0 || cnt_tt == 0 || cnt_wz == 0)
        throw std::runtime_error("genEventSumw = 0 for some sample!");

    // --------------------------
    // Define histogram binning
    // --------------------------
    const int NBINS = 99;
    const double XMIN = 60;
    const double XMAX = 120;

    TH1D *hDY = new TH1D("hDY","",NBINS,XMIN,XMAX);
    TH1D *hTT = new TH1D("hTT","",NBINS,XMIN,XMAX);
    TH1D *hWZ = new TH1D("hWZ","",NBINS,XMIN,XMAX);
    TH1D *hDT = new TH1D("hDT","",NBINS,XMIN,XMAX);

    // ------------------------------------------------------------
    // RDataFrames
    // ------------------------------------------------------------
    RDataFrame d_dy("Events", fname_dy);
    RDataFrame d_tt("Events", fname_tt);
    RDataFrame d_wz("Events", fname_wz);
    RDataFrame d_dt("Events", fname_dt);

    // ------------------------------------------------------------
    // Define Z1 quantities and GoodZ1 selection
    auto defZ = [&](RDataFrame &df) {
        return df
            .Define("Z1_mass" , pickZ1, {"bestZIdx","ZCand_mass","ZCand_flav"})
            .Define("Z1_flav" , pickFlav, {"bestZIdx","ZCand_flav"})
            .Define("GoodZ1", [&](int flav){ return selectChannel(flav, CHANNEL); }, {"Z1_flav"});
    };

    // DY
    auto D_dy = defZ(d_dy)
        .Filter("GoodZ1")
        .Define("w", [=](float w){ return w*1000.*lumi/cnt_dy; }, {"overallEventWeight"});

    // TT
    auto D_tt = defZ(d_tt)
        .Filter("GoodZ1")
        .Define("w", [=](float w){ return w*1000.*lumi/cnt_tt; }, {"overallEventWeight"});

    // WZ
    auto D_wz = defZ(d_wz)
        .Filter("GoodZ1")
        .Define("w", [=](float w){ return w*1000.*lumi/cnt_wz; }, {"overallEventWeight"});

    // Data
    auto D_dt = defZ(d_dt).Filter("GoodZ1");

    // ------------------------------------------------------------
    // Fill histograms
    // ------------------------------------------------------------
    D_dy.Histo1D({"hDY","",NBINS,XMIN,XMAX},"Z1_mass","w")->Copy(*hDY);
    D_tt.Histo1D({"hTT","",NBINS,XMIN,XMAX},"Z1_mass","w")->Copy(*hTT);
    D_wz.Histo1D({"hWZ","",NBINS,XMIN,XMAX},"Z1_mass","w")->Copy(*hWZ);
    D_dt.Histo1D({"hDT","",NBINS,XMIN,XMAX},"Z1_mass")->Copy(*hDT);

    // Total MC
    TH1D *hMC = (TH1D*)hDY->Clone("hMC");
    hMC->Add(hTT);
    hMC->Add(hWZ);

    // Normalize MC to Data
    double scale = hDT->Integral() / hMC->Integral();
    hMC->Scale(scale);

    // =====================================================================
    // PLOTS
    // =====================================================================
    TCanvas *c = new TCanvas("c","",1400,1000);
    c->Divide(1,2);

    // Top pad
    c->cd(1);
    gPad->SetPad(0,0.3,1,1);
    gPad->SetBottomMargin(0.05);
    gPad->SetLeftMargin(0.15);


    // Convert #1f77b4 to RGB 0-1 scale
    Float_t r = 31.0/255.0;
    Float_t g = 119.0/255.0;
    Float_t b = 180.0/255.0;

    // Define a new ROOT color (ID = 1001, pick an unused ID)
    Int_t myBlue = TColor::GetFreeColorIndex();  
    new TColor(myBlue, r, g, b);

    // Now set it
    hMC->SetTitle("");
    hMC->SetStats(0); 
    hMC->SetLineColor(myBlue);
    hMC->GetXaxis()->SetLabelSize(0); 

    //hMC->SetLineColor(kBlue);      // blue MC line
    hMC->SetLineWidth(2);
    hMC->SetFillStyle(0);          // no fill
    hMC->Draw("HIST");

    hDT->SetMarkerStyle(20);       // black points
    hDT->SetMarkerSize(1.0);
    hDT->SetLineColor(kBlack);
    hDT->Draw("E SAME");
    hMC->GetYaxis()->SetTitle("Events/bin width");
    //hMC->GetYaxis()->CenterTitle(true);
    hMC->GetYaxis()->SetTitleOffset(0.5);
    hMC->GetYaxis()->SetTitleSize(0.07);
    //hMC->GetYaxis()->SetLabelSize(0.08);
    TLegend *leg = new TLegend(0.65,0.7,0.88,0.88);
    leg->AddEntry(hMC,"MC (DY+t#bar{t}+WZ)","l");
    leg->AddEntry(hDT,"Data","lp");
    leg->SetLineColor(kWhite);
    leg->Draw();

    CMS_lumi cmsLabel;
    cmsLabel.set_lumi((TPad*)gPad, lumi);

    // Ratio pad
    c->cd(2);
    gPad->SetPad(0,0.0,1,0.3);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.25);
    gPad->SetLeftMargin(0.15);

    TH1D *hRatio = (TH1D*)hMC->Clone("hRatio");
    hRatio->Divide(hDT);
    hRatio->SetTitle("");
    hRatio->SetStats(0);
    hRatio->SetMinimum(0.5);
    hRatio->SetMaximum(1.5);
    // Turn off automatic labels
    //hRatio->GetYaxis()->SetNdivisions(0);

    // Manually set only two labels
    //hRatio->GetYaxis()->SetBinLabel(1, "0.5");
    //hRatio->GetYaxis()->SetBinLabel(2, "1.5");
    hRatio->GetYaxis()->SetTitle("MC/Data");
    hRatio->GetYaxis()->CenterTitle(true);
    hRatio->GetYaxis()->SetTitleOffset(0.2);
    hRatio->GetYaxis()->SetTitleSize(0.15);
    hRatio->GetYaxis()->SetLabelSize(0.08);
    //hRatio->GetYaxis()->SetLabelOffset(0.1);

    hRatio->GetXaxis()->SetTitle("m_{ll} [GeV]");
    hRatio->GetXaxis()->SetTitleSize(0.12);
    hRatio->GetXaxis()->SetLabelSize(0.10);
    hRatio->SetMarkerStyle(8);       // black points
    hRatio->SetMarkerSize(1.0);
    hRatio->SetLineColor(kBlack);
    hRatio->Draw("P");

    c->SaveAs(Form("plot_%s_%s_RDF.png", CHANNEL.c_str(), YEAR.c_str()));
    c->SaveAs(Form("plot_%s_%s_RDF.pdf", CHANNEL.c_str(), YEAR.c_str()));

    std::cout << "\nSaved plot: plot_" << CHANNEL << "_" << YEAR << "_RDF\n";
}
