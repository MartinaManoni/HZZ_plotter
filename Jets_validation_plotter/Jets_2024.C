// Martina Manoni plotter for Jet control plots
#include "CMS_lumi.C"
#include <tuple>
#include <vector>
#include <filesystem>
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>

// -----------------------------------------------------------------------------
// Draw Ratio Plot
void DrawRatioPlot(string name, string var_name, TCanvas *c,
                   TH1D *data, TH1D *MC, double _lumi)
{
    std::cout << "[INFO] Drawing plot (CMS-TDR palette): " << name << std::endl;

    // --------------------------------------------
    // Style
    // --------------------------------------------
    gStyle->SetOptStat(0);

    // Normalize MC to data
    double dataIntegral = data->Integral();
    double mcIntegral = MC->Integral();
    double scaleFactor = dataIntegral / mcIntegral;
    //double scaleFactor = 1.0;

    // Normalize MC (optional)
    MC->Scale(scaleFactor);

    // Canvas pads
    c->Divide(1,2);
    TPad *pad1 = (TPad*)c->cd(1);
    TPad *pad2 = (TPad*)c->cd(2);

    // Top pad formatting
    pad1->SetPad(0,0.33,1,1);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.12);
    pad1->SetRightMargin(0.05);

    // Bottom pad formatting
    pad2->SetPad(0,0,1,0.33);
    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.30);
    pad2->SetLeftMargin(0.12);
    pad2->SetRightMargin(0.05);

    // --------------------------------------------
    // TOP PAD
    // --------------------------------------------
    pad1->cd();

    // CMS-TDR inspired colors
    MC->SetFillColor(TColor::GetColor("#6baed6"));   // soft blue
    MC->SetLineColor(TColor::GetColor("#2171b5"));   // darker blue outline
    MC->SetLineWidth(2);

    data->SetMarkerStyle(20);
    data->SetMarkerColor(kBlack);
    data->SetLineColor(kBlack);

    // Axis label: Events / binWidth
    double binWidth = MC->GetXaxis()->GetBinWidth(1);
    TString yLabel = TString::Format("Events / %.2f", binWidth);
    double max = 0.;
        for (int bin = 0; bin < MC->GetSize() - 2; bin++){
          if ( MC->GetBinContent(bin + 1) > max) max = MC->GetBinContent(bin + 1);
        }

    double maxMC = MC->GetMaximum();
    double maxData = data->GetMaximum();

    MC->SetMaximum(1.35 * std::max(maxMC, maxData));

    MC->GetYaxis()->SetTitle(yLabel);
    MC->GetYaxis()->SetTitleSize(0.06);
    MC->GetYaxis()->SetLabelSize(0.05);
    MC->GetYaxis()->SetTitleOffset(0.95);

    MC->GetXaxis()->SetLabelSize(0); // hide x labels on top

    // Draw MC
    MC->Draw("HIST");

    // Draw data
    data->Draw("P E1 SAME");

    // Legend
    TLegend *leg = new TLegend(0.70,0.60,0.92,0.87);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.045);
    leg->AddEntry(MC,  "DY + t#bar{t} MC", "f");
    leg->AddEntry(data,"Data",            "p");
    leg->Draw();

    // CMS lumi text
    //CMS_lumi(pad1, 0, 0);

    // --------------------------------------------
    // RATIO PAD
    // --------------------------------------------
    pad2->cd();

    // Build ratio = data / MC
    TH1D *ratio = (TH1D*)data->Clone("ratio");
    ratio->Divide(MC);

    ratio->SetMarkerStyle(20);
    ratio->SetMarkerSize(0.7);
    ratio->SetMarkerColor(kBlack);
    ratio->SetLineColor(kBlack);

    ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetTitle("Data/MC");
    ratio->GetYaxis()->SetTitleSize(0.11);
    ratio->GetYaxis()->SetLabelSize(0.10);
    ratio->GetYaxis()->SetTitleOffset(0.5);

    ratio->GetXaxis()->SetTitle(var_name.c_str());
    ratio->GetXaxis()->SetTitleSize(0.14);
    ratio->GetXaxis()->SetLabelSize(0.12);

    ratio->Draw("P E1");

    // Unity line
    TLine *line = new TLine(MC->GetXaxis()->GetXmin(),1,
                            MC->GetXaxis()->GetXmax(),1);
    line->SetLineWidth(2);
    line->SetLineStyle(2);
    line->Draw("SAME");

    // --------------------------------------------
    // Save
    // --------------------------------------------
    c->SaveAs((name + ".pdf").c_str());
    c->SaveAs((name + ".png").c_str());
    c->SaveAs((name + ".C").c_str());

    c->Clear();
}


// -----------------------------------------------------------------------------
// Compute leading jet
// -----------------------------------------------------------------------------
double computeLeadJet(ROOT::RVec<float> &eta, ROOT::RVec<float> &pt)
{
    double highest = 0;
    int idx = -1;
    for (int i=0; i<pt.size(); i++)
        if (pt[i] > highest && pt[i] > 20) { highest = pt[i]; idx = i; }

    return (idx >= 0 ? eta[idx] : 0);
}

// -----------------------------------------------------------------------------
// Sum genEventSumw over multiple files
// -----------------------------------------------------------------------------
double get_genEventSumw(const std::vector<std::string>& files)
{
    ROOT::RDataFrame rdf("Runs", files);
    return *(rdf.Sum("genEventSumw"));
}

// -----------------------------------------------------------------------------
// Histo for MC (vector of files)
// -----------------------------------------------------------------------------
auto Histo_MC(const std::vector<std::string>& files, double _lumi,
              std::string var, int bins, double xmin, double xmax)
{
    double gen_sumWeights = get_genEventSumw(files);
    ROOT::RDataFrame rdf("Events", files);

    auto sorted = rdf
        .Define("Horn_Cut", "((abs(Jet_eta) < 2.5 || abs(Jet_eta) > 3) || Jet_pt >50)")
        .Define("Jet_Mask", "(Jet_pt > 30) && Horn_Cut && (Jet_jetId == 6 || Jet_jetId == 2) && Jet_ZZMask == false ") //&& Jet_ZZMask == false
        .Define("Jet_Mask_smearUp",   "(Jet_smearUp_pt > 30) && Horn_Cut && (Jet_jetId == 6 || Jet_jetId == 2) && Jet_ZZMask == false ") //&& Jet_ZZMask == false
        .Define("Jet_Mask_smearDown", "(Jet_smearDn_pt > 30) && Horn_Cut && (Jet_jetId == 6 || Jet_jetId == 2) && Jet_ZZMask == false") //&& Jet_ZZMask == false

        .Define("FilteredJet_eta", "Jet_eta[Jet_Mask]")
        .Define("FilteredJet_pt",  "Jet_pt[Jet_Mask]")

        .Define("FilteredJet_smearUp_pt",   "Jet_smearUp_pt[Jet_Mask]")
        .Define("FilteredJet_smearDn_pt",   "Jet_smearDn_pt[Jet_Mask]")

        .Define("SortedJet_eta", "FilteredJet_eta[Reverse(Argsort(FilteredJet_pt))]")
        .Define("SortedJet_pt",  "FilteredJet_pt[Reverse(Argsort(FilteredJet_pt))]")

        .Define("SortedJet_smearUp_pt", "FilteredJet_smearUp_pt[Reverse(Argsort(FilteredJet_pt))]")
        .Define("SortedJet_smearDn_pt", "FilteredJet_smearDn_pt[Reverse(Argsort(FilteredJet_pt))]");

    auto skim = sorted
        .Filter("nZCand >= 0 && SortedJet_eta.size() > 0 && Flag_JetVetoed == 0")
        .Define("weight", Form("1000 * overallEventWeight * %f / %f", _lumi, gen_sumWeights))
        .Define("leadEta",  "SortedJet_eta.at(0)")
        .Define("leadPt",   "SortedJet_pt.at(0)")
        .Define("leadEtaUp",   computeLeadJet, {"SortedJet_eta", "SortedJet_smearUp_pt"})
        .Define("leadEtaDown", computeLeadJet, {"SortedJet_eta", "SortedJet_smearDn_pt"})
        .Define("leadPtUp",    computeLeadJet, {"SortedJet_pt", "SortedJet_smearUp_pt"})
        .Define("leadPtDown",  computeLeadJet, {"SortedJet_pt", "SortedJet_smearDn_pt"});

    auto h_nom = skim.Histo1D({(var+"_MC").c_str(), (var+"_MC").c_str(), bins, xmin, xmax}, var, "weight");
    auto h_up  = skim.Histo1D({(var+"Up_MC").c_str(), (var+"Up_MC").c_str(), bins, xmin, xmax}, var+"Up", "weight");
    auto h_dn  = skim.Histo1D({(var+"Down_MC").c_str(), (var+"Down_MC").c_str(), bins, xmin, xmax}, var+"Down", "weight");

    return std::make_tuple(h_nom, h_up, h_dn);
}

// -----------------------------------------------------------------------------
// DATA HISTO (multiple files)
// -----------------------------------------------------------------------------
ROOT::RDF::RResultPtr<TH1D>
Histo_Data(const std::vector<std::string>& files,
           std::string var, int bins, double xmin, double xmax)
{
    ROOT::RDataFrame rdf("Events", files);

    auto sorted = rdf
        .Define("Horn_Cut", "((abs(Jet_eta) < 2.5 || abs(Jet_eta) > 3) || Jet_pt > 50)")
        .Define("Jet_Mask", "(Jet_pt >30) && Horn_Cut && (Jet_jetId == 6 || Jet_jetId == 2) && Jet_ZZMask == false") //&& Jet_ZZMask == false
        .Define("FilteredJet_eta", "Jet_eta[Jet_Mask]")
        .Define("FilteredJet_pt",  "Jet_pt[Jet_Mask]")
        .Define("SortedJet_eta", "FilteredJet_eta[Reverse(Argsort(FilteredJet_pt))]")
        .Define("SortedJet_pt",  "FilteredJet_pt[Reverse(Argsort(FilteredJet_pt))]");

    auto skim = sorted
        .Filter("nZCand >= 0 && SortedJet_eta.size() > 0 && Flag_JetVetoed == 0")
        .Define("leadEta", "SortedJet_eta.at(0)")
        .Define("leadPt",  "SortedJet_pt.at(0)");

    return skim.Histo1D({(var+"_data").c_str(), (var+"_data").c_str(), bins, xmin, xmax}, var);
}

// -----------------------------------------------------------------------------
// MAIN
// -----------------------------------------------------------------------------
void Jets_2024()
{
    ROOT::EnableImplicitMT();

    std::string var = "leadPt";
    std::string var_name = "Leading #p_T";
    std::string period = "2024";

    std::map<std::string, std::tuple<int,double,double>> var_binning;
    var_binning["leadEta"] = {47, -4.7, 4.7};
    var_binning["leadPt"] = {40, 30, 200};
    var_binning["leadNjets"] = {8, 2, 10};

    int n_bins;
    double xmin, xmax;
    std::tie(n_bins,xmin,xmax) = var_binning[var];

    // -----------------------------------------------------------
    // INPUT MULTIPLE FILES HERE
    // -----------------------------------------------------------

    std::vector<std::string> DY_files = {
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Tau/DYJetsTo2Tau_0/DYJetsTo2Tau/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Tau/DYJetsTo2Tau_1/DYJetsTo2Tau/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Tau/DYJetsTo2Tau_2/DYJetsTo2Tau/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Tau/DYJetsTo2Tau_3/DYJetsTo2Tau/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Tau/DYJetsTo2Tau_4/DYJetsTo2Tau/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Tau/DYJetsTo2Tau_5/DYJetsTo2Tau/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Tau/DYJetsTo2Tau_6/DYJetsTo2Tau/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Tau/DYJetsTo2Tau_7/DYJetsTo2Tau/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Tau/DYJetsTo2Tau_8/DYJetsTo2Tau/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Tau/DYJetsTo2Tau_9/DYJetsTo2Tau/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Tau/DYJetsTo2Tau_10/DYJetsTo2Tau/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Tau/DYJetsTo2Tau_11/DYJetsTo2Tau/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Tau/DYJetsTo2Tau_12/DYJetsTo2Tau/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Tau/DYJetsTo2Tau_13/DYJetsTo2Tau/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Tau/DYJetsTo2Tau_14/DYJetsTo2Tau/ZZ4lAnalysis.root",

        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Mu/DYJetsTo2Mu_0/DYJetsTo2Mu/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Mu/DYJetsTo2Mu_1/DYJetsTo2Mu/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Mu/DYJetsTo2Mu_2/DYJetsTo2Mu/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Mu/DYJetsTo2Mu_3/DYJetsTo2Mu/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Mu/DYJetsTo2Mu_4/DYJetsTo2Mu/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Mu/DYJetsTo2Mu_5/DYJetsTo2Mu/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Mu/DYJetsTo2Mu_6/DYJetsTo2Mu/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Mu/DYJetsTo2Mu_7/DYJetsTo2Mu/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Mu/DYJetsTo2Mu_8/DYJetsTo2Mu/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Mu/DYJetsTo2Mu_9/DYJetsTo2Mu/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Mu/DYJetsTo2Mu_10/DYJetsTo2Mu/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Mu/DYJetsTo2Mu_11/DYJetsTo2Mu/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Mu/DYJetsTo2Mu_12/DYJetsTo2Mu/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Mu/DYJetsTo2Mu_13/DYJetsTo2Mu/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2Mu/DYJetsTo2Mu_14/DYJetsTo2Mu/ZZ4lAnalysis.root",

        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2E/DYJetsTo2E_0/DYJetsTo2E/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2E/DYJetsTo2E_1/DYJetsTo2E/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2E/DYJetsTo2E_2/DYJetsTo2E/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2E/DYJetsTo2E_3/DYJetsTo2E/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2E/DYJetsTo2E_4/DYJetsTo2E/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2E/DYJetsTo2E_5/DYJetsTo2E/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2E/DYJetsTo2E_6/DYJetsTo2E/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2E/DYJetsTo2E_7/DYJetsTo2E/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2E/DYJetsTo2E_8/DYJetsTo2E/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2E/DYJetsTo2E_9/DYJetsTo2E/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2E/DYJetsTo2E_10/DYJetsTo2E/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2E/DYJetsTo2E_11/DYJetsTo2E/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2E/DYJetsTo2E_12/DYJetsTo2E/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/DYJetsTo2E/DYJetsTo2E_13/DYJetsTo2E/ZZ4lAnalysis.root",
        
    };

    std::vector<std::string> TT_files = {
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/TTto2L2Nu/TTto2L2Nu_0/TTto2L2Nu/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/TTto2L2Nu/TTto2L2Nu_1/TTto2L2Nu/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/TTto2L2Nu/TTto2L2Nu_2/TTto2L2Nu/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/TTto2L2Nu/TTto2L2Nu_3/TTto2L2Nu/ZZ4lAnalysis.root"
    };

    std::vector<std::string> data_files = {
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/EGamma02024Cv1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/EGamma02024Dv1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/EGamma02024Ev1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/EGamma02024Fv1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/EGamma02024Gv2/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/EGamma02024Hv2/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/EGamma02024Iv1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/EGamma02024Iv2/ZZ4lAnalysis.root",

    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/EGamma12024Cv1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/EGamma12024Dv1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/EGamma12024Ev1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/EGamma12024Fv1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/EGamma12024Gv2/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/EGamma12024Hv1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/EGamma12024Iv1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/EGamma12024Iv2/ZZ4lAnalysis.root",

    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/Muon02024Cv1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/Muon02024Dv1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/Muon02024Ev1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/Muon02024Fv1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/Muon02024Gv1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/Muon02024Hv1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/Muon02024Iv1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/Muon02024Iv2/ZZ4lAnalysis.root",

    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/Muon12024Cv1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/Muon12024Dv1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/Muon12024Ev1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/Muon12024Fv1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/Muon12024Gv2/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/Muon12024Hv2/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/Muon12024Iv1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/Muon12024Iv2/ZZ4lAnalysis.root",

    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/MuonEG2024Cv1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/MuonEG2024Dv1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/MuonEG2024Ev1/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/MuonEG2024Fv2/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/MuonEG2024Gv3/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/MuonEG2024Hv2/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/MuonEG2024Iv2/ZZ4lAnalysis.root",
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/MuonEG2024Iv2v2/ZZ4lAnalysis.root"
    };


    double lumi = 108.8;

    // MC
    auto [h_DY, h_DY_up, h_DY_dn] = Histo_MC(DY_files, lumi, var, n_bins, xmin, xmax);
    auto [h_TT, h_TT_up, h_TT_dn] = Histo_MC(TT_files, lumi, var, n_bins, xmin, xmax);

    // Combine DY + TT
    h_DY->Scale(3);
    h_DY->Add(h_TT.GetPtr());
    h_DY_up->Add(h_TT_up.GetPtr());
    h_DY_dn->Add(h_TT_dn.GetPtr());

    // DATA
    auto h_data = Histo_Data(data_files, var, n_bins, xmin, xmax);

    // Save ROOT file
    std::filesystem::create_directories("2024/");
    TFile outfile(("2024/" + var + "_histos.root").c_str(),"RECREATE");
    h_data->Write();
    h_DY->Write();
    h_DY_up->Write();
    h_DY_dn->Write();
    outfile.Close();

    //h_DY_up.GetPtr(),
    //h_DY_dn.GetPtr(),

    // Draw plots
    TCanvas *c = new TCanvas();
    DrawRatioPlot("2024/" + var + "_plot", var_name, c,
                  h_data.GetPtr(),
                  h_DY.GetPtr(),
                  lumi);
}

