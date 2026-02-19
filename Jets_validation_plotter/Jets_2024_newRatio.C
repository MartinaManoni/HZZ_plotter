// Martina Manoni plotter for Jet control plots, reads directly each chunck 
#include "CMS_lumi.C"
#include <tuple>
#include <vector>
#include <filesystem>
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>



#include <filesystem>
#include <TFile.h>

namespace fs = std::filesystem;

std::vector<std::string> collect_root_files(const std::string& baseDir)
{
    std::vector<std::string> files;

    for (const auto& entry : fs::recursive_directory_iterator(baseDir)) {
        if (!entry.is_regular_file()) continue;

        const auto& path = entry.path();
        if (path.filename() == "ZZ4lAnalysis.root") {
            files.push_back(path.string());
        }
    }
    return files;
}


// -----------------------------------------------------------------------------
void DrawRatioPlot(string name, string var_name, TCanvas *c,
                   TH1D *data,
                   TH1D *MC,
                   TH1D *MCUp,
                   TH1D *MCDn,
                   TH1D *MCUpscale,
                   TH1D *MCDnscale,
                   double _lumi)
{
    gStyle->SetOptStat(0);

    // =========================================================
    // 1) Normalize MC FIRST (important!)
    // =========================================================
    double scaleFactor = data->Integral() / MC->Integral();

    MC->Scale(scaleFactor);
    MCUp->Scale(scaleFactor);
    MCDn->Scale(scaleFactor);
    MCUpscale->Scale(scaleFactor);
    MCDnscale->Scale(scaleFactor);

    // =========================================================
    // 2) Build total MC uncertainty envelopes
    // =========================================================
    TH1D *MCUpFinal = (TH1D*) MC->Clone("MCUpFinal");
    TH1D *MCDnFinal = (TH1D*) MC->Clone("MCDnFinal");

    MCUpFinal->Reset();
    MCDnFinal->Reset();

    for (int bin = 1; bin <= MC->GetNbinsX(); bin++) {

        double nom = MC->GetBinContent(bin);
        double mc_stat = MC->GetBinError(bin);  // sqrt(sum w^2)

        double dJEC_up   = MCUp->GetBinContent(bin)      - nom;
        double dJER_up   = MCUpscale->GetBinContent(bin) - nom;
        double dJEC_down = nom - MCDn->GetBinContent(bin);
        double dJER_down = nom - MCDnscale->GetBinContent(bin);

        double total_up = std::sqrt(mc_stat*mc_stat +
                                    dJEC_up*dJEC_up +
                                    dJER_up*dJER_up);

        double total_dn = std::sqrt(mc_stat*mc_stat +
                                    dJEC_down*dJEC_down +
                                    dJER_down*dJER_down);

        MCUpFinal->SetBinContent(bin, nom + total_up);
        MCDnFinal->SetBinContent(bin, std::max(0.0, nom - total_dn));
    }

    // =========================================================
    // 3) Canvas & top pad (unchanged logic)
    // =========================================================
    c->Divide(1,2);

    c->cd(1);
    gPad->SetBottomMargin(0.02);
    gPad->SetTopMargin(0.18);
    gPad->SetLeftMargin(0.10);

    MC->SetFillColor(TColor::GetColor("#6baed6"));
    MC->SetLineColor(TColor::GetColor("#2171b5"));
    MC->SetLineWidth(2);
    MC->Draw("HIST");

    data->SetMarkerStyle(20);
    data->SetMarkerSize(0.8);
    data->Draw("P SAME");

    CMS_lumi cms;
    cms.set_lumi((TPad*)gPad, _lumi);

    // =========================================================
    // 4) Ratio pad — EXACT Python equivalent
    // =========================================================
    c->cd(2);
    gPad->SetBottomMargin(0.20);
    gPad->SetTopMargin(0.02);
    gPad->SetLeftMargin(0.10);

    // --- Data / MC ---
    TH1D *Ratio = (TH1D*) data->Clone("Ratio");
    Ratio->Divide(MC);
    Ratio->SetMarkerStyle(20);
    Ratio->SetMarkerSize(0.8);

    // --- Ratio uncertainty band ---
    TH1D *ratio_band_up = (TH1D*) MCUpFinal->Clone("ratio_band_up");
    TH1D *ratio_band_dn = (TH1D*) MCDnFinal->Clone("ratio_band_dn");

    for (int bin = 1; bin <= MC->GetNbinsX(); bin++) {

        double mc = MC->GetBinContent(bin);
        if (mc <= 0) {
            ratio_band_up->SetBinContent(bin, 1.0);
            ratio_band_dn->SetBinContent(bin, 1.0);
            continue;
        }

        ratio_band_up->SetBinContent(bin, MCUpFinal->GetBinContent(bin) / mc);
        ratio_band_dn->SetBinContent(bin, MCDnFinal->GetBinContent(bin) / mc);
    }

    ratio_band_up->SetFillColor(kGray+1);
    ratio_band_up->SetLineColor(kGray+1);
    ratio_band_up->SetMinimum(0.5);
    ratio_band_up->SetMaximum(1.5);
    ratio_band_up->GetYaxis()->SetTitle("Data / MC");
    ratio_band_up->Draw("HIST");

    ratio_band_dn->SetFillColor(kGray);
    ratio_band_dn->SetLineColor(kGray);
    ratio_band_dn->Draw("HIST SAME");

    Ratio->Draw("P SAME");

    TLine *line = new TLine(
        Ratio->GetXaxis()->GetXmin(), 1.0,
        Ratio->GetXaxis()->GetXmax(), 1.0
    );
    line->SetLineStyle(2);
    line->Draw();

    // =========================================================
    // 5) Legend
    // =========================================================
    c->cd(1);
    TLegend *leg = new TLegend(0.64,0.45,0.85,0.75);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(MC, "DY + t#bar{t} MC", "f");
    leg->AddEntry(data, "Data", "p");
    leg->AddEntry(ratio_band_up, "Stat. ⊕ Syst.", "f");
    leg->Draw();

    c->SaveAs((TString)name + ".pdf");
    c->SaveAs((TString)name + ".png");
}



// -----------------------------------------------------------------------------
// Compute leading jet
// -----------------------------------------------------------------------------
double computeLeadJet(ROOT::RVec<float> &eta, ROOT::RVec<float> &pt)
{
    double highest = 0;
    int idx = -1;
    for (int i=0; i<pt.size(); i++)
        if (pt[i] > highest && pt[i] > 30) { highest = pt[i]; idx = i; }

    return (idx >= 0 ? eta[idx] : 0);
}

// -----------------------------------------------------------------------------
// Sum genEventSumw over multiple files
// -----------------------------------------------------------------------------
double get_genEventSumw(const std::vector<std::string>& files)
{
    ROOT::EnableImplicitMT();
    ROOT::RDF::RSnapshotOptions opts;
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
    ROOT::EnableImplicitMT();
    ROOT::RDF::RSnapshotOptions opts;
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
        .Define("FilteredJet_scaleUp_pt", "Jet_scaleUp_pt[Jet_Mask]")
        .Define("FilteredJet_scaleDn_pt", "Jet_scaleDn_pt[Jet_Mask]")

        .Define("SortedJet_eta", "FilteredJet_eta[Reverse(Argsort(FilteredJet_pt))]")
        .Define("SortedJet_pt",  "FilteredJet_pt[Reverse(Argsort(FilteredJet_pt))]")

        .Define("SortedJet_smearUp_pt", "FilteredJet_smearUp_pt[Reverse(Argsort(FilteredJet_pt))]")
        .Define("SortedJet_smearDn_pt", "FilteredJet_smearDn_pt[Reverse(Argsort(FilteredJet_pt))]")
        .Define("SortedJet_scaleUp_pt", "FilteredJet_scaleUp_pt[Reverse(Argsort(FilteredJet_pt))]")
        .Define("SortedJet_scaleDn_pt", "FilteredJet_scaleDn_pt[Reverse(Argsort(FilteredJet_pt))]");

    auto skim = sorted
        .Filter("nZCand >= 0 && SortedJet_eta.size() > 0 && Flag_JetVetoed == 0")
        .Define("weight", Form("1000 * overallEventWeight * %f / %f", _lumi, gen_sumWeights))
        .Define("leadEta",  "SortedJet_eta.at(0)")
        .Define("leadPt",   "SortedJet_pt.at(0)")
        .Define("leadEtaUp",   computeLeadJet, {"SortedJet_eta", "SortedJet_smearUp_pt"})
        .Define("leadEtaDown", computeLeadJet, {"SortedJet_eta", "SortedJet_smearDn_pt"})
        .Define("leadEtaDownscale", computeLeadJet, {"SortedJet_eta", "SortedJet_scaleDn_pt"})
        .Define("leadEtaUpscale", computeLeadJet, {"SortedJet_eta", "SortedJet_scaleUp_pt"})
        .Define("leadPtUp",    computeLeadJet, {"SortedJet_pt", "SortedJet_smearUp_pt"})
        .Define("leadPtDown",  computeLeadJet, {"SortedJet_pt", "SortedJet_smearDn_pt"})
        .Define("leadPtUpscale", computeLeadJet, {"SortedJet_pt", "SortedJet_scaleUp_pt"})
        .Define("leadPtDownscale", computeLeadJet, {"SortedJet_pt", "SortedJet_scaleDn_pt"});

    auto h_nom = skim.Histo1D({(var+"_MC").c_str(), (var+"_MC").c_str(), bins, xmin, xmax}, var, "weight");
    h_nom->Sumw2();
    auto h_up  = skim.Histo1D({(var+"Up_MC").c_str(), (var+"Up_MC").c_str(), bins, xmin, xmax}, var+"Up", "weight");
    h_up->Sumw2();
    auto h_dn  = skim.Histo1D({(var+"Down_MC").c_str(), (var+"Down_MC").c_str(), bins, xmin, xmax}, var+"Down", "weight");
    h_dn->Sumw2();
    auto hist_up_scale = skim.Histo1D({(var + "Up_MC_scale").c_str(), (var + "Up_MC_scale").c_str(), bins, xmin, xmax}, var + "Upscale", "weight");
    hist_up_scale->Sumw2();
    auto hist_dn_scale = skim.Histo1D({(var + "Down_MC_scale").c_str(), (var + "Down_MC_scale").c_str(), bins, xmin, xmax}, var + "Downscale", "weight"); 
    hist_dn_scale->Sumw2();
    return std::make_tuple(h_nom, h_up, h_dn, hist_up_scale, hist_dn_scale);
}

// -----------------------------------------------------------------------------
// DATA HISTO (multiple files)
// -----------------------------------------------------------------------------
ROOT::RDF::RResultPtr<TH1D>
Histo_Data(const std::vector<std::string>& files,
           std::string var, int bins, double xmin, double xmax)
{

    ROOT::EnableImplicitMT();
    ROOT::RDF::RSnapshotOptions opts;
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
void Jets_2024_newRatio()
{
    gEnv->SetValue("TTreeProcessorMT.SkipBadFiles", 1);
    gEnv->SetValue("TTreeProcessorMT.MaxBadFiles", -1); // no limit
    ROOT::EnableImplicitMT();
    //ROOT::Internal::RDF::SetSkipBadFiles(true);

    std::string var = "leadEta";
    std::string var_name = "Leading #eta";
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

    auto DY_files = collect_root_files(
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/270126/MC/2024/DYJetsToLL"
    );

    std::cout << "[INFO] DY files loaded: " << DY_files.size() << std::endl;


    auto TT_files = collect_root_files(
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/270126/MC/2024/TTto2L2Nu"
    );

    std::cout << "[INFO] TT files loaded: " << TT_files.size() << std::endl;

    auto data_files = collect_root_files(
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/270126/Data/2024"
    );

    std::cout << "[INFO] Data files loaded: " << data_files.size() << std::endl;



    // Filter out corrupted files
    auto filterGood = [](const std::vector<std::string>& files){
        std::vector<std::string> good;
        for (auto &f : files) {
            TFile *file = TFile::Open(f.c_str());
            if (file && !file->IsZombie()) good.push_back(f);
            if(file) file->Close();
        }
        return good;
    };

    // Apply filter
    DY_files   = filterGood(DY_files);
    TT_files   = filterGood(TT_files);
    data_files = filterGood(data_files);

    std::cout << "[INFO] DY files loaded: " << DY_files.size() << std::endl;
    std::cout << "[INFO] TT files loaded: " << TT_files.size() << std::endl;
    std::cout << "[INFO] Data files loaded: " << data_files.size() << std::endl;

    double lumi = 108.8;

    // MC
    auto [h_DY, h_DY_up, h_DY_dn, h_DY_up_scale, h_DY_dn_scale] = Histo_MC(DY_files, lumi, var, n_bins, xmin, xmax);
    auto [h_TT, h_TT_up, h_TT_dn, h_TT_up_scale, h_TT_dn_scale] = Histo_MC(TT_files, lumi, var, n_bins, xmin, xmax);

    // Combine DY + TT
    h_DY->Scale(3);
    h_DY->Add(h_TT.GetPtr());
    h_DY_up->Scale(3);
    h_DY_dn->Scale(3);
    h_DY_up_scale->Scale(3);
    h_DY_dn_scale->Scale(3);

    h_DY_up->Add(h_TT_up.GetPtr());
    h_DY_dn->Add(h_TT_dn.GetPtr());
    h_DY_up_scale -> Add(h_TT_up_scale.GetPtr());
    h_DY_dn_scale -> Add(h_TT_dn_scale.GetPtr());

    // DATA
    auto h_data = Histo_Data(data_files, var, n_bins, xmin, xmax);

    // Save ROOT file
    std::filesystem::create_directories("2024_newRatio/");
    TFile outfile(("2024_newRatio/" + var + "_histos.root").c_str(),"RECREATE");
    h_data->Write();
    h_DY->Write();
    h_DY_up->Write();
    h_DY_dn->Write();
    h_DY_up_scale -> Write();
    h_DY_dn_scale -> Write();
    outfile.Close();

    //h_DY_up.GetPtr(),
    //h_DY_dn.GetPtr(),

    // Draw plots
    TCanvas *c = new TCanvas();
    DrawRatioPlot("2024_newRatio/" + var + "_plot", var_name, c,
                  h_data.GetPtr(),
                  h_DY.GetPtr(),
                  h_DY_up.GetPtr(),
                  h_DY_dn.GetPtr(), 
                  h_DY_up_scale.GetPtr(), 
                  h_DY_dn_scale.GetPtr(),
                  lumi);
}