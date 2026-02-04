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
void DrawRatioPlot(string name, string var_name, TCanvas *c, TH1D *data, TH1D *MC, TH1D*MCUp, TH1D*MCDn, TH1D*MCUpscale, TH1D*MCDnscale, double _lumi){
        std::cout << "[INFO] Starting DrawRatioPlot function for plot: " << name << std::endl;
        std::cout << "[INFO] Lumi: " << _lumi << std::endl;

        gStyle->SetOptStat(0);

        std::cout << "Scale and smearing calculation start "<<std::endl;
        // Sum uncertainties in quadrature
        TH1D *MCUpFinal = new TH1D("MCUpFinal", "MCUpFinal", MC->GetSize() - 2, MC->GetXaxis()->GetXmin(), MC->GetXaxis()->GetXmax());
        TH1D *MCDnFinal = new TH1D("MCDnFinal", "MCDnFinal", MC->GetSize() - 2, MC->GetXaxis()->GetXmin(), MC->GetXaxis()->GetXmax());

        for (int bin = 0; bin < MC->GetSize() - 2; bin++) {
          // --- Nominal ---
          double nom = MC->GetBinContent(bin + 1);

          // --- Compute shifts (signed deviations) ---
          double dJECup     = MCUp->GetBinContent(bin + 1)     - nom;
          double dJECdn     = MCDn->GetBinContent(bin + 1)     - nom;
          double dScaleUp   = MCUpscale->GetBinContent(bin + 1) - nom;
          double dScaleDn   = MCDnscale->GetBinContent(bin + 1) - nom;

          // --- Quadrature accumulators ---
          double up2 = 0.0;
          double dn2 = 0.0;

          // Helper to accumulate by sign
          auto accumulate = [&](double d) {
              if (d > 0) up2 += d*d;
              else       dn2 += d*d;
          };

          // --- Add *all* signed variations ---
          accumulate(dJECup);
          accumulate(dJECdn);
          accumulate(dScaleUp);
          accumulate(dScaleDn);

          // --- Include MC statistical uncertainty ---
          double mc_stat = MC->GetBinError(bin + 1);
          up2 += mc_stat * mc_stat;
          dn2 += mc_stat * mc_stat;

          // --- Final uncertainties ---
          double total_up = sqrt(up2);
          double total_dn = sqrt(dn2);

          // --- Set final MC uncertainty envelopes ---
          MCUpFinal->SetBinContent(bin + 1, nom + total_up);
          MCDnFinal->SetBinContent(bin + 1, nom - total_dn);

          // --- Print debug ---
          std::cout << "Bin " << bin+1
                    << ": nom=" << nom
                    << "  dJECup="   << dJECup
                    << "  dJECdn="   << dJECdn
                    << "  dScaleUp=" << dScaleUp
                    << "  dScaleDn=" << dScaleDn
                    << "  total_up=" << total_up
                    << "  total_dn=" << total_dn
                    << std::endl;
        }

        // Normalize MC to data
        double dataIntegral = data->Integral();
        double mcIntegral = MC->Integral();
        double scaleFactor = dataIntegral / mcIntegral;
        //double scaleFactor = 1.0;
        MC->Scale(scaleFactor);
        MCUp->Scale(scaleFactor);
        MCDn->Scale(scaleFactor);
        MCUpscale->Scale(scaleFactor);
        MCDnscale->Scale(scaleFactor);
        MCUpFinal->Scale(scaleFactor);
        MCDnFinal->Scale(scaleFactor);

        //Splits the canvas into two vertical pads
        c->Divide(0,2,0,0);
        // TPad * pad1 = new TPad("pad1","pad1", 0, 0.3, 1, 1.0);
        // pad1->Draw();
        c->cd(1);


        gPad->SetBottomMargin(0.02);
        gPad->SetTopMargin(0.18);
        gPad->SetLeftMargin(0.10);
        
        //***Main Histogram Plot (Top Panel)***


        //MC->SetTitle("");
        //data->SetTitle("");
        //MCUpFinal->SetTitle("");
        //MCDnFinal->SetTitle("");

        double max = 0.;
        for (int bin = 0; bin < MC->GetSize() - 2; bin++){
          if ( MC->GetBinContent(bin + 1) > max) max = MC->GetBinContent(bin + 1);
        }
        MC->SetMaximum(1.4*max); //Scales the y-axis maximum to 140% of the largest bin for better visibility.

        
        // CMS-TDR inspired colors
        MC->SetFillColor(TColor::GetColor("#6baed6"));   // soft blue
        MC->SetLineColor(TColor::GetColor("#2171b5"));   // darker blue outline
        MC->SetLineWidth(2);
        double binWidth = MC->GetXaxis()->GetBinWidth(1); 
        TString yAxisLabel = TString::Format("Events/%.2f", binWidth);
        MC->GetYaxis()->SetTitle(yAxisLabel);
        MC->Draw("HIST");
        MC->GetXaxis()->SetLabelSize(0);
        MC->GetYaxis()->SetTitleSize(0.07);
        MC->GetYaxis()->SetLabelSize(0.07);
        MC->GetYaxis()->SetTitleOffset(10.);

        MC->GetYaxis()->SetNoExponent(kFALSE);  // allow 10^x
        MC->GetYaxis()->SetMoreLogLabels(kTRUE);
        MC->GetYaxis()->SetMaxDigits(2);

        MC->GetXaxis()->SetLabelSize(0); // hide x labels on top


        std::cout << "Scale and smearing calculation end "<<std::endl;
        data->SetMarkerStyle(20);
        data->SetMarkerColor(kBlack);
        data->SetLineColor(kBlack);
        data->SetMarkerSize(0.8);
        data->Draw("P SAME");


        // ---- CMS label ----
        CMS_lumi cms;
        cms.set_lumi((TPad*)gPad, _lumi);


        MC->SetTitle("");
        data->SetTitle("");
        MCUpFinal->SetTitle("");
        MCDnFinal->SetTitle("");
        
        //***Ratio Plot (Bottom Panel)***
        c->cd(2);

        gPad->SetBottomMargin(0.20);
        gPad->SetTopMargin(0.02);
        gPad->SetLeftMargin(0.10);

        TH1F * Ratio = new TH1F("Ratio","Ratio", MC->GetSize() - 2, MC->GetXaxis()->GetXmin(), MC->GetXaxis()->GetXmax());
        TH1F * sigma_up = new TH1F("sigma_up","sigma_up", MC->GetSize() - 2, MC->GetXaxis()->GetXmin(), MC->GetXaxis()->GetXmax());
        TH1F * sigma_dn = new TH1F("sigma_dn","sigma_dn", MC->GetSize() - 2, MC->GetXaxis()->GetXmin(), MC->GetXaxis()->GetXmax());

        for(int bin = 1; bin <= MC->GetNbinsX(); bin++){
            double mcVal   = MC->GetBinContent(bin);
            double dataVal = data->GetBinContent(bin);
            double dataErr = data->GetBinError(bin);

            if(mcVal < 1e-6){
                Ratio->SetBinContent(bin, 0);
                Ratio->SetBinError(bin, 0);
                sigma_up->SetBinContent(bin, 0);
                sigma_dn->SetBinContent(bin, 0);
            } else {
                // Fractional ratio and data error
                Ratio->SetBinContent(bin, dataVal / mcVal - 1);
                Ratio->SetBinError(bin, 0);

                // Fractional MC systematic uncertainties
                sigma_up->SetBinContent(bin, MCUpFinal->GetBinContent(bin) / mcVal - 1);
                sigma_dn->SetBinContent(bin, MCDnFinal->GetBinContent(bin) / mcVal - 1);
            }
        }

        std::cout << "[INFO] Fractional MC uncertainties (sigma_up, sigma_dn):" << std::endl;
        for(int bin = 1; bin <= MC->GetNbinsX(); bin++){
            double sigmaUpVal = sigma_up->GetBinContent(bin);
            double sigmaDnVal = sigma_dn->GetBinContent(bin);
            std::cout << "Bin " << bin 
                      << ": sigma_up = " << sigmaUpVal 
                      << ", sigma_dn = " << sigmaDnVal << std::endl;
        }

        Ratio->SetTitle("");
        sigma_up->SetTitle("");
        sigma_dn->SetTitle("");

        sigma_up->SetTitle("");
        sigma_up->SetFillColor(kGray+2);
        sigma_up->SetLineColor(kGray+2);
        sigma_up->SetMaximum(1.0);
        sigma_up->SetMinimum(-1.0);
        std::string ytitle_plots = "Leading Jet " + var_name;
        sigma_up->GetXaxis()->SetTitle(ytitle_plots.c_str());
        sigma_up->GetYaxis()->SetTitle("Data/MC");
        sigma_up->GetYaxis()->SetTitleOffset(10.);
        sigma_up->Draw("HIST");
        sigma_up->GetXaxis()->SetLabelSize(0.07);
        sigma_up->GetXaxis()->SetTitleSize(0.07);
        sigma_up->GetYaxis()->SetLabelSize(0.07);
        sigma_up->GetYaxis()->SetTitleSize(0.07);

        sigma_dn->SetFillColor(kGray);
        sigma_dn->SetLineColor(kGray);
        sigma_dn->Draw("HIST SAME");

        gPad->RedrawAxis();

        Ratio->SetMarkerStyle(20);
        Ratio->SetMarkerSize(0.8);
        Ratio->Draw("P SAME");


        double x_min = Ratio->GetXaxis()->GetXmin();
        double x_max = Ratio->GetXaxis()->GetXmax();
        TLine *line0 = new TLine(x_min, 0, x_max, 0);
        line0->SetLineColor(kBlack);  // or kBlack
        line0->SetLineStyle(2);     // dashed line
        line0->SetLineWidth(2);
        line0->Draw("SAME");

        c->cd(1);
        TLegend * leg = new TLegend(0.64,0.45,0.85,0.75);
        leg->SetBorderSize(0);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.07);
        leg->AddEntry(MC, "DY + t#bar{t}  MC", "f" );
        leg->AddEntry(data, "Data", "p");
        //leg->AddEntry(UncertaintyBand, "Uncertainty Band", "f");
        leg->AddEntry(sigma_up, "Stat.+Syst Unc.", "f");
        leg->Draw();

        c->Update();

        c->SaveAs((TString)name+".pdf");
        c->SaveAs((TString)name+".png");
        c->SaveAs((TString)name+".C");

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
    auto h_up  = skim.Histo1D({(var+"Up_MC").c_str(), (var+"Up_MC").c_str(), bins, xmin, xmax}, var+"Up", "weight");
    auto h_dn  = skim.Histo1D({(var+"Down_MC").c_str(), (var+"Down_MC").c_str(), bins, xmin, xmax}, var+"Down", "weight");
    auto hist_up_scale = skim.Histo1D({(var + "Up_MC_scale").c_str(), (var + "Up_MC_scale").c_str(), bins, xmin, xmax}, var + "Upscale", "weight");
    auto hist_dn_scale = skim.Histo1D({(var + "Down_MC_scale").c_str(), (var + "Down_MC_scale").c_str(), bins, xmin, xmax}, var + "Downscale", "weight"); 

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