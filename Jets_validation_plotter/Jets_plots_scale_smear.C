//Martina Manoni plotter for Jet control plots
#include "CMS_lumi.C"
#include <tuple>
#include <vector>
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>

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

            if(mcVal == 0){
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

        /*
        double temp;
        for(int bin = 0; bin < Ratio->GetSize() - 2; bin++){
          temp = Ratio->GetBinContent(bin + 1);
          Ratio->SetBinContent( bin + 1, temp);
        }

        for(int bin = 0; bin < sigma_up->GetSize() - 2; bin++){
          temp = sigma_up->GetBinContent(bin + 1);
          sigma_up->SetBinContent( bin + 1, temp);
        }

        for(int bin = 0; bin < sigma_dn->GetSize() - 2; bin++){
          temp = sigma_dn->GetBinContent(bin + 1);
          sigma_dn->SetBinContent( bin + 1, temp);
        }
        */

        Ratio->SetTitle("");
        sigma_up->SetTitle("");
        sigma_dn->SetTitle("");

        sigma_up->SetTitle("");
        sigma_up->SetFillColor(kGray);
        sigma_up->SetLineColor(kGray);
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

//This function calculates the eta of the leading jet after applying a downward or upward variation to the jet energy
double computeLeadJet(ROOT::RVec<float> &_JetEta, ROOT::RVec<float> &_JetPt_JER){
  double highest = 0;
  int i_up = 0;
  for(int i=0; i<_JetEta.size(); i++){
    if (_JetPt_JER.at(i) > highest && _JetPt_JER.at(i) > 30){
            highest = _JetPt_JER.at(i);
            i_up = i;
    }
  }
  return _JetEta.at(i_up);
}


// Function to get genEventSumw using RDataFrame
double get_genEventSumw(const std::string& fpath) {
    ROOT::RDataFrame rdf_runs("Runs", fpath); // Define RDataFrame for the "Runs" tree
    auto gen_sumWeights = *(rdf_runs.Sum("genEventSumw")); // Sum over genEventSumw column
    //std::cout << "[INFO] genEventSumw for file " << fpath << ": " << gen_sumWeights << std::endl;
    return gen_sumWeights;
}


// Histo_MC
tuple<ROOT::RDF::RResultPtr<TH1D>, ROOT::RDF::RResultPtr<TH1D>, ROOT::RDF::RResultPtr<TH1D>, ROOT::RDF::RResultPtr<TH1D>, ROOT::RDF::RResultPtr<TH1D>> Histo_MC(
    string fpath, double _lumi, string var , int bins, double xmin, double xmax) {
    double gen_sumWeights = get_genEventSumw(fpath); // Get the sum of genEventSumw from Runs tree
    std::cout << "[INFO] genEventSumw for file " << fpath << ": " << gen_sumWeights << std::endl;

    ROOT::RDataFrame rdf("Events", fpath);
    // Step 1: Apply Jets filters
    auto sorted_rdf = rdf
        .Define("Horn_Cut", 
        "((abs(Jet_eta) < 2.5 || abs(Jet_eta) > 3) || Jet_pt > 50)")

        .Define(
        "HF_Cut",
        "(abs(Jet_eta) < 3 || abs(Jet_eta) > 5 || Jet_pt > 50)")
    
        .Define("Jet_Mask", 
        "(Jet_pt > 30) && Horn_Cut && HF_Cut && Jet_jetId==6 && Jet_ZZMask == false")
        
        .Define("Jet_Mask_scaleUp", 
        "(Jet_scaleUp_pt > 30) && Horn_Cut && HF_Cut && Jet_jetId==6 && Jet_ZZMask == false")
        
        .Define("Jet_Mask_scaleDown", 
        "(Jet_scaleDn_pt > 30) && Horn_Cut && HF_Cut && Jet_jetId==6 && Jet_ZZMask == false") 
        
        .Define("Jet_Mask_smearUp", 
        "(Jet_smearUp_pt > 30) && Horn_Cut && HF_Cut && Jet_jetId==6 && Jet_ZZMask == false") 
        
        .Define("Jet_Mask_smearDown", 
        "(Jet_smearDn_pt > 30) && Horn_Cut && HF_Cut && Jet_jetId==6 && Jet_ZZMask == false") 
        
        .Define("FilteredJet_eta", "Jet_eta[Jet_Mask]")
        .Define("FilteredJet_pt", "Jet_pt[Jet_Mask]")
        .Define("FilteredJet_Njets", "FilteredJet_pt.size()")

        .Define("FilteredJet_scaleUp_Njets", "Jet_pt[Jet_Mask_scaleUp].size()")
        .Define("FilteredJet_scaleDn_Njets", "Jet_pt[Jet_Mask_scaleDown].size()")
        .Define("FilteredJet_smearUp_Njets", "Jet_pt[Jet_Mask_smearUp].size()")
        .Define("FilteredJet_smearDn_Njets", "Jet_pt[Jet_Mask_smearDown].size()")

        .Define("FilteredJet_smearUp_pt", "Jet_smearUp_pt[Jet_Mask]")
        .Define("FilteredJet_smearDn_pt", "Jet_smearDn_pt[Jet_Mask]")
        .Define("FilteredJet_scaleUp_pt", "Jet_scaleUp_pt[Jet_Mask]")
        .Define("FilteredJet_scaleDn_pt", "Jet_scaleDn_pt[Jet_Mask]")

        .Define("SortedJet_eta", "FilteredJet_eta[Reverse(Argsort(FilteredJet_pt))]")
        .Define("SortedJet_pt", "FilteredJet_pt[Reverse(Argsort(FilteredJet_pt))]")
        
        .Define("SortedJet_smearUp_pt", "FilteredJet_smearUp_pt[Reverse(Argsort(FilteredJet_pt))]")
        .Define("SortedJet_smearDn_pt", "FilteredJet_smearDn_pt[Reverse(Argsort(FilteredJet_pt))]")
        .Define("SortedJet_scaleUp_pt", "FilteredJet_scaleUp_pt[Reverse(Argsort(FilteredJet_pt))]")
        .Define("SortedJet_scaleDn_pt", "FilteredJet_scaleDn_pt[Reverse(Argsort(FilteredJet_pt))]");
        
    // Step 2: Apply the event-level filter
    auto skim_rdf = sorted_rdf.Filter("nZCand >= 0 && SortedJet_eta.size() > 0 && Flag_JetVetoed == 0") //Z control region used (nZCand)
        .Define("weight", Form("1000 * overallEventWeight * %f / %f", _lumi, gen_sumWeights))
        .Define("leadEta", "SortedJet_eta.at(0)")
        .Define("leadPt", "SortedJet_pt.at(0)")
        .Define("leadNjets", "FilteredJet_Njets")

        .Define("leadEtaUp", computeLeadJet, {"SortedJet_eta", "SortedJet_smearUp_pt"})
        .Define("leadEtaDown", computeLeadJet, {"SortedJet_eta", "SortedJet_smearDn_pt"})
        .Define("leadEtaDownscale", computeLeadJet, {"SortedJet_eta", "SortedJet_scaleDn_pt"})
        .Define("leadEtaUpscale", computeLeadJet, {"SortedJet_eta", "SortedJet_scaleUp_pt"})


        .Define("leadPtUp", computeLeadJet, {"SortedJet_pt", "SortedJet_smearUp_pt"})
        .Define("leadPtDown", computeLeadJet, {"SortedJet_pt", "SortedJet_smearDn_pt"})
        .Define("leadPtUpscale", computeLeadJet, {"SortedJet_pt", "SortedJet_scaleUp_pt"})
        .Define("leadPtDownscale", computeLeadJet, {"SortedJet_pt", "SortedJet_scaleDn_pt"})

        .Define("leadNjetsUp", "FilteredJet_smearUp_Njets")
        .Define("leadNjetsDown", "FilteredJet_smearDn_Njets" )
        .Define("leadNjetsUpscale", "FilteredJet_scaleUp_Njets")
        .Define("leadNjetsDownscale", "FilteredJet_scaleDn_Njets");

        std::cout << " HISTO MC: create histograms "<<std::endl;
        auto hist_nominal = skim_rdf.Histo1D({(var + "_MC").c_str(), (var + "_MC").c_str(), bins, xmin, xmax}, var, "weight");
        auto hist_up = skim_rdf.Histo1D({(var + "Up_MC").c_str(), (var + "Up_MC").c_str(), bins, xmin, xmax}, var + "Up", "weight");
        auto hist_dn = skim_rdf.Histo1D({(var + "Down_MC").c_str(), (var + "Down_MC").c_str(), bins, xmin, xmax}, var + "Down", "weight"); 
        auto hist_up_scale = skim_rdf.Histo1D({(var + "Up_MC_scale").c_str(), (var + "Up_MC_scale").c_str(), bins, xmin, xmax}, var + "Upscale", "weight");
        auto hist_dn_scale = skim_rdf.Histo1D({(var + "Down_MC_scale").c_str(), (var + "Down_MC_scale").c_str(), bins, xmin, xmax}, var + "Downscale", "weight"); 

        std::cout << " HISTO MC: finished creating histograms"<<std::endl;

        return std::make_tuple(hist_nominal, hist_up, hist_dn, hist_up_scale, hist_dn_scale);
}


// Function to create histograms for data samples
ROOT::RDF::RResultPtr<TH1D> Histo_Data(string fpath, string var, int bins, double xmin, double xmax) {
    ROOT::RDataFrame rdf("Events", fpath);

    auto sorted_data = rdf
        .Define("Horn_Cut", 
        "((abs(Jet_eta) < 2.5 || abs(Jet_eta) > 3) || Jet_pt > 50)") 
        .Define(
        "HF_Cut",
        "(abs(Jet_eta) < 3 || abs(Jet_eta) > 5 || Jet_pt > 50)")
        .Define("Jet_Mask", 
        "(Jet_pt > 30) && Horn_Cut && HF_Cut && Jet_jetId==6 && Jet_ZZMask == false")
        .Define("FilteredJet_eta", "Jet_eta[Jet_Mask]")
        .Define("FilteredJet_pt", "Jet_pt[Jet_Mask]")
        .Define("FilteredJet_Njets", "FilteredJet_pt.size()")

        .Define("SortedJet_eta", "FilteredJet_eta[Reverse(Argsort(FilteredJet_pt))]")
        .Define("SortedJet_pt", "FilteredJet_pt[Reverse(Argsort(FilteredJet_pt))]");

    // Apply event-level filter and define leadEta and leadPt
    auto skim_data = sorted_data.Filter("nZCand >= 0 && SortedJet_eta.size() > 0 && Flag_JetVetoed == 0") //ZLCand_lepIdx nZCand
        .Define("leadEta", "SortedJet_eta.at(0)")
        .Define("leadPt", "SortedJet_pt.at(0)")
        .Define("leadNjets", "FilteredJet_Njets");

    // Create histogram for the specified variable
    std::cout << "HISTO DATA: finished creating histograms "<<std::endl;
    auto hist_data = skim_data.Histo1D({(var + "_data").c_str(), (var + "_data").c_str(), bins, xmin, xmax}, var);

    return hist_data;
}


void Jets_plots_scale_smear(){
  std::cout<<"Jets_plots_scale_smear"<<std::endl;
  ROOT::EnableImplicitMT();
  std::string var = "leadEta"; //"leadNjets";//"leadEta";//leadPt //leadchEmEF //leadEta //leadneEmEF //leadNJets
  std::string var_name = "#eta"; //"p_{T}";//"#eta"; ////Neutral Em Energy Fraction //charged Em Energy Fraction //Number of jets
  std::string period = "postBPix2023"; //"preEE2022" "postBPix2023" "preBPix2023"

  // Define a map that associates 'var' with the corresponding binning
  std::map<std::string, std::tuple<int, double, double>> var_binning;
  
  // binning configurations
  var_binning["leadEta"] = std::make_tuple(47, -4.7, +4.7);
  var_binning["leadPt"] = std::make_tuple(40, 30, 200.0);
  var_binning["leadNjets"] = std::make_tuple(8, 2, 10);

  int n_bins;
  double min_val, max_val;
  std::tie(n_bins, min_val, max_val) = var_binning[var];

  /*std::cout<<"period preEE2022"<<std::endl;
  std::string fDY   = "/eos/user/m/mmanoni/HZZ_prod_241125_NoEleMuo_YesJetCorr/2022_MC/PROD_samplesNano_2022_MC_cc84ce40/DYJetsToLL/ZZ4lAnalysis.root";
  std::string fTT   = "/eos/user/m/mmanoni/HZZ_prod_241125_NoEleMuo_YesJetCorr/2022_MC/PROD_samplesNano_2022_MC_cc84ce40/TTto2L2Nu/ZZ4lAnalysis.root";
  std::string fdata = "/eos/user/m/mmanoni/HZZ_prod_241125_NoEleMuo_YesJetCorr/2022_Data/PROD_samplesNano_2022_Data_cc84ce40/Data_eraCD_preEE.root";
  double lumi  = 7.98;*/

  /*std::cout<<"period postEE2022"<<std::endl;
  string fDY   = "/eos/user/m/mmanoni/HZZ_prod_241125_NoEleMuo_YesJetCorr/2022EE_MC/PROD_samplesNano_2022EE_MC_cc84ce40/DYJetsToLL/ZZ4lAnalysis.root";
  string fTT   = "/eos/user/m/mmanoni/HZZ_prod_241125_NoEleMuo_YesJetCorr/2022EE_MC/PROD_samplesNano_2022EE_MC_cc84ce40/TTto2L2Nu/ZZ4lAnalysis.root"; 
  string fdata = "/eos/user/m/mmanoni/HZZ_prod_241125_NoEleMuo_YesJetCorr/2022_Data/PROD_samplesNano_2022_Data_cc84ce40/Data_eraEFG_postEE.root";
  double lumi  = 26.67;*/

  /*std::cout<<"period preBPix2023"<<std::endl;
  string fDY   = "/eos/user/m/mmanoni/HZZ_prod_241125_NoEleMuo_YesJetCorr/2023preBPix_MC/PROD_samplesNano_2023preBPix_MC_cc84ce40/DYJetsToLL/ZZ4lAnalysis.root";
  string fTT   = "/eos/user/m/mmanoni/HZZ_prod_241125_NoEleMuo_YesJetCorr/2023preBPix_MC/PROD_samplesNano_2023preBPix_MC_cc84ce40/TTto2L2Nu/ZZ4lAnalysis.root";
  string fdata = "/eos/user/m/mmanoni/HZZ_prod_241125_NoEleMuo_YesJetCorr/2023_Data/PROD_samplesNano_2023_Data_cc84ce40/Data_eraC_preBPix.root";
  double lumi  = 18.06;*/

  std::cout<<"period postBPix2023"<<std::endl;
  string fDY   = "/eos/user/m/mmanoni/HZZ_prod_241125_NoEleMuo_YesJetCorr/2023postBPix_MC/PROD_samplesNano_2023postBPix_MC_cc84ce40/DYJetsToLL/ZZ4lAnalysis.root";
  string fTT   = "/eos/user/m/mmanoni/HZZ_prod_241125_NoEleMuo_YesJetCorr/2023postBPix_MC/PROD_samplesNano_2023postBPix_MC_cc84ce40/TTto2L2Nu/ZZ4lAnalysis.root";
  string fdata = "/eos/user/m/mmanoni/HZZ_prod_241125_NoEleMuo_YesJetCorr/2023_Data/PROD_samplesNano_2023_Data_cc84ce40/Data_eraD_postBPix.root";
  double lumi  = 9.69;


  ROOT::RDF::RResultPtr<TH1D> hist_DY;
  ROOT::RDF::RResultPtr<TH1D> hist_DY_up;
  ROOT::RDF::RResultPtr<TH1D> hist_DY_dn;
  ROOT::RDF::RResultPtr<TH1D> hist_DY_up_scale;
  ROOT::RDF::RResultPtr<TH1D> hist_DY_dn_scale;
  tie(hist_DY,hist_DY_up,hist_DY_dn, hist_DY_up_scale, hist_DY_dn_scale) = Histo_MC(fDY, lumi, var, n_bins, min_val, max_val); //Histo_MC is expected to return a tuple of three values, Each of these variables (hist_DY, hist_DY_up, hist_DY_dn) will be assigned one element from the tuple returned by Histo_MC

  ROOT::RDF::RResultPtr<TH1D> hist_ttbar;
  ROOT::RDF::RResultPtr<TH1D> hist_ttbar_up;
  ROOT::RDF::RResultPtr<TH1D> hist_ttbar_dn;
  ROOT::RDF::RResultPtr<TH1D> hist_ttbar_up_scale;
  ROOT::RDF::RResultPtr<TH1D> hist_ttbar_dn_scale;
  tie(hist_ttbar,hist_ttbar_up,hist_ttbar_dn,hist_ttbar_up_scale,hist_ttbar_dn_scale) = Histo_MC(fTT, lumi, var, n_bins, min_val, max_val);

  // Process data histogram
  ROOT::RDF::RResultPtr<TH1D> hist_data = Histo_Data(fdata, var, n_bins, min_val, max_val);

  hist_DY    -> Add(hist_ttbar.GetPtr());
  hist_DY_up -> Add(hist_ttbar_up.GetPtr());
  hist_DY_dn -> Add(hist_ttbar_dn.GetPtr());
  hist_DY_up_scale -> Add(hist_ttbar_up_scale.GetPtr());
  hist_DY_dn_scale -> Add(hist_ttbar_dn_scale.GetPtr());

  //std::cout << "Saving distributions into root file ..." << std::endl;
  //TFile * outfile = new TFile("{var}_histos.root", "RECREATE");

  std::string subdir = "2023_postBPix_HF_final/"; 
  std::filesystem::create_directories(subdir);


  std::cout << "Saving distributions into root file ..." << std::endl;
  std::string filename = subdir + var +"_"+ period + "_scale_smear_Norm.root";  // Concatenate the value of 'var' with the string
  TFile *outfile = new TFile(filename.c_str(), "RECREATE");

  outfile    -> cd();
  hist_data  -> Write();
  hist_DY    -> Write();
  hist_DY_up -> Write();
  hist_DY_dn -> Write();
  hist_DY_up_scale -> Write();
  hist_DY_dn_scale -> Write();
  outfile    -> Close();
  

  TCanvas * canvas = new TCanvas("canvas","canvas",1200,900);
  std::string filename_plots = subdir + var +"_"+ period + "_plots_scale_smear_JetID6";
  DrawRatioPlot(filename_plots, var_name, canvas,
                hist_data.GetPtr(),
                hist_DY.GetPtr(), hist_DY_up.GetPtr(), hist_DY_dn.GetPtr(), hist_DY_up_scale.GetPtr(), hist_DY_dn_scale.GetPtr(),
                lumi);

}