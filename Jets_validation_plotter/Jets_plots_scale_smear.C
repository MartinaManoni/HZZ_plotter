//Martina Manoni plotter for Jet control plots
#include "CMS_lumi.C"
#include <tuple>
#include <vector>
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>

void DrawRatioPlot(string name, string var_name, TCanvas *c, TH1D *data, TH1D *MC, TH1D*MCUp, TH1D*MCDn, TH1D*MCUpscale, TH1D*MCDnscale, double _lumi){
        std::cout << "[INFO] Starting DrawRatioPlot function for plot: " << name << std::endl;
        std::cout << "[INFO] Lumi: " << _lumi << std::endl;

        // Normalize MC to data
        double dataIntegral = data->Integral();
        double mcIntegral = MC->Integral();
        double scaleFactor = dataIntegral / mcIntegral;
        MC->Scale(scaleFactor);
        MCUp->Scale(scaleFactor);
        MCDn->Scale(scaleFactor);
        MCUpscale->Scale(scaleFactor);
        MCDnscale->Scale(scaleFactor);

        //Splits the canvas into two vertical pads
        c->Divide(0,2,0,0);
        // TPad * pad1 = new TPad("pad1","pad1", 0, 0.3, 1, 1.0);
        // pad1->Draw();
        c->cd(1);

        gPad->SetBottomMargin(0.02);
        gPad->SetTopMargin(0.18);
        gPad->SetLeftMargin(0.10);
        
        //***Main Histogram Plot (Top Panel)***
        double max = 0.;
        for (int bin = 0; bin < MC->GetSize() - 2; bin++){
          if ( MC->GetBinContent(bin + 1) > max) max = MC->GetBinContent(bin + 1);
        }
        MC->SetMaximum(1.4*max); //Scales the y-axis maximum to 140% of the largest bin for better visibility.

        MC->SetFillColor(kOrange + 1);
        double binWidth = MC->GetXaxis()->GetBinWidth(1); 
        TString yAxisLabel = TString::Format("Events/%.2f", binWidth);
        MC->GetYaxis()->SetTitle(yAxisLabel);
        MC->Draw("HIST");
        MC->GetXaxis()->SetLabelSize(0);
        MC->GetYaxis()->SetTitleSize(0.07);
        MC->GetYaxis()->SetLabelSize(0.07);
        MC->GetYaxis()->SetTitleOffset(0.7);

        std::cout << "Scale and smearing calculation start "<<std::endl;
        // Sum uncertainties in quadrature
        TH1D *MCUpFinal = new TH1D("MCUpFinal", "MCUpFinal", MC->GetSize() - 2, MC->GetXaxis()->GetXmin(), MC->GetXaxis()->GetXmax());
        TH1D *MCDnFinal = new TH1D("MCDnFinal", "MCDnFinal", MC->GetSize() - 2, MC->GetXaxis()->GetXmin(), MC->GetXaxis()->GetXmax());

        for (int bin = 0; bin < MC->GetSize() - 2; bin++) {
            double up_var = 0.0, dn_var = 0.0;

            double var1 = MCUp->GetBinContent(bin + 1) - MC->GetBinContent(bin + 1);
            double var2 = MCUpscale->GetBinContent(bin + 1) - MC->GetBinContent(bin + 1);
            double var3 = MCDn->GetBinContent(bin + 1) - MC->GetBinContent(bin + 1);
            double var4 = MCDnscale->GetBinContent(bin + 1) - MC->GetBinContent(bin + 1);
            
            if (var1 > 0) up_var += pow(var1, 2);
            else dn_var += pow(var1, 2);

            if (var2 > 0) up_var += pow(var2, 2);
            else dn_var += pow(var2, 2);

            if (var3 > 0) up_var += pow(var3, 2);
            else dn_var += pow(var3, 2);

            if (var4 > 0) up_var += pow(var4, 2);
            else dn_var += pow(var4, 2);

            MCUpFinal->SetBinContent(bin + 1, sqrt(up_var));
            MCDnFinal->SetBinContent(bin + 1, sqrt(dn_var));
        }
        std::cout << "Scale and smearing calculation end "<<std::endl;
        data->SetLineColor(kBlack);
        data->SetMarkerStyle(20);
        data->SetMarkerSize(0.6);
        data->Draw("p E1 X0 SAME");

        //***Ratio Plot (Bottom Panel)***
        c->cd(2);

        gPad->SetBottomMargin(0.20);
        gPad->SetTopMargin(0.02);
        gPad->SetLeftMargin(0.10);

        TH1F * Ratio = new TH1F("Ratio","Ratio", MC->GetSize() - 2, MC->GetXaxis()->GetXmin(), MC->GetXaxis()->GetXmax());
        TH1F * sigma_up = new TH1F("sigma_up","sigma_up", MC->GetSize() - 2, MC->GetXaxis()->GetXmin(), MC->GetXaxis()->GetXmax());
        TH1F * sigma_dn = new TH1F("sigma_dn","sigma_dn", MC->GetSize() - 2, MC->GetXaxis()->GetXmin(), MC->GetXaxis()->GetXmax());

        for(int bin = 0; bin < MC->GetSize() - 2; bin++){
          if (MC->GetBinContent(bin + 1) == 0){
                  Ratio   ->SetBinContent(bin + 1, 0);
                  Ratio   ->SetBinError(bin + 1, 0);
                  sigma_up->SetBinContent(bin + 1, 0);
                  sigma_dn->SetBinContent(bin + 1, 0);
          }else{
                  Ratio   ->SetBinContent(bin + 1, float(data->GetBinContent(bin + 1))/MC->GetBinContent(bin + 1));
                  Ratio   ->SetBinError(bin + 1, 1./pow((data->GetBinContent(bin + 1)),2));
                  sigma_up->SetBinContent(bin + 1, ((float(MCUpFinal->GetBinContent(bin + 1) + MC->GetBinContent(bin + 1)))/MC->GetBinContent(bin + 1)));
                  sigma_dn->SetBinContent(bin + 1, ((float((- MCDnFinal->GetBinContent(bin + 1)) + MC->GetBinContent(bin + 1)))/MC->GetBinContent(bin + 1)));
          }
        }

        double temp;
        for(int bin = 0; bin < Ratio->GetSize() - 2; bin++){
          temp = Ratio->GetBinContent(bin + 1);
          Ratio->SetBinContent( bin + 1, temp - 1);
        }

        for(int bin = 0; bin < sigma_up->GetSize() - 2; bin++){
          temp = sigma_up->GetBinContent(bin + 1);
          sigma_up->SetBinContent( bin + 1, temp - 1);
        }

        for(int bin = 0; bin < sigma_dn->GetSize() - 2; bin++){
          temp = sigma_dn->GetBinContent(bin + 1);
          sigma_dn->SetBinContent( bin + 1, temp - 1);
        }

        sigma_up->SetTitle("");
        sigma_up->SetFillColor(kGray);
        sigma_up->SetLineColor(kGray);
        sigma_up->SetMaximum(1.0);
        sigma_up->SetMinimum(-1.0);
        std::string ytitle_plots = "Leading Jet " + var_name;
        sigma_up->GetXaxis()->SetTitle(ytitle_plots.c_str());
        sigma_up->GetYaxis()->SetTitle("(Data/MC)-1");
        sigma_up->GetYaxis()->SetTitleOffset(0.6);
        sigma_up->Draw("HIST");
        sigma_up->GetXaxis()->SetLabelSize(0.07);
        sigma_up->GetXaxis()->SetTitleSize(0.07);
        sigma_up->GetYaxis()->SetLabelSize(0.07);
        sigma_up->GetYaxis()->SetTitleSize(0.07);

        sigma_dn->SetFillColor(kGray+2);
        sigma_dn->SetLineColor(kGray+2);
        sigma_dn->Draw("HIST SAME");

        gPad->RedrawAxis();

        Ratio->SetMarkerStyle(20);
        Ratio->SetMarkerSize(0.6);
        Ratio->Draw("p E1 X0 SAME");

        c->cd(1);
        TLegend * leg = new TLegend(0.74,0.45,0.95,0.75);
        leg->SetBorderSize(0);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.07);
        leg->AddEntry(MC, "DY + t#bar{t}  MC", "f" );
        leg->AddEntry(data, "Data", "p");
        //leg->AddEntry(UncertaintyBand, "Uncertainty Band", "f");
        leg->AddEntry(sigma_up, "JEC Uncertainty", "f");
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
        .Define("Jet_Mask", 
        "(Jet_pt > 30) && (Jet_jetId == 6 || Jet_jetId == 2) && Jet_ZZMask == false")
        
        .Define("Jet_Mask_scaleUp", 
        "(Jet_scaleUp_pt > 30) && (Jet_jetId == 6 || Jet_jetId == 2) && ((abs(Jet_eta) < 2.5 || abs(Jet_eta) > 3) || Jet_scaleUp_pt > 50) && Jet_ZZMask == false")
        
        .Define("Jet_Mask_scaleDown", 
        "(Jet_scaleDn_pt > 30) && (Jet_jetId == 6 || Jet_jetId == 2) && ((abs(Jet_eta) < 2.5 || abs(Jet_eta) > 3) || Jet_scaleDn_pt > 50) && Jet_ZZMask == false") 
        
        .Define("Jet_Mask_smearUp", 
        "(Jet_smearUp_pt > 30) && (Jet_jetId == 6 || Jet_jetId == 2) && ((abs(Jet_eta) < 2.5 || abs(Jet_eta) > 3) || Jet_smearUp_pt > 50) && Jet_ZZMask == false") 
        
        .Define("Jet_Mask_smearDown", 
        "(Jet_smearDn_pt > 30) && (Jet_jetId == 6 || Jet_jetId == 2) && ((abs(Jet_eta) < 2.5 || abs(Jet_eta) > 3) || Jet_smearDn_pt > 50) && Jet_ZZMask == false") 
        
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
        //(abs(Jet_eta) < 2.5 || abs(Jet_eta) > 3) || Jet_pt > 50) //&& Jet_ZZMask == false// && Jet_ZZMask == false //|| Jet_jetId == 2 && Jet_ZZMask == false//((abs(Jet_eta) < 2.5 || abs(Jet_eta) > 3) || Jet_pt > 50) &&
        //((abs(Jet_eta) < 2.5 || abs(Jet_eta) > 3) || Jet_pt > 50) &&
        
    // Step 2: Apply the event-level filter
    auto skim_rdf = sorted_rdf.Filter("ZLCand_lepIdx >= 0 && SortedJet_eta.size() > 0 && Flag_JetVetoed == 0") //Z control region used (nZCand)
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
        .Define("Jet_Mask", 
        "(Jet_pt > 30) && (Jet_jetId == 6 || Jet_jetId == 2) && Jet_ZZMask == false") // ((abs(Jet_eta) < 2.5 || abs(Jet_eta) > 3) || Jet_pt > 50)//&& Jet_ZZMask == false//|| Jet_jetId == 2 && Jet_ZZMask == false //((abs(Jet_eta) < 2.5 || abs(Jet_eta) > 3) || Jet_pt > 50) && 
        .Define("FilteredJet_eta", "Jet_eta[Jet_Mask]")
        .Define("FilteredJet_pt", "Jet_pt[Jet_Mask]")
        .Define("FilteredJet_Njets", "FilteredJet_pt.size()")

        .Define("SortedJet_eta", "FilteredJet_eta[Reverse(Argsort(FilteredJet_pt))]")
        .Define("SortedJet_pt", "FilteredJet_pt[Reverse(Argsort(FilteredJet_pt))]");

    // Apply event-level filter and define leadEta and leadPt
    auto skim_data = sorted_data.Filter("ZLCand_lepIdx >= 0 && SortedJet_eta.size() > 0 && Flag_JetVetoed == 0") //ZLCand_lepIdx nZCand
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
  std::string var = "leadEta";//"leadNjets";//"leadEta";//leadPt //leadchEmEF //leadEta //leadneEmEF //leadNJets
  std::string var_name = "#eta"; //"p_{T}";//"#eta"; ////Neutral Em Energy Fraction //charged Em Energy Fraction //Number of jets
  std::string period = "preEE2022"; //"preEE2022" "postBPix2023" "preBPix2023"

  // Define a map that associates 'var' with the corresponding binning
  std::map<std::string, std::tuple<int, double, double>> var_binning;
  
  // binning configurations
  var_binning["leadEta"] = std::make_tuple(47, -4.7, +4.7);
  var_binning["leadPt"] = std::make_tuple(40, 30, 200.0);
  var_binning["leadNjets"] = std::make_tuple(8, 2, 10);

  int n_bins;
  double min_val, max_val;
  std::tie(n_bins, min_val, max_val) = var_binning[var];

  std::cout<<"period preEE2022"<<std::endl;
  std::string fTT   = "/eos/user/m/mmanoni/HZZ_samples22_newprod/MC_JetClean_Jes/PROD_samplesNano_2022_MC_3e5e00ec/TTto2L2Nu/ZZ4lAnalysis.root";
  std::string fDY   = "/eos/user/m/mmanoni/HZZ_samples22_newprod/MC_JetClean_Jes/PROD_samplesNano_2022_MC_3e5e00ec/DYJetsToLL/ZZ4lAnalysis.root";
  std::string fdata = "/eos/user/m/mmanoni/HZZ_samples22_newprod/Data_JetClean_Jes/PROD_samplesNano_2022_Data_3e5e00ec/Data_eraCD_preEE.root";
  double lumi  = 7.98;

  /*std::cout<<"period postEE2022"<<std::endl;
  string fDY   = "/eos/user/m/mmanoni/HZZ_samples22_newprod/MC_JetClean_Jes/PROD_samplesNano_2022EE_MC_3e5e00ec/DYJetsToLL/ZZ4lAnalysis.root";
  string fTT   = "/eos/user/m/mmanoni/HZZ_samples22_newprod/MC_JetClean_Jes/PROD_samplesNano_2022EE_MC_3e5e00ec/TTto2L2Nu/ZZ4lAnalysis.root"; 
  string fdata = "/eos/user/m/mmanoni/HZZ_samples22_newprod/Data_JetClean_Jes/PROD_samplesNano_2022_Data_3e5e00ec/Data_eraEFG_postEE.root";
  double lumi  = 26.67;*/

  /*std::cout<<"period preBPix2023"<<std::endl;
  string fDY   = "/eos/user/m/mmanoni/HZZ_samples23_newprod/MC_JetClean_Jes/PROD_samplesNano_2023preBPix_MC_3e5e00ec/DYJetsToLL/ZZ4lAnalysis.root";
  string fTT   = "/eos/user/m/mmanoni/HZZ_samples23_newprod/MC_JetClean_Jes/PROD_samplesNano_2023preBPix_MC_3e5e00ec/TTto2L2Nu/ZZ4lAnalysis.root";
  string fdata = "/eos/user/m/mmanoni/HZZ_samples23_newprod/Data_JetClean_Jes/PROD_samplesNano_2023_Data_3e5e00ec/Data_eraC_preBPix.root";
  double lumi  = 17.8;*/

  /*std::cout<<"period preBPix2023"<<std::endl;
  string fDY   = "/eos/user/m/mmanoni/HZZ_samples23_newprod/MC_updated_Trig_newJetID_oldEleSS/PROD_samplesNano_2023preBPix_MC_03abdca3/DYJetsToLL/ZZ4lAnalysis.root";
  string fTT   = "/eos/user/m/mmanoni/HZZ_samples23_newprod/MC_updated_Trig_newJetID_oldEleSS/PROD_samplesNano_2023preBPix_MC_03abdca3/TTto2L2Nu/ZZ4lAnalysis.root";
  string fdata = "/eos/user/m/mmanoni/HZZ_samples23_newprod/Data_updated_Trig_newJetID_oldEleSS/PROD_samplesNano_2023_Data_03abdca3/Data_eraC_preBPix.root";
  double lumi  = 17.8;*/

  /*std::cout<<"period postBPix2023"<<std::endl;
  string fDY   = "/eos/user/m/mmanoni/HZZ_samples23_newprod/MC_JetClean_Jes/PROD_samplesNano_2023postBPix_MC_3e5e00ec/DYJetsToLL/ZZ4lAnalysis.root";
  string fTT   = "/eos/user/m/mmanoni/HZZ_samples23_newprod/MC_JetClean_Jes/PROD_samplesNano_2023postBPix_MC_3e5e00ec/TTto2L2Nu/ZZ4lAnalysis.root";
  string fdata = "/eos/user/m/mmanoni/HZZ_samples23_newprod/Data_JetClean_Jes/PROD_samplesNano_2023_Data_3e5e00ec/Data_eraD_postBPix.root";
  double lumi  = 9.5;*/


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

  std::string subdir = "2022preEE/"; 
  std::filesystem::create_directories(subdir);


  std::cout << "Saving distributions into root file ..." << std::endl;
  std::string filename = subdir + var +"_"+ period + "_scale_smear_Clean_NOHorncut_Norm_ZLCand_lepIdx.root";  // Concatenate the value of 'var' with the string
  TFile *outfile = new TFile(filename.c_str(), "RECREATE");

  outfile    -> cd();
  hist_data  -> Write();
  hist_DY    -> Write();
  hist_DY_up -> Write();
  hist_DY_dn -> Write();
  hist_DY_up_scale -> Write();
  hist_DY_dn_scale -> Write();
  outfile    -> Close();
  

  TCanvas * canvas = new TCanvas();
  std::string filename_plots = subdir + var +"_"+ period + "_plots_scale_smear_Clean_NOHorncut_Norm_ZLCand_lepIdx";
  DrawRatioPlot(filename_plots, var_name, canvas,
                hist_data.GetPtr(),
                hist_DY.GetPtr(), hist_DY_up.GetPtr(), hist_DY_dn.GetPtr(), hist_DY_up_scale.GetPtr(), hist_DY_dn_scale.GetPtr(),
                lumi);

}