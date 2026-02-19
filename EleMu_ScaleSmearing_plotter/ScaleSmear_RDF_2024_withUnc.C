// ============================================================================
// ScaleSmear_RDF_2024_withUnc.C
// ROOT plotter with MC uncertainty band (stat + scale/smear)
// ============================================================================
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>
#include <TLine.h>
#include <iostream>
#include <vector>
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
double getGenEventSumw(const std::vector<std::string> &fnames) {
    double genEventSumw = 0;
    for (auto &fname : fnames) {
        TFile f(fname.c_str());
        if (f.IsZombie()) { std::cerr << "Cannot open: " << fname << std::endl; continue; }
        TTree *Runs = (TTree*)f.Get("Runs");
        if (!Runs) { std::cerr << "No Runs tree in: " << fname << std::endl; continue; }
        double sumw; Long64_t cnt;
        Runs->SetBranchAddress("genEventSumw", &sumw);
        Runs->SetBranchAddress("genEventCount", &cnt);
        for (Long64_t i=0; i<Runs->GetEntries(); i++) {
            Runs->GetEntry(i);
            genEventSumw += sumw;
        }
        f.Close();
    }
    return genEventSumw;
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
void ScaleSmear_RDF_2024_withUnc() {

    // --------------------------
    // Config
    // --------------------------
    std::string YEAR    = "2024";
    std::string CHANNEL = "mu";  // "mu" or "e"
    double lumi         = 108.8;  // fb^-1

    // --------------------------
    // Bin settings
    // --------------------------
    const int NBINS = 99;
    const double XMIN = 60;
    const double XMAX = 120;

    // --------------------------
    // Set scale/smear uncertainties
    // --------------------------
    double scale_unc, smear_unc;
    if (CHANNEL == "mu") { scale_unc = 0.0005; smear_unc = 0.025; } // 0.05%, 2.5%
    if (CHANNEL == "e")  { scale_unc = 0.0025; smear_unc = 0.07; }  // 0.25%, 7%

    // --------------------------
    // Files (fill DY_files, TT_files, data_files)
    // --------------------------
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
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/TTto2L2Nu/TTto2L2Nu_3/TTto2L2Nu/ZZ4lAnalysis.root",
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
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/MuonEG2024Iv2v2/ZZ4lAnalysis.root",
    };

    // --------------------------
    // Get genEventSumw
    // --------------------------
    double cnt_dy = getGenEventSumw(DY_files);
    double cnt_tt = getGenEventSumw(TT_files);
    if (cnt_dy == 0 || cnt_tt == 0)
        throw std::runtime_error("genEventSumw = 0 for some sample!");

    // --------------------------
    // RDataFrames
    // --------------------------
    RDataFrame d_dy("Events", DY_files);
    RDataFrame d_tt("Events", TT_files);
    RDataFrame d_dt("Events", data_files);

    auto defZ = [&](RDataFrame &df) {
        return df
            .Define("Z1_mass" , pickZ1, {"bestZIdx","ZCand_mass","ZCand_flav"})
            .Define("Z1_flav" , pickFlav, {"bestZIdx","ZCand_flav"})
            .Define("GoodZ1", [&](int flav){ return selectChannel(flav, CHANNEL); }, {"Z1_flav"});
    };

    auto D_dy = defZ(d_dy)
        .Filter("GoodZ1")
        .Define("w", [=](float w){ return w*1000.*lumi/cnt_dy; }, {"overallEventWeight"});
    auto D_tt = defZ(d_tt)
        .Filter("GoodZ1")
        .Define("w", [=](float w){ return w*1000.*lumi/cnt_tt; }, {"overallEventWeight"});
    auto D_dt = defZ(d_dt).Filter("GoodZ1");

    // --------------------------
    // Fill histograms
    // --------------------------
    TH1D *hDY = new TH1D("hDY","",NBINS,XMIN,XMAX);
    TH1D *hTT = new TH1D("hTT","",NBINS,XMIN,XMAX);
    TH1D *hDT = new TH1D("hDT","",NBINS,XMIN,XMAX);

    D_dy.Histo1D({"hDY","",NBINS,XMIN,XMAX},"Z1_mass","w")->Copy(*hDY);
    D_tt.Histo1D({"hTT","",NBINS,XMIN,XMAX},"Z1_mass","w")->Copy(*hTT);
    D_dt.Histo1D({"hDT","",NBINS,XMIN,XMAX},"Z1_mass")->Copy(*hDT);

    // Total MC
    TH1D *hMC = (TH1D*)hDY->Clone("hMC");
    hMC->Add(hTT);

    // Normalize MC to Data
    double scale = hDT->Integral() / hMC->Integral();
    hMC->Scale(scale);

    // --------------------------
    // Create MC uncertainty band (stat + scale/smear)
    // --------------------------
    TH1D *hMC_up   = (TH1D*)hMC->Clone("hMC_up");
    TH1D *hMC_down = (TH1D*)hMC->Clone("hMC_down");

    for (int i=1; i <= hMC->GetNbinsX(); i++) {
        double nom = hMC->GetBinContent(i);
        double stat = hMC->GetBinError(i); // MC stat
        double syst = nom * sqrt(scale_unc*scale_unc + smear_unc*smear_unc);
        double total = sqrt(stat*stat + syst*syst);
        hMC_up->SetBinContent(i, nom + total);
        hMC_down->SetBinContent(i, std::max(0., nom - total));
    }

    // --------------------------
    // Plotting
    // --------------------------
    TCanvas *c = new TCanvas("c","",900,800);
    c->Divide(1,2);

    // Top pad
    c->cd(1);
    gPad->SetPad(0,0.3,1,1);
    gPad->SetBottomMargin(0.05);

    hMC->SetLineColor(kBlue);      
    hMC->SetLineWidth(2);
    hMC->SetFillStyle(0);
    hMC->Draw("HIST");

    // Uncertainty band
    hMC_up->SetFillColor(kBlue);
    hMC_up->SetFillStyle(3004); // hatched
    hMC_up->Draw("E2 SAME");

    hDT->SetMarkerStyle(20);
    hDT->SetMarkerSize(1.0);
    hDT->SetLineColor(kBlack);
    hDT->Draw("E SAME");

    TLegend *leg = new TLegend(0.65,0.7,0.88,0.88);
    leg->AddEntry(hDT,"Data","lp");
    leg->AddEntry(hMC,"MC (DY+tt)","l");
    leg->AddEntry(hMC_up,"MC unc.","f");
    leg->Draw();

    CMS_lumi cmsLabel;
    cmsLabel.set_lumi((TPad*)gPad, lumi);

    // Ratio pad
    c->cd(2);
    gPad->SetPad(0,0.0,1,0.3);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.25);

    TH1D *hRatio = (TH1D*)hMC->Clone("hRatio");
    hRatio->Divide(hDT);
    hRatio->SetMinimum(0.5);
    hRatio->SetMaximum(1.5);
    hRatio->GetYaxis()->SetTitle("MC/Data");
    hRatio->GetYaxis()->SetTitleSize(0.10);
    hRatio->GetYaxis()->SetLabelSize(0.08);
    hRatio->GetXaxis()->SetTitle("m_{ll} [GeV]");
    hRatio->GetXaxis()->SetTitleSize(0.12);
    hRatio->GetXaxis()->SetLabelSize(0.10);
    hRatio->Draw("E");

    // Uncertainty band in ratio
    TH1D *hRatio_up   = (TH1D*)hMC_up->Clone("hRatio_up");
    TH1D *hRatio_down = (TH1D*)hMC_down->Clone("hRatio_down");
    hRatio_up->Divide(hDT);
    hRatio_down->Divide(hDT);
    hRatio_up->SetFillColor(kBlue);
    hRatio_up->SetFillStyle(3004);
    hRatio_up->Draw("E2 SAME");

    c->SaveAs(Form("plot_%s_%s_RDF_unc.png", CHANNEL.c_str(), YEAR.c_str()));
    c->SaveAs(Form("plot_%s_%s_RDF_unc.pdf", CHANNEL.c_str(), YEAR.c_str()));

    std::cout << "\nSaved plot with MC uncertainties: plot_" << CHANNEL << "_" << YEAR << "_RDF_unc\n";
}
