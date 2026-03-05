#!/bin/env python3
# PLOTS FILLER - Clean Version with Data and ZX

import ROOT
import argparse
import os
import array

ROOT.EnableImplicitMT()

# =========================
# PATHS & LUMI
# =========================
MC_PATHS = {
    "2022": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2022_MC/",
    "2022EE": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2022EE_MC/",
    "2023preBPix": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2023preBPix_MC/",
    "2023postBPix": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2023postBPix_MC/",
    "2024": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_MC/",
}

DATA_PATHS = {
    "2022": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2022_Data/Data_eraCD_preEE.root",
    "2022EE": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2022_Data/Data_eraEFG_postEE.root",
    "2023preBPix": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2023_Data/Data_eraC_preBPix.root",
    "2023postBPix": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2023_Data/Data_eraD_postBPix.root",
    "2024": [
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Cv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Dv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Ev1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Fv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Gv2/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Hv2/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Iv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Iv2/ZZ4lAnalysis.root",

            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Cv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Dv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Ev1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Fv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Gv2/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Hv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Iv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Iv2/ZZ4lAnalysis.root",

            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Cv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Dv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Ev1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Fv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Gv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Hv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Iv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Iv2/ZZ4lAnalysis.root",

            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Cv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Dv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Ev1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Fv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Gv2/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Hv2/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Iv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Iv2/ZZ4lAnalysis.root",

            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Cv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Dv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Ev1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Fv2/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Gv3/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Hv2/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Iv2/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Iv2v2/ZZ4lAnalysis.root"
            ]
}

ZX_PATHS = {
    "2022": "/afs/cern.ch/user/m/mmanoni/HZZ_plotter_CMSSW14/CMSSW_14_1_6/src/HZZ_plotter/HZZ_mass4l_plotter/ZX_results_2022.root",
    "2022EE": "/afs/cern.ch/user/m/mmanoni/HZZ_plotter_CMSSW14/CMSSW_14_1_6/src/HZZ_plotter/HZZ_mass4l_plotter/ZX_results_2022EE.root",
    "2023preBPix": "/afs/cern.ch/user/m/mmanoni/HZZ_plotter_CMSSW14/CMSSW_14_1_6/src/HZZ_plotter/HZZ_mass4l_plotter/ZX_results_2023preBPix.root",
    "2023postBPix": "/afs/cern.ch/user/m/mmanoni/HZZ_plotter_CMSSW14/CMSSW_14_1_6/src/HZZ_plotter/HZZ_mass4l_plotter/ZX_results_2023postBPix.root",
    "2024": "/afs/cern.ch/user/m/mmanoni/HZZ_plotter_CMSSW14/CMSSW_14_1_6/src/HZZ_plotter/HZZ_mass4l_plotter/ZX_results_2024.root",
}

LUMI = {
    "2022": 7.98,
    "2022EE": 26.67,
    "2023preBPix": 18.06,
    "2023postBPix": 9.69,
    "2024": 108.82,
}

# =========================
# SUM OF WEIGHTS
# =========================
def get_genEventSumw(file):
    hCounters = file.Get("Counters")
    return hCounters.GetBinContent(40)

# =========================
# HISTOGRAM DEFINITIONS
# =========================
def define_histograms(df_full, df_window, samplename, isMC):
    histos = {}
    weight_col = "weight" if isMC else None

    def book(model, var):
        if weight_col:
            return df_full.Histo1D(model, var, weight_col)
        else:
            return df_full.Histo1D(model, var)

    bins_costhetaZ1 = array.array('d', [-1.0, -0.75, -0.50, -0.25, 0.0, 0.25, 0.50, 0.75, 1.0])

    histos[f"costhetaZ1_{samplename}_FULL"] = book(
        ROOT.RDF.TH1DModel(f"costhetaZ1_{samplename}_FULL",
                           "cos(#theta_{Z1});cos(#theta_{Z1});Events",
                           len(bins_costhetaZ1) - 1, bins_costhetaZ1),
        "costheta1"
    )

    histos[f"costhetaZ1_{samplename}_105to160"] = book(
        ROOT.RDF.TH1DModel(f"costhetaZ1_{samplename}_105to160",
                           "cos(#theta_{Z1});cos(#theta_{Z1});Events",
                           len(bins_costhetaZ1) - 1, bins_costhetaZ1),
        "costheta1"
    )

    return histos

# =========================
# RUN SINGLE SAMPLE
# =========================
def run_sample(sample_name, filepath, output_file, period, isMC=True):
    df_full = ROOT.RDataFrame("ZZTree/candTree", filepath)
    df_window = df_full.Filter("ZZMass > 105 && ZZMass < 160")

    if isMC:
        f = ROOT.TFile.Open(filepath)
        genEventSumw = get_genEventSumw(f)
        f.Close()

        lumi = LUMI[period]
        weight_expr = f"overallEventWeight * {lumi} * 1000 * dataMCWeight/{genEventSumw}"
        df_full = df_full.Define("weight", weight_expr)
        df_window = df_window.Define("weight", weight_expr)

    histos = define_histograms(df_full, df_window, sample_name, isMC)

    output_file.cd()
    for h in histos.values():
        h.GetValue().Write()

# =========================
# RUN MC
# =========================
def run_mc(period):
    path = MC_PATHS.get(period)
    if not path:
        raise ValueError(f"Unknown MC period: {period}")

    samples = ["WWZ", "WZZ", "ZZZ", "ggTo4mu_Contin_MCFM701", "ggTo4e_Contin_MCFM701",
               "ggTo4tau_Contin_MCFM701", "ggTo2e2mu_Contin_MCFM701", "ggTo2e2tau_Contin_MCFM701",
               "ggTo2mu2tau_Contin_MCFM701", "ZZTo4l", "VBFH125", "ggH125", "WplusH125",
               "WminusH125", "ZH125", "ttH125"]

    fout = ROOT.TFile.Open(f"H4l_MC_{period}_DIFF.root", "RECREATE")
    for sample in samples:
        filename = os.path.join(path, sample, "ZZ4lAnalysis_SKIMMED.root")
        run_sample(sample, filename, fout, period, isMC=True)
    fout.Close()

# =========================
# RUN DATA
# =========================
def run_data(period):
    paths = DATA_PATHS.get(period)
    if not paths:
        raise ValueError(f"Unknown DATA period: {period}")

    fout = ROOT.TFile.Open(f"H4l_Data_{period}_DIFF.root", "RECREATE")
    files = paths if isinstance(paths, list) else [paths]
    for fpath in files:
        run_sample("DATA", fpath, fout, period, isMC=False)
    fout.Close()

# =========================
# RUN ZX (special treatment)
# =========================
def run_zx(period):
    path = ZX_PATHS.get(period)
    if not path or not os.path.exists(path):
        print(f"[WARNING] ZX file not found for period {period}: {path}")
        return

    df = ROOT.RDataFrame("candTree", path)

    df = df.Define("m4l", "ZZMass") \
           .Define("costhetaZ1", "costheta1") \
           .Define("Z1mass", "Z1Mass") \
           .Define("Z2mass", "Z2Mass") \
           .Define("Z1flav", "Z1Flav") \
           .Define("Z2flav", "Z2Flav") \
           .Define("weight", "weight1")  # ZX weight

    df_SR = df.Filter("m4l > 105 && m4l < 160")

    fout = ROOT.TFile.Open(f"H4l_ZX_{period}_DIFF.root", "RECREATE")
    histos = define_histograms(df, df_SR, "ZX", isMC=True)

    for hname, h in histos.items():
        h_clone = ROOT.TH1F(hname, hname, h.GetValue().GetNbinsX(),
                             h.GetValue().GetXaxis().GetXmin(), h.GetValue().GetXaxis().GetXmax())
        for i in range(1, h_clone.GetNbinsX() + 1):
            h_clone.SetBinContent(i, h.GetValue().GetBinContent(i))
            h_clone.SetBinError(i, h.GetValue().GetBinError(i))
        h_clone.Write()
    fout.Close()

# =========================
# MAIN
# =========================
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--period", default="2022")
    parser.add_argument("--data", action="store_true")
    parser.add_argument("--zx", action="store_true")
    args = parser.parse_args()

    if args.data:
        run_data(args.period)
    elif args.zx:
        run_zx(args.period)
    else:
        run_mc(args.period)