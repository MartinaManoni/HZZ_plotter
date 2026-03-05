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
    "2022": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2022_Data/Data_eraCD_preEE_SKIMMED.root",
    "2022EE": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2022_Data/Data_eraEFG_postEE_SKIMMED.root",
    "2023preBPix": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2023_Data/Data_eraC_preBPix_SKIMMED.root",
    "2023postBPix": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2023_Data/Data_eraD_postBPix_SKIMMED.root",
    "2024": [
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Cv1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Dv1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Ev1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Fv1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Gv2/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Hv2/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Iv1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Iv2/ZZ4lAnalysis_SKIMMED.root",

            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Cv1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Dv1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Ev1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Fv1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Gv2/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Hv1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Iv1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Iv2/ZZ4lAnalysis_SKIMMED.root",

            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Cv1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Dv1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Ev1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Fv1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Gv1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Hv1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Iv1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Iv2/ZZ4lAnalysis_SKIMMED.root",

            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Cv1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Dv1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Ev1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Fv1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Gv2/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Hv2/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Iv1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Iv2/ZZ4lAnalysis_SKIMMED.root",

            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Cv1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Dv1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Ev1/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Fv2/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Gv3/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Hv2/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Iv2/ZZ4lAnalysis_SKIMMED.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Iv2v2/ZZ4lAnalysis_SKIMMED.root"
            ]
}

ZX_PATHS = {
    "2022": "/afs/cern.ch/user/m/mmanoni/HZZ_plotter_CMSSW14/CMSSW_14_1_6/src/HZZ_plotter/HZZ_mass4l_plotter/ZX_results_2022_obs.root",
    "2022EE": "/afs/cern.ch/user/m/mmanoni/HZZ_plotter_CMSSW14/CMSSW_14_1_6/src/HZZ_plotter/HZZ_mass4l_plotter/ZX_results_2022EE_obs.root",
    "2023preBPix": "/afs/cern.ch/user/m/mmanoni/HZZ_plotter_CMSSW14/CMSSW_14_1_6/src/HZZ_plotter/HZZ_mass4l_plotter/ZX_results_2023preBPix_obs.root",
    "2023postBPix": "/afs/cern.ch/user/m/mmanoni/HZZ_plotter_CMSSW14/CMSSW_14_1_6/src/HZZ_plotter/HZZ_mass4l_plotter/ZX_results_2023postBPix_obs.root",
    "2024": "/afs/cern.ch/user/m/mmanoni/HZZ_plotter_CMSSW14/CMSSW_14_1_6/src/HZZ_plotter/HZZ_mass4l_plotter/ZX_results_2024_obs.root",
}

LUMI = {
    "2022": 7.98,
    "2022EE": 26.67,
    "2023preBPix": 18.06,
    "2023postBPix": 9.69,
    "2024": 108.82,
}

OBSERVABLES = {

    "costhetaZ1": {
        "var": "costheta1",
        "title": "cos(#theta_{Z1});cos(#theta_{Z1});Events",
        "bins": [-1.0,-0.75,-0.50,-0.25,0.0,0.25,0.50,0.75,1.0]
    },

    "pT4l": {
        "var": "ZZPt",
        "title": "p_{T}^{4l};p_{T}^{4l} [GeV];Events",
        "bins": [0,10,16,22,28,36,46,60,80,106,146,200]
    },

    "rapidity4l": {
        "var": "ZZy",
        "title": "|y_{4l}|;|y_{4l}|;Events",
        "bins": [0,0.12,0.24,0.36,0.5,0.6,0.75,0.9,1.1,1.3,2.5]
    }

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

    def book(df, model, var):
        if weight_col:
            return df.Histo1D(model, var, weight_col)
        else:
            return df.Histo1D(model, var)

    for obs, cfg in OBSERVABLES.items():

        bins = array.array('d', cfg["bins"])
        nbins = len(bins) - 1
        var = cfg["var"]
        title = cfg["title"]

        # FULL histogram
        histos[f"{obs}_{samplename}_FULL"] = book(
            df_full,
            ROOT.RDF.TH1DModel(
                f"{obs}_{samplename}_FULL",
                title,
                nbins,
                bins
            ),
            var
        )

        # Signal window histogram
        histos[f"{obs}_{samplename}_105to160"] = book(
            df_window,
            ROOT.RDF.TH1DModel(
                f"{obs}_{samplename}_105to160",
                title,
                nbins,
                bins
            ),
            var
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

    files = paths if isinstance(paths, list) else [paths]

    # Build ONE dataframe with ALL files
    df_full = ROOT.RDataFrame("ZZTree/candTree", files)

    df_window = df_full.Filter("ZZMass > 105 && ZZMass < 160")

    histos = define_histograms(df_full, df_window, "DATA", isMC=False)

    fout = ROOT.TFile.Open(f"H4l_Data_{period}_DIFF.root", "RECREATE")

    for h in histos.values():
        h.GetValue().Write()

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

    # Map columns to OBSERVABLES
    df = df.Define("costhetaZ1", "costheta1") \
           .Define("pT4l","ZZPt")\
           .Define("rapidity4l", "ZZy")\
           .Define("weight", "weight1")  # ZX weight

    df_SR = df.Filter("ZZMass > 105 && ZZMass < 160")

    fout = ROOT.TFile.Open(f"H4l_ZX_{period}_DIFF.root", "RECREATE")

    for obs, cfg in OBSERVABLES.items():
        bins = array.array('d', cfg["bins"])
        for suffix, rdf in [("_FULL", df), ("_105to160", df_SR)]:
            hname = f"{obs}_ZX{suffix}"
            # Directly create histogram with correct binning and weights
            h = rdf.Histo1D(
                ROOT.RDF.TH1DModel(hname, hname, len(bins)-1, bins),
                cfg["var"],
                "weight"
            ).GetValue()
            h.Write()

    fout.Close()
    print(f"[INFO] ZX histograms for {period} written to H4l_ZX_{period}_DIFF.root")
# =========================
# MAIN
# =========================
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--periods", nargs="+", default=["2022", "2022EE", "2023preBPix", "2023postBPix", "2024"], help="List of periods/years to run over")
    parser.add_argument("--mc", action="store_true")
    parser.add_argument("--data", action="store_true")
    parser.add_argument("--zx", action="store_true")
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()

    for period in args.periods:
        print(f"Running for period: {period}")
        if args.data:
            run_data(period)
        elif args.zx:
            run_zx(period)
        elif args.mc:
            run_mc(period)
        elif args.all:
            run_mc(period)
            run_data(period)
            run_zx(period)
