#!/bin/env python3
# PLOTS FILLER - Clean Version with Data and ZX

import ROOT
import argparse
import os
import array

ROOT.EnableImplicitMT()
ROOT.gROOT.SetBatch(True)

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
    "2022": "/eos/user/m/mmanoni/test_cleaning/2022_Data/Data_eraCD_preEE_SKIMMED.root",
    "2022EE": "/eos/user/m/mmanoni/test_cleaning/2022_Data/Data_eraEFG_postEE_SKIMMED.root",
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

FINAL_STATES_OS = {
    "4e":    "Z1Flav==-121 && Z2Flav==-121",
    "4mu":   "Z1Flav==-169 && Z2Flav==-169",
    "2e2mu": "Z1Flav==-121 && Z2Flav==-169",
    "2mu2e": "Z1Flav==-169 && Z2Flav==-121",
    "2e2mu_com":"(Z1Flav==-121 && Z2Flav==-169) || (Z1Flav==-169 && Z2Flav==-121)",
}

FINAL_STATES_SS = {
    "4e":    "(abs(Z1Flav)==121) && (Z2Flav==121)",
    "4mu":   "(abs(Z1Flav)==169) && (Z2Flav==169)",
    "2e2mu": "(abs(Z1Flav)==121) && (Z2Flav==169)",
    "2mu2e": "(abs(Z1Flav)==169) && (Z2Flav==121)",
    "2e2mu_com":"((abs(Z1Flav)==121) && (Z2Flav==169)) || ((abs(Z1Flav)==169) && (Z2Flav==121))",
}

OBSERVABLES = {
    "Nj": {
        "var": "Nj_clipped",
        "title": "N_{jets};N_{jets};Events",
        "bins": [0,1,2,3,4, 5]
    },

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
def run_sample(sample, filepath, output_file, period, isMC=True):

    #DATAFRAME SR
    if sample == "DYJetsToLL" and period == "2024":
        print("--------using DY 2024 logic")
        sr_filepath = "/eos/user/m/mmanoni/DY_HADDED_MARTINA.root"
        
        df_full = ROOT.RDataFrame("ZZTree/candTree", sr_filepath)
        df_window = df_full.Filter("ZZMass > 105 && ZZMass < 160")
    else:
        print("--------using other logic for sample", sample)
        df_full = ROOT.RDataFrame("ZZTree/candTree", filepath)
        df_window = df_full.Filter("ZZMass > 105 && ZZMass < 160")

    #df_window = df_full.Filter("ZZMass > 105 && ZZMass < 160")
    df_CR = ROOT.RDataFrame("CRZLLTree/candTree", filepath)


    df_CRZL = ROOT.RDataFrame("CRZLTree/candTree", filepath)

    if isMC:
        f = ROOT.TFile.Open(filepath)
        genEventSumw = get_genEventSumw(f)
        f.Close()

        lumi = LUMI[period]
        weight_expr = f"overallEventWeight * {lumi} * 1000 / {genEventSumw}"
        #DATAFRAME SR
        df_full = df_full.Define("weight", weight_expr)
        df_window = df_window.Define("weight", weight_expr)
        
        #DATAFRAME CR ZLL
        df_CR = df_CR.Define("weight", weight_expr)

        #DATAFRAME CR ZL
        df_CRZL = df_CRZL.Define("weight", weight_expr)

    regions = {
        "3P1F": df_CR.Filter("CRflag == 8388608"),
        "2P2F": df_CR.Filter("CRflag == 4194304"),
        "SS":   df_CR.Filter("CRflag == 2097152"),
        "SIP":  df_CR.Filter("CRflag == 21"),
        }


    histos = {}


    df_full = df_full.Define("Nj_clipped", "Nj >= 4 ? 4 : Nj")
    df_window = df_window.Define("Nj_clipped", "Nj >= 4 ? 4 : Nj")

    df_CRZL = df_CRZL.Define("Nj_clipped", "Nj >= 4 ? 4 : Nj")
    df_CRZL_window = df_CRZL.Filter("ZZMass > 105 && ZZMass < 160")


    # -------------------------
    # Control Region histograms
    # -------------------------
    for region, df_cr_region in regions.items():

        #select correct final state selection depending on SS or OS regions
        fs_dict = FINAL_STATES_SS if region in ["SS", "SIP"] else FINAL_STATES_OS
        #--------------------

        df_cr_region = df_cr_region.Define("Nj_clipped", "Nj >= 4 ? 4 : Nj")

        df_cr_window = df_cr_region.Filter("ZZMass > 105 && ZZMass < 160")

        # INCLUSIVE CR
        histos_cr_inc = define_histograms(
            df_cr_region,
            df_cr_window,
            f"{sample}_{region}",
            isMC
        )
        histos.update(histos_cr_inc)

        # CR + FINAL STATES
        for fs, cut in fs_dict.items():

            df_cr_fs = df_cr_region.Filter(cut)
            df_cr_fs_window = df_cr_fs.Filter("ZZMass > 105 && ZZMass < 160")

            histos_cr_fs = define_histograms(
                df_cr_fs,
                df_cr_fs_window,
                f"{sample}_{region}_{fs}",
                isMC
            )

            histos.update(histos_cr_fs)

    # -------------------------
    # ZL CR Incluive
    # -------------------------
    histos_CRZL = define_histograms(
        df_CRZL,
        df_CRZL_window,
        f"{sample}_CRZL",
        isMC
    )

    histos.update(histos_CRZL)



    # -------------------------
    # Inclusive histograms
    # -------------------------
    histos_inc = define_histograms(
        df_full,
        df_window,
        sample,
        isMC
    )

    histos.update(histos_inc)

    # -------------------------
    # SR by final state
    # -------------------------

    for fs, cut in FINAL_STATES_OS.items():

        df_fs_full = df_full.Filter(cut)
        df_fs_window = df_window.Filter(cut)

        histos_fs = define_histograms(
            df_fs_full,
            df_fs_window,
            f"{sample}_{fs}",
            isMC
        )

        histos.update(histos_fs)
    
    output_file.cd()

    for h in histos.values():
        hist = h.GetValue()
        hist.Scale(1.0, "width")
        hist.Write()

# =========================
# RUN MC
# =========================
def run_mc(period):
    path = MC_PATHS.get(period)
    if not path:
        raise ValueError(f"Unknown MC period: {period}")

    samples = ["DYJetsToLL","TTto2L2Nu", "WZto3LNu"]

    fout = ROOT.TFile.Open(f"H4l_MC_{period}_CHECKS_ZX.root", "RECREATE")
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

    histos = {}


    # -------------------------
    # Inclusive histograms
    # -------------------------
    histos_inc = define_histograms(
        df_full,
        df_window,
        "DATA",
        isMC=False
    )

    histos.update(histos_inc)

    # -------------------------
    # Final state histograms
    # -------------------------

    for fs, cut in FINAL_STATES.items():

        df_fs_full = df_full.Filter(cut)
        df_fs_window = df_window.Filter(cut)

        histos_fs = define_histograms(
            df_fs_full,
            df_fs_window,
            f"DATA_{fs}",
            isMC=False
        )

        histos.update(histos_fs)

    fout = ROOT.TFile.Open(f"H4l_Data_{period}_CHECKS_ZX.root", "RECREATE")
    for h in histos.values():
        hist = h.GetValue()
        hist.Scale(1.0, "width")
        hist.Write()

    fout.Close()

# =========================
# RUN ZX (special treatment)
# =========================

# =========================
# RUN ZX (corrected for final states)
# =========================
def run_zx(period):
    path = ZX_PATHS.get(period)
    if not path or not os.path.exists(path):
        print(f"[WARNING] ZX file not found for period {period}: {path}")
        return

    df = ROOT.RDataFrame("candTree", path)

    # Map columns to standard names for observables
    df = df.Define("pT4l", "ZZPt") \
           .Define("rapidity4l", "ZZy") \
           .Define("massZ1", "Z1Mass") \
           .Define("massZ2", "Z2Mass") \
           .Define("weight", "weight1")  # ZX weight

    # -------------------------
    # Inclusive histograms
    # -------------------------
    df_SR = df.Filter("ZZMass > 105 && ZZMass < 160")

    fout = ROOT.TFile.Open(f"H4l_ZX_{period}_CHECKS_ZX.root", "RECREATE")

    for obs, cfg in OBSERVABLES.items():
        bins = array.array('d', cfg["bins"])
        for suffix, rdf in [("_FULL", df), ("_105to160", df_SR)]:
            hname = f"{obs}_ZX{suffix}"
            h = rdf.Histo1D(
                ROOT.RDF.TH1DModel(hname, hname, len(bins)-1, bins),
                cfg["var"],
                "weight"
            ).GetValue()
            h.Scale(1.0, "width")
            h.Write()

    # -------------------------
    # Final state histograms (corrected)
    # -------------------------
    FINAL_STATES_ZX = {
        "4e":    "(abs(Z1Flav)==121) && (Z2Flav==121)",
        "4mu":   "(abs(Z1Flav)==169) && (Z2Flav==169)",
        "2e2mu": "(abs(Z1Flav)==121) && (Z2Flav==169)",
        "2mu2e": "(abs(Z1Flav)==169) && (Z2Flav==121)",
        "2e2mu_com": "((abs(Z1Flav)==121) && (Z2Flav==169)) || ((abs(Z1Flav)==169) && (Z2Flav==121))",
    }

    for fs, cut in FINAL_STATES_ZX.items():
        df_fs = df.Filter(cut)
        df_fs_SR = df_fs.Filter("ZZMass > 105 && ZZMass < 160")

        for obs, cfg in OBSERVABLES.items():
            bins = array.array('d', cfg["bins"])
            for suffix, rdf in [("_FULL", df_fs), ("_105to160", df_fs_SR)]:
                hname = f"{obs}_ZX_{fs}{suffix}"
                h = rdf.Histo1D(
                    ROOT.RDF.TH1DModel(hname, hname, len(bins)-1, bins),
                    cfg["var"],
                    "weight"
                ).GetValue()
                h.Scale(1.0, "width")
                h.Write()

    fout.Close()
    print(f"[INFO] ZX histograms for {period} written to H4l_ZX_{period}_DIFF.root")
# =========================
# MAIN
# =========================
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--periods", nargs="+", default=["2024"], help="List of periods/years to run over")#, "2022", "2022EE", "2023preBPix", "2023postBPix", "2024"
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
