#!/bin/env python3

import ROOT
import os

# =========================================================
# ROOT GLOBAL SETTINGS
# =========================================================
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2()

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadLeftMargin(0.14)
ROOT.gStyle.SetPadBottomMargin(0.13)
ROOT.gStyle.SetPadRightMargin(0.05)
ROOT.gStyle.SetPadTopMargin(0.08)
ROOT.gStyle.SetTitleSize(0.05, "XYZ")
ROOT.gStyle.SetLabelSize(0.045, "XYZ")
ROOT.gStyle.SetTitleOffset(1.1, "X")
ROOT.gStyle.SetTitleOffset(1.4, "Y")

canvasSizeX = 900
canvasSizeY = 700


OBSERVABLES = {

    "costhetaZ1": {
        "title": "cos(#theta_{Z1});cos(#theta_{Z1});Events/bin width"
    },

    "costhetaZ2": {
        "title": "cos(#theta_{Z2});cos(#theta_{Z2});Events/bin width"
    },

    "costhetastar": {
        "title": "cos(#theta*);cos(#theta*);Events/bin width"
    },

    "phi": {
        "title": "#phi;#phi;Events/bin width"
    },

    "phi1": {
        "title": "#phi_{1};#phi_{1};Events/bin width"
    },

    "pT4l": {
        "title": "p_{T}^{4l};p_{T}^{4l} [GeV];Events/bin width"
    },

    "rapidity4l": {
        "title": "|y_{4l}|;|y_{4l}|;Events/bin width"
    },

    "Nj": {
        "title": "N_{jets};N_{jets};Events/bin width"
    },

    "dphijj": {
        "title": "#Delta#phi_{jj};#Delta#phi_{jj};Events/bin width"
    },

    "massZ1": {
        "title": "m_{Z1};m_{Z1} [GeV];Events/bin width"
    },

    "massZ2": {
        "title": "m_{Z2};m_{Z2} [GeV];Events/bin width"
    },

    "pTj1": {
        "title": "p_{T}^{j1};p_{T}^{j1} [GeV];Events/bin width"
    },

    "pTj2": {
        "title": "p_{T}^{j2};p_{T}^{j2} [GeV];Events/bin width"
    },

    "mjj": {
        "title": "m_{jj};m_{jj} [GeV];Events/bin width"
    },

    "absdetajj": {
        "title": "|#Delta#eta_{jj}|;|#Delta#eta_{jj}|;Events/bin width"
    },

    "pTHj": {
        "title": "p_{T}^{Hj};p_{T}^{Hj} [GeV];Events/bin width"
    },

    "TCjMax": {
        "title": "T_{Cj};T_{Cj} [GeV];Event/bin width"
    },
    
    "TBjMax": {
        "title": "T_{Bj};T_{Bj} [GeV];Events/bin width"
    },
}

# =========================================================
# COLORS
# =========================================================
category_colors = {
    "ZX": "#669966",
    "EW": "#0331B9",
    "qqZZ": "#99ccff",
    "ggZZ": "#4b78ff",
    "H125": "#ff9999",
}

FINAL_STATES = [
    "inclusive",
    "4e",
    "4mu",
    "2e2mu",
    "2mu2e",
]

# =========================================================
# SAMPLE → CATEGORY MAP (MC ONLY)
# =========================================================
sample_map = {
    "WWZ": "EW",
    "WZZ": "EW",
    "ZZZ": "EW",

    "ZZTo4l": "qqZZ",

    "ggTo4mu_Contin_MCFM701": "ggZZ",
    "ggTo4e_Contin_MCFM701": "ggZZ",
    "ggTo4tau_Contin_MCFM701": "ggZZ",
    "ggTo2e2mu_Contin_MCFM701": "ggZZ",
    "ggTo2e2tau_Contin_MCFM701": "ggZZ",
    "ggTo2mu2tau_Contin_MCFM701": "ggZZ",

    "VBFH125": "H125",
    "ggH125": "H125",
    "WplusH125": "H125",
    "WminusH125": "H125",
    "ZH125": "H125",
    "ttH125": "H125",
}

# =========================================================
# BUILD MC HISTOGRAMS
# =========================================================
def build_mc(files_dict, variable, suffix, fs):

    summed = {"EW": None, "ggZZ": None, "qqZZ": None, "H125": None}

    for period, rootfile in files_dict.items():

        print(f"[INFO] Adding MC period: {period}")

        for sample, category in sample_map.items():

            #hist_name = f"{variable}_{sample}_{suffix}"
            if fs == "inclusive":
                hist_name = f"{variable}_{sample}_{suffix}"
            else:
                hist_name = f"{variable}_{sample}_{fs}_{suffix}"

            hist = rootfile.Get(hist_name)

            if not hist:
                continue

            clone = hist.Clone(f"{hist_name}_{period}")

            if summed[category] is None:
                summed[category] = clone
            else:
                summed[category].Add(clone)

    # Style
    for cat, hist in summed.items():
        if hist:
            hist.SetFillColor(ROOT.TColor.GetColor(category_colors[cat]))
            hist.SetLineColor(ROOT.kBlack)

    return summed


# =========================================================
# SUM ZX BACKGROUND
# =========================================================
def sum_zx(files_dict, variable, suffix, fs):

    hsum = None

    for period, rootfile in files_dict.items():

        print(f"[INFO] Adding ZX period: {period}")

        if fs == "inclusive":
            hist_name = f"{variable}_ZX_{suffix}"
        else:
            hist_name = f"{variable}_ZX_{fs}_{suffix}"

        hist = rootfile.Get(hist_name)

        if not hist:
            continue

        clone = hist.Clone(f"{hist_name}_{period}")

        if hsum is None:
            hsum = clone
        else:
            hsum.Add(clone)

    if hsum:
        hsum.SetFillColor(ROOT.TColor.GetColor(category_colors["ZX"]))
        hsum.SetLineColor(ROOT.kBlack)

    return hsum


# =========================================================
# SUM DATA
# =========================================================
def sum_data(files_dict, variable, suffix, fs):

    hsum = None

    for period, rootfile in files_dict.items():

        print(f"[INFO] Adding DATA period: {period}")

        if fs == "inclusive":
            hist_name = f"{variable}_DATA_{suffix}"
        else:
            hist_name = f"{variable}_DATA_{fs}_{suffix}"
        hist = rootfile.Get(hist_name)

        if not hist:
            continue

        clone = hist.Clone(f"{hist_name}_{period}")

        if hsum is None:
            hsum = clone
        else:
            hsum.Add(clone)

    if hsum:
        hsum.SetMarkerStyle(20)
        hsum.SetMarkerSize(1.2)
        hsum.SetLineColor(ROOT.kBlack)

    return hsum


# =========================================================
# DRAW PLOT
# =========================================================
def draw_plot(mc_hists, hzx, hdata, variable, suffix, fs, outdir):

    title = OBSERVABLES[variable]["title"]

    stack = ROOT.THStack(
        f"stack_{variable}_{suffix}",
        title
    )

    if hzx:
        stack.Add(hzx)

    # stack order
    order = ["EW", "ggZZ", "qqZZ"]

    for cat in order:
        if mc_hists[cat]:
            stack.Add(mc_hists[cat])

    if mc_hists["H125"]:
        stack.Add(mc_hists["H125"])

    canvas = ROOT.TCanvas(
        f"{variable}_{suffix}",
        f"{variable}_{suffix}",
        canvasSizeX,
        canvasSizeY
    )

    canvas.SetTicks()

    stack.Draw("HIST")
    stack.SetMaximum(stack.GetMaximum() * 1.5)

    if hdata:
        hdata.Draw("E SAME")

    legend = ROOT.TLegend(0.65, 0.65, 0.88, 0.88)

    if mc_hists["H125"]:
        legend.AddEntry(mc_hists["H125"], "H(125)", "f")

    if hzx:
        legend.AddEntry(hzx, "Z+X", "f")

    if mc_hists["qqZZ"]:
        legend.AddEntry(mc_hists["qqZZ"], "q#bar{q} → ZZ", "f")

    if mc_hists["ggZZ"]:
        legend.AddEntry(mc_hists["ggZZ"], "gg → ZZ", "f")

    if mc_hists["EW"]:
        legend.AddEntry(mc_hists["EW"], "EW", "f")

    if hdata:
        legend.AddEntry(hdata, "Data", "pe")

    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.Draw()

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)

    latex.DrawLatex(0.15, 0.92, "CMS Preliminary")
    latex.DrawLatex(0.58, 0.92, f"{suffix}  (13.6 TeV)")

    canvas.SaveAs(f"{outdir}/{variable}_{fs}_{suffix}.png")

    print(f"[INFO] Saved {outdir}/{variable}_{suffix}.png")


# =========================================================
# MAIN
# =========================================================
if __name__ == "__main__":

    periods = ["2024"]
    #"2022", "2022EE", "2023preBPix","2023postBPix", "2024"

    # MC files
    mc_files = {}
    for p in periods:

        fname = f"H4l_MC_{p}_DIFF.root"

        if os.path.exists(fname):
            mc_files[p] = ROOT.TFile.Open(fname)

    # ZX files
    zx_files = {}
    for p in periods:

        fname = f"H4l_ZX_{p}_DIFF.root"

        if os.path.exists(fname):
            zx_files[p] = ROOT.TFile.Open(fname)

    # DATA files
    data_files = {}
    for p in periods:

        fname = f"H4l_Data_{p}_DIFF.root"

        if os.path.exists(fname):
            data_files[p] = ROOT.TFile.Open(fname)

    outdir = "plots_Obs_byfinal_state_2024"
    os.makedirs(outdir, exist_ok=True)

    for variable in OBSERVABLES:

        print("\n====================================")
        print(f"[INFO] Observable: {variable}")
        print("====================================")

        for fs in FINAL_STATES:

            print(f"\n[INFO] Final state: {fs}")

            for suffix in ["FULL", "105to160"]:

                print(f"[INFO] Region: {suffix}")

                mc_hists = build_mc(mc_files, variable, suffix, fs)

                hzx = sum_zx(zx_files, variable, suffix, fs)

                hdata = sum_data(data_files, variable, suffix, fs)

                draw_plot(mc_hists, hzx, hdata, variable, suffix, fs, outdir)