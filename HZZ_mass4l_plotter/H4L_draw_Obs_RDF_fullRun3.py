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
def build_mc(files_dict, variable, suffix):

    summed = {"EW": None, "ggZZ": None, "qqZZ": None, "H125": None}

    for period, rootfile in files_dict.items():

        print(f"[INFO] Adding MC period: {period}")

        for sample, category in sample_map.items():

            hist_name = f"{variable}_{sample}_{suffix}"
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
def sum_zx(files_dict, variable, suffix):

    hsum = None

    for period, rootfile in files_dict.items():

        print(f"[INFO] Adding ZX period: {period}")

        hist_name = f"{variable}_ZX_{suffix}"
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
def sum_data(files_dict, variable, suffix):

    hsum = None

    for period, rootfile in files_dict.items():

        print(f"[INFO] Adding DATA period: {period}")

        hist_name = f"{variable}_DATA_{suffix}"
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
def draw_plot(mc_hists, hzx, hdata, variable, suffix, outdir):

    stack = ROOT.THStack(
        f"stack_{variable}_{suffix}",
        f";rapidity;Events"
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

    canvas.SaveAs(f"{outdir}/{variable}_{suffix}.png")

    print(f"[INFO] Saved {outdir}/{variable}_{suffix}.png")


# =========================================================
# MAIN
# =========================================================
if __name__ == "__main__":

    periods = ["2022", "2022EE"]
    variable = "rapidity4l"

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

    outdir = "plots_allPeriods"
    os.makedirs(outdir, exist_ok=True)

    for suffix in ["FULL", "105to160"]:

        print("\n==============================")
        print(f"[INFO] Making plot: {suffix}")
        print("==============================")

        mc_hists = build_mc(mc_files, variable, suffix)

        hzx = sum_zx(zx_files, variable, suffix)

        hdata = sum_data(data_files, variable, suffix)

        draw_plot(mc_hists, hzx, hdata, variable, suffix, outdir)