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
    "EW": "#0331B9",
    "qqZZ": "#99ccff",
    "ggZZ": "#4b78ff",
    "H125": "#ff9999",
}

# =========================================================
# SAMPLE → CATEGORY MAP
# =========================================================
sample_map = {
    "WWZ": "EW", "WZZ": "EW", "ZZZ": "EW",
    "ZZTo4l": "qqZZ",
    "ggTo4mu_Contin_MCFM701": "ggZZ",
    "ggTo4e_Contin_MCFM701": "ggZZ",
    "ggTo4tau_Contin_MCFM701": "ggZZ",
    "ggTo2e2mu_Contin_MCFM701": "ggZZ",
    "ggTo2e2tau_Contin_MCFM701": "ggZZ",
    "ggTo2mu2tau_Contin_MCFM701": "ggZZ",
    "VBFH125": "H125", "ggH125": "H125",
    "WplusH125": "H125", "WminusH125": "H125",
    "ZH125": "H125", "ttH125": "H125",
}

# =========================================================
# BUILD MC STACK
# =========================================================
def build_stack(files_dict, variable, suffix):

    summed = {cat: None for cat in category_colors.keys()}

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

    # Set appearance
    for cat, hist in summed.items():
        if hist:
            hist.SetFillColor(ROOT.TColor.GetColor(category_colors[cat]))
            hist.SetLineColor(ROOT.kBlack)

    stack = ROOT.THStack(
        f"stack_{variable}_{suffix}",
        f";cos(#theta_{{Z1}});Events"
    )

    for cat in ["EW", "ggZZ", "qqZZ", "H125"]:
        if summed[cat]:
            stack.Add(summed[cat])

    return stack, [summed[cat] for cat in ["EW", "ggZZ", "qqZZ", "H125"]]

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
def draw_plot(stack, hlist, hdata, variable, suffix, outdir):

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
    legend.AddEntry(hlist[3], "H(125)", "f")
    legend.AddEntry(hlist[2], "q#bar{q} #rightarrow ZZ", "f")
    legend.AddEntry(hlist[1], "gg #rightarrow ZZ", "f")
    legend.AddEntry(hlist[0], "EW", "f")
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

    periods = ["2022", "2022EE", "2023preBPix", "2023postBPix", "2024"]
    variable = "costhetaZ1"

    # Open MC files
    mc_files = {}
    for p in periods:
        fname = f"H4l_MC_{p}_DIFF.root"
        if os.path.exists(fname):
            mc_files[p] = ROOT.TFile.Open(fname)

    # Open DATA files
    data_files = {}
    for p in periods:
        fname = f"H4l_Data_{p}_DIFF.root"
        if os.path.exists(fname):
            data_files[p] = ROOT.TFile.Open(fname)

    outdir = "plots_allPeriods"
    os.makedirs(outdir, exist_ok=True)

    # Produce BOTH versions
    for suffix in ["FULL", "105to160"]:

        print(f"\n==============================")
        print(f"[INFO] Making plot: {suffix}")
        print(f"==============================")

        stack, hlist = build_stack(mc_files, variable, suffix)
        hdata = sum_data(data_files, variable, suffix)

        draw_plot(stack, hlist, hdata, variable, suffix, outdir)