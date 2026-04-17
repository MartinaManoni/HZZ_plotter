#!/bin/env python3

import ROOT
import os

ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2()

# =========================
# CONFIG
# =========================
variable = "Z2"
suffixes = ["FULL", "105to160"]

final_states = [
    "inclusive",
    "4e",
    "4mu",
    "2e2mu",
    "2mu2e",
    "2e2mu_com",
]

control_regions = ["3P1F", "2P2F", "SS", "SIP"]

samples = ["DYJetsToLL", "TTto2L2Nu", "WZto3LNu"]

year = "2024"

colors = {
    "DYJetsToLL": ROOT.kAzure+1,
    "TTto2L2Nu": ROOT.kRed+1,
    "WZto3LNu": ROOT.kGreen+2,
}

# =========================
# LOAD MC FILES
# =========================
periods = ["2024"]

mc_files = {
    p: ROOT.TFile.Open(f"H4l_MC_{p}_CHECKS_ZX.root")
    for p in periods
    if os.path.exists(f"H4l_MC_{p}_CHECKS_ZX.root")
}

# =========================
# HELPERS
# =========================
def get_hist_sum(hist_name):
    hsum = None

    for f in mc_files.values():
        h = f.Get(hist_name)
        if not h:
            continue

        h = h.Clone()

        if hsum is None:
            hsum = h
        else:
            hsum.Add(h)

    return hsum


# =========================
# BUILD STACK
# =========================
def build_stack(prefix):

    stack = ROOT.THStack(f"stack_{prefix}", "N_{jets};N_{jets};Events")
    hist_list = []

    for sample in samples:

        hist_name = f"{variable}_{sample}_{prefix}"
        h = get_hist_sum(hist_name)

        if not h:
            continue

        h.SetFillColor(colors[sample])
        h.SetLineColor(ROOT.kBlack)

        stack.Add(h)
        hist_list.append(h)

    return stack, hist_list


# =========================
# PARSE LABELS
# =========================
def parse_prefix(prefix):

    region = "SR"
    fs = "inclusive"
    window = "FULL"

    tokens = prefix.split("_")

    # detect CR / CRZL / etc
    if tokens[0] in control_regions + ["CRZL"]:
        region = tokens[0]
        tokens = tokens[1:]

    # detect final state
    for t in tokens:
        if t in final_states and t != "inclusive":
            fs = t

    # detect window
    if "105to160" in tokens:
        window = "105-160"

    return region, fs, window


# =========================
# DRAW FUNCTION
# =========================
def draw(prefix, outname):

    c = ROOT.TCanvas("c", "c", 900, 700)

    stack, hist_list = build_stack(prefix)

    if not hist_list:
        print(f"[WARNING] No histograms for {prefix}")
        return

    stack.Draw("HIST")
    stack.SetMaximum(stack.GetMaximum() * 1.5)

    # =========================
    # LEGEND
    # =========================
    leg = ROOT.TLegend(0.65, 0.65, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    for h in hist_list:
        sample_name = h.GetName().split("_")[1]
        leg.AddEntry(h, sample_name, "f")

    leg.Draw()

    # =========================
    # LABEL (Region / FS / Year)
    # =========================
    region, fs, window = parse_prefix(prefix)

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.035)

    latex.DrawLatex(0.15, 0.85, f"{region} region")
    latex.DrawLatex(0.15, 0.80, f"{fs}")
    latex.DrawLatex(0.15, 0.75, f"{window}")
    latex.DrawLatex(0.15, 0.70, f"{year}")

    c.SaveAs(outname)
    print(f"[INFO] Saved {outname}")


# =========================
# OUTPUT DIR
# =========================
outdir = f"plots_{variable}_MC_only_2024"
os.makedirs(outdir, exist_ok=True)


# =========================================================
# 1) SIGNAL REGION (SR)
# =========================================================
for fs in final_states:
    for suffix in suffixes:

        if fs == "inclusive":
            prefix = suffix
        else:
            prefix = f"{fs}_{suffix}"

        draw(prefix, f"{outdir}/{variable}_SR_{prefix}.png")


# =========================================================
# 2) CONTROL REGIONS (CR)
# =========================================================
for cr in control_regions:

    # inclusive CR
    for suffix in suffixes:
        prefix = f"{cr}_{suffix}"
        draw(prefix, f"{outdir}/{variable}_{prefix}.png")

    # CR + final states
    for fs in final_states:
        if fs == "inclusive":
            continue

        for suffix in suffixes:
            prefix = f"{cr}_{fs}_{suffix}"
            draw(prefix, f"{outdir}/{variable}_{prefix}.png")


# =========================================================
# 3) CRZL (inclusive only)
# =========================================================
for suffix in suffixes:
    prefix = f"CRZL_{suffix}"
    draw(prefix, f"{outdir}/{variable}_CRZL_{suffix}.png")