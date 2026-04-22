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
    p: ROOT.TFile.Open(f"H4l_MCplusDATA_{year}.root")
    for p in periods
    if os.path.exists(f"H4l_MCplusDATA_{year}.root")
}

# =========================
# LOAD DATA FILE  (NEW)
# =========================
data_file = None
data_path = f"H4l_MCplusDATA_{year}.root"

if os.path.exists(data_path):
    data_file = ROOT.TFile.Open(data_path)
    print(f"[INFO] Loaded data file: {data_path}")
else:
    print(f"[WARNING] Data file not found: {data_path}")

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
# DATA HELPER (NEW)
# =========================
def get_data_hist(prefix):
    if not data_file:
        return None

    # try direct match first
    name = f"{variable}_DATA_{prefix}"
    h = data_file.Get(name)

    if h:
        h = h.Clone()
    else:
        # fallback: sometimes suffix issues happen
        keys = [k.GetName() for k in data_file.GetListOfKeys()]
        match = next((k for k in keys if k.startswith(f"{variable}_DATA_{prefix}")), None)
        if not match:
            return None
        h = data_file.Get(match).Clone()

    h.SetMarkerStyle(20)
    h.SetMarkerSize(1.1)
    h.SetLineColor(ROOT.kBlack)
    h.SetMarkerColor(ROOT.kBlack)
    return h


# =========================
# BUILD STACK
# =========================
def build_stack(prefix):

    stack = ROOT.THStack(f"stack_{prefix}", "N_{jets};N_{jets};Events")
    hist_list = []

    for sample in samples:

        hist_name = f"{variable}_{sample}_{prefix}"
        print("hist_name", hist_name)
        h = get_hist_sum(hist_name)

        print(f"Got histogram for {sample} with name {hist_name}: {h}")

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

    if tokens[0] in control_regions + ["CRZL"]:
        region = tokens[0]
        tokens = tokens[1:]

    for t in tokens:
        if t in final_states and t != "inclusive":
            fs = t

        # ADD THIS
        if t in ["Ze", "Zmu"]:
            fs = t

    if "105to160" in tokens:
        window = "105-160"

    return region, fs, window


# =========================
# DRAW FUNCTION (MODIFIED ONLY HERE)
# =========================
def draw(prefix, outname):

    print("prefix:", prefix)
    print("outname:", outname)

    c = ROOT.TCanvas("c", "c", 900, 700)

    # build histos
    
    stack, hist_list = build_stack(prefix)
    print("hist_list", hist_list)
    hdata = get_data_hist(prefix)

    if not hist_list:
        print(f"[WARNING] No histograms for {prefix}")
        return

    # =========================
    # COMPUTE MAX (MC + DATA)
    # =========================
    mc_max = max([h.GetMaximum() for h in hist_list]) if hist_list else 0
    data_max = hdata.GetMaximum() if hdata else 0

    ymax = max(mc_max, data_max) * 1.3

    # =========================
    # DRAW
    # =========================
    stack.Draw("HIST")
    stack.SetMaximum(ymax)

    if hdata:
        hdata.Draw("E1 SAME")

    # =========================
    # LEGEND
    # =========================
    leg = ROOT.TLegend(0.65, 0.65, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    for h in hist_list:
        sample_name = h.GetName().split("_")[1]
        leg.AddEntry(h, sample_name, "f")

    if hdata:
        leg.AddEntry(hdata, "Data", "lep")

    leg.Draw()

    # =========================
    # LABEL
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
outdir = f"plots_{variable}_MCplusData_only_2024"
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

    for suffix in suffixes:
        prefix = f"{cr}_{suffix}"
        draw(prefix, f"{outdir}/{variable}_{prefix}.png")

    for fs in final_states:
        if fs == "inclusive":
            continue

        for suffix in suffixes:
            prefix = f"{cr}_{fs}_{suffix}"
            draw(prefix, f"{outdir}/{variable}_{prefix}.png")


# =========================================================
# 3) CRZL
# =========================================================
for suffix in ["FULL"]:
    
    print("Drawing CRZL FULL only...")

    prefix = "CRZL_FULL"
    draw(prefix, f"{outdir}/{variable}_CRZL_FULL.png")


    prefix_ze = f"CRZL_Ze_{suffix}"
    draw(prefix_ze, f"{outdir}/{variable}_CRZL_Ze_{suffix}.png")

    prefix_zmu = f"CRZL_Zmu_{suffix}"
    draw(prefix_zmu, f"{outdir}/{variable}_CRZL_Zmu_{suffix}.png")