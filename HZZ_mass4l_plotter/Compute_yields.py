#!/bin/env python3
import ROOT
import os

ROOT.EnableImplicitMT()

# -----------------------------
# Paths and constants
# -----------------------------
MC_PATHS = {
    "2022": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2022_MC/",
    "2022EE": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2022EE_MC/",
    "2023preBPix": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2023preBPix_MC/",
    "2023postBPix": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2023postBPix_MC/",
    "2024": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/",
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

FINAL_STATES = ["4mu", "4e", "2e2mu", "2mu2e", "4l"]

QQZZ_SAMPLES = ["ZZTo4l"]
GGZZ_SAMPLES = [
    "ggTo4mu_Contin_MCFM701",
    "ggTo4e_Contin_MCFM701",
    "ggTo2e2mu_Contin_MCFM701",
    "ggTo4tau_Contin_MCFM701",
    "ggTo2e2tau_Contin_MCFM701",
    "ggTo2mu2tau_Contin_MCFM701",
]
SIGNAL_SAMPLES = ["ggH125", "VBFH125", "WplusH125", "WminusH125", "ZH125", "ttH125"]

# Flavor mapping for 2e2mu <-> 2mu2e
Z_FLAVORS = {
    "4mu": (-169, -169),
    "4e": (-121, -121),
    "2e2mu": [(-169, -121)],
    "2mu2e": [(-121, -169)],
}

# -----------------------------
# Helpers
# -----------------------------
def get_genEventSumw(f):
    runs = f.Runs
    nRuns = runs.GetEntries()
    genEventSumw = 0.
    for iRun in range(nRuns):
        runs.GetEntry(iRun)
        genEventSumw += runs.genEventSumw
    return genEventSumw

def compute_yields(mc_files, year, zx_file=None, m4l_min=70, m4l_max=3000):
    yields = {}
    lumi = LUMI[year]

    for sample, filepath in mc_files.items():
        if not os.path.exists(filepath):
            continue

        df = ROOT.RDataFrame("Events", filepath)
        df = df.Filter("bestCandIdx != -1").Filter("HLT_passZZ4l")
        df = df.Define("m4l", "ZZCand_mass[bestCandIdx]") \
               .Define("Z1flav", "ZZCand_Z1flav[bestCandIdx]") \
               .Define("Z2flav", "ZZCand_Z2flav[bestCandIdx]")

        f = ROOT.TFile.Open(filepath)
        genEventSumw = get_genEventSumw(f)
        f.Close()

        df = df.Define("weight", f"overallEventWeight * ZZCand_dataMCWeight * {lumi} * 1000. / {genEventSumw}")
        df = df.Filter(f"m4l > {m4l_min} && m4l < {m4l_max}")

        sample_yields = {}
        # Per final state
        for fs in FINAL_STATES:
            if fs == "4l":
                sample_yields[fs] = df.Sum("weight").GetValue()
            else:
                flavs = Z_FLAVORS[fs]
                if isinstance(flavs[0], tuple):
                    selection = " || ".join([f"(Z1flav == {f[0]} && Z2flav == {f[1]})" for f in flavs])
                else:
                    selection = f"(Z1flav == {flavs[0]} && Z2flav == {flavs[1]})"
                sample_yields[fs] = df.Filter(selection).Sum("weight").GetValue()

        yields[sample] = sample_yields

    # ZX
    if zx_file and os.path.exists(zx_file):
        df_zx = ROOT.RDataFrame("candTree", zx_file)
        df_zx = df_zx.Define("m4l", "ZZMass").Define("finState", "FinState").Define("weight", "weight1")
        df_zx = df_zx.Filter(f"m4l > {m4l_min} && m4l < {m4l_max}")
        zx_yields = {}
        for fs, val in zip(["4mu","4e","2e2mu","2mu2e","4l"], [0,1,2,3,None]):
            if val is None:
                zx_yields[fs] = df_zx.Sum("weight").GetValue()
            else:
                zx_yields[fs] = df_zx.Filter(f"finState == {val}").Sum("weight").GetValue()
        yields["ZX"] = zx_yields

    return yields

# -----------------------------
# Run per year and aggregate
# -----------------------------
ALL_YIELDS = {}
YEARS = ["2022","2022EE","2023preBPix","2023postBPix","2024"]

for year in YEARS:
    print(f"\n====================== {year} ======================")
    path = MC_PATHS[year]
    zx_file = ZX_PATHS.get(year)

    samples = [
        "WWZ", "WZZ", "ZZZ",
        "ggTo4mu_Contin_MCFM701", "ggTo4e_Contin_MCFM701",
        "ggTo4tau_Contin_MCFM701", "ggTo2e2mu_Contin_MCFM701",
        "ggTo2e2tau_Contin_MCFM701", "ggTo2mu2tau_Contin_MCFM701",
        "ZZTo4l",
        "VBFH125", "ggH125", "WplusH125", "WminusH125",
        "ZH125", "ttH125"
    ]

    mc_files = {}
    for s in samples:
        filename = os.path.join(path, s, "ZZ4lAnalysis.root")
        if os.path.exists(filename):
            mc_files[s] = filename

    yields = compute_yields(mc_files, year, zx_file)
    ALL_YIELDS[year] = yields

    # Print per year per final state
    print("Sample           4mu      4e       2e2mu    2mu2e    4l")
    for sample, ydict in yields.items():
        print(f"{sample:15s} " + " ".join(f"{ydict[fs]:8.2f}" for fs in FINAL_STATES))

# -----------------------------
# Aggregated results
# -----------------------------
def aggregate_years(year_list):
    agg = {fs:{"qqZZ":0,"ggZZ":0,"ZX":0,"Signal":0,"Bkg":0} for fs in FINAL_STATES}
    for year in year_list:
        ydict = ALL_YIELDS.get(year, {})
        for fs in FINAL_STATES:
            qqZZ_total = sum(ydict.get(s, {}).get(fs,0) for s in QQZZ_SAMPLES)
            ggZZ_total = sum(ydict.get(s, {}).get(fs,0) for s in GGZZ_SAMPLES)
            signal_total = sum(ydict.get(s, {}).get(fs,0) for s in SIGNAL_SAMPLES)
            zx_total = ydict.get("ZX", {}).get(fs,0)
            bkg_total = qqZZ_total + ggZZ_total + zx_total
            agg[fs]["qqZZ"] += qqZZ_total
            agg[fs]["ggZZ"] += ggZZ_total
            agg[fs]["ZX"] += zx_total
            agg[fs]["Signal"] += signal_total
            agg[fs]["Bkg"] += bkg_total
    return agg

# Print custom aggregates
CUSTOM_AGG = {
    "2022+2022EE": ["2022","2022EE"],
    "2022+2022EE+2023preBPix+2023postBPix": ["2022","2022EE","2023preBPix","2023postBPix"]
}

for label, years_list in CUSTOM_AGG.items():
    agg = aggregate_years(years_list)
    print(f"\n===== Aggregated Yields: {label} =====")
    print("FS      qqZZ    ggZZ     ZX       Bkg      Signal")
    for fs in FINAL_STATES:
        print(f"{fs:6s} " + " ".join(f"{agg[fs][cat]:8.2f}" for cat in ["qqZZ","ggZZ","ZX","Bkg","Signal"]))

# -----------------------------
# Total over all years
# -----------------------------
agg_total = aggregate_years(YEARS)
print(f"\n===== Total over all years =====")
print("FS      qqZZ    ggZZ     ZX       Bkg      Signal")
for fs in FINAL_STATES:
    print(f"{fs:6s} " + " ".join(f"{agg_total[fs][cat]:8.2f}" for cat in ["qqZZ","ggZZ","ZX","Bkg","Signal"]))
