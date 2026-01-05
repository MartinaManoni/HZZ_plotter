import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from collections import defaultdict
from scipy.stats import chi2
import ROOT

plt.style.use(hep.style.CMS)

YEAR = "2024"
CHANNEL = "e"  # "e", "mu", or "incl"
lumi = 108.8  # 1/fb

# ----------------------------
# MC and Data files
# ----------------------------
DY_files = [
    # DY -> 2Tau
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
]

TT_files = [
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/TTto2L2Nu/TTto2L2Nu_0/TTto2L2Nu/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/TTto2L2Nu/TTto2L2Nu_1/TTto2L2Nu/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/TTto2L2Nu/TTto2L2Nu_2/TTto2L2Nu/ZZ4lAnalysis.root",
        "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_MC/TTto2L2Nu/TTto2L2Nu_3/TTto2L2Nu/ZZ4lAnalysis.root"
]

data_files = [
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
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/2024_Data/MuonEG2024Iv2v2/ZZ4lAnalysis.root"
]

# ----------------------------
# Helper functions
# ----------------------------
def select_channel(zflav, channel):
    zf = abs(zflav)
    if channel == "mu": sel = (zf == 169)
    elif channel == "e": sel = (zf == 121)
    elif channel == "incl": sel = (zf == 169) | (zf == 121)
    else: raise ValueError("CHANNEL must be 'mu', 'e', or 'incl'")
    return ak.flatten(sel)

def poisson_interval(k, alpha=0.317):
    low, high = chi2.ppf(alpha/2, 2*k)/2, chi2.ppf(1-alpha/2, 2*k+2)/2
    return k-low, high-k

def get_genEventSumw(filenames):
    total_sumw = 0.0
    for f in filenames:
        tf = ROOT.TFile.Open(f)
        runs = tf.Runs
        for iRun in range(runs.GetEntries()):
            runs.GetEntry(iRun)
            total_sumw += runs.genEventSumw
        tf.Close()
    return total_sumw

def read_mc_arrays(filenames, branches):
    arrays = defaultdict(list)
    # Use a dict to explicitly point to the TTree
    file_tree_map = {fn: "Events" for fn in filenames}  # all files, read 'Events' TTree
    for chunk in uproot.iterate(file_tree_map, branches=branches, step_size="100 MB"):
        for br in branches:
            arrays[br].append(chunk[br])
    for br in arrays:
        arrays[br] = ak.concatenate(arrays[br])
    return arrays

# ----------------------------
# Read MC arrays
# ----------------------------
DY_branches = ["bestZIdx", "ZCand_mass", "ZCand_flav", "overallEventWeight"]
TT_branches = DY_branches

dy = read_mc_arrays(DY_files, DY_branches)
tt = read_mc_arrays(TT_files, TT_branches)

# Data arrays (merge multiple files if needed)
dt_zc_list, dt_zm_list, dt_z1f_list = [], [], []
for f in data_files:
    with uproot.open(f) as tf:
        dt_zc_list.append(tf['Events/bestZIdx'].array())
        dt_zm_list.append(tf['Events/ZCand_mass'].array())
        dt_z1f_list.append(tf['Events/ZCand_flav'].array())
dt_zc = ak.concatenate(dt_zc_list)
dt_zm = ak.concatenate(dt_zm_list)
dt_z1f = ak.concatenate(dt_z1f_list)

# ----------------------------
# Sum genEventSumw
# ----------------------------
dy_cnt = get_genEventSumw(DY_files)
tt_cnt = get_genEventSumw(TT_files)

# ----------------------------
# Prepare MC sample
# ----------------------------
def prepare_sample(arrays, cnt):
    m = arrays["ZCand_mass"][(ak.singletons(arrays["bestZIdx"]))]
    zf = arrays["ZCand_flav"][(ak.singletons(arrays["bestZIdx"]))]
    w  = arrays["overallEventWeight"]
    sel = select_channel(zf, CHANNEL)
    return m[sel], w[sel] * 1000 * lumi / cnt

m_dy, w_dy = prepare_sample(dy, dy_cnt)
m_tt, w_tt = prepare_sample(tt, tt_cnt)

# Data selection
dt_z1mass = dt_zm[(ak.singletons(dt_zc))]
dt_z1mass = dt_z1mass[select_channel(dt_z1f[(ak.singletons(dt_zc))], CHANNEL)]

# ----------------------------
# Histogramming
# ----------------------------
CP_BINNING = np.linspace(60, 120, 100)
n_dy, bins = np.histogram(ak.flatten(m_dy), weights=w_dy, bins=CP_BINNING)
n_tt, _ = np.histogram(ak.flatten(m_tt), weights=w_tt, bins=CP_BINNING)
n_dt, _ = np.histogram(ak.flatten(dt_z1mass), bins=CP_BINNING)
mc_tot = n_dy + n_tt
binsc = 0.5*(bins[1:] + bins[:-1])

# ----------------------------
# Plot
# ----------------------------
fig, ax = plt.subplots(2,1, figsize=(14,10), gridspec_kw={"height_ratios":(3,1)}, sharex=True)
ax[0].step(binsc, mc_tot*(sum(n_dt)/sum(mc_tot)), where='mid', label=r'MC (DY + TT)')
ax[0].errorbar(binsc, n_dt, yerr=poisson_interval(n_dt), fmt='.', color='k', label='Data')
ax[0].set_ylabel('Events / bin width')
ax[0].legend()

ax[1].errorbar(binsc, mc_tot*(sum(n_dt)/sum(mc_tot))/n_dt, fmt='.', color='k')
ax[1].set_ylim(0.5,1.5)
xlabel = r"$m_{ee}$" if CHANNEL=="e" else r"$m_{\mu\mu}$" if CHANNEL=="mu" else r"$m_{\ell\ell}$"
ax[1].set_xlabel(xlabel)
ax[1].set_ylabel('MC/Data')

hep.cms.label(llabel="Preliminary", data=True, lumi=round(lumi,2), com="13.6", ax=ax[0])

plt.savefig(f"plot_{CHANNEL}_2024.png")
plt.savefig(f"plot_{CHANNEL}_2024.pdf")
plt.close()
