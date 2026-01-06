#by Matteo Bonanomi
import uproot
import numpy as np
import awkward as ak
from tqdm import tqdm
from collections import defaultdict
import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.patches import Patch
import matplotlib.patches as patches
from scipy.stats import chi2
import ROOT

plt.style.use(hep.style.CMS)

def get_files_for_year(YEAR):

    base = " /eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/031125/" #with corrections
    #base = "/eos/user/m/mmanoni/HZZ_prod_241125_NoEleMuo_YesJetCorr/" #without corrections
    CHANNEL = "mu"  # "mu", "e", "incl"
    # ----------------------------
    # 2022 preEE
    # ----------------------------
    if YEAR == "2022preEE":
        lumi = 7.98
        #mc_dir   = f"{base}/2022_MC/PROD_samplesNano_2022_MC_cc84ce40/"
        mc_dir   = f"{base}/2022_MC/"
        #data_dir = f"{base}/2022_Data/PROD_samplesNano_2022_Data_cc84ce40/"
        data_dir = f"{base}/2022_Data/"
        datafile = f"{data_dir}/Data_eraCD_preEE.root"

    # ----------------------------
    # 2022 postEE
    # ----------------------------
    elif YEAR == "2022postEE":
        lumi = 26.67
        #mc_dir   = f"{base}/2022EE_MC/PROD_samplesNano_2022EE_MC_cc84ce40/"
        #data_dir = f"{base}/2022_Data/PROD_samplesNano_2022_Data_cc84ce40/"
        mc_dir   = f"{base}/2022EE_MC/"
        data_dir = f"{base}/2022_Data/"
        datafile = f"{data_dir}/Data_eraEFG_postEE.root"

    # ----------------------------
    # 2023 preBPix
    # ----------------------------
    elif YEAR == "2023preBPix":
        lumi = 18.06
        mc_dir   = f"{base}/2023preBPix_MC/"
        data_dir = f"{base}/2023_Data/"
        #mc_dir   = f"{base}/2023preBPix_MC/PROD_samplesNano_2023preBPix_MC_cc84ce40/"
        #data_dir = f"{base}/2023_Data/PROD_samplesNano_2023_Data_cc84ce40/"
        datafile = f"{data_dir}/Data_eraC_preBPix.root"

    # ----------------------------
    # 2023 postBPix
    # ----------------------------
    elif YEAR == "2023postBPix":
        lumi = 9.69
        mc_dir   = f"{base}/2023postBPix_MC/"
        data_dir = f"{base}/2023_Data/"
        #mc_dir   = f"{base}/2023postBPix_MC/PROD_samplesNano_2023postBPix_MC_cc84ce40/"
        #data_dir = f"{base}/2023_Data/PROD_samplesNano_2023_Data_cc84ce40"
        datafile = f"{data_dir}/Data_eraD_postBPix.root"

    # ----------------------------
    # Add more years here...
    # ----------------------------
    else:
        raise ValueError(f"Unknown YEAR={YEAR}")

    # Common MC samples
    fname_dy = f"{mc_dir}/DYJetsToLL/ZZ4lAnalysis.root"
    fname_tt = f"{mc_dir}/TTto2L2Nu/ZZ4lAnalysis.root"
    fname_wz = f"{mc_dir}/WZto3LNu/ZZ4lAnalysis.root"

    return CHANNEL, fname_dy, fname_tt, fname_wz, datafile, lumi


def select_channel(zflav, CHANNEL):
    zf = abs(zflav)
    
    if CHANNEL == "mu":          # muons
        sel = (zf == 169)

    elif CHANNEL == "e":         # electrons
        sel = (zf == 121)

    elif CHANNEL == "incl":      # muons + electrons
        sel = (zf == 169) | (zf == 121)

    else:
        raise ValueError("CHANNEL must be 'mu', 'e', or 'incl'")
    
    return ak.flatten(sel)

def poisson_interval(k, alpha=0.317): 
    a = alpha
    low, high = (chi2.ppf(a/2, 2*k) / 2, chi2.ppf(1-a/2, 2*k + 2) / 2)
    return k-low, high-k

def get_genEventSumw(input_file, maxEntriesPerSample=None):
    f = input_file
    runs = f.Runs
    event = f.Events
    nRuns = runs.GetEntries()
    nEntries = event.GetEntries()

    iRun = 0
    genEventCount = 0
    genEventSumw = 0.

    while iRun < nRuns and runs.GetEntry(iRun):
        genEventCount += runs.genEventCount
        genEventSumw += runs.genEventSumw
        iRun += 1
    
    print("gen=", genEventCount, "sumw=", genEventSumw)

    print("Total genEventCount=", genEventCount, "Total sumw=", genEventSumw)
    print(f"Entries in dataset: {nEntries}")

    if genEventCount == 0:
        print("WARNING: genEventCount is 0. Possible issue with input file or dataset.")

    if maxEntriesPerSample is not None:
        print(f"Scaling to {maxEntriesPerSample} entries")
        if nEntries > maxEntriesPerSample:
            genEventSumw = genEventSumw * maxEntriesPerSample / nEntries
            nEntries = maxEntriesPerSample
        print("    scaled to:", nEntries, "sumw=", genEventSumw)

    return genEventSumw


def plot_ZCand_samples(bestZIdx, ZCand_mass, ZCand_flav, sample_name):
    fig, ax = plt.subplots(3, 1, figsize=(10, 12))
    
    ax[0].hist(ak.to_numpy(bestZIdx), bins=50, color='blue', alpha=0.7)
    ax[0].set_xlabel('bestZIdx')
    ax[0].set_ylabel('Counts')
    ax[0].set_title(f'{sample_name} - bestZIdx')
    
    ax[1].hist(ak.flatten(ZCand_mass), bins=50, color='green', alpha=0.7)
    ax[1].set_xlabel('ZCand_mass')
    ax[1].set_ylabel('Counts')
    ax[1].set_title(f'{sample_name} - ZCand_mass')
    
    ax[2].hist(ak.flatten(ZCand_flav), bins=50, color='red', alpha=0.7)
    ax[2].set_xlabel('ZCand_flav')
    ax[2].set_ylabel('Counts')
    ax[2].set_title(f'{sample_name} - ZCand_flav')
    
    plt.tight_layout()
    plt.savefig(f'plot_{sample_name}.png')
    plt.close()

YEAR = '2022postEE'

CHANNEL, fname_dy, fname_tt, fname_wz, fname_dt, lumi = get_files_for_year(YEAR)

# Read files and get data
with uproot.open(fname_dy) as f:
    dy_zc = f['Events/bestZIdx'].array()
    dy_zm = f['Events/ZCand_mass'].array()
    dy_z1f = f['Events/ZCand_flav'].array()
    dy_ow = f['Events/overallEventWeight'].array()
    dy_cnt = get_genEventSumw(ROOT.TFile.Open(fname_dy))
    
with uproot.open(fname_tt) as f:
    tt_zc = f['Events/bestZIdx'].array()
    tt_zm = f['Events/ZCand_mass'].array()
    tt_z1f = f['Events/ZCand_flav'].array()
    tt_ow = f['Events/overallEventWeight'].array()
    tt_cnt = get_genEventSumw(ROOT.TFile.Open(fname_tt))
    
with uproot.open(fname_wz) as f:
    wz_zc = f['Events/bestZIdx'].array()
    wz_zm = f['Events/ZCand_mass'].array()
    wz_z1f = f['Events/ZCand_flav'].array()
    wz_ow = f['Events/overallEventWeight'].array()
    wz_cnt = get_genEventSumw(ROOT.TFile.Open(fname_wz))

with uproot.open(fname_dt) as f:
    dt_zc = f['Events/bestZIdx'].array()
    dt_zm = f['Events/ZCand_mass'].array()
    dt_z1f = f['Events/ZCand_flav'].array()

#plot_ZCand_samples(dy_zc, dy_zm, dy_z1f, "DYnew")
#plot_ZCand_samples(tt_zc, tt_zm, tt_z1f, "TTnew")
#plot_ZCand_samples(dt_zc, dt_zm, dt_z1f, "Datanew")

if dy_cnt == 0 or tt_cnt == 0 or wz_cnt == 0:
    raise ValueError("One or more genEventSumw values are zero, leading to division by zero.")

dt_z1mass = dt_zm[(ak.singletons(dt_zc))]
dt_z1flav = dt_z1f[(ak.singletons(dt_zc))]

m_dy = dy_zm[(ak.singletons(dy_zc))]
zf_dy = dy_z1f[(ak.singletons(dy_zc))]
w_dy = dy_ow

m_tt = tt_zm[(ak.singletons(tt_zc))]
zf_tt = tt_z1f[(ak.singletons(tt_zc))]
w_tt = tt_ow

m_wz = wz_zm[(ak.singletons(wz_zc))]
zf_wz = wz_z1f[(ak.singletons(wz_zc))]
w_wz = wz_ow

# Selections
#169 muons
#121 ele
#(abs(zf_dy) == 169) | (abs(zf_dy) == 121) inclusive

sel_dy = select_channel(zf_dy, CHANNEL)
m_dy = m_dy[sel_dy]
w_dy = w_dy[sel_dy]

sel_tt = select_channel(zf_tt, CHANNEL)
m_tt = m_tt[sel_tt]
w_tt = w_tt[sel_tt]

sel_wz = select_channel(zf_wz, CHANNEL)
m_wz = m_wz[sel_wz]
w_wz = w_wz[sel_wz]

sel_dt = select_channel(dt_z1flav, CHANNEL)
dt_z1mass = dt_z1mass[sel_dt]


m_tot = np.concatenate([ak.flatten(m_dy), ak.flatten(m_tt), ak.flatten(m_wz)])

w_tot_dy = w_dy * 1000 * lumi / dy_cnt
w_tot_tt = w_tt * 1000 * lumi / tt_cnt
w_tot_wz = w_wz * 1000 * lumi / wz_cnt

w_tot = np.concatenate([w_tot_dy, w_tot_tt, w_tot_wz])

CP_BINNING = np.linspace(60, 120, 100)


# The scaling of DY and TT xsecs are there because we spotted some wrong values in the csv used for the prod
# Doesn't really matter if then we normalize MC to Data. OTOH TT is off by one order of magnitude, so better to correct
n_dy, bins = np.histogram(ak.flatten(m_dy), weights=w_dy * 1000 * lumi / dy_cnt , bins=CP_BINNING) #* (6225.4 / 5558.0)
n_tt, bins = np.histogram(ak.flatten(m_tt), weights=w_tt * 1000 * lumi / tt_cnt , bins=CP_BINNING) #* (87.3 / 762.1)
n_wz, bins = np.histogram(ak.flatten(m_wz), weights=w_wz * 1000 * lumi / wz_cnt, bins=CP_BINNING)

n_dt, bins_dt = np.histogram(ak.flatten(dt_z1mass), bins=CP_BINNING)
mc_tot = n_dy + n_tt + n_wz

# Plot and normalize to Data
binsc = 0.5 * (bins[1:] + bins[:-1])
fig, ax = plt.subplots(
    nrows=2, 
    ncols=1,
    figsize=(14, 10),
    gridspec_kw={"height_ratios": (3, 1)},
    sharex=True
)

binsc_dt = 0.5 * (bins_dt[1:] + bins_dt[:-1])

ax[0].step(binsc, mc_tot * (sum(n_dt) / sum(mc_tot)), where='mid', label=r'MC (DY + t$\bar{t}$ + WZ)')
ax[0].errorbar(binsc_dt, n_dt, poisson_interval(n_dt), marker='.', color='k', linestyle='None', label='Data')
ax[0].set_ylabel('Events / bin width')

ax[1].errorbar(binsc, mc_tot * (sum(n_dt) / sum(mc_tot)) / n_dt, marker='.', color='k', linestyle='None')

# -----------------------------------------------------------
# Horizontal reference line at MC/Data = 1
# -----------------------------------------------------------
ax[1].axhline(1.0, color="black", linestyle="--", linewidth=1)


# ===========================================================
#   MC STATISTICAL UNCERTAINTY (per bin)
# ===========================================================
dy_w2, _ = np.histogram(ak.flatten(m_dy), weights=(w_tot_dy)**2, bins=CP_BINNING)
tt_w2, _ = np.histogram(ak.flatten(m_tt), weights=(w_tot_tt)**2, bins=CP_BINNING)
wz_w2, _ = np.histogram(ak.flatten(m_wz), weights=(w_tot_wz)**2, bins=CP_BINNING)

mc_stat_err = np.sqrt(dy_w2 + tt_w2 + wz_w2)
print("mc_stat_err", mc_stat_err)

# Normalize MC to data
scale_norm = sum(n_dt) / sum(mc_tot)
mc_scaled = mc_tot * scale_norm
mc_stat_scaled = mc_stat_err * scale_norm

ratio = mc_scaled / n_dt
ratio_stat_err = mc_stat_scaled / n_dt

# ===========================================================
# CONSTANT SYSTEMATIC UNCERTAINTIES (YOU PROVIDE THESE TWO)
# ===========================================================
# Example: scale = 1%, resolution = 4%
# Replace these two after you send the numbers

#2022 ele
#scale_value = 0.003
#res_value   = 0.12

#2023 ele
#scale_value = 0.002
#res_value   = 0.10

#2022 mu
scale_value = 0.0015
res_value   = 0.05

#2023 mu
#scale_value = 0.0015
#res_value   = 0.05


scale_syst = scale_value * np.ones_like(binsc)
res_syst   = res_value   * np.ones_like(binsc)

ratio_scale_err = ratio * scale_syst
ratio_res_err   = ratio * res_syst


# ===========================================================
#   TOTAL UNCERTAINTY (stat ⊕ scale ⊕ resolution)
# ===========================================================
ratio_tot_err = np.sqrt(
    ratio_stat_err**2 +
    ratio_scale_err**2 +
    ratio_res_err**2
)


# ===========================================================
#   DRAW TOTAL UNCERTAINTY BAND
# ===========================================================
ax[1].fill_between(
    binsc,
    1 - ratio_tot_err,
    1 + ratio_tot_err,
    color="blue",
    alpha=0.5,
    label="Total unc. (stat + scale + res)"
)

ax[1].legend()











#ax[1].set_xlabel(r'$m_{ee}$')
#ax[1].set_xlabel(r'$m_{\mu\mu}$')

if CHANNEL == "e":
    xlabel = r"$m_{ee}$"
elif CHANNEL == "mu":
    xlabel = r"$m_{\mu\mu}$"
elif CHANNEL == "incl":
    xlabel = r"$m_{\ell\ell}$"  # for inclusive, use generic lepton notation
else:
    xlabel = ""

ax[1].set_xlabel(xlabel)
ax[1].set_ylabel('MC/Data')
ax[1].set_ylim(0.5, 1.5)
# Add legend
ax[0].legend()

hep.cms.label(
    llabel="Preliminary",
    data=True,
    lumi=round(lumi, 2),
    com="13.6",  # <--- change this
    ax=ax[0]
)

plt.savefig(f"plot_{CHANNEL}_{YEAR}_Corr.png")
plt.savefig(f"plot_{CHANNEL}_{YEAR}_Corr.pdf")