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

    if genEventSumw == 0:
        print("WARNING: genEventSumw is 0. This will cause division errors in weight calculations.")


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

# File paths
#MATTEOS
#fname_dy = "/eos/user/m/mbonanom/nanoProd_Run3_2022EE_ControlPlots/samplesNanoDY_2022EE_MC_87a15506/DYJetsToLL/ZZ4lAnalysis.root"
#fname_tt = "/eos/user/m/mbonanom/nanoProd_Run3_2022EE_ControlPlots/samplesNanoTT_2022EE_MC_87a15506/TTto2L2Nu/ZZ4lAnalysis.root"
#fname_wz = "/eos/user/m/mbonanom/nanoProd_Run3_2022EE_ControlPlots/samplesNanoWZ_2022EE_MC_87a15506/WZto3LNu/ZZ4lAnalysis.root"
#fname_dt = "/eos/user/m/mbonanom/nanoProd_Run3_Data_ControlPlots/samplesNano_2022_Data_70a8ade8/post_EE/ZZ4lAnalysis.root"

#MARTINAS 2023 preBPix (ERA C)
#NO SS
'''fname_dy = "/eos/user/m/mmanoni/HZZ_samples23_newprod/MC_updated_Trig_noEleMuCorr/PROD_samplesNano_2023preBPix_MC_03abdca3/DYJetsToLL/ZZ4lAnalysis.root"
fname_tt = "/eos/user/m/mmanoni/HZZ_samples23_newprod/MC_updated_Trig_noEleMuCorr/PROD_samplesNano_2023preBPix_MC_03abdca3/TTto2L2Nu/ZZ4lAnalysis.root"
fname_wz = "/eos/user/m/mmanoni/HZZ_samples23_newprod/MC_updated_Trig_noEleMuCorr/PROD_samplesNano_2023preBPix_MC_03abdca3/WZto3LNu/ZZ4lAnalysis.root"
fname_dt = "/eos/user/m/mmanoni/HZZ_samples23_newprod/Data_updated_Trig_noEleMuCorr/PROD_samplesNano_2023_Data_03abdca3/Data_eraC_preBPix.root"'''

#WITH SS
'''fname_dy = "/eos/user/m/mmanoni/HZZ_samples23_newprod/MC_updated_Trig/PROD_samplesNano_2023preBPix_MC_03abdca3/DYJetsToLL/ZZ4lAnalysis.root"
fname_tt = "/eos/user/m/mmanoni/HZZ_samples23_newprod/MC_updated_Trig/PROD_samplesNano_2023preBPix_MC_03abdca3/TTto2L2Nu/ZZ4lAnalysis.root"
fname_wz = "/eos/user/m/mmanoni/HZZ_samples23_newprod/MC_updated_Trig/PROD_samplesNano_2023preBPix_MC_03abdca3/WZto3LNu/ZZ4lAnalysis.root"
fname_dt = "/eos/user/m/mmanoni/HZZ_samples23_newprod/Data_updated_Trig/PROD_samplesNano_2023_Data_03abdca3/Data_eraC_preBPix.root"'''


#Lourdes production preBPix
'''fname_dy = "/eos/user/l/lurda/CMS/HZZ/XS_analysis/250226/2023preBPix_MC/DYJetsToLL/ZZ4lAnalysis.root"
fname_tt = "/eos/user/l/lurda/CMS/HZZ/XS_analysis/250226/2023preBPix_MC/TTto2L2Nu/ZZ4lAnalysis.root"
fname_wz = "/eos/user/l/lurda/CMS/HZZ/XS_analysis/250226/2023preBPix_MC/WZto3LNu/ZZ4lAnalysis.root"
fname_dt = "/eos/user/l/lurda/CMS/HZZ/XS_analysis/250226/2023_Data/Data_eraC_preBPix.root"'''

#MARTINAS 2023 postBPix (ERA D) 
#no SS
'''fname_dy = "/eos/user/m/mmanoni/HZZ_samples23_newprod/MC_updated_Trig_noEleMuCorr/PROD_samplesNano_2023postBPix_MC_03abdca3/DYJetsToLL/ZZ4lAnalysis.root"
fname_tt = "/eos/user/m/mmanoni/HZZ_samples23_newprod/MC_updated_Trig_noEleMuCorr/PROD_samplesNano_2023postBPix_MC_03abdca3/TTto2L2Nu/ZZ4lAnalysis.root"
fname_wz = "/eos/user/m/mmanoni/HZZ_samples23_newprod/MC_updated_Trig_noEleMuCorr/PROD_samplesNano_2023postBPix_MC_03abdca3/WZto3LNu/ZZ4lAnalysis.root"
fname_dt = "/eos/user/m/mmanoni/HZZ_samples23_newprod/Data_updated_Trig_noEleMuCorr/PROD_samplesNano_2023_Data_03abdca3/Data_eraD_postBPix.root"'''
#with SS
'''fname_dy = "/eos/user/m/mmanoni/HZZ_samples23_newprod/MC_updated_Trig_newJetID_oldEleSS/PROD_samplesNano_2023postBPix_MC_03abdca3/DYJetsToLL/ZZ4lAnalysis.root"
fname_tt = "/eos/user/m/mmanoni/HZZ_samples23_newprod/MC_updated_Trig_newJetID_oldEleSS/PROD_samplesNano_2023postBPix_MC_03abdca3/TTto2L2Nu/ZZ4lAnalysis.root"
fname_wz = "/eos/user/m/mmanoni/HZZ_samples23_newprod/MC_updated_Trig_newJetID_oldEleSS/PROD_samplesNano_2023postBPix_MC_03abdca3/WZto3LNu/ZZ4lAnalysis.root"
fname_dt = "/eos/user/m/mmanoni/HZZ_samples23_newprod/Data_updated_Trig_newJetID_oldEleSS/PROD_samplesNano_2023_Data_03abdca3/Data_eraD_postBPix.root"'''

fname_dy = "/eos/user/m/mmanoni/HZZ_samples23_newprod/MC_newEleSS_2403/PROD_samplesNano_2023preBPix_MC/DYJetsToLL/ZZ4lAnalysis.root"
fname_tt = "/eos/user/m/mmanoni/HZZ_samples23_newprod/MC_newEleSS_2403/PROD_samplesNano_2023preBPix_MC/TTto2L2Nu/ZZ4lAnalysis.root"
fname_wz = "/eos/user/m/mmanoni/HZZ_samples23_newprod/MC_newEleSS_2403/PROD_samplesNano_2023preBPix_MC/WZto3LNu/ZZ4lAnalysis.root"
fname_dt = "/eos/user/m/mmanoni/HZZ_samples23_newprod/Data_newEleSS_2403/PROD_samplesNano_2023_Data/Data_eraC_preBPix.root"

# Determine luminosity
#17.794	Run3Summer23 era C pre BPix
#9.451	Run3Summer23BPix era D post BPix 
if "preBPix" in fname_dt:
    lumi = 17.794
else:
    lumi = 9.451

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

print("dy_zc", ak.type(dy_zc))
print("dy_zc", dy_zc)

plot_ZCand_samples(dy_zc, dy_zm, dy_z1f, "DYnew")
plot_ZCand_samples(tt_zc, tt_zm, tt_z1f, "TTnew")
plot_ZCand_samples(dt_zc, dt_zm, dt_z1f, "Datanew")

if dy_cnt == 0 or tt_cnt == 0 or wz_cnt == 0:
    raise ValueError("One or more genEventSumw values are zero, leading to division by zero.")

print("ak.type(dt_zc", ak.type(dt_zc))  # Check if it's a flat array or nested
print("ak.type(dt_zm", ak.type(dt_zm))

print("dt_zc", dt_zc)
print("dt_zm", dt_zm)
print("dt_z1f", dt_z1f)

print("ak.singletons(dt_zc))", ak.singletons(dt_zc))
dt_z1mass = dt_zm[(ak.singletons(dt_zc))]
print("dt_z1mass", dt_z1mass)
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

print(f"DY weights - Min: {np.min(w_dy)}, Max: {np.max(w_dy)}")
print(f"TT weights - Min: {np.min(w_tt)}, Max: {np.max(w_tt)}")
print(f"WZ weights - Min: {np.min(w_wz)}, Max: {np.max(w_wz)}")

print(f"dy_cnt: {dy_cnt}, tt_cnt: {tt_cnt}, wz_cnt: {wz_cnt}")

# Selections
#169 muons
#121 ele
#(abs(zf_dy) == 169) | (abs(zf_dy) == 121) inclusive

print("Ciao 0")
sel_dy = (abs(zf_dy) ==121)
sel_dy = ak.flatten(sel_dy)
m_dy = m_dy[sel_dy]
w_dy = w_dy[sel_dy]

sel_tt = (abs(zf_tt) == 121)
sel_tt = ak.flatten(sel_tt)
m_tt = m_tt[sel_tt]
w_tt = w_tt[sel_tt]

sel_wz = (abs(zf_wz) == 121)
sel_wz = ak.flatten(sel_wz)
m_wz = m_wz[sel_wz]
w_wz = w_wz[sel_wz]

sel_dt = (abs(dt_z1flav) == 121)
sel_dt = ak.flatten(sel_dt)
dt_z1mass = dt_z1mass[sel_dt]

print("dt_z1mass", dt_z1mass)
print("Ciao 1")

m_tot = np.concatenate([ak.flatten(m_dy), ak.flatten(m_tt), ak.flatten(m_wz)])

print("Checking for NaN or Inf in weights prima...")
print("NaN in w_dy:", np.any(np.isnan(w_dy)), "Inf in w_dy:", np.any(np.isinf(w_dy)))
print("NaN in w_tt:", np.any(np.isnan(w_tt)), "Inf in w_tt:", np.any(np.isinf(w_tt)))
print("NaN in w_wz:", np.any(np.isnan(w_wz)), "Inf in w_wz:", np.any(np.isinf(w_wz)))

w_tot_dy = w_dy * 1000 * lumi / dy_cnt
w_tot_tt = w_tt * 1000 * lumi / tt_cnt
w_tot_wz = w_wz * 1000 * lumi / wz_cnt

print("Ciao 2")
w_tot = np.concatenate([w_tot_dy, w_tot_tt, w_tot_wz])

CP_BINNING = np.linspace(60, 120, 100)

print("Checking for NaN or Inf in weights...")
print("NaN in w_dy:", np.any(np.isnan(w_dy)), "Inf in w_dy:", np.any(np.isinf(w_dy)))
print("NaN in w_tt:", np.any(np.isnan(w_tt)), "Inf in w_tt:", np.any(np.isinf(w_tt)))
print("NaN in w_wz:", np.any(np.isnan(w_wz)), "Inf in w_wz:", np.any(np.isinf(w_wz)))


# The scaling of DY and TT xsecs are there because we spotted some wrong values in the csv used for the prod
# Doesn't really matter if then we normalize MC to Data. OTOH TT is off by one order of magnitude, so better to correct
print("Ciao 3")
n_dy, bins = np.histogram(ak.flatten(m_dy), weights=w_dy * 1000 * lumi / dy_cnt , bins=CP_BINNING) #* (6225.4 / 5558.0)
n_tt, bins = np.histogram(ak.flatten(m_tt), weights=w_tt * 1000 * lumi / tt_cnt , bins=CP_BINNING) #* (87.3 / 762.1)
n_wz, bins = np.histogram(ak.flatten(m_wz), weights=w_wz * 1000 * lumi / wz_cnt, bins=CP_BINNING)

print("Ciao 4")
print("dt_z1mass", dt_z1mass)
print("ak.flatten(dt_z1mass)", ak.flatten(dt_z1mass))
n_dt, bins_dt = np.histogram(ak.flatten(dt_z1mass), bins=CP_BINNING)
print("n_dt", n_dt)
print("bins_dt", bins_dt)
mc_tot = n_dy + n_tt + n_wz

# Plot and normalize to Data
binsc = 0.5 * (bins[1:] + bins[:-1])
fig, ax = plt.subplots(
    nrows=2, 
    ncols=1,
    figsize=(10, 8),
    gridspec_kw={"height_ratios": (3, 1)},
    sharex=True
)

binsc_dt = 0.5 * (bins_dt[1:] + bins_dt[:-1])

print("Ciao 5")

ax[0].step(binsc, mc_tot * (sum(n_dt) / sum(mc_tot)), where='mid', label='MC')
ax[0].errorbar(binsc_dt, n_dt, poisson_interval(n_dt), marker='.', color='k', linestyle='None', label='Data')
ax[0].set_ylabel('Events / bin width')
print("Ciao 6")
print("binsc", binsc)
print("mc_tot", mc_tot)
print("sum(n_dt)",sum(n_dt))
print("sum(mc_tot)",sum(mc_tot))

ax[1].errorbar(binsc, mc_tot * (sum(n_dt) / sum(mc_tot)) / n_dt, marker='.', color='k', linestyle='None')
#ax[1].set_xlabel(r'$m_{ee}$')
#ax[1].set_xlabel(r'$m_{\mu\mu}$')
ax[1].set_xlabel(r'$m_{ee}$')
ax[1].set_ylabel('MC/Data')
ax[1].set_ylim(0.5, 2)
# Add legend
ax[0].legend()
print("Ciao 7")

hep.cms.label(
    llabel="Preliminary",
    data=True,
    lumi=round(lumi, 2),
    ax=ax[0]
)

plt.savefig("plot_ee_newSS_29Mar.png")