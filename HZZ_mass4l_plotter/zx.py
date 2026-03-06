import numpy as np
import pandas as pd
import uproot
from math import sqrt, log
import ROOT

# Configuration
years = ["2022", "2022EE", "2023preBPix", "2023postBPix", "2024"]
eos_path_template = '/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/FAKERATES/{}/'
branches_ZX = [
    "ZZMass",
    "Z1Flav",
    "Z2Flav",
    "LepLepId",
    "LepEta",
    "LepPt",
    "Z1Mass",
    "Z2Mass",
    "costheta1",
    "costheta2",
    "costhetastar",
    "Phi",
    "Phi1",
    "ZZPt",
    "ZZy",
    "Nj",
    "mjj",
    "absdetajj",
    "dphijj",
    "pTj1",
    "pTj2",
    "pTHj",
    "pTHjj",
    "mHj",
    "TBjMax",
    "TCjMax",
]


def FindFinalState(z1_flav, z2_flav):
    if z1_flav == -121:
        if z2_flav == 121:
            return 0  # 4e
        if z2_flav == 169:
            return 2  # 2e2mu
    if z1_flav == -169:
        if z2_flav == 121:
            return 3  # 2mu2e
        if z2_flav == 169:
            return 1  # 4mu


def GetFakeRate(lep_Pt, lep_eta, lep_ID, g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE):
    my_lep_Pt = min(lep_Pt, 79.)
    my_lep_ID = abs(lep_ID)

    if 5 < my_lep_Pt <= 7:
        bin = 0
    elif 7 < my_lep_Pt <= 10:
        bin = 1
    elif 10 < my_lep_Pt <= 20:
        bin = 2
    elif 20 < my_lep_Pt <= 30:
        bin = 3
    elif 30 < my_lep_Pt <= 40:
        bin = 4
    elif 40 < my_lep_Pt <= 50:
        bin = 5
    elif 50 < my_lep_Pt <= 80:
        bin = 6
    else:
        bin = 6  # fallback

    if my_lep_ID == 11:
        bin = bin - 1  # Electron FR has no bin for [5, 7]
        if abs(lep_eta) < 1.479:
            return g_FR_e_EB.GetY()[bin]
        else:
            return g_FR_e_EE.GetY()[bin]

    if my_lep_ID == 13:
        if abs(lep_eta) < 1.2:
            return g_FR_mu_EB.GetY()[bin]
        else:
            return g_FR_mu_EE.GetY()[bin]


def openFR(year):
    print("OpenFR")
    eos_path = eos_path_template.format(year)
    fnameFR = eos_path + f'FakeRates_SS_{year}.root'

    # Open file with ROOT to access TGraphErrors
    input_file_FR = ROOT.TFile(fnameFR)
    g_FR_mu_EB = input_file_FR.Get("FR_SS_muon_EB")
    g_FR_mu_EE = input_file_FR.Get("FR_SS_muon_EE")
    g_FR_e_EB = input_file_FR.Get("FR_SS_electron_EB")
    g_FR_e_EE = input_file_FR.Get("FR_SS_electron_EE")

    return g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE


def findFSZX(df):
    print("findFSZX")
    df['FinState'] = [FindFinalState(x, y) for x, y in zip(df['Z1Flav'], df['Z2Flav'])]
    return df


def comb(year):
    # cb_SS
    if year == "2022":
        return np.array([
            1.239, 
            1.093, 
            1.057, 
            1.254,
            ])  # 2022preEE OK
    if year == "2022EE":
        return np.array([
            1.067, # 4e
            1.015, # 4mu
            1.049, # 2e2mu
            0.905, # 2mu2e
            ])  # 2022EE OK
    if year == "2023preBPix":
        return np.array([
            1.116, 
            1.036, 
            0.989, 
            1.141,
            ])  # 2023preBPix OK
    if year == "2023postBPix":
        return np.array([
            0.795, # 4e
            1.025, # 4mu
            1.074, # 2e2mu
            1.078, # 2mu2e
            ])  # 2023postBPix OK
    if year == "2024":
        return np.array([
            0.787, # 4e
            0.960, # 4mu
            0.958, # 2e2mu
            0.749, # 2mu2e
            ])  # 2024 OK


def ratio(year):
    # fs_ROS_SS
    if year == "2022":
        return np.array([
            1.030, 
            1.165, 
            1.057,
            1.254,
            ])  # 2022preEE
    if year == "2022EE":
        return np.array([
            0.990,   # 4e
            0.997,  # 4mu
            1.039,   # 2e2mu
            1.016,  # 2mu2e
            ])  # 2022EE
    if year == "2023preBPix":
        return np.array([
            0.992, 
            1.024, 
            1.102, 
            1.024,
            ])  # 2023preBPix OK
    if year == "2023postBPix":
        return np.array([
            1.006,   # 4e
            1.040,  # 4mu
            1.078,   # 2e2mu
            1.025,  # 2mu2e
            ])  # 2023postBPix OK
    if year == "2024":
        return np.array([
            0.997,   # 4e
            1.028,  # 4mu
            1.051,   # 2e2mu
            1.024,  # 2mu2e
            ])  # 2024 OK

def ZXYield(df, year, g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE):
    print("ZXYield")
    cb_SS = comb(year)
    fs_ROS_SS = ratio(year)

    #print("df", df)

    #vec = df.to_numpy()
    #Yield = np.zeros(len(vec), float)

    #print(len(vec))
    #print("vec", vec)

    Yield = np.zeros(len(df), float)

    for i in range(len(df)):
        row = df.iloc[i]
        finSt = row['FinState']
        lepPt = row['LepPt']   # This is a list/array
        lepEta = row['LepEta'] # list/array
        lepID = row['LepLepId']  # list/array of lepton IDs

        #print("lepPt", lepPt)
        #print("lepPt2", lepPt[2])
        #print("lepEta", lepEta)
        #print("lepID", lepID)

        fr1 = GetFakeRate(lepPt[2], lepEta[2], lepID[2], g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE)
        fr2 = GetFakeRate(lepPt[3], lepEta[3], lepID[3], g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE)

        Yield[i] = cb_SS[finSt] * fs_ROS_SS[finSt] * fr1 * fr2

    return Yield

def doZX(year, g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE, data_files):
    dir_name = "CRZLLTree"  # TDirectoryFile
    tree_name = "candTree"   # Actual TTree inside the directory

    dfs = []  # List to hold DataFrames from each file

    for data in data_files:
        print(f"Reading data file: {data}")
        with uproot.open(data) as file:
            if dir_name not in file:
                raise KeyError(f"Directory '{dir_name}' not found in file {data}.")
            subdir = file[dir_name]
            if tree_name not in subdir:
                raise KeyError(f"TTree '{tree_name}' not found in '{dir_name}' in file {data}.")
            ttreeZX = subdir[tree_name]
            dfZX = ttreeZX.arrays(branches_ZX, library="pd")
            dfs.append(dfZX)

    # Concatenate all files into one DataFrame
    dfZX = pd.concat(dfs, ignore_index=True)

    # Keep only same-sign events
    dfZX = dfZX[dfZX.Z2Flav > 0]

    # Determine final state
    dfZX = findFSZX(dfZX)

    # Compute weights
    dfZX['weight1'] = ZXYield(dfZX, year, g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE)

    return dfZX

def save_to_root(df, output_path, tree_name="candTree"):
    # uproot requires arrays, convert DataFrame columns to dict of arrays
    arrays = {col: df[col].values for col in df.columns}
    
    with uproot.recreate(output_path) as root_file:
        root_file[tree_name] = arrays
    print(f"Saved ROOT file to {output_path}")

# ------------------------- Main -------------------------
def zx():
    d_ZX = {}

    for year in years:
        g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE = openFR(year)

        # Define data files per year
        if year != "2024":
            data_files_dict = {
                "2022": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2022_Data/Data_eraCD_preEE_SKIMMED.root",
                "2022EE": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2022_Data/Data_eraEFG_postEE_SKIMMED.root",
                "2023preBPix": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2023_Data/Data_eraC_preBPix_SKIMMED.root",
                "2023postBPix": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2023_Data/Data_eraD_postBPix_SKIMMED.root",
            }
            data_files = [data_files_dict[year]] 
        else:
            # Replace this with your list of 2024 files
            data_files = [
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

        d_ZX[year] = doZX(year, g_FR_mu_EB, g_FR_mu_EE, g_FR_e_EB, g_FR_e_EE, data_files)

        # Save the DataFrame with weights to ROOT file
        output_file = f"ZX_results_{year}_obs.root"
        save_to_root(d_ZX[year], output_file)

        print(year, 'done')

    return d_ZX


# If script is run directly
if __name__ == "__main__":
    zx()
