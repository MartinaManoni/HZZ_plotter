#!/bin/env python3
# PLOTS FILLER by Martina 12/12/2025
import ROOT
import argparse
import os

ROOT.EnableImplicitMT()

# Import utility functions
#from ZZAnalysis.NanoAnalysis.tools import get_genEventSumw
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

# Define paths
MC_PATHS = {
    "2022": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2022_MC/",
    "2022EE": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2022EE_MC/",
    "2023preBPix": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2023preBPix_MC/",
    "2023postBPix": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2023postBPix_MC/",
    "2024": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_MC/",
}

DATA_PATHS = {
    "2022": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2022_Data/Data_eraCD_preEE.root",
    "2022EE": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2022_Data/Data_eraEFG_postEE.root",
    "2023preBPix": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2023_Data/Data_eraC_preBPix.root",
    "2023postBPix": "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2023_Data/Data_eraD_postBPix.root",
    "2024": [
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Cv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Dv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Ev1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Fv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Gv2/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Hv2/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Iv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma02024Iv2/ZZ4lAnalysis.root",

            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Cv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Dv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Ev1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Fv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Gv2/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Hv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Iv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/EGamma12024Iv2/ZZ4lAnalysis.root",

            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Cv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Dv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Ev1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Fv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Gv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Hv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Iv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon02024Iv2/ZZ4lAnalysis.root",

            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Cv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Dv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Ev1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Fv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Gv2/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Hv2/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Iv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/Muon12024Iv2/ZZ4lAnalysis.root",

            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Cv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Dv1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Ev1/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Fv2/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Gv3/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Hv2/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Iv2/ZZ4lAnalysis.root",
            "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/MuonEG2024Iv2v2/ZZ4lAnalysis.root"
            ]
}

ZX_PATHS = {
    "2022": "/afs/cern.ch/user/m/mmanoni/HZZ_plotter_CMSSW14/CMSSW_14_1_6/src/HZZ_plotter/HZZ_mass4l_plotter/ZX_results_2022_Njets_FIX.root",
    "2022EE": "/afs/cern.ch/user/m/mmanoni/HZZ_plotter_CMSSW14/CMSSW_14_1_6/src/HZZ_plotter/HZZ_mass4l_plotter/ZX_results_2022EE_Njets_FIX.root",
    "2023preBPix": "/afs/cern.ch/user/m/mmanoni/HZZ_plotter_CMSSW14/CMSSW_14_1_6/src/HZZ_plotter/HZZ_mass4l_plotter/ZX_results_2023preBPix_Njets_FIX.root",
    "2023postBPix": "/afs/cern.ch/user/m/mmanoni/HZZ_plotter_CMSSW14/CMSSW_14_1_6/src/HZZ_plotter/HZZ_mass4l_plotter/ZX_results_2023postBPix_Njets_FIX.root",
    "2024": "/afs/cern.ch/user/m/mmanoni/HZZ_plotter_CMSSW14/CMSSW_14_1_6/src/HZZ_plotter/HZZ_mass4l_plotter/ZX_results_2024_Njets_FIX.root",
}

Z_FLAVORS = {
    "4mu": (-169, -169),
    "4e": (-121, -121),
    "2e2mu": [(-169, -121), (-121, -169)],
}

LUMI = {
    "2022": 7.98, # fb^-1
    "2022EE": 26.67,
    "2023preBPix": 18.06,
    "2023postBPix": 9.69,
    "2024": 108.82,
}


# Add this near the top, after importing ROOT
ROOT.gInterpreter.Declare("""
#include <vector>
#include <cmath>

bool pass_dR(float eta1, float phi1, float eta2, float phi2) {
    float dphi = std::fabs(phi1 - phi2);
    if (dphi > M_PI) dphi = 2*M_PI - dphi;
    float deta = eta1 - eta2;
    return std::sqrt(deta*deta + dphi*dphi) > 0.4; // typical jet-lepton cleaning
}

ROOT::VecOps::RVec<int> get_jet_idx(const std::string& year,
                                    const ROOT::VecOps::RVec<unsigned char>& jetId,
                                    const ROOT::VecOps::RVec<bool>& zzMask,
                                    const ROOT::VecOps::RVec<float>& jet_pt,
                                    const ROOT::VecOps::RVec<float>& jet_eta,
                                    const ROOT::VecOps::RVec<float>& jet_phi,
                                    const ROOT::VecOps::RVec<bool>& lep_zzfullsel,
                                    const ROOT::VecOps::RVec<float>& lep_eta,
                                    const ROOT::VecOps::RVec<float>& lep_phi,
                                    int requiredId = 6,
                                    float minPt = 30.0,
                                    float maxAbsEta = 5.0)
{
    ROOT::VecOps::RVec<int> jet_indices_final;
    const std::size_t njet = jet_eta.size();
    const std::size_t nlep = lep_eta.size();

    for (std::size_t i = 0; i < njet; ++i) {
        if ((int)jetId[i] != requiredId) continue;
        if (std::abs(jet_eta[i]) >= maxAbsEta) continue;

        bool passdR = true;
        for (std::size_t j = 0; j < nlep; ++j) {
            if (!lep_zzfullsel[j]) continue;
            if (!pass_dR(jet_eta[i], jet_phi[i], lep_eta[j], lep_phi[j])) {
                passdR = false;
                break;
            }
        }
        if (!passdR) continue;

        float pt_cut = minPt;
        if ((std::abs(jet_eta[i]) >= 2.5f && std::abs(jet_eta[i]) < 3.0f) ||
            (std::abs(jet_eta[i]) >= 3.0f && (year == "2022" || year == "2022EE" || year == "2023" || year == "2023BPix")))
        {
            pt_cut = 50.0f;
        }

        if (jet_pt[i] > pt_cut) jet_indices_final.push_back(i);
    }

    return jet_indices_final;
}
""")

ROOT.gInterpreter.Declare("""
ROOT::RVec<float> concatenate(ROOT::RVec<float> &A, ROOT::RVec<float> &B){
    int sizeA = A.size();
    int sizeB = B.size();
    int sizeC = sizeA + sizeB;
    ROOT::RVec<float> C(sizeC);

    for (int i = 0; i < sizeA; ++i) {
        C[i] = A[i];
    }
    for (int i = 0; i < sizeB; ++i) {
        C[sizeA+i] = B[i];
    }

    return C;
}
""")

ROOT.gInterpreter.Declare("""
ROOT::RVec<short> concatenate(ROOT::RVec<int> &A, ROOT::RVec<int> &B){
    int sizeA = A.size();
    int sizeB = B.size();
    int sizeC = sizeA + sizeB;
    ROOT::RVec<float> C(sizeC);

    for (int i = 0; i < sizeA; ++i) {
        C[i] = A[i];
    }
    for (int i = 0; i < sizeB; ++i) {
        C[sizeA+i] = B[i];
    }

    return C;
}
""")

ROOT.gInterpreter.Declare("""
ROOT::RVec<bool> concatenate(ROOT::RVec<bool> &A, ROOT::RVec<bool> &B){
    int sizeA = A.size();
    int sizeB = B.size();
    int sizeC = sizeA + sizeB;
    ROOT::RVec<float> C(sizeC);

    for (int i = 0; i < sizeA; ++i) {
        C[i] = A[i];
    }
    for (int i = 0; i < sizeB; ++i) {
        C[sizeA+i] = B[i];
    }

    return C;
}
""")

ROOT.gInterpreter.Declare("""
ROOT::RVec<unsigned char> concatenate(ROOT::RVec<unsigned char> &A, ROOT::RVec<unsigned char> &B){
    int sizeA = A.size();
    int sizeB = B.size();
    int sizeC = sizeA + sizeB;
    ROOT::RVec<float> C(sizeC);

    for (int i = 0; i < sizeA; ++i) {
        C[i] = A[i];
    }
    for (int i = 0; i < sizeB; ++i) {
        C[sizeA+i] = B[i];
    }

    return C;
}
""")



def get_genEventSumw(input_file, maxEntriesPerSample=None):
    '''
       Util function to get the sum of weights per event.
       Returns the sum of weights, similarly to what we
       stored in Counters->GetBinContent(40) in the miniAODs.
    '''
    f = input_file

    runs  = f.Runs
    event = f.Events
    nRuns = runs.GetEntries()
    nEntries = event.GetEntries()

    iRun = 0
    genEventCount = 0
    genEventSumw = 0.

    while iRun < nRuns and runs.GetEntry(iRun) :
        genEventCount += runs.genEventCount
        genEventSumw += runs.genEventSumw
        iRun +=1
    print ("gen=", genEventCount, "sumw=", genEventSumw)

    if maxEntriesPerSample is not None:
        print(f"Scaling to {maxEntriesPerSample} entries")
        if nEntries>maxEntriesPerSample :
            genEventSumw = genEventSumw*maxEntriesPerSample/nEntries
            nEntries=maxEntriesPerSample
        print("    scaled to:", nEntries, "sumw=", genEventSumw)

    return genEventSumw

def define_histograms(df, df_SR, samplename, isMC):
    histos = {}
    weight_col = "weight" if isMC else None

    def book(df_input, hname, model, var):
        return df_input.Histo1D(model, var, weight_col) if weight_col else df_input.Histo1D(model, var)

    # General histograms (full range)

    histos[f"ZZMass_1GeV_{samplename}"] = book(df,
        f"ZZMass_1GeV_{samplename}",
        ROOT.RDF.TH1DModel("ZZMass_1GeV_" + samplename, "ZZMass_1GeV_" + samplename, 65, 105., 170.),
        "m4l")

    histos[f"ZZMass_2GeV_{samplename}"] = book(df,
        f"ZZMass_2GeV_{samplename}",
        ROOT.RDF.TH1DModel("ZZMass_2GeV_" + samplename, "ZZMass_2GeV_" + samplename, 65, 70., 200.),
        "m4l")

    '''histos[f"ZZMass_4GeV_{samplename}"] = book(df,
        f"ZZMass_4GeV_{samplename}",
        ROOT.RDF.TH1DModel("ZZMass_4GeV_" + samplename, "ZZMass_4GeV_" + samplename, 233, 70., 1002.),
        "m4l")'''

    histos[f"ZZMass_5GeV_{samplename}"] = book(df,
        f"ZZMass_5GeV_{samplename}",
        ROOT.RDF.TH1DModel("ZZMass_5GeV_" + samplename, "ZZMass_5GeV_" + samplename, 66, 70., 400.),
        "m4l")

    '''histos[f"ZZMass_10GeV_{samplename}"] = book(df,
        f"ZZMass_10GeV_{samplename}",
        ROOT.RDF.TH1DModel("ZZMass_10GeV_" + samplename, "ZZMass_10GeV_" + samplename, 93, 70., 1000.),
        "m4l")'''

    histos[f"Z1Mass_{samplename}"] = book(df,
        f"Z1Mass_{samplename}",
        ROOT.RDF.TH1DModel("Z1Mass_" + samplename, "Z1Mass_" + samplename, 40, 40., 120.),
        "Z1mass")

    histos[f"Z2Mass_{samplename}"] = book(df,
        f"Z2Mass_{samplename}",
        ROOT.RDF.TH1DModel("Z2Mass_" + samplename, "Z2Mass_" + samplename, 54, 12., 120.),
        "Z2mass")

    # Signal region only
    histos[f"Z1Mass_SR_{samplename}"] = book(df_SR,
        f"Z1Mass_SR_{samplename}",
        ROOT.RDF.TH1DModel("Z1Mass_SR_" + samplename, "Z1Mass_SR_" + samplename, 40, 40., 120.),
        "Z1mass")

    histos[f"Z2Mass_SR_{samplename}"] = book(df_SR,
        f"Z2Mass_SR_{samplename}",
        ROOT.RDF.TH1DModel("Z2Mass_SR_" + samplename, "Z2Mass_SR_" + samplename, 40, 0., 80.),
        "Z2mass")


    CHANNEL_TO_FINSTATE = {
        "4e": 0,
        "4mu": 1,
        "2e2mu": [2, 3],
    }

    for ch in ["4mu", "4e", "2e2mu"]:
        if "ZX" in samplename:
            # Use FinState instead of Z1/Z2 flav
            finstates = CHANNEL_TO_FINSTATE[ch]
            if isinstance(finstates, list):
                selection = " || ".join([f"(finState == {fs})" for fs in finstates])
            else:
                selection = f"(finState == {finstates})"
        else:
            flav = Z_FLAVORS[ch]
            if isinstance(flav, tuple):
                selection = f"(Z1flav == {flav[0]} && Z2flav == {flav[1]})"
            else:
                selection = " || ".join([f"(Z1flav == {f[0]} && Z2flav == {f[1]})" for f in flav])

        
        df_ch = df.Filter(selection)
        df_SR_ch = df_SR.Filter(selection)

        histos[f"ZZMass_1GeV_{ch}_{samplename}"] = book(df_ch,
            f"ZZMass_1GeV_{ch}_{samplename}",
            ROOT.RDF.TH1DModel(f"ZZMass_1GeV_{ch}_{samplename}", "", 65, 105., 170.),
            "m4l")

        histos[f"ZZMass_2GeV_{ch}_{samplename}"] = book(df_ch,
            f"ZZMass_2GeV_{ch}_{samplename}",
            ROOT.RDF.TH1DModel(f"ZZMass_2GeV_{ch}_{samplename}", "", 65, 70., 200.),
            "m4l")

        '''histos[f"ZZMass_4GeV_{ch}_{samplename}"] = book(df_ch,
            f"ZZMass_4GeV_{ch}_{samplename}",
            ROOT.RDF.TH1DModel(f"ZZMass_4GeV_{ch}_{samplename}", "", 233, 70., 1002.),
            "m4l")'''

        histos[f"ZZMass_5GeV_{ch}_{samplename}"] = book(df_ch,
            f"ZZMass_5GeV_{ch}_{samplename}",
            ROOT.RDF.TH1DModel(f"ZZMass_5GeV_{ch}_{samplename}", "", 66, 70., 400.),
            "m4l")
        
        '''histos[f"ZZMass_10GeV_{ch}_{samplename}"] = book(df_ch,
            f"ZZMass_10GeV_{ch}_{samplename}",
            ROOT.RDF.TH1DModel(f"ZZMass_10GeV_{ch}_{samplename}", "",93, 70., 1000.),
            "m4l")'''

        histos[f"Z1Mass_{ch}_{samplename}"] = book(df_ch,
            f"Z1Mass_{ch}_{samplename}",
            ROOT.RDF.TH1DModel(f"Z1Mass_{ch}_{samplename}", "", 40, 40., 120.),
            "Z1mass")

        histos[f"Z2Mass_{ch}_{samplename}"] = book(df_ch,
            f"Z2Mass_{ch}_{samplename}",
            ROOT.RDF.TH1DModel(f"Z2Mass_{ch}_{samplename}", "", 54, 12., 120.),
            "Z2mass")

        histos[f"Z1Mass_SR_{ch}_{samplename}"] = book(df_SR_ch,
            f"Z1Mass_SR_{ch}_{samplename}",
            ROOT.RDF.TH1DModel(f"Z1Mass_SR_{ch}_{samplename}", "", 40, 40., 120.),
            "Z1mass")

        histos[f"Z2Mass_SR_{ch}_{samplename}"] = book(df_SR_ch,
            f"Z2Mass_SR_{ch}_{samplename}",
            ROOT.RDF.TH1DModel(f"Z2Mass_SR_{ch}_{samplename}", "", 40, 0., 80.),
            "Z2mass")

    return histos

def run_sample(sample_name, filepaths, output_file, isMC=True):
    """
    Process a sample (MC or Data) using ROOT RDataFrame.
    filepaths can be a single string or a list of ROOT files.
    """
    # Ensure filepaths is a list
    if isinstance(filepaths, str):
        filepaths = [filepaths]

    # Check each file individually
    valid_files = [f for f in filepaths if os.path.exists(f)]
    if not valid_files:
        print(f"[WARNING] No valid files found for {sample_name}")
        return

    print(f"[INFO] Processing {sample_name} from {len(valid_files)} files")
    
    # Create RDataFrame from the list of valid files
    df = ROOT.RDataFrame("Events", valid_files)

    # Basic cuts
    df = df.Filter("bestCandIdx != -1").Filter("HLT_passZZ4l")

    if isMC:
        # Use the first file to get genEventSumw
        f = ROOT.TFile.Open(filepaths[0])
        genEventSumw = get_genEventSumw(f, 1e12)
        f.Close()

        df = df.Define("genEventSumw", str(genEventSumw))
        df = df.Define("weight", "overallEventWeight * ZZCand_dataMCWeight / genEventSumw")


        #genEventSumw = get_genEventSumw(ROOT.TFile.Open(filepaths), 1e12)
        #df = df.Define("genEventSumw", str(genEventSumw))
        #df = df.Define("weight", "overallEventWeight * ZZCand_dataMCWeight / genEventSumw")

    df = df.Define("m4l", "ZZCand_mass[bestCandIdx]") \
           .Define("Z1mass", "ZZCand_Z1mass[bestCandIdx]") \
           .Define("Z2mass", "ZZCand_Z2mass[bestCandIdx]") \
           .Define("Z1flav", "ZZCand_Z1flav[bestCandIdx]") \
           .Define("Z2flav", "ZZCand_Z2flav[bestCandIdx]")

    df_SR = df.Filter("m4l > 118 && m4l < 130")

    histos = define_histograms(df, df_SR, sample_name, isMC)
    # ----------------------------
    # Add Nj>=2 histograms
    # ----------------------------
    df = df.Define('electron_ZZFullSel','Electron_ZZFullSel') \
        .Define('muon_ZZFullSel','Muon_ZZFullSel') \
        .Define('Leptons_ZZFullSel', "concatenate(electron_ZZFullSel,muon_ZZFullSel)") \
        .Define('lep_eta', "concatenate(Electron_eta, Muon_eta)") \
        .Define('lep_phi', "concatenate(Electron_phi, Muon_phi)") \
        .Define("jet_indices", f'get_jet_idx("{sample_name}", Jet_jetId, Jet_ZZMask, Jet_pt, Jet_eta, Jet_phi, Leptons_ZZFullSel, lep_eta, lep_phi)') \
        .Define("nJets", "jet_indices.size()")

    df_Nj2 = df.Filter("nJets >= 2")
    df_SR_Nj2 = df_Nj2.Filter("m4l > 118 && m4l < 130")

    histos_Nj2 = define_histograms(df_Nj2, df_SR_Nj2, f"{sample_name}_Nj2", isMC)
    histos.update(histos_Nj2)
    print(f"[INFO] Added histograms for {sample_name} with Nj>=2")


    output_file.cd()
    for hname, h in histos.items():
        h_clone = ROOT.TH1F(hname, hname, h.GetValue().GetNbinsX(), h.GetValue().GetXaxis().GetXmin(), h.GetValue().GetXaxis().GetXmax())
        for i in range(1, h_clone.GetNbinsX() + 1):
            h_clone.SetBinContent(i, h.GetValue().GetBinContent(i))
            h_clone.SetBinError(i, h.GetValue().GetBinError(i))
        h_clone.Write()

def run_zx(period):

    path = ZX_PATHS.get(period)
    if not path:
        raise ValueError(f"Unknown ZX period: {period}")

    if not os.path.exists(path):
        print(f"[WARNING] ZX file not found: {path}")
        return

    df = ROOT.RDataFrame("candTree", path)

    # Define necessary branches
    df = df.Define("m4l", "ZZMass") \
           .Define("Z1mass", "Z1Mass") \
           .Define("Z2mass", "Z2Mass") \
           .Define("Z1flav", "Z1Flav") \
           .Define("Z2flav", "Z2Flav") \
           .Define("finState", "FinState") \
           .Define("weight", "weight1")

    df_SR = df.Filter("m4l > 118 && m4l < 130")

    fout = ROOT.TFile.Open(f"H4l_ZX_{period}_Nj2_FIX.root", "RECREATE")
    histos = define_histograms(df, df_SR, "ZX", isMC=True)

    print("df", df)

    df = df.Define("nJets", 'Nj')

    df_Nj2 = df.Filter("nJets >= 2")
    df_SR_Nj2 = df_Nj2.Filter("m4l > 118 && m4l < 130")

    histos_Nj2 = define_histograms(df_Nj2, df_SR_Nj2, "ZX_Nj2", isMC=True)
    histos.update(histos_Nj2)

    for hname, h in histos.items():
        h_clone = ROOT.TH1F(hname, hname, h.GetValue().GetNbinsX(), h.GetValue().GetXaxis().GetXmin(), h.GetValue().GetXaxis().GetXmax())
        for i in range(1, h_clone.GetNbinsX() + 1):
            h_clone.SetBinContent(i, h.GetValue().GetBinContent(i))
            h_clone.SetBinError(i, h.GetValue().GetBinError(i))
        h_clone.Write()
    fout.Close()

def run_data(period):
    """
    Process Data for a given period and fill histograms using RDataFrame.

    Handles both single ROOT files and lists of ROOT files.
    """
    paths = DATA_PATHS.get(period)
    if not paths:
        raise ValueError(f"Unknown data period: {period}")

    # Ensure we have a list of files
    if isinstance(paths, str):
        file_list = [paths]
    else:
        file_list = paths

    # Filter out non-existent files
    valid_files = [f for f in file_list if os.path.exists(f)]
    if not valid_files:
        print(f"[WARNING] No valid data files found for {period}")
        return

    print(f"[INFO] Processing Data for period {period} from {len(valid_files)} files")

    # Open output file
    fout = ROOT.TFile.Open(f"H4l_Data_{period}_Nj2_FIX.root", "RECREATE")

    # Call run_sample with the list of files
    run_sample("Data", valid_files, fout, isMC=False)

    fout.Close()
    print(f"[INFO] Finished processing Data for period {period}")

def run_mc(period):
    path = MC_PATHS.get(period)

    if not path:
        raise ValueError(f"Unknown MC period: {period}")

    samples = [
        "WWZ", "WZZ", "ZZZ",
        "ggTo4mu_Contin_MCFM701", "ggTo4e_Contin_MCFM701",
        "ggTo4tau_Contin_MCFM701", "ggTo2e2mu_Contin_MCFM701",
        "ggTo2e2tau_Contin_MCFM701", "ggTo2mu2tau_Contin_MCFM701",
        "ZZTo4l",
        "VBFH125", "ggH125", "WplusH125", "WminusH125",
        "ZH125", "ttH125"
    ]

    fout = ROOT.TFile.Open(f"H4l_MC_{period}_Nj2_FIX.root", "RECREATE")

    # ------------------------------------------------
    # 1) Collect MC files
    # ------------------------------------------------
    mc_files = {}
    for sample in samples:
        filename = os.path.join(path, sample, "ZZ4lAnalysis.root")
        if os.path.exists(filename):
            mc_files[sample] = filename
        else:
            print(f"[WARNING] File not found for sample {sample}")

    # ------------------------------------------------
    # 3) Run histogram production
    # ------------------------------------------------
    for sample, filename in mc_files.items():
        print(f"[INFO] Processing {sample} from 1 file")
        run_sample(sample, filename, fout, isMC=True)

    fout.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["data", "mc", "zx", "both"], default="all")
    parser.add_argument("--period", choices=["2022", "2022EE", "2023preBPix", "2023postBPix", "2024"], default="2022")
    args = parser.parse_args()

    if args.mode == "data":
        run_data(args.period)
    elif args.mode == "mc":
        run_mc(args.period)
    elif args.mode == "zx":
        run_zx(args.period)
    elif args.mode == "all":
        run_mc(args.period)
        run_data(args.period)
        run_zx(args.period)