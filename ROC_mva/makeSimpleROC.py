import ROOT
import matplotlib.pyplot as plt

# Open ROOT file
DYNanoV12 = ROOT.TFile.Open(
    "/grid_mnt/data__data.polcms/cms/adewit/TauAnalysis/DY_2022_nanov12.root"
)
evts = DYNanoV12.Get("Events")

# Dummy histograms for counting
ele_prompt_all = ROOT.TH1F("ele_prompt_all", "ele_prompt_all", 10, 0, 5000)
ele_nonprompt_all = ROOT.TH1F("ele_nonprompt_all", "ele_nonprompt_all", 10, 0, 5000)
ele_prompt_pass = ROOT.TH1F("ele_prompt_pass", "ele_prompt_pass", 10, 0, 5000)
ele_nonprompt_pass = ROOT.TH1F("ele_nonprompt_pass", "ele_nonprompt_pass", 10, 0, 5000)

# --- Total prompt electrons ---
evts.Draw(
    "Electron_pt>>ele_prompt_all",
    "(Electron_genPartFlav==1 || Electron_genPartFlav==15) "
    "&& Electron_pt<10 && abs(Electron_eta)<0.8",
    "goff"
)
neles_allprompt = ele_prompt_all.Integral()

# --- Total nonprompt electrons ---
evts.Draw(
    "Electron_pt>>ele_nonprompt_all",
    "!(Electron_genPartFlav==1 || Electron_genPartFlav==15) "
    "&& Electron_pt<10 && abs(Electron_eta)<0.8",
    "goff"
)
neles_allnonprompt = ele_nonprompt_all.Integral()

# MVA cut scan values
cut_values = [
    -1, -0.9, -0.8, -0.7, -0.5, -0.4, -0.3, -0.2, -0.1, 0,
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
    0.82, 0.84, 0.86, 0.88, 0.9, 0.91, 0.92, 0.93,
    0.94, 0.95, 0.96, 0.97, 0.98, 0.99,
    0.991, 0.992, 0.993, 0.994, 0.995, 0.996, 0.997, 0.998, 0.999
]

signal_eff = []
background_eff = []

for cut in cut_values:
    # Reset histograms before reuse
    ele_prompt_pass.Reset()
    ele_nonprompt_pass.Reset()

    # Prompt passing cut
    evts.Draw(
        "Electron_pt>>ele_prompt_pass",
        f"(Electron_genPartFlav==1 || Electron_genPartFlav==15) "
        f"&& Electron_pt<10 && abs(Electron_eta)<0.8 "
        f"&& Electron_mvaHZZIso>{cut}",
        "goff"
    )
    neles_promptcut = ele_prompt_pass.Integral()

    # Nonprompt passing cut
    evts.Draw(
        "Electron_pt>>ele_nonprompt_pass",
        f"!(Electron_genPartFlav==1 || Electron_genPartFlav==15) "
        f"&& Electron_pt<10 && abs(Electron_eta)<0.8 "
        f"&& Electron_mvaHZZIso>{cut}",
        "goff"
    )
    neles_nonpromptcut = ele_nonprompt_pass.Integral()

    # Compute efficiencies (protect against division by zero)
    sig_eff = neles_promptcut / neles_allprompt if neles_allprompt > 0 else 0
    bkg_rej = 1 - (neles_nonpromptcut / neles_allnonprompt) if neles_allnonprompt > 0 else 0

    signal_eff.append(sig_eff)
    background_eff.append(bkg_rej)

# Plot ROC
fig, ax1 = plt.subplots(figsize=(8, 8), constrained_layout=True)

ax1.plot(signal_eff, background_eff, marker='o', linestyle='-', label="2018 training")

ax1.set_xlabel("Signal efficiency")
ax1.set_ylabel("Background rejection")
ax1.legend()

plt.savefig("roc_barrel_lowpt.pdf")
plt.close()