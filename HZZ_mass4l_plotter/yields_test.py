import ROOT

periods = ["2022", "2022EE", "2023preBPix", "2023postBPix", "2024"]

signals = [
    "ggH125",
    "VBFH125",
    "WplusH125",
    "WminusH125",
    "ZH125",
    "ttH125"
]

observables = [
    "costhetaZ1","costhetaZ2","costhetastar","phi","phi1",
    "pT4l","rapidity4l","Nj","dphijj","massZ1","massZ2",
    "pTj1","pTj2","mjj","absdetajj","pTHj","mHj","TCjMax","TBjMax"
]

def get_signal_yield(file, obs):
    total = 0
    for sig in signals:
        h = file.Get(f"{obs}_{sig}_105to160")
        if h:
            total += h.Integral("width")
    return total

def get_zx_yield(file, obs):
    h = file.Get(f"{obs}_ZX_105to160")
    return h.Integral("width") if h else 0

# Dictionary to store per-period yields per observable
yields = {obs: {} for obs in observables}

for period in periods:
    mc_file = ROOT.TFile.Open(f"H4l_MC_{period}_DIFF.root")
    zx_file = ROOT.TFile.Open(f"H4l_ZX_{period}_DIFF.root")

    for obs in observables:
        sig_y = get_signal_yield(mc_file, obs)
        zx_y  = get_zx_yield(zx_file, obs)

        yields[obs][period] = {"Signal": sig_y, "ZX": zx_y}

    mc_file.Close()
    zx_file.Close()

# Compute total Run-3 yields per observable
for obs in observables:
    total_sig = sum(yields[obs][p]["Signal"] for p in periods)
    total_zx  = sum(yields[obs][p]["ZX"] for p in periods)
    yields[obs]["Run3"] = {"Signal": total_sig, "ZX": total_zx}

# Print nicely
for obs in observables:
    print(f"\n================ YIELDS for {obs} (105 < m4l < 160) ================\n")
    print(f"{'Period':12s} {'Signal':>12s} {'ZX':>12s}")
    print("------------------------------------------------------")
    for period in periods + ["Run3"]:
        print(f"{period:12s} {yields[obs][period]['Signal']:12.3f} {yields[obs][period]['ZX']:12.3f}")
    print("======================================================\n")