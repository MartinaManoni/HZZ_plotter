import ROOT
import array 

# Enable multithreading (VERY important for speed)
ROOT.EnableImplicitMT()

# =========================
# Input
# =========================
fname = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_MC/DYJetsToLL/ZZ4lAnalysis_SKIMMED.root"

df = ROOT.RDataFrame("CRZLTree/candTree", fname)

# =========================
# CRZL SELECTION (exact translation of your C++)
# =========================
df_sel = df.Filter("Z1Mass > 40 && Z1Mass < 120") \
           .Filter("""
                !((LepPt[0] > LepPt[1] && (LepPt[0] < 20 || LepPt[1] < 10)) ||
                  (LepPt[1] > LepPt[0] && (LepPt[1] < 20 || LepPt[0] < 10)))
           """) \
           .Filter("abs(LepEta[2]) <= 2.5") \
           .Filter("!(LepSIP[2] > 4 || Lepdxy[2] > 0.5 || Lepdz[2] > 1.0)") \
           .Filter("MET <= 25") \
           .Filter("Nj == 1")

# =========================
# Define lepton-3 variables
# =========================
df_sel = df_sel.Define("pt", "LepPt[2]") \
               .Define("eta", "abs(LepEta[2])") \
               .Define("flav", "abs(LepLepId[2]) == 11 ? 0 : 1")  # 0=ele, 1=mu

# =========================
# Loose selection
# =========================
df_sel = df_sel.Define("isLoose", "pt > 5")

# =========================
# Tight selection (YOUR REAL LOGIC)
# =========================
df_sel = df_sel.Define(
    "isTight",
    """
    LepisID[2] &&
    (
        abs(LepLepId[2]) == 11
        ? true
        : (LepCombRelIsoPF[2] < 0.35)
    )
    """
)

# =========================
# Split numerator / denominator
# =========================
df_den_ele = df_sel.Filter("flav == 0 && isLoose")
df_num_ele = df_sel.Filter("flav == 0 && isLoose && isTight")

df_den_mu  = df_sel.Filter("flav == 1 && isLoose")
df_num_mu  = df_sel.Filter("flav == 1 && isLoose && isTight")

# =========================
# Binning
# =========================
pt_bins = array.array('d', [5, 7, 10, 20, 30, 40, 50, 80])

# =========================
# HISTOGRAMS (FIXED CORRECTLY)
# =========================
h_den_ele = df_den_ele.Histo1D(
    ("h_den_ele", "Electron Denominator; pT; Events", len(pt_bins)-1, pt_bins),
    "pt"
)

h_num_ele = df_num_ele.Histo1D(
    ("h_num_ele", "Electron Numerator; pT; Events", len(pt_bins)-1, pt_bins),
    "pt"
)

h_den_mu = df_den_mu.Histo1D(
    ("h_den_mu", "Muon Denominator; pT; Events", len(pt_bins)-1, pt_bins),
    "pt"
)

h_num_mu = df_num_mu.Histo1D(
    ("h_num_mu", "Muon Numerator; pT; Events", len(pt_bins)-1, pt_bins),
    "pt"
)

# =========================
# FORCE COMPUTATION
# =========================
den_e = h_den_ele.GetValue()
num_e = h_num_ele.GetValue()
den_m = h_den_mu.GetValue()
num_m = h_num_mu.GetValue()

# =========================
# COMPUTE FAKE RATES
# =========================
print("\n=== ELECTRON FAKE RATE ===")
for i in range(1, den_e.GetNbinsX()+1):
    fr = num_e.GetBinContent(i) / den_e.GetBinContent(i) if den_e.GetBinContent(i) > 0 else 0
    print(f"{den_e.GetBinLowEdge(i)}-{den_e.GetBinLowEdge(i+1)}: {fr:.4f}")

print("\n=== MUON FAKE RATE ===")
for i in range(1, den_m.GetNbinsX()+1):
    fr = num_m.GetBinContent(i) / den_m.GetBinContent(i) if den_m.GetBinContent(i) > 0 else 0
    print(f"{den_m.GetBinLowEdge(i)}-{den_m.GetBinLowEdge(i+1)}: {fr:.4f}")


import ROOT
import array

# =========================
# Get histograms
# =========================
den_e = h_den_ele.GetValue()
num_e = h_num_ele.GetValue()

den_m = h_den_mu.GetValue()
num_m = h_num_mu.GetValue()

# =========================
# Fake rate arrays
# =========================
x = array.array('d')
y_e = array.array('d')
y_m = array.array('d')

ex = array.array('d')
ey_e = array.array('d')
ey_m = array.array('d')

# =========================
# CANVAS for histograms
# =========================
c1 = ROOT.TCanvas("c1", "Num/Den distributions", 1200, 500)
c1.Divide(2, 1)

# =========================
# ELECTRONS: plot num + den
# =========================
c1.cd(1)
den_e.SetLineColor(ROOT.kBlue)
num_e.SetLineColor(ROOT.kRed)

den_e.Draw("HIST")
num_e.Draw("HIST SAME")

legend = ROOT.TLegend(0.6, 0.7, 0.88, 0.88)
legend.AddEntry(den_e, "Denominator", "l")
legend.AddEntry(num_e, "Numerator", "l")
legend.Draw()

# =========================
# MUONS: plot num + den
# =========================
c1.cd(2)
den_m.SetLineColor(ROOT.kBlue)
num_m.SetLineColor(ROOT.kRed)

den_m.Draw("HIST")
num_m.Draw("HIST SAME")

legend2 = ROOT.TLegend(0.6, 0.7, 0.88, 0.88)
legend2.AddEntry(den_m, "Denominator", "l")
legend2.AddEntry(num_m, "Numerator", "l")
legend2.Draw()

c1.SaveAs("num_den_distributions.png")

# =========================
# BUILD FAKE RATE GRAPH
# =========================
nbins = den_e.GetNbinsX()

for i in range(1, nbins + 1):

    # bin center
    xval = den_e.GetBinCenter(i)
    x.append(xval)
    ex.append(0)

    # electron FR
    d = den_e.GetBinContent(i)
    n = num_e.GetBinContent(i)

    fr = n / d if d > 0 else 0
    y_e.append(fr)

    # binomial error (safe)
    err = ROOT.TMath.Sqrt(fr * (1 - fr) / d) if d > 0 else 0
    ey_e.append(err)

    # muon FR
    d2 = den_m.GetBinContent(i)
    n2 = num_m.GetBinContent(i)

    fr2 = n2 / d2 if d2 > 0 else 0
    y_m.append(fr2)

    err2 = ROOT.TMath.Sqrt(fr2 * (1 - fr2) / d2) if d2 > 0 else 0
    ey_m.append(err2)

# =========================
# GRAPHS
# =========================
g_ele = ROOT.TGraphErrors(nbins, x, y_e, ex, ey_e)
g_mu  = ROOT.TGraphErrors(nbins, x, y_m, ex, ey_m)

g_ele.SetTitle("Electron Fake Rate; p_{T} [GeV]; Fake Rate")
g_mu.SetTitle("Muon Fake Rate; p_{T} [GeV]; Fake Rate")

g_ele.SetMarkerStyle(20)
g_mu.SetMarkerStyle(21)

g_ele.SetLineColor(ROOT.kRed)
g_mu.SetLineColor(ROOT.kBlue)

# =========================
# DRAW FAKE RATE PLOT
# =========================
c2 = ROOT.TCanvas("c2", "Fake Rates", 800, 600)

g_ele.Draw("AP")
g_mu.Draw("P SAME")

legend3 = ROOT.TLegend(0.6, 0.7, 0.88, 0.88)
legend3.SetFillStyle(0) 
legend3.SetBorderSize(0)
legend3.AddEntry(g_ele, "Electrons", "p")
legend3.AddEntry(g_mu, "Muons", "p")
legend3.Draw()

c2.SaveAs("fake_rates.png")