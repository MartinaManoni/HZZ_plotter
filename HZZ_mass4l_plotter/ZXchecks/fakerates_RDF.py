import ROOT
import array

#TO DOs
# Add 0 JETS
# Print out number of numerator and denominator in each category
# Run it for all backgrounds
# confirm again that the original code applies the same cuts

# Enable multithreading
ROOT.EnableImplicitMT()

# =========================
# Input
# =========================
fname = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_MC/WZto3LNu/ZZ4lAnalysis_SKIMMED.root"

df = ROOT.RDataFrame("CRZLTree/candTree", fname)

# =========================
# BASE SELECTION (NO Nj CUT HERE)
# =========================
df_sel = df.Filter("Z1Mass > 40 && Z1Mass < 120") \
           .Filter("""
                !((LepPt[0] > LepPt[1] && (LepPt[0] < 20 || LepPt[1] < 10)) ||
                  (LepPt[1] > LepPt[0] && (LepPt[1] < 20 || LepPt[0] < 10)))
           """) \
           .Filter("abs(LepEta[2]) <= 2.5") \
           .Filter("!(LepSIP[2] > 4 || Lepdxy[2] > 0.5 || Lepdz[2] > 1.0)") \
           .Filter("MET <= 25")

# =========================
# Define variables
# =========================
df_sel = df_sel.Define("pt", "LepPt[2]") \
               .Define("eta", "abs(LepEta[2])") \
               .Define("flav", "abs(LepLepId[2]) == 11 ? 0 : 1") \
               .Define("isLoose", "pt > 5") \
               .Define(
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
# Binning
# =========================
pt_bins = array.array('d', [5, 7, 10, 20, 30, 40, 50, 80])

# =========================
# FUNCTION TO BUILD FR GRAPHS
# =========================
def make_fake_rate_graphs(df_base, nj_cut, color_ele, color_mu):

    df_cut = df_base.Filter(nj_cut)

    # Split
    df_den_ele = df_cut.Filter("flav == 0 && isLoose")
    df_num_ele = df_cut.Filter("flav == 0 && isLoose && isTight")

    df_den_mu  = df_cut.Filter("flav == 1 && isLoose")
    df_num_mu  = df_cut.Filter("flav == 1 && isLoose && isTight")

    # Histos
    h_den_ele = df_den_ele.Histo1D(("hde_"+nj_cut, "", len(pt_bins)-1, pt_bins), "pt")
    h_num_ele = df_num_ele.Histo1D(("hne_"+nj_cut, "", len(pt_bins)-1, pt_bins), "pt")

    h_den_mu = df_den_mu.Histo1D(("hdm_"+nj_cut, "", len(pt_bins)-1, pt_bins), "pt")
    h_num_mu = df_num_mu.Histo1D(("hnm_"+nj_cut, "", len(pt_bins)-1, pt_bins), "pt")

    den_e = h_den_ele.GetValue()
    num_e = h_num_ele.GetValue()
    den_m = h_den_mu.GetValue()
    num_m = h_num_mu.GetValue()

    nbins = den_e.GetNbinsX()

    x = array.array('d')
    y_e = array.array('d')
    y_m = array.array('d')
    ex = array.array('d')
    ey_e = array.array('d')
    ey_m = array.array('d')

    print(f"\n=== Category: {nj_cut} ===")

    for i in range(1, nbins + 1):

        xval = den_e.GetBinCenter(i)
        x.append(xval)
        ex.append(0)

        # electrons
        d = den_e.GetBinContent(i)
        n = num_e.GetBinContent(i)
        fr = n/d if d > 0 else 0
        err_e = ROOT.TMath.Sqrt(fr*(1-fr)/d) if d > 0 else 0

        y_e.append(fr)
        ey_e.append(err_e)

        # muons
        d2 = den_m.GetBinContent(i)
        n2 = num_m.GetBinContent(i)
        fr2 = n2/d2 if d2 > 0 else 0
        err_m = ROOT.TMath.Sqrt(fr2*(1-fr2)/d2) if d2 > 0 else 0
        
        y_m.append(fr2)
        ey_m.append(err_m)

        # PRINT
        print(
            f"Bin {i} (pt ~ {xval:.1f}): "
            f"e: {n}/{d} = {fr:.4f} ± {err_e:.4f}, "
            f"mu: {n2}/{d2} = {fr2:.4f} ± {err_m:.4f}"
    )

    g_ele = ROOT.TGraphErrors(nbins, x, y_e, ex, ey_e)
    g_mu  = ROOT.TGraphErrors(nbins, x, y_m, ex, ey_m)

    # styles
    g_ele.SetMarkerStyle(20)
    g_mu.SetMarkerStyle(24)

    g_ele.SetMarkerColor(color_ele)
    g_ele.SetLineColor(color_ele)

    g_mu.SetMarkerColor(color_mu)
    g_mu.SetLineColor(color_mu)

    print("Electrons:")
    print("  Denominator:", den_e.Integral())
    print("  Numerator  :", num_e.Integral())

    print("Muons:")
    print("  Denominator:", den_m.Integral())
    print("  Numerator  :", num_m.Integral())

    return g_ele, g_mu

# =========================
# BUILD ALL CATEGORIES
# =========================
g_ele_0j, g_mu_0j = make_fake_rate_graphs(df_sel, "Nj == 0", ROOT.kGreen, ROOT.kGreen+2)
g_ele_1j, g_mu_1j = make_fake_rate_graphs(df_sel, "Nj == 1", ROOT.kRed, ROOT.kRed+2)
g_ele_2j, g_mu_2j = make_fake_rate_graphs(df_sel, "Nj >= 2", ROOT.kBlue, ROOT.kBlue+2)
g_ele_inc, g_mu_inc = make_fake_rate_graphs(df_sel, "Nj >= 0", ROOT.kBlack, ROOT.kGray+2)

# =========================
# PLOT
# =========================
c = ROOT.TCanvas("c", "Fake Rates vs Nj", 800, 600)

frame = c.DrawFrame(5, 0, 80, 1)
frame.SetTitle("Fake Rate; p_{T} [GeV]; Fake Rate")

# electrons
g_ele_0j.Draw("P SAME")
g_ele_1j.Draw("P SAME")
g_ele_2j.Draw("P SAME")
g_ele_inc.Draw("P SAME")

# muons
g_mu_0j.Draw("P SAME")
g_mu_1j.Draw("P SAME")
g_mu_2j.Draw("P SAME")
g_mu_inc.Draw("P SAME")

# legend
leg = ROOT.TLegend(0.55, 0.6, 0.88, 0.88)
leg.SetFillStyle(0)
leg.SetBorderSize(0)

leg.AddEntry(g_ele_0j, "e, N_{j} = 0", "p")
leg.AddEntry(g_ele_1j, "e, N_{j} = 1", "p")
leg.AddEntry(g_ele_2j, "e, N_{j} #geq 2", "p")
leg.AddEntry(g_ele_inc, "e, inclusive", "p")

leg.AddEntry(g_mu_0j, "#mu, N_{j} = 0", "p")
leg.AddEntry(g_mu_1j, "#mu, N_{j} = 1", "p")
leg.AddEntry(g_mu_2j, "#mu, N_{j} #geq 2", "p")
leg.AddEntry(g_mu_inc, "#mu, inclusive", "p")

leg.Draw()

c.SaveAs("fake_rates_vs_Nj_WZ_final.png")