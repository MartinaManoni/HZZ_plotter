### Draw and decorate plots produced with H4l_fill.py.
# This is just a quick example, for illustration purposes only!!!
# It lacks all the plots and features of the full miniAOD plotter.
#
# run 
# python3 H4l_draw_mZZ_periods2022.py

from __future__ import print_function
import glob
import optparse
import os
import os.path as osp
import sys
from datetime import date
import math
import ctypes
import ROOT
ROOT.gROOT.SetBatch(True)
import CMSGraphics, CMS_lumi
import numpy as np
from array import array
ROOT.PyConfig.IgnoreCommandLineOptions = True

inFilenameMC2024     = 'H4l_MC_2024_RDF.root'
inFilenameMC2024   = 'H4l_MC_2024_RDF.root'
inFilenameData2024   = "H4l_Data_2024_RDF.root"
inFilenameData2024 = "H4l_Data_2024_RDF.root"
inFilenameZX2024   = "H4l_ZX_2024_RDF.root"
inFilenameZX2024 = "H4l_ZX_2024_RDF.root"

lumi_24 = 108.82

outFilename = "Plots_2024_RDF.root"

## output directory
today = date.today()
print('Creating output dir...')
out_dir = str(today)+'_plots_2024'
os.makedirs(out_dir, exist_ok=True) #check if output dir exist


# plots options
blindPlots = True
blindHLow = 105.
blindHHi  = 140.
blindHM   = 500.
epsilon=0.1
addEmptyBins = True

# Set style matching the one used for HZZ plots
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetErrorX(0)
ROOT.gStyle.SetPadTopMargin(0.05)  
ROOT.gStyle.SetPadBottomMargin(0.13)
ROOT.gStyle.SetPadLeftMargin(0.16) 
ROOT.gStyle.SetPadRightMargin(0.03)
ROOT.gStyle.SetLabelOffset(0.008, "XYZ")
ROOT.gStyle.SetLabelSize(0.04, "XYZ")
ROOT.gStyle.SetAxisColor(1, "XYZ")
ROOT.gStyle.SetStripDecimals(True)
ROOT.gStyle.SetTickLength(0.03, "XYZ")
ROOT.gStyle.SetNdivisions(510, "XYZ")
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetTitleSize(0.05, "XYZ")
ROOT.gStyle.SetTitleOffset(1.00, "X")
ROOT.gStyle.SetTitleOffset(1.25, "Y")
ROOT.gStyle.SetLabelOffset(0.008, "XYZ")
ROOT.gStyle.SetLabelSize(0.04, "XYZ")

canvasSizeX=910
canvasSizeY=700

#####################
def printCanvases(type="png", path=".") :
    canvases = ROOT.gROOT.GetListOfCanvases()
    for c in canvases :
        c.Print(path+"/"+c.GetTitle()+"."+type)

def printCanvas(c, type="png", name=None, path="." ) :
    if name == None : name = c.GetTitle()
    name=name.replace(">","")
    name=name.replace("<","")
    name=name.replace(" ","_")
    c.Print(path+"/"+name+"."+type)


######################
def Stack(f2024, f2024ZX, variable="ZZmass_", version = "", finalState = 'fs_4l'):
    print("Stacking")
    # final state
    if(finalState == 'fs_4e'):
        fs_string = '4e_'
    elif(finalState == 'fs_4mu'):
        fs_string = '4mu_'
    elif(finalState == 'fs_2e2mu'):
        fs_string = '2e2mu_'
    elif(finalState == 'fs_4l'):
        fs_string = ''
    else:
        raise ValueError('Error: wrong final state!')

    # define histo name
    name = variable + (version if variable == "ZZMass_" else "") + fs_string
    print("variable", variable)
    print("version", version)
    print("fs_string", fs_string)
    print(f"[Stack] Using histogram base name: {name}")

    
    #------------EW------------------#
    #---------------------#
    # 2024
    WWZ24  = f2024.Get(name+"WWZ")
    WZZ24  = f2024.Get(name+"WZZ")
    ZZZ24  = f2024.Get(name+"ZZZ")
    EW24 = WWZ24.Clone("h_EW24")
    for i in [WZZ24, ZZZ24]:
        EW24.Add(i,1.)
    EW24.Scale(lumi_24*1000.)
    
    EW = EW24.Clone("h_EW")

    EW.SetLineColor(ROOT.TColor.GetColor("#000099"))
    EW.SetFillColor(ROOT.TColor.GetColor("#0331B9"))

    
    #-----------qqZZ---------------#
    # 2024
    ZZTo4l24 = f2024.Get(name+"ZZTo4l").Clone()
    ZZTo4l24.Scale(lumi_24*1000.)

    ZZTo4l = ZZTo4l24.Clone("h_ZZTo4l")

    ZZTo4l.SetLineColor(ROOT.TColor.GetColor("#000099"))
    ZZTo4l.SetFillColor(ROOT.TColor.GetColor("#99ccff"))
    
    #-----------signal------------#
    # 2024
    VBF125_24     = f2024.Get(name+"VBFH125")
    ggH125_24     = f2024.Get(name+"ggH125") 
    WplusH125_24  = f2024.Get(name+"WplusH125")
    WminusH125_24 = f2024.Get(name+"WminusH125")
    ZH125_24      = f2024.Get(name+"ZH125")
    ttH125_24     = f2024.Get(name+"ttH125")

    signal24 = VBF125_24.Clone("h_signal24")
    for i in [ggH125_24, WplusH125_24, WminusH125_24, ZH125_24, ttH125_24]:
        signal24.Add(i,1.)
    signal24.Scale(lumi_24*1000.)

    signal = signal24.Clone("h_signal")
      
    signal.SetLineColor(ROOT.TColor.GetColor("#cc0000"))
    signal.SetFillColor(ROOT.TColor.GetColor("#ff9b9b"))

    
    #------------ggTo-----------------#
    # 2024
    ggTo4mu_24     = f2024.Get(name+"ggTo4mu_Contin_MCFM701") 
    ggTo4e_24      = f2024.Get(name+"ggTo4e_Contin_MCFM701")
    ggTo4tau_24    = f2024.Get(name+"ggTo4tau_Contin_MCFM701")
    ggTo2e2mu_24   = f2024.Get(name+"ggTo2e2mu_Contin_MCFM701")
    ggTo2e2tau_24  = f2024.Get(name+"ggTo2e2tau_Contin_MCFM701")
    ggTo2mu2tau_24 = f2024.Get(name+"ggTo2mu2tau_Contin_MCFM701")

    ggToZZ24 = ggTo4mu_24.Clone("h_ggToZZ24")
    for i in [ggTo4e_24, ggTo4tau_24, ggTo2e2mu_24, ggTo2e2tau_24, ggTo2mu2tau_24]:
        ggToZZ24.Add(i,1.)
    ggToZZ24.Scale(lumi_24*1000.)
    ggToZZ = ggToZZ24.Clone("h_ggToZZ")

    # Set color/style
    ggToZZ.SetLineColor(ROOT.TColor.GetColor("#000099"))
    ggToZZ.SetFillColor(ROOT.TColor.GetColor("#4b78ff"))
    
    ##############################
    ##############################

    # 2024
    hzx_24 = f2024ZX.Get(name+"ZX").Clone()
    # hzx_24.Scale(lumi_24*1000.)
    hzx = hzx_24.Clone("h_ZX")

    # Set color/style
    hzx.SetLineColor(ROOT.TColor.GetColor("#003300"))
    hzx.SetFillColor(ROOT.TColor.GetColor("#669966"))
    
    #------------------Stack----------#
    if variable == "ZZmass_":
        if version=="4GeV_" :
            hs = ROOT.THStack("Stack_4GeV", "; m_{#it{4l}} (GeV) ; Events / 4 GeV" )
        elif version=="5GeV_" :
            hs = ROOT.THStack("Stack_5GeV", "; m_{#it{4l}} (GeV) ; Events / 5 GeV" )
        elif version=="10GeV_" :
            hs = ROOT.THStack("Stack_10GeV", "; m_{#it{4l}} (GeV) ; Events / 10 GeV" )
        else:
            hs = ROOT.THStack("Stack_2GeV", "; m_{#it{4l}} (GeV) ; Events / 2 GeV" )
    else:
        hs = ROOT.THStack(f"Stack_{variable}", f"; {variable} (GeV); Events")

    hs.Add(hzx,"HISTO")
    hs.Add(EW,"HISTO")
    hs.Add(ggToZZ,"HISTO")
    hs.Add(ZZTo4l,"HISTO")
    hs.Add(signal,"HISTO")
    
    return hs, [hzx, EW, ggToZZ, ZZTo4l, signal]


### Get a TGraph for data, blinded if required
def dataGraph (f1, variable="ZZMass_", version = "", finalState = 'fs_4l', blind = True):

    # final state
    if(finalState == 'fs_4e'):
        fs_string = '4e_'
    elif(finalState == 'fs_4mu'):
        fs_string = '4mu_'
    elif(finalState == 'fs_2e2mu'):
        fs_string = '2e2mu_'
    elif(finalState == 'fs_4l'):
        fs_string = ''
    else:
        raise ValueError('Error: wrong final state!')

    # define histo name
    name = variable + (version if variable == "ZZMass_" else "") + fs_string
    print(name)

    hd1 = f1.Get(name+"Data")
  
    # Combine all into one
    hd = hd1.Clone('h_data')
    
    nbinsIn = hd.GetNbinsX()
    nbins = 0
    
    x = np.array([0.]*nbinsIn, dtype='double')
    y = np.array([0.]*nbinsIn, dtype='double')
    errX = np.array([0.]*nbinsIn, dtype='double')  
    UpErr = np.array([0.]*nbinsIn, dtype='double')
    LowErr = np.array([0.]*nbinsIn, dtype='double')

    for i in range (1, nbinsIn):

        # New: blind 300–500 as well
        blindLow2 = 300.
        blindHigh2 = 400.

        if blind and (
                (blindHLow <= hd.GetBinCenter(i) <= blindHHi) or
                (blindLow2 <= hd.GetBinCenter(i) <= blindHigh2) or
                (hd.GetBinCenter(i) >= blindHM)
            ):
            continue

        if not addEmptyBins and hd.GetBinContent(i) == 0 : continue 
        x[nbins]      = hd.GetBinCenter(i)
        y[nbins]      = hd.GetBinContent(i)
        UpErr[nbins]  = hd.GetBinErrorUp(i)
        LowErr[nbins] = hd.GetBinErrorLow(i)
        nbins += 1

    Data = ROOT.TGraphAsymmErrors(nbins,x,y,errX,errX,LowErr,UpErr)
    Data.SetMarkerStyle(20)
    Data.SetLineColor(ROOT.kBlack)
    Data.SetMarkerSize(0.9)
    return Data



## --------------------------------
if __name__ == "__main__" :
    
    # 2024
    fMC2024       = ROOT.TFile.Open(inFilenameMC2024, "READ")
    fData2024     = ROOT.TFile.Open(inFilenameData2024, "READ")
    fZX2024       = ROOT.TFile.Open(inFilenameZX2024, "READ")

    of = ROOT.TFile.Open(outFilename, "recreate")

    # Labels for log plots
    xlabelsv = [80, 100, 150, 200, 250, 300, 350, 400, 450, 500]
    label_margin = -0.1
    xlabels = [None]*len(xlabelsv)
    for i, label in enumerate(xlabelsv): 
        xlabels[i] = ROOT.TLatex(label, label_margin , str(label))
        xlabels[i].SetTextAlign(23)
        xlabels[i].SetTextFont(42)
        xlabels[i].SetTextSize(0.04)

    # --- plots
    finalStates = ['fs_4e', 'fs_4mu', 'fs_2e2mu', 'fs_4l']
    variables = ['ZZMass_', 'Z1Mass_', 'Z2Mass_', 'Z1Mass_SR_', 'Z2Mass_SR_']

    # Lumi dictionary per year
    lumi_dict = {
        '2024': lumi_24  # adjust as needed
    }

    for fs in finalStates: 
        print(f"\nProcessing final state: {fs}")

        for variable in variables:
            versions = ["2GeV_", "4GeV_", "5GeV_", "10GeV_"] if variable == "ZZMass_" else [""]

            for version in versions:
                print(f"  Variable: {variable}, version: {version if version else 'no version'}")

                # Load histograms for all years
                HStack_var, h_list_var = Stack(
                    fMC2024, fZX2024,
                    variable=variable, version=version, finalState=fs
                )

                # Combine data from all years
                HData_var = dataGraph(fData2024,
                    variable=variable, version=version, finalState=fs,
                    blind=(blindPlots and variable=="ZZMass_")
                )

                # Create canvas
                canvas_name = f"{variable}{version}_full_allYears_{fs}".replace("__", "_")
                Canvas = ROOT.TCanvas(canvas_name, canvas_name, canvasSizeX, canvasSizeY)
                Canvas.SetTicks()

                # Set maximum
                xmin = ctypes.c_double(0.)
                ymin = ctypes.c_double(0.)
                xmax = ctypes.c_double(0.)
                ymax = ctypes.c_double(0.)
                HData_var.ComputeRange(xmin, ymin, xmax, ymax)
                HStack_var.SetMaximum(math.ceil(max(HStack_var.GetMaximum(), ymax.value)))

                # === Set axis labels on the first histogram in stack ===
                first_histo = h_list_var[0]  # Z+X is usually first
                if variable == "ZZMass_":
                    first_histo.GetXaxis().SetRangeUser(70., 300.)
                    first_histo.GetXaxis().SetTitle("m_{4l} (GeV)")

                    # Custom ticks every 50
                    for tick in range(100, 401, 50):
                        line = ROOT.TLine(tick, 0, tick, 0)  # Only labels, grid optional
                        # TLatex labels
                        lbl = ROOT.TLatex(tick, -0.05*HStack_var.GetMaximum(), str(tick))
                        lbl.SetTextAlign(23)
                        lbl.SetTextFont(42)
                        lbl.SetTextSize(0.04)
                        lbl.Draw()

                elif variable in ["Z1Mass", "Z2Mass", "Z1Mass_SR", "Z2Mass_SR"]:
                    first_histo.GetXaxis().SetRangeUser(0., 120.)
                    first_histo.GetXaxis().SetTitle(variable + " (GeV)")

                first_histo.GetXaxis().SetTitleSize(0.05)
                first_histo.GetXaxis().SetTitleOffset(1.1)
                first_histo.GetYaxis().SetTitle("Events")
                first_histo.GetYaxis().SetTitleSize(0.05)
                first_histo.GetYaxis().SetTitleOffset(1.25)

                # Draw stack
                HStack_var.Draw("HISTO")

                # Blinding for ZZMass
                if blindPlots and variable == "ZZMass_":
                    bblind = ROOT.TBox(blindHLow, 0, blindHHi, HStack_var.GetMaximum() - epsilon)
                    bblind.SetFillColor(ROOT.kGray)
                    bblind.SetFillStyle(3002)
                    bblind.Draw("same")

                    # New blinding box: 300–500 GeV
                    bblind2 = ROOT.TBox(300., 0, 400., HStack_var.GetMaximum() - epsilon)
                    bblind2.SetFillColor(ROOT.kGray)
                    bblind2.SetFillStyle(3002)
                    bblind2.Draw("same")

                # Draw data points on top
                HData_var.Draw("samePE1")

                # X-axis labels tweaks for ZZMass
                if variable == "ZZMass_":
                    HStack_var.GetXaxis().SetLabelSize(0)
                    for label in xlabels:
                        label.Draw()

                ROOT.gPad.RedrawAxis()

                # === Legend ===
                legend = ROOT.TLegend(0.72, 0.70, 0.94, 0.92)
                legend.AddEntry(HData_var, "Data", "p")
                legend.AddEntry(h_list_var[4], "H(125)", "f")
                legend.AddEntry(h_list_var[3], "q#bar{q}#rightarrow ZZ", "f")
                legend.AddEntry(h_list_var[2], "gg#rightarrow ZZ", "f")
                legend.AddEntry(h_list_var[1], "EW", "f")
                legend.AddEntry(h_list_var[0], "Z+X", "f")
                legend.SetFillColor(ROOT.kWhite)
                legend.SetLineColor(ROOT.kWhite)
                legend.SetTextFont(43)
                legend.SetTextSize(20)
                legend.Draw()

                # === CMS / Lumi ===
                CMS_lumi.writeExtraText = True
                CMS_lumi.extraText      = "Preliminary"
                CMS_lumi.lumi_sqrtS     = r"108.82 fb^{-1} (13.6 TeV)"
                CMS_lumi.cmsTextSize    = 1
                CMS_lumi.lumiTextSize   = 0.7
                CMS_lumi.extraOverCmsTextSize = 0.80
                CMS_lumi.relPosX = 0.12
                CMS_lumi.CMS_lumi(Canvas, 0, 0)

                # Save
                Canvas.Update()
                Canvas.Draw()
                Canvas.Write()
                printCanvas(Canvas, "png", path=out_dir)