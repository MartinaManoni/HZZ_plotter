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

inFilenameMC2022     = 'H4l_MC_2022_RDF.root'
inFilenameMC2022EE   = 'H4l_MC_2022EE_RDF.root'
inFilenameData2022   = "H4l_Data_2022_RDF.root"
inFilenameData2022EE = "H4l_Data_2022EE_RDF.root"
inFilenameZX2022   = "H4l_ZX_2022_RDF.root"
inFilenameZX2022EE = "H4l_ZX_2022EE_RDF.root"

lumi_22   = 7.98 # 1/fb
lumi_22EE  = 26.67 # 1/fb

inFilenameMC2023pre     = 'H4l_MC_2023preBPix_RDF.root'
inFilenameMC2023pre   = 'H4l_MC_2023preBPix_RDF.root'
inFilenameData2023pre   = "H4l_Data_2023preBPix_RDF.root"
inFilenameData2023pre = "H4l_Data_2023preBPix_RDF.root"
inFilenameZX2023pre   = "H4l_ZX_2023preBPix_RDF.root"
inFilenameZX2023pre = "H4l_ZX_2023preBPix_RDF.root"

lumi_23pre   = 18.06

inFilenameMC2023post     = 'H4l_MC_2023postBPix_RDF.root'
inFilenameMC2023post   = 'H4l_MC_2023postBPix_RDF.root'
inFilenameData2023post   = "H4l_Data_2023postBPix_RDF.root"
inFilenameData2023post = "H4l_Data_2023postBPix_RDF.root"
inFilenameZX2023post   = "H4l_ZX_2023postBPix_RDF.root"
inFilenameZX2023post = "H4l_ZX_2023postBPix_RDF.root"

lumi_23post  = 9.69

inFilenameMC2024     = 'H4l_MC_2024_RDF.root'
inFilenameMC2024   = 'H4l_MC_2024_RDF.root'
inFilenameData2024   = "H4l_Data_2024_RDF.root"
inFilenameData2024 = "H4l_Data_2024_RDF.root"
inFilenameZX2024   = "H4l_ZX_2024_RDF.root"
inFilenameZX2024 = "H4l_ZX_2024_RDF.root"

lumi_24 = 108.82

outFilename = "Plots_fullRun3_RDF.root"

## output directory
today = date.today()
print('Creating output dir...')
out_dir = str(today)+'_plots_fullRun3'
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

def get_m4l_label(finalState):
    if finalState == "fs_4e":
        return "m_{4e} (GeV)"
    elif finalState == "fs_4mu":
        return "m_{4\\mu} (GeV)"
    elif finalState == "fs_2e2mu":
        return "m_{2e2\\mu} (GeV)"
    elif finalState == "fs_4l":
        return "m_{4\\ell} (GeV)"
    else:
        return "m_{4l} (GeV)"

######################
def Stack(f2022, f2022EE, f2022ZX, f2022EEZX, f2023, f2023post, f2023ZX, f2023postZX, f2024, f2024ZX, variable="ZZmass_", version = "", finalState = 'fs_4l'):
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
    # 2022
    WWZ  = f2022.Get(name+"WWZ")
    print("WWZ",name)
    WZZ  = f2022.Get(name+"WZZ")
    ZZZ  = f2022.Get(name+"ZZZ")
    EWSamples = [WWZ, WZZ, ZZZ]
    EW = WWZ.Clone("h_EW")
    for i in EWSamples:
        EW.Add(i,1.)
    EW.Scale(lumi_22*1000.)

    #---------------------#
    # 2022EE
    WWZee  = f2022EE.Get(name+"WWZ")
    WZZee  = f2022EE.Get(name+"WZZ")
    ZZZee  = f2022EE.Get(name+"ZZZ")
    EWSamplesee = [WZZee, ZZZee]
    EWee = WWZee.Clone("h_EWee")
    for i in EWSamplesee:
        EWee.Add(i,1.)
    EWee.Scale(lumi_22EE*1000.)

    # 2023preBPix
    WWZ23pre  = f2023.Get(name+"WWZ")
    WZZ23pre  = f2023.Get(name+"WZZ")
    ZZZ23pre  = f2023.Get(name+"ZZZ")
    EW23pre = WWZ23pre.Clone("h_EW23pre")
    for i in [WZZ23pre, ZZZ23pre]:
        EW23pre.Add(i,1.)
    EW23pre.Scale(lumi_23pre*1000.)
    EW.Add(EW23pre,1.) # add 2023preBPix

    # 2023postBPix
    WWZ23post  = f2023post.Get(name+"WWZ")
    WZZ23post  = f2023post.Get(name+"WZZ")
    ZZZ23post  = f2023post.Get(name+"ZZZ")
    EW23post = WWZ23post.Clone("h_EW23post")
    for i in [WZZ23post, ZZZ23post]:
        EW23post.Add(i,1.)
    EW23post.Scale(lumi_23post*1000.)
    EW.Add(EW23post,1.) # add 2023postBPix

    # 2024
    WWZ24  = f2024.Get(name+"WWZ")
    WZZ24  = f2024.Get(name+"WZZ")
    ZZZ24  = f2024.Get(name+"ZZZ")
    EW24 = WWZ24.Clone("h_EW24")
    for i in [WZZ24, ZZZ24]:
        EW24.Add(i,1.)
    EW24.Scale(lumi_24*1000.)
    
    EW.Add(EW24,1.) # add 2024

    EW.SetLineColor(ROOT.TColor.GetColor("#000099"))
    EW.SetFillColor(ROOT.TColor.GetColor("#0331B9"))

    
    #-----------qqZZ---------------#
    # 2022
    ZZTo4l = f2022.Get(name+"ZZTo4l").Clone()
    ZZTo4l.Scale(lumi_22*1000.)

    # 2022EE
    ZZTo4lee = f2022EE.Get(name+"ZZTo4l").Clone()
    ZZTo4lee.Scale(lumi_22EE*1000.)
    ZZTo4l.Add(ZZTo4lee,1.) # full 2022

    # 2023preBPix
    ZZTo4l23pre = f2023.Get(name+"ZZTo4l").Clone()
    ZZTo4l23pre.Scale(lumi_23pre*1000.)
    ZZTo4l.Add(ZZTo4l23pre,1.) # add 2023preBPix

    # 2023postBPix
    ZZTo4l23post = f2023post.Get(name+"ZZTo4l").Clone()
    ZZTo4l23post.Scale(lumi_23post*1000.)
    ZZTo4l.Add(ZZTo4l23post,1.) # add 2023postBPix

    # 2024
    ZZTo4l24 = f2024.Get(name+"ZZTo4l").Clone()
    ZZTo4l24.Scale(lumi_24*1000.)
    ZZTo4l.Add(ZZTo4l24,1.) # add 2024

    ZZTo4l.SetLineColor(ROOT.TColor.GetColor("#000099"))
    ZZTo4l.SetFillColor(ROOT.TColor.GetColor("#99ccff"))
    
    #-----------signal------------#
    # 2022
    VBF125     = f2022.Get(name+"VBFH125")
    ggH125     = f2022.Get(name+"ggH125") 
    WplusH125  = f2022.Get(name+"WplusH125")
    WminusH125 = f2022.Get(name+"WminusH125")
    ZH125      = f2022.Get(name+"ZH125")
    ttH125     = f2022.Get(name+"ttH125")

    signalSamples = [ggH125, WplusH125, WminusH125, ZH125, ttH125]
    signal = VBF125.Clone("h_signal")
    for i in signalSamples:
        signal.Add(i, 1.)
    signal.Scale(lumi_22*1000.) 

    # 2022EE
    VBF125ee     = f2022EE.Get(name+"VBFH125")
    ggH125ee     = f2022EE.Get(name+"ggH125")
    WplusH125ee  = f2022EE.Get(name+"WplusH125")
    WminusH125ee = f2022EE.Get(name+"WminusH125")
    ZH125ee      = f2022EE.Get(name+"ZH125")
    ttH125ee     = f2022EE.Get(name+"ttH125")

    signalSamplesee = [ggH125ee, WplusH125ee, WminusH125ee, ZH125ee, ttH125ee]
    signalee = VBF125ee.Clone("h_signalee")
    for i in signalSamplesee:
        signalee.Add(i, 1.)
    signalee.Scale(lumi_22EE*1000.)
    signal.Add(signalee,1.) # full 2022

    # 2023preBPix
    VBF125_23pre     = f2023.Get(name+"VBFH125")
    ggH125_23pre     = f2023.Get(name+"ggH125") 
    WplusH125_23pre  = f2023.Get(name+"WplusH125")
    WminusH125_23pre = f2023.Get(name+"WminusH125")
    ZH125_23pre      = f2023.Get(name+"ZH125")
    ttH125_23pre     = f2023.Get(name+"ttH125")

    signal23pre = VBF125_23pre.Clone("h_signal23pre")
    for i in [ggH125_23pre, WplusH125_23pre, WminusH125_23pre, ZH125_23pre, ttH125_23pre]:
        signal23pre.Add(i,1.)
    signal23pre.Scale(lumi_23pre*1000.)
    signal.Add(signal23pre,1.) # add 2023preBPix

    # 2023postBPix
    VBF125_23post     = f2023post.Get(name+"VBFH125")
    ggH125_23post     = f2023post.Get(name+"ggH125") 
    WplusH125_23post  = f2023post.Get(name+"WplusH125")
    WminusH125_23post = f2023post.Get(name+"WminusH125")
    ZH125_23post      = f2023post.Get(name+"ZH125")
    ttH125_23post     = f2023post.Get(name+"ttH125")

    signal23post = VBF125_23post.Clone("h_signal23post")
    for i in [ggH125_23post, WplusH125_23post, WminusH125_23post, ZH125_23post, ttH125_23post]:
        signal23post.Add(i,1.)
    signal23post.Scale(lumi_23post*1000.)
    signal.Add(signal23post,1.) # add 2023postBPix

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
    signal.Add(signal24,1.) # add 2024
      
    signal.SetLineColor(ROOT.TColor.GetColor("#cc0000"))
    signal.SetFillColor(ROOT.TColor.GetColor("#ff9b9b"))

    
    #------------ggTo-----------------#
    # 2022
    ggTo4mu     = f2022.Get(name+"ggTo4mu_Contin_MCFM701") 
    ggTo4e      = f2022.Get(name+"ggTo4e_Contin_MCFM701")
    ggTo4tau    = f2022.Get(name+"ggTo4tau_Contin_MCFM701")
    ggTo2e2mu   = f2022.Get(name+"ggTo2e2mu_Contin_MCFM701")
    ggTo2e2tau  = f2022.Get(name+"ggTo2e2tau_Contin_MCFM701")
    ggTo2mu2tau = f2022.Get(name+"ggTo2mu2tau_Contin_MCFM701")

    ggZZSamples = [ ggTo4e, ggTo4tau, ggTo2e2mu, ggTo2e2tau, ggTo2mu2tau]
    ggToZZ = ggTo4mu.Clone("h_ggTo")
    for i in ggZZSamples:
        ggToZZ.Add(i,1.)
    ggToZZ.Scale(lumi_22*1000.)

    # 2022EE
    ggTo4muee     = f2022EE.Get(name+"ggTo4mu_Contin_MCFM701") 
    ggTo4eee     = f2022EE.Get(name+"ggTo4e_Contin_MCFM701")
    ggTo4tauee    = f2022EE.Get(name+"ggTo4tau_Contin_MCFM701")
    ggTo2e2muee  = f2022EE.Get(name+"ggTo2e2mu_Contin_MCFM701")
    ggTo2e2tauee  = f2022EE.Get(name+"ggTo2e2tau_Contin_MCFM701")
    ggTo2mu2tauee = f2022EE.Get(name+"ggTo2mu2tau_Contin_MCFM701")

    ggZZSamplesee= [ ggTo4eee, ggTo4tauee, ggTo2e2muee, ggTo2e2tauee, ggTo2mu2tauee]
    ggToZZee = ggTo4muee.Clone("h_ggToee")
    for i in ggZZSamplesee:
        ggToZZee.Add(i,1.)
    ggToZZee.Scale(lumi_22EE*1000.)
    ggToZZ.Add(ggToZZee, 1.)

    # 2023preBPix
    ggTo4mu_23pre     = f2023.Get(name+"ggTo4mu_Contin_MCFM701") 
    ggTo4e_23pre      = f2023.Get(name+"ggTo4e_Contin_MCFM701")
    ggTo4tau_23pre    = f2023.Get(name+"ggTo4tau_Contin_MCFM701")
    ggTo2e2mu_23pre   = f2023.Get(name+"ggTo2e2mu_Contin_MCFM701")
    ggTo2e2tau_23pre  = f2023.Get(name+"ggTo2e2tau_Contin_MCFM701")
    ggTo2mu2tau_23pre = f2023.Get(name+"ggTo2mu2tau_Contin_MCFM701")

    ggToZZ23pre = ggTo4mu_23pre.Clone("h_ggTo23pre")
    for i in [ggTo4e_23pre, ggTo4tau_23pre, ggTo2e2mu_23pre, ggTo2e2tau_23pre, ggTo2mu2tau_23pre]:
        ggToZZ23pre.Add(i,1.)
    ggToZZ23pre.Scale(lumi_23pre*1000.)
    ggToZZ.Add(ggToZZ23pre,1.) # add 2023preBPix

    # 2023postBPix
    ggTo4mu_23post     = f2023post.Get(name+"ggTo4mu_Contin_MCFM701") 
    ggTo4e_23post      = f2023post.Get(name+"ggTo4e_Contin_MCFM701")
    ggTo4tau_23post    = f2023post.Get(name+"ggTo4tau_Contin_MCFM701")
    ggTo2e2mu_23post   = f2023post.Get(name+"ggTo2e2mu_Contin_MCFM701")
    ggTo2e2tau_23post  = f2023post.Get(name+"ggTo2e2tau_Contin_MCFM701")
    ggTo2mu2tau_23post = f2023post.Get(name+"ggTo2mu2tau_Contin_MCFM701")

    ggToZZ23post = ggTo4mu_23post.Clone("h_ggTo23post")
    for i in [ggTo4e_23post, ggTo4tau_23post, ggTo2e2mu_23post, ggTo2e2tau_23post, ggTo2mu2tau_23post]:
        ggToZZ23post.Add(i,1.)
    ggToZZ23post.Scale(lumi_23post*1000.)
    ggToZZ.Add(ggToZZ23post,1.) # add 2023postBPix

    # 2024
    ggTo4mu_24     = f2024.Get(name+"ggTo4mu_Contin_MCFM701") 
    ggTo4e_24      = f2024.Get(name+"ggTo4e_Contin_MCFM701")
    ggTo4tau_24    = f2024.Get(name+"ggTo4tau_Contin_MCFM701")
    ggTo2e2mu_24   = f2024.Get(name+"ggTo2e2mu_Contin_MCFM701")
    ggTo2e2tau_24  = f2024.Get(name+"ggTo2e2tau_Contin_MCFM701")
    ggTo2mu2tau_24 = f2024.Get(name+"ggTo2mu2tau_Contin_MCFM701")

    ggToZZ24 = ggTo4mu_24.Clone("h_ggTo24")
    for i in [ggTo4e_24, ggTo4tau_24, ggTo2e2mu_24, ggTo2e2tau_24, ggTo2mu2tau_24]:
        ggToZZ24.Add(i,1.)
    ggToZZ24.Scale(lumi_24*1000.)
    ggToZZ.Add(ggToZZ24,1.) # add 2024

    # Set color/style
    ggToZZ.SetLineColor(ROOT.TColor.GetColor("#000099"))
    ggToZZ.SetFillColor(ROOT.TColor.GetColor("#4b78ff"))
    
    ##############################
    ##############################
    # ZX
    # 2022
    hzx = f2022ZX.Get(name+"ZX").Clone()
    # hzx.Scale(lumi_22*1000.)

    hzxee = f2022EEZX.Get(name+"ZX").Clone()
    # hzxee.Scale(lumi_22EE*1000.)
    hzx.Add(hzxee,1.) # full 2022

    # 2023preBPix
    hzx_23pre = f2023ZX.Get(name+"ZX").Clone()
    # hzx_23pre.Scale(lumi_23pre*1000.)
    hzx.Add(hzx_23pre,1.) # add 2023preBPix

    # 2023postBPix
    hzx_23post = f2023postZX.Get(name+"ZX").Clone()
    # hzx_23post.Scale(lumi_23post*1000.)
    hzx.Add(hzx_23post,1.) # add 2023postBPix

    # 2024
    hzx_24 = f2024ZX.Get(name+"ZX").Clone()
    # hzx_24.Scale(lumi_24*1000.)
    hzx.Add(hzx_24,1.) # add 2024

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
def dataGraph (f1, f2, f3, f4, f5, variable="ZZMass_", version = "", finalState = 'fs_4l', blind = True):

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
    hd2 = f2.Get(name+"Data")
    hd3 = f3.Get(name+"Data")
    hd4 = f4.Get(name+"Data")
    hd5 = f5.Get(name+"Data")
  
    # Combine all into one
    hd = hd1.Clone('h_data')
    hd.Add(hd2, 1.)
    hd.Add(hd3, 1.)
    hd.Add(hd4, 1.)
    hd.Add(hd5, 1.)
    
    nbinsIn = hd.GetNbinsX()
    nbins = 0
    
    x = np.array([0.]*nbinsIn, dtype='double')
    y = np.array([0.]*nbinsIn, dtype='double')
    errX = np.array([0.]*nbinsIn, dtype='double')  
    UpErr = np.array([0.]*nbinsIn, dtype='double')
    LowErr = np.array([0.]*nbinsIn, dtype='double')

    for i in range (1, nbinsIn):
        #if blind and ((hd.GetBinCenter(i)>=blindHLow and hd.GetBinCenter(i)<=blindHHi) or hd.GetBinCenter(i)>=blindHM) : continue

        # New: blind 300–500 as well
        blindLow2 = 300.
        blindHigh2 = 500.

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
    

    fMC2022     = ROOT.TFile.Open(inFilenameMC2022,"READ")
    fData2022   = ROOT.TFile.Open(inFilenameData2022,"READ")
    fZX2022   = ROOT.TFile.Open(inFilenameZX2022,"READ")

    fMC2022EE   = ROOT.TFile.Open(inFilenameMC2022EE,"READ")
    fData2022EE = ROOT.TFile.Open(inFilenameData2022EE,"READ")
    fZX2022EE = ROOT.TFile.Open(inFilenameZX2022EE,"READ")

    # 2023preBPix
    fMC2023pre    = ROOT.TFile.Open(inFilenameMC2023pre, "READ")
    fData2023pre  = ROOT.TFile.Open(inFilenameData2023pre, "READ")
    fZX2023pre    = ROOT.TFile.Open(inFilenameZX2023pre, "READ")

    # 2023postBPix
    fMC2023post   = ROOT.TFile.Open(inFilenameMC2023post, "READ")
    fData2023post = ROOT.TFile.Open(inFilenameData2023post, "READ")
    fZX2023post   = ROOT.TFile.Open(inFilenameZX2023post, "READ")

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
        '2022': lumi_22,
        '2022EE': lumi_22EE,
        '2023pre': lumi_23pre,
        '2023post': lumi_23post,
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
                    fMC2022, fMC2022EE, fZX2022, fZX2022EE,
                    fMC2023pre, fMC2023post, fZX2023pre, fZX2023post,
                    fMC2024, fZX2024,
                    variable=variable, version=version, finalState=fs
                )

                # Combine data from all years
                HData_var = dataGraph(
                    fData2022, fData2022EE,
                    fData2023pre, fData2023post, fData2024,
                    variable=variable, version=version, finalState=fs,
                    blind=(blindPlots and variable=="ZZMass_")
                )

                # Create canvas
                canvas_name = f"{variable}{version}_full_allYears_{fs}".replace("__", "_")
                Canvas = ROOT.TCanvas(canvas_name, canvas_name, canvasSizeX, canvasSizeY)
                Canvas.SetTicks()

                print("CIAO")
                Canvas.GetListOfPrimitives().Print()

                # Set maximum
                xmin = ctypes.c_double(0.)
                ymin = ctypes.c_double(0.)
                xmax = ctypes.c_double(0.)
                ymax = ctypes.c_double(0.)
                HData_var.ComputeRange(xmin, ymin, xmax, ymax)
                HStack_var.SetMaximum(math.ceil(max(HStack_var.GetMaximum(), ymax.value)))

                # === Set axis labels on the first histogram in stack ===
                

                # Draw stack
                HStack_var.Draw("HISTO")

                first_histo = h_list_var[0]  # must be non-null
                first_histo.GetXaxis().SetRangeUser(70., 500.)

                # Blinding for ZZMass
                if blindPlots and variable == "ZZMass_":
                    bblind = ROOT.TBox(105, 0, 140, HStack_var.GetMaximum() - epsilon)
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
                CMS_lumi.lumi_sqrtS     = r"171 fb^{-1} (13.6 TeV)"
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
