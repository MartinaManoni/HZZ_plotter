# HZZ_plotter

Collection of plotting scripts for H→ZZ→4ℓ control and validation plots in Run 3.

## HZZ4l Mass Plotter

Generates m4l, Z1, and Z2 distributions for the full Run 3 H→ZZ→4ℓ analysis.

### Workflow

1. **ZX Background**  
   `zx.py` computes Z+X (fake-lepton) background yields for 2022–2024 from ROOT files.

2. **Histogram Filling**  
   `H4l_fill_RDF.py` reads MC, Data, and ZX files, applies event selection, computes weights, fills histograms (m4l, Z1, Z2) for all channels and the Higgs signal region, and saves ROOT files.

3. **Plotting**  
   `H4l_draw_RDF_fullRun3.py` combines yearly histograms, creates stacked plots, applies blinding in the Higgs mass region, and saves plots as ROOT files and images.

## Z Peak / Scale & Smearing Plots

Scripts: `HZZ_Zpeak_plot_2023.py` (2022–2023), `HZZ_Zpeak_plot_2024.py` (2024)

- Reads MC (DY, t̄t, WZ) and Data ROOT files.  
- Selects muon, electron, or inclusive channels.  
- Normalizes MC to data and includes statistical and systematic uncertainties.  
- Produces MC vs. Data plots with ratio and total uncertainty bands.  
- Saves results as PNG and PDF.

## Jets Control Region Plots

Scripts: `Jets_plots_scale_smear.C` (2022–2023), `Jets_2024.C` (2024)  

- Plots leading jet distributions (`η`, `pT`, `Njets`) for MC and data using ROOT.  
- Applies jet selections, computes event weights, and includes scale & smearing uncertainties.  
- Produces two-panel plots with data/MC comparison and ratio.  
- Saves ROOT histograms and plots (PDF/PNG) for specified Run 3 periods.


