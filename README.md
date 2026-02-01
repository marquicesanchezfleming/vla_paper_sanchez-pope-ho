Hello! If you are reading this, then you are perhaps interested in the work that we (Sanchez-Fleming, Pope, Ho) have done at Cornell University, between June 2024 and May 2026. In this repository, you will find 3 directories, titled "VLA_summer25", "Stage2-Plotting", and "Optical Properties". There general descriptions for each directory is as follows:
  
  VLA_summer25: This directory contains the bulk of the data analysis that was done during the summer of 2025. The main working scripts can all be found in 'main.ipynb', which is a pattern that will reoccur throughout the other directories. Inside, the bulk of what you can find is the code pertaining to the plotting of the radio luminosity light curves for the 17 radio-bright objects that we recieved extedned (S/C/X) band observations for, through project code VLA25A/386. Various other plots, including frequency vs flux, power law fits, and spectral energy distributions (SED's). 

  Stage2-Plotting: In this directory, we delve more into the actual plotting, as well as putting things together. As before, the main function is found in 'main.ipynb'. Inside, you will find code that can be used to create a LaTeX table of radio observations, histograms of SNe type in our BTS sample, archival observations of radio bright SNe, and radio detections from the Australian Square Kilometre Array Pathfinder (ASKAP) radio telescope. 

  Optical Properties: This directory is by far the most graphical, containing most of all the code neccessary to create the figures seen in our paper. Inside, you will find code required to construct the following plots:

    Optical Light Curves 
    Spectras
    Luminosity vs time plots (differentiates by H-rich and H-poor SNe)
    Temporal Index plots
    GRB constarint plots
    Time-Ordered Spectra
    Reorganized Stacking Plots

Remember, when in doubt, follow the 'main.ipynb'! There may (and likely is!) old csv's that are not particularly useful, but if there is something that can be created, then when I save a plot, I will typically save it as a PDF. 
