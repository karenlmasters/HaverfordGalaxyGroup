# HaverfordGalaxyGroup
Code collaboration for the Galaxy Group at Haverford College. 

Warning - this README is not automatic, so sometimes there's code here which is not documented/sorted corrected. 

* 1ADVICE.txt - advice to please read this README file! 

## Example Codes Either Written to help, or output of student papers

### Basic Statistics/Coding Techniques
* `BinnedMedianExample.ipynb` - functions to apply binned statistics to data with an example use on some random data. By Emmy Wisz with some modifications from Karen Masters. https://github.com/karenlmasters/HaverfordGalaxyGroup/blob/main/BasicCode/BinnedMedianExample.ipynb

* `Fitsinput_simpleplots.ipynb` -- Example code for how to read in a fits file, slice the data on selection from a column, make basic histograms and scatter plots. Written for mangaHIall.fits file... but should be easily modified to any .fits table.  https://github.com/karenlmasters/HaverfordGalaxyGroup/blob/main/BasicCode/Fitsinput_simpleplots.ipynb

* `basic_stats.py` -- Various useful statistics functions from Brooke Simmons. Includes a function to calculate proper errors on fractions, and to match two distributions (e.g. by mass). https://github.com/karenlmasters/HaverfordGalaxyGroup/blob/main/BasicCode/basic_stats.py

* `mpl_style.py` -- Style file for Matplotlib which may be useful (from Coleman Krawcyzk) https://github.com/karenlmasters/HaverfordGalaxyGroup/blob/main/BasicCode/mpl_style.py

* `SurvivalAnalysisTutorial.ipynb` -- Work in progress tutorial on Survival Analysis. Will likely use HI-MaNGA data as an example at some point. https://github.com/karenlmasters/HaverfordGalaxyGroup/blob/main/BasicCode/SurvivalAnalysisTutorial.ipynb


### Working with MaNGA Data
* `FourPanelGalaxyPlot.ipynb` -- Code to make a four panel plot for a MaNGA galaxy, showing the gri image and three other MaNGA maps. 
https://github.com/karenlmasters/HaverfordGalaxyGroup/blob/main/MaNGA/FourPanelGalaxyPlot.ipynb

* `MaNGAGalaxyofTheDay.ipynb` -- Code by Karen Masters which mades a four panel plot of any input MaNGA galaxy. Hoping to use for a Bluesky bot for 'MaNGA galaxy of the day'. There are enough galaxies in the MaNGA sample to do this for almost 30 years. https://github.com/karenlmasters/HaverfordGalaxyGroup/blob/main/MaNGA/MaNGAGalaxyofTheDay.ipynb

* `Pipe3D_access.ipynb` -- example notebook to access some Pipe3D analysis of MaNGA data (by KLM) https://github.com/karenlmasters/HaverfordGalaxyGroup/blob/main/MaNGA/Pipe3d_access.ipynb

* `RadialGradientsMaNGAEditedMar20.ipynb` -- create a radial metallicity trend for one or multiple galaxies. Originally from KLM, minor changes by Emmy Wisz. https://github.com/karenlmasters/HaverfordGalaxyGroup/blob/main/MaNGA/RadialGradientsMaNGAEditedMar20.ipynb

* `Metallicity_BinnedStatisticOverSamples.ipynb` -- sample selection for spiral galaxies binned into arm winding level, loops over all galaxies in a sample and saves a metallicity trend for each. From Emmy Wisz BMC'23 senior thesis https://github.com/karenlmasters/HaverfordGalaxyGroup/blob/main/MaNGA/Metallicity_BinnedStatisticOverSamples.ipynb

* `Metallicity_MassBinning_LinearRegressionFits_FitsofFits.ipynb` -- mass binning within subsamples, fitting a straight line to trends. From Emmy Wisz '23 senior thesis. https://github.com/karenlmasters/HaverfordGalaxyGroup/blob/main/MaNGA/Metallicity_MassBinning_LinearRegressionFits_FitsofFits.ipynb

* `Radial_Metallicity_Rainbow_ex.ipynb` - the rainbow plot of metallicity trends in bins of galaxy stellar mass. From Emmy Wisz. https://github.com/karenlmasters/HaverfordGalaxyGroup/blob/main/MaNGA/Radial_Metallicity_Rainbow_ex.ipynb


#### Working with Galaxy Zoo 3D Masks (and MaNGA Data)
*  `VAC_GZ3D_tutorial_KLMedits.ipynb` -- edits to the GZ3D VAC tutorial in Marvin so that it works on Sciserver.org (by KLM) https://github.com/karenlmasters/HaverfordGalaxyGroup/blob/main/GZ3D/VAC_GZ3D_tutorial_KLMedits.ipynb

*  `Arm_Interarm_Tutorial.ipynb` -- tutorial to learn how to calculate excess in stellar mass density of spiral galaxy arms. By Maša Kilibarda '26 with edits by KLM https://github.com/karenlmasters/HaverfordGalaxyGroup/blob/main/GZ3D/Arm_Interarm_Tutorial.ipynb 

* `KLM_RRL_PAandLength_Notebook.ipynb` -- Rachel Langgin BMC'23 uploaded her code to measure bar position angles from GZ:3D masks  https://github.com/karenlmasters/HaverfordGalaxyGroup/blob/main/GZ3D/KLM_RRL_PAandLength_Notebook.ipynb


### Working with HI-MaNGA Data/HI Data
* `HI-MassFraction-StellarMassPlot.ipynb` -- Simplified from `HIMaNGA_First_Paper_EmilyHarrington_KLM.ipynb` - just makes the HI mass fraction against stellar mass plot. https://github.com/karenlmasters/HaverfordGalaxyGroup/blob/main/HIMaNGA/HI-MassFraction-StellarMassPlot.ipynb 

*  `HIMaNGA_First_Paper_EmilyHarrington_KLM.ipynb` -- Code to make plots in Masters et al. 2019 on HI-MaNGA. Mostly written by Emily Harrington BMC '20. https://github.com/karenlmasters/HaverfordGalaxyGroup/blob/main/HIMaNGA/HIMaNGA_First_Paper_EmilyHarrington_KLM.ipynb

* `HI-Estimate Code Template.ipynb`-- Code by Elisabeth Brann BMC'26 which applied different techniques (model estimates and survival analysis) to HI scaling relations. https://github.com/karenlmasters/HaverfordGalaxyGroup/blob/main/HIMaNGA/HI-Estimate%20Code%20Template.ipynb

* `Model_HI_GlobalProfiles.ipynb` - Toy model for a HI global profile based on a model rotation curve and HI radial profile (by Karen Masters, with some help from Jessica Washington KNAC student summer 2021) https://github.com/karenlmasters/HaverfordGalaxyGroup/blob/main/HIMaNGA/Model_HI_GlobalProfiles.ipynb

* `SurvivalAnalysisTutorial.ipynb` -- Work in progress tutorial on Survival Analysis. Will likely use HI-MaNGA data as an example at some point. See below for instructions on how to set up a conda environment for this notebook. https://github.com/karenlmasters/HaverfordGalaxyGroup/blob/main/BasicCode/SurvivalAnalysisTutorial.ipynb



## Other Useful Repositories

ASTR352 (Extragalactic Data Science Class) https://github.com/karenlmasters/ASTR352JupyterActivities

Masha's Spiral Arm Project Respository: https://github.com/mashakilibarda/Spiral_Arms_Mass_Project 

Research Group Members also might like to check out Group Google Drive, which has some more private code sharing stuff. 

David Stark's survival analysis code (https://github.com/dvstark/survival/tree/main)

## Instructions for setting up an Conda environment for survival analysis codes

Note: These steps were performed using an Apply Macbook Pro with M1 processor and Sequoia 15.5

0) Install R from here:
	https://cran.rstudio.com/

1) This may be optional, but here is how conda was set up in this working example:

	-Install miniconda. Anaconda probably fine too, just takes longer. 

	-Add conda-forge to the conda channels:
		>> conda config --add channels conda-forge

	-Install a few packages that improve conda's speed when solving environments, and improve how it handles packages installed via pip
		>> conda install -n base conda-libmamba-solver
		>> conda config --set solver libmamba
		>> conda config --set pip_interop_enabled true 

2) Setup new conda environment w/ some standard packages. I called my environment rpy2. We're not actually installing rpy2 yet though. Conda-forge installs an older version that seems to force use of an older version of R that's incompatible with packages we need…

	>> conda create -n rpy2 python=3.12 jupyterlab numpy matplotlib lifelines fortranformat

3) Activate the environment

	>> conda activate rpy2

4) Now, install rpy2. Set to the latest version

	>> pip install rpy2==3.6.1

5) Check the version of R. You can do this just by opening R and seeing what it prints. This is what I get:
	>> R
	R version 4.5.1 (2025-06-13) -- "Great Square Root"
	 Copyright (C) 2025 The R Foundation for Statistical Computing
	 Platform: aarch64-apple-darwin20


Did not install R packages properly within python…trying inside R

6) open python and import rpy2. Make sure there are no errors

	>> import rpy2
