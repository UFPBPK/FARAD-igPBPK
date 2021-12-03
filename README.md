# FARAD-igPBPK
# This is a ReadMe file.
# This repository contains all raw data and all code for the project entitled "An Interactive Generic Physiologically Based Pharmacokinetic (igPBPK) Modelling Platform to Predict Drug Withdrawal Intervals in Cattle and Swine: A Case Study on Flunixin, Florfenicol and Penicillin G".
# This manuscript will be submitted to a peer-reviewed journal soon.

1. "Code_App": This folder is used to run APP
  a. R file for "Server" and "ui": These are R Shiny files that are used to run the "App".
  b. R file for "GenPBPK": This is a generic PBPK model code based on the Mrgsolve package.
  c. R file for "Function": The R file includes required functions when the user runs the APP.
  d. RDS file for "Fit_....rds": There are sixe rds files, which provide the optimized parameters.
  e. other: All the folders and files are used to support the App.


2. "Code_Artwork": This folder contains all code that is needed to reproduce all the figures present in the manuscript.
  a. R file for "Figure_2, Figure_3....": There are six R code files to reproduce all the artworks in the manuscript.
  b. R file for "CirPlot": The R file is a function to produe the circle plot in Figure 4.
  c. R file for "Pars": The R file is a function to run MC simulation.
  d. Excel files for "Data_...": The excel files provide the observed data extracted from the literature
  e. Rds files: All the rds files are used to support reprducing the artworks

3. "Code_GetPC": This folder is used to estimate the partition coefficients
  a. R file for "Get PC": The code is a function to get the partition coefficient from different methods.
  b. R file for "Kpfun_": There are six files to provide six different methods to estimate partition coefficients.

4. "Code_ModFit": Thif folder contains code files that are needed to do model calibration.
  a. Excel files for "Data_...": The excel files provide the observed data extracted from the literature.
  b. Excel files for "Tiss_": The tissue composition data which are needed to estimate partition coefficients.
  c. R file for "Modfit_..": The code for model calibration.
  d. R file for "Pars": The R file is a function to run MC simulations.
  e. R file for "Get PC": The code is a function to get the partition coefficients from different methods.
  f. R file for "Kpfun_": There are six files to provide six different methods to estimate partition coefficients.
  g. R file for "GenPBPK": This is a generic PBPK model code based on the Mrgsolve package.  
