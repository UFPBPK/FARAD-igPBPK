# FARAD-igPBPK 
# This is a ReadMe file (Tutorial for the project code files).
# This repository contains all raw data and all code files for the project entitled "An Interactive Generic Physiologically Based Pharmacokinetic (igPBPK) Modeling Platform to Predict Drug Withdrawal Intervals in Cattle and Swine: A Case Study on Flunixin, Florfenicol and Penicillin G". Four folders are included in this project, including Code_App, Code_Artwork, Code_GetPC and Code_ModFit. 
# This manuscript is acceptable pending revision by Toxicological Sciences.

1.	Project folders

•	"Code_App": This folder is used to run APP.  

a)	R file for "Server" and "ui": These are R Shiny files that are used to launch the "App".  

b)	R file for "GenPBPK": This is the generic PBPK model code based on the Mrgsolve package.  

c)	R file for "Function": The R file includes required functions when the user runs the APP.  

d)	RDS file for "Fit_....rds": There are six rds files, which provide the optimized parameters.  

e)	Other files in this folder: All other files and the subfolders are used to support the App. 

•	"Code_Artwork": This folder contains all code files that are needed to reproduce results presented in all the figures and tables in the manuscript.  
a)	R files for "Figure_2, Figure_3, Figure_4, Figure_5_6, Figure_S1, Figure_S2": There are the R code files to reproduce all the artworks in the manuscript.  
b)	R file for “Table S4-7”: This is the R file used to generate results presented in Tables S4-S7 in the Supplementary Data file.  

c)	R file for "CirPlot": This R file is a function to produce the circle plot in Figure 4.  

d)	R file for "Pars": This R file is a function to run Monte Carlo simulations.  

e)	Excel files for "Data_FLO, Data_FLU, Data_PG": These Excel files provide the observed data extracted from the literature. 

f)	Rds files: All the rds files contain the simulation results and are used to support reproducing the artworks. 

g)	R file for "GenPBPK": This is the generic PBPK model code based on the Mrgsolve package.  

h)	R file for "GenPBPK": This is the generic PBPK model code based on the Mrgsolve package.  

•	"Code_GetPC": This folder is used to estimate the partition coefficients.  

a)	R file for "GetPC": This file is a function to get/estimate partition coefficients using different methods (listed below).  

b)	R file for "Kpfun_Berez, Kpfun_PT, Kpfun_Qppr, Kpfun_RR, Kpfun_Schmitt": These five files provide five different methods to estimate partition coefficients (references listed in the manuscript).

c)	Tiss_Cattle, Tiss_Pig, Tiss_Rat: These files contain tissue composition data in cattle, pigs, and rats, respectively. 

•	"Code_ModFit": This folder contains code files that are needed to do model calibration.  

a)	Excel files for "Data_FLO, Data_FLU, Data_PG": These Excel files provide the observed data extracted from the literature.  

b)	Excel files for "Tiss_Cattle, Tiss_Pig, Tiss_Rat": These files contain tissue composition data in cattle, pigs, and rats, respectively, which are needed to estimate partition coefficients.  

c)	R file for "Modfit_FLO_Cattle, Modfit_FLO_Swine, Modfit_FLU_Cattle, Modfit_FLU_Swine, Modfit_PG_Cattle, Modfit_PG_Swine": These files are for model calibration for each drug in each species.   

d)	R file for "Pars": This R file is a function to run Monte Carlo simulations.  

e)	R file for "GetPC": This file is a function to get/estimate partition coefficients using different methods.  

f)	R file for "Kpfun_Berez, Kpfun_PT, Kpfun_Qppr, Kpfun_RR, Kpfun_Schmitt": These five files provide five different methods to estimate partition coefficients (references listed in the manuscript).  

g)	R file for "GenPBPK": This is a generic PBPK model code based on the Mrgsolve package.  

2.	Instructions on the model code 

•	This instruction can be separated into three parts, including Part I: Model fitting; Part II: reproducing all results presented in the figures and tables in the manuscript and Part III: Run APP  

•	Before you run all the codes, you need to check (1) R version (make sure the R version > 4.0) and (2) install some required packages: “FME”, “ggplot2”, “mrgsolve”, “dplyr”, “tidyr”, “patchwork”, “ggprism”. 

2.1 Part I: Model fitting 

•	Download the ‘Code_ModFit” folder   

•	Open the RStudio and set working directory as the directory where the folder ‘Code_ModFit’ locates.

•	Run the R code for model fitting: Open the R file “Modfit_FLO_Cattle” (you can choose different one with the name starting with “Modfit_..) and run all the code.  

•	Get the results: After running all the code in the selected R file, the results of model fitting have been saved into R objects. For example, you can get the best parameters list if you run the code “exp(Fit_Cattle_FLO$par)”. You also can check the adjusted R-square by running the code “summary(fit)”. You also can check the Goodness-of-Fit plot for calibration data and evaluation data by run the code “p1” and “p2”.  

2.2 Part II: reproducing all results presented in the figures and tables in the manuscript 

•	Download the ‘Code_Artwork” folder   

•	Open the RStudio and set the working directory as the directory where the folder “Code_Artwork” locates.  

•	Open the R file “Figure_2.R” (you can choose any one you want to reproduce) and run all the code in the R file. 

•	Get the results: After running all the code in selected “Figure_2.R” file, the figure will saved in the your working directory. For example, if you run the code “Figure_2.R”, the “Fig 2a. tiff”, “Fig 2b1.tiff”, “Fig 2b2.tiff” and “Fig 2b3.tiff” will save in your working directory. 

2.3 Part III: Run App 

•	Download the ‘Code_App” folder  

•	Open the RStudio and set working directory as the directory where the folder ‘Code_App” locates.  

•	Double click the ‘ui.R’ (or ‘sever.R’, either one is ok) in the app file to open the app in RStuido.   

•	Once you open the ‘ui.R’ file, you can click the “Run App” button at the upper right of window to run the app on your local computer. Note that if this is your first time to launch the App, it may ask you to install necessary packages, such as “extrafont” and “ggh4x”. To run the app in the web browser to avoid abnormal issues, please click the ‘Open in Browser’ button to run the app with the http format (it still runs locally). Please also refer to the Supplementary Data file for a step-by-step tutorial for this App. 

