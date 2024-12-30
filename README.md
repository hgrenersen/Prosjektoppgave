# Prosjektoppgave
Code written and data used for the work with writing my specialization project (TMA4500) at NTNU. 

Scripts containing functions to generate data, implement kriging, different CV (LOO, K-fold and block) estimators and estimating importance weights (Voronoi, KMM & KLIEP) can be found in the folder "Scripts". Scripts that actually generate data, calculate estimates and store these things, are placed at the top level of the repository. The Data folder is sorted based on the value of the LGCP parameter, with each folder containing the data used, the estimated IWs, the true IWs and quantities needed to calculate the classical CV estimates. Unfortunately, the underlying realization which is sampled to generate the data, and the intensity field, are not available on Github, since the files are too large.

Here is also a link to a visualization of the estimated IWs for some datasets: https://ao3gaf-hgrenersen.shinyapps.io/ImportanceWeights/

Explanations of the files
- LOOCV_Gridded.R: Script to generate gridded data, and calculate LOOCV estimates as well for evaluating their performance.
- generateDataAndDoCV.R: The script generates the data used in all simulation studies, i.e. for values of the LGCP parameter beta in {0, 0.5, 1, 1.5, 2}. It also calculates classical CV estimates as LOO, K-fold and block for all of the datasets.
- IW_timing.R: Estimates (and times this) the IWs for all types of IWs considered, across all different datasets
- calcTrueQuants.R: Calculates the main quantity of interest, the average ECRPS over a 200 by 200 grid, which the CV estimators might be thought of as estimators of.
- calcAllFuncs.R: Includes additional functions needed, as for tuning the bandwidth to be used in KMM. 
- Tables.Rmd: Gathers results from CV and IWCV estimators and compares them to the average ECRPS in terms of MSE. 
