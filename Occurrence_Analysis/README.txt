Code for the Bayesian hierarchical logistic regression models analyzing contaminant occurrence in smallmouth bass ovary and juvenile full-body tissue across various river sites in the Chesapeake Bay Watershed. 

All necessary data sets to run the occurrence analyses are found in the Occurrence Analysis folder and are as follows:

dat.final.p1.juv.csv is the data for the occurrence analysis that evaluates contaminant occurrence in only juvenile tissue. 

dat.final.p1.ovary.csv is the data for the occurrence analysis that evaluates contaminant occurrence in only ovary tissue. 

dat.final.p2.csv is the data for the occurrence analysis that compares the contaminant occurrence probability between ovary tissue and juvenile tissue.

Additional unedited data is found in the Data folder. 



Description of R Scripts:

Occurrence_Analysis.R includes the code for all three contaminant occurrence models (i.e., juvenile-only, ovary-only, and both juvenile and ovary combined). Each model is contained within two for-loops that loop over all the contaminants and then each land-use type. Code for pulling models estimates from the output for plotting and table creation is included within the script and within the for-loop for each analysis. 
