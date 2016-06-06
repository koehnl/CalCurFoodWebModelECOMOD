# CalCurFoodWebModelECOMOD
# Code and data for running analyses from Koehn et al. 2016 "Developing a high taxonomic resolution food web model to assess the #functional role of forage fish in the California Current ecosystem" Ecological Modelling

# Data
# Data for the Koehn et al. 2016 Ecopath model is contained in 3 files
# 1. Koehn.et.al.2016_groupinfo.csv has the basic Ecopath parameters 
# 2. Koehn.et.al.2016_parameters.csv again has the basic Ecopath parameters along with CV values for each parameter based on pedigrees for data quality (see Koehn et al. 2016)
# 3. Koehn.et.al.2016_dietmatrix.csv has the diet matrix for the Ecopath model

# Code
# The repository includes 2 Rcode files:
# 1. The Ecopath_Rcode_Koehn.et.al.2016.R contains the function "runecopath" for running an Ecopath model in R
# 2. Koehn.et.al.2016_MonteCarloCode.R has functions for running the Monte Carlo routine from Koehn et al. 2016. This includes finding probability distributions for each Ecopath parameters based on the specified CVs or specified pedigree for diets. This code also contains a function to make a set number of draws from the distributions for each parameter, including the diet matrix, and finding which set of parameters produces a balanced Ecopath model.
