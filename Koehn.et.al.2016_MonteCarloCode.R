# Koehn et al. 2016 Monte Carlo routine
# Code to read in parameter values and CVs and diet and pedigree, create distributions for parameters, and pull out new Ecopath values
# Also has a function to make some set amount of Monte Carlo draws and run a set amount of Ecopath models to see 
# which models are balanced.

# READ IN DATA (Example with Koehn et al. 2016 California Current Model)

workdir <- # set directory as location of Koehn.et.al.2016_groupinfo, Koehn.et.al.2016_dietmatrix, and Koehn.et.al.2016_parameters
setwd(workdir)

BasicParams <- read.csv(file = "Koehn.et.al.2016_groupinfo.csv", header = TRUE)
Type = BasicParams$GroupType
BA = BasicParams$BA
# Read in Biomass, PB, QB, EE, and CV values for each
params <- read.csv(file = "Koehn.et.al.2016_parameters.csv", header = TRUE) 
func.groups = params$Group
density.vec = params$density; CV.vec = params$CV_bio # initial densities and density CVs
PB.vec = params$PB; PB_CV.vec = params$PB_CV # initial PB values
QB.vec = params$QB; QB_CV.vec = params$QB_CV # initial QB values
EE.vec = params$EE  # initial ecotrophic efficiency values
yield.vec = params$Yield; Y_CV.vec = params$Y_CV


######################################################################
####### Diet ##########
#####################################################################
diet <- read.csv(file = "Koehn.et.al.2016_dietmatrix.csv", header = TRUE)
pedigree = diet$pedigree[1:92] # scores for the quality of each diet
diet_use = diet[,2:94] # removes group names
pedigree_multiplier = as.numeric(diet$pedigree_mult[1:94]) # values to multiply diet pedigree by to produce ranges of values
# in the dirichlet distribution

library(MCMCpack)
pedigree.convert <- function(pedigree) {
  multiplier = vector(mode = 'numeric', length = length(pedigree))
  for(i in 1:length(pedigree)) {
    if(pedigree[i] == 1) {
      multiplier[i] = 2000 # 900, 1200
    } else if (pedigree[i] == 0.8) {
      multiplier[i] = 1500 # 600, 900
    } else if (pedigree[i] == 0.6) {
      multiplier[i] = 1000 #300 , 600
    } else if (pedigree[i] == 0.4) {
      multiplier[i] = 500 # 200 300
    } else if (pedigree[i] == 0.2) {
      multiplier[i] = 300 # 100 200
    } else {
      multiplier[i] = 150 #50 100
    }
  }
  return(multiplier)
}

# function to pull 1 new diet matrix given original diet matrix and the pedigree
LEKdirichlet <- function(diet.matrix, pedigree_multiplier) {
  n = ncol(diet.matrix)
  # need 1 empty diet matrix
  matrixnew = matrix(NA, nrow = nrow(diet.matrix), ncol = n)
  #matrixnew = diet.matrix
  for (i in 2:(n-1)) {
    diet =  diet.matrix[,i]*pedigree_multiplier[i] # get diet for one group
    run = rdirichlet(1, diet)
    matrixnew[,i] = run
  } 
  matrixnew[,1] = 0 
  matrixnew[,93] = 0 # for detritus
  return(matrixnew)
}


#############################################################################################################
# Probability distributions needed

lognormal <- function(mean, CV) {
  meanlog = log(mean)-((CV^2)/2)  # mu = log(E[X]) - (sigma^2)/2
  return(params = list(meanlog = meanlog, sdlog = CV))
}

#############################################################################################################
# New density values from probability distribution

density.values.noloop <- function(density.vec, CV.vec) {
  index = which(density.vec > -1)
  n = length(index)
  new.density = density.vec 
  mean = density.vec 
  CV = CV.vec
  temp = lognormal(mean[index], CV[index])
  new.density[index]<-rlnorm(n,temp$meanlog, temp$sdlog)
  new.density[length(density.vec)] = 10
  return(new.density)
}


################################################################################################################################
#Function to pull new PB and QB values from probability distributions constructed from the original PB or QB value and a CV. 
# THIS FUNCTION ASSUMES THAT THERE IS A VALUE FOR PB and QB FOR ALL FUNCTIONAL GROUPS - NOT SOLVING FOR PB or QB*****

PBQB.values.noloop <- function(PB.vec, PB_CV.vec, QB.vec, QB_CV.vec) {
  indexPB = which(PB.vec > 0); indexQB = which(QB.vec > 0)
  n = length(indexPB); n2 = length(indexQB)
  new.PB = PB.vec; new.QB = QB.vec
  #mean = PB.vec; CV = PB_CV.vec
  tempPB = lognormal(PB.vec[indexPB], PB_CV.vec[indexPB])
  new.PB[indexPB] = rlnorm(n, tempPB$meanlog, tempPB$sdlog)
  tempQB = lognormal(QB.vec[indexQB], QB_CV.vec[indexQB])
  new.QB[indexQB] = rlnorm(n2, tempQB$meanlog, tempQB$sdlog)
  PB_QB = cbind(new.PB, new.QB)
  return(PB_QB)
}


################################################################################################################################
###########  EE ##################
# Uniform distribution function for EE
# All organisms where EE is entered (not solved for as the last parameter) have EE values ~= 0.8
# (except for phytoplankton with an EE of ~0.4, and dungeness crab ~0.64, and detritus/offal pools)
# All organisms with EE = 0.8, and dungeness crab, have EE's possibly ranging from 0.5-0.9
# phytoplankton range is from 0.25-0.75, and detritus/offal from 0-0.2
# Need to pull from different ranges, depending on what organism
# Phytoplankton is ALWAYS first 


EE.uniform.noloop <- function(EE.vec) {
  index = which(EE.vec > -1)
  n = length(index[-1])
  index = index[-1]
  new.EE = EE.vec
  new.EE[1] = runif(1, min = 0.25, max = 0.75)
  new.EE[index] = runif(n, min = 0.5, max = 0.9)
  return(new.EE)
}

############################################################################################
# Probability distribution for Yield values

Yield_noloop <- function(yield.vec, Y_CV.vec) {
  index = which(yield.vec > 0)
  n = length(index)
  new.yield = yield.vec
  mean = yield.vec 
  CV = Y_CV.vec
  temp = lognormal(mean[index], CV[index])
  new.yield[index]<-rlnorm(n,temp$meanlog, temp$sdlog)
  return(new.yield)
}

##########################################################################################
# Put each draw of Biomass, PB, QB, EE, and Yield into a different matrix of parameters
new_params <- function(density.vec, CV.vec, PB.vec, PB_CV.vec, QB.vec, QB_CV.vec, EE.vec, yield.vec, Y_CV.vec) {
  param_matrix = matrix(nrow = length(density.vec) , ncol = 5) 
  param_matrix[,1] = density.values.noloop(density.vec, CV.vec)
  temp = PBQB.values.noloop(PB.vec, PB_CV.vec, QB.vec, QB_CV.vec)
  param_matrix[,2] = temp[,1]; param_matrix[,3] = temp[,2]
  param_matrix[,4] = EE.uniform.noloop(EE.vec)
  param_matrix[,5] = Yield_noloop(yield.vec, Y_CV.vec)
  return(param_matrix)
}

# Run multiple ecopath models and see which are balanced - saves parameter and diet matrix files for balanced runs
# NOTE: Run "runecopath" function in file "Ecopath_Rcode_Koehn.et.al.2016.R" before running 
run_many_ecopath2 <- function(density.vec, CV.vec, PB.vec, PB_CV.vec, QB.vec, QB_CV.vec, EE.vec, yield.vec, Y_CV.vec, diet.matrix, pedigree_multiplier, size) {
  time1=Sys.time()
  Group = func.groups
  n = length(Group)
  BA = rep(0, length = n)
  test_balanced = matrix(nrow = n, ncol = size)
  for (i in 1:size) {
    param_matrix = new_params(density.vec, CV.vec, PB.vec, PB_CV.vec, QB.vec, QB_CV.vec, EE.vec, yield.vec, Y_CV.vec)
    matrixnew = LEKdirichlet(diet.matrix, pedigree_multiplier)
    ecopathinput = list(Group = Group, Biomass = param_matrix[,1], PB = param_matrix[,2],
                        QB = param_matrix[,3], EE = param_matrix[,4], Yield = param_matrix[,5], BA = BA, Type = Type, DC = matrixnew)
    ecopathout=runecopath(ecopathinput) 
    test_balanced[,i] = ecopathout$EEout
    index = which(ecopathout$EEout[-93] > 1)
    if(length(index) == 0) {
      paramfile = paste('params',i,'.csv', sep = '')
      dietfile = paste('diet',i, '.csv', sep = '')
      write.csv(param_matrix, file = paramfile)
      write.csv(matrixnew, file = dietfile)
    }
  }
  timediff = Sys.time() - time1 # total time to run
  print(timediff)
  return(test_balanced)
}

set.seed(6) # Need to set.seed if running multiple groups of draws or else will pull the same parameters for each set of runs

test = run_many_ecopath2(density.vec, CV.vec, PB.vec, PB_CV.vec, QB.vec, QB_CV.vec, EE.vec, yield.vec, Y_CV.vec, diet.matrix = diet_use, pedigree_multiplier, 10)

