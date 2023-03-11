# SECTION 02. This script imputes trait data using mean/mode replacement, k-nearest neighbour (KNN), missForest (RF), and multivariate imputation using chained equations (MICE) with and without phylogeny (in the form of phylogenetic eigenvectors derived from a phylogenetic tree). Given a complete-case dataset, it simulates data missing at random (MAR) for different taxonomic orders (using avian orders as an example here). It then uses different imputation methods to fill in the missing values using the known observations. The outputs for this script are measures of imputation performance for each trait in the dataset.

# Acknowledgments. ----

# All code/functions used in analyses adapted from: May et al. (2023). A real data-based simulation strategy to select an imputation method for mixed-type trait data. PLOS Computational Biology. In Press. 
# GitHub: https://github.com/jmay29/sim-imp-strategy

# This script makes use of the following imputation functions:

# kNN() function in the "VIM" package for imputation.
# Citations: Kowarik A, Templ M (2016). “Imputation with the R Package VIM.” Journal of Statistical Software, 74(7), 1–16. doi: 10.18637/jss.v074.i07.
# https://cran.r-project.org/web/packages/VIM/VIM.pdf

# missForest() function in the "missForest" package for imputation.
# Citations: Daniel J. Stekhoven (2013). missForest: Nonparametric Missing Value Imputation using Random Forest. R package version 1.4.
# Stekhoven D. J., & Buehlmann, P. (2012). MissForest - non-parametric missing value imputation for mixed-type data. Bioinformatics, 28(1), 112-118.

# mice() function in the "mice" package for imputation.
# Citations: van Buuren S, Groothuis-Oudshoorn K (2011). “mice: Multivariate Imputation by Chained Equations in R.” Journal of Statistical Software, 45(3), 1-67. https://www.jstatsoft.org/v45/i03/.
# R package version 3.13.0. https://cran.r-project.org/web/packages/mice/mice.pdf

### 1. Load libraries and functions. ----

library(ape)
library(data.table)
library(doParallel)
library(doRNG)
library(Information)
library(lmtest)
library(matrixStats)
library(mclust)
library(Metrics)
library(mice)
library(missForest)
library(mltools)
library(naniar)
library(nnet)
library(phytools)
library(plotrix)
library(MPSEM)
library(rsample)
library(simputation)
library(tidyverse)
library(VIM)
source("R/Functions/DataHandling_Functions.R")
source("R/Functions/Imputation_Functions.R")
source("R/Functions/Phylo_Functions.R")
source("R/Functions/Simpute_Functions.R")
source("R/Functions/ParlMICEWrapper.R") ## had to modify the parlmice function from the mice package due to a bug & to run in parallel.

### 2. Data loading and variable assignment. ----

# Trees. ---
# Read in trees.
treeFiles <- list.files(path = "Data/Trees/", pattern = "_ultra.tre")
l_trees <- lapply(treeFiles, function(x) read.tree(paste("Data/Trees/", x, sep = "")))
# Clean up tree file names.
treeNames <- gsub("_ultra.tre", "", treeFiles)
# Name l_trees according to file name as this contains both order and gene name information.
names(l_trees) <- treeNames
# Subset to bird examples.
l_trees <- l_trees[names(l_trees) %in% c("Apodiformes_COI", "Charadriiformes_COI", "Galliformes_COI")]

# Complete-cases datasets. ---
# NOTE: Upon creation of candidate datasets, code from May et al. (2023) was used to create near complete-case datasets to allow for some degree of missingness (to expand sample sizes).
# Link to script: https://github.com/jmay29/sim-imp-strategy/blob/main/R/Scripts/Optional_NearlyCompleteCaseCreation.R

# To quickly recreate these datasets, read in cleaned original trait datasets from previous script.
rawFiles <- c("Apodiformes.csv", "Charadriiformes.csv", "Galliformes.csv")
l_dfRaw <- lapply(rawFiles, fread, data.table = F)
# Subset dataframes according to traits that were selected criteria using May et al. (2023) strategy:
birdTraits <- list(c("latitude", "temp_range", "annual_temp", "precip_range", "annual_precip", "migration_1", "habitat",	"island",	"range_size",	"hwi", "for_aerial",	"body_mass_log"), c("territoriality",	"migration_3",	"diet",	"island",	"range_size",	"hwi",	"body_mass_log"), c("latitude",	"precip_range",	"annual_precip",	"temp_range",	"annual_temp",	"habitat",	"territoriality",	"migration_1",	"island",	"range_size",	"hwi",	"body_mass_log"))
# Subset dataframes in dfRaw according to traits.
l_dfRaw <- mapply(function(x, y){
  x <- x[, c("species_name", y)]
  }, x = l_dfRaw, y = birdTraits)
# Match dataframes to names in trees.
l_match <- mapply(DropAndMatch, l_trees, l_dfRaw, SIMPLIFY = F)
# Extract updated datasets.
l_dfCC <- lapply(l_match, function(x) x[[2]])
# Create vector of order names.
orderNames <- c("Apodiformes", "Charadriiformes", "Galliformes")
# Name dataframes according to order.
names(l_dfCC) <- orderNames
# Ensure blanks are NAs in all dataframes.
l_dfCC <- lapply(l_dfCC, function(x){
  x[x == ""] <- NA
  return(x)
})
# Write these to file.
mapply(fwrite, l_dfCC, file = paste(orderNames, "_NC.csv", sep = ""))

# Raw datasets. ---
# Name l_dfRaw datasets according to order.
names(l_dfRaw) <- orderNames
# Ensure blanks are NAs in all dataframes.
l_dfRaw <- lapply(l_dfRaw, function(x){
  x[x == ""] <- NA
  return(x)
})


# Traits. ---
# Extract trait names for each order.
l_traits <- lapply(l_dfCC, function(x) setdiff(colnames(x), "species_name"))
# Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical for each order.
l_traitTypes <- mapply(BreakIntoTypes, data = l_dfCC, traitCols = l_traits, SIMPLIFY = F)

# Extract numerical traits for each order.
l_contTraits <- lapply(l_traitTypes, function(x) x$Numerical)
# Extract categorical traits for each order.
l_catTraits <- lapply(l_traitTypes, function(x) x$Categorical)

# Create an integer for the number of replicates to use for missingness simulations and imputation. Default here is 100 replicates.
r <- 100
# Get a list of the functions loaded into the global environment.
l_functions <- mget(lsf.str())
# Index to get the simulation/imputation (simpute) function names.
index <- grep("Simpute", names(l_functions))
l_simputeFunctions <- l_functions[index]

### 3. Simulation/imputation. ----

# Set a seed to ensure results can be replicated.
set.seed(1591)
# Create empty list to hold results for each order.
l_l_simputeResults <- CreateNamedList(listLength = length(orderNames), elementNames = names(orderNames))

# For every order.
for(o in 1:length(orderNames)){
  
  # Take oth order.
  ord <- orderNames[[o]]
  # Get complete-case dataset for the order.
  dfCC <- l_dfCC[[grep(ord, names(l_dfCC))]]
  # Get raw dataset for the order.
  dfRaw <- l_dfRaw[[grep(ord, names(l_dfRaw))]]
  # Get trees for the order.
  l_ordTrees <- l_trees[grep(ord, names(l_trees))]
  # Update names for tree so only gene names are included.
  ordTreeNames <- gsub(paste(ord, "_", sep = ""), "", names(l_ordTrees))
  # Get trait names for the order.
  ordTraits <- l_traits[[grep(ord, names(l_traits))]]
  # Get numeric traits for the order.
  ordContTraits <- l_contTraits[[grep(ord, names(l_contTraits))]]
  # Get categorical traits for the order.
  ordCatTraits <- l_catTraits[[grep(ord, names(l_catTraits))]]
  
  # Create empty list to hold results for each simpute function.
  l_simputeResults <- CreateNamedList(listLength = length(l_simputeFunctions), elementNames = names(l_simputeFunctions))
  
  # For every simpute function..
  for(f in 1:length(l_simputeResults)){
    
    # Take fth function.
    simputeFunc <- l_simputeFunctions[[f]]
    # Get name of function.
    funcName <- names(l_simputeFunctions)[[f]]
    # Get name of imputation method.
    impMethod <- gsub(x = funcName, pattern = "Simpute", replacement = "")
    
    # Simpute without phylogeny.
    result <- simputeFunc(data = dfCC, raw = dfRaw, vars = ordTraits, int = r)
    
    # If the imputation method is MeanMode..
    if(impMethod == "MeanMode"){
      
      # Apply AverageErrors() function to write average error rates to file.
      AverageErrors(results = result$errorRates, data = dfCC, cont = ordContTraits, cat = ordCatTraits, method = impMethod, taxa = ord)
      
      # Append result to l_simputeResults.
      l_simputeResults[[f]] <- result
    } else {
      
      # Apply AverageErrors() function to write average error rates to file. Set paramTrack = T so we can track parameter values as well.
      res <- AverageErrors(results = result$errorRates, data = dfCC, cont = ordContTraits, cat = ordCatTraits, method = impMethod, taxa = ord, paramTrack = T)
    
      # If numerical traits were simulated..
      if(length(setdiff(ordContTraits, result$traitsNotSimulated)) > 0){
        # Apply AverageCorr() function to write average correlation coefficients to file. Set paramTrack = T so we can track parameter values as well.
      AverageCorr(results = result$corrCoef, param = res, data = dfCC, method = impMethod, taxa = ord)
      }
      
      # Simpute with phylogeny.
      phyResult <- mapply(simputeFunc, tree = l_ordTrees, MoreArgs = list(data = dfCC, raw = dfRaw, vars = ordTraits, int = r, phyImp = T), SIMPLIFY = F)
      # Extract error rates from result object.
      errors <- lapply(phyResult, function(x) x$errorRates)
      # Write average error rates to file.
      l_res <- mapply(AverageErrors, results = errors, treeName = ordTreeNames, MoreArgs = list(data = dfCC, cont = ordContTraits, cat = ordCatTraits, method = impMethod, taxa = ord, paramTrack = T), SIMPLIFY = F)
      
      # If numerical traits were simulated..
      if(length(setdiff(ordContTraits, result$traitsNotSimulated)) > 0){
      # Extract correlation coefficients from result object.
      corrs <- lapply(phyResult, function(x) x$corrCoef)
      # Write average error rates to file.
      mapply(AverageCorr, results = corrs, param = l_res, treeName = ordTreeNames, MoreArgs = list(data = dfCC, method = impMethod, taxa = ord), SIMPLIFY = F)
      }
      
      # Combine results.
      l_results <- list(withoutPhy = result, withPhy = phyResult) 
      # Append l_results to l_simputeResults.
      l_simputeResults[[f]] <- l_results
    }
    
  }

  # Append results to l_l_simputeResults.
  l_l_simputeResults[[o]] <- l_simputeResults
}
