# SECTION 03. This script imputes the original dataset with missing values based on the results of simulations/imputations(Section 02).

# Required input for this script:
# dfCC = Cleaned complete-case dataset.
# dfRaw = Cleaned original dataset (with missing values).
# l_trees = Phylogenetic trees in Newick format (if using phylogenetic imputation).
# Error rate csv files for each trait/missing pattern.
# Best parameter csv files for each trait/missing pattern.

# Output:
# "DescriptiveStats.csv" = descriptive statistics of complete, original, and imputed numerical trait data.
# "ImputedData.csv" = final imputed dataset using best method

### Acknowledgments. ----

# Script adapted from https://github.com/jmay29/sim-imp-strategy/blob/main/R/Scripts/Section_3_FinalImputation.R.

### 1. Load libraries and functions. ----

library(ape)
library(data.table)
library(lmtest)
library(Information)
library(matrixStats)
library(Metrics)
library(mice)
library(missForest)
library(mltools)
library(naniar)
library(nnet)
library(pastecs)
library(phytools)
library(plotrix)
library(MPSEM)
library(rsample)
library(simputation)
library(tidyverse)
library(VIM)
library(viridis)
source("R/Functions/DataHandling_Functions.R")
source("R/Functions/Imputation_Functions.R")
source("R/Functions/Phylo_Functions.R")
source("R/Functions/Simpute_Functions.R")

### 2. Data loading and variable assignment. ----
# Set the seed.
set.seed(547)

# Read in the cleaned original trait dataset.
dfRaw <- fread("Apodiformes.csv", data.table = F)
# Ensure blanks are NAs.
dfRaw[dfRaw == ""] <- NA
# Read in cleaned complete-case trait dataset (using the PlotCompleteCase traits that were simulated MAR).
dfCC <- fread("Apodiformes_NC.csv", data.table = F)
# Ensure blanks are NAs.
dfCC[dfCC == ""] <- NA

# Subset to only include traits that were successfully simulated MAR (according to error rate files).
traitsSIM <- c("temp_range",	"annual_temp",	"precip_range",	"annual_precip",	"migration_1",	"habitat",	"for_aerial")
dfCC <- dfCC[, c("species_name", traitsSIM)]
dfRaw <- dfRaw[, c("species_name", traitsSIM)]
# Check dataframes to ensure they only have trait and taxonomy information!
colnames(dfCC)
colnames(dfRaw)

# Variable assignment. ---
# Assign order name.
orderName <- "Apodiformes"
# Look at column names of dfCC.
colnames(dfCC)
# Create vector of column names that contain taxonomic information.
taxCols <- c("species_name")
# Extract trait names from dfRaw.
traits <- setdiff(colnames(dfRaw), taxCols)
# Extract trait names from dfCC (these were the traits that were simulated MAR).
traitsMAR <- setdiff(colnames(dfCC), taxCols)
# Convert habitat to character.
dfCC$habitat <- as.character(dfCC$habitat)
dfRaw$habitat <- as.character(dfRaw$habitat)
# Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
l_traits <- BreakIntoTypes(dfRaw, traits)
l_traitsMAR <- BreakIntoTypes(dfCC, traitsMAR)
# Extract numerical traits.
contTraits <- l_traits[[1]]
contTraitsMAR <- l_traitsMAR[[1]]
# Extract categorical traits.
catTraits <- l_traits[[2]]
catTraitsMAR <- l_traitsMAR[[2]] ## None for this order
# Ensure correct format for catTraits.
dfRaw[, catTraits] <- lapply(dfRaw[, catTraits], as.factor)
dfCC[, catTraitsMAR] <- lapply(dfCC[, catTraitsMAR], as.factor)

### 3. Handling error rate and parameter files. ----

# Error rate file handling. ---
# Section 2 wrote error files for each trait and imputation method. Now we must read those files into R. For example, I have relocated the error files to a folder called Results/ErrorRates/All/.
# Read in the error rate files as dataframes.
errorFiles <- list.files(path = "Results/ErrorRates/", pattern = paste(orderName, ".+", "ErrorRates.csv", sep = ""))
l_dfErrors <- lapply(paste("Results/ErrorRates/", errorFiles, sep = ""), fread, data.table = F)
# Name according to file.
names(l_dfErrors) <- errorFiles

# Create a list to hold all of the error rate dataframes for each trait simulated MAR.
l_dfTraitAll <- CreateNamedList(listLength = length(traitsMAR), elementNames = traitsMAR)

# For every trait..
for(t in 1:length(traitsMAR)){
  # Take tth method.
  trait <- traitsMAR[[t]]
  # Index to get the corresponding dfErrors for that method.
  index <- grep(pattern = trait, x = names(l_dfErrors))
  l_dfTraitErrors <- l_dfErrors[index]
  names(l_dfTraitErrors)
  # Replace every 3rd column name with "error_rate" and 4th column name with SE so we can easily merge categorical/numerical traits together. Also add name o 
  l_dfTraitErrors <- lapply(l_dfTraitErrors, function(x){
    colnames(x)[3] <- "error_rate"
    colnames(x)[4] <- "SE"
    return(x)
  })
  # Rbind dataframes using rbindlist. Setting idcol = T so we can add a trait column as well.
  dfTraitError <- rbindlist(l_dfTraitErrors, idcol = T)
  # Append to l_dfTraitAll.
  l_dfTraitAll[[t]] <- dfTraitError
  
}

# Parameter file handling. ---
# Repeat same steps to handle the parameter files.
# Read in the parameter files as dataframes.
paramFiles <- list.files(path = "Results/ErrorRates/", pattern = paste(orderName, ".+", "_Parameters.csv", sep = ""))
l_dfParams <- lapply(paste("Results/ErrorRates/", paramFiles, sep = ""), fread, data.table = F)
# Name according to file.
names(l_dfParams) <- paramFiles
# Create a list to hold parameter dataframes for each trait.
l_dfParamsAll <- CreateNamedList(listLength = length(traitsMAR), elementNames = traitsMAR)

# For every trait..
for(t in 1:length(traitsMAR)){
  # Take tth method.
  trait <- traitsMAR[[t]]
  # Index to get the corresponding l_dfParamsAll for that trait.
  index <- grep(pattern = trait, x = names(l_dfParams))
  l_dfTraitParams <- l_dfParams[index]
  names(l_dfTraitParams)
  
  # Rbind dataframes using rbindlist. Setting idcol = T so we can add a trait column as well.
  dfTraitParams <- rbindlist(l_dfTraitParams, idcol = T)
  # Append to l_dfParamsAll.
  l_dfParamsAll[[t]] <- dfTraitParams
  
}


### 4. Determining best imputation method and parameter values for the dataset. ----

# Create lists to hold winner votes for each trait.
l_winningImp <- CreateNamedList(listLength = length(traitsMAR), elementNames = traitsMAR)

# For every trait..
for(t in 1:length(traitsMAR)){
  # Take tth trait.
  trait <- traitsMAR[[t]]
  # Index to get the corresponding error dataframes for the trait.
  dfTraitErrorAll <- l_dfTraitAll[[grep(trait, names(l_dfTraitAll))]]
  # Create trait columns.
  dfTraitErrorAll$trait <- trait
  # Create method column by removing trait, method, and file extension from .id column.
  dfTraitErrorAll$method <- gsub(pattern = trait, replacement = "", x = dfTraitErrorAll$.id)
  dfTraitErrorAll$method <- gsub(pattern = paste(orderName, "__", sep = ""), replacement = "", x = dfTraitErrorAll$method)
  dfTraitErrorAll$method <- gsub(pattern = "_ErrorRates.csv", replacement = "", x = dfTraitErrorAll$method)
  
  # If trait is numerical..
  if(trait %in% contTraitsMAR){
    # Determine the method that resulted in the lowest error rate.
    l_winningImp[[t]] <- dfTraitErrorAll$method[which.min(dfTraitErrorAll$error_rate)] 
  } else if (trait %in% catTraitsMAR){
    # Determine the method that resulted in the highest ARI.
    l_winningImp[[t]] <- dfTraitErrorAll$method[which.max(dfTraitErrorAll$error_rate)]
  }
}

# Determine the winning method!
BEST <- names(sort(table(unlist(l_winningImp)), decreasing = T))[1]

# Let's find the parameter values that performed best for our optimal method.

# For every dataframe in l_dfParamsAll..
for(t in 1:length(l_dfParamsAll)){
  # Take tth trait.
  trait <- traitsMAR[[t]]
  # Extract parameter dataframe.
  dfParamsMethod <- l_dfParamsAll[[t]]
  # Create trait column.
  dfParamsMethod$trait <- trait
  # Add a method column by cleaning up .id column.
  dfParamsMethod$method <- gsub(pattern = trait, replacement = "", x = dfParamsMethod$.id)
  dfParamsMethod$method <- gsub(pattern = paste(orderName, "__", sep = ""), replacement = "", x = dfParamsMethod$method)
  dfParamsMethod$method <- gsub(pattern = "_Parameters.csv", replacement = "", x = dfParamsMethod$method) 
  # Replace in l_dfParamsAll.
  l_dfParamsAll[[t]] <- dfParamsMethod
}

# Bind all dataframes in l_dfParamsAll.
dfParamsAll <- bind_rows(l_dfParamsAll)

# Subset to only those parameters used for our optimal method.
dfParamsOpt <- dfParamsAll[dfParamsAll$method == BEST, ]
# Extract parameter values.
bestParams <- dfParamsOpt$param
# Name according to traits.
names(bestParams) <- dfParamsOpt$trait

### 5. Imputation prep. ----

# Make copy of dfRaw.
dfRawSubset <- dfRaw

# Let's check the distributions of the data and for class imbalances in categorical data.
# Descriptive info for numerical data.
contRes <- lapply(dfRawSubset[, contTraits], GetNumericalInfo)
contRes
# Descriptive info for categorical data:
catRes <- lapply(dfRawSubset[, catTraits], GetCategoricalInfo)

# For every categorical trait..
for(cat in 1:length(catRes)){
  # Print the name of the trait.
  cat("Trait:", names((catRes))[[cat]])
  # Extract the result.
  result <- catRes[[cat]]
  # Print the results.
  print(result[1:4])
  # Convert frequency count into format amenable to plotting.
  dfCount  <- data.frame(result[[3]])
  # Plot barplot.
  barplot(height = dfCount$Freq, names.arg = dfCount$vector, col = "skyblue", main = paste(names((catRes))[[cat]], sep = " "))
}

# # Apodiformes EDITS:
# Remove open (3) for habitat because there only a few obs for these groups (less than 10%).
table(dfRawSubset$habitat)
dfRawSubset <- dfRawSubset[-which(dfRawSubset$habitat == 3), ]
table(dfRawSubset$habitat)
# Drop unused levels.
dfRawSubset$habitat <- droplevels(dfRawSubset$habitat)

# Selecting predictors. ---
# Here, we select those traits that have significant correlations to use as predictors for imputation (this will vary for each trait).
# Apply the SelectPredictors() function to obtain a list of predictors for each trait.
l_predictors <- SelectPredictors(dfRawSubset[, traits])

# Normalize numeric data. ---
# Make copy of data so we can normalize the numeric data while still retaining original data.
dfNorm <- dfRawSubset
# Use lapply to scale numerical traits.
dfNorm[, contTraits] <- lapply(dfNorm[, contTraits], scale)
dfNorm[, contTraits] <- lapply(dfNorm[, contTraits], as.numeric)

# Phylogenetic imputation prep. ---
# Set phylogenetic imputation to true if proceeding (default is F):
phyImp <- F

### 6. Final imputation. ----

# Replicate dfNorm so we can impute the missing values in a new dataset.
dfNormImp <- dfNorm

# For each trait..
for(t in 1:length(traits)) {
  # Take the name of the tth trait.
  trait <- traits[[t]]
  # If phylogenetic imputation was selected...
  if(phyImp == T){
    # Take the corresponding predictors with appended eigenvectors for that trait.
    preds <- l_EVPredictors[[grep(trait, names(l_EVPredictors))]]
  } else {
    # Take the corresponding predictors for that trait.
    preds <- l_predictors[[grep(trait, names(l_predictors))]]
  }
  # Take the best parameter for the trait.
  param <- bestParams[[grep(trait, names(bestParams))]]
  # If phylogenetic imputation was selected...
  if(phyImp == T){
    # Take the tth dfNormEV.
    dfNormMiss <- l_dfNormEV[[grep(trait, names(l_EVPredictors))]] 
  } else {
    # Make a copy of dfNorm to impute.
    dfNormMiss <- dfNorm
  }
  # If optimal method is KNN..
  if(grepl("KNN", BEST)){
    # Impute using kNN and optimal value of k.
    dfImputed <- kNN(dfNormMiss, variable = trait, dist_var = preds, k = param)
    # If optimal method is RF..
  } else if(grepl("RF", BEST)){
    # Impute the dataset (which contains the trait in question and its predictors) using missForest and optimal value of ntree.
    imputedRF <- missForest(as.data.frame(dfNormMiss[, c(trait, preds)]), ntree = param)
    # Access the imputed dataframe.
    dfImputed <- imputedRF$ximp
    # Add species_name back.
    dfImputed$species_name <- dfNormMiss$species_name
    # Round the imputed count data.
    #dfImputed[, intTraits] <- lapply(dfImputed[, intTraits], function(x) as.integer(round(x)))
    # If optimal method is MICE..
  } else if(grepl("MICE", BEST)){
    
    # If phylogenetic imputation was selected..
    if(phyImp == T){
      # MICE requires a predictor matrix for imputation. So, we will first merge all of the dataframe in l_dfNormEV into one dataframe. This is so we can identify the eigenvector columns needed for imputation of each trait and update this in the predictor matrix. This should also reduce computation time since we are just imputing one dataframe and not a list of dataframes. 
      dfMissEV <- l_dfNormEV[[1]]
      # Identify the eigenvector columns and remove those from dfMissEV.
      evIndex <- grep(pattern = "V_", colnames(dfMissEV))
      dfMissEV <- dfMissEV[, -evIndex]
      # Now, let's merge dfMissEV with the corresponding eigenvector columns in each dataframe in l_dfNormEV.
      for(e in 1:length(l_dfNormEV)){
        # Take the eth l_dfNormEV.
        dfMisseth <- l_dfNormEV[[e]]
        # Identify eigenvector columns in eth dataframe.
        index <- grep("V_", colnames(dfMisseth))
        evCols <- colnames(dfMisseth)[index]
        # Merge these columns with dfMisseth (including species_name).
        dfMissEV <- merge(dfMissEV, dfMisseth[, c("species_name", evCols)], by = "species_name")
      }
      # Create predictor matrix for use in MICE imputation with appended eigenvectors. We can indicate the column names of dfMissEV (excluding species_name) as names of variables to consider in imputation process. May take a while depending on size of dataset.
      matPredictors <- CreatePredictorMatrix(dfMissing = dfMissEV, cols = colnames(dfMissEV)[-1], predictors = l_EVPredictors)
      # Ensure all of the rows for the eigenvectors to 0 as they themselves will not be imputed.
      eigens <- grep(pattern = "V_", rownames(matPredictors))
      matPredictors[eigens, ] <- 0
      # Impute the datasets using optimal m and matPredictors (can take a while depending on size of dataset). May take a while depending on size of dataset.
      imputedMICE <- mice(dfMissEV[, c(colnames(matPredictors))], predictorMatrix = matPredictors, m = param, maxit = 10, print = FALSE)
      # Combine the imputed data into one dataframe.
      dfImputed <- CombineMIDataframes(midata = imputedMICE, method = "MICE", m = param, contVars = contTraits, catVars = catTraits)
    } else {
      # Create predictor matrix for use in MICE imputation (based on results of our trait predictor screening).
      matPredictors <- CreatePredictorMatrix(cols = traits, predictors = l_predictors, dfMissing = dfNormMiss)
      # Impute the datasets using optimal m and matPredictors (can take a while depending on size of dataset).
      imputedMICE <- mice(dfNormMiss[, c(colnames(matPredictors))], predictorMatrix = matPredictors, m = param, maxit = 10, print = FALSE)
      # Combine the imputed data into one dataframe.
      dfImputed <- CombineMIDataframes(midata = imputedMICE, method = "MICE", m = param, contVars = contTraits, catVars = catTraits)
    }
  }
  # Ensure dfNormImp and dfImputed are in same order.
  dfNormImp <- dfNormImp[order(dfNormImp$species_name), ]
  dfImputed <- dfImputed[order(dfImputed$species_name), ]
  # Replace the data for the trait in dfNormImp with the newly imputed data in dfImputed.
  dfNormImp[, trait] <- dfImputed[, trait]
}

# Ensure dataframes are in the same order.
dfCC <- dfCC[c("species_name", traitsMAR)]
dfNorm <- dfNorm[c("species_name", traitsMAR)]
dfNormImp <- dfNormImp[c("species_name", traitsMAR)]
# Combine dataframes into list.
l_dfAll <- list(dfCC, dfNorm, dfNormImp)
# Name the list according to df.
names(l_dfAll) <- c("dfCC", "dfNorm", "dfNormImp")
# Now we want to edit the column names in each dataframe so we can easily identify which dataset the trait data are from. First, create list of strings to append to each col name (in order of l_dfAll).
colStrings <- list("cc", "raw", "raw_imp")
# Now apply PasteColNames() using mapply. Because we wrapped traits in a list() it is recycled for each iteration (our constants). Specified SIMPLIFY = F so it doesn't reduce results to a matrix.
l_dfAll <- mapply(PasteColNames, data = l_dfAll, colNames = list(traitsMAR), string = colStrings, sep = "_", SIMPLIFY = F)
# Merge all the dataframes by species name.
dfAll <- Reduce(function(...) merge(..., by = "species_name", all = T), l_dfAll)
# Remove outgroup from dfAll.
dfAll <- dfAll[!dfAll$species_name == "Nyctidromus_albicollis", ]

# Backtransform contTraitsMAR data back so we can log-transform instead for visualization purposes.
# Calculate mean of each numeric trait in temp_range.
# TODO: Convert following lines to function.
avg <- mean(dfRawSubset$temp_range, na.rm = T)
# Calculate SD of each numeric trait in original dfRawSubset.
SD <- sd(dfRawSubset$temp_range, na.rm = T)
# Backtransform normalized and imputed normalized data using Unscale function.
dfAll$temp_range_raw <- Unscale(vector = dfAll$temp_range_raw, scale = SD, centre = avg)
dfAll$temp_range_raw_imp <- Unscale(vector = dfAll$temp_range_raw_imp, scale = SD, centre = avg)

# Calculate mean of each numeric trait in annual_temp.
avg <- mean(dfRawSubset$annual_temp, na.rm = T)
# Calculate SD of each numeric trait in original dfRawSubset.
SD <- sd(dfRawSubset$annual_temp, na.rm = T)
# Backtransform normalized and imputed normalized data using Unscale function.
dfAll$annual_temp_raw <- Unscale(vector = dfAll$annual_temp_raw, scale = SD, centre = avg)
dfAll$annual_temp_raw_imp <- Unscale(vector = dfAll$annual_temp_raw_imp, scale = SD, centre = avg)

# Calculate mean of each numeric trait in precip_range.
avg <- mean(dfRawSubset$precip_range, na.rm = T)
# Calculate SD of each numeric trait in original dfRawSubset.
SD <- sd(dfRawSubset$precip_range, na.rm = T)
# Backtransform normalized and imputed normalized data using Unscale function.
dfAll$precip_range_raw <- Unscale(vector = dfAll$precip_range_raw, scale = SD, centre = avg)
dfAll$precip_range_raw_imp <- Unscale(vector = dfAll$precip_range_raw_imp, scale = SD, centre = avg)

# Calculate mean of each numeric trait in annual_precip.
avg <- mean(dfRawSubset$annual_precip, na.rm = T)
# Calculate SD of each numeric trait in original dfRawSubset.
SD <- sd(dfRawSubset$annual_precip, na.rm = T)
# Backtransform normalized and imputed normalized data using Unscale function.
dfAll$annual_precip_raw <- Unscale(vector = dfAll$annual_precip_raw, scale = SD, centre = avg)
dfAll$annual_precip_raw_imp <- Unscale(vector = dfAll$annual_precip_raw_imp, scale = SD, centre = avg)

### 8. Descriptive statistics. ----
contCols <- sapply(dfAll, is.numeric)
# Remove categorical as we are obtaining descriptive stats for numerical traits.
dfContSubset <- dfAll[, contCols]
# Apply stat.desc to dfContSubset to get some descriptive stats for each trait and dataset.
dfDescCont <- stat.desc(dfContSubset)
# Order alphabetically.
dfDescCont <- dfDescCont[, order(colnames(dfDescCont))]
# Write to file.
write.csv(dfDescCont, "ApodiformesDescriptiveStats.csv")

# Order dfAll alphabetically and write to file.
dfAll <- dfAll[, order(colnames(dfAll))]
write.csv(dfAll, "ApodiformesFinalDataset.csv")
