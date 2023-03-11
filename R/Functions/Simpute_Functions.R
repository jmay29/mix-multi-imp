# Functions for simulating missingness at random (MAR) and imputing values in complete-case datasets. Returns lists of error rates for both numerical and categorical variables.

MeanModeSimpute <- function(data, raw, vars, int = 100) {
  
  # Given a complete-case dataset, this function simulates missingness at random (MAR) and imputes values using mean/mode replacement. Determines error rates for numerical (mean squared error - MSE) and categorical variables (proportion falsely classified - PFC).
  # data = complete-case dataset containing trait data. Must also contain a column with species name information = "species_name"
  # raw = original dataset containing trait data with missing values. Used to build logistic regression models and simulate data that are MAR in the complete-case dataset. Must also contain a column with species name information = "species_name"
  # vars = names of columns for which to simulate missing data
  # int = number of iterations (missingness replicates). Default is 100
  
  # Trait preparation. ---
  # Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
  l_traits <- BreakIntoTypes(data, vars)
  # Extract numerical traits.
  contTraits <- l_traits[[1]]
  # Extract categorical traits.
  catTraits <- l_traits[[2]]
  # Ensure categorical traits are factor type for logistic regression model building (data class required for glm function).
  # If there is more than one categorical trait..
  if(length(catTraits) > 1){
    # Use lapply to convert categorical traits to factor class.
    data[, catTraits] <- lapply(data[, catTraits], as.factor)
    raw[, catTraits] <- lapply(raw[, catTraits], as.factor)
    # If there is only one categorical trait..
  } else if(length(catTraits) == 1) {
    # Convert single trait to factor class.
    data[, catTraits] <- as.factor(data[, catTraits])
    raw[, catTraits] <- as.factor(raw[, catTraits])
  }
  # Identify integer (count) traits, if any.
  intTraits <- GetTraitNames(data = data[, vars], class = "integer")

  # Determine original sample size for each trait.
  l_sampleSizes <- lapply(data[, vars], function(x) length(na.omit(x)))
  
  # Logistic regression fitting. ---
  # Apply FitLogReg() function to raw data to identify significant predictors of missingness for each trait.
  l_models <- FitLogReg(data = raw, cols = vars)
  # Now, let's identify traits that CAN be simulated MAR versus those that cannot (i.e. they must be simulated MCAR). Traits that must be simulated MCAR either 1) have no missing values in the original data or 2) do not have any predictors of missingness in the dataset upon fitting of the logistic regression models. IDMissPattern() also refines the logistic regression models (i.e. drops insignificant terms) through use of the RefineModels() function.
  l_missPatt <- IDMissPattern(data = raw, vars = vars, models = l_models)
  # Extract MCAR traits.
  traitsNOSIM <- l_missPatt[[1]]
  # Extract MAR traits.
  traitsMAR <- l_missPatt[[2]]
  # If traitsMAR is empty (no traits can be simulated MAR)..
  if(length(traitsMAR) == 0){
    stop("No traits can be simulated MAR!") 
  }
  # Extract final MAR models.
  l_MARfinalModels <- l_missPatt[[3]]
  
  # Imputation prep. ---
  # Create list to hold the dataframes with simulated missingness for each iteration.
  l_dfMiss <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the missingness proportion in each dfMiss.
  l_l_missingness <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the imputed dataframes for each iteration.
  l_dfImp <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create lists to hold the error rates for each iteration.
  l_Error <- CreateNamedList(listLength = int, elementNames = 1:int)
  
  # For every iteration...
  for(i in 1:int) {
    
    # set.seed.
    set.seed(i)
    # MAR simulation. ---
    # Create copy of complete-case dataframe (untransformed) to introduce NAs into.
    dfMissOrig <- data
    # Create a list to hold the missingness proportion for each trait to be simulated MAR.
    l_missingness <- CreateNamedList(listLength = length(traitsMAR), elementNames = traitsMAR)
    
    # For every trait that can be simulated MAR..
    for(t in 1:length(traitsMAR)){
      # Take tth trait.
      trait <- traitsMAR[[t]]
      # Get the corresponding MAR model.
      index <- grep(pattern = trait, x = names(l_MARfinalModels))
      modelMAR <- l_MARfinalModels[[index]]
      # Simulate MAR data in complete-case dataset.
      res <- SimMAR(model = modelMAR, data = data)
      # Get original sample size of trait. ++
      n <- l_sampleSizes[[grep(trait, names(l_sampleSizes))]]
      # Subtract data originally missing in data from the number missing in dfMissOrig after MAR simulation to determine actual number of NAs introduced. ++
      marN <- sum(is.na(res)) - sum(is.na(data[[trait]]))
      # If missingness is less than 0.80 (enough data for imputation)..++
      if(marN/n < 0.80){
        # Replace complete-case column in dfMissOrig with res. ++
        dfMissOrig[, trait] <- res 
        # Divide by n to determine missingness percentage for trait and append to l_missingness. ++
        l_missingness[[t]] <- marN/n
      } else {
        # Assign NA as missingness proportion value.
        l_missingness[[t]] <- NA
      }
    }
    
    # Finally, append traits that could not be simulated MAR (no NAs introduced based on fitted logistic regression models) or that exceed 0.80 missingness (hard time getting accurate imputations) to traitsNOSIM.
    traitsNOSIM <- c(traitsNOSIM, names(which(l_missingness == 0)), names(which(is.na(l_missingness))))
    # Update traitsMAR based on this result.  
    traitsMAR <- setdiff(traitsMAR, traitsNOSIM)
    # If traitsMAR is empty (no traits can be simulated MAR)..
    if(length(traitsMAR) == 0){
      stop("No traits can be simulated MAR!") 
    }
    
    # Make copy of data and rename. $$$
    dfNorm <- data
    # Normalize numerical traits prior to imputation. $$$
    # If there is more than one numerical trait..
    if(length(contTraits) > 1){
      # Use lapply to scale numerical traits.
      dfNorm[, contTraits] <- lapply(dfNorm[, contTraits], scale)
      dfNorm[, contTraits] <- lapply(dfNorm[, contTraits], as.numeric)
      # If there is only one numerical trait..
    } else if(length(contTraits) == 1) {
      # Scale single trait. ==
      dfNorm[, contTraits] <- as.numeric(scale(dfNorm[, contTraits]))
    }
    # Introduce missing values from dfMissOrig into dfNorm.
    dfMiss <- as.data.frame(mapply(function(x, y) ifelse(is.na(x), x, y), x = dfMissOrig[, c("species_name", contTraits)], y = dfNorm[, c("species_name", contTraits)], SIMPLIFY = F))
    # Add categorical traits back.
    dfMiss <- merge(dfMiss, dfMissOrig[, c("species_name", catTraits)], by = "species_name")
    # Bind missingness indicator columns and reorganize the dataframe.
    dfMiss <- BindAndOrganize(dfMiss, vars)
    # Make sure dfMiss and dfNorm (the original data) are both ordered alphabetically.
    dfMiss <- dfMiss[order(dfMiss$species_name), ]
    dfNorm <- dfNorm[order(dfNorm$species_name), ]
    
    # Imputation. ---
    # Ensure categorical traits are character type for mode imputation, including binary variables because we will be finding most frequent category (mode) value later on.
    # If there is more than one categorical trait..
    if(length(catTraits) > 1){
      # Use lapply to convert categorical traits to character class.
      dfNorm[, catTraits] <- lapply(dfNorm[, catTraits], as.character)
      dfMiss[, catTraits] <- lapply(dfMiss[, catTraits], as.character)
      # If there is only one categorical trait..
    } else if(length(catTraits) == 1) {
      # Convert single trait to character class.
      dfNorm[, catTraits] <- as.character(dfNorm[, catTraits])
      dfMiss[, catTraits] <- as.character(dfMiss[, catTraits])
    }
    # The ImputeMeanMode() function imputes values using mean/mode replacement for continuous and categorical variables respectively, and returns a list of error rates (MSE for continuous traits and PFC for categorical traits).
    impResult <- ImputeMeanMode(dfTrue = dfNorm, dfMissing = dfMiss, cols = vars, cont = contTraits, cat = catTraits)
    
    # Result handling. ---
    # Extract dfImp.
    dfImp <- impResult[[1]]
    # Extract errorRates.
    errorRates <- impResult[[2]]
    
    # Back-transforming data. ---
    # First, match data species (original species in complete-case) to species now in dfNorm in case any were removed (e.g. outgroups).
    dfOrig <- data[data$species_name %in% dfNorm$species_name, ]
    # Ensure dfOrig and dfImp are in same order. $$$
    dfOrig <- dfOrig[order(dfOrig$species_name), ]
    dfImp <- dfImp[order(dfImp$species_name), ]
    # Apply the BackTransform function to the continuous traits. $$$
    dfImp <- BackTransform(origData = dfOrig, tfData = dfImp, cols = contTraits)
    
    # Append to ith elements of lists.
    l_dfMiss[[i]] <- dfMissOrig
    l_l_missingness[[i]] <- l_missingness
    l_dfImp[[i]] <- dfImp
    l_Error[[i]] <- errorRates
    
  }
  
  # Final trait check. ---
  # If there are any traits that could not be simulated MAR..
  if(length(traitsNOSIM) > 0){
    # If the trait was identified in a previous iteration to be an trait that could not be simulated MAR, it is possible missing values were introduced in other iterations. These will likely be very low numbers of NAs and not enough for error rate analyses. So here we will subset l_l_Error to ensure it only contains values for traits that were simulated MAR.
    l_Error <- lapply(l_Error, function(x) {
      # Identify error rates associated with traitsMAR.
      index <- which(names(x) %in% traitsMAR)
      # Remove traits that were identified as MCAR in later iterations.
      x <- x[index]
    })
    print("The following traits could not be simulated MAR:")
    print(unique(traitsNOSIM))
  }
  
  # Average out the missingness for each trait. ---
  # Unlist missingness proportions.
  missingness <- unlist(l_l_missingness)
  # Calculate the average missingness proportion for each trait.
  avgMiss <- lapply(traitsMAR, function(x) {
    # Identify missingness proportions associated with the trait.
    index <- grep(pattern = x, names(missingness))
    # Take the mean.
    average <- mean(missingness[index], na.rm = T)
    # Name average according to the trait.
    names(average) <- x
    # Return the average missingness proportion.
    return(average)
  })
  
  # Subset to sample sizes for traits that could be simulated MAR.
  l_MARn <- l_sampleSizes[names(l_sampleSizes) %in% traitsMAR]
  # Order l_MARn by avgMiss.
  l_MARn <- l_MARn[names(unlist(avgMiss))]
  # Calculate the average number of NAs introduced for each trait.
  avgNAs <- mapply(function(x, y) x * y, y = l_MARn, x = avgMiss)
  # Rename according to l_MARn.
  names(avgNAs) <- names(l_MARn)
  
  # Create list to hold the results.
  l_results <- list(completeCaseData = data, rawData = raw, traits = l_traits, numeric = contTraits, categorical = catTraits, integer = intTraits, sampleSizes = l_sampleSizes, traitsNotSimulated = traitsNOSIM, traitsMAR = traitsMAR, MARfinalModels = l_MARfinalModels, reps = int, missingData = l_dfMiss, missingness = l_l_missingness, averageMissingness = avgMiss, averageNAs = avgNAs, imputedData = l_dfImp, errorRates = l_Error)
  # Return l_results.
  return(l_results)
  
}

KNNSimpute <- function(data, raw, vars, int = 100, k = 50, phyImp = F, tree = NULL) {
  
  # Given a complete-case dataset, this function simulates missingness at random (MAR) and imputes values using the kNN() function in the "VIM" package for imputation. Determines error rates for numerical (MSE) and categorical variables (ARI). +++
  # Citations: Kowarik A, Templ M (2016). “Imputation with the R Package VIM.” Journal of Statistical Software, 74(7), 1–16. doi: 10.18637/jss.v074.i07.
  # https://cran.r-project.org/web/packages/VIM/VIM.pdf
  
  # data = complete-case dataset containing trait data. Must also contain a column containing species name data = "species_name".
  # raw* = original dataset containing trait data with missing values. Used to build logistic regression models and simulate data that are MAR in the complete-case dataset. Must also contain a column with species name information = "species_name"
  # vars = names of columns containing trait data
  # int = number of iterations (missingness replicates)
  # k = maximum number of nearest neighbours to test
  # phyImp = whether to include phylogenetic information in the imputation process
  # tree = phylo object to be decomposed into phylogenetic eigenvectors if phyImp = T
  
  # Data matching. ---
  # If phylogenetic imputation was chosen..
  if(phyImp == T) {
    # Make sure the trait data and tree tips match.
    l_matched <- DropAndMatch(tree, data)
    # Extract updated tree.
    tree <- l_matched[[1]]
    # Extract updated dataframe.
    data <- l_matched[[2]]
  }
  
  # Trait preparation. ---
  # Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
  l_traits <- BreakIntoTypes(data, vars)
  # Extract numerical traits.
  contTraits <- l_traits[[1]]
  # Extract categorical traits.
  catTraits <- l_traits[[2]]
  # Ensure categorical traits are factor type for logistic regression model building (data class required for glm function).
  # If there is more than one categorical trait..
  if(length(catTraits) > 1){
    # Use lapply to convert categorical traits to factor class.
    data[, catTraits] <- lapply(data[, catTraits], as.factor)
    raw[, catTraits] <- lapply(raw[, catTraits], as.factor)
    # If there is only one categorical trait..
  } else if(length(catTraits) == 1) {
    # Convert single trait to factor class.
    data[, catTraits] <- as.factor(data[, catTraits])
    raw[, catTraits] <- as.factor(raw[, catTraits])
  }
  # Determine original sample size for each trait. &&
  l_sampleSizes <- lapply(data[, vars], function(x) length(na.omit(x)))
  
  # Predictor selection. ---
  # Here, we select those traits that have significant correlations to use as predictors for imputation (this will vary for each trait). Apply the SelectPredictors() function to obtain a list of predictors for each trait.
  l_predictors <- SelectPredictors(data[, vars])
  
  # Logistic regression fitting. ---
  # Apply FitLogReg() function to raw data to identify significant predictors of missingness for each trait.
  l_models <- FitLogReg(data = raw, cols = vars)
  # Now, let's identify traits that CAN be simulated MAR versus those that cannot. Traits that cannot be simulated MAR either 1) have no missing values in the original data or 2) do not have any predictors of missingness in the dataset upon fitting of the logistic regression models. IDMissPattern() also refines the logistic regression models (i.e. drops insignificant terms) through use of the RefineModels() function. +++
  l_missPatt <- IDMissPattern(data = raw, vars = vars, models = l_models)
  # Extract names of traits that could not be simulated MAR. +++
  traitsNOSIM <- l_missPatt[[1]]
  # Extract MAR traits.
  traitsMAR <- l_missPatt[[2]]
  # If traitsMAR is empty (no traits can be simulated MAR)..
  if(length(traitsMAR) == 0){
    stop("No traits can be simulated MAR!") 
  }
  # Extract final MAR models.
  l_MARfinalModels <- l_missPatt[[3]]
  
  # Imputation prep. ---
  # Create lists to hold the error rates for each iteration.
  l_l_Error <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create lists to hold the correlation coefficients for each iteration. +++
  l_l_Corr <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the dataframes with simulated missingness from each iteration.
  l_dfMissOrig <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the missingness proportion in each dfMiss.
  l_l_missingness <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the imputed dataframes from each iteration.
  l_l_dfImp <- CreateNamedList(listLength = int, elementNames = 1:int)
  # If phylogenetic imputation was chosen..
  if(phyImp == T) {
    # Create list to hold the dataframes and predictors with appended eigenvectors.
    l_l_evs <- CreateNamedList(listLength = int, elementNames = 1:int)
  }
  
  # For every iteration...
  for(i in 1:int) {
    
    # set.seed.
    set.seed(i)
    # MAR simulation. ---
    # Create copy of complete-case dataframe (untransformed) to introduce NAs into.
    dfMissOrig <- data
    # Create a list to hold the missingness proportion for each trait to be simulated MAR.
    l_missingness <- CreateNamedList(listLength = length(traitsMAR), elementNames = traitsMAR)
    
    # For every trait that can be simulated MAR..
    for(t in 1:length(traitsMAR)){
      # Take tth trait.
      trait <- traitsMAR[[t]]
      # Get the corresponding MAR model.
      index <- grep(pattern = trait, x = names(l_MARfinalModels))
      modelMAR <- l_MARfinalModels[[index]]
      # Simulate MAR data in complete-case dataset.
      res <- SimMAR(model = modelMAR, data = data)
      # Get original sample size of trait. ++
      n <- l_sampleSizes[[grep(trait, names(l_sampleSizes))]]
      # Subtract data originally missing in data from the number missing in dfMissOrig after MAR simulation to determine actual number of NAs introduced. ++
      marN <- sum(is.na(res)) - sum(is.na(data[[trait]]))
      # If missingness is less than 0.80 (enough data for imputation)..++
      if(marN/n < 0.80){
        # Replace complete-case column in dfMissOrig with res. ++
        dfMissOrig[, trait] <- res 
        # Divide by n to determine missingness percentage for trait and append to l_missingness. ++
        l_missingness[[t]] <- marN/n
      } else {
        # Assign NA as missingness proportion value.
        l_missingness[[t]] <- NA
      }
    }
    
    # Append traits that could not be simulated MAR (no NAs introduced based on fitted logistic regression models) OR traits that exceed 0.80 missingness to traitsNOSIM. +++
    traitsNOSIM <- c(traitsNOSIM, names(which(l_missingness == 0)), names(which(is.na(l_missingness))))
    # Update traitsMAR based on this result.  
    traitsMAR <- setdiff(traitsMAR, traitsNOSIM)
    # If traitsMAR is empty (no traits can be simulated MAR)..
    if(length(traitsMAR) == 0){
      stop("No traits can be simulated MAR!") 
    }
    
    # Dataframe organization. ---
    # Make copy of data and rename. $$$
    dfNorm <- data
    # Normalize numerical traits prior to imputation.
    # If there is more than one numerical trait..
    if(length(contTraits) > 1){
      # Use lapply to scale numerical traits.
      dfNorm[, contTraits] <- lapply(dfNorm[, contTraits], scale)
      dfNorm[, contTraits] <- lapply(dfNorm[, contTraits], as.numeric)
      # If there is only one numerical trait..
    } else if(length(contTraits) == 1) {
      # Scale single trait. ==
      dfNorm[, contTraits] <- as.numeric(scale(dfNorm[, contTraits]))
    }
    # Introduce missing values from dfMissOrig into dfNorm.
    dfMiss <- as.data.frame(mapply(function(x, y) ifelse(is.na(x), x, y), x = dfMissOrig[, c("species_name", contTraits)], y = dfNorm[, c("species_name", contTraits)], SIMPLIFY = F))
    # Add categorical traits back.
    dfMiss <- merge(dfMiss, dfMissOrig[, c("species_name", catTraits)], by = "species_name")
    # Bind missingness indicator columns and reorganize the dataframe.
    dfMiss <- BindAndOrganize(dfMiss, vars)
    # Make sure dfMiss and dfNorm (the original data) are both ordered alphabetically.
    dfMiss <- dfMiss[order(dfMiss$species_name), ]
    dfNorm <- dfNorm[order(dfNorm$species_name), ]
    # Ensure categorical traits are factor type.
    # If there is more than one categorical trait..
    if(length(catTraits) > 1){
      # Use lapply to convert categorical traits to factor class.
      dfNorm[, catTraits] <- lapply(dfNorm[, catTraits], as.factor)
      dfMiss[, catTraits] <- lapply(dfMiss[, catTraits], as.factor)
      # If there is only one categorical trait..
    } else if(length(catTraits) == 1) {
      # Convert single trait to factor class.
      dfNorm[, catTraits] <- as.factor(dfNorm[, catTraits])
      dfMiss[, catTraits] <- as.factor(dfMiss[, catTraits])
    }
    
    
    # Phylogenetic eigenvector decomposition. ---
    if(phyImp == T) {
      # Append eigenvectors to dfMiss and list of predictors. Each trait will have a corresponding dataframe and list of predictors including the eigenvectors.
      l_evs <- AppendEigenvectors(data = dfMiss, vars = traitsMAR, tree = tree, predictors = l_predictors)
      # Extract list of dfMiss.
      l_dfMissEV <- l_evs[[1]]
      # Extract updated list of predictors.
      l_EVPredictors <- l_evs[[2]]
    }
    
    # Identify outgroup(s) in trait datasets (this is the species that contains no trait data dfNorm). &&
    outgroup <- dfNorm$species_name[apply(dfNorm[, vars], 1, function(x) all(is.na(x)))]
    # If found in dataframe.. &&
    if(length(outgroup) > 0){
      # Remove outgroup from dataframes.
      dfNorm <- dfNorm[!dfNorm$species_name %in% outgroup, ]
      dfMiss <- dfMiss[!dfMiss$species_name %in% outgroup, ]
      if(phyImp == T){
        # Also remove outgroup from l_dfMissEV.
        l_dfMissEV <- lapply(l_dfMissEV, function(x) x[!x$species_name %in% outgroup, ])
      }
    }
    
    # Imputation. ---   
    # The ImputeKNN() function entails a loop that uses different values of k (the number of nearest neighbours to use in the kNN algorithm) and returns a list of error rates (MSE for continuous traits and ARI for categorical traits) for each parameter value tested. +++
    # If phylogenetic imputation was selected..
    if(phyImp == T) {
      # Impute data using eigenvectors.
      impResult <- ImputeKNN(dfTrue = dfNorm, dfMissing = dfMiss, cols = traitsMAR, cont = contTraits, cat = catTraits, k = 50, predictors = l_EVPredictors, phyImp = T, l_dfMissing = l_dfMissEV)
    } else if(phyImp == F){
      # Impute data only using trait data.
      impResult <- ImputeKNN(dfTrue = dfNorm, dfMissing = dfMiss, cols = traitsMAR, cont = contTraits, cat = catTraits, k = 50, predictors = l_predictors)
    }
    
    # Result handling. ---  
    # Extract l_dfImp.  
    l_dfImp <- impResult[[1]]
    # Extract errorRates  .
    l_Error <- impResult[[2]]
    # Extract corrCoefs. +++
    l_Corr <- impResult[[3]]
    
    # Back-transforming data. ---
    # First, match data species (original species in complete-case) to species now in dfNorm in case any were removed (e.g. outgroups).
    dfOrig <- data[data$species_name %in% dfNorm$species_name, ]
    # Ensure dfOrig and dfImp are in same order. $$$
    dfOrig <- dfOrig[order(dfOrig$species_name), ]
    l_dfImp <- lapply(l_dfImp, function(x) x[order(x$species_name), ])
    
    # For every imputed dataframe..
    for(d in 1:length(l_dfImp)){
      # Take the dth dfImp.
      dfImp <- l_dfImp[[d]]
      # Apply BackTransform function to contTraits in dfImp. $$$
      dfImp <- BackTransform(origData = dfOrig, tfData = dfImp, cols = traitsMAR[traitsMAR %in% contTraits])
      # Replace dfImp in l_dfImp with newly backtransformed dataset.
      l_dfImp[[d]] <- dfImp
    }
    
    # Append results to ith element of lists.
    l_dfMissOrig[[i]] <- dfMissOrig
    l_l_missingness[[i]] <- l_missingness
    l_l_dfImp[[i]] <- l_dfImp
    l_l_Error[[i]] <- l_Error
    l_l_Corr[[i]] <- l_Corr # +++
    # If phylogenetic imputation was selected..
    if(phyImp == T) {
      l_l_evs[[i]] <- l_evs
    }
    
  } ## i
  
  # Final trait check. --- +++
  # If there are any traits that could not be simulated MAR..
  if(length(traitsNOSIM) > 0){
    # If the trait was identified in a previous iteration to be a trait that couldn't be simulated MAR, it is possible missing values were introduced in other iterations. These will likely be very low numbers of NAs and not enough for error rate analyses. So here we will subset l_l_Error to ensure it only contains values for varsMAR.
    l_l_Error <- lapply(l_l_Error, function(x) {
      # Identify error rates associated with traitsMAR.
      index <- which(names(x) %in% traitsMAR)
      # Remove traits that were identified as NOSIM in later iterations.
      x <- x[index]
    })
    
    l_l_Corr <- lapply(l_l_Corr, function(x) {
      # Identify error rates associated with traitsMAR.
      index <- which(names(x) %in% traitsMAR)
      # Remove traits that were identified as NOSIM in later iterations.
      x <- x[index]
    })
    
    print("The following traits could not be simulated MAR:")
    print(unique(traitsNOSIM))
  }
  
  # Average out the missingness for each trait. ---
  # Unlist missingness proportions.
  missingness <- unlist(l_l_missingness)
  # Calculate the average missingness proportion for each trait.
  avgMiss <- lapply(traitsMAR, function(x) {
    # Identify missingness proportions associated with the trait.
    index <- grep(pattern = x, names(missingness))
    # Take the mean.
    average <- mean(missingness[index], na.rm = T)
    # Name average according to the trait.
    names(average) <- x
    # Return the average missingness proportion.
    return(average)
  })
  # Subset to sample sizes for traits that could be simulated MAR.
  l_MARn <- l_sampleSizes[names(l_sampleSizes) %in% traitsMAR]
  # Order l_MARn by avgMiss. +++
  l_MARn <- l_MARn[names(unlist(avgMiss))]
  # Calculate the average number of NAs introduced for each trait.
  avgNAs <- mapply(function(x, y) x * y, y = l_MARn, x = avgMiss)
  # Rename according to l_MARn.
  names(avgNAs) <- names(l_MARn)
  
  # If phylogenetic imputation was selected..
  if(phyImp == T) { 
    # Create list to hold the results. +++
    l_results <- list(completeCaseData = data, rawData = raw, traits = l_traits, numeric = contTraits, categorical = catTraits, sampleSizes = l_sampleSizes, traitsNotSimulated = traitsNOSIM, traitsSimulated = traitsMAR, MARfinalModels = l_MARfinalModels, k = k, reps = int, predictors = l_EVPredictors, tree = tree, eigenvectors = l_l_evs, missingData = l_dfMissOrig, missingness = l_l_missingness, averageMissingness = avgMiss, averageNAs = avgNAs, imputedData = l_l_dfImp, errorRates = l_l_Error, corrCoef = l_l_Corr)
  } else if(phyImp == F){
    # Create list to hold the results. +++
    l_results <- list(completeCaseData = data, rawData = raw, traits = l_traits, numeric = contTraits, categorical = catTraits, sampleSizes = l_sampleSizes, traitsNotSimulated = traitsNOSIM, traitsSimulated = traitsMAR, MARfinalModels = l_MARfinalModels, k = k, reps = int, predictors = l_predictors, missingData = l_dfMissOrig, missingness = l_l_missingness, averageMissingness = avgMiss, averageNAs = avgNAs, imputedData = l_l_dfImp, errorRates = l_l_Error, corrCoef = l_l_Corr)
  }
  # Return results.  
  return(l_results)
  
}

MICESimpute <- function(data, raw, vars, int = 100, mSets = c(5, 10, 40), phyImp = F, tree = NULL) {
  
  # Given a complete-case dataset, this function simulates missingness at random (MAR) and imputes values using the mice() function in the "mice" package. Determines error rates for numerical (MSE) and categorical variables (ARI). +++
  # Citations: van Buuren S, Groothuis-Oudshoorn K (2011). “mice: Multivariate Imputation by Chained Equations in R.” Journal of Statistical Software, 45(3), 1-67. https://www.jstatsoft.org/v45/i03/.
  # R package version 3.13.0. https://cran.r-project.org/web/packages/mice/mice.pdf
  # Vignettes consulted: Gerko Vink and Stef van Buuren. miceVignettes. https://www.gerkovink.com/miceVignettes/
  # Rianne Schouten and Gerko Vink. Wrapper function parlMICE. https://www.gerkovink.com/parlMICE/Vignette_parlMICE.html
  
  # data = complete-case dataset containing trait data. Must also contain a column containing species name called "species_name".
  # raw = original dataset containing trait data with missing values. Used to build logistic regression models and simulate data that are MAR in the complete-case dataset. Must also contain a column with species name information = "species_name"  
  # vars = names of columns containing trait data
  # int = number of iterations (missingness replicates)
  # mSets = vector containing values of m to test (number of multiply imputed dataframes)
  # phyImp = whether to include phylogenetic information in the imputation process
  # tree = phylo object to be decomposed into phylogenetic eigenvectors if phyImp = T
  
  # Data matching. ---
  # If phylogenetic imputation was chosen..
  if(phyImp == T) {
    # Make sure the trait data and tree tips match.
    l_matched <- DropAndMatch(tree, data)
    # Extract updated tree.
    tree <- l_matched[[1]]
    # Extract updated dataframe.
    data <- l_matched[[2]]
  }
  
  # Trait preparation. ---
  # Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
  l_traits <- BreakIntoTypes(data, vars)
  # Extract numerical traits.
  contTraits <- l_traits[[1]]
  # Extract categorical traits.
  catTraits <- l_traits[[2]]
  # Ensure categorical traits are factor type for logistic regression model building (data class required for glm function).
  # If there is more than one categorical trait..
  if(length(catTraits) > 1){
    # Use lapply to convert categorical traits to factor class.
    data[, catTraits] <- lapply(data[, catTraits], as.factor)
    raw[, catTraits] <- lapply(raw[, catTraits], as.factor)
    # If there is only one categorical trait..
  } else if(length(catTraits) == 1) {
    # Convert single trait to factor class.
    data[, catTraits] <- as.factor(data[, catTraits])
    raw[, catTraits] <- as.factor(raw[, catTraits])
  }
  # Identify integer (count) traits, if any.   
  intTraits <- GetTraitNames(data = data[, vars], class = "integer")
  # Determine original sample size for each trait. &&
  l_sampleSizes <- lapply(data[, vars], function(x) length(na.omit(x)))
  
  # Predictor selection. ---
  # Here, we select those traits that have significant correlations to use as predictors for imputation (this will vary for each trait). Apply the SelectPredictors() function to obtain a list of predictors for each trait.
  l_predictors <- SelectPredictors(data[, vars])
  
  # Logistic regression fitting. ---  
  # Apply FitLogReg() function to raw data to identify significant predictors of missingness for each trait.
  l_models <- FitLogReg(data = raw, cols = vars)
  # Now, let's identify traits that CAN be simulated MAR versus those that cannot. Traits that cannot be simulated MAR either 1) have no missing values in the original data or 2) do not have any predictors of missingness in the dataset upon fitting of the logistic regression models. IDMissPattern() also refines the logistic regression models (i.e. drops insignificant terms) through use of the RefineModels() function. +++
  l_missPatt <- IDMissPattern(data = raw, vars = vars, models = l_models)
  # Extract names of traits that could not be simulated MAR. +++
  traitsNOSIM <- l_missPatt[[1]]
  # Extract MAR traits.
  traitsMAR <- l_missPatt[[2]]
  # If traitsMAR is empty (no traits can be simulated MAR)..
  if(length(traitsMAR) == 0){
    stop("No traits can be simulated MAR!") 
  }
  # Extract final MAR models.
  l_MARfinalModels <- l_missPatt[[3]]
  
  # Imputation prep. ---
  # Create lists to hold the error rates for each iteration.  
  l_l_Error <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the dataframes with simulated missingness from each iteration.    
  # Create lists to hold the correlation coefficients for each iteration. +++
  l_l_Corr <- CreateNamedList(listLength = int, elementNames = 1:int)
  l_dfMissOrig <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the missingness proportion in each dfMiss.    
  l_l_missingness <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the imputed dataframes from each iteration.    
  l_l_dfImp <- CreateNamedList(listLength = int, elementNames = 1:int)
  # If phylogenetic imputation was chosen..
  if(phyImp == T) {
    # Create list to hold the dataframes and predictors with appended eigenvectors.    
    l_l_evs <- CreateNamedList(listLength = int, elementNames = 1:int)
    l_dfMissEV <- CreateNamedList(listLength = int, elementNames = 1:int)
  }
  
  # For every iteration...
  for(i in 1:int) {
    
    # set.seed.
    set.seed(i)
    # MAR simulation. ---  
    # Create copy of complete-case dataframe (untransformed) to introduce NAs into.
    dfMissOrig <- data
    # Create a list to hold the missingness proportion for each trait to be simulated MAR.  
    l_missingness <- CreateNamedList(listLength = length(traitsMAR), elementNames = traitsMAR)
    # For every trait that can be simulated MAR..
    for(t in 1:length(traitsMAR)){
      # Take tth trait.  
      trait <- traitsMAR[[t]]
      # Get the corresponding MAR model.
      index <- grep(pattern = trait, x = names(l_MARfinalModels))
      modelMAR <- l_MARfinalModels[[index]]
      # Simulate MAR data in complete-case dataset.
      res <- SimMAR(model = modelMAR, data = data)
      # Get original sample size of trait. ++
      n <- l_sampleSizes[[grep(trait, names(l_sampleSizes))]]
      # Subtract data originally missing in data from the number missing in dfMissOrig after MAR simulation to determine actual number of NAs introduced. ++
      marN <- sum(is.na(res)) - sum(is.na(data[[trait]]))
      # If missingness is less than 0.80 (enough data for imputation)..++
      if(marN/n < 0.80){
        # Replace complete-case column in dfMissOrig with res. ++
        dfMissOrig[, trait] <- res 
        # Divide by n to determine missingness percentage for trait and append to l_missingness. ++
        l_missingness[[t]] <- marN/n
      } else {
        # Assign NA as missingness proportion value.
        l_missingness[[t]] <- NA
      }
    }
    
    # Finally, append traits that could not be simulated MAR (no NAs introduced based on fitted logistic regression models) or that exceed 0.80 missingness (hard time getting accurate imputations) to traitsNOSIM.
    traitsNOSIM <- c(traitsNOSIM, names(which(l_missingness == 0)), names(which(is.na(l_missingness))))
    # Update traitsMAR based on this result.  
    traitsMAR <- setdiff(traitsMAR, traitsNOSIM)
    # If traitsMAR is empty (no traits can be simulated MAR)..
    if(length(traitsMAR) == 0){
      stop("No traits can be simulated MAR!") 
    }
    
    # Dataframe organization. ---
    # Make copy of data and rename. $$$
    dfNorm <- data
    # Normalize numerical traits prior to imputation. $$$
    # If there is more than one numerical trait..
    if(length(contTraits) > 1){
      # Use lapply to scale numerical traits.
      dfNorm[, contTraits] <- lapply(dfNorm[, contTraits], scale)
      dfNorm[, contTraits] <- lapply(dfNorm[, contTraits], as.numeric)
      # If there is only one numerical trait..
    } else if(length(contTraits) == 1) {
      # Scale single trait. ==
      dfNorm[, contTraits] <- as.numeric(scale(dfNorm[, contTraits]))
    }
    # Introduce missing values from dfMissOrig into dfNorm. !!!
    dfMiss <- as.data.frame(mapply(function(x, y) ifelse(is.na(x), x, y), x = dfMissOrig[, c("species_name", contTraits)], y = dfNorm[, c("species_name", contTraits)], SIMPLIFY = F))
    # Add categorical traits back. !!!
    dfMiss <- merge(dfMiss, dfMissOrig[, c("species_name", catTraits)], by = "species_name")
    # Bind missingness indicator columns and reorganize the dataframe.
    dfMiss <- BindAndOrganize(dfMiss, vars)
    # Make sure dfMiss and dfNorm (the original data) are both ordered alphabetically.
    dfMiss <- dfMiss[order(dfMiss$species_name), ]
    dfNorm <- dfNorm[order(dfNorm$species_name), ]
    # If there is more than one categorical trait..
    if(length(catTraits) > 1){
      # Use lapply to convert categorical traits to factor class.
      dfNorm[, catTraits] <- lapply(dfNorm[, catTraits], as.factor)
      dfMiss[, catTraits] <- lapply(dfMiss[, catTraits], as.factor)
      # If there is only one categorical trait..
    } else if(length(catTraits) == 1) {
      # Convert single trait to factor class.
      dfNorm[, catTraits] <- as.factor(dfNorm[, catTraits])
      dfMiss[, catTraits] <- as.factor(dfMiss[, catTraits])
    }
    
    # Phylogenetic eigenvector decomposition. ---
    if(phyImp == T) {
      # Append eigenvectors to dfMiss and list of predictors. Each trait will have a corresponding dataframe and list of predictors including the eigenvectors.
      l_evs <- AppendEigenvectors(data = dfMiss, vars = vars, tree = tree, predictors = l_predictors)
      # Extract list of dfMiss.
      l_dfMissEV <- l_evs[[1]]
      # Extract updated list of predictors.
      l_EVPredictors <- l_evs[[2]]
      # Dataframe list collapse. ---   
      # MICE requires a predictor matrix for imputation. So, we will first merge all of the dataframe in l_dfMissEV into one dataframe. This is so we can identify the eigenvector columns needed for imputation of each trait and update this in the predictor matrix. This should also same computation time since we are just imputing one dataframe and not a list of dataframes. 
      dfMissEV <- l_dfMissEV[[1]]
      # Identify the eigenvector columns and remove those dfMiss.
      evIndex <- grep(pattern = "V_", colnames(dfMissEV))
      dfMissEV <- dfMissEV[, -evIndex]
      # Also remove the missingness indicator columns.
      missIndex <- grep(pattern = "_NA", colnames(dfMissEV))
      dfMissEV <- dfMissEV[, -missIndex]
      # Now, let's merge dfMissEV with the corresponding eigenvector columns in each dataframe in l_dfMissEV.
      for(e in 1:length(l_dfMissEV)){
        # Take the eth l_dfMissEV.
        dfMisseth <- l_dfMissEV[[e]]
        # Identify eigenvector columns in eth dataframe.
        index <- grep("V_", colnames(dfMisseth))
        evCols <- colnames(dfMisseth)[index]
        # Merge these columns with dfMisseth (including species_name).
        dfMissEV <- merge(dfMissEV, dfMisseth[, c("species_name", evCols)], by = "species_name")
      }
    }
    # Identify outgroup(s) in trait datasets (this is the species that contains no trait data dfNorm). &&
    outgroup <- dfNorm$species_name[apply(dfNorm[, vars], 1, function(x) all(is.na(x)))]
    # If found in dataframe.. &&
    if(length(outgroup) > 0){
      # Remove outgroup from dataframes.
      dfNorm <- dfNorm[!dfNorm$species_name %in% outgroup, ]
      dfMiss <- dfMiss[!dfMiss$species_name %in% outgroup, ]
      # If phylogenetic imputation was selected.. &&   
      if(phyImp == T) {
        l_dfMissEV <- lapply(l_dfMissEV, function(x) x[!x$species_name %in% outgroup, ])
        dfMissEV <- dfMissEV[!dfMissEV$species_name %in% outgroup, ]
      }
    }
    
    # Predictor matrix creation. ---   
    # If phylogenetic imputation wasn't chosen..
    if(phyImp == F){
      # Create predictor matrix for use in MICE imputation (based on results of our trait predictor screening).
      matPredictors <- CreatePredictorMatrix(cols = vars, predictors = l_predictors, dfMissing = dfMiss)
      # If phylogenetic imputation was chosen..
    } else if(phyImp == T){
      # Create predictor matrix for use in MICE imputation with appended eigenvectors. We can indicate the column names of dfMissEV (excluding species_name) as names of variables to consider in imputation process.
      matPredictors <- CreatePredictorMatrix(dfMissing = dfMissEV, cols = colnames(dfMissEV)[-1], predictors = l_EVPredictors)
      # We can set all of the rows for the eigenvectors to 0 as they themselves will not be imputed.
      eigens <- grep(pattern = "V_", rownames(matPredictors))
      matPredictors[eigens, ] <- 0
    }
    
    # Imputation. ---   
    # The ImputeMICE() function entails a loop that uses different values of m (the number of multiply imputed datasets) and returns a list of error rates (MSE for continuous traits and ARI for categorical traits) for each parameter value tested.
    # If phylogenetic imputation wasn't selected..  
    if(phyImp == F) {
      # Impute data only using trait data. +
      impResult <- ImputeMICE(dfTrue = dfNorm, dfMissing = dfMiss, dfOrig = data, cols = traitsMAR, cont = contTraits, cat = catTraits, inter = intTraits, mSets = mSets, matPredictors = matPredictors, seed = i)
      # If phylogenetic imputation was selected..
    } else if(phyImp == T){
      # Bind back the missing columns.   
      dfMissEV <- bind_shadow(dfMissEV, only_miss = T)
      # Impute data using eigenvectors. +
      impResult <- ImputeMICE(dfTrue = dfNorm, dfMissing = dfMissEV, dfOrig = data, cols = traitsMAR, cont = contTraits, cat = catTraits, inter = intTraits, mSets = mSets, matPredictors = matPredictors, seed = i)
    }
    # Result handling. ---  
    # Extract l_dfImp.  
    l_dfImp <- impResult[[1]]
    # Extract errorRates  .
    l_Error <- impResult[[2]]
    # Extract corrCoefs. +++
    l_Corr <- impResult[[3]]
    
    # Back-transforming data. ---  
    # First, match data species (original species in complete-case) to species now in dfNorm in case any were removed (e.g. outgroups).
    dfOrig <- data[data$species_name %in% dfNorm$species_name, ]
    # Ensure dfOrig and dfImp are in same order. $$$
    dfOrig <- dfOrig[order(dfOrig$species_name), ]
    l_dfImp <- lapply(l_dfImp, function(x) x[order(x$species_name), ])
    # For every imputed dataframe..
    for(d in 1:length(l_dfImp)){
      # Take the dth dfImp.
      dfImp <- l_dfImp[[d]]
      # Apply BackTransform function to contTraits in dfImp. $$$
      dfImp <- BackTransform(origData = dfOrig, tfData = dfImp, cols = traitsMAR[traitsMAR %in% contTraits])
      # If there are any integer traits..
      if(length(intTraits) > 0){
        # Round to nearest whole number.
        dfImp[, intTraits] <- lapply(dfImp[, intTraits], function(x) as.integer(round(x)))
      }
      # Replace dfImp in l_dfImp with newly backtransformed dataset.
      l_dfImp[[d]] <- dfImp
    }
    
    # Append results to ith element of lists.    
    l_dfMissOrig[[i]] <- dfMissOrig
    l_l_missingness[[i]] <- l_missingness      
    l_l_dfImp[[i]] <- l_dfImp
    l_l_Error[[i]] <- l_Error
    l_l_Corr[[i]] <- l_Corr # +++
    # If phylogenetic imputation was selected..
    if(phyImp == T) {
      l_l_evs[[i]] <- l_evs
      l_dfMissEV[[i]] <- dfMissEV
    }
  } ## i
  
  # Final trait check. ---  
  # If there are any traits that could not be simulated MAR..
  if(length(traitsNOSIM) > 0){
    # If the trait was identified in a previous iteration to be a trait that couldn't be simulated MAR, it is possible missing values were introduced in other iterations. These will likely be very low numbers of NAs and not enough for error rate analyses. So here we will subset l_l_Error to ensure it only contains values for varsMAR.
    l_l_Error <- lapply(l_l_Error, function(x) {
      lapply(x, function(e){
        # Identify error rates associated with traitsMAR.
        index <- which(names(e) %in% traitsMAR)
        # Remove traits that were identified as NOSIM in later iterations.
        e <- e[index]
      })
    })
    print("The following traits could not be simulated MAR:")
    print(unique(traitsNOSIM))
    
    l_l_Corr <- lapply(l_l_Corr, function(x) {
      lapply(x, function(e){
        # Identify error rates associated with traitsMAR.
        index <- which(names(e) %in% traitsMAR)
        # Remove traits that were identified as NOSIM in later iterations.
        e <- e[index]
      })
    })
    
  }
  
  # Average out the missingness for each trait. ---    
  # Unlist missingness proportions.
  missingness <- unlist(l_l_missingness)
  # Calculate the average missingness proportion for each trait.
  avgMiss <- lapply(traitsMAR, function(x) {
    # Identify missingness proportions associated with the trait.
    index <- grep(pattern = x, names(missingness))
    # Take the mean.
    average <- mean(missingness[index], na.rm = T)
    # Name average according to the trait.
    names(average) <- x
    # Return the average missingness proportion.
    return(average)
  })
  # Subset to sample sizes for traits that could be simulated MAR.
  l_MARn <- l_sampleSizes[names(l_sampleSizes) %in% traitsMAR]
  # Order l_MARn by avgMiss.
  l_MARn <- l_MARn[names(unlist(avgMiss))]
  # Calculate the average number of NAs introduced for each trait.
  avgNAs <- mapply(function(x, y) x * y, y = l_MARn, x = avgMiss)
  # Rename according to l_MARn.
  names(avgNAs) <- names(l_MARn)
  
  # If phylogenetic imputation was selected..
  if(phyImp == T) { 
    # Create list to hold the results. +++
    l_results <- list(completeCaseData = data, rawData = raw, traits = l_traits, numeric = contTraits, categorical = catTraits, integer = intTraits, traitsNotSimulated = traitsNOSIM, traitsSimulated = traitsMAR, MARfinalModels = l_MARfinalModels, m = mSets, predictorMatrix = matPredictors, reps = int, predictors = l_EVPredictors, tree = tree, eigenvectors = l_l_evs, missingData = l_dfMissOrig, eigenTraitData = l_dfMissEV, missingness = l_l_missingness, averageMissingness = avgMiss, averageNAs = avgNAs, imputedData = l_l_dfImp, errorRates = l_l_Error, corrCoef = l_l_Corr)
  } else if(phyImp == F){
    # Create list to hold the results. +++
    l_results <- list(completeCaseData = data, rawData = raw, traits = l_traits, numeric = contTraits, categorical = catTraits, integer = intTraits, traitsNotSimulated = traitsNOSIM, traitsSimulated = traitsMAR, MARfinalModels = l_MARfinalModels, m = mSets, predictorMatrix = matPredictors, reps = int, predictors = l_predictors, missingData = l_dfMissOrig, missingness = l_l_missingness, averageMissingness = avgMiss, averageNAs = avgNAs, imputedData = l_l_dfImp, errorRates = l_l_Error, corrCoef = l_l_Corr)
  }
  # Return l_results.  
  return(l_results)
  
}

RFSimpute <- function(data, raw, vars, int = 100, ntrees = c(100, 1000), phyImp = F, tree = NULL) {
  
  # Given a complete-case dataset, this function simulates missingness at random (MAR) and imputes values using the missForest() function in the "missForest" package. Determines error rates for numerical (MSE) and categorical variables (ARI).
  # Citations: Daniel J. Stekhoven (2013). missForest: Nonparametric Missing Value Imputation using Random Forest. R package version 1.4.
  # Stekhoven D. J., & Buehlmann, P. (2012). MissForest - non-parametric missing value imputation for mixed-type data. Bioinformatics, 28(1), 112-118.
  
  # data = complete-case dataset containing trait data. Must also contain a column containing species name called "species_name".
  # raw* = original dataset containing trait data with missing values. Used to build logistic regression models and simulate data that are MAR in the complete-case dataset. Must also contain a column with species name information = "species_name"
  # vars = names of columns containing trait data
  # int = number of iterations (missingness replicates)
  # ntrees = vector containing values of ntree to test
  # phyImp = whether to include phylogenetic information in the imputation process
  # tree = phylo object to be decomposed into phylogenetic eigenvectors if phyImp = T
  
  # Data matching. ---
  # If phylogenetic imputation was chosen..
  if(phyImp == T) {
    # Make sure the trait data and tree tips match.
    l_matched <- DropAndMatch(tree, data)
    # Extract updated tree.
    tree <- l_matched[[1]]
    # Extract updated dataframe.
    data <- l_matched[[2]]
  }
  
  # Trait preparation. ---
  # Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
  l_traits <- BreakIntoTypes(data, vars)
  # Extract numerical traits.
  contTraits <- l_traits[[1]]
  # Extract categorical traits.
  catTraits <- l_traits[[2]]
  # Ensure categorical traits are factor type for logistic regression model building (data class required for glm function).
  # If there is more than one categorical trait..
  if(length(catTraits) > 1){
    # Use lapply to convert categorical traits to factor class.
    data[, catTraits] <- lapply(data[, catTraits], as.factor)
    raw[, catTraits] <- lapply(raw[, catTraits], as.factor)
    # If there is only one categorical trait..
  } else if(length(catTraits) == 1) {
    # Convert single trait to factor class.
    data[, catTraits] <- as.factor(data[, catTraits])
    raw[, catTraits] <- as.factor(raw[, catTraits])
  }
  # Identify integer (count) traits, if any.
  intTraits <- GetTraitNames(data = data[, vars], class = "integer")
  # Determine original sample size for each trait. &&
  l_sampleSizes <- lapply(data[, vars], function(x) length(na.omit(x)))
  
  # Predictor selection. ---
  # Here, we select those traits that have significant correlations to use as predictors for imputation (this will vary for each trait). Apply the SelectPredictors() function to obtain a list of predictors for each trait.
  l_predictors <- SelectPredictors(data[, vars])
  
  # Logistic regression fitting. ---
  # Apply FitLogReg() function to raw data to identify significant predictors of missingness for each trait.
  l_models <- FitLogReg(data = raw, cols = vars)
  # Now, let's identify traits that CAN be simulated MAR versus those that cannot. Traits that cannot be simulated MAR either 1) have no missing values in the original data or 2) do not have any predictors of missingness in the dataset upon fitting of the logistic regression models. IDMissPattern() also refines the logistic regression models (i.e. drops insignificant terms) through use of the RefineModels() function. +++
  l_missPatt <- IDMissPattern(data = raw, vars = vars, models = l_models)
  # Extract names of traits that could not be simulated MAR. +++
  traitsNOSIM <- l_missPatt[[1]]
  # Extract MAR traits.
  traitsMAR <- l_missPatt[[2]]
  # If traitsMAR is empty (no traits can be simulated MAR)..
  if(length(traitsMAR) == 0){
    stop("No traits can be simulated MAR!") 
  }
  # Extract final MAR models.
  l_MARfinalModels <- l_missPatt[[3]]
  
  # Imputation prep. ---
  # Create lists to hold the error rates for each iteration.  
  l_l_Error <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create lists to hold the correlation coefficients for each iteration. +++
  l_l_Corr <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the dataframes with simulated missingness from each iteration.    
  l_dfMissOrig <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the missingness proportion in each dfMiss.    
  l_l_missingness <- CreateNamedList(listLength = int, elementNames = 1:int)
  # Create list to hold the imputed dataframes from each iteration.    
  l_l_dfImp <- CreateNamedList(listLength = int, elementNames = 1:int)
  # If phylogenetic imputation was chosen..
  if(phyImp == T) {
    # Create list to hold the dataframes and predictors with appended eigenvectors.    
    l_l_evs <- CreateNamedList(listLength = int, elementNames = 1:int)
  }
  
  # For every iteration...
  for(i in 1:int) {
    
    # set.seed.
    set.seed(i)
    # MAR simulation. ---  
    # Create copy of complete-case dataframe (untransformed) to introduce NAs into.
    dfMissOrig <- data 
    # Create a list to hold the missingness proportion for each trait to be simulated MAR.  
    l_missingness <- CreateNamedList(listLength = length(traitsMAR), elementNames = traitsMAR)
    # For every trait that can be simulated MAR..
    for(t in 1:length(traitsMAR)){
      # Take tth trait.  
      trait <- traitsMAR[[t]]
      # Get the corresponding MAR model.
      index <- grep(pattern = trait, x = names(l_MARfinalModels))
      modelMAR <- l_MARfinalModels[[index]]
      # Simulate MAR data in complete-case dataset.
      res <- SimMAR(model = modelMAR, data = data)
      # Get original sample size of trait. ++
      n <- l_sampleSizes[[grep(trait, names(l_sampleSizes))]]
      # Subtract data originally missing in data from the number missing in dfMissOrig after MAR simulation to determine actual number of NAs introduced. ++
      marN <- sum(is.na(res)) - sum(is.na(data[[trait]]))
      # If missingness is less than 0.80 (enough data for imputation)..++
      if(marN/n < 0.80){
        # Replace complete-case column in dfMissOrig with res. ++
        dfMissOrig[, trait] <- res 
        # Divide by n to determine missingness percentage for trait and append to l_missingness. ++
        l_missingness[[t]] <- marN/n
      } else {
        # Assign NA as missingness proportion value.
        l_missingness[[t]] <- NA
      }
    }
    
    # Finally, append traits that could not be simulated MAR (no NAs introduced based on fitted logistic regression models) or that exceed 0.80 missingness (hard time getting accurate imputations) to traitsNOSIM.
    traitsNOSIM <- c(traitsNOSIM, names(which(l_missingness == 0)), names(which(is.na(l_missingness))))
    # Update traitsMAR based on this result.  
    traitsMAR <- setdiff(traitsMAR, traitsNOSIM)
    # If traitsMAR is empty (no traits can be simulated MAR)..
    if(length(traitsMAR) == 0){
      stop("No traits can be simulated MAR!") 
    }
    
    # Dataframe organization. ---
    # Make copy of data and rename. $$$
    dfNorm <- data
    # Normalize numerical traits prior to imputation. $$$
    # If there is more than one numerical trait..
    if(length(contTraits) > 1){
      # Use lapply to scale numerical traits.
      dfNorm[, contTraits] <- lapply(dfNorm[, contTraits], scale)
      dfNorm[, contTraits] <- lapply(dfNorm[, contTraits], as.numeric)
      # If there is only one numerical trait..
    } else if(length(contTraits) == 1) {
      # Scale single trait. ==
      dfNorm[, contTraits] <- as.numeric(scale(dfNorm[, contTraits]))
    }
    # Introduce missing values from dfMissOrig into dfNorm. !!!
    dfMiss <- as.data.frame(mapply(function(x, y) ifelse(is.na(x), x, y), x = dfMissOrig[, c("species_name", contTraits)], y = dfNorm[, c("species_name", contTraits)], SIMPLIFY = F))
    # Add categorical traits back. !!!
    dfMiss <- merge(dfMiss, dfMissOrig[, c("species_name", catTraits)], by = "species_name")
    # Bind missingness indicator columns and reorganize the dataframe.
    dfMiss <- BindAndOrganize(dfMiss, vars)
    # Make sure dfMiss and dfNorm (the original data) are both ordered alphabetically.
    dfMiss <- dfMiss[order(dfMiss$species_name), ]
    dfNorm <- dfNorm[order(dfNorm$species_name), ]
    # If there is more than one categorical trait..
    if(length(catTraits) > 1){
      # Use lapply to convert categorical traits to factor class.
      dfNorm[, catTraits] <- lapply(dfNorm[, catTraits], as.factor)
      dfMiss[, catTraits] <- lapply(dfMiss[, catTraits], as.factor)
      # If there is only one categorical trait..
    } else if(length(catTraits) == 1) {
      # Convert single trait to factor class.
      dfNorm[, catTraits] <- as.factor(dfNorm[, catTraits])
      dfMiss[, catTraits] <- as.factor(dfMiss[, catTraits])
    }
    
    # Phylogenetic eigenvector decomposition. ---
    if(phyImp == T) {
      # Append eigenvectors to dfMiss and list of predictors. Each trait will have a corresponding dataframe and list of predictors including the eigenvectors.
      l_evs <- AppendEigenvectors(data = dfMiss, vars = traitsMAR, tree = tree, predictors = l_predictors)
      # Extract list of dfMiss.
      l_dfMissEV <- l_evs[[1]]
      # Extract updated list of predictors.
      l_EVPredictors <- l_evs[[2]]
    }
    # Identify outgroup(s) in trait datasets (this is the species that contains no trait data dfNorm). &&
    outgroup <- dfNorm$species_name[apply(dfNorm[, vars], 1, function(x) all(is.na(x)))]
    # If found in dataframe.. &&
    if(length(outgroup) > 0){
      # Remove outgroup from dataframes.
      dfNorm <- dfNorm[!dfNorm$species_name %in% outgroup, ]
      dfMiss <- dfMiss[!dfMiss$species_name %in% outgroup, ]
      if(phyImp == T){
        # Also remove outgroup from l_dfMissEV.
        l_dfMissEV <- lapply(l_dfMissEV, function(x) x[!x$species_name %in% outgroup, ])
      }
    }
    
    # Imputation. ---   
    # The ImputeRF() function entails a loop that uses different values of ntree (the number of trees grown in the forest) and returns a list of error rates (MSE for continuous traits and ARI for categorical traits) for each parameter value tested. +++
    # If phylogenetic imputation was selected..  
    if(phyImp == T) {
      # Impute data using eigenvectors. +
      impResult <- ImputeRF(dfTrue = dfNorm, dfMissing = dfMiss, dfOrig = data, cols = traitsMAR, cont = contTraits, cat = catTraits, inter = intTraits, ntrees = c(100, 1000), predictors = l_EVPredictors, phyImp = T, l_dfMissing = l_dfMissEV, seed = i)
    } else if(phyImp == F){
      # Impute data only using trait data. +
      impResult <- ImputeRF(dfTrue = dfNorm, dfMissing = dfMiss, dfOrig = data, cols = traitsMAR, cont = contTraits, cat = catTraits, inter = intTraits, ntrees = c(100, 1000), predictors = l_predictors, seed = i)
    }
    
    # Result handling. ---  
    # Extract l_dfImp.  
    l_dfImp <- impResult[[1]]
    # Extract errorRates  .
    l_Error <- impResult[[2]]
    # Extract corrCoefs. +++
    l_Corr <- impResult[[3]]
    
    # Back-transforming data. ---  
    # First, match data species (original species in complete-case) to species now in dfNorm in case any were removed (e.g. outgroups).
    dfOrig <- data[data$species_name %in% dfNorm$species_name, ]
    # Ensure dfOrig and dfImp are in same order. $$$
    dfOrig <- dfOrig[order(dfOrig$species_name), ]
    l_dfImp <- lapply(l_dfImp, function(x) x[order(x$species_name), ])
    
    # For every imputed dataframe..
    for(d in 1:length(l_dfImp)){
      # Take the dth dfImp.
      dfImp <- l_dfImp[[d]]
      # Apply BackTransform function to contTraits in dfImp. $$$
      dfImp <- BackTransform(origData = dfOrig, tfData = dfImp, cols = traitsMAR[traitsMAR %in% contTraits])
      # If there are any integer traits..
      if(length(intTraits) > 0){
        # Round to nearest whole number.
        dfImp[, intTraits] <- lapply(dfImp[, intTraits], function(x) as.integer(round(x)))
      }
      # Replace dfImp in l_dfImp with newly backtransformed dataset.
      l_dfImp[[d]] <- dfImp
    }
    
    # Append results to ith element of lists.    
    l_dfMissOrig[[i]] <- dfMissOrig
    l_l_missingness[[i]] <- l_missingness      
    l_l_dfImp[[i]] <- l_dfImp
    l_l_Error[[i]] <- l_Error
    l_l_Corr[[i]] <- l_Corr # +++
    # If phylogenetic imputation was selected..
    if(phyImp == T) {
      l_l_evs[[i]] <- l_evs
    }
  } ## i
  
  # Final trait check. --- +++
  # If there are any traits that could not be simulated MAR..
  if(length(traitsNOSIM) > 0){
    # If the trait was identified in a previous iteration to be a trait that couldn't be simulated MAR, it is possible missing values were introduced in other iterations. These will likely be very low numbers of NAs and not enough for error rate analyses. So here we will subset l_l_Error to ensure it only contains values for varsMAR.
    l_l_Error <- lapply(l_l_Error, function(x) {
      # Identify error rates associated with traitsMAR.
      index <- which(names(x) %in% traitsMAR)
      # Remove traits that were identified as MCAR in later iterations.
      x <- x[index]
    })
    print("The following traits could not be simulated MAR:")
    print(unique(traitsNOSIM))
    
    l_l_Corr <- lapply(l_l_Corr, function(x) {
      # Identify error rates associated with traitsMAR.
      index <- which(names(x) %in% traitsMAR)
      # Remove traits that were identified as MCAR in later iterations.
      x <- x[index]
    })
    
  }
  
  # Average out the missingness for each trait. ---    
  # Unlist missingness proportions.
  missingness <- unlist(l_l_missingness)
  # Calculate the average missingness proportion for each trait.
  avgMiss <- lapply(traitsMAR, function(x) {
    # Identify missingness proportions associated with the trait.
    index <- grep(pattern = x, names(missingness))
    # Take the mean.
    average <- mean(missingness[index], na.rm = T)
    # Name average according to the trait.
    names(average) <- x
    # Return the average missingness proportion.
    return(average)
  })
  # Subset to sample sizes for traits that could be simulated MAR. &&
  l_MARn <- l_sampleSizes[names(l_sampleSizes) %in% traitsMAR]
  # Order l_MARn by avgMiss.
  l_MARn <- l_MARn[names(unlist(avgMiss))]
  # Calculate the average number of NAs introduced for each trait. &&
  avgNAs <- mapply(function(x, y) x * y, y = l_MARn, x = avgMiss)
  # Rename according to l_MARn.
  names(avgNAs) <- names(l_MARn)
  
  # If phylogenetic imputation was selected..    
  if(phyImp == T) { 
    # Create list to hold the results.
    l_results <- list(completeCaseData = data, rawData = raw, traits = l_traits, numeric = contTraits, categorical = catTraits, integer = intTraits, traitsNotSimulated = traitsNOSIM, traitsSimulated = traitsMAR, MARfinalModels = l_MARfinalModels, ntrees = ntrees, reps = int, predictors = l_EVPredictors, tree = tree, eigenvectors = l_l_evs, missingData = l_dfMissOrig, missingness = l_l_missingness, averageMissingness = avgMiss, averageNAs = avgNAs, imputedData = l_l_dfImp, errorRates = l_l_Error, corrCoef = l_l_Corr)
  } else if(phyImp == F){
    # Create list to hold the results.
    l_results <- list(completeCaseData = data, rawData = raw, traits = l_traits, numeric = contTraits, categorical = catTraits, integer = intTraits, traitsNotSimulated = traitsNOSIM, traitsSimulated = traitsMAR, MARfinalModels = l_MARfinalModels, ntrees = ntrees, reps = int, predictors = l_predictors, missingData = l_dfMissOrig, missingness = l_l_missingness, averageMissingness = avgMiss, averageNAs = avgNAs, imputedData = l_l_dfImp, errorRates = l_l_Error, corrCoef = l_l_Corr)
  }
  # Return l_results.  
  return(l_results)
  
}