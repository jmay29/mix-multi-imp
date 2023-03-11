# Data handling functions for preparing trait and sequence data for imputation runs.

BreakIntoTypes <- function(data, traitCols) {
  
  # Function for identifying trait types of columns in a dataframe.
  
  # data = dataframe containing trait and taxonomic data
  # traitsCols = columns containing trait data
  
  # Identify categorical traits. Specifying drop = F for single traits so we can preserve dataframe format.
  catTraits <- GetTraitNames(data = data[, traitCols, drop = F], class = "character")
  # Identify binary traits and append them to catTraits as they will also be treated as factors. Omitting NAs so NA isn't considered its own category.
  binTraits <- GetTraitNames(data = na.omit(data[, traitCols, drop = F]), class = "binary")
  catTraits <- c(binTraits, catTraits)
  # Convert categorical traits to factors to differentiate them from numeric traits.
  # If there is more than 1 categorical trait..
  if(length(catTraits) > 1){
    # Use lapply.
    data[, catTraits] <- lapply(data[, catTraits], as.factor)
    # if there's only one categorical trait..
  } else if(length(catTraits) == 1){
    data[[catTraits]] <- as.factor(data[[catTraits]])
    # If there aren't any categorical traits..
  } else if(length(catTraits) == 0){
    # Assign NA to catTraits.
    catTraits <- NA
  }
  
  # Identify numerical traits.
  contTraits <- GetTraitNames(data = data[, traitCols, drop = F], class = "numeric")
  # If there aren't any numerical traits..
  if(length(contTraits) == 0){
    # Assign NA to contTraits.
    contTraits <- NA
  }
  
  # Combine contTraits and catTraits into a list.
  l_traits <- list(contTraits, catTraits)
  # Name according to trait type.
  names(l_traits) <- c("Numerical", "Categorical")
  
  # Return l_traits.
  return(l_traits)
  
}

CleanColumns <- function(data) {
  
  # Function for printing column names of a dataframe and cleaning the names using the "janitor" package.
  # data = dataframe
  
  # Print original column names.
  print("Original column names:")
  print(colnames(data))
  # Convert column names of data into snake case (e.g. var_one).
  dfClean <- janitor::clean_names(data)
  # Print new column names.
  print("New column names:")
  print(colnames(dfClean))
  
  # Return dataframe with cleaned names.
  return(dfClean)
  
}

CleanDataframe <- function(data, orderName) {
  
  # Function for subsetting dataframe for a certain order and to only retain those rows with species name information.
  # data = dataframe with taxonomy and sequence information
  # orderName = character vector containing name of order on which to subset the data
  
  # Subsetting the data by order and selecting only relevant columns.
  dfCleaned <- subset(data, order == orderName, select = c(species_name, processid, COI))
  # Removing NAs.
  dfCleaned <- na.omit(dfCleaned)
  # Remove rownames.
  rownames(dfCleaned) <- NULL
  
  return(dfCleaned)
  
}

CountByGroup <- function(data, orderName){
  
  # Tidyverse function to tally number of species per family and return a dataframe of n/proportions.
  
  # data = dataframe
  # orderName = order name to append to a column called order_name
  
  dfCount <- data %>%
    # Group by family.
    group_by(family) %>%
    # Tally number of species per family.
    tally() %>%
    # Add a column called prop that is the number of species/total.
    mutate(prop = n/sum(n)) %>%
    # Add a column called order_name.
    mutate(order_name = orderName)
  
  return(dfCount)
  
}

CountVotes <- function(vector, string){
  
  # Function for counting number of instances of a string in a dataframe column.
  # vector = character vector
  # string = character string to count
  
  # Find pattern in vector using grep. Get length of result.
  votes <- length(grep(string, vector))
  
  return(votes)
  
}

CreateNamedList <- function(listLength, elementNames) {
  
  # Function for creating a named list.
  # listLength = length of list
  # elementNames = vector element containing names for the list
  
  # Create empty list with length = listLength.
  l_list <- vector(mode = "list", length = listLength)
  # Name the list with provided vector.
  names(l_list) <- elementNames
  
  # Return named list.
  return(l_list)
  
}

CreateTraitSubset <- function(df, traits, taxonomy, speciesColumn, critNumber) {
  
  # This function identifies the traits that maximize the number of complete cases and returns a complete-case dataframe.
  
  # df = dataframe containing trait information
  # traits = list of column names to consider (should be names of traits)
  # taxonomy = column names that contain taxonomy information
  # speciesColumn = character vector that contains name of species column in dataframe
  # critNumber = the minimum sample size required for the complete-case subset
  
  # Get the order name for the taxon.
  orderName <- unique(df$order)
  # Use the IdentifyTraits to identify how many complete cases we get for each trait subset.
  suppressWarnings(traitNum <- IdentifyTraits(df, traits))
  # How many of the traits create subsets greater than critNumber?
  goodTraits <- names(which(traitNum > critNumber))
  
  # If goodTraits contains more than 3 traits...
  if(length(goodTraits) >= 3) {
    
    # Add the taxonomy info back.
    goodTraits <- c(taxonomy, goodTraits)
    # Get the complete-case dataset using the TakeCompleteCases function.
    dfComplete <- TakeCompleteCases(df, goodTraits)
    # Getting it ready for merging with sequence data. Renaming Species column and adding underscores to species names.
    speciesIndex <- which(colnames(dfComplete) == speciesColumn)
    colnames(dfComplete)[speciesIndex] <- "species_name"
    dfComplete$species_name <- gsub(" ", "_", dfComplete$species_name)
    
    # Return the complete-case subset.
    return(dfComplete)
    
  } else if(length(goodTraits) < 3) {
    
    print(paste("Not enough trait data for", orderName, "at this sample size:", critNumber, sep = " "))
    
  }
  
}

CritNumberLoop <- function(critNumbers =  seq(100, 5000, by = 100), l_df, traitCols, taxCols) {
  
  # Function to create complete-case datasets for trait data, considering different criterion for the required sample size.
  
  # numbers = numeric vector, minimum must be 100
  # l_df = list of dataframes containing trait data
  # traitCols = column names containing trait data
  # taxCols = column names containing taxonomy data
  
  # For each value of critNumber...
  for(i in 1:length(critNumbers)) {
    
    # Take the ith critNumber.
    num <- critNumbers[[i]]
    # Create complete-case subsets for each taxon using the CreateTraitSubset() function.
    l_dfCompleteCases <- lapply(l_df, CreateTraitSubset, traits = traitCols, taxonomy = taxCols, speciesColumn = "species_name", critNumber = num)
    # Subset the list to dataframe objects using the SubsetListByClass() function.
    l_dfCompleteCases <- SubsetListByClass(l_dfCompleteCases, "data.frame")
    # For those taxa that meet the sample size criterion, write the dataframes to file using the WriteCompleteCases() function.
    lapply(l_dfCompleteCases, WriteCompleteCases)
    
  }
  
}

ExtractTraitDataset <- function(l_df, n = 100, orderName) {
  
  # Function for extracting a dataframe from a list of dataframes and writing it to file, if its sample size (nrow) is greater than n.
  # l_df = list of dataframes
  # n = number of rows in the dataframe. Default is 100
  # orderName = order to filter for
  
  # Extract the largest dataset from the list, if it has a sample size greater than 100.
  sampleSizes <- sapply(l_candidateFiles, nrow)
  # If any of the sample sizes are greater than n...
  if (max(sampleSizes) > n) {
    
    # What's the index of the largest sample size?
    candidate <- which.max(sampleSizes)
    # Subset this dataframe from the list of candidate dataframes.
    dfCandidate <- l_candidateFiles[[candidate]]
    # Create a file name for the dataframe.
    fileName <- paste(orderName, "COI_Trait.csv", sep = c("_"))
    # Write it to file.
    write.csv(dfCandidate, fileName)
    
    print(paste("Wrote", fileName, "to file", sep = " "))
    
    # Return the dataframe with the largest sample size.
    return(dfCandidate)
    
  } else {
    
    print(paste("None of the dataframes have a sample size greater than", n, sep = " "))
    
  }
  
}

GetCategoricalInfo <- function(vector) {
  
  # Function for providing descriptive information for a categorical (character) vector.
  # vector = character vector
  
  # Create empty list to hold results.
  result <- vector(mode = "list", length = 5)
  # Name the list.
  names(result) <- c("n", "NumberOfCategories", "Frequencies", "Proportions")
  
  # Get the sample size (n) of the data.
  result[[1]] <- length(na.omit(vector))
  # Get the number of categories (excluding NAs).
  result[[2]] <- length(unique(na.omit(vector)))
  # Get the categorical frequencies including NAs.
  result[[3]] <- table(vector, useNA = "ifany")
  # Get the categorical proportions including NAs.
  result[[4]] <- table(vector, useNA = "ifany")/length(vector)
  
  # Return the list of results.
  return(result)
  
}

GetCounts <- function(data, var, taxCol){
  
  # Function for determining range, min, and max counts of vector. Also returns order names with min/max counts for associated vector
  # data = dataframe
  # var = column name on which to perform calculations
  # taxCol = column name containing taxonomy information
  
  # Ensure dataframe format.
  data <- as.data.frame(data)
  
  # Range of numerical traits.
  print(paste("Range of", var, ":", sep = " "))
  print(range(data[, var]))
  # Minimum value:
  print(paste("Max of", var, ":", max(data[, var]), sep = " "))
  print(data[, taxCol][data[, var] == max(data[, var])])
  # Maximum value:
  print(paste("Min of", var, ":", min(data[, var]), sep = " "))
  print(data[, taxCol][data[, var] == min(data[, var])])
  
}

GetNumericalInfo <- function(vector) {
  
  # Function for providing descriptive information for a numerical vector.
  # vector = numerical vector
  
  # Create empty list to hold results.
  result <- vector(mode = "list", length = 7)
  # Name the list.
  names(result) <- c("n", "Mean", "Median", "Range", "NumberNAs", "ProportionNAs", "Plot")
  
  # Get the sample size (n) of the data.
  result[[1]] <- length(na.omit(vector))
  # Get the mean of the data.
  result[[2]] <- mean(vector, na.rm = T)
  # Get the median of the data.
  result[[3]] <- median(vector, na.rm = T)
  # Get the range of the data.
  result[[4]] <- range(vector, na.rm = T)
  # Get the numer of NAs in the data.
  result[[5]] <- sum(is.na(vector))
  # Get the proportion of NAs in the data.
  result[[6]] <- sum(is.na(vector))/length(vector)
  # Plot the data.
  result[[7]] <- hist(vector, col = "skyblue")
  
  # Return the list of results.
  return(result)
  
}

GetTraitNames <- function(data, class = "numeric") { 
  
  ## TODO: Fix warning message if single trait provided.
  
  # Function for extracting trait names from a dataframe that are of a certain class.
  
  # data = dataframe containing trait data
  # class = type of class for which to extract names of trait columns (either "numeric", "binary", or "character").
  
  if(class == "numeric") {
    # Apply is.numeric() function across columns.
    cons <- sapply(data, is.numeric)
    # Get the names of the columns that are TRUE.
    traits <- names(which(cons == TRUE))
  } else if(class == "binary") {
    # Apply is.binary() function from the "Information" package across columns.
    bin <- sapply(data, Information::is.binary)
    # Get the names of the columns that are TRUE.
    traits <- names(which(bin == TRUE))
  } else if(class == "integer") {
    # Apply is.character() function across columns.
    ints <- sapply(data, is.integer)
    # Get the names of the columns that are TRUE.
    traits <- names(which(ints == TRUE))
  } else if(class == "character") {
    # Apply is.character() function across columns.
    chars <- sapply(data, is.character)
    # Get the names of the columns that are TRUE.
    traits <- names(which(chars == TRUE))
  }
  
  # Return names of the traits.
  return(traits)
  
}

FixTraits <- function(data, columns) {
  
  # Function for fixing binary/multicategorical trait columns. Flips NAs to 0s.
  # data = dataframe containing trait information
  # columns = column numbers
  
  # Replace all of the NAs in these columns with 0s.
  data[, columns] <- lapply(data[, columns], function(x) ifelse(is.na(x) == T, 0, 1))
  # Identify species with 0s for all columns and switch to NAs.
  missing <- apply(data[, columns], 1, function(x) all(x == 0))
  missing <- which(missing == TRUE)
  data[missing, columns] <- NA
  
  return(data)
  
}

FixSpeciesColumn <- function(data, columnNum) {
  
  # Function for renaming column that contains species name information and inserting an underscore.
  
  # data = dataframe that contains species name information
  # columnNum = number of column that contains species name information
  
  # Rename the column to "species_name".
  colnames(data)[columnNum] <- "species_name"
  # Insert an underscore.
  data$species_name <- gsub(" ", "_", data$species_name)
  
  return(data)
  
}

IdentifyTaxa <- function(list, traits, numSpecies = 100, numTraits = 3) {
  
  # Function for identifying taxa with complete information available.
  # list = list of dataframes containing trait information. Each dataframe should correspond to a different order and be named by that order.
  # trait = column names that correspond to trait data
  # numSpecies = minimum number of observations (species) that should have complete trait information
  # numTraits = minimum number of traits that should have complete observations
  
  # Apply the IdentifyTraits function over the list of dataframes. This will provide information about which traits will maximize the complete cases in the dataset.
  completeTraits <- lapply(list, IdentifyTraits, variables = traits)
  # Which orders have the minimum criteria?
  # numSpecies (rows)
  traits100 <- lapply(completeTraits, function(x) sum(x > numSpecies))
  # numTraits (columns)
  goodTaxa <- names(which(traits100 > numTraits))
  
  # Return the names of the taxa that met the minimum criteria.
  return(goodTaxa)
  
}

IdentifyTraits <- function(data, variables) {
  # Function for identifying the traits that maximize the number of complete-cases.
  
  # data = dataframe containing species and trait data
  # variables = vector of trait names (should match column names in data)
  
  # Convert data to data.table format.
  dfMultivariable <- as.data.table(data[, variables])
  # Order the columns by the amount of missing data (NA values).
  suppressWarnings(dfTraitsNA <- sort(dfMultivariable[, lapply(.SD, function(x) sum(is.na(x)))]))
  # Reorder the dataset so that the columns with the least amount of NA values are now first.
  setcolorder(dfMultivariable, names(dfTraitsNA))
  
  # Now I want to loop through the traits, removing one column (trait) at a time and count the number of complete cases. This will provide us some information as to which traits would provide an adequate sample size for downstream analysis.
  # Take the number of columns in dfMultivariable.
  len <- ncol(dfMultivariable)
  # Create a numeric vector to hold the results of the loop.
  all.cc <- NULL
  # Start the loop:
  for (i in 1:len) {
    # Works best if you set dfMultivariable back to a dataframe.
    x <- as.data.frame(dfMultivariable)
    # x is the placeholder dataframe in the loop.
    x <- x[, 1:len]
    # Determine which rows are "complete" using the "len" subset of traits.
    x <- complete.cases(x)
    # Complete rows of data will be "TRUE".
    x <- which(x == "TRUE")
    # Find the number of complete cases.
    x <- length(x)
    # Add it to the all.cc variable that's holding all of the results of the loop.
    all.cc[i] <- x
    # Minus 1 from tempLen so we can check the next subset of traits (we started at the last column because the columns were ordered by number of NA values).
    len <- len - 1
  }
  # Now, decide where to cut the datatable. (i.e. pick an adequate subset of traits that maximize sample size).
  names(all.cc) <- rev(colnames(dfMultivariable))
  # Look at the results.
  all.cc
  
  # Return the vector that contains information about how to maximize complete cases.
  return(all.cc)
  
}

OutlierCheck <- function(data, col){
  
  # Function for identifying outliers based on the interquartile range (IQR) of the data. Returns a subsetted dataframe of the outliers so user can check the values for potential error.
  # data = dataframe. Make sure dataframe also contains a column called "species_name" with species name information.
  # col = column name for which to identify outliers
  
  # Code for calculating thresholds adapted from: 
  # Author: https://stackoverflow.com/users/1312519/by0.
  # https://stackoverflow.com/questions/12866189/calculating-the-outliers-in-r.
  
  # Determine the 25% quantile.
  lowerQuantile <- quantile(data[, col], na.rm = T)[2]
  # Determine the 75% quantile.
  upperQuantile <- quantile(data[, col], na.rm = T)[4]
  # Calculate the IQR.
  iqr <- upperQuantile - lowerQuantile
  # Calculate upper threshold ((3 x the IQR) + upperQuantile).
  upperThreshold <- (iqr * 3) + upperQuantile
  # Calculate lower threshold (lowerQuantile - (3 x the IQR)).
  lowerThreshold <- lowerQuantile - (iqr * 3)
  
  # Identify outliers based on whether they exceed the upper or lower thresholds.
  outliers <- which(data[, col] > upperThreshold | data[, col] < lowerThreshold)
  
  # Subset the outliers. Ensure species name information is also kept.
  dfOutliers <- data[outliers, c("species_name", col)]
  
  if(nrow(dfOutliers) == 0){
    print(paste("No outliers detected for:", col, sep = " "))
  } else {
    print(paste("Outliers detected for:", col, sep = " "))
    # Return the dataframe of outliers and their corresponding extreme values.
    return(dfOutliers)
  }
  
}

OverHeadSampleSizeCheck <- function(orderName, traitFileName, dfSeq) {
  
  # Function for performing sample size check (SampleSizeCheck.R) given the order name and name of file containing trait data.
  # orderName = name of order
  # traitFileName = name of file containing trait data for the order
  # dfSeq = dataframe containing COI sequence data
  
  # Subset the dataframe for the order.
  dfOrderCOI <- CleanDataframe(dfSeq, orderName)
  # Read in the trait data.
  dfTraits <- read.csv(traitFileName)
  # Clean up the trait dataframe.
  dfTraits$X <- NULL
  
  # Merge gene and trait datasets together.
  dfCOI_Trait <- SampleSizeCheck(dfOrderCOI, dfTraits, variable = "species_name", 100)
  
  # For now, remove the sequence data column because in our next script we are finding centroid sequences.
  dfCOI_Trait$COI <- NULL
  
  return(dfCOI_Trait)
  
}

MergeAndClean <- function(df1, df2, colNames, newNames = colNames) {
  
  # Function for merging two dataframes together and rearranging columns.
  # df1 = first dataframe to merge
  # df2 = second dataframe to merge
  # colNames = column names to keep in final merged dataframe
  # newNames = new columns names for final merged dataframe. Should be same length as colNames. Default is to keep original colnames.
  
  df1 <- as.data.frame(df1)
  df2 <- as.data.frame(df2)
  
  # Merge the dataframes.
  dfMerge <- merge(df1, df2, by = "species_name")
  # Subset by colNames.
  dfMerge <- dfMerge[, colNames]
  # Rename columns.
  names(dfMerge) <- newNames
  
  return(dfMerge)
  
}

MergeAndUpdate <- function(df1, df2, orderCol1, orderCol2) {
  
  # Function for merging two dataframes together and updating order-level information in the merged dataframe.
  # df1 = dataframe 1
  # df2 = dataframe 2
  # orderCol1 = column name with order-level information in df1
  # orderCol2 = column name with order-level information in df2
  
  # Merge df1 and df2.
  dfMerge <- merge(df1, df2, by = "species_name", all = T)
  # Update order-level information in dfMerge.
  dfMerge[, orderCol1] <- ifelse(is.na(dfMerge[, orderCol1]) == T, dfMerge[, orderCol2], dfMerge[, orderCol1])
  
  # Return the merged dataframe with updated order-level information.
  return(dfMerge)
  
}

PasteColNames <- function(data, colNames, string, ...) {
  
  # Function for pasting a string onto a subset of column names in a dataframe.
  # data = dataframe
  # colNames = column names to paste string onto
  # string = string to paste onto column names
  # ... = arguments for paste
  
  # Get trait col indices.
  index <- which(colnames(data) %in% colNames)
  # Paste col names so we can identify complete data.
  colnames(data)[index] <- paste(colnames(data[, colNames]), string, ...)
  
  # Return dataframe with new column names.
  return(data)
  
}

ReplaceValues <- function(data, col, patterns, replacements){
  
  # Function to replace values in a specific column in dataframe.
  # data = dataframe
  # col = column name in which to replace values
  # patterns = vector of pattern names to replace
  # replacements = replacement patterns for corresponding elements in patterns
  
  # For every pattern..
  for(p in 1:length(patterns)){
    # Replace pattern with replacement in target column.
    data[, col] <- gsub(pattern = patterns[p], replacements[p], x = data[, col])
  }
  
  # Return dataframe with updated column.
  return(data)
  
}

RmNAAll <- function(data, cols){
  
  # Function to remove rows with all NAs in target columns.
  # data = dataframe
  # cols = column names to consider when identifying NAs
  
  # Identify rows with all NAs.
  index <- apply(data[, cols], MARGIN = 1, function(x) all(is.na(x)))
  # Remove from dataframe.
  dfNew <- data[!index, ]
  
  return(dfNew)
  
}

SampleSizeCheck <- function(df1, df2, variable, number) {
  
  # Function for determining the number of complete cases when merging two dataframes together.
  # df1 = first dataframe
  # df2 = second dataframe
  # variable = column name to merge by
  # number = threshold sample size must meet
  
  
  # Merge dataframes together.
  dfMerge <- merge(df1, df2, by = variable)
  
  if(nrow(dfMerge) > number) {
    
    print("PERFECT!")
    print(nrow(dfMerge))
    
  } else if (nrow(dfMerge) < number) {
    
    print("INADEQUATE!")
    print(nrow(dfMerge))
    
    
  } else {
    
    print("This is NOT a dataframe... ?")
    
  }
  
  return(dfMerge)
  
}

ScreenCategories <- function(variable, threshold = 0.10) {
  
  # Function for determining whether a variable contains categories with observations below a certain threshold.
  
  # variable = character/factor type vector
  # varName = names of variable
  # threshold = numeric threshold (e.g. 0.10)
  
  # Count the number of observations per category.
  counts <- plyr::count(variable)
  # Calculate the proportions for each category.
  proportions <- counts$freq/sum(counts$freq)
  # Name the proportions according to category.
  names(proportions) <- counts$x
  # Create an empty list to hold categories that meet the required threshold.
  goodCat <- vector(mode = "numeric")
  
  print("Results:")
  
  # # For every category...
  for(i in 1:length(proportions)) {
    # Take the ith category.
    category <- proportions[i]
    # If the category comprises less than 10% of the dataset...
    if(category < threshold) {
      print(paste(category, "is less than", threshold, "of dataset. Remove", names(category), "!"))
    } else if(category > threshold) {
      print(paste(category, "is greater than", threshold, "of dataset. Keep", names(category), "!"))
      # Append to goodCat.
      goodCat[[i]] <- category
      # Name according to category.
      names(goodCat)[[i]] <- names(category)
    }
  }
  
  # Remove NAs from goodCat.
  goodCat <- na.omit(goodCat)
  
  # Check if there is more than one level in the variable. Otherwise the whole trait should be removed because it is invariant.
  if(length(goodCat) == 1) {
    print("This variable only contains one category. Removal from dataset is recommended.")
  } else {
    print(paste("This variable contains", length(goodCat), "categories. It can remain in the dataset."))
  }
  
  # Return the proportions for each category.
  return(goodCat)
  
}

SpeciesNameCheck <- function(data) {
  
  # Function for ensuring there are no observations with missing species name information and no species duplicates.
  # data = dataframe with "species_name" column
  
  # Replace all blanks with NAs.
  data[data == ""] <- NA
  # Remove all rows where species_name info is missing.
  data <- data[!is.na(data$species_name) == T, ]
  
  # Check that are no rows with missing species info.
  if(sum(is.na(data$species_name)) == 0) {
    print("No rows with missing species info! :)")
  } else if (sum(is.na(data$species_name)) > 0) {
    print("There are rows with missing species info! :(")
  }
  # Check that are no duplicate species names.
  if(sum(duplicated(data$species_name)) == 0) {
    print("No species duplicates! :)")
  } else if (sum(duplicated(data$species_name)) > 0) {
    print("There are species duplicates! :(")
  }
  # Return data.
  return(data)
  
}

SplitAndBreak <- function(data, colName, name = NULL) {
  
  # Function for splitting a dataframe into a list, grouped by a particular variable (column). The list is then broken up into its respective components and added to the global environment.
  # data = dataframe to split
  # colName = column (grouping variable) to split dataframe by
  # name = character string to append to names of list elements
  
  # Split the dataframe into a list, grouped by colName.
  l_df <- split(data, data[, colName])
  # Append name to the names of the list elements.
  names(l_df) <- paste(name, names(l_df), sep = "")
  # list2env is a handy function for breaking down a list into its respective components and into the global environment. They are named according to the list elements.
  list2env(l_df, .GlobalEnv)
  
  # Return list of dataframes.
  return(l_df)
  
}

SubsetAndMerge <- function(list, ...) {
  
  # Function for subsetting a list of dataframes by order name (which is contained in the species_name column pertaining to the mammal alignment) and merging the list of dataframes.
  
  # Subsetting for order name.
  l_dfSubset <- lapply(list, function(x) SubsetByOrder(df = x, ...))
  # Merge all of the dataframes using Reduce().
  dfMerge <- Reduce(function(...) merge(..., by = "species_name", all = T), l_dfSubset)
  
  return(dfMerge)
  
}

SubsetByOrder <- function(df, string) {
  
  # Function for subsetting a dataframe based on string match in seq.name column (messy record title column). Also cleaning up the name a bit.
  
  # df = dataframe that contains a seq.name column (but in this case it is messy record title for a specific alignment)
  # string = string match that I want to subset by in that column
  
  # Subset by string using grepl() function.
  dfSubset <- df[grepl(string, df$seq.name), ]
  # Replacing underscores with spaces.
  dfSubset$seq.name <- gsub(pattern = "_", replacement = " ", x = dfSubset$seq.name)
  # Taking the first two words of the string, which correspond to the species names for which the record was obtained.
  dfSubset$seq.name <- word(dfSubset$seq.name, 1, 2)
  # Renaming seq.name to species_name.
  colnames(dfSubset)[1] <- "species_name"
  # Add the underscores back into species_name to match trait data.
  dfSubset$species_name <- gsub(pattern = " ", replacement = "_", x = dfSubset$species_name)
  
  return(dfSubset)
  
}

SubsetList <- function(list, names) {
  
  # Function for subsetting a list.
  # list = named list
  # names = names of the elements you want to subset
  
  # Identify the index positions in the list.
  index <- which(names(list) %in% names)
  # Subset the list based on index positions.
  l_df <- list[index]
  
  return(l_df)
  
}

SubsetListByClass <- function(list, targetClass) {
  
  # Function for subsetting a list by object class.
  
  # list = list object
  # targetClass = class to subset by
  
  # Get the classes of the list elements.
  classes <- lapply(list, class)
  # Get the indices of the dataframes.
  indices <- which(classes == targetClass)
  # Subset out the dataframes.
  subset <- list[indices]
  
  # Return the subsetted list.
  return(subset)
  
}

TakeCompleteCases <- function(data, variables) {
  # Function for subsetting dataset based on column names providing and removing NAs.
  # data = dataframe containing species and trait information
  # variables = column names to subset
  
  # Subset the column names in variables.
  dfSubset <- subset(data, select = variables)
  # Take only the complete-cases.
  dfSubset <- na.omit(dfSubset)
  
  return(dfSubset)
  
}

WriteCompleteCases <- function(data) {
  
  # Function for writing complete case dataset to file.
  # df = dataframe containing complete case trait information. Should be at the order level.
  
  # Create file name.
  fileName <- paste(unique(data$order), nrow(data), ".csv", sep = "")
  
  # Write the dataframe to file.
  write.csv(data, fileName)
  print(paste("Wrote records to ", fileName, sep = ""))
  
}

WriteDataframeToFile <- function(df, orderCol, fileName) {
  
  # Function for writing dataframe to file.
  # df = dataframe containing complete case trait information. Should be at the order level.
  # orderCol = column that contains order level info
  # fileName = name to give to csv file
  
  # Create file name.
  orderName <- names(sort(table(df[[orderCol]]), decreasing = T))[1]
  fileName <- paste(orderName, fileName, ".csv", sep = "")
  
  # Write the dataframe to file.
  write.csv(df, fileName)
  print(paste("Wrote records to ", fileName, sep = ""))
  
}
