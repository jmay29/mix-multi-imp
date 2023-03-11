# Functions related to phylogenetic analyses.

AppendEigenvectors <- function(data, vars, tree, predictors){
  
  # Function for appending eigenvectors to dataframe containing trait data with missing values. Returns a list of dataframes with corresponding eigenvectors for each trait.
  
  # data = dataframe with missing values
  # vars = names of traits (columns in data)
  # tree = phylo object
  # predictors = list of predictors for each trait
  
  # Create a list to hold dataframes of predictors and eigenvectors for each trait.
  l_dfMiss <- CreateNamedList(listLength = length(vars), elementNames = vars)
  # Create a list to hold names of updated predictors.
  l_newPreds <- CreateNamedList(listLength = length(vars), elementNames = vars)
  
  # For every trait..
  for(t in 1:length(vars)){
    # Take the tth trait.
    trait <- vars[[t]]
    # Make a complete-case subset of trait data.
    dfTrait <- na.omit(data[, c("species_name", trait)])
    # Decompose phylogenetic distrance matrix into orthogonal vectors (phylogenetic eigenvectors).
    PEMtrait <- DecomposeTree(tree = tree, data = dfTrait, trait = trait, method = "PEM", evselection = "cutoff", cutoff = 75)
    # Update the PEM using derived parameters and extract the eigenvectors so we also have them for species with missing data (cite: Johnson et al. 2021).
    dfTraitEVs <- UpdatePEM(tree, alpha = PEMtrait$a, rate = PEMtrait$psi)
    # Subset to those eigenvectors selected according to the cutoff threshold.
    dfTraitEVs <- dfTraitEVs[, names(dfTraitEVs) %in% names(PEMtrait[[1]])]
    # Identify column numbers of the eigenvectors.
    index <- grep("V_", colnames(dfTraitEVs))
    # Append name of trait to eigenvector names.
    colnames(dfTraitEVs)[index] <- paste(colnames(dfTraitEVs)[index], trait, sep = "_")
    # Assign eigenvector names to a variable (without species_name).
    eigenvectors <- colnames(dfTraitEVs[, -ncol(dfTraitEVs)])
    # Identify the predictors for the trait in question.
    index <- grep(trait, names(predictors))
    # Update the names of the predictors.
    newPreds <- c(predictors[[index]], eigenvectors)
    # Append the eigenvectors to data.
    l_dfMiss[[t]] <- merge(data, dfTraitEVs, by = "species_name")
    # Append the updated predictors.
    l_newPreds[[t]] <- newPreds
  }
  # Return list of dataframes and modified predictors.
  return(list(l_dfMiss, l_newPreds))
  
}

Calibrate <- function(phylo, tips = "root", minAge = 1, maxAge = 1, fileName, ...) {
  
  # Function for calibrating tree given known fossil dates (calibration points).
  # Uses apes "chronos" function (Paradis).
  # phylo = phylogenetic tree with branch lengths representing expected number of subs/per.
  # ageMin = 
  # tips = character vector containing names of taxa for which fossil data info is available for calibration, default is to calibrate the root
  # minAge = estimated minimum age of node, default is 1
  # maxAge = estimated maximum age of node, default is 1
  # fileName = name of file to write tree to in Newick format
  
  # Use makeChronosCalib function to create a dataframe of calibration points.
  calib <- makeChronosCalib(phylo, node = "root", age.min = minAge, age.max = maxAge)
  # Use chronos function to estimate a chronogram. Default parameters are lambda = 1 (smoothing parameter) and model = "correlated" clock. Setting control = chronos.control uses default parameters for tol, iter.max, eval.max, nb.rate.cat, dual.iter.max, and epsilon. See ?chronos for further details.
  chrono <- chronos(phylo, calibration = calib, control = chronos.control(), ...)
  # Write chronogram to Newick format.
  write.tree(chrono, fileName)
  print(paste("Wrote file to", fileName))
  
  # Return chronogram of class "chronos" and "phylo".
  return(chrono)
  
}

DecomposeTree <- function(tree, method = "PEM", evselection = "cutoff", cutoff = 75, data, trait = NULL) {
  
  # Function for decomposing phylogenetic distance matrix into orthogonal vectors using the PVRdecomp (package "PVR") or PEM.fitSimple (package "MPSEM") functions.
  
  # Code adapted from:
  # PVR documentation: Santos T. 2018. https://cran.r-project.org/web/packages/PVR/PVR.pdf
  # MPSEM tutorial: G Geunard. 2019. A phylogenetic modelling tutorial using Phylogenetic Eigenvector Maps (PEM) as implemented in R package MPSEM (0.3-6).
  
  # tree = object of class "phylo". Branch lengths should be in expected # of subs/site/unit of time (ultrametric tree).
  # method = one of "PVR" 'phylogenetic eigenvector regression' or "PEM" phylogenetic eigenvector mapping
  # evselection = one of "cutoff" or "stepwise" methods for selecting eigenvector number. If cutoff is selected, a threshold variability method is used; if "stepwise" is selected, stepwise regressiod is used (function lmforwardsequentialAICc from the MPSEM package)
  # cutoff = number indicating the variability cutoff for the eigenvalues for eigenvector selection (default is 65)
  # data = dataframe containing trait data
  # trait = when PEM is specified, indicates target trait (response variable)
  
  # Make sure the trait data and tree tips match.
  tree <- drop.tip(phy = tree, tip = tree$tip.label[!tree$tip.label %in% data$species_name])
  # Make sure the order matches between the phylogeny and dataframe.
  data <- data[match(tree$tip.label, data$species_name), ]
  # Ensure data is in data.frame format.
  data <- as.data.frame(data)
  
  if(method == "PVR") {
    # Decomposing phylogenetic distance matrix into orthogonal vectors.
    matrixPhylo <- PVRdecomp(tree)
    # Extract the eigenvalues so we can determine how many eigenvectors we should be using.
    dfEigenvalues <- as.data.frame(matrixPhylo@Eigen$values) ## Place eigenvalues into dataframe.
    # How much variability does the eigenvalue explain?
    dfEigenvalues$variability <- dfEigenvalues$`matrixPhylo@Eigen$values`/sum(dfEigenvalues$`matrixPhylo@Eigen$values`)
    # Get the percentage.
    dfEigenvalues$percentage <- dfEigenvalues$variability*100
    # Get the cumulative percentage for each eigenvalues.
    dfEigenvalues$cumsum <- cumsum(dfEigenvalues$percentage)
    # Plot the cumulative percentages.
    plot(dfEigenvalues$cumsum)
    # Add a horizontal line at cutoff.
    abline(h = cutoff, col = "red")
    # How many eigenvalues explain the variability?
    evs <- which(dfEigenvalues$cumsum <= cutoff)
    EVN <- length(evs)
    # Now we need to extract the eigenvectors for each species.
    # Extract the eigenvalues so we can determine how many eigenvectors we should be using.
    dfEigenvectors <- as.data.frame(matrixPhylo@Eigen$vectors) ## Place eigenvalues into dataframe.
    # Add the species names to this dataframe so we can merge with the trait data later on.
    dfEigenvectors$species_name <- matrixPhylo@phylo$tip.label
    # Subset according to the number of eigenvalues we determined previously.
    dfEigenvectors <- dfEigenvectors[, c(1:EVN, ncol(dfEigenvectors))]
    # Return dataframe of eigenvectors.
    return(dfEigenvectors)
    
  } else if(method == "PEM") {
    # Convert phylogenetic tree to graph-class object needed to build phylogenetic eigenvector map (PEM).
    graphPhylo <- Phylo2DirectedGraph(tree)
    # Calculate phylogenetic eigenvector maps (PEMs). Specifying PEM.fitSimple allows to estimate values of a (steepness parameter) and psi (relative evolution rate along edges).
    PEM <- PEM.fitSimple(
      y = data[, trait],
      x = NULL,
      w = graphPhylo)
    # If the cutoff method for eigenvector selection was selected..
    if(evselection == "cutoff"){
      # Extract the singular values so we can determine how many eigenvectors we should be using. 
      dfEigenvalues <- as.data.frame(PEM$d) ## Place singular values into dataframe.
      # Square to convert to obtains eigenvalues.
      dfEigenvalues$`PEM$d` <- dfEigenvalues$`PEM$d`*dfEigenvalues$`PEM$d`
      # How much variability does the eigenvalue explain?
      dfEigenvalues$variability <- dfEigenvalues$`PEM$d`/sum(dfEigenvalues$`PEM$d`)
      # Get the percentage.
      dfEigenvalues$percentage <- dfEigenvalues$variability*100
      # Get the cumulative percentage for each eigenvalues.
      dfEigenvalues$cumsum <- cumsum(dfEigenvalues$percentage)
      # Plot the cumulative percentages.
      plot(dfEigenvalues$cumsum)
      # Add a horizontal line at cutoff.
      abline(h = cutoff, col = "red")
      # How many eigenvalues explain the variability?
      evs <- which(dfEigenvalues$cumsum <= cutoff)
      EVN <- length(evs)
      # If number of eigenvectors is greater than 20..
      if(EVN > 20){
        # Limit to 20 eigenvectors so we can reduce chance of overfitting imputation models.
        EVN <- 20
      }
      # Now we need to extract the eigenvectors for each species.
      dfEigenvectors <- as.data.frame(PEM)
      dfEigenvectors$species_name <- rownames(dfEigenvectors)
      # Subset according to the number of eigenvalues we determined previously using PVR.
      dfEigenvectors <- dfEigenvectors[, c(1:EVN, ncol(dfEigenvectors))]
      rownames(dfEigenvectors) <- NULL
      
      # If the stepwise method was selected..
    } else if(evselection == "stepwise"){
      # Apply lmforwardsequentialAICc function.
      lm1 <- lmforwardsequentialAICc(y = data[, trait], object = PEM)
      # Extract eigenvectors.
      EVN <- names(lm1$coefficients)[-1]
      # Now we need to extract the eigenvectors for each species.
      dfEigenvectors <- as.data.frame(PEM)
      dfEigenvectors$species_name <- rownames(dfEigenvectors)
      # Rearrange dataframe.
      dfEigenvectors <- dfEigenvectors[, c(EVN, "species_name")]
      # Remove row names.
      rownames(dfEigenvectors) <- NULL
    }
    # Combine dfEigenvectors into a list with parameters estimated for PEM graph (a = steepness parameter; psi = relative rate of evolution along branch lengths)
    result <- list(eigenvectors = dfEigenvectors, a = unique(PEM$a), psi = unique(PEM$psi))
    # Return result.
    return(result)
  }
  
}

DropAndMatch <- function(tree, data){
  
  # Function for matching tree tip labels to a dataframe. Also ensures the order is the same.
  # tree = phylo object with tips as species names
  # data = dataframe with column entitled "species_name"
  
  # Make sure the trait data and tree tips match.
  tree <- drop.tip(phy = tree, tip = tree$tip.label[!tree$tip.label %in% data$species_name])
  # Make sure the order matches between the phylogeny and dataframe.
  data <- data[match(tree$tip.label, data$species_name), ]
  
  # Return dataframe and tree in list format.
  return(list(phylo = tree, dataframe = data))
  
}

GapNFilter <- function(df) {
  
  # Function for removing sequences with N/gap content above a certain threshold.
  # df = dataframe with sequence data
  # threshold = threshold of N/gap content to filter by
  
  # Determine the number of positions where an *internal* N or gap is found for each sequence.
  df$internal_gapN <- str_count(df$COI, c("[N-]"))
  # Remove sequence if the number of Ns or gaps is greater than 1% (0.01) of the total length of the sequence.
  df$percentage_gapN <- df$internal_gapN/nchar(df$COI)
  df <- df[!df$percentage_gapN > 0.01, ]
  # Remove these columns as they are no longer needed.
  df <- subset(df, select = -c(internal_gapN, percentage_gapN))
  
  return(df)
  
}

PickBestLambda <- function(cvResults) {
  
  # Function for plotting the D2 values of different chronograms. D2 signifies "the influence of each observation on overall date estimates" when CV is set to TRUE for chronopl (ape package, Paradis). See ?chronopl for further detail. The lambda value that results in the lowest D2 value is considered the best lambda value.
  
  CV <- sapply(cvResults, function(x) sum(attr(x, "D2")))
  names(CV) <- 10^(-4:1)
  # Extract the lambda value that minimizes D2.
  bestLambda <- as.numeric(names(which.min(CV)))
  
  # Return the plot.
  return(bestLambda)
  
}

OutlierCheckDM <- function(mat, threshold){
  
  # Function for identifying outliers in a distance matrix based on the interquartile range (IQR) of the data. Returns a subsetted dataframe of the outliers so user can check the values for potential error.
  # data = distance matrix
  
  # Code for calculating thresholds adapted from: 
  # Author: https://stackoverflow.com/users/1312519/by0.
  # https://stackoverflow.com/questions/12866189/calculating-the-outliers-in-r.
  
  # Use the upper threshold of the IQR to detect outliers.
  lowerQuantile <- quantile(mat, na.rm = T)[2]
  upperQuantile <- quantile(mat, na.rm = T)[4]
  iqr <- upperQuantile - lowerQuantile
  upperThreshold <- (iqr * 1.5) + upperQuantile
  # Remove 0 values so that these are not considered (when a species is compared to itself - the diagonal values).
  mat[mat == 0] <- NA
  # Convert to datatable.
  dfOutliers <- as.data.table(mat, keep.rownames = T)
  # Rename the "rn" column (row names).
  setnames(dfOutliers, "rn", "species_name")
  # Identify BINs with no relatives within "typical" range of genetic divergence (i.e. all of their genetic distances are greater than 1.5 x IQR upper threshold.)
  dfOutliers <- dfOutliers[, outlier := apply(.SD, 1, function(x) all(x > upperThreshold, na.rm = T))][outlier == TRUE]
  
  # If no outliers were found..
  if(nrow(dfOutliers) == 0){
    
    print("No outliers found!")
    
    # If outliers were found..
  } else if(nrow(dfOutliers) > 0){
    
    print("Outliers found!")
    
    # Return dataframe of outliers.
    return(dfOutliers)
    
  }
  
}

ReadAndMerge <- function(..., col, newName) {
  
  # Function for reading in and merging sequence data into one dataframe. Also renames the first column to "species_name"
  
  # col = column name to merge by
  # newName = name to rename first column
  
  # Read in the alignments. Only retaining alignments with n > 100 for now.
  l_df <- ReadInAlignments(...)
  # Merge all of the dataframes using Reduce().
  dfMerge <- Reduce(function(...) merge(..., by = col, all = T), l_df)
  # Rename species_name column.
  colnames(dfMerge)[1] <- newName
  
  return(dfMerge)
  
}

ReadInAlignments <- function(filePattern, fileType) {
  
  # Function for reading in alignments and returning a list of dataframes with record title and sequence data.
  
  # filePattern = Character vector containing common pattern in fasta file names
  # fileType = should be one of "fasta" or "phy"
  
  if(fileType == "fasta") {
    # Read the fasta files in.
    fastaFiles <- list.files(pattern = filePattern)
    l_fastaFiles <- lapply(fastaFiles, readDNAStringSet)
    # Convert them into dataframes.
    l_dfFastaFiles <- lapply(l_fastaFiles, function(x) data.frame(Title = names(x), Sequence = paste(x) ))
    return(l_dfFastaFiles)
    
  } else if(fileType == "phy") {
    
    # Read the phy files in.
    phyFiles <- list.files(pattern = filePattern)
    l_dfPhyFiles <- lapply(phyFiles, read.phylip)
    # Name the list of dataframes.
    names(l_dfPhyFiles) <- phyFiles
    
    for (i in 1:length(phyFiles)) {
      
      # Take the ith dataframe in the list.
      df <- l_dfPhyFiles[[i]]
      # Take the ith file name in the vector.
      info <- phyFiles[[i]]
      
      # Name the columns in each of the dataframes according to the name of the file.
      # Append gene name to 2nd column name (column containing seq info).
      colnames(df)[2] <- paste(info, colnames(df)[2], sep = "_")
      
      # Replace dataframe.
      l_dfPhyFiles[[i]] <- df
      
    }
    
    return(l_dfPhyFiles)
    
  } else {
    
    print("Something is wrong!")
    
  }
  
}

TestPhyloSig <- function(df, colName, phylo, sig = "lambda") {
  
  # Function that estimates phylogenetic signal of the trait given a phylogenetic tree.
  # df = dataframe containing trait information. Must also contain column called "species_name" with species name info
  # colName = column name of trait of interest
  # phylo = tree of class phylo
  # sig = one of "K" (Blomberg's K) or "lambda" (Pagel's lambda)
  
  # Subset the trait from the dataframe.
  trait <- df[, colName]
  # If the trait is numeric...
  if(is.numeric(trait)) {
    # Ensure df is data.frame format.
    df <- as.data.frame(df)
    # Ensure that the order of the data matches the tree.
    df <- df[match(phylo$tip.label, df$species_name), ]
    # Name it according to the tips of the provided tree.
    names(trait) <- phylo$tip.label
    # Estimate Blomberg's K (Blomberg, 2003) or Pagel's lambda using the phylosig function from the "phytools" package. A test of significance is also performed by setting test = TRUE.
    sig <- phylosig(phylo, trait, test = TRUE, method = sig)
    # If the trait is factor (categorical)...
  } else if(is.factor(trait)) {
    # Ensure df is data.frame format.
    df <- as.data.frame(df)
    # First prune the tree so that only tips in the dataframe are present. This is necessary for the phylo.d function to work properly.
    prunedPhylo <- drop.tip(phy = phylo, tip = phylo$tip.label[!phylo$tip.label %in% df$species_name])
    # Make sure the dataframe matches the order of the tree's tip labels.
    df <- df[match(prunedPhylo$tip.label, df$species_name), ]
    # Estimate the D metric (measure of phylogenetic signal for binary traits) by Fritz and Purvis (2010) using the phylo.d function from the "caper" package.
    # Temporarily rename column as "colName" so that phylo.d function works (as it takes colName as literal name of column).
    names(df)[names(df) == colName] <- "colName"
    sig <- phylo.d(df, prunedPhylo, names.col = species_name, binvar = colName) 
  }
  # Return the estimate of phylogenetic signal.
  return(sig)
  
}

SelectCentroidSeqsDECIPHER <- function(df) {
  
  # Function for selecting centroid sequence given a dataframe of sequences.
  # Acknowlegements: This function adapted from code provided by Matt Orton (https://github.com/m-orton/R-Scripts) - also licensed under GNU GPL v3.
  # df = dataframe containing species names, sequence data, and sequence counts per species.
  
  # Set to data.table format.
  setDT(df)
  # Subset dataframe to find BINs with more than one sequence.
  dfLargeBins <- df[seq_count > 1]
  # Remove gaps from the sequences.
  df[, COI := gsub("-", "", COI)] 
  # Subset out the species with more than 1 sequence.
  dfCentroidSeqs <- df[species_name %in% dfLargeBins$species_name]
  # We also have to create another separate dataframe with species that only have one sequence, called dfSingletons.
  dfSingletons <- df[!species_name %in% dfLargeBins$species_name]
  # We then take the dfCentroidSeqs sequences and group them by groupCol.
  largeBinList <- split(dfCentroidSeqs, by = "species_name")
  # Convert all the sequences in largeBinList to DNAStringSet format for the multiple sequence alignment.
  DNAStringSetList <- lapply(largeBinList, function(x) DNAStringSet(x$COI))
  # Name DNAStringSetList using the processids.
  for (i in seq(from = 1, to = length(unique(dfCentroidSeqs$species_name)), by = 1)) {
    names(DNAStringSetList[[i]]) <- largeBinList[[i]]$processid
  }
  # Align the sequences in each BIN using DECIPHER.
  alignmentList <- lapply(DNAStringSetList, function(x) AlignTranslation(x, geneticCode = getGeneticCode("SGC1")))
  # Convert each BIN alignment to DNAbin format.
  alignmentList <- lapply(alignmentList, function(x) as.DNAbin(x))
  # Estimates the genetic distance between sequences in each BIN with the TN93 model.
  distanceMatrixList <- lapply(alignmentList, function(x) dist.dna(x, model = "TN93", as.matrix = TRUE))
  # Find the centroid sequence using the genetic distance matrix. It is the sequence in a BIN with minimum average pairwise distance to all other sequences in the BIN.
  centroidSeqs <- sapply(distanceMatrixList, function(x) which.min(rowSums(x)))
  centroidSeqs <- names(centroidSeqs)
  centroidSeqs <- gsub("^.*\\.", "", centroidSeqs)
  # Subset dfCentroidSeqs by the recordIDs of the centroid sequences.
  dfCentroidSeqs <- dfCentroidSeqs[dfCentroidSeqs$processid %in% centroidSeqs]
  # Combine the singletons and centroid sequences into a new dataframe. Now each BIN has a representative sequence.
  dfCentroidSeqs <- rbind(dfCentroidSeqs, dfSingletons)
  
  return(dfCentroidSeqs)
  
}

SelectLambdaLik <- function(phylo, minAge, maxAge, lambdas = c(-4:1)) {
  
  # Function for selecting lambda value based on log likelihood.
  # phylo = phylogenetic tree with branch lengths as expected num subs/per site
  # minAge = estimated minimum age of node
  # maxAge = estimated maximum age of node
  # lambdas = list of lambda values to consider (numbers correspond to scientific notation e.g. -4:1 would be 10^-4 ... 10^1).
  
  # Create vector to hold log likelihood values of trees.
  l_ploglik <- vector(mode = "list", length = length(lambdas))
  
  for(i in 1:length(l_ploglik)) {
    
    # Take ith lambda exponent.
    e <- lambdas[[i]]
    
    # Calibrate the tree using the provided fossil dates and the corresponding model (one of three).
    tree <- Calibrate(phylo, minAge = minAge, maxAge = maxAge, lambda = 10^e, fileName = paste(10^e, ".tre", sep = "_"))
    # Extract log likelihood value.
    likelihood <- attr(tree, "ploglik")
    # Name the likelihood according to the lambda value.
    names(likelihood) <- 10^e
    # Append calibrated tree to l_ploglik.
    l_ploglik[[i]] <- likelihood
    
  }
  
  # Which lambda resulted in tree with the largest log likelihood value?
  bestLambda <- as.numeric(names(which.max(unlist(l_ploglik))))
  # Return list of log likelihood values for each chronogram.
  return(bestLambda)
  
}

SelectLambdaValue <- function(phylo, tips = "root", minAge, maxAge, lambdas = c(-4:1)) {
  
  # Function for selecting lambda value using cross-validation.
  # Adapted from Paradis, E. (2012). Analysis of Phylogenetics and Evolution with R. Second Edition. New York, New York: Springer. Pgs 184-185, 198-199.
  # phylo = phylogenetic tree with branch lengths as expected num subs/per site
  # lambdas = list of lambda values to consider (numbers correspond to scientific notation e.g. -4:1 would be 10^-4 ... 10^1).
  # tips = character vector containing names of taxa for which fossil data info is available for calibration, default is the root
  # minAge = estimated minimum age of node
  # maxAge = estimated maximum age of node
  
  # Get the node number of the most recent common ancestor for the provided taxa.
  #nodeNum <- getMRCA(phylo, tips)
  
  
  # Create empty list to store CV results. We are considering 6 different lambda values.
  cvResults <- vector("list", length = length(lambdas))
  # Iterate through the different lambda values.
  for (i in 1:length(lambdas)) {
    
    # Take ith lambda exponent.
    e <- lambdas[[i]]
    # Apply chronopl function for cross-validation of corresponding lambda value (CV = T).
    cvResults[[i]] <- chronopl(phylo, lambda = 10^e, age.min = minAge, age.max = maxAge, node = tips, CV = TRUE)
    
  }
  
  # Return list of chronograms that have been smoothed using different lambda values.
  return(cvResults)
  
}

SelectOutgroup <- function(type, object, dfGeneName, OGtaxa, geneCol) {
  
  # Function for appending outgroup sequence data to mitochondrial or nuclear pre-centroid dataframes.
  # type = one of "amphNuc", "mammNuc", or "mt"
  # object = either a list of dataframes containing sequence data (amphNuc) or singular dataframe (mammNuc)
  # dfGeneName = name of dataframe that contains seq data for outgroup
  # OGtaxa = name of outgroup taxa
  # geneCol = name of column in dfGeneName that contains seq data.
  
  if(type == "amphNuc") {
    
    # Which dataframe in the list of dataframes contains the necessary seq data?
    geneSub <- which(names(object) == dfGeneName)
    # Extract that dataframe from the list of dataframes.
    dfGene <- object[[geneSub]]
    # We need taxonomic information for these species so let's merge with the original COI subset which has this information.
    dfGene <- merge(dfGene, dfCOI[, c("species_name", "order")], by.x = "seq.name", by.y = "species_name")
    # Subset for outgroup order.
    dfOutgroup <- subset(dfGene, order == OGtaxa)
    # Set to data.table format.
    setDT(dfOutgroup)
    # Find the number of sequences per species.
    dfOutgroup[, seq_count := length(get(geneCol)), keyby = seq.name]
    # Select species with highest number of seqs (means it is well sampled comparatively and best pick for outgroup species of that order).
    mostSampled <- dfOutgroup[which(dfOutgroup$seq_count == max(dfOutgroup$seq_count)), ]
    # If more than one, just take the first species name.
    bestOG <- unique(mostSampled$seq.name)[1]
    # Subset dataframe based on that species name.
    dfOutgroup <- subset(dfOutgroup, seq.name == bestOG)
    # Remove duplicates (multiple COI seqs but only one nuclear seq so we can do this).
    dfOutgroup <- dfOutgroup[!duplicated(dfOutgroup$seq.name)]
    # Rename seq.name to species_name.
    colnames(dfOutgroup)[1] <- "species_name"
    # Remove number of seqs col.
    dfOutgroup$seq_count <- NULL
    
    return(dfOutgroup)
    
  } else if(type == "mammNuc") {
    
    # Subset for the columns we need for matching.
    dfGene <- object[, c("species_name", geneCol)]
    # We need taxonomic information for these species so let's merge with the original COI subset which has this information.
    dfOutgroup <- merge(dfGene, dfCOI[, c("species_name", "order")], by = "species_name")
    # Remove NAs.
    dfOutgroup <- na.omit(dfOutgroup)
    # Set to data.table format.
    setDT(dfOutgroup)
    # Find the number of sequences per species.
    dfOutgroup[, seq_count := length(get(geneCol)), keyby = species_name]
    # Select species with highest number of seqs (means it is well sampled comparatively and best pick for outgroup species of that order).
    mostSampled <- dfOutgroup[which(dfOutgroup$seq_count == max(dfOutgroup$seq_count)), ]
    # If more than one, just take the first species name.
    bestOG <- unique(mostSampled$species_name)[1]
    # Subset dataframe based on that species name.
    dfOutgroup <- subset(dfOutgroup, species_name == bestOG)
    # Remove duplicates (multiple COI seqs but only one nuclear seq so we can do this).
    dfOutgroup <- dfOutgroup[!duplicated(dfOutgroup$species_name)]
    # Remove number of seqs col.
    dfOutgroup$seq_count <- NULL
    
    return(dfOutgroup)
    
  } else if(type == "mt") {
    
    # Subsetting the data by order and selecting only relevant columns.
    dfOutgroup <- subset(dfCOI, order == OGtaxa, select = c(order, species_name, processid, COI))
    # Removing NAs.
    dfOutgroup <- na.omit(dfOutgroup)
    setDT(dfOutgroup)
    # Find the number of sequences per species.
    dfOutgroup[, seq_count := length(COI), keyby = species_name]
    # Select species with highest number of seqs (means it is well sampled comparatively and best pick for outgroup species of that order).
    mostSampled <- dfOutgroup[which(dfOutgroup$seq_count == max(dfOutgroup$seq_count)), ]
    # If more than one, just take the first species name.
    bestOG <- unique(mostSampled$species_name)[1]
    # Subset dataframe based on that species name.
    dfOutgroup <- subset(dfOutgroup, species_name == bestOG)
    # Remove rownames.
    rownames(dfOutgroup) <- NULL
    # Remove number of seqs col.
    dfOutgroup$seq_count <- NULL
    
    return(dfOutgroup)
    
  } else {
    
    "ERROR!"
    
  }
}

SelectRateModel <- function(phylo, tips = "root", minAge, maxAge) {
  
  # Function for model of substitution rate variation among branches as offered by the chronos function in the ape package (Paradis). See ?chronos for further detail.
  # phylo = phylogenetic tree with branch lengths as expected num subs/per site
  # tips = character vector containing names of taxa for which fossil data info is available for calibration, default is the root
  # minAge = estimated minimum age of node
  # maxAge = estimated maximum age of node
  
  # Create a vector of the different model options.
  models <- c("correlated", "relaxed")
  
  # Get the node number of the most recent common ancestor for the provided taxa.
  #nodeNum <- getMRCA(phylo, tips)
  
  l_ploglik <- vector(mode = "list", length = 2)
  
  for(i in 1:length(l_ploglik)) {
    
    # Calibrate the tree using the provided fossil dates and the corresponding model (one of three).
    tree <- Calibrate(phylo, minAge = minAge, maxAge = maxAge, model = models[[i]], fileName = paste(models[[i]], ".tre", sep = "_"))
    # Extract log likelihood value.
    likelihood <- attr(tree, "ploglik")
    # Name the likelihood according to the model.
    names(likelihood) <- models[[i]]
    # Append calibrated tree to l_chrono.
    l_ploglik[[i]] <- likelihood
    
  }
  
  # Which model has the largest log likelihood value?
  bestModel <- names(which.max(unlist(l_ploglik)))
  # Return list of log likelihood values for each chronogram.
  return(bestModel)
  
}

SubsetAndMergeDT <- function(list, ...) {
  
  # Data.table version: Function for subsetting a list of dataframes by order name (which is contained in the species_name column pertaining to the alignment) and merging the list of dataframes.
  
  # Subsetting for order name.
  l_dfSubset <- lapply(list, function(x) SubsetByOrderDT(df = x, ...))
  # Merge all of the dataframes using Reduce().
  dfMerge <- Reduce(function(...) merge(..., by = "species_name", all = T), l_dfSubset)
  
  return(dfMerge)
  
}

SubsetByOrderDT <- function(df, string) {
  
  # Data.table version: Function for subsetting a dataframe based on string match in seq.name column (messy record title column). Also cleaning up the name a bit.
  
  # df = dataframe that contains a seq.name column (but in this case it is messy record title for a specific alignment)
  # string = string match that I want to subset by in that column
  
  # Convert to data.table
  setDT(df)
  # Subset by string using grepl() function.
  dfSubset <- df[grepl(string, seq.name)]
  # Replacing underscores with spaces.
  dfSubset[, seq.name := gsub(pattern = "_", replacement = " ", x = seq.name)]
  # Taking the last two words of the string, which correspond to the species names for which the record was obtained.
  dfSubset[, seq.name := word(seq.name, start = -2, end = -1)]
  # Remove duplicate species names (resulting from gene names in sequence record titles).
  dfSubset <- dfSubset[!duplicated(dfSubset$seq.name), ]
  # Renaming seq.name to species_name.
  colnames(dfSubset)[1] <- "species_name"
  # Add the underscores back into species_name to match trait data.
  dfSubset[, species_name := gsub(pattern = " ", replacement = "_", x = species_name)]
  # Convert back to dataframe.
  dfSubset <- as.data.frame(dfSubset)
  
  return(dfSubset)
  
}

UpdatePEM <- function(tree, alpha, rate){
  
  # Function for building phylogenetic eigenvector map (PEM) using empirically derived parameters.
  # tree = object of class "phylo". Branch lengths should be in expected # of subs/site/unit of time (ultrametric tree).
  # alpha = alpha (steepness) parameter
  # rate = relative evolution rate of the trait
  
  # Code adapted from:
  # MPSEM tutorial: G Geunard. 2019. A phylogenetic modelling tutorial using Phylogenetic Eigenvector Maps (PEM) as implemented in R package MPSEM (0.3-6).
  
  # Convert phylogenetic tree to graph-class object needed to build phylogenetic eigenvector map (PEM).
  graphPhylo <- Phylo2DirectedGraph(tree)
  # Create new PEM using derived parameter values.
  updatedPEM <- PEM.build(graphPhylo, a = alpha, psi = rate)
  # Extract the eigenvectors into a dataframe.
  dfEigenvectors <- as.data.frame(updatedPEM)
  # Create species_name column.
  dfEigenvectors$species_name <- rownames(dfEigenvectors)
  # Remove rownames.
  rownames(dfEigenvectors) <- NULL
  # Return dataframe of updated eigenvectors.
  return(dfEigenvectors)
  
}
