# Script for obtaining correlation coefficients and plotting violin and stacked bar plots for complete-case, original, and imputed datasets.

# Acknowledgments:
# https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
# https://drsimonj.svbtle.com/pretty-scatter-plots-with-ggplot2
# http://www.sthda.com/english/wiki/ggplot2-colors-how-to-change-colors-automatically-and-manually
# https://ggplot2.tidyverse.org/reference/aes_colour_fill_alpha.html

# For modifying correlation matrix:
# Author: akrun. Link: https://stackoverflow.com/questions/56193208/sort-for-top-matrix-correlations-and-remove-reverse-duplicates-without-apply

# Violin plots:
# https://r-graph-gallery.com/violin_and_boxplot_ggplot2.html
# https://rpkgs.datanovia.com/ggpubr/reference/ggviolin.html

# Grouped barplots:
# https://r-graph-gallery.com/48-grouped-barplot-with-ggplot2.html

# Load libraries and functions. ----
library(colorspace)
library(corrplot)
library(data.table)
library(Hmisc)
library(psych)
library(fields)
library(MASS)
library(wesanderson)
library(ggridges)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(hrbrthemes)
library(viridis)
source("R/Functions/DataHandling_Functions.R")
source("R/Functions/Imputation_Functions.R")
source("R/Functions/Phylo_Functions.R")
source("R/Functions/Simpute_Functions.R")

# Load data. ----
# Read in dataframe.
dfAll <- fread("ApodiformesFinalDataset.csv", data.table = F)
# Ensure blanks are NAs.
dfAll[dfAll == ""] <- NA

# Updating dataset to remove species in dfCC that had missing values in any trait.
dfNC <- read.csv("Apodiformes_NC.csv")
# Removing rows with any missing values so it is truly CC for plotting comparison.
dfUpdate <- na.omit(dfNC)
# Create vector of species names to keep in dfCC.
goodSpecies <- dfUpdate$species_name
# Get index of species names not in goodSpecies.
badSpecies <- setdiff(dfNC$species_name, goodSpecies)
# Identify cc columns in dfAll. ++
ccCols <- grep("_cc$", colnames(dfAll))
# Update these species with NAs in dfAll.
dfAll[dfAll$species_name %in% badSpecies, ccCols] <- NA

# Variable assignment. ---
# Assign order name.
orderName <- "Apodiformes"
# Look at column names of dfCC.
colnames(dfAll)
# Remove misc col.
dfAll$V1 <- NULL
# Create vector of column names that contain taxonomic information.
taxCols <- c("species_name")
# What column is species name?
specCol <- grep("species_name", colnames(dfAll))

# Split up according to complete-case, original, and imputed.
# Complete-case:
dfCC <- dfAll[, c(grep("_cc$", colnames(dfAll)), specCol)]
colnames(dfCC)
# Fix names.
colnames(dfCC)[1:(ncol(dfCC)-1)] <- c("annual_precip", "annual_temp", "aerial_foraging", "habitat", "migration", "precip_range", "temp_range")
# Original:
dfRaw <- dfAll[, c(grep("_raw$", colnames(dfAll)), specCol)]
# Fix names.
colnames(dfRaw)
colnames(dfRaw)[1:(ncol(dfCC)-1)] <- c("annual_precip", "annual_temp", "aerial_foraging", "habitat", "migration", "precip_range", "temp_range")
# Imputed:
dfImp <- dfAll[, c(grep("_raw_imp$", colnames(dfAll)), specCol)]
# Fix names.
colnames(dfImp)
colnames(dfImp)[1:(ncol(dfCC)-1)] <- c("annual_precip", "annual_temp", "aerial_foraging", "habitat", "migration", "precip_range", "temp_range")

# Combine dataframes into list.
l_dfTraits <- list(dfCC, dfRaw, dfImp)
# Name accordingly.
names(l_dfTraits) <- c("cc", "original", "imputed")

# Plots. ----

# Create list to hold correlation coefficient information for each dataset.
l_dfCorr <- CreateNamedList(listLength = length(l_dfTraits), elementNames = names(l_dfTraits))

# For every dataframe..
for(i in 1:length(l_dfTraits)){
  
  # Take the ith dataframe.
  dfTraits <- l_dfTraits[[i]]
  # Get name of dataframe.
  nameDf <- names(l_dfTraits)[i]
  # Extract trait names from dfTraits.
  traits <- setdiff(colnames(dfTraits), taxCols)
  # Apply BreakIntoTypes() function to identify which traits are numerical and which are categorical.
  l_traits <- BreakIntoTypes(dfTraits, traits)
  # Extract numerical traits.
  contTraits <- l_traits[[1]]
  # Extract categorical traits.
  catTraits <- l_traits[[2]]
  # Ensure correct format for catTraits.
  if(length(catTraits) > 1) {
    # If more than one cat trait, use lapply.
    dfTraits[, catTraits] <- lapply(dfTraits[, catTraits], as.factor) 
  } else {
    dfTraits[, catTraits] <- as.factor(dfTraits[, catTraits]) 
  }
  # Log-transform contTraits to improve visualization.
  dfTraits[, contTraits] <- LogTransform(dfTraits[, contTraits], contTraits)
  
  # Calculate correlation coefficients (pearson). ---
  l_corrs <- rcorr(as.matrix(dfTraits[, contTraits]))
  # Reformat all dataframes within l_corrs so we can readily merge dataframes together.
  l_corrs <- lapply(l_corrs, function(x) {
    # Modify format of dataframe amendable to plotting.
    x <- as.data.frame(x)
    # Replace lower triangle of dataset with NAs.
    x[lower.tri(x, diag = TRUE)] <- NA
    # Create a column to hold trait-trait info.
    x$row <- rownames(x)
    # Pivot dataframe.
    x <- pivot_longer(x, cols = !row)
    # Remove NAs.
    x <- na.omit(x)
    # Add trait_combo column.
    x$trait_combo <- paste(x$row, x$name, sep = "/")
    # Remove columns we don't need anymore.
    x[, c("row", "name")] <- NULL
    
    return(x)
  })
  
  # Merge all dataframes in l_corrs together.
  dfCorr <- Reduce(function(...) merge(..., by = c("trait_combo")), l_corrs)
  # Rename columns.
  colnames(dfCorr)[2:4] <- c("r", "n", "p")
  # Ensure dataframe format.
  dfCorr <- as.data.frame(dfCorr)
  
  # Append to l_dfCorr.
  l_dfCorr[[i]] <- dfCorr
  
}

# Combine dataframes in l_dfCorr.
dfCorrAll <- do.call("cbind", l_dfCorr)
# Remove row names.
rownames(dfCorrAll) <- NULL
# Add order_name column.
dfCorrAll$order_name <- orderName
# Remove redundant cols.
dfCorrAll$original.trait_combo <- NULL
dfCorrAll$imputed.trait_combo <- NULL
# Rename first column.
colnames(dfCorrAll)[1] <- "trait_combo"

# Write dfCorrAll to file.
fwrite(dfCorrAll, paste(orderName, "dfCorrAll.csv", sep = "_"))

# Violin plots. ----

# Separate dfAll by numerical and categorical traits.
colnames(dfAll)
dfAllCont <- dfAll[, c(1:6, 16:18, 20:22)]
dfAllCat <- dfAll[, c(7:15)]

# Pivot dataframes.
dfAllContPivot <- pivot_longer(dfAllCont, colnames(dfAllCont))
# Create dataset column.
dfAllContPivot$dataset <- case_when(grepl("_cc$", dfAllContPivot$name) ~ "cc",
                                    grepl("_raw$", dfAllContPivot$name) ~ "original",
                                    grepl("_raw_imp$", dfAllContPivot$name) ~ "imputed")
# Remove NAs.
dfAllContPivot <- na.omit(dfAllContPivot)
# Clean up name column to remove dataset IDs.
dfAllContPivot$name <- gsub(pattern = "_cc$|_raw$|_raw_imp$", replacement = "", x = dfAllContPivot$name)
# Update trait names.
contTraits <- unique(dfAllContPivot$name)
# Split up dfAllContPivot by trait.
l_dfAllContPivot <- split(dfAllContPivot, dfAllContPivot$name)
# Set colour palette.
cols <- c("#F3E79A", "#ED7C97", "#9F7FCD")

# For every trait in l_dfAllContPivot..
for(t in 1:length(l_dfAllContPivot)){
  
  # Take the tth dataframe.
  dfAllCont <- l_dfAllContPivot[[t]]
  # Get tth trait.
  trait <- unique(dfAllCont$name)

  # if(trait == "mass_g.x"){
  #   # Remove the huge capybera! to improve visualization.
  #   #dfAllCont <- dfAllCont[!dfAllCont$value > 3000.0000, ]
  #   dfAllCont$value <- log(dfAllCont$value)
  # }
  
  # Modify label to add sample sizes to dataset variables.
  sampleSizes <- plyr::count(dfAllCont$dataset)
  # Sort by decreasing so imputed is first and cc last.
  sampleSizes <- sort(sampleSizes$freq, decreasing = T)
  # Create vector to append to dataset labels.
  sampleSizesMOD <- paste(c("Imputed \n(n = ", "Original \n(n = ", "Complete-case \n(n = "), sampleSizes, ")", sep = "")
  
  # Append sample sizes to dataset names so we can include those in the plot.
  dfAllCont$dataset <- case_when(grepl("cc", dfAllCont$dataset) ~ sampleSizesMOD[[3]],
                                 grepl("original", dfAllCont$dataset) ~ sampleSizesMOD[[2]],
                                 grepl("imputed", dfAllCont$dataset) ~ sampleSizesMOD[[1]])
  
  # Create plot title.
  plotTitle <- paste(orderName, trait, "ViolinPlot.tiff", sep = "_")

  tiff(plotTitle, units = "in", width = 10, height = 8, res = 600)
  plot <- ggviolin(dfAllCont, x = "dataset", y = "value", fill = "dataset", orientation = "horiz", add = c("boxplot"), palette = cols, color = NA, width = 0.9, order = sampleSizesMOD, add.params = list(fill = "white", alpha = 0.7, color = "black"), title = trait, trim = F, font.xtickslab = c(21, "bold", "black"), font.ytickslab = c(20, "bold", "black"), xlab = "", ylab = "", legend = "none", ggtheme = theme_minimal())
  print(plot)
  dev.off()
  
}

# Barplots for categorical traits. ----

# Convert to factor type.
dfAllCat[,] <- lapply(dfAllCat[,], as.factor)
# Pivot dataframes.
dfAllCatPivot <- pivot_longer(dfAllCat, colnames(dfAllCat))
# Create dataset column.
dfAllCatPivot$dataset <- case_when(grepl("_cc$", dfAllCatPivot$name) ~ "cc",
                                    grepl("_raw$", dfAllCatPivot$name) ~ "original",
                                    grepl("_raw_imp$", dfAllCatPivot$name) ~ "imputed")
# Remove NAs.
dfAllCatPivot <- na.omit(dfAllCatPivot)
# Clean up name column to remove dataset IDs.
dfAllCatPivot$name <- gsub(pattern = "_cc$|_raw$|_raw_imp$", replacement = "", x = dfAllCatPivot$name)
# Update trait names.
catTraits <- unique(dfAllCatPivot$name)
# Split up dfAllCatPivot by trait.
l_dfAllCatPivot <- split(dfAllCatPivot, dfAllCatPivot$name)

cols <- c("#F3E79A", "#9F7FCD", "#F7A086", "#CF63A6")

# For every trait in l_dfAllCatPivot..
for(t in 1:length(l_dfAllCatPivot)){
  
  # Take the tth dataframe.
  dfPivot <- l_dfAllCatPivot[[t]]
  # Get tth trait.
  trait <- unique(dfPivot$name)
  # Group by name and value.
  dfCount <- group_by(na.omit(dfPivot), dataset, value)
  # Get counts.
  dfCount <- dplyr::summarise(dfCount, count = n())
  # Create proportion column.
  setDT(dfCount)
  dfCount[, groupN := sum(count), by = dataset]
  dfCount[, prop := count/groupN]
  
  # Modify label to add sample sizes to dataset variables.
  sampleSizes <- plyr::count(dfPivot$dataset)
  # Sort by decreasing so imputed is first and cc last.
  sampleSizes <- sort(sampleSizes$freq, decreasing = T)
  # Create vector to append to dataset labels.
  sampleSizesMOD <- paste(c("Imputed \n(n = ", "Original \n(n = ", "Complete-case \n(n = "), sampleSizes, ")", sep = "")
  
  # Append sample sizes to dataset names so we can include those in the plot.
  dfCount$dataset <- case_when(grepl("cc", dfCount$dataset) ~ sampleSizesMOD[[3]],
                                 grepl("original", dfCount$dataset) ~ sampleSizesMOD[[2]],
                                 grepl("imputed", dfCount$dataset) ~ sampleSizesMOD[[1]])

  # Create plot title.
  plotTitle <- paste(orderName, trait, "StackedBarPlot.tiff", sep = "_")
  
  # Reorder dataset according to how we want it plotted.
  dfCount$dataset <- factor(dfCount$dataset, levels = c(sampleSizesMOD[[3]], sampleSizesMOD[[2]], sampleSizesMOD[[1]]))
  # Ensure trait is factor type.
  dfCount$value <- as.factor(dfCount$value)
  
  # Plot.
  tiff(plotTitle, units = "in", width = 13, height = 14, res = 600)
  
  plot <- ggplot(data = dfCount, mapping = aes(x = dataset, y = prop, fill = value)) + 
    geom_bar(stat="identity", width = 0.95) + 
    theme_light() +
    scale_fill_manual(values = cols) +
    geom_text(aes(label = scales::percent(prop, accuracy = 0.1)), vjust = 0.4, size = 11, position = position_stack(vjust = 0.5), fontface = "bold") +
    theme(axis.text.y=element_text(size = 31, face = "bold"), axis.text.x=element_text(size = 29, face = "bold"), legend.text = element_text(face = 'bold', size = 25)) +
    labs(title = trait,
         y = "", x = "") 
  
  print(plot)
  dev.off()
  
}
