# Figure. Heat map of error metrics. ----

# Acknowledgements. ----
# Tutorials consulted:
# https://medium.com/analytics-vidhya/create-heatmap-in-r-using-ggplot2-d4c02dccec28
# Stackoverflow discussion consulted:
# Link: https://stackoverflow.com/questions/57814024/order-pheatmap-by-annotation. Author: GordonShumway.

# Library and function loading. ----
library(data.table)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(reshape)
library(tidyverse)
library(viridis)
source("R/Functions/DataHandling_Functions.R")
source("R/Functions/Imputation_Functions.R")

# 1. Data loading and variable assignment. ----

# Traits. ---
# Create vector of column names that contain Taxonomic information.
taxCols <- c("species_name")
# Create list of numerical traits for each order.
l_contTraits <- list(Anguilliformes = c("depth_range_deep"), Characiformes = c("length"), Clupeiformes = c("southernmost", "northernmost"), Cypriniformes = c("length"), Gadiformes = c("depth_range_deep"), Perciformes = NA, Pleuronectiformes = NA, Scorpaeniformes = c("depth_range_shallow", "depth_range_deep"), Siluriformes = c("length"), Anura = c("litter_size_max_n", "litter_size_min_n", "body_size_mm"), Caudata = c("body_size_mm"), Apodiformes = c("temp_range", "annual_temp", "precip_range", "annual_precip"), Charadriiformes = NA, Galliformes = c("latitude", "precip_range", "annual_precip", "temp_range", "annual_temp"), Passeriformes = c("temp_range", "annual_temp", "precip_range", "annual_precip"), Piciformes = c("latitude", "precip_range", "annual_precip", "temp_range", "annual_temp"), Psittaciformes = c("latitude", "temp_range", "annual_temp", "precip_range", "annual_precip"), Carnivora = NA, Chiroptera = c("mass_g.x", "adult_forearm_len_mm", "pet_mean_mm", "aet_mean_mm", "temp_mean_01deg_c", "precip_mean_mm"), Rodentia = c("mass_g.x", "temp_mean_01deg_c", "precip_mean_mm", "gr_area_km2"), Squamata = c("latitude_centroid_from_roll_et_al_2017", "smallest_clutch", "largest_clutch", "female_svl" ))
# Create list of categorical traits for each order.
l_catTraits <- list(Anguilliformes = c("resilience"), Characiformes = NA, Clupeiformes = c("body_shape_ii", "striking_features", "posof_mouth"), Cypriniformes = c("operculum_present", "body_shape_i"), Gadiformes = NA, Perciformes = c("notched", "operculum_present", "body_shape_i", "resilience"), Pleuronectiformes = c("resilience"), Scorpaeniformes = NA, Siluriformes = c("resilience"), Anura = c("noc", "diu", "arb", "aqu", "fos"), Caudata = c("ter", "fos", "breeding"), Apodiformes = c("for_aerial", "migration_1", "habitat"), Charadriiformes = c("territoriality", "migration_3"), Galliformes = NA, Passeriformes = c("for_forest", "for_ground", "habitat", "territoriality"), Piciformes = c("for_ground"), Psittaciformes = NA, Carnivora = c("activity_diurnal", "activity_crepuscular", "activity_nocturnal"), Chiroptera = c("eats_plants", "eats_animals", "for_strat_value"), Rodentia = c("eats_animals", "activity_crepuscular"), Squamata = c("activity_time"))
# Combine l_contTraits and l_catTraits into one list.
l_traits <- mapply(function(x, y){
  z <- c(x, y)
  z <- na.omit(z)
}, x = l_contTraits, y = l_catTraits)
# Create vector of order names.
orderNames <- names(l_catTraits)

# Read in the error rate files as dataframes.
# Identify file names.
errorFiles <- list.files(path = "Results/ErrorRates/", pattern = "_ErrorRates.csv")
# Read in error files.
l_dfErrors <- lapply(paste("Results/ErrorRates/", errorFiles, sep = ""), fread, data.table = F)
# Name according to file.
names(l_dfErrors) <- errorFiles

# 2. Error rate heat maps (COI). ----

# Create list to hold error rate dataframes for each order, based on numerical and categorical traits.
l_dfAllCont <- CreateNamedList(listLength = length(orderNames), elementNames = orderNames)
l_dfAllCat <- CreateNamedList(listLength = length(orderNames), elementNames = orderNames)

# For every order..
for(o in 1:length(orderNames)){
  
  # Extract order.
  orderName <- orderNames[[o]]
  # Extract trait names for the Taxon.
  traits <- l_traits[[o]]
  # Extract numerical traits for the order.
  contTraits <- l_contTraits[[o]]
  # Extract numerical traits for the order.
  catTraits <- l_catTraits[[o]]
  
  # Get error rate results for the order.
  index <- grep(orderName, x = names(l_dfErrors))
  l_dfOrderErrors <- l_dfErrors[index]
  
  # Make a list to hold the dataframes for each trait.
  l_dfTraitError <- CreateNamedList(length(traits), elementNames = traits)
  
  # For each trait...
  for(t in 1:length(traits)) {
    
    # Take the tth trait.
    trait <- traits[[t]]
    
    # Index to get the corresponding dfErrors for the trait.
    index <- grep(trait, names(l_dfOrderErrors))
    # Subset l_dfErrors.
    l_dfSubset <- l_dfOrderErrors[index]
    names(l_dfSubset)
    
    for(a in 1:length(l_dfSubset)) {
      
      # Take the ath l_dfSubset.
      dfSubset <- l_dfSubset[[a]]
      # Add new column name called "trait".
      dfSubset$trait <- paste(orderName, trait)
      # Remove extraneous cols (V1, missingness_type, SE) for now.
      dfSubset <- dfSubset [, -c(1:2, 4)]
      # Rename column 1.
      colnames(dfSubset)[1] <- "error_rate"
      
      # Extract title for df.
      title <- names(l_dfSubset)[[a]]
      # Use gsub to parse out method name from title:
      method <- gsub(paste(orderName, "_", trait, "_", sep = ""), "", title) ## Order name and trait
      method <- gsub("_ErrorRates.csv", "", method) ## File extension
      # Add new column name called "method".
      dfSubset$method <- method
      
      # Replace original dfSubset.
      l_dfSubset[[a]] <- dfSubset
      
    }
    
    # Rbind all dfSubsets together.
    dfTraitError <- do.call(rbind, l_dfSubset)
    # Append to l_dfTraitError.
    l_dfTraitError[[t]] <- dfTraitError
    
  }
  
  # Rbind all dfTraitErrors together.
  dfAllOrder <- do.call(rbind, l_dfTraitError)
  rownames(dfAllOrder) <- NULL
  
  # Split dfAllOrder based on numerical and categorical traits.
  # Identify rows with numerical trait data.
  index <- sapply(contTraits, grep, x = dfAllOrder$trait)
  # Subset based on index.
  dfAllCont <- dfAllOrder[index, ]
  # Identify rows with categorical trait data.
  index <- sapply(catTraits, grep, x = dfAllOrder$trait)
  # Subset based on index.
  dfAllCat <- dfAllOrder[index, ]
  
  # Append to l_dfAllCont and l_dfAllCat.
  l_dfAllCont[[o]] <- dfAllCont
  l_dfAllCat[[o]] <- dfAllCat
  
}

# Bind dataframes in l_dfAllCont and l_dfAllCat.
dfFinalCont <- do.call(rbind, l_dfAllCont)
dfFinalCat <- do.call(rbind, l_dfAllCat)
# Remove row names.
rownames(dfFinalCont) <- NULL
rownames(dfFinalCat) <- NULL
# Get rid of NAs.
dfFinalCont <- na.omit(dfFinalCont)

dfFinalCat <- na.omit(dfFinalCat)
# Replace MeanMode with Mean and Mode for numerical traits.
dfFinalCont[dfFinalCont == "MeanMode"] <- "Mean"
dfFinalCat[dfFinalCat == "MeanMode"] <- "Mode"
# Pivot dataframe so there is a single row per trait.
dfContPivot <- pivot_wider(dfFinalCont, names_from = method, values_from = error_rate)
dfCatPivot <- pivot_wider(dfFinalCat, names_from = method, values_from = error_rate)
# To get COI cols only.
dfContPivotCOI <- dfContPivot[, 1:8]
dfCatPivotCOI <- dfCatPivot[, 1:8]

# TODO: Update trait names in error rate files.

# Convert error rate values to matrix format.
matAllCont <- as.matrix(dfContPivotCOI[, c(2:8)])
matAllCat <- as.matrix(dfCatPivotCOI[, c(2:8)])

# Create dataframe to hold Taxon and trait information (used in heatmap for grouping purposes).
dfTaxonCont <- data.frame("trait" = dfContPivotCOI$trait)
dfTaxonCat <- data.frame("trait" = dfCatPivotCOI$trait)
# Create Taxon column (class).
dfTaxonCont$trait
dfTaxonCont$Taxon <- c(rep("Actinopterygii", 9), rep("Amphibia", 4), rep("Aves", 23),  rep("Mammalia", 10), rep("Reptilia", 4))
dfTaxonCat$Taxon <- c(rep("Actinopterygii", 12), rep("Amphibia", 8), rep("Aves", 10), rep("Mammalia", 8), "Reptilia")
# Add row names based on trait.
rownames(dfTaxonCont) <- dfContPivotCOI$trait
rownames(dfTaxonCat) <- dfCatPivotCOI$trait

# Create type column to include info on biological vs. environmental.
dfTaxonCont$Type <- factor(c("Environmental", "Biological", rep("Environmental", 2), "Biological", rep("Environmental", 3), rep("Biological", 5),  rep("Environmental", 23), rep("Biological", 2), rep("Environmental", 4), "Biological", rep("Environmental", 2), "Biological", "Environmental", rep("Biological", 3)), levels = c("Biological", "Environmental"), ordered = T)
dfTaxonCat$Type <- c(rep("Biological", 14), rep("Environmental", 5), rep("Biological", 3), "Environmental", rep("Biological", 4), "Environmental", rep("Biological", 11))
# Remove trait column.
dfTaxonCont$trait <- NULL
dfTaxonCat$trait <- NULL
# Match row names in matrix to names in dfTaxon.
rownames(matAllCont) <- rownames(dfTaxonCont)
rownames(matAllCat) <- rownames(dfTaxonCat)

# Convert to factors and ordered = T so we can add order to the annotations (by Type and then Taxon).
dfTaxonCont$Type <- factor(dfTaxonCont$Type, ordered = T)
dfTaxonCont$Taxon <- factor(dfTaxonCont$Taxon, ordered = T)
dfTaxonCont <- dfTaxonCont[order(dfTaxonCont$Type, dfTaxonCont$Taxon, rownames(dfTaxonCont)), ]
dfTaxonCat <- dfTaxonCat[order(dfTaxonCat$Type, dfTaxonCat$Taxon, rownames(dfTaxonCat)), ]

# Reorder rows based on annotation order.
matAllCont <- matAllCont[rownames(dfTaxonCont), ]
matAllCat <- matAllCat[rownames(dfTaxonCat), ]
# Reorder columns so we can more readily see the patterns in the heatmaps.
colnames(matAllCont)
matAllCont <- matAllCont[, c("KNN", "MICE", "RF", "Mean", "KNN_COI", "MICE_COI", "RF_COI")]
colnames(matAllCat)
matAllCat <- matAllCat[, c("KNN", "MICE", "RF", "Mode", "KNN_COI", "MICE_COI", "RF_COI")]

# Create colour annotations for numerical traits.
ann_colors_cont = list(Type = c(Environmental = "cornflowerblue", Biological = "indianred1"), Taxon = c(Actinopterygii = "lightskyblue2", Amphibia = "lightgreen", Aves = "lightgoldenrod", Mammalia = "coral", Reptilia = "rosybrown2"))
# Create colour annotations for categorical traits.
ann_colors_cat = list(Type = c(Environmental = "cornflowerblue", Biological = "indianred1"), Taxon = c(Actinopterygii = "lightskyblue2", Amphibia = "lightgreen", Aves = "lightgoldenrod", Mammalia = "coral", Reptilia = "rosybrown2"))

# Colours. 
COLS <- colorRampPalette(brewer.pal(9, "GnBu"))(10)
COLS2 <- rev(COLS)
COLS3 <- c(COLS2, "white")

# Assign number of rows per group so we can add white space between trait Types (easier to read):
dfTaxonCont$Taxon
gapsCont <- c(3, 7, 11, 14, 20, 43, 49)
dfTaxonCat$Taxon
gapsCat <- c(12, 15, 23, 31, 32, 37)
# Assign number of cols per group so we can add white space between method Types (easier to read):
colnames(matAllCont)
gapsCOL <- c(3, 4)

# Plot numerical trait results.
pheatmap(matAllCont, color = COLS3, breaks = c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.5), cluster_rows = F, scale = "none", cluster_cols = F, gaps_row = gapsCont, gaps_col = gapsCOL, annotation_row = dfTaxonCont, fontsize = 8, cellwidth = 75, angle_col = "45", annotation_colors = ann_colors_cont, main = "Numerical traits: NMSE", filename = "COI_NUM.png")

# Plot categorical trait results.
pheatmap(matAllCat, color = colorRampPalette(brewer.pal(8, "RdPu"))(25), cluster_rows = F, scale = "none", gaps_col = gapsCOL, cluster_cols = F, gaps_row = gapsCat, annotation_row = dfTaxonCat, fontsize = 8, legend = T, cellwidth = 75, angle_col = "45",  annotation_colors = ann_colors_cat, main = "Categorical traits: ARI", filename = "COI_CAT.png")


# 3. Error rate heat maps (RAG1). ----

# To get RAG1 cols only.
dfContPivotRAG1 <- dfContPivot[, 1:11]
dfCatPivotRAG1 <- dfCatPivot[, 1:11]
# Remove NAs so we are only looking at traits imputed using this gene.
dfContPivotRAG1 <- na.omit(dfContPivotRAG1)
dfCatPivotRAG1 <- na.omit(dfCatPivotRAG1)

# Convert error rate values to matrix format.
matAllContRAG1 <- as.matrix(dfContPivotRAG1[, c(2:11)])
matAllCatRAG1 <- as.matrix(dfCatPivotRAG1[, c(2:11)])
# Assign row names based on trait to matrices.
rownames(matAllContRAG1) <- dfContPivotRAG1$trait
rownames(matAllCatRAG1) <-  dfCatPivotRAG1$trait

# Create dataframe to hold Taxon and trait information (used in heatmap for grouping purposes).
dfTaxonContRAG1 <- data.frame("trait" = dfContPivotRAG1$trait)
dfTaxonCatRAG1 <- data.frame("trait" = dfCatPivotRAG1$trait)
# Add row names based on trait.
rownames(dfTaxonContRAG1) <- dfContPivotRAG1$trait
rownames(dfTaxonCatRAG1) <- dfCatPivotRAG1$trait

# Merge with dfTaxonCont and dfTaxonCat to get Taxonomy and trait Type info.
dfTaxonContRAG1 <- merge(dfTaxonContRAG1, dfTaxonCont, by = "row.names")
dfTaxonCatRAG1 <-  merge(dfTaxonCatRAG1, dfTaxonCat, by = "row.names")

# Replace row names with Row.names column (created when merging).
rownames(dfTaxonContRAG1) <- dfTaxonContRAG1$Row.names
rownames(dfTaxonCatRAG1) <- dfTaxonCatRAG1$Row.names
# Remove trait and Row.names columns.
dfTaxonContRAG1$trait <- NULL
dfTaxonCatRAG1$trait <- NULL
dfTaxonContRAG1$Row.names <- NULL
dfTaxonCatRAG1$Row.names <- NULL

# Convert to factors and ordered = T so we can add order to the annotations (by Type and then Taxon). ***
dfTaxonContRAG1$Type <- factor(dfTaxonContRAG1$Type, ordered = T)
dfTaxonContRAG1$Taxon <- factor(dfTaxonContRAG1$Taxon, ordered = T)
dfTaxonContRAG1 <- dfTaxonContRAG1[order(dfTaxonContRAG1$Type, dfTaxonContRAG1$Taxon, rownames(dfTaxonContRAG1)), ]
dfTaxonCatRAG1$Type <- factor(dfTaxonCatRAG1$Type, ordered = T)
dfTaxonCatRAG1$Taxon <- factor(dfTaxonCatRAG1$Taxon, ordered = T)
dfTaxonCatRAG1 <- dfTaxonCatRAG1[order(dfTaxonCatRAG1$Type, dfTaxonCatRAG1$Taxon, rownames(dfTaxonCatRAG1)), ]

# Reorder rows based on annotation order. ***
matAllContRAG1 <- matAllContRAG1[rownames(dfTaxonContRAG1), ]
matAllCatRAG1 <- matAllCatRAG1[rownames(dfTaxonCatRAG1), ]
# Reorder columns so we can more readily see the patterns in the heatmaps.
colnames(matAllContRAG1)
matAllContRAG1 <- matAllContRAG1[, c("KNN", "MICE", "RF", "Mean", "KNN_COI", "MICE_COI", "RF_COI", "KNN_RAG1", "MICE_RAG1", "RF_RAG1")]
colnames(matAllCatRAG1)
matAllCatRAG1 <- matAllCatRAG1[, c("KNN", "MICE", "RF", "Mode", "KNN_COI", "MICE_COI", "RF_COI", "KNN_RAG1", "MICE_RAG1", "RF_RAG1")]

# Annotations for numerical traits.
unique(dfTaxonContRAG1$Taxon)
table(dfTaxonContRAG1$Type)
ann_colors_cont = list(
  Taxon = c(Actinopterygii = "lightskyblue2", Aves = "lightgoldenrod", Reptilia = "rosybrown2"), Type = c(Environmental = "cornflowerblue", Biological = "indianred1"))

# Annotations for categorical traits.
unique(dfTaxonCatRAG1$Taxon)
table(dfTaxonCatRAG1$Type)
ann_colors_cat = list(
  Taxon = c(Actinopterygii = "lightskyblue2", Aves = "lightgoldenrod",  Reptilia = "rosybrown2"), Type = c(Environmental = "cornflowerblue", Biological = "indianred1"))

# Colours.
COLS <- colorRampPalette(brewer.pal(8, "GnBu"))(9)
COLS2 <- rev(COLS)
COLS3 <- c(COLS2)

# Assign number of rows per group so we can add white space between trait Types (easier to read):
dfTaxonContRAG1$Taxon
gapsCont <- c(1, 4, 8)
dfTaxonCatRAG1$Taxon
gapsCat <- c(6, 9, 10)

# Assign number of cols per group so we can add white space between method Types (easier to read):
colnames(matAllContRAG1)
gapsCOL <- c(3, 4, 7)

pheatmap(matAllContRAG1, color = COLS3, breaks = c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.5), cluster_rows = F, scale = "none", cluster_cols = F, gaps_row = gapsCont, gaps_col = gapsCOL, annotation_row = dfTaxonContRAG1, fontsize = 12, cellwidth = 60, angle_col = "45", annotation_colors = ann_colors_cont, main = "Numerical traits: NMSE", filename = "RAG1_NUMERICAL.png", width = 14)

pheatmap(matAllCatRAG1, color = colorRampPalette(brewer.pal(8, "RdPu"))(25), cluster_rows = F, scale = "none", cluster_cols = F, gaps_row = gapsCat, gaps_col = gapsCOL, annotation_row = dfTaxonCatRAG1, fontsize = 12, cutree_cols = 5, legend = T, cellwidth = 60,  angle_col = "45",  annotation_colors = ann_colors_cat, main = "Categorical traits: ARI", filename = "RAG1_CATEGORICAL.png", width = 14)
