# SECTION 01. This scripts aggregates trait data from different sources using avian trait data as an example. The goal here is to maximize the number of complete cases, considering relevant trait data for each order. The output of this script is a number of candidate trait datasets for different orders and a dataset of original data with missing values.

# Acknowledgments. ----

# The R packages "traitdataform" and "rdataretriever" were used to obtain trait data from a number of different databases.
# Citations: Schneider FD (2020). Package 'traitdataform' version 0.6.1. https://ecologicaltraitdata.github.io/traitdataform/.
# McGlinn D, Sharma P, Harris DJ, Senyondo H, Ye H, Taylor S, Pandey A, Bansal H, Pohlman M, White R (2020). Package "rdatatriever" version 3.0.0. https://cran.r-project.org/web/packages/rdataretriever/rdataretriever.pdf.

# The Open Traits Network provides links and information for many of these datasets: https://opentraits.org/datasets.html

# Data were obtained from the following databases for this script:

# AnAge.
# Citation: de Magalhães JP, Costa J. 2009. A database of vertebrate longevity records and their relation to other life-history traits. J Evol Biol. 22:1770–1774.
# Tacutu R, Thornton D, Johnson E, Budovsky A, Barardo D, Craig T, Diana E, Lehmann G, Toren D, Wang J, et al. 2018. Human Ageing Genomic Resources: new and updated databases. Nucleic Acids Res. 46(D1):D1083–D1090. doi:10.1093/nar/gkx1042.
# URL for download: https://genomics.senescence.info/download.html#anage
# Licensed under CC BY 3.0. https://creativecommons.org/licenses/by/3.0/
# See section 2 for alterations made.

# Amniota life history database.
# Downloaded via the "traitdataform" package on March 9th, 2021.
# Citation: Myhrvold NP, Baldridge E, Chan B, Sivam D, Freeman DL, Ernest SKM. 2015. An amniote life-history database to perform comparative analyses with birds, mammals, and reptiles. Ecology. 96(11):3109.
# Myhrvold NP, Baldridge E, Chan B, Sivam D, Freeman DL, Ernest SKM. 2016. Data from: An amniote life-history database to perform comparative analyses with birds, mammals, and reptiles. https://doi.org/10.6084/m9.figshare.3563457.v1.
# Licensed under CC0 1.0. https://creativecommons.org/publicdomain/zero/1.0/
# See section 3 for alterations made.

# Vertebrate Home Range Sizes.
# Downloaded via the "rdataretriever" package on March 9th, 2021.
# Citations: Tamburello N, Côté IM, Dulvy NK. 2015. Energy and the Scaling of Animal Space Use. Am Nat. 186(2):196–211. doi:10.1086/682070.
# Tamburello N, Côté IM, Dulvy NK. 2015. Data from: Energy and the scaling of animal space use. Dryad Dataset. doi:https://doi.org/10.5061/dryad.q5j65.
# Licensed under CC0 1.0. https://creativecommons.org/publicdomain/zero/1.0/
# See section 4 for alterations made.

# Avian hand-wing index.
# Citation: Sheard C, Neate-Clegg MHC, Alioravainen N, Jones SEI, Vincent C, MacGregor HEA, Bregman TP, Claramunt S, Tobias JA. 2020. Ecological drivers of global gradients in avian dispersal inferred from wing morphology. Nat Commun. 11(1):2463. doi:10.1038/s41467-020-16313-6.
# Sheard C. 2020. catherinesheard/Global-HWI v1.1 (v1.1). Zenodo. doi:https://doi.org/10.5281/zenodo.3832215.
# Licensed under Creative Commons BY 4.0. https://creativecommons.org/licenses/by/4.0/
# See section 5 for alterations made.

# Elton traits.
# Citations: Wilman H, Belmaker J, Simpson J, de la Rosa C, Rivadeneira MM, Jetz W. 2014. EltonTraits 1.0: Species-level foraging attributes of the world’s birds and mammals. Ecology. 95(7):2027–2027. doi:10.1890/13-1917.1.
# Wilman H, Belmaker J, Simpson J, de la Rosa C, Rivadeneira MM, Jetz W. 2014. Data from: EltonTraits 1.0: Species-level foraging attributes of the world’s birds and mammals. https://esapubs.org/archive/ecol/E095/178/.
# Licensed under Creative Commons BY 4.0. https://creativecommons.org/licenses/by/4.0/
# See section 6 for alterations made.

### 1. Load libraries and functions. ----

library(data.table)
library(janitor)
library(readxl)
library(rdataretriever)
library(tidyverse)
library(traitdataform)
source("R/Functions/DataHandling_Functions.R")

### 2. Splitting AnAge. ----

# Load data from the AnAge dataset.
anage <- read.delim("AnAge.txt")
# Clean column names.
anage <- CleanColumns(anage)
# Subset by column name.
dfAnAge <- subset(anage, select = -c(hagrid, kingdom, phylum, common_name, source, sample_size, data_quality, references))

# How many species are there in each class?
table(dfAnAge$class)
# Create a species name column by pasting together the genus and species columns.
dfAnAge$species_name <- paste(dfAnAge$genus, dfAnAge$species, sep = "_")
# Let's create separate traits for maximum age depending on whether the specimen is wild or raised in captivity (See J. P. DE MAGALHÃES  J. COSTA, 2009 for rationale).
dfAnAge$max_longevity_wild <- ifelse(dfAnAge$specimen_origin == "wild", dfAnAge$maximum_longevity_yrs, NA)
dfAnAge$max_longevity_captive <- ifelse(dfAnAge$specimen_origin == "captivity", dfAnAge$maximum_longevity_yrs, NA)
# Remove columns we don't need anymore.
dfAnAge <- subset(dfAnAge, select = -c(specimen_origin, maximum_longevity_yrs))

# Split the dataframe by class.
SplitAndBreak(dfAnAge, colName = "class", name = "AnAge")

### 3. Splitting Amniota. ----

# Load data from the amniote life-history database using the pulldata() function from the "traitdataform" package.
pulldata("amniota")
# In this dataset, -999 represents NA values. Here, I am replacing all -999s with NAs.
amniota[amniota == -999] <- NA
# Clean column names.
amniota <- CleanColumns(amniota)
# Subset by column name.
dfAmniotes <- subset(amniota, select = -c(subspecies, common_name, no_sex_svl_cm, no_sex_maturity_d, no_sex_body_mass_g))

# How many species are there in each class?
table(dfAmniotes$class)
# Create a species name column by pasting together the genus and species columns.
dfAmniotes$species_name <- paste(dfAmniotes$genus, dfAmniotes$species, sep = "_")

# Split the amniota dataframe by class, as we have bird, mammal, and reptile information in this dataset.
amniotaSplit <- split(dfAmniotes, dfAmniotes$class)
names(amniotaSplit)
# list2env is a handy function for breaking down a list into its respective components and into the global environment. They are named according to the list elements.
list2env(amniotaSplit, .GlobalEnv)
# Split the dataframe by class.
SplitAndBreak(dfAmniotes, colName = "class", name = "Amni")

### 4. Splitting HomeRanges. ----

# Save the data to file.
homeranges <- fread("Tamburelloetal_HomeRangeDatabase.csv", data.table = F)
# Replace all blanks with NAs.
homeranges[homeranges == ""] <- NA
# Clean column names.
homeranges <- CleanColumns(homeranges)
# Subset by column name.
dfHomeRanges <- subset(homeranges, select = -c(taxon, common_name, primarymethod, n, log10_mass, alternative_mass_reference, log10_hra, hra_reference, realm, dimension, log10_preymass, prey_size_reference))
# Create a vector that contains taxonomy info.
taxCols <- c("class", "order", "family", "genus")
# Make the first word upper case for the taxonomy columns.
dfHomeRanges[, taxCols] <- lapply(dfHomeRanges[, taxCols], str_to_title)

# How many species are there in each class?
table(dfHomeRanges$class)
# Create a species_name column.
dfHomeRanges$species_name <- paste(dfHomeRanges$genus, dfHomeRanges$species, sep = "_")

# Split the dataframe by class.
SplitAndBreak(dfHomeRanges, colName = "class", name = "HR")

### 5. Avian hand-wind index dataset handling: Aves. ----

# Load data from the Sheard Avian hand-wing index dataset.
aves <- read_xlsx("Dataset HWI 2020-04-10.xlsx", sheet = 1)
# Replace all blanks with NAs.
aves[aves == ""] <- NA
aves[aves == "NA"] <- NA
# Clean column names.
aves <- CleanColumns(aves)
# Subset by column name.
dfAves <- subset(aves, select = -c(species_id, species_name, iucn_name, sample_size, synonym, notes))
# Rename species column.
colnames(dfAves)[2] <- "species_name"

### 6. Elton handling: Aves. ----

# Load data from the Elton 1.0 dataset.
eltonaves <- read.delim("BirdFuncDat.txt")
# Replace all blanks with NAs.
eltonaves[eltonaves == ""] <- NA
# Clean column names.
eltonaves <- CleanColumns(eltonaves)
# Subset by column name. Note: Using diet_5cat because this basically summarizes the other diet columns for us.
dfEltonAves <- subset(eltonaves, select = -c(spec_id, pass_non_pass, bl_family_latin, bl_family_english, bl_fam_sequ_id, taxo, english, diet_source, diet_inv, diet_vend, diet_vect, diet_vfish, diet_vunk, diet_scav, diet_fruit, diet_nect, diet_seed, diet_plant_o, diet_entered_by, for_strat_source, for_strat_entered_by, body_mass_comment, record_comment))
# Rename species column to species_name to keep it consistent with other dataframes and add underscores to species names.
colnames(dfEltonAves)[2] <- "species_name"
# Add underscores to species_name column for matching purposes later on.
dfEltonAves$species_name <- gsub(" ", "_", dfEltonAves$species_name)

# Here, we perform some checks on data records and simplify some of the traits for this dataset.

# TRAIT: diet_5cat.
# We will only keep "A", "B", and "C" observations since these assignments are for the species level (and not inferred from genus or family).
table(dfEltonAves$diet_certainty)
dfEltonAves$diet_5cat[dfEltonAves$diet_certainty != "A" & dfEltonAves$diet_certainty != "B" & dfEltonAves$diet_certainty != "C"] <- NA
# Remove the metadata column.
dfEltonAves$diet_certainty <- NULL

# TRAITS: for_ (foraging).
# To simplify the number of foraging categories, we are re-coding the foraging traits to "for_water", "for_ground", "for_forest", and "for_aerial". These traits currently contain the percentage of foraging behaviour in certain habitats. So, for example, if the percentage is greater than 0 for "for_strat_wataroundsurf", it would be assigned a 1 for "for_water".
# for_water
dfEltonAves$for_water <- ifelse(dfEltonAves$for_strat_watbelowsurf > 0 | dfEltonAves$for_strat_wataroundsurf > 0, 1, 0)
# for_ground
dfEltonAves$for_ground <- ifelse(dfEltonAves$for_strat_ground > 0, 1, 0)
# for_forest
dfEltonAves$for_forest <- ifelse(dfEltonAves$for_strat_understory > 0 | dfEltonAves$for_strat_midhigh > 0 | dfEltonAves$for_strat_canopy > 0, 1, 0)
# for_aerial
dfEltonAves$for_aerial <- ifelse(dfEltonAves$for_strat_aerial > 0, 1, 0)
# Remove the columns we don't need anymore.
dfEltonAves <- subset(dfEltonAves, select = -c(for_strat_watbelowsurf, for_strat_wataroundsurf, for_strat_ground, for_strat_understory, for_strat_midhigh, for_strat_canopy, for_strat_aerial))
# We are only keeping values based on species-level estimates (for_strat_spec_level == 1).
table(dfEltonAves$for_strat_spec_level)
# Create a vector of the column names for the foraging traits.
foragingTraits <- c("for_water", "for_ground", "for_forest", "for_aerial", "pelagic_specialist")
# Remove observations for foraging traits with values of 0 in the for_strat_spec_level column.
dfEltonAves[, foragingTraits] <- lapply(dfEltonAves[, foragingTraits], function(x) ifelse(dfEltonAves$for_strat_spec_level == 1, x, NA))
# Remove the metadata column.
dfEltonAves$for_strat_spec_level <- NULL

# TRAIT: body_mass_value.
# Only keeping observations with "Dunning08" in body_mass_source column (others are inferred values).
table(dfEltonAves$body_mass_source)
dfEltonAves$body_mass_value[dfEltonAves$body_mass_source != "Dunning08"] <- NA
# In addition, we are only keeping values based on species-level estimates (1) in "body_mass_spec_level" column.
table(dfEltonAves$body_mass_spec_level)
dfEltonAves$body_mass_value[dfEltonAves$body_mass_spec_level != 1] <- NA
# Remove the metadata columns.
dfEltonAves <- subset(dfEltonAves, select = -c(body_mass_source, body_mass_spec_level))

### 7. Complete-case dataset creation: Aves. ----

# Merge amniota and anage data.
sort(colnames(AmniAves))
sort(colnames(AnAgeAves))
Aves <- MergeAndUpdate(AmniAves, AnAgeAves, "order.x", "order.y")
colnames(Aves)
# Update columns.
Aves <- subset(Aves, select = -c(class.x, class.y, order.y, family.y, genus.y, species.y))
colnames(Aves)[2:5] <- c("order", "family", "genus", "species")

# Merge with homeranges data.
colnames(HRAves)
Aves <- MergeAndUpdate(Aves, HRAves, "order.x", "order.y")
colnames(Aves)
# Update columns.
Aves <- subset(Aves, select = -c(class, order.y, family.y, genus.y, species.y))
colnames(Aves)[2:5] <- c("order", "family", "genus", "species")

# Merge with Elton data.
colnames(dfEltonAves)
Aves <- MergeAndUpdate(Aves, dfEltonAves, "order", "ioc_order")
colnames(Aves)
# Update columns.
Aves$ioc_order <- NULL

# Merge with Sheard data.
colnames(dfAves)
Aves <- MergeAndUpdate(Aves, dfAves, "order.x", "order.y")
colnames(Aves)
# Update columns.
Aves$order.y <- NULL
colnames(Aves)[2] <- "order"

# Ensure there are no rows with missing species_name info or duplicated species names.
Aves <- SpeciesNameCheck(Aves)

# How many species are there in each order?
sort(table(Aves$order))
# Split the dataframe by order.
avesSplit <- split(Aves, Aves$order)
# Check the names of the list.
names(avesSplit)

# How many orders?
length(unique(Aves$order))

# Which columns contain trait data?
colnames(Aves)
# Assign the column names for trait data to a vector.
avesTraits <- colnames(Aves)
avesTraits <- avesTraits[!avesTraits %in% c("species_name", "class", "order", "family", "genus", "species")]
avesTraits
# Assign the column names that contain taxonomy data to a vector.
avesTax <- setdiff(colnames(Aves), avesTraits)
avesTax

# Identify taxa that have complete observations for at least 100 species and 3 traits. 
goodAves <- IdentifyTaxa(avesSplit, avesTraits, numSpecies = 100)
# Let's subset the list to only contain the taxa that meet these criteria.
l_goodAves <- SubsetList(avesSplit, goodAves)
# Write these dataframes to file.
mapply(fwrite, l_goodAves, paste(names(l_goodAves), ".csv", sep = ""))

# Let's create complete-case subsets for this taxa.
# Apply the CreateTraitSubset() function across our subsetted list of dataframes to create a list of complete-case dataframes.  Here, critNumber is an argument that can be tweaked to maximize sample size for each order and we do so using the CritNumberLoop() function.
CritNumberLoop(l_df = l_goodAves, traitCols = avesTraits, taxCols = avesTax)
