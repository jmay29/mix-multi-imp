# Script for creating scatterplots between correlation coefficients.

# Acknowlegments:

# Consulted this tutorial for creating a pretty scatterplot in R: 
# https://drsimonj.svbtle.com/pretty-scatter-plots-with-ggplot2
# Author: Simon Jackson

# Load libraries and functions. ----
library(colorspace)
library(corrplot)
library(data.table)
library(Hmisc)
library(psych)
library(fields)
library(MASS)
library(ggpubr)
library(wesanderson)
library(ggridges)
library(ggpmisc)
library(tidyverse)
source("R/Functions/DataHandling_Functions.R")
source("R/Functions/Imputation_Functions.R")
source("R/Functions/Phylo_Functions.R")
source("R/Functions/Simpute_Functions.R")

# Load data. ----
# Read in dataframes containing pairwise correlation coefficients.
corrFiles <- list.files(pattern = "_dfCorrAll.csv", path = "Results/FinalImputations/")
l_dfPlot <- lapply(corrFiles, function(x) fread(paste("Results/FinalImputations/", x, sep = ""), data.table = F))
# Bind the dataframes.
dfPlot <- rbindlist(l_dfPlot)

# First examine which relationships have too similar sample sizes between either original and cc or original and imputed (aka it doesn't tell us much).
dfPlot <- dfPlot[-c(12, 21, 26, 17, 24, 27), ]
# Set thresholds for scale.
lowerThreshold <- -0.85
upperThreshold <- 1

# Scatterplots. ----

pal <- wes_palette("GrandBudapest2", type = c("continuous"))
colnames(dfPlot)

# Original vs. complete-case.
tiff("original_cc_CorrScatter.tiff", units = "in", width = 10, height = 7, res = 700)

plot <- ggplot(dfPlot, aes(x = original.r, cc.r)) +
  geom_point(aes(shape = order_name, color = order_name), size = 10, show.legend = T, alpha = 0.7) +
  theme_minimal() +
  scale_color_manual(values = pal) +
  scale_shape_manual(values = c(17, 15, 16, 18)) +
  xlab("Original dataset (r)") +
  ylab("Complete-case dataset (r)") +
  ylim(lowerThreshold, upperThreshold) +
  geom_smooth(method = lm, se = F, show.legend = F, col = "skyblue") +
  theme(axis.title.x = element_text(vjust = 0.6, size = 20, face = "bold"), axis.title.y = element_text(vjust = 0.6, size = 20, face = "bold"), axis.text.y = element_text(vjust = 0.6, size = 18, face = "bold"), axis.text.x = element_text(vjust = 0.6, size = 18, face = "bold")) + theme(legend.text = element_text(face = "bold")) + guides(shape = guide_legend(override.aes = list(size = 5)))

print(plot)

dev.off()

# Original vs. imputed.

# Create formula for the plot.
formula <- original.r ~ imputed.r

tiff("original_imputed_CorrScatter.tiff", units="in", width=10, height=7, res=700)
plot <- ggplot(dfPlot, aes(x = original.r, imputed.r)) +
  geom_point(aes(shape = order_name, color = order_name), size = 10, show.legend = T, alpha = 0.7) +
  theme_minimal() +
  scale_color_manual(values = pal) +
  scale_shape_manual(values=c(17, 15, 16, 18)) +
  xlab("Original dataset (r)") +
  ylab("Imputed dataset (r)") +
  ylim(lowerThreshold, upperThreshold) + 
  geom_smooth(method = lm, se = F, show.legend = F, col = "skyblue") +
  theme(axis.title.x = element_text(vjust = 0.6, size = 20, face = "bold"), axis.title.y = element_text(vjust = 0.6, size = 20, face = "bold"), axis.text.y = element_text(vjust = 0.6, size = 18, face = "bold"), axis.text.x = element_text(vjust = 0.6, size = 18, face = "bold")) + theme(legend.text = element_text(face = "bold")) + guides(shape = guide_legend(override.aes = list(size = 5)))

print(plot)

dev.off()
