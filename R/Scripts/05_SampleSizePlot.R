# Script for creating sample size plots (complete-case vs. original).

# Load libaries and functions. ----
library(data.table)
library(ggplot2)
library(tidyverse)

# Sample size plot. ----

# Create vector of sample sizes.
sampleSizesCC <- c(Anguilliformes = 215, Characiformes = 193, Clupeiformes = 119, Cypriniformes = 117, Gadiformes = 164, Perciformes = 228, Pleuronectiformes = 206, Scorpaeniformes = 447, Siluriformes = 225, Anura = 417, Caudata = 156, Apodiformes = 176, Charadriiformes = 183, Galliformes = 110, Passeriformes = 387, Piciformes = 141, Psittaciformes = 105, Carnivora = 162, Chiroptera = 241, Rodentia = 177, Squamata = 152)

# Create a new dataframe containing sample size info.
dfSS <- data.frame(cc = sampleSizesCC)
# Add order_name column.
dfSS$order_name <- rownames(dfSS)
# Remove row names.
rownames(dfSS) <- NULL
# Create a grouping variable so we can add colour to the plot.
dfSS$order_name
dfSS$group <- c(rep("Actinopterygii", 9), rep("Amphibia", 2), rep("Aves", 6), rep("Mammalia", 3), "Reptilia")
# Convert to factor type (this is necessary for ggplot).
dfSS$order_name <- factor(dfSS$order_name, levels = dfSS$order_name[order(dfSS$group)])

# Plots.
# Complete-case sample sizes AFTER matching with sequence data AND allowing for some missingness (aka nearly-complete). ---
# Set classic theme.
theme_set(theme_classic())
# Write plot to current working directory.
tiff("Complete-case_Barplot.tiff", units="in", width=10, height=6, res=300)

plot <- ggplot(data = dfSS, aes(order_name, cc, fill = group)) + 
  geom_bar(stat="identity", width=0.5) + 
  scale_fill_manual("Legend", values = c("lightskyblue2", "lightgreen", "lightgoldenrod", "coral", "rosybrown2")) + ylab("Sample Size") + xlab("Order") +
  labs(title="", 
       subtitle="", 
       caption="") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6, size = 11)) + coord_flip()
print(plot)

dev.off()
