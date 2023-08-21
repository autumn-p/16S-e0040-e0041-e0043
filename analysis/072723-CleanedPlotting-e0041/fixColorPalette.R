#Import libraries
library(gridExtra)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)
library(VennDiagram) 
library(patchwork)
library(gridExtra)
library(foreach)

#make the file paths relative*****

# Set up file paths
outPath = "analysis/out-cleaned" #out
dataframePath = "data/ps_all.txt.gz" #raw data
appendCol_path = "data/metadatae0041.tsv" #metadata
KCHpalette_path = "config/KCHcolors-Silva-partial.txt" #color palette
referenceASV = "config/referenceASVs-e0026.txt" #for appending ASV codes

# Read in dataframe
datae0041raw <- read.table(dataframePath, header=TRUE, stringsAsFactors = FALSE)
# Remove columns that are no longer necessary for analysis.
datae0041 <- datae0041raw %>%
  dplyr::select(plate, well, OTU, Kingdom, Phylum, Class, Order, Family, Genus, count, relAbundance)
# Adding "subject" column based on the "plate" columnn
datae0041 <- datae0041 %>% 
  mutate(subject = ifelse(plate == "e0041-A-5", "XEA", ifelse(plate == "e0041-B-5", "XBA", NA)))


# Import metadata table
appendCol <- read.table(appendCol_path, header = TRUE)
# Join the metadata table to the original data frame
datae0041meta <- left_join(datae0041, appendCol, by=c("well", "subject"))
# Filter the data frame to include only rows with relAbundance greater than 0.1%
datae0041meta <- datae0041meta %>% filter(relAbundance > 0.001)
# Raw alpha diversity and alpha diversity w/ limit of detection by well
alpha_diversity_e0041 <- datae0041meta %>%
  group_by(subject, well) %>%
  summarize(alpha_diversity_e0041 = sum(count > 0))
# Join alpha_diversity to the original table
datae0041meta <- datae0041meta %>%
  left_join(alpha_diversity_e0041, by = c('subject', 'well'))


# Import the color palette
KCHpalette <- read.table(KCHpalette_path, header = TRUE) %>%
  mutate(taxashort=gsub("*.*\\.","",taxa))
# Create a stripped-down color palette containing only the families present in the dataset
KCHpalettee0041 <- KCHpalette %>%
  filter(taxashort %in% sort(unique(datae0041$Family)))
# Make a named list
KCHpalettee0041vector <- KCHpalettee0041$hex
names(KCHpalettee0041vector) <- KCHpalettee0041$taxashort

datae0041meta %>%
  filter(well=="B1", plate=="e0041-A-5") %>%
  ggplot() +
  geom_bar(aes(x=well, y=relAbundance, fill=Family), stat="identity") +
  scale_fill_manual(values=KCHpalettee0041vector)

# Debug gray families in community plot.
# Identify the families that are not in the palette vector.
problemFamilies <- unique(datae0041 %>% filter(!(Family %in% KCHpalette$taxaShort)) %>%
                            pull(Family))

