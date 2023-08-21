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


# Set up file paths
outPath = "analysis/out-poster" #out
dataframePath = "data/ps_all.txt.gz" #raw data
appendCol_path = "data/metadatae0041.tsv" #metadata
KCHpalette_path = "config/KCHcolors-Silva-partial.txt" #color palette
referenceASV = "config/referenceASVs-e0026.txt" #for appending ASV codes
source("config/plotDefaults.R")

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

# Save metadata table to out-cleaned folder
write.table(datae0041meta, paste0(outPath, "/datae0041meta.txt"), row.names = FALSE, quote = FALSE)


# Creating dataframe with each well of each plate having one row
alpha_diversity_table <- datae0041meta %>%
  group_by(subject, well) %>%
  distinct(alpha_diversity_e0041, .keep_all = TRUE)
# Based off naming convention community (ex. "XBB-029") or "blank"
alpha_diversity_table <- alpha_diversity_table %>%
  mutate(community_mixture = case_when(
    # super+recipient portion
    donor == "super-community" & recipient != "blank" ~ "super+recipient",
    # recipient_only portion
    donor == "blank" & recipient != "blank" ~ "recipient_only",
    # donor_only portion
    donor != "blank" & recipient == "blank" ~ "donor_only",
    # donor+recipient portion excluding super-community
    donor != "blank" & donor != "super-community" & recipient != "blank" ~ "donor+recipient",
    # everything else
    TRUE ~ NA_character_
  ))

# Add the community mixture ID for each well to the metadata table
datae0041meta_mixID <- datae0041meta %>%
  left_join(alpha_diversity_table %>% select(subject, well, community_mixture), by = c('subject', 'well'))

# Read in reference ASV table
referenceASV_table <- read.table(referenceASV, header = TRUE, stringsAsFactors = FALSE)
# Convert OTU to upper case
referenceASV_table$OTU <- toupper(referenceASV_table$OTU)
# Join the referenceASV table column ASVnum to the metadata table based on the OTU column
datae0041meta_mixID <- datae0041meta_mixID %>%
  left_join(referenceASV_table %>% select(OTU, ASVnum), by = 'OTU')

# Filter data for only XEA
datae0041meta_mixID <- datae0041meta_mixID %>% filter(subject == "XEA", !is.na(ASVnum))


# Filter data for the specified donor, recipient, and mixture
# use foreach to merge output
# don't use square brackets
# Loop through each donor+recipient community and add colonization status
dataASVorigin <- foreach(i=unique(datae0041meta_mixID %>% filter(community_mixture=="donor+recipient") %>% pull(well)), .combine="rbind") %do% {
  #i = "B1"
  # Subset the data for the current donor+recipient community
  subset_data <- datae0041meta_mixID %>% filter(well == i)
  i_donor = unique(subset_data$donor)
  i_recipient = unique(subset_data$recipient)
  i_replicate = unique(subset_data$replicate)
  # Identify the corresponding donor_only and recipient_only wells
  #corresponding_donor_well <- subset_data$well[subset_data$community_mixture == "donor_only"]
  corresponding_donor_well <- datae0041meta_mixID %>%
    filter(donor == i_donor & community_mixture == "donor_only" & replicate == i_replicate)
  corresponding_recipient_well <- datae0041meta_mixID %>%
    filter(recipient == i_recipient & community_mixture == "recipient_only" & replicate == i_replicate)
  
  # Identify ASVnums for each well
  #asvs_in_donor_recipient <- subset_data$ASVnum[subset_data$community_mixture == "donor+recipient"]
  asvs_in_donor_recipient <- unique(subset_data$ASVnum)
  #asvs_in_donor <- subset_data$ASVnum[subset_data$well == corresponding_donor_well]
  asvs_in_donor <- unique(corresponding_donor_well$ASVnum)
  #asvs_in_recipient <- subset_data$ASVnum[subset_data$well == corresponding_recipient_well]
  asvs_in_recipient <- unique(corresponding_recipient_well$ASVnum)
  
  # Colonization Status Update
  subset_data <- subset_data %>% filter(community_mixture == "donor+recipient") %>%
    mutate(colonization_status = case_when(
      asvs_in_donor_recipient %in% asvs_in_donor & !(asvs_in_donor_recipient %in% asvs_in_recipient) ~ "Colonizer_Donor",
      asvs_in_donor_recipient %in% asvs_in_recipient & !(asvs_in_donor_recipient %in% asvs_in_donor) ~ "Native_Recipient",
      asvs_in_donor_recipient %in% asvs_in_donor & asvs_in_donor_recipient %in% asvs_in_recipient ~ "Hybrid_D+R",
      TRUE ~ "Weirdo_Neither"
    ))
  # Append the result to the list
  #result_list[[i]] <- subset_data
  
  return(subset_data)
  
}

# POSTER 
# recipient stacked plot
# Remove legends
# Fix gray
# Titles are small
# Theme change
# Generate a stacked bar plot for all of the mixtures involving one recipient community
# x-axis is the donor, y-axis is the relative abundance, facets are replicates of the recipient community
community_abundance_bar_recipient_wrapped <- datae0041meta_mixID %>% 
  filter(community_mixture == "recipient_only") %>% 
  ggplot() +
  geom_bar(aes(x=recipient, y=relAbundance, fill=factor(Family)), color = "black", stat= "identity") +
  scale_fill_manual(values=KCHpalettee0041vector) +
  xlab("Recipient") +
  ylab("Relative Abundance") +
  ggtitle("Relative Abundance Distribution of Recipient Communities") +
  facet_wrap(~ well ~recipient, scales = "free_x", ncol = 4) +
  DEFAULTS.THEME_PRES +
  scale_x_discrete(breaks = NULL) 
community_abundance_bar_recipient_wrapped
# Save plot
save_plot(paste0(outPath, "/community_abundance_bar_recipient_facet_wrapped.png"), community_abundance_bar_recipient_wrapped, base_width = 10, base_height = 15)

# POSTER
# Donors Stacked Bar Plot
community_abundance_bar_donor_wrapped <- datae0041meta_mixID %>% 
  filter(community_mixture == "donor_only") %>% 
  ggplot() +
  geom_bar(aes(x=donor, y=relAbundance, fill=factor(Family)), color = "black", stat= "identity") +
  scale_fill_manual(values=KCHpalettee0041vector) +
  xlab("Donor") +
  ylab("Relative Abundance") +
  ggtitle("Relative Abundance Distribution of Donor Communities") +
  facet_wrap(~ well ~donor, scales = "free_x", ncol = 4) +
  DEFAULTS.THEME_PRES +
  scale_x_discrete(breaks = NULL) 
community_abundance_bar_donor_wrapped
# Save plot
save_plot(paste0(outPath, "/community_abundance_bar_donor_wrapped.png"), community_abundance_bar_donor_wrapped, base_width = 10, base_height = 15)


# POSTER - change to recipient community type
# Remove legends
# Fix gray
# Titles are small
# Theme change
# Define the desired order of community mixtures

# Define the desired order of community mixtures
desired_order <- c("donor_only", "recipient_only", "donor+recipient")
# X-Axis Labels
desired_labels <- c("Donor Only", "Recipient Only", "Donor + Recipient")

datae0041meta_mixID <- datae0041meta_mixID %>%
  mutate(recipient = ifelse(recipient == "blank", "no_recipient", recipient))

# Alpha Diversity for D, R, Mix
alpha_diversity_summary <- datae0041meta_mixID %>%
  mutate(community_mixture = ifelse(
    community_mixture == "super+recipient", "donor+recipient", community_mixture)) %>%
  mutate(community_mixture = factor(community_mixture, levels = desired_order)) %>%
  group_by(well) %>%
  summarize(alpha_diversity = first(alpha_diversity_e0041),
            community_mixture = first(community_mixture),
            recipient_community = ifelse(community_mixture == "donor+recipient" | community_mixture == "recipient_only" , recipient, "no_recipient")) %>%
  filter(!is.na(community_mixture))
  

# Create the scatter plot
alpha_diversity_scatterPlot <- alpha_diversity_summary %>%
  ggplot() +
  geom_boxplot(aes(x = community_mixture, y = alpha_diversity)) +
  geom_point(aes(x = community_mixture, y = alpha_diversity, 
                 color = recipient_community),
             position = position_jitter(width = 0.2)) +
  xlab("Community Mixture") +
  ylab("Number of Species") +
  scale_color_brewer(palette="Dark2") +
  ggtitle("Alpha Diversity by Community") +
  DEFAULTS.THEME_PRES +
  scale_x_discrete(labels = desired_labels) +
  scale_fill_discrete(labels = c("no_recipient" = "No Recipient"))
alpha_diversity_scatterPlot

# Save plot
save_plot(paste0(outPath, "/alpha_diversity_summary.png"), alpha_diversity_scatterPlot, base_width = 10, base_height = 10)



# POSTER
# Remove legends
# Paste all together
# Fix gray
# Titles are small
# Theme change
# Pre-Abx D+R=Mixture
unique_preAbx_donor_asv_count <- datae0041meta_mixID %>%
  filter(well=="A1") %>%
  distinct(ASVnum) %>%
  nrow()
# Donor stacked bar plot
donor_abundance_sum_preAbx <- datae0041meta_mixID %>% 
  filter(well=="A1") %>% 
  ggplot() +
  geom_bar(aes(x = community_mixture, y = relAbundance, fill = factor(Family)),
           color = "black", stat = "identity") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Relative Abundance Distribution of Donor Community") +
  geom_text(aes(x = 1, y = 1.05, label = paste(unique_preAbx_donor_asv_count, "Species")),
            hjust = 0.5, size = 4) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRES +
  theme_nothing()
donor_abundance_sum_preAbx
# Save plot
save_plot(paste0(outPath, "/donor_abundance_sum_preAbx.png"), donor_abundance_sum_preAbx, base_width = 10, base_height = 10)

# Recipient Abundance Plot - plot to visualize the distribution of recipients 
# Calculate the number of unique asv values
unique_preAbx_recipient_asv_count <- datae0041meta_mixID %>%
  filter(well=="B12") %>%
  distinct(ASVnum) %>%
  nrow()
# Stacked bar plot 
recipient_abundance_sum_preAbx <- datae0041meta_mixID %>% 
  filter(well == "B12") %>% 
  ggplot() +
  geom_bar(aes(x = community_mixture, y = relAbundance, fill = factor(Family)),
           color = "black", stat = "identity") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Relative Abundance Distribution of Recipient Community") +
  geom_text(aes(x = 1, y = 1.05, label = paste(unique_preAbx_recipient_asv_count, "Species")),
            hjust = 0.5, size = 4) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRES +
  theme_nothing()
recipient_abundance_sum_preAbx
# Save plot
save_plot(paste0(outPath, "/recipient_abundance_sum_preAbx.png"), recipient_abundance_sum_preAbx, base_width = 10, base_height = 10)


# Full opacity plot
# Calculate the number of unique asv values
unique_preAbx_mixture_asv_count <- datae0041meta_mixID %>%
  filter(well=="B1") %>%
  distinct(ASVnum) %>%
  nrow()
# Stacked Bar Plot
full_opacity_mix_preAbx_sum <- datae0041meta_mixID %>% 
  filter(well=="B1") %>% 
  ggplot() +
  geom_bar(aes(x = community_mixture, y = relAbundance, fill = factor(Family)),
           color = "black", stat = "identity") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Mixture") +
  geom_text(aes(x = 1, y = 1.05, label = paste(unique_preAbx_mixture_asv_count, "Species")),
            hjust = 0.5, size = 4) +
  guides(fill = "none") +
  DEFAULTS.THEME_PRES +
  theme_nothing()
full_opacity_mix_preAbx_sum

# Modulated opacity plot
modulated_opacity_preAbx_sum <- dataASVorigin %>% 
  filter(well=="B1") %>% 
  ggplot() +
  geom_bar(aes(x = community_mixture, y = relAbundance, fill = factor(Family),
               alpha = ifelse(colonization_status == "Colonizer_Donor", 1, 0)),
           color = "black", stat = "identity") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  scale_alpha_identity() +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Colonizers") +
  geom_text(aes(x = 1, y = 1.05, label = paste(unique_preAbx_mixture_asv_count, "Species")),
                      hjust = 0.5, size = 4) +
  guides(fill = "none")  +
  guides(alpha = "none") +
  DEFAULTS.THEME_PRES +
  theme_nothing()
modulated_opacity_preAbx_sum

# Arrange plots side by side using gridExtra
arranged_colonizer_abundance_stacked_preAbx <- grid.arrange(donor_abundance_sum_preAbx, recipient_abundance_sum_preAbx, full_opacity_mix_preAbx_sum, modulated_opacity_preAbx_sum, ncol = 4)

# Save plot
save_plot(paste0(outPath, "/arranged_colonizer_abundance_stacked_preAbx.png"), arranged_colonizer_abundance_stacked_preAbx, base_width = 20, base_height = 10)





#POSTER
#PostAbx V1 D+R=Mix
# POSTER
# Remove legends
# Paste all together
# Fix gray
# Titles are small
# Theme change

#Donor unique ASV
unique_postAbxV1_donor_asv_count <- datae0041meta_mixID %>%
  filter(well=="A2") %>%
  distinct(ASVnum) %>%
  nrow()
# Donor stacked bar plot
donor_abundance_sum_postAbxV1 <- datae0041meta_mixID %>% 
  filter(well=="A2") %>% 
  ggplot() +
  geom_bar(aes(x = community_mixture, y = relAbundance, fill = factor(Family)),
           color = "black", stat = "identity") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Relative Abundance Distribution of Donor Community") +
  geom_text(aes(x = 1, y = 1.05, label = paste(unique_postAbxV1_donor_asv_count, "Species")),
            hjust = 0.5, size = 4) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRES +
  theme_nothing()
donor_abundance_sum_postAbxV1

# Recipient Abundance Plot - plot to visualize the distribution of recipients 
# Calculate the number of unique asv values
unique_postAbxV1_recipient_asv_count <- datae0041meta_mixID %>%
  filter(well=="D12") %>%
  distinct(ASVnum) %>%
  nrow()
# Stacked bar plot 
recipient_abundance_sum_postAbxV1 <- datae0041meta_mixID %>% 
  filter(well == "D12") %>% 
  ggplot() +
  geom_bar(aes(x = community_mixture, y = relAbundance, fill = factor(Family)),
           color = "black", stat = "identity") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Relative Abundance Distribution of Recipient Community") +
  geom_text(aes(x = 1, y = 1.05, label = paste(unique_postAbxV1_recipient_asv_count, "Species")),
            hjust = 0.5, size = 4) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRES +
  theme_nothing()
recipient_abundance_sum_postAbxV1

# Full opacity plot
# Calculate the number of unique asv values
unique_postAbxV1_mixture_asv_count <- datae0041meta_mixID %>%
  filter(well=="D2") %>%
  distinct(ASVnum) %>%
  nrow()
# Stacked Bar Plot
full_opacity_mix_postAbxV1_sum <- datae0041meta_mixID %>% 
  filter(well=="D2") %>% 
  ggplot() +
  geom_bar(aes(x = community_mixture, y = relAbundance, fill = factor(Family)),
           color = "black", stat = "identity") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Mixture") +
  geom_text(aes(x = 1, y = 1.05, label = paste(unique_postAbxV1_mixture_asv_count, "Species")),
            hjust = 0.5, size = 4) +
  guides(fill = "none") +
  DEFAULTS.THEME_PRES +
  theme_nothing()
full_opacity_mix_postAbxV1_sum

# Modulated opacity plot
modulated_opacity_postAbxV1_sum <- dataASVorigin %>% 
  filter(well=="D2") %>% 
  ggplot() +
  geom_bar(aes(x = community_mixture, y = relAbundance, fill = factor(Family),
               alpha = ifelse(colonization_status == "Colonizer_Donor", 1, 0)),
           color = "black", stat = "identity") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  scale_alpha_identity() +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Colonizers") +
  geom_text(aes(x = 1, y = 1.05, label = paste(unique_postAbxV1_mixture_asv_count, "Species")),
            hjust = 0.5, size = 4) +
  guides(fill = "none")  +
  guides(alpha = "none") +
  DEFAULTS.THEME_PRES +
  theme_nothing()
modulated_opacity_postAbxV1_sum

# Arrange plots side by side using gridExtra
arranged_colonizer_abundance_stacked_postAbxV1 <- grid.arrange(donor_abundance_sum_postAbxV1, recipient_abundance_sum_postAbxV1, full_opacity_mix_postAbxV1_sum, modulated_opacity_postAbxV1_sum, ncol = 4)

# Save plot
save_plot(paste0(outPath, "/arranged_colonizer_abundance_stacked_postAbxV1.png"), arranged_colonizer_abundance_stacked_postAbxV1, base_width = 20, base_height = 10)


#POSTER
#PostAbx V2 D+R=Mix
# POSTER
# Remove legends
# Paste all together
# Fix gray
# Titles are small
# Theme change

#Donor unique ASV
unique_postAbxV2_donor_asv_count <- datae0041meta_mixID %>%
  filter(well=="A5") %>%
  distinct(ASVnum) %>%
  nrow()
# Donor stacked bar plot
donor_abundance_sum_postAbxV2 <- datae0041meta_mixID %>% 
  filter(well=="A5") %>% 
  ggplot() +
  geom_bar(aes(x = community_mixture, y = relAbundance, fill = factor(Family)),
           color = "black", stat = "identity") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Relative Abundance Distribution of Donor Community") +
  geom_text(aes(x = 1, y = 1.05, label = paste(unique_postAbxV2_donor_asv_count, "Species")),
            hjust = 0.5, size = 4) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRES +
  theme_nothing()
donor_abundance_sum_postAbxV2

# Recipient Abundance Plot - plot to visualize the distribution of recipients 
# Calculate the number of unique asv values
unique_postAbxV2_recipient_asv_count <- datae0041meta_mixID %>%
  filter(well=="F12") %>%
  distinct(ASVnum) %>%
  nrow()
# Stacked bar plot 
recipient_abundance_sum_postAbxV2 <- datae0041meta_mixID %>% 
  filter(well == "F12") %>% 
  ggplot() +
  geom_bar(aes(x = community_mixture, y = relAbundance, fill = factor(Family)),
           color = "black", stat = "identity") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Relative Abundance Distribution of Recipient Community") +
  geom_text(aes(x = 1, y = 1.05, label = paste(unique_postAbxV2_recipient_asv_count, "Species")),
            hjust = 0.5, size = 4) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRES +
  theme_nothing()
recipient_abundance_sum_postAbxV2

# Full opacity plot
# Calculate the number of unique asv values
unique_postAbxV2_mixture_asv_count <- datae0041meta_mixID %>%
  filter(well=="F5") %>%
  distinct(ASVnum) %>%
  nrow()
# Stacked Bar Plot
full_opacity_mix_postAbxV2_sum <- datae0041meta_mixID %>% 
  filter(well=="F5") %>% 
  ggplot() +
  geom_bar(aes(x = community_mixture, y = relAbundance, fill = factor(Family)),
           color = "black", stat = "identity") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Mixture") +
  geom_text(aes(x = 1, y = 1.05, label = paste(unique_postAbxV2_mixture_asv_count, "Species")),
            hjust = 0.5, size = 4) +
  guides(fill = "none") +
  DEFAULTS.THEME_PRES +
  theme_nothing()
full_opacity_mix_postAbxV2_sum

# Modulated opacity plot
modulated_opacity_postAbxV2_sum <- dataASVorigin %>% 
  filter(well=="F5") %>% 
  ggplot() +
  geom_bar(aes(x = community_mixture, y = relAbundance, fill = factor(Family),
               alpha = ifelse(colonization_status == "Colonizer_Donor", 1, 0)),
           color = "black", stat = "identity") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  scale_alpha_identity() +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Colonizers") +
  geom_text(aes(x = 1, y = 1.05, label = paste(unique_postAbxV2_mixture_asv_count, "Species")),
            hjust = 0.5, size = 4) +
  guides(fill = "none")  +
  guides(alpha = "none") +
  DEFAULTS.THEME_PRES +
  theme_nothing()
modulated_opacity_postAbxV2_sum

# Arrange plots side by side using gridExtra
arranged_colonizer_abundance_postAbxV2 <- grid.arrange(donor_abundance_sum_postAbxV2, recipient_abundance_sum_postAbxV2, full_opacity_mix_postAbxV2_sum, modulated_opacity_postAbxV2_sum, ncol = 4)

# Save plot
save_plot(paste0(outPath, "/arranged_colonizer_abundance_postAbxV2.png"), arranged_colonizer_abundance_postAbxV2, base_width = 20, base_height = 10)



# POSTER
# Colonizer Rel Abundance Percent

# Define the desired order of community mixtures
desired_order <- c("XEA-pre", "XEA-post-V1", "XEA-post-V2")
# X-Axis Labels
desired_labels <- c("Pre-Antibiotics", "Post-Antibiotics V1", "Post-Antibiotics V2")

# Calculate the percentage relative abundance of colonizers within each well
# Create the scatterplot
colonizer_percent_scatterplot <- dataASVorigin %>%
  mutate(recipient = factor(recipient, levels = desired_order)) %>%
  group_by(well, recipient) %>%
  summarize(colonizer_percent = sum(colonization_status == "Colonizer_Donor") / n() * 100) %>%
  filter(colonizer_percent!=0, colonizer_percent!=100) %>%
  ggplot() +
  geom_point(aes(x = recipient, y = colonizer_percent), position = position_jitter(width = 0.2)) +
  xlab("Recipient Community") +
  ylab("Percentage Colonizers (%)") +
  ggtitle("Percentage Colonizers") +
  DEFAULTS.THEME_PRES +
  scale_x_discrete(labels = desired_labels) + ylim(0,100)
colonizer_percent_scatterplot
# Save plot
save_plot(paste0(outPath, "/colonizer_percent_scatterplot.png"), colonizer_percent_scatterplot, base_width = 5, base_height = 5)

# POSTER
# Number of Colonizers per bacterial family
# Calculate the number of colonizers within each bacterial family
# Create the scatterplot - PRE
pre_colonizer_family_scatterplot <- dataASVorigin %>%
  #  & !Family=="Lachnospiraceae"
  filter(colonization_status == "Colonizer_Donor" & recipient == "XEA-pre") %>%
  group_by(Family) %>%
  summarize(colonizer_count = n()) %>%
  ggplot() +
  geom_point(aes(x = reorder(Family, colonizer_count), y = colonizer_count), stat = "identity", fill = "dodgerblue") +
  xlab("Bacterial Family") +
  ylab("Number of Colonizers") +
  ggtitle("Pre-Antibiotics") +
  DEFAULTS.THEME_PRES +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
pre_colonizer_family_scatterplot
# Save plot
save_plot(paste0(outPath, "/pre_colonizer_family_scatterplot.png"), pre_colonizer_family_scatterplot, base_width = 5, base_height = 5)

# Create the scatterplot - POST V1
postV1_colonizer_family_scatterplot <- dataASVorigin %>%
  #  & !Family=="Lachnospiraceae"
  filter(colonization_status == "Colonizer_Donor" & recipient == "XEA-post-V1") %>%
  group_by(Family) %>%
  summarize(colonizer_count = n()) %>%
  ggplot() +
  geom_point(aes(x = reorder(Family, colonizer_count), y = colonizer_count), stat = "identity", fill = "dodgerblue") +
  xlab("Bacterial Family") +
  ylab("Number of Colonizers") +
  ggtitle("Post-Antibiotics V1") +
  DEFAULTS.THEME_PRES +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
postV1_colonizer_family_scatterplot
# Save plot
save_plot(paste0(outPath, "/postV1_colonizer_family_scatterplot.png"), postV1_colonizer_family_scatterplot, base_width = 5, base_height = 5)

# Create the scatterplot - POST V2
postV2_colonizer_family_scatterplot <- dataASVorigin %>%
  #  & !Family=="Lachnospiraceae"
  filter(colonization_status == "Colonizer_Donor" & recipient == "XEA-post-V2") %>%
  group_by(Family) %>%
  summarize(colonizer_count = n()) %>%
  ggplot() +
  geom_point(aes(x = reorder(Family, colonizer_count), y = colonizer_count), stat = "identity", fill = "dodgerblue") +
  xlab("Bacterial Family") +
  ylab("Number of Colonizers") +
  ggtitle("Post-Antibiotics V2") +
  DEFAULTS.THEME_PRES +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
postV2_colonizer_family_scatterplot
# Save plot
save_plot(paste0(outPath, "/postV2_colonizer_family_scatterplot.png"), postV2_colonizer_family_scatterplot, base_width = 5, base_height = 5)
