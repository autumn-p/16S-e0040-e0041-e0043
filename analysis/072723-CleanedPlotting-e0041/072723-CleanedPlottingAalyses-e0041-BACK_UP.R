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
library(purrr)

#make the file paths relative*****

# Set up file paths
outPath = "C:/Users/aparr/16S-e0041-e0043/analysis/out-cleaned" #out
dataframePath = "C:/Users/aparr/16S-e0041-e0043/data/ps_all.txt.gz" #raw data
appendCol_path = "C:/Users/aparr/16S-e0041-e0043/data/metadatae0041.tsv" #metadata
KCHpalette_path = "C:/Users/aparr/16S-e0041-e0043/config/KCHcolors-Silva-partial.txt" #color palette
referenceASV = "C:/Users/aparr/16S-e0041-e0043/config/referenceASVs-e0026.txt" #for appending ASV codes

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
KCHpalette <- read.table(KCHpalette_path, header = TRUE)
# Create a stripped-down color palette containing only the families present in the dataset
datae0041meta <- datae0041meta %>%
  mutate(fullSilvataxonomy=paste(Kingdom,Phylum,Class,Order,Family, sep="."))
KCHpalettee0041 <- KCHpalette %>%
  mutate(taxaNoDashes=gsub("-",".",taxa)) %>%
  filter(taxa %in% sort(unique(datae0041meta$fullSilvataxonomy)) | 
           taxaNoDashes %in% (sort(unique(datae0041meta$fullSilvataxonomy)))) %>%
  mutate(taxashort=ifelse(taxashort=="", gsub(".*\\.","",taxa), taxashort))
# Make a named list
KCHpalettee0041vector <- KCHpalettee0041$hex
names(KCHpalettee0041vector) <- KCHpalettee0041$taxashort

# Save metadata table to out-cleaned folder
write.table(datae0041meta, paste0(outPath, "/datae0041meta.txt"), row.names = FALSE, quote = FALSE)


# Creating dataframe with each well of each plate having one row
alpha_diversity_table <- datae0041meta %>%
  group_by(subject, well) %>%
  distinct(alpha_diversity_e0041, .keep_all = TRUE)

# Plotting replicates by community group
# Scatter Plot assessing reliability based on replicates alpha diversity
alpha_diversity_scatterPlot_all <- alpha_diversity_table %>% 
  group_by(donor, recipient) %>%
  ggplot() +
  geom_point(aes(x = interaction(donor, recipient), y = alpha_diversity_e0041, fill = interaction(donor, recipient)), stat = "identity", color = "black") +
  xlab("Community Type") +
  ylab("Alpha Diversity") +
  ggtitle("Alpha Diversity by Replicate") +
  theme_bw() +
  theme(legend.position="none", axis.text.x = element_text(angle = 60, hjust = 1))
alpha_diversity_scatterPlot_all

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

# Alpha Diversity boxplot on Scatter plot with replicates grouped by community mixture type
alpha_diversity_scatterBoxPlot_grouped <- alpha_diversity_table %>% 
  group_by(community_mixture) %>%
  ggplot() +
  geom_jitter(aes(x = community_mixture, y = alpha_diversity_e0041, fill = interaction(community_mixture)), color = "black") +
  geom_boxplot(aes(x = community_mixture, y = alpha_diversity_e0041, fill = interaction(community_mixture)), color = "black") +
  xlab("Community Mixture") +
  ylab("Alpha Diversity") +
  ggtitle("Alpha Diversity by Community") +
  theme_bw() +
  theme(legend.position="none", axis.text.x = element_text(angle = 60, hjust = 1))
alpha_diversity_scatterBoxPlot_grouped
# Save plot
save_plot(paste0(outPath, "/replicateMixGrouped_scatterBoxPlot.png"), alpha_diversity_scatterBoxPlot_grouped, base_width = 10, base_height = 10)


# Generate a stacked bar plot for all of the mixtures involving one recipient community
# x-axis is the donor, y-axis is the relative abundance, facets are replicates of the recipient community
community_abundance_bar_recipient_wrapped <- datae0041meta_mixID %>% 
  filter(community_mixture == "recipient_only") %>% 
  ggplot() +
  geom_bar(aes(x=community_mixture, y=relAbundance, fill=factor(Family)), color = "black", stat= "identity") +
  scale_fill_manual(values=KCHpalettee0041vector) +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Relative Abundance Distribution of Families based on Community Types") +
  facet_wrap(~ well + subject, ncol = 4) 
community_abundance_bar_recipient_wrapped
# Save plot
save_plot(paste0(outPath, "/community_abundance_bar_recipient_facet_wrapped.png"), community_abundance_bar_recipient_wrapped, base_width = 25, base_height = 25)

# Debug gray families in community plot.
# Find the families that are not in the palette vector.
unique(datae0041meta_mixID %>% filter(!(Family %in% names(KCHpalettee0041vector))) %>%
  pull(Family))

# Read in reference ASV table
referenceASV_table <- read.table(referenceASV, header = TRUE, stringsAsFactors = FALSE)
# Convert OTU to upper case
referenceASV_table$OTU <- toupper(referenceASV_table$OTU)
# Join the referenceASV table column ASVnum to the metadata table based on the OTU column
datae0041meta_mixID <- datae0041meta_mixID %>%
  left_join(referenceASV_table %>% select(OTU, ASVnum), by = 'OTU')

# Filter data for only XEA
datae0041meta_mixID <- datae0041meta_mixID %>% filter(subject == "XEA" & !is.na(ASVnum))
  

# Initialize an empty list to store results
#result_list <- list()


# Filter data for the specified donor, recipient, and mixture
# use foreach to merge output
# don't use square brackets
# Loop through each donor+recipient community and add colonization status
dataASVorigin <- foreach(i=unique(datae0041meta_mixID %>% filter(community_mixture == "donor+recipient" | community_mixture == "super+recipient") %>% pull(well)), .combine="rbind") %do% {
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
  subset_data <- subset_data %>% filter(community_mixture == "donor+recipient" | community_mixture =="super+recipient") %>%
    mutate(colonization_status = case_when(
      asvs_in_donor_recipient %in% asvs_in_donor & !(asvs_in_donor_recipient %in% asvs_in_recipient) ~ "Colonizer_Donor",
      asvs_in_donor_recipient %in% asvs_in_recipient & !(asvs_in_donor_recipient %in% asvs_in_donor) ~ "Native_Recipient",
      asvs_in_donor_recipient %in% asvs_in_donor & asvs_in_donor_recipient %in% asvs_in_recipient ~ "Hybrid_D+R",
      TRUE ~ "Weirdo_Neither"
    ))
  
  # Append the result to the list
  #result_list[[i]] <- subset_data
  
  # Return the modified subset_data
  return(subset_data)
  
}

# Combine the results into a single data frame
#final_result <- do.call(rbind, result_list)

# plot the number of new colonizers in each community mixture
plot_colonizers <- dataASVorigin %>%
  filter(colonization_status == "Colonizer_Donor") %>%
  group_by(recipient) %>%
  summarise(numColonizers = n()) %>%
  ggplot() +
  geom_bar(aes(x = recipient, y = numColonizers), stat = "identity") +
  xlab("Recipient Communities") +
  ylab("Count of Colonizers") +
  ggtitle("Number of Donor Colonizers") +
  theme_minimal()
plot_colonizers
# Save Plot
ggsave(filename = paste0(outPath, "/plot_colonizers_per_recipient_comm.png"), plot = plot_colonizers, width = 10, height = 10, units = "cm")


# Get number of colonizers, grouped by donor and recipient
numColonizers <- dataASVorigin %>%
  filter(colonization_status == "Colonizer_Donor") %>%
  group_by(donor, recipient) %>%
  summarise(numColonizers = n())
# Plot number of colonizer
plotColonizersByDonorRecipient <- ggplot(numColonizers, aes(x = donor, y = numColonizers, fill = recipient)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Donor") +
  ylab("Number of Colonizers") +
  ggtitle("Number of Colonizers in Mixture Communities") +
  theme_minimal() +
  theme(legend.position = "top") +
  guides(fill = guide_legend(title = "Recipient"))
plotColonizersByDonorRecipient
# Save plot
save_plot(paste0(outPath, "/plotColonizersByDonorRecipient.png"), plotColonizersByDonorRecipient, base_width = 15, base_height = 10)


# plot the families of the new colonizers in each community mixture
# Filter data to include only rows with "Colonizer_Donor" colonization status
colonizers_per_family <- dataASVorigin %>%
  filter(colonization_status == "Colonizer_Donor") %>%
  count(Family, name = "Count_of_Colonizers")
# Create the bar plot for families of new colonizers
plot_colonizer_families <- colonizers_per_family %>%
  ggplot() +
  geom_bar(aes(x = Family, y = Count_of_Colonizers, fill = Family), stat = "identity") +
  xlab("Families of Colonizer Donors") +
  ylab("Count of Colonizers") +
  ggtitle("Number of Colonizers per Family") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  
plot_colonizer_families
# Save Plot
save_plot(paste0(outPath, "/plot_colonizer_families.png"), plot_colonizer_families, base_width = 15, base_height = 10)


# Colonization Success Plot - plot #/rel abundance of colonizers in each recipient community 
# to compare which communities were the most vulnerable to colonization
# x-axis: recipient community, y-axis: # of colonizers, total % rel abundance of colonizers, 
# plot each mixture as a point
colonization_success_plot <- dataASVorigin %>%
  group_by(well, recipient) %>%
  summarize(total_relAbundance_colonizers = sum(relAbundance[colonization_status == "Colonizer_Donor"], na.rm = TRUE)) %>%
  ggplot() +
  geom_jitter(aes(x = recipient, y = total_relAbundance_colonizers)) +
  xlab("Recipient") +
  ylab("Relative Abundance of Colonizers") +
  ggtitle("Colonization Success By Relative Abundance") +
  theme_minimal() +
  theme(legend.position = "right")
colonization_success_plot
# Save Plot
save_plot(paste0(outPath, "/colonization_success.png"), colonization_success_plot, base_width = 10, base_height = 10)



# Rank Abundance Plot - plot to visualize the distribution of colonizers in the recipient 
# community to identify how dominant the colonizers are in the community
# calculate the total relative abundance of colonizers in each mixture
# show the abundance of colonizers in a stacked bar plot
# (try changing the alpha on the plot based on whether it’s a colonizer)
# try making a stacked bar plot that only has the colonizers
rank_abundance_plot <- dataASVorigin %>% 
  ggplot() +
  geom_bar(aes(x = well, y = relAbundance, fill = factor(Family),
               alpha = ifelse(colonization_status == "Colonizer_Donor", 1, .5)),
           color = "black", stat = "identity") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  xlab("Well") +
  ylab("Relative Abundance") +
  ggtitle("Relative Abundance Distribution of Families based on Community Types") +
  facet_wrap(~ recipient, ncol = 1, scales = "free_x") +
  guides(alpha = "none") +
  theme(strip.text = element_blank()) +
  theme(strip.text = element_text(size = 12, face = "bold", color = "black"),  
        strip.background = element_rect(fill = "lightgray")) +
  labs(strip_text = "Recipient: {strip_variable}") 
rank_abundance_plot
# Save plot
save_plot(paste0(outPath, "/rank_abundance_plot.png"), rank_abundance_plot, base_width = 23, base_height = 15)





# Full opacity plot
full_opacity_plot <- datae0041meta_mixID %>% 
  filter(community_mixture == "donor+recipient") %>% 
  ggplot() +
  geom_bar(aes(x = community_mixture, y = relAbundance, fill = factor(Family)),
           color = "black", stat = "identity") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Full Opacity") +
  facet_wrap(~ well + subject, ncol = 4) +
  guides(fill = "none")
# Modulated opacity plot
modulated_opacity_plot <- datae0041meta_mixID %>% 
  filter(community_mixture == "donor+recipient") %>% 
  ggplot() +
  geom_bar(aes(x = community_mixture, y = relAbundance, fill = factor(Family),
               alpha = ifelse(colonization_status == "Colonizer_Donor", 0.5, 1)),
           color = "black", stat = "identity") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Modulated Opacity") +
  facet_wrap(~ well + subject, ncol = 4) +
  guides(alpha = "none")  # To remove alpha from the legend
# Arrange plots side by side using gridExtra
arranged_colonizer_abundance_stacked <- grid.arrange(full_opacity_plot, modulated_opacity_plot, ncol = 2)
# Save plot
save_plot(paste0(outPath, "/arranged_colonizer_abundance_stacked.png"), arranged_colonizer_abundance_stacked, base_width = 20, base_height = 10)


# Shared Colonizers - a plot to compare the overlap of colonizers between different recipient 
# communities to show which colonizers are consistently successful across various communities
# bar plot - # of times each ASV appears in the list of successful colonizers, compared to the 
# of times it appears
# could also try this at the family level
data_shared_colonizers <- dataASVorigin %>%
  filter(colonization_status == "Colonizer_Donor") %>%
  group_by(ASVnum, recipient) %>%
  summarize(count = n()) %>%
  ungroup()
# Create the Shared Colonizers plot
shared_colonizers_plot <- data_shared_colonizers %>% 
  ggplot() +
  geom_bar(aes(x = ASVnum, y = count, fill = recipient), stat = "identity", position = "dodge") +
  xlab("ASV") +
  ylab("Count of Appearances") +
  ggtitle("Shared Colonizers - Overlap of Successful Colonizers") +
  theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
shared_colonizers_plot
# Save plot
save_plot(paste0(outPath, "/shared_colonizers_plot.png"), shared_colonizers_plot, base_width = 20, base_height = 10)




# Differential Abundance Analysis - identify ASVs that are significantly enriched or depleted 
# in the recipient community after colonization
# maybe we can start by ID-ing ASVs in the recipient communities that disappear in the 
# community mixtures? do they consistently disappear?

# Identify ASVs that disappear going from preAbx to postAbxV1
disappeared_asvs_V1 <- dataASVorigin %>%
  filter(recipient == "XEA-pre" & !(ASVnum %in% dataASVorigin[dataASVorigin$recipient == "XEA-post-V1", "ASVnum"]))

# Identify ASVs that disappear going from preAbx to postAbxV2
disappeared_asvs_V2 <- dataASVorigin %>%
  filter(recipient == "XEA-pre" & !(ASVnum %in% dataASVorigin[dataASVorigin$recipient == "XEA-post-V2", "ASVnum"]))



# Indicator Species Analysis - identify indicator species that are strongly associated with 
# successful colonization in different recipient communities
# maybe we can also try annotating each species in the donor communities and seeing if they 
# can colonize each of the three recipient communities?

# Create new columns for ability to colonize pre/postv1/postv2
abilityToColonizePreV1V2 <- dataASVorigin %>%
  filter(colonization_status == "Colonizer_Donor") %>%
  mutate(
    can_colonize_preAbx = ifelse(
      ASVnum %in%  & recipient == "XEA-pre",
      "Yes",
      "No"
    ),
    can_colonize_postAbxV1 = ifelse(
      ASVnum %in% asvs_in_donor & recipient == "XEA-post-V1",
      "Yes",
      "No"
    ),
    can_colonize_postAbxV2 = ifelse(
      ASVnum %in% asvs_in_donor & recipient == "XEA-post-V2",
      "Yes",
      "No"
    )
  )





