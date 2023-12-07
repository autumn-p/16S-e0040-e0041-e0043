# Import libraries
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
outPath <- "analysis/out-cleaned"   # Output directory
dataframePath <- "data/ps_all.txt.gz"  # Raw data
appendCol_path <- "data/metadatae0041.tsv"  # Metadata
KCHpalette_path <- "config/KCHcolors-Silva-partial.txt"  # Color palette
referenceASV <- "config/referenceASVs-e0026.txt"  # Reference ASVs

# Read in dataframe
datae0041raw <- read.table(dataframePath, header=TRUE, stringsAsFactors = FALSE)
# Remove unnecessary columns for analysis
datae0041 <- datae0041raw %>%
  dplyr::select(plate, well, OTU, Kingdom, Phylum, Class, Order, Family, Genus, count, relAbundance)
# Adding "subject" column based on the "plate" column
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
  mutate(fullSilvataxonomy=paste(Kingdom, Phylum, Class, Order, Family, sep="."))
KCHpalettee0041 <- KCHpalette %>%
  mutate(taxaNoDashes=gsub("-",".",taxa)) %>%
  filter(taxa %in% sort(unique(datae0041meta$fullSilvataxonomy)) | 
           taxaNoDashes %in% (sort(unique(datae0041meta$fullSilvataxonomy)))) %>%
  mutate(taxashort=ifelse(taxashort=="", gsub(".*\\.","",taxa), taxashort))
# Make a named list
KCHpalettee0041vector <- KCHpalettee0041$hex
names(KCHpalettee0041vector) <- KCHpalettee0041$taxashort


# Create an empty data frame to store the results
result_data <- tibble()

# #Annotate every donor ASV with its ability to colonize the recipient in the mixture community
# dataColonizationSuccess <- foreach(i=unique(datae0041meta %>% filter(donor != "blank" & recipient == "blank") %>% pull(well)), .combine="rbind") %do% {
#   print(i)
#   # Subset the data for the current donor community
#   subset_data <- datae0041meta %>% filter(well == i)
#   i_donor = unique(subset_data$donor)
#   i_replicate = unique(subset_data$replicate)
# 
#   # Finding corresponding wells based on the donor well
#   corresponding_preabx_well <- datae0041meta %>%
#     filter(donor == i_donor & recipient == "XEA-pre" & replicate == i_replicate)
#   corresponding_postabxV1_well <- datae0041meta %>%
#     filter(donor == i_donor & recipient == "XEA-post-V1" & replicate == i_replicate)
#   corresponding_preabx_only_well <- datae0041meta %>%
#     filter(donor == "blank" & recipient == "XEA-pre" & replicate == i_replicate)
#   corresponding_postabxV1_only_well <- datae0041meta %>%
#     filter(donor == "blank" & recipient == "XEA-post-V1" & replicate == i_replicate)
# 
#   # Identify ASVnums for each well
#   asvs_in_donor <- unique(subset_data$OTU)
#   asvs_in_preabx_mix <- unique(corresponding_preabx_well$OTU)
#   asvs_in_postabxV1_mix <- unique(corresponding_postabxV1_well$OTU)
# 
#   # two columns - does it colonize pre successfully and second does it colonize post successfully
#   # third annotation - uncertain (in both d and r)
#   # fourth column - no colonization
#   # ah snap, what about d & r but not in mix???
# 
#   # Colonization preAbx column should have true, false, NA value
#   # Check colonization status
#   colonization_status <- tibble(
#     well = rep(i, length(asvs_in_donor)),
#     OTU = asvs_in_donor,
#     Colonization_Preabx = asvs_in_donor %in% asvs_in_preabx_mix & !(asvs_in_donor %in% corresponding_preabx_only_well$OTU),
#     Colonization_PostabxV1 = asvs_in_donor %in% asvs_in_postabxV1_mix & !(asvs_in_donor %in% corresponding_postabxV1_only_well$OTU),
#     Uncertain_Colonization = (asvs_in_donor %in% asvs_in_preabx_mix & asvs_in_donor %in% corresponding_preabx_only_well$OTU) |
#       (asvs_in_donor %in% asvs_in_postabxV1_mix & asvs_in_donor %in% corresponding_postabxV1_only_well$OTU),
#     Not_Colonizing = !(asvs_in_donor %in% asvs_in_preabx_mix) & !(asvs_in_donor %in% asvs_in_postabxV1_mix)
#   )
# 
#   print(colonization_status)
# 
#   # Append the result to the overall data frame
#   result_data <- bind_rows(result_data, colonization_status)
# }

# Different column version - few less didn't colonize and a few more uncertain colonizers
# Annotate every donor ASV with its ability to colonize the recipient in the mixture community
# dataColonizationSuccess <- foreach(i=unique(datae0041meta %>% filter(donor != "blank" & recipient == "blank") %>% pull(well)), .combine="rbind") %do% {
#   print(i)
#   # Subset the data for the current donor community
#   subset_data <- datae0041meta %>% filter(well == i)
#   i_donor = unique(subset_data$donor)
#   i_replicate = unique(subset_data$replicate)
#   
#   # Finding corresponding wells based on the donor well
#   corresponding_preabx_well <- datae0041meta %>%
#     filter(donor == i_donor & recipient == "XEA-pre" & replicate == i_replicate)
#   corresponding_postabxV1_well <- datae0041meta %>%
#     filter(donor == i_donor & recipient == "XEA-post-V1" & replicate == i_replicate)
#   corresponding_preabx_only_well <- datae0041meta %>%
#     filter(donor == "blank" & recipient == "XEA-pre" & replicate == i_replicate)
#   corresponding_postabxV1_only_well <- datae0041meta %>%
#     filter(donor == "blank" & recipient == "XEA-post-V1" & replicate == i_replicate)
#   
#   # Identify ASVnums for each well
#   asvs_in_donor <- unique(subset_data$OTU)
#   asvs_in_preabx_mix <- unique(corresponding_preabx_well$OTU)
#   asvs_in_postabxV1_mix <- unique(corresponding_postabxV1_well$OTU)
#   asvs_in_preabx_only <- unique(corresponding_preabx_only_well$OTU)
#   asvs_in_postabxV1_only <- unique(corresponding_postabxV1_only_well$OTU)
#   
#   # two columns - does it colonize pre successfully and second does it colonize post successfully
#   # third annotation - uncertain (in both d and r)
#   # fourth column - no colonization
#   # ah snap, what about d & r but not in mix???
#   
#   # Colonization preAbx column should have true, false, NA value
#   # Check colonization status
#   colonization_status <- tibble(
#     well = rep(i, length(asvs_in_donor)),
#     OTU = asvs_in_donor,
#     Colonization_Preabx = asvs_in_donor %in% asvs_in_preabx_mix & !(asvs_in_donor %in% corresponding_preabx_only_well$OTU),
#     Colonization_PostabxV1 = asvs_in_donor %in% asvs_in_postabxV1_mix & !(asvs_in_donor %in% corresponding_postabxV1_only_well$OTU),
#     Uncertain_Colonization = (asvs_in_donor %in% asvs_in_preabx_only) | (asvs_in_donor %in% asvs_in_postabxV1_only),
#     Not_Colonizing = !(asvs_in_donor %in% asvs_in_preabx_mix) | !(asvs_in_donor %in% asvs_in_postabxV1_mix)
#   )
#   
#   print(colonization_status)
#   
#   # Append the result to the overall data frame
#   result_data <- bind_rows(result_data, colonization_status)
# }
# 
# # Merge the result_data with the original datae0041meta dataframe based on well and OTU columns
# datae0041meta <- left_join(datae0041meta, result_data, by = c("well", "OTU"))
# 
# # Create the donor_colonizer_type column
# datae0041meta <- datae0041meta %>%
#   mutate(donor_colonizer_type = case_when(
#     Colonization_Preabx & Colonization_PostabxV1  ~ "Universal_Colonizer",
#     (Colonization_Preabx & !Colonization_PostabxV1) | (!Colonization_Preabx & Colonization_PostabxV1) ~ "Conditional_Colonizer",
#     Uncertain_Colonization ~ "Uncertain_Colonizer",
#     Not_Colonizing ~ "Did_Not_Colonize",
#     TRUE ~ NA_character_
#   ))

# Annotate every donor ASV with its ability to colonize the recipient in the mixture community
dataColonizationSuccess <- foreach(i = unique(datae0041meta %>% filter(donor != "blank" & recipient == "blank") %>% pull(well)), .combine = "rbind") %do% {
  print(i)
  # Subset the data for the current donor community
  subset_data <- datae0041meta %>% filter(well == i)
  i_donor <- unique(subset_data$donor)
  i_replicate <- unique(subset_data$replicate)
  
  # Finding corresponding wells based on the donor well
  corresponding_preabx_well <- datae0041meta %>%
    filter(donor == i_donor & recipient == "XEA-pre" & replicate == i_replicate)
  corresponding_postabxV1_well <- datae0041meta %>%
    filter(donor == i_donor & recipient == "XEA-post-V1" & replicate == i_replicate)
  corresponding_preabx_only_well <- datae0041meta %>%
    filter(donor == "blank" & recipient == "XEA-pre" & replicate == i_replicate)
  corresponding_postabxV1_only_well <- datae0041meta %>%
    filter(donor == "blank" & recipient == "XEA-post-V1" & replicate == i_replicate)
  
  # Incorporate logic directly into subset_data using mutate
  subset_data <- subset_data %>% 
    mutate(
      Colonization_Preabx = OTU %in% unique(corresponding_preabx_well$OTU) & !(OTU %in% unique(corresponding_preabx_only_well$OTU)),
      Colonization_PostabxV1 = OTU %in% unique(corresponding_postabxV1_well$OTU) & !(OTU %in% unique(corresponding_postabxV1_only_well$OTU)),
      Uncertain_Colonization = (OTU %in% unique(corresponding_preabx_only_well$OTU)) | (OTU %in% unique(corresponding_postabxV1_only_well$OTU)),
      Not_Colonizing = !(OTU %in% unique(corresponding_preabx_well$OTU)) | !(OTU %in% unique(corresponding_postabxV1_well$OTU))
    )
  
  print(subset_data)
  
  # Append the result to the overall data frame
  result_data <- bind_rows(result_data, subset_data)
}

# Create the donor_colonizer_type column
result_data <- result_data %>%
  mutate(donor_colonizer_type = case_when(
    Colonization_Preabx & Colonization_PostabxV1  ~ "Universal_Colonizer",
    (Colonization_Preabx & !Colonization_PostabxV1) | (!Colonization_Preabx & Colonization_PostabxV1) ~ "Conditional_Colonizer",
    Uncertain_Colonization ~ "Uncertain_Colonizer",
    Not_Colonizing ~ "Did_Not_Colonize",
    TRUE ~ NA_character_
  ))


# PLAN FOR NEXT ANALYSES
## How many colonizers are in each category, across all of the donors (foreach)...build from there
### Across families? Plot the data in as many ways as possible to do sanity checks
#### Brainstorm plots to check where things may have gone wrong
#### Ex. Number of colonizers that colonizer pre vs post and check with intuition.

# Count the number of colonizers in each category across all donors
colonizer_counts <- result_data %>%
  group_by(donor_colonizer_type) %>%
  summarize(count = n())
# Print results
print(colonizer_counts)

# Bar plot of types of colonizers
colonizer_type_bars <- colonizer_counts%>% ggplot() +
  geom_bar(aes(x = donor_colonizer_type, y = count, fill = donor_colonizer_type), stat = "identity") +
  labs(title = "Colonizer Counts Across Donors",
       x = "Colonizer Type",
       y = "Count") +
  theme_minimal()
colonizer_type_bars
# Save plot
save_plot(paste0(outPath, "/colonizer_type_bars.png"), colonizer_type_bars, base_width = 15, base_height = 10)
# Too many NAs?

# OoOoO same thing in a pie chart (showing the proportion of colonizers in each category)
pie_chart <- datae0041meta %>% ggplot() +
  geom_bar(aes(x = 1, fill = donor_colonizer_type), width = 1) +
  coord_polar(theta = "y") +
  labs(title = "Proportion of Colonizers in Each Category") +
  theme_void()
pie_chart
# Save plot
save_plot(paste0(outPath, "/pie_chart.png"), pie_chart, base_width = 10, base_height = 10)

# Bar plot across families
Colonization_across_families_plot <- result_data %>% ggplot() +
  geom_bar(aes(x = donor_colonizer_type, fill = Family)) +
  #facet_wrap(~donor_colonizer_type) +
  labs(title = "Colonization Status Across Families",
       x = "Donor Colonizer Type",
       y = "Count") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
Colonization_across_families_plot
# Save plot
save_plot(paste0(outPath, "/Colonization_across_families_plot.png"), Colonization_across_families_plot, base_width = 20, base_height = 10)

# Line plot to show how similar overall the replicate 1s are to the replicate 2s in each colonizer type category
line_plot_time <- datae0041meta %>% ggplot() +
  geom_line(aes(x = replicate, y = count, color = donor_colonizer_type)) +
  labs(title = "Colonization Status Over Time",
       x = "Replicate",
       y = "Count") +
  theme_minimal()
line_plot_time
# Save plot
save_plot(paste0(outPath, "/line_plot_time.png"), line_plot_time, base_width = 15, base_height = 10)


# Number of colonizers that colonized pre vs post
pre_post_colonization_plot <- datae0041meta %>% ggplot() +
  geom_bar(aes(x = "Colonization_Preabx", fill = Colonization_Preabx), position = "dodge", show.legend = TRUE) +
  geom_bar(aes(x = "Colonization_PostabxV1", fill = Colonization_PostabxV1), position = "dodge", show.legend = TRUE) +
  labs(title = "Number of Colonizers that Colonized Pre vs Post",
       x = "Colonization Status",
       y = "Count",
       fill = "Colonization Status") +
  theme_minimal()
pre_post_colonization_plot
# Save plot
save_plot(paste0(outPath, "/pre_post_colonization_plot.png"), pre_post_colonization_plot, base_width = 15, base_height = 10)

# Distribution of donor_colonizer_type per donor and corresponding wells
donor_colonizer_distribution_plot <- datae0041meta %>% ggplot() +
  geom_bar(aes(x = donor_colonizer_type, fill = donor_colonizer_type), position = "dodge", show.legend = FALSE) +
  facet_wrap(~ donor + well, scales = "free") +
  labs(title = "Distribution of Donor Colonizer Type per Donor and Corresponding Wells",
       x = "Donor Colonizer Type",
       y = "Count") +
  theme_minimal()
donor_colonizer_distribution_plot
# Save plot
save_plot(paste0(outPath, "/donor_colonizer_distribution_plot.png"), donor_colonizer_distribution_plot, base_width = 30, base_height = 30)

# Replicate comparison plot
replicate_comparison_plot <- datae0041meta %>%
  ggplot() +
  geom_bar(aes(x = Colonization_Preabx, fill = replicate), position = "dodge", show.legend = TRUE) +
  labs(title = "Replicate Comparison",
       x = "Colonization Status (Pre-antibiotics)",
       y = "Count",
       fill = "Replicate") +
  theme_minimal()
replicate_comparison_plot
# Save plot
save_plot(paste0(outPath, "/replicate_comparison_plot.png"), replicate_comparison_plot, base_width = 15, base_height = 10)
# work on this one

#work on this one too
# Well vs. Corresponding Replicate facet-wrapped plots
well_vs_replicate_plot <- datae0041meta %>%
  ggplot() +
  geom_bar(aes(x = Colonization_Preabx, fill = replicate), position = "dodge", show.legend = TRUE) +
  facet_wrap(~well) +
  labs(title = "Well vs. Corresponding Replicate",
       x = "Colonization Status (Pre-antibiotics)",
       y = "Count",
       fill = "Replicate") +
  theme_minimal()
well_vs_replicate_plot
# Save plot
save_plot(paste0(outPath, "/replicate_comparison_plot.png"), replicate_comparison_plot, base_width = 15, base_height = 10)



