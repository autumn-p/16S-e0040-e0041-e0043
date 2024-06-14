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

# Save result_data table to out-cleaned folder
write.table(result_data, paste0(outPath, "/result_data.txt"), row.names = FALSE, quote = FALSE)

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
pie_chart <- result_data %>% ggplot() +
  geom_bar(aes(x = 1, fill = donor_colonizer_type), width = 1) +
  coord_polar(theta = "y") +
  labs(title = "Proportion of Colonizers in Each Category") +
  theme_void()
pie_chart
# Save plot
save_plot(paste0(outPath, "/pie_chart.png"), pie_chart, base_width = 10, base_height = 10)

# Bar plot across families
#<<<<<<< Updated upstream
Colonization_across_families_plot <- result_data %>% ggplot() +
  geom_bar(aes(x = donor_colonizer_type, fill = Family)) +
  #facet_wrap(~donor_colonizer_type) +
  labs(title = "Colonization Status Across Families",
       x = "Donor Colonizer Type",
       y = "Count") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#=======
Colonization_across_families_plot <- datae0041meta %>% filter(donor_colonizer_type == "Conditional_Colonizer") %>% ggplot() +
  geom_bar(aes(x = Family, fill = Family)) +
  labs(title = "Colonization Status Across Families",
       x = "Donor Colonizer Type",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1))
#>>>>>>> Stashed changes
Colonization_across_families_plot
# Save plot
save_plot(paste0(outPath, "/Colonization_across_families_plot.png"), Colonization_across_families_plot, base_width = 20, base_height = 10)

# Bar plot across faceted families
Colonization_across_families_faceted_plot <- result_data %>% ggplot() +
  geom_bar(aes(x = donor_colonizer_type, fill = Family)) +
  facet_wrap(~Family) +
  labs(title = "Colonization Status Across Families",
       x = "Donor Colonizer Type",
       y = "Count") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
Colonization_across_families_faceted_plot
# Save plot
save_plot(paste0(outPath, "/Colonization_across_families_faceted_plot.png"), Colonization_across_families_faceted_plot, base_width = 20, base_height = 10)

colonizer_families <- result_data %>% ggplot() +
  geom_bar(aes(x = Family, fill = Family)) +
  facet_wrap(~donor_colonizer_type, ncol=1) +
  labs(title = "Colonization Status Across Families",
       x = "Donor Colonizer Type",
       y = "Count") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
colonizer_families
# Save plot
save_plot(paste0(outPath, "/colonizer_families.png"), colonizer_families, base_width = 20, base_height = 10)

# Conditional colonizers plotted with distribution of pre vs post
result_data_counts <- result_data %>%
  group_by(donor_colonizer_type, Family) %>%
  summarise(Preabx_count = sum(Colonization_Preabx),
            Postabx_count = sum(Colonization_PostabxV1)) %>%
  ungroup()
Colonization_conditional_pre_post_plot <- result_data_counts %>% filter(donor_colonizer_type == "Conditional_Colonizer") %>% ggplot() +
  geom_bar(aes(x = donor_colonizer_type, y = Preabx_count, fill = "Preabx"), stat = "identity", position = "dodge") +
  geom_bar(aes(x = donor_colonizer_type, y = Postabx_count, fill = "Postabx"), stat = "identity", position = "dodge") +
  facet_wrap(~Family) +
  labs(title = "Colonization Status Across Families",
       x = "Donor Colonizer Type",
       y = "Count") +
  scale_fill_manual(values = c("Preabx" = "blue", "Postabx" = "red")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
Colonization_conditional_pre_post_plot
# Save plot
save_plot(paste0(outPath, "/Colonization_conditional_pre_post_plot.png"), Colonization_conditional_pre_post_plot, base_width = 20, base_height = 10)


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

# Reshape the data to long format for pre and post conditions
result_data_long <- datae0041meta %>%
  pivot_longer(cols = c(Colonization_Preabx, Colonization_PostabxV1, Colonization_PostabxV2),
               names_to = "colonization_stage",
               values_to = "colonized") %>%
  filter(colonized) %>%
  group_by(Family, colonization_stage) %>%
  summarize(number_of_colonizers = n(), .groups = 'drop')

# Update the factor levels for colonization_stage
result_data_long$colonization_stage <- factor(result_data_long$colonization_stage, 
                                              levels = c("Colonization_Preabx", "Colonization_PostabxV1", "Colonization_PostabxV2"),
                                              labels = c("Pre-Abx", "Post-Abx V1", "Post-Abx V2"))

# Create the plot
pre_post_colonization_plot <- result_data_long %>%
  ggplot(aes(x = Family, y = number_of_colonizers, fill = Family)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ colonization_stage) +
  labs(title = "Number of Colonizers Pre and Post Antibiotic Treatment",
       x = "Family",
       y = "Number of Colonizers") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Display the plot
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

# Progress since last meeting with Kat
# Normalizes the counts within each family based on the total number of rows for each 
# combination of Family and donor_colonizer_type and makes all total to 1
result_data_normalized <- result_data %>%
  group_by(Family, donor_colonizer_type) %>%
  mutate(Normalized_Count = 1 / n()) %>%
  ungroup()

# Bar plot across faceted families with normalized counts
Colonization_combo_across_families_faceted_normalized_plot <- result_data_normalized %>% ggplot() +
  geom_bar(aes(x = donor_colonizer_type, fill = Family, y = Normalized_Count), stat = "identity", position = "dodge") +
  facet_wrap(~Family) +
  labs(title = "Normalized Colonization Status Across Families",
       x = "Donor Colonizer Type",
       y = "Normalized Count") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
Colonization_combo_across_families_faceted_normalized_plot
# Save plot
save_plot(paste0(outPath, "/Colonization_combo_across_families_faceted_normalized_plot.png"), Colonization_combo_across_families_faceted_normalized_plot, base_width = 15, base_height = 10)

#*******
#bar plot with # of colonizers that colonized the pre vs post abx comm; x-axis families and y-axis 
#number of colonizers ; facets for pre, post v1, and v2
# Reshape the data to long format for faceting by pre and post conditions
result_data_long <- result_data_normalized %>%
  pivot_longer(cols = c(Colonization_Preabx, Colonization_PostabxV1),
               names_to = "colonization_stage",
               values_to = "colonized") %>%
  filter(colonized) %>%
  group_by(Family, colonization_stage) %>%
  summarize(number_of_colonizers = n(), .groups = 'drop')

# Update the factor levels for colonization_stage
result_data_long$colonization_stage <- factor(result_data_long$colonization_stage, 
                                              levels = c("Colonization_Preabx", "Colonization_PostabxV1"),
                                              labels = c("Pre-Abx", "Post-Abx"))

# Create the plot
Colonization_normalized_families_comm_plot <- result_data_long %>%
  ggplot(aes(x = Family, y = number_of_colonizers, fill = Family)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ colonization_stage) +
  labs(title = "Number of Colonizers Pre and Post Antibiotic Treatment",
       x = "Family",
       y = "Number of Colonizers") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
Colonization_normalized_families_comm_plot
# Save plot
save_plot(paste0(outPath, "/Colonization_normalized_families_comm_plot.png"), Colonization_normalized_families_comm_plot, base_width = 15, base_height = 10)

# Plotting difference in colonizers
# Reshape the data to long format for pre and post conditions
result_data_long <- result_data_normalized %>%
  pivot_longer(cols = c(Colonization_Preabx, Colonization_PostabxV1),
               names_to = "colonization_stage",
               values_to = "colonized") %>%
  filter(colonized) %>%
  group_by(Family, colonization_stage) %>%
  summarize(number_of_colonizers = n(), .groups = 'drop')

# Spread the data to wide format to calculate differences
result_data_wide <- result_data_long %>%
  pivot_wider(names_from = colonization_stage, values_from = number_of_colonizers, values_fill = list(number_of_colonizers = 0)) %>%
  mutate(difference_in_colonizers = `Colonization_PostabxV1` - `Colonization_Preabx`)

# Update the factor levels for consistency
result_data_wide$Family <- factor(result_data_wide$Family)

# Create the plot for the difference in colonizers
Colonization_difference_families_comm_plot <- result_data_wide %>%
  ggplot(aes(x = Family, y = difference_in_colonizers, fill = Family)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Difference in Number of Colonizers Pre and Post Antibiotic Treatment",
       x = "Family",
       y = "Difference in Number of Colonizers (Post - Pre)") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
Colonization_difference_families_comm_plot

# Save plot
save_plot(paste0(outPath, "/Colonization_difference_families_comm_plot.png"), Colonization_difference_families_comm_plot, base_width = 15, base_height = 10)

#new version with just conditional
# Plot the filtered data
Colonization_conditional_faceted_normalized_plot <- result_data_normalized %>% filter(donor_colonizer_type == "Conditional_Colonizer") %>% ggplot() +
  geom_bar(aes(x = donor_colonizer_type, fill = Family, y = Normalized_Count), stat = "identity", position = "dodge") +
  facet_wrap(~Family) +
  labs(title = "Normalized Colonization Status Across Families",
       x = "Donor Colonizer Type",
       y = "Normalized Count") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
Colonization_conditional_faceted_normalized_plot
# Save plot
save_plot(paste0(outPath, "/Colonization_conditional_faceted_normalized_plot.png"), Colonization_conditional_faceted_normalized_plot, base_width = 15, base_height = 10)

# ORRR
# Nevermind this is wronggg
# Normalize counts within each family
result_data_normalized <- result_data %>%
  group_by(Family) %>%
  mutate(Normalized_Count = sum(!is.na(donor_colonizer_type))) %>%
  ungroup()

# Bar plot across faceted families with normalized counts
Colonization_across_families_faceted_normalized_plot <- result_data_normalized %>% ggplot() +
  geom_bar(aes(x = donor_colonizer_type, fill = Family, y = Normalized_Count), stat = "identity", position = "dodge") +
  facet_wrap(~Family) +
  labs(title = "Normalized Colonization Status Across Families",
       x = "Donor Colonizer Type",
       y = "Normalized Count") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
Colonization_across_families_faceted_normalized_plot
# Save plot
save_plot(paste0(outPath, "/Colonization_across_families_faceted_normalized_plot.png"), Colonization_across_families_faceted_normalized_plot, base_width = 15, base_height = 10)

# Distribution of donor_colonizer_type per donor and corresponding wells
cleaned_donor_colonizer_distribution_plot <- result_data %>%
  ggplot(aes(x = donor_colonizer_type, fill = donor_colonizer_type)) +
  geom_bar(position = "dodge") +
  facet_wrap(~ donor + well, scales = "free") +
  labs(title = "Distribution of Donor Colonizer Type per Donor and Corresponding Wells",
       x = "Donor Colonizer Type",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
cleaned_donor_colonizer_distribution_plot
# Save plot
save_plot(paste0(outPath, "/cleaned_donor_colonizer_distribution_plot.png"), cleaned_donor_colonizer_distribution_plot, base_width = 15, base_height = 20)

library(ggplot2)

# Distribution of donor_colonizer_type per donor and corresponding wells
replicate_distribution_plot <- result_data %>%
  ggplot(aes(x = donor_colonizer_type, fill = donor_colonizer_type)) +
  geom_bar(position = "dodge") +
  facet_grid(replicate ~ donor + well, scales = "free") +
  labs(title = "Distribution of Donor Colonizer Type per Replicate and Corresponding Wells",
       x = "Donor Colonizer Type",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  
replicate_distribution_plot
# Save plot
save_plot(paste0(outPath, "/replicate_distribution_plot.png"), replicate_distribution_plot, base_width = 15, base_height = 10)

# Distribution of donor_colonizer_type per donor and corresponding wells
replicate_adjusted_plot <- result_data %>%
  ggplot(aes(x = donor_colonizer_type, fill = donor_colonizer_type)) +
  geom_bar(position = "dodge") +
  facet_grid(donor ~ replicate + well, scales = "free") +
  labs(title = "Distribution of Donor Colonizer Type per Donor, Replicate, and Corresponding Wells",
       x = "Donor Colonizer Type",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  
replicate_adjusted_plot
# Save plot
save_plot(paste0(outPath, "/replicate_adjusted_plot.png"), replicate_adjusted_plot, base_width = 15, base_height = 10)

# Subset data for Sutterellaceae
sutterellaceae_data <- result_data %>% filter(Family == "Sutterellaceae")

# Bar plot of colonization types within Sutterellaceae
colonization_sutte_plot <- sutterellaceae_data %>% ggplot() +
  geom_bar(aes(x = donor_colonizer_type, fill = donor_colonizer_type), position = "dodge") +
  labs(title = "Colonization Status Within Sutterellaceae",
       x = "Donor Colonizer Type",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
colonization_sutte_plot

# Subset data for Sutterellaceae in pre-antibiotic and post-antibioticV1
sutterellaceae_pre <- result_data %>% filter(Family == "Sutterellaceae" & Colonization_Preabx == TRUE)
sutterellaceae_postV1 <- result_data %>% filter(Family == "Sutterellaceae" & Colonization_PostabxV1 == TRUE)

# Bar plot for colonization types in Sutterellaceae pre and post-antibioticV1
sutte_colonization_pre_post_plot <- bind_rows(mutate(sutterellaceae_pre, AntibioticStatus = "Pre"),
                                        mutate(sutterellaceae_postV1, AntibioticStatus = "PostV1")) %>%
  ggplot() +
  geom_bar(aes(x = donor_colonizer_type, fill = AntibioticStatus), position = "dodge") +
  labs(title = "Colonization Status in Sutterellaceae - Pre vs. Post-AntibioticV1",
       x = "Donor Colonizer Type",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
sutte_colonization_pre_post_plot
# Save plot
save_plot(paste0(outPath, "/sutte_colonization_pre_post_plot.png"), sutte_colonization_pre_post_plot, base_width = 15, base_height = 10)

# Subset data for Rikenellaceae in pre-antibiotic and post-antibioticV1
rikenellaceae_pre <- result_data %>%
  filter(Family == "Rikenellaceae" & Colonization_Preabx == TRUE)

rikenellaceae_postV1 <- result_data %>%
  filter(Family == "Rikenellaceae" & Colonization_PostabxV1 == TRUE)

# Bar plot for colonization types in Rikenellaceae pre and post-antibioticV1
rikenellaceae_pre_post_plot <- bind_rows(mutate(rikenellaceae_pre, AntibioticStatus = "Pre"),
                                         mutate(rikenellaceae_postV1, AntibioticStatus = "PostV1")) %>%
  ggplot() +
  geom_bar(aes(x = donor_colonizer_type, fill = AntibioticStatus), position = "dodge") +
  labs(title = "Colonization Status in Rikenellaceae - Pre vs. Post-AntibioticV1",
       x = "Donor Colonizer Type",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
rikenellaceae_pre_post_plot
# Save plot
save_plot(paste0(outPath, "/rikenellaceae_pre_post_plot.png"), rikenellaceae_pre_post_plot, base_width = 15, base_height = 10)

# Subset data for Erysipelotrichaceae in pre-antibiotic and post-antibioticV1
erysipelotrichaceae_pre <- result_data %>%
  filter(Family == "Erysipelotrichaceae" & Colonization_Preabx == TRUE)

erysipelotrichaceae_postV1 <- result_data %>%
  filter(Family == "Erysipelotrichaceae" & Colonization_PostabxV1 == TRUE)

# Bar plot for colonization types in Erysipelotrichaceae pre and post-antibioticV1
erysipelotrichaceae_pre_post_plot <- bind_rows(mutate(erysipelotrichaceae_pre, AntibioticStatus = "Pre"),
                                               mutate(erysipelotrichaceae_postV1, AntibioticStatus = "PostV1")) %>%
  ggplot() +
  geom_bar(aes(x = donor_colonizer_type, fill = AntibioticStatus), position = "dodge") +
  labs(title = "Colonization Status in Erysipelotrichaceae - Pre vs. Post-AntibioticV1",
       x = "Donor Colonizer Type",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
erysipelotrichaceae_pre_post_plot
# Save plot
save_plot(paste0(outPath, "/erysipelotrichaceae_pre_post_plot.png"), erysipelotrichaceae_pre_post_plot, base_width = 15, base_height = 10)

# Initialize an empty data frame for postV2_result_data
postV2_result_data <- data.frame()

# Annotate every donor ASV with its ability to colonize the recipient in the mixture community
postV2dataColonizationSuccess <- foreach(i = unique(datae0041meta %>% filter(donor != "blank" & recipient == "blank") %>% pull(well)), .combine = "rbind") %do% {
  print(i) 
  # Subset the data for the current donor community
  postV2_subset_data <- datae0041meta %>% filter(well == i)
  i_donor <- unique(postV2_subset_data$donor)
  i_replicate <- unique(postV2_subset_data$replicate)
  
  # Finding corresponding wells based on the donor well
  corresponding_preabx_well <- datae0041meta %>%
    filter(donor == i_donor & recipient == "XEA-pre" & replicate == i_replicate)
  corresponding_postabxV1_well <- datae0041meta %>%
    filter(donor == i_donor & recipient == "XEA-post-V1" & replicate == i_replicate)
  corresponding_postabxV2_well <- datae0041meta %>%
    filter(donor == i_donor & recipient == "XEA-post-V2" & replicate == i_replicate)
  corresponding_preabx_only_well <- datae0041meta %>%
    filter(donor == "blank" & recipient == "XEA-pre" & replicate == i_replicate)
  corresponding_postabxV1_only_well <- datae0041meta %>%
    filter(donor == "blank" & recipient == "XEA-post-V1" & replicate == i_replicate)
  corresponding_postabxV2_only_well <- datae0041meta %>%
    filter(donor == "blank" & recipient == "XEA-post-V2" & replicate == i_replicate)
  
  # Incorporate logic directly into subset_data using mutate
  postV2_subset_data <- postV2_subset_data %>% 
    mutate(
      Colonization_Preabx = OTU %in% unique(corresponding_preabx_well$OTU) & !(OTU %in% unique(corresponding_preabx_only_well$OTU)),
      Colonization_PostabxV1 = OTU %in% unique(corresponding_postabxV1_well$OTU) & !(OTU %in% unique(corresponding_postabxV1_only_well$OTU)),
      Colonization_PostabxV2 = OTU %in% unique(corresponding_postabxV2_well$OTU) & !(OTU %in% unique(corresponding_postabxV2_only_well$OTU)),
      Uncertain_Colonization = (OTU %in% unique(corresponding_preabx_only_well$OTU)) | (OTU %in% unique(corresponding_postabxV1_only_well$OTU)) | (OTU %in% unique(corresponding_postabxV2_only_well$OTU)),
      Not_Colonizing = !(OTU %in% unique(corresponding_preabx_well$OTU)) | !(OTU %in% unique(corresponding_postabxV1_well$OTU)) | !(OTU %in% unique(corresponding_postabxV2_well$OTU))
    )
  
  print(postV2_subset_data)
  
  # Append the result to the overall data frame
  postV2_result_data <- bind_rows(postV2_result_data, postV2_subset_data)
}

# Create the donor_colonizer_type column
postV2_result_data <- postV2_result_data %>%
  mutate(donor_colonizer_type = case_when(
    Colonization_Preabx & Colonization_PostabxV1  & Colonization_PostabxV2 ~ "Universal_Colonizer",
    (Colonization_Preabx & !Colonization_PostabxV1 & !Colonization_PostabxV2) | (!Colonization_Preabx & Colonization_PostabxV1 & !Colonization_PostabxV2) | (!Colonization_Preabx & !Colonization_PostabxV1 & Colonization_PostabxV2)~ "Conditional_Colonizer",
    Uncertain_Colonization ~ "Uncertain_Colonizer",
    Not_Colonizing ~ "Did_Not_Colonize",
    TRUE ~ NA_character_
  ))

# Subset data for each family in postV1 and postV2
family_postV1 <- postV2_result_data %>% filter(Colonization_PostabxV1 == TRUE)
family_postV2 <- postV2_result_data %>% filter(Colonization_PostabxV2 == TRUE)

# Combine data for plotting
family_comparison <- bind_rows(mutate(family_postV1, AntibioticStatus = "PostV1"),
                               mutate(family_postV2, AntibioticStatus = "PostV2"))

# Create the bar plot
family_colonization_comparison_plot <- family_comparison %>%
  ggplot(aes(x = Family, fill = AntibioticStatus)) +
  geom_bar(position = "dodge") +
  labs(title = "Colonization Status Comparison between PostV1 and PostV2 by Family",
       x = "Family",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
family_colonization_comparison_plot
# Save plot
save_plot(paste0(outPath, "/family_colonization_comparison_plot.png"), family_colonization_comparison_plot, base_width = 15, base_height = 10)

# Kat's step 1
statistics_data <- result_data %>%
  filter(donor_colonizer_type=="Conditional_Colonizer") %>%
  select(OTU, Family, donor_colonizer_type)

# Kat's step 2 broken down

count_observed <- function(result_data) {
  # Step 1: Filter data for conditional colonizers
  conditional_colonizers <- result_data %>%
    filter(donor_colonizer_type == "Conditional_Colonizer" & plate=="e0041-A-5")
  
  # Step 2: Identify families that are different between pre-abx and post-abx V1
  pre_abx_families <- conditional_colonizers %>%
    filter(Colonization_Preabx == "TRUE") %>%
    select(Family) %>%
    distinct()
  
  post_abx_v1_families <- conditional_colonizers %>%
    filter(Colonization_PostabxV1 == "TRUE") %>%
    select(Family) %>%
    distinct()
  
  # Either conditional colonizer family to pre-abx or post-abx-V1
  unique_to_pre_abx <- setdiff(unique(pre_abx_families$Family), unique(post_abx_v1_families$Family))
  unique_to_post_abx_v1 <- setdiff(unique(post_abx_v1_families$Family), unique(pre_abx_families$Family))
  
  # Combine the unique families from both pre-abx and post-abx V1
  different_families <- union(unique_to_pre_abx, unique_to_post_abx_v1)
  
  # Step 3: Count the number of conditional colonizers in different families
  observed_count <- conditional_colonizers %>%
    filter(Family %in% different_families) %>%
    select(OTU, Family, donor_colonizer_type)
  
  return(observed_count)
}

# Probably do this if I actually used the function
#observed_count <- count_observed(result_data)
#print(observed_count)

# Kat's Step 3
# Function to shuffle data
shuffled_data <- result_data %>%
  select(OTU, Family, donor_colonizer_type) %>%
  slice(sample(n()))
  
# Kat's Step 4 count
# Function to count differential colonizers in families of interest
observed_count_permuted <- shuffled_data %>%
  filter(Family %in% different_families) %>%
  count(Family, name = "Count")


# Kat's step 4 permute
# Function to generate permuted values
generate_permuted_values <- function(shuffled_data, different_families, num_permutations = 100) {
  permuted_values <- replicate(num_permutations, {
    count_differential_colonizers(shuffled_data, different_families)
  })
  
  return(permuted_values)
}




# Number of permutations
num_permutations <- 100

# Initialize a list to store permuted values
permuted_values <- vector("list", length = num_permutations)

# Perform permutations
for (i in 1:num_permutations) {
  # Shuffle data
  shuffled_data <- result_data %>%
    select(OTU, Family, donor_colonizer_type) %>%
    slice(sample(n()))
  
  # Count differential colonizers in families of interest
  observed_count_permuted <- shuffled_data %>%
    filter(Family %in% different_families) %>%
    count(Family, name = "Count")
  
  # Store the result in the list
  permuted_values[[i]] <- observed_count_permuted
}
