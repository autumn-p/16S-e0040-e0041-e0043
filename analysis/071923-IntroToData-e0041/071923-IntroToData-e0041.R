#Import libraries
library(gridExtra)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)
library(VennDiagram) 
library(patchwork)
library(gridExtra)

# Set up file paths
outPath = "/analysis/out" # out
dataframePath = "/data/ps_all.txt.gz" #raw data
appendCol_path = "data/metadatae0041.tsv" #metadata
KCHpalette_path = "/config/KCHcolors-Silva-partial.txt" #color palette
referenceASV = "/config/referenceASVs-e0026.txt"

# Read in dataframe
datae0041raw <- read.table(dataframePath, header=TRUE, stringsAsFactors = FALSE)
# Remove columns that are no longer necessary for analysis.
datae0041 <- datae0041raw %>%
  dplyr::select(plate, well, OTU, Kingdom, Phylum, Class, Order, Family, Genus, count, relAbundance)


# Add the "plate" column to the data frame
#datae0041$plate <- plate_values

# Adding "subject" column based on the "plate" columnn
datae0041 <- datae0041 %>% 
  mutate(subject = ifelse(plate == "e0041-A-5", "XEA", ifelse(plate == "e0041-B-5", "XBA", NA)))

# Import metadata table
appendCol <- read.table(appendCol_path, header = TRUE)

# Join the metadata table to the original data frame
datae0041meta <- left_join(datae0041, appendCol, by=c("well", "subject"))

# Filter the data frame to include only rows with relAbundance greater than 0.1%
datae0041meta <- datae0041meta %>% filter(relAbundance > 0.001)

# Raw alpha diversity and alpha diversity w/ limit of detection**
# Get alpha diversity by well
alpha_diversity_e0041 <- datae0041meta %>%
  group_by(subject, well) %>%
  summarize(alpha_diversity_e0041 = sum(count > 0))

# Join alpha_diversity to the original table**
datae0041meta <- datae0041meta %>%
  left_join(alpha_diversity_e0041, by = c('subject', 'well'))

# Import the color palette
KCHpalette <- read.table(KCHpalette_path, header = TRUE)

# Create a stripped-down color palette containing only the families
# present in the dataset.
datae0041meta <- datae0041meta %>%
  mutate(fullSilvataxonomy=paste(Kingdom,Phylum,Class,Order,Family, sep="."))
KCHpalettee0041 <- KCHpalette %>%
  filter(taxa %in% sort(unique(datae0041meta$fullSilvataxonomy))) %>%
  mutate(taxashort=ifelse(taxashort=="", gsub(".*\\.","",taxa), taxashort))

# Make a named list
KCHpalettee0041vector <- KCHpalettee0041$hex
names(KCHpalettee0041vector) <- KCHpalettee0041$taxashort

# Save table with metadata to out folder
write.table(datae0041meta, paste0(outPath, "/datae0041meta.txt"), row.names = FALSE, quote = FALSE)




# Bar plot showing the relative abundance of each community  
#   filter(recipient %in% c("XEA-pre", "XEA-post-V1", "XEA-post-V2", "XBA-pre-abx", "XBA-post-abx", "XBA-drop-out")) %>% 
# Create a bar plot using ggplot2
communityAbundanceBarByCommunity <- datae0041meta %>% 
  filter(well=="B1" & subject=="XBA") %>% 
  ggplot() +
  geom_bar(aes(x=replicate, y=relAbundance, fill=factor(Family)), color = "black", stat= "identity") +
  scale_fill_manual(values=KCHpalettee0041vector) +
  facet_wrap(donor~recipient) +
  xlab("Recipient") +
  ylab("Relative Abundance") +
  ggtitle("Distribution of Family Relative Abundance")
communityAbundanceBarByCommunity

# Save plot
save_plot(paste0(outPath, "/communityAbundanceBarByComm2.png"), communityAbundanceBarByCommunity, base_width = 15, base_height = 15)

# Creating dataframe with each well of each plate having one row
alpha_diversity_table <- datae0041meta %>%
  group_by(subject, well) %>%
  distinct(alpha_diversity_e0041, .keep_all = TRUE)



# Bar plot to compare the distribution of alpha diversity values among a strainOnly vs donor + recipient
# Subset the data by community type (donor and recipient)

# Create a bar plot using ggplot2
alpha_diversity_barPlot <- alpha_diversity_table %>%
  filter(well == "B11" & subject == "XEA") %>%
  ggplot() +
  geom_bar(aes(x = paste("donor: ", donor, " recipient: ", recipient), y = alpha_diversity_e0041), stat = "identity", color = "black", fill = "blue") +
  xlab("Community Type") +
  ylab("Alpha Diversity") +
  ggtitle("Distribution of Alpha Diversity by Community Type") +
  theme_bw()
alpha_diversity_barPlot

# Save plot
save_plot(paste0(outPath, "/alphaDiversity_XEApre_barPlot.png"), alpha_diversity_barPlot, base_width = 5, base_height = 5)






# Create a bar plot using ggplot2 with average alpha diversity for each community type
#donor only vs donor + recipient
#B11 vs B12
#XEA pre + super community vs XEA pre
#D11 vs D12
#XEA post V1 + super community vs XEA post V1
#F11 vs F12
#XEA post V2 + super community vs XEA post V2

XEA_alpha_diversity_barPlot <- alpha_diversity_table %>% 
  filter(well %in% c("B11", "B12", "D11", "D12", "F11", "F12") & subject == "XEA") %>%
  group_by(donor, recipient) %>%
  ggplot() +
  geom_bar(aes(x = interaction(donor, recipient), y = alpha_diversity_e0041, fill = interaction(donor, recipient)), stat = "identity", color = "black") +
  xlab("Community Type") +
  ylab("Alpha Diversity") +
  ggtitle("Distribution of Alpha Diversity by Community Type") +
  theme_bw()
XEA_alpha_diversity_barPlot

# Save plot
save_plot(paste0(outPath, "/alphaDiversity_XEAprepostV1V2_barPlot.png"), XEA_alpha_diversity_barPlot, base_width = 15, base_height = 10)


#B11 vs B12
#XBA pre + super community vs XBA pre
#D11 vs D12
#XBA post + super community vs XBA post
#F11 vs F12
#XBA drop out + super community vs XBA drop out
XBA_alpha_diversity_barPlot <- alpha_diversity_table %>%
  filter(well %in% c("B11", "B12", "D11", "D12", "F11", "F12") & subject == "XBA") %>%
  group_by(donor, recipient) %>%
  ggplot() +
  geom_bar(aes(x = interaction(donor, recipient), y = alpha_diversity_e0041, fill = interaction(donor, recipient)), stat = "identity", color = "black") +
  xlab("Community Type") +
  ylab("Alpha Diversity") +
  ggtitle("Distribution of Alpha Diversity by Community Type") +
  theme_bw()
XBA_alpha_diversity_barPlot

# Save plot
save_plot(paste0(outPath, "/alphaDiversity_XBAprepostdrop_barPlot.png"), XBA_alpha_diversity_barPlot, base_width = 15, base_height = 10)





# Just replicate 1 & 2 under same conditions
# Scatter Plot assessing reliability based on replicates alpha diversity
alpha_diversity_scatterPlot <- alpha_diversity_table %>% 
  filter(well %in% c("B1", "C1") & subject == "XEA") %>%
  group_by(donor, recipient) %>%
  ggplot() +
  geom_point(aes(x = interaction(donor, recipient), y = alpha_diversity_e0041, fill = interaction(donor, recipient)), stat = "identity", color = "black") +
  xlab("Community Type") +
  ylab("Alpha Diversity") +
  ggtitle("Alpha Diversity by Replicate") +
  theme_bw()
alpha_diversity_scatterPlot

# Save plot
save_plot(paste0(outPath, "/single_replicateReliability_scatterplot.png"), alpha_diversity_scatterPlot)



# New variable "community_type" with the combination of "donor" and "recipient"
#datae0041meta <- datae0041meta %>%
#  mutate(community_type = paste("donor:", donor, "recipient:", recipient))



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

# Save plot
save_plot(paste0(outPath, "/all_replicateReliability_scatterplot.png"), alpha_diversity_scatterPlot_all, base_width = 20, base_height = 15)



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

# SUMMARY PLOT
# Alpha Diversity Scatter plot with replicates grouped by community mixture type
alpha_diversity_scatterPlot_grouped <- alpha_diversity_table %>% 
  group_by(community_mixture) %>%
  ggplot() +
  geom_jitter(aes(x = community_mixture, y = alpha_diversity_e0041, fill = interaction(community_mixture)), color = "black") +
  xlab("Community Mixture") +
  ylab("Alpha Diversity") +
  ggtitle("Alpha Diversity by Community") +
  theme_bw() +
  theme(legend.position="none", axis.text.x = element_text(angle = 60, hjust = 1))
alpha_diversity_scatterPlot_grouped

# Save plot
save_plot(paste0(outPath, "/replicateMixGrouped_scatterplot.png"), alpha_diversity_scatterPlot_grouped, base_width = 10, base_height = 10)

# SUMMARY PLOT
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




datae0041meta_mixID <- datae0041meta %>%
  left_join(alpha_diversity_table %>% select(subject, well, community_mixture), by = c('subject', 'well'))


# SUMMARY PLOT
# Create a stacked bar plot to visualize the taxonomic composition of the communities
# Create a bar plot using ggplot2
#one bar
community_abundance_bar_summary <- datae0041meta_mixID %>% 
  filter(donor == 'XBB-029' & recipient == 'XEA-pre' & replicate == '1') %>% 
  #filter(community_mixture %in% c("donor_only", "recipient_only", "super+recipient", "donor+recipient")) %>% 
  #filter(well == "B1" & subject == "XEA") %>% 
  ggplot() +
  geom_bar(aes(x=community_mixture, y=relAbundance, fill=factor(Family)), color = "black", stat= "identity") +
  scale_fill_manual(values=KCHpalettee0041vector) +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Relative Abundance Distribution of Families based on Community Types")
community_abundance_bar_summary

save_plot(paste0(outPath, "/community_abundance_bar_summary.png"), community_abundance_bar_summary, base_width = 10, base_height = 10)

#many bars
community_abundance_bar_summary_grouped <- datae0041meta_mixID %>% 
  #filter(donor == 'XBB-029' & recipient == 'XEA-pre' & replicate == '1') %>% 
  filter(community_mixture %in% c("donor_only", "recipient_only", "super+recipient", "donor+recipient") & subject == 'XEA' & replicate == '1') %>% 
  #filter(well == "B1" & subject == "XEA") %>% 
  ggplot() +
  geom_bar(aes(x=community_mixture, y=relAbundance, fill=factor(Family)), color = "black", stat= "identity") +
  scale_fill_manual(values=KCHpalettee0041vector) +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Relative Abundance Distribution of Families based on Community Types") +
  facet_wrap(~ donor + recipient, ncol = 4) 
community_abundance_bar_summary_grouped

save_plot(paste0(outPath, "/community_abundance_bar_summary_grouped.png"), community_abundance_bar_summary_grouped, base_width = 25, base_height = 25)

#Make a stacked bar plot for one donor community
community_abundance_bar_donorOnly <- datae0041meta_mixID %>% 
  filter(community_mixture == "donorOnly" & well == "A1" & subject == "XEA") %>% 
  ggplot() +
  geom_bar(aes(x=community_mixture, y=relAbundance, fill=factor(Family)), color = "black", stat= "identity") +
  scale_fill_manual(values=KCHpalettee0041vector) +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Relative Abundance Distribution of Families based on Community Types") +
  #facet_wrap(~ well + subject, ncol = 4)  
community_abundance_bar_donorOnly

save_plot(paste0(outPath, "/community_abundance_bar_donorOnly.png"), community_abundance_bar_donorOnly, base_width = 25, base_height = 25)


# Make a stacked bar plot for one recipient community
community_abundance_bar_recipientOnly <- datae0041meta_mixID %>% 
  filter(community_mixture == "recipientOnly") %>% 
  ggplot() +
  geom_bar(aes(x=community_mixture, y=relAbundance, fill=factor(Family)), color = "black", stat= "identity") +
  scale_fill_manual(values=KCHpalettee0041vector) +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Relative Abundance Distribution of Families based on Community Types") +
  facet_wrap(~ well + subject, ncol = 4)  
community_abundance_bar_recipientOnly

save_plot(paste0(outPath, "/community_abundance_bar_recipientOnly.png"), community_abundance_bar_recipientOnly, base_width = 25, base_height = 25)



# Make a stacked bar plot for that donor + recipient community
community_abundance_bar_donor_recipient <- datae0041meta_mixID %>% 
  filter(community_mixture == "donor+recipient") %>% 
  ggplot() +
  geom_bar(aes(x=community_mixture, y=relAbundance, fill=factor(Family)), color = "black", stat= "identity") +
  scale_fill_manual(values=KCHpalettee0041vector) +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Relative Abundance Distribution of Families based on Community Types") +
  facet_wrap(~ well + subject, ncol = 4)
community_abundance_bar_donor_recipient

save_plot(paste0(outPath, "/community_abundance_bar_donor_recipient.png"), community_abundance_bar_donor_recipient, base_width = 25, base_height = 25)





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

save_plot(paste0(outPath, "/community_abundance_bar_recipient_facet_wrapped.png"), community_abundance_bar_recipient_wrapped, base_width = 25, base_height = 25)


# append the ASV names according to the current reference database
# referenceASV database e0032 - use the file titled config/referenceASVs-e0026.txt 
# the column we’re interested in is called ASVnum, and you can join based on the OTU column

# Read in reference ASV table
referenceASV_table <- read.table(referenceASV, header = TRUE, stringsAsFactors = FALSE)

# Convert OTU to upper case
referenceASV_table$OTU <- toupper(referenceASV_table$OTU)

# Join the referenceASV table column ASVnum to the metadata table based on the OTU column
datae0041meta_mixID <- datae0041meta_mixID %>%
  left_join(referenceASV_table %>% select(OTU, ASVnum), by = 'OTU')


# write a function to take a donor community, recipient community, and community mixture, 
# and return the list of ASVs from the donor community that successfully colonize the 
# recipient community

# Doesn't previous limit of detection take care of having to do it again here?

# Filter data for the specified donor, recipient, and mixture
# Get ASVs in the donor-only community (well A1)
asvs_in_donor <- datae0041meta_mixID %>% 
  filter(subject == "XEA" & well == "A1" & community_mixture == "donor_only") %>%
  pull(ASVnum)
  
# Get ASVs in the recipient-only community (well B12)
asvs_in_recipient <- datae0041meta_mixID %>% 
  filter(subject == "XEA" & well == "B12" & community_mixture == "recipient_only") %>%
  pull(ASVnum)

#  Get ASVs in the donor + recipient community mixture (well B1)
asvs_in_donor_recipient <- datae0041meta_mixID %>% 
  filter(subject == "XEA" & well == "B1" & community_mixture == "donor+recipient") %>%
  pull(ASVnum)

# Append a new column to the donor_data dataframe specifying colonization status

beyonce_filtered <- datae0041meta_mixID %>%
  filter(well %in% c("A1", "B1", "B12") & subject == "XEA")

beyonce_B1_filtered$colonization_status <- NA

# use OTU not ASVnum (if there are overlapping ASVs)
# change to case_when()

beyonce_B1_filtered <- beyonce_B1_filtered %>% filter(community_mixture == "donor+recipient") %>%
  mutate(colonization_status = case_when(
    asvs_in_donor_recipient %in% asvs_in_donor & !(asvs_in_donor_recipient %in% asvs_in_recipient) ~ "Colonizer_Donor",
    asvs_in_donor_recipient %in% asvs_in_recipient & !(asvs_in_donor_recipient %in% asvs_in_donor) ~ "Native_Recipient",
    asvs_in_donor_recipient %in% asvs_in_donor & asvs_in_donor_recipient %in% asvs_in_recipient ~ "Hybrid_D+R",
    TRUE ~ "Weirdo_Neither"
  ))
  

# plot the number of new colonizers in each community mixture

# Filter data to include only rows with "Colonizer_Donor" colonization status
count_colonizers <- beyonce_B1_filtered %>%
  filter(colonization_status == "Colonizer_Donor") %>%
  group_by(well) %>% summarize(numColonizers=n())

# Create a data frame with the count of colonizers
df_colonizers <- data.frame(Label = "Colonizer Donors", Count = count_colonizers)

# Create the bar plot
plot_colonizers <- count_colonizers %>%
  ggplot() +
  geom_bar(aes(x = well, y = numColonizers), stat = "identity") +
  xlab("Colonizer Donors") +
  ylab("Count of Colonizers") +
  ggtitle("Number of Donor Colonizers in Well B1 (XBA)") +
  theme_minimal()
plot_colonizers

# Save Plot
#save_plot(paste0(outPath, "/plot_colonizers_B1.png"), plot_colonizers_B1, base_width = 10, base_height = 10)

ggsave(filename = paste0(outPath, "/plot_colonizers_B1.png"), plot = plot_colonizers, width = 10, height = 10, units = "cm")


# plot the families of the new colonizers in each community mixture
# Filter data to include only rows with "Colonizer_Donor" colonization status
colonizers_per_family <- beyonce_B1_filtered %>%
  filter(colonization_status == "Colonizer_Donor") %>%
  count(Family, name = "Count_of_Colonizers")

# Create the bar plot for families of new colonizers
plot_colonizer_families_B1 <- colonizers_per_family %>%
  ggplot() +
  geom_bar(aes(x = Family, y = Count_of_Colonizers, fill = Family), stat = "identity") +
  xlab("Families of Colonizer Donors") +
  ylab("Count of Colonizers") +
  ggtitle("Number of Colonizers per Family in Well B1") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  
plot_colonizer_families_B1

# Save Plot
save_plot(paste0(outPath, "/plot_colonizer_families_B1.png"), plot_colonizer_families_B1, base_width = 15, base_height = 10)


# compare the number of new colonizers for different recipient communities

# I would change the fill to recipient right? or group by before i filter in pip in

# Colonization Success Plot - plot #/rel abundance of colonizers in each recipient community 
# to compare which communities were the most vulnerable to colonization
# x-axis: recipient community, y-axis: # of colonizers, total % rel abundance of colonizers, 
# plot each mixture as a point

colonization_success_plot <- beyonce_B1_filtered %>%
  group_by(subject, well) %>%
  summarize(total_colonizers = sum(colonization_status == "Colonizer_Donor", na.rm = TRUE),
            total_relAbundance_colonizers = sum(relAbundance[colonization_status == "Colonizer_Donor"], na.rm = TRUE)) %>%
  ggplot() +
  geom_point(aes(x = well, y = total_colonizers, size = total_relAbundance_colonizers)) +
  xlab("Well") +
  ylab("# of Colonizers") +
  ggtitle("Colonization Success Plot") +
  theme_minimal() +
  theme(legend.position = "right")

colonization_success_plot


# Save Plot
save_plot(paste0(outPath, "/colonization_success_B1.png"), colonization_success_plot, base_width = 15, base_height = 10)

# Relative abundance colonization success

colonization_success_plot_relAbundance <- beyonce_B1_filtered %>%
  group_by(subject, well) %>%
  summarize(total_relAbundance_colonizers = sum(as.numeric(relAbundance[colonization_status == "Colonizer_Donor"]), na.rm = TRUE)) %>%
  ggplot() +
  geom_point(aes(x = well, y = total_relAbundance_colonizers)) +
  xlab("Well") +
  ylab("Relative Abundance of Colonizers") +
  ggtitle("Colonization Success Plot - Relative Abundance") +
  theme_minimal() +
  theme(legend.position = "right")
colonization_success_plot_relAbundance

# Save plot
save_plot(paste0(outPath, "/colonization_success_B1_relAbundance.png"), colonization_success_plot_relAbundance, base_width = 15, base_height = 10)




# Rank Abundance Plot - plot to visualize the distribution of colonizers in the recipient 
# community to identify how dominant the colonizers are in the community
# calculate the total relative abundance of colonizers in each mixture
# show the abundance of colonizers in a stacked bar plot
# (try changing the alpha on the plot based on whether it’s a colonizer)
# try making a stacked bar plot that only has the colonizers

rank_abundance_plot <- beyonce_B1_filtered %>% 
  filter(community_mixture == "donor+recipient") %>% 
  ggplot() +
  geom_bar(aes(x = community_mixture, y = relAbundance, fill = factor(Family),
               alpha = ifelse(colonization_status == "Colonizer_Donor", 1, .5)),
           color = "black", stat = "identity") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Relative Abundance Distribution of Families based on Community Types") +
  facet_wrap(~ well + subject, ncol = 4) +
  guides(alpha = "none")  # To remove alpha from the legend
rank_abundance_plot

# Save plot
save_plot(paste0(outPath, "/rank_abundance_plot_B1.png"), rank_abundance_plot, base_width = 15, base_height = 10)

# Donor Abundance Plot - plot to visualize the distribution of donors 
# Calculate the number of unique asv values
unique_donor_asv_count <- beyonce_filtered %>%
  filter(community_mixture == "donor_only") %>%
  distinct(ASVnum) %>%
  nrow()
# Stacked bar plot
donor_abundance_plot <- beyonce_filtered %>% 
  filter(community_mixture == "donor_only") %>% 
  ggplot() +
  geom_bar(aes(x = community_mixture, y = relAbundance, fill = factor(Family)),
           color = "black", stat = "identity") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Relative Abundance Distribution of Donor Community") +
  geom_text(aes(x = 1, y = -0.1, label = paste(unique_donor_asv_count, "Species")),
            vjust = -1, hjust = 0.5, size = 4) 
donor_abundance_plot
# Save plot
save_plot(paste0(outPath, "/donor_abundance_plot_A1_v2.png"), donor_abundance_plot, base_width = 15, base_height = 10)


# Recipient Abundance Plot - plot to visualize the distribution of recipients 
# Calculate the number of unique asv values
unique_recipient_asv_count <- beyonce_filtered %>%
  filter(community_mixture == "recipient_only") %>%
  distinct(ASVnum) %>%
  nrow()
# Stacked bar plot 
recipient_abundance_plot <- beyonce_filtered %>% 
  filter(community_mixture == "recipient_only") %>% 
  ggplot() +
  geom_bar(aes(x = community_mixture, y = relAbundance, fill = factor(Family)),
           color = "black", stat = "identity") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Relative Abundance Distribution of Recipient Community") +
  geom_text(aes(x = 1, y = -0.1, label = paste(unique_recipient_asv_count, "Species")),
            vjust = -1, hjust = 0.5, size = 4) 
recipient_abundance_plot
# Save plot
save_plot(paste0(outPath, "/recipient_abundance_plot_B12.png"), recipient_abundance_plot, base_width = 15, base_height = 10)


# Full opacity plot
# Calculate the number of unique asv values
unique_mixture_asv_count <- beyonce_filtered %>%
  filter(community_mixture == "donor+recipient") %>%
  distinct(ASVnum) %>%
  nrow()
# Stacked Bar Plot
full_opacity_plot <- beyonce_B1_filtered %>% 
  filter(community_mixture == "donor+recipient") %>% 
  ggplot() +
  geom_bar(aes(x = community_mixture, y = relAbundance, fill = factor(Family)),
           color = "black", stat = "identity") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Mixture") +
  geom_text(aes(x = 1, y = -0.1, label = paste(unique_mixture_asv_count, "Species")),
            vjust = -1, hjust = 0.5, size = 4) +
  facet_wrap(~ well + subject, ncol = 4) +
  guides(fill = "none")

# Modulated opacity plot
modulated_opacity_plot <- beyonce_B1_filtered %>% 
  filter(community_mixture == "donor+recipient") %>% 
  ggplot() +
  geom_bar(aes(x = community_mixture, y = relAbundance, fill = factor(Family),
               alpha = ifelse(colonization_status == "Colonizer_Donor", 0.5, 1)),
           color = "black", stat = "identity") +
  scale_fill_manual(values = KCHpalettee0041vector) +
  xlab("Community Mixture Type") +
  ylab("Relative Abundance") +
  ggtitle("Colonizers") +
  facet_wrap(~ well + subject, ncol = 4) +
  guides(alpha = "none")  # To remove alpha from the legend

# Arrange plots side by side using gridExtra
arranged_colonizer_abundance_stacked <- grid.arrange(full_opacity_plot, modulated_opacity_plot, ncol = 2)

# Save plot
save_plot(paste0(outPath, "/arranged_colonizer_abundance_stacked_V2.png"), arranged_colonizer_abundance_stacked, base_width = 20, base_height = 10)
ggsave(
  filename = paste0(outPath, "/arranged_colonizer_abundance_stacked_V2.png"),
  plot = arranged_colonizer_abundance_stacked,
  width = 20,
  height = 10
)

# Scatter Plot



# Shared Colonizers - a plot to compare the overlap of colonizers between different recipient 
# communities to show which colonizers are consistently successful across various communities
# bar plot - # of times each ASV appears in the list of successful colonizers, compared to the 
# of times it appears
# could also try this at the family level

data_shared_colonizers <- beyonce_B1_filtered %>%
  filter(colonization_status == "Colonizer_Donor") %>%
  group_by(ASVnum, community_mixture) %>%
  summarize(count = n()) %>%
  ungroup()

# Create the Shared Colonizers plot
shared_colonizers_plot <- data_shared_colonizers %>% 
  ggplot() +
  geom_bar(aes(x = ASVnum, y = count, fill = community_mixture), stat = "identity", position = "dodge") +
  xlab("ASV") +
  ylab("Count of Appearances") +
  ggtitle("Shared Colonizers - Overlap of Successful Colonizers") +
  theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank())
shared_colonizers_plot




# Differential Abundance Analysis - identify ASVs that are significantly enriched or depleted 
# in the recipient community after colonization
# maybe we can start by ID-ing ASVs in the recipient communities that disappear in the 
# community mixtures? do they consistently disappear?


# Identify ASVs that disappear going from preAbx to postAbxV1
disappeared_asvs_V1 <- beyonce_B1_filtered %>% anti_join(recipient == "XEA-preAbx", recipient == "XEA-postAbxV1", by = "ASVnum")

# Identify ASVs that disappear going from preAbx to postAbxV2
disappeared_asvs_V2 <- beyonce_B1_filtered %>% anti_join(recipient == "XEA-preAbx", recipient == "XEA-postAbxV2", by = "ASVnum")









# Indicator Species Analysis - identify indicator species that are strongly associated with 
# successful colonization in different recipient communities
# maybe we can also try annotating each species in the donor communities and seeing if they 
# can colonize each of the three recipient communities?

