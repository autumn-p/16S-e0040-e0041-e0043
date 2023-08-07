#Import libraries
library(gridExtra)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)
library(VennDiagram) 

# Set up file paths
outPath = "C:/Users/aparr/16S-e0041-e0043/analysis/out" # out
dataframePath = "C:/Users/aparr/16S-e0041-e0043/data/ps_all.txt.gz" #raw data
appendCol_path = "C:/Users/aparr/16S-e0041-e0043/data/metadatae0041.tsv" #metadata
KCHpalette_path = "C:/Users/aparr/16S-e0041-e0043/config/KCHcolors-Silva-partial.txt" #color palette

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



# Based off naming convention community (ex. "XBB-029") or "blank
alpha_diversity_table <- alpha_diversity_table %>%
  mutate(community_mixture = case_when(
    # super+recipient portion
    donor == "super-community" & recipient != "blank" ~ "super+recipient",
    # recipient_only portion
    donor == "blank" & recipient != "blank" ~ "recipient_only",
    # donor_only, should I be excluding super-community****
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


#Box Plot of Alpha Diversity by Subject:
#Create a box plot to visualize the distribution of alpha diversity for each subject


