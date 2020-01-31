#############################################
### Threshold to distinguish human-associated
### vs. sewer-associated ASVs
### Lou LaMartina, finalized Jan 30, 2020
#############################################


###############
### PURPOSE ###
###############

# We aligned 16S rRNA v4v5 sequences from wastewater treatment plant 
# influent datasets (71 cities, 5-year Milwaukee time series, and 
# residential sewers) to those from the Human Microbiome Project. 
#
# Any othat were 100% matches were considered derived from humans. 
#
# However, many potential "environmental" organisms
# occurred in the human microbiome.
#
# So we designed a threshold to filter them out:
#
# Among ASVs that were shared between WWTP influnet and HMP:
#
# if the minimum abundance (5-percentile) of a wastewater ASV 
# exceeds the maximum abundance (95-percentile) of that ASV in the human microbiome,
# it was reclassified as a sewer-associated sequence.
#
# min(wastewater) > max(HMP) -> sewer ASV
#
# if the minimum abundance (5-percentile) of a wastewater ASV 
# is less than the maximum abundance (95-percentile) of that ASV in the human microbiome,
# it was considered a human-associated sequence.
#
# min(wastewater) < max(HMP) -> human ASV


setwd("~/Desktop/TimeSeries_final")
library(dplyr)
library(reshape2)
library(phyloseq)
library(matrixStats)




#################
### load data ###
#################

# (shortcut to threshold section)
#load("./RData/HMP/TimeSeries_threshold_data.RData")


# phyloseq objects
HMP_object <- readRDS("./RData/HMP/HMP_phyloseq_object.RData") # 17.38879 secs
TimeSeries_object <- readRDS("./RData/TimeSeries_phyloseq_object.RData")
Cities_object <- readRDS("./RData/Cities_phyloseq_object.RData")


# combine all wastwater samples into one PS object
WWTP_object <- merge_phyloseq(TimeSeries_object, Cities_object)


# load results of alignment - wastewater ASVs that perfectly matched HMP sequences
HMP_align <- read.csv("./RData/HMP/HMP_alignments.csv")


# subset HMP align so it only has ASVs that are in the HMP object
HMP_align <- HMP_align[HMP_align$SeqID %in% taxa_names(HMP_object),]

# extract taxonomy table
Taxonomy_all <- data.frame(WWTP_object@tax_table@.Data)
Taxonomy_all$ASV <- rownames(Taxonomy_all)


# create biomes variable
Biomes <- c("Stool", "Skin", "Oral", "Vaginal")




#######################################
### subset phyloseq object by biome ###
#######################################

# subset by (simplified) body sites
Stool_object <- subset_samples(HMP_object, Biome == "Stool")
Skin_object <- subset_samples(HMP_object, Biome == "Skin")
Oral_object <- subset_samples(HMP_object, Biome == "Oral")
Vaginal_object <- subset_samples(HMP_object, Biome == "Vaginal")


# remove original to save memory
rm(HMP_object)


# remove ASVs that did not occur in these samples, or samples with no ASVs
Stool_object <- subset_taxa(Stool_object, taxa_sums(Stool_object) > 0)
Skin_object <- subset_taxa(Skin_object, taxa_sums(Skin_object) > 0)
Oral_object <- subset_taxa(Oral_object, taxa_sums(Oral_object) > 0)
Vaginal_object <- subset_taxa(Vaginal_object, taxa_sums(Vaginal_object) > 0)


# create a list of HMP objects
HMP_objects.ls <- list(Stool_object, Skin_object, Oral_object, Vaginal_object)
names(HMP_objects.ls) <- Biomes


# convert to relative abundances (takes a minute)
HMP_relabun_objects.ls <- lapply(HMP_objects.ls, function(i) 
  transform_sample_counts(i, function(x) x / sum(x)))

WWTP_relabun_object <- transform_sample_counts(WWTP_object, function(x) x / sum(x))


# subset to only sequences that were found in influent
# (for some reason subset_taxa doesn't like working with lapply)
Stool_relabun_object_aligned <- subset_taxa(HMP_relabun_objects.ls[["Stool"]], 
                                            taxa_names(HMP_relabun_objects.ls[["Stool"]]) %in% HMP_align$SeqID)

Skin_relabun_object_aligned <- subset_taxa(HMP_relabun_objects.ls[["Skin"]], 
                                            taxa_names(HMP_relabun_objects.ls[["Skin"]]) %in% HMP_align$SeqID)

Oral_relabun_object_aligned <- subset_taxa(HMP_relabun_objects.ls[["Oral"]], 
                                            taxa_names(HMP_relabun_objects.ls[["Oral"]]) %in% HMP_align$SeqID)

Vaginal_relabun_object_aligned <- subset_taxa(HMP_relabun_objects.ls[["Vaginal"]], 
                                            taxa_names(HMP_relabun_objects.ls[["Vaginal"]]) %in% HMP_align$SeqID)


# make new list
HMP_relabuns_aligned_objects.ls <- list(Stool_relabun_object_aligned, Skin_relabun_object_aligned, 
                                  Oral_relabun_object_aligned, Vaginal_relabun_object_aligned)

names(HMP_relabuns_aligned_objects.ls) <- Biomes




########################
### define threshold ###
########################


# find mean relative abundances of unique sequences in HMP
HMP_aligned_means.ls <- lapply(HMP_relabuns_aligned_objects.ls, function(i) 
  taxa_sums(i) / length(sample_names(i)))


# change HMP sequence IDs to the ASV IDs that they matched to.
# doing this a long convoluted way to make sure i'm not getting them mixed up:
# 1. convert to data frame, 2. make SeqID variable, 3. merge with HMP_align to add ASV variable, 
# 4. change row names to those ASVs, 5. remove seqID and ASV columns
for(i in Biomes){
  HMP_aligned_means.ls[[i]] <- data.frame(Mean_in_HMP = HMP_aligned_means.ls[[i]])
  HMP_aligned_means.ls[[i]]$SeqID <- rownames(HMP_aligned_means.ls[[i]])
  HMP_aligned_means.ls[[i]] <- merge(HMP_aligned_means.ls[[i]], HMP_align, by = "SeqID")
  rownames(HMP_aligned_means.ls[[i]]) <- HMP_aligned_means.ls[[i]]$ASV
  HMP_aligned_means.ls[[i]] <- HMP_aligned_means.ls[[i]][-1]
}


# mean relative abundances in WWTP
WWTP_means <- data.frame(Mean_in_WWTP = taxa_sums(WWTP_relabun_object) / length(sample_names(WWTP_relabun_object)),
                         ASV = taxa_names(WWTP_relabun_object))


# extract abundance matrices into data frames
HMP_relabuns_aligned_dfs.ls <- lapply(HMP_relabuns_aligned_objects.ls, function(i)
  data.frame(i@otu_table@.Data))

WWTP_relabun <- data.frame(WWTP_relabun_object@otu_table@.Data)


# change HMP sequence IDs to the ASV IDs that they matched to (same as above)
for(i in Biomes){
  HMP_relabuns_aligned_dfs.ls[[i]]$SeqID <- rownames(HMP_relabuns_aligned_dfs.ls[[i]])
  HMP_relabuns_aligned_dfs.ls[[i]] <- merge(HMP_relabuns_aligned_dfs.ls[[i]], HMP_align, by = "SeqID")
  rownames(HMP_relabuns_aligned_dfs.ls[[i]]) <- HMP_relabuns_aligned_dfs.ls[[i]]$ASV
  HMP_relabuns_aligned_dfs.ls[[i]] <- HMP_relabuns_aligned_dfs.ls[[i]][-c(1,ncol(HMP_relabuns_aligned_dfs.ls[[i]]))]
}


# convert to data matrices
HMP_relabuns_aligned_dms.ls <- lapply(HMP_relabuns_aligned_dfs.ls, function(i) as.matrix(i))


# find 95-percentiles of human sequences in microbiome samples
HMP_95pct.ls <- lapply(HMP_relabuns_aligned_dms.ls, function(i)
  rowQuantiles(i, probs = 0.95))


# find 5-percentiles of all ASVs in wastewater datasets
WWTP_5pct <- rowQuantiles(t(WWTP_relabun), probs = 0.05)


# extract shared values: i only want to see the percentiles of ASVs that are in both wastewater and human microbiome
Shared_5pct_ASVs.ls <- list()
for(i in Biomes){
  Shared_5pct_ASVs.ls[[i]] <- data.frame(WWTP = subset(WWTP_5pct, names(WWTP_5pct) %in% names(HMP_95pct.ls[[i]])),
                                         ASV = names(subset(WWTP_5pct, names(WWTP_5pct) %in% names(HMP_95pct.ls[[i]]))))
}
names(Shared_5pct_ASVs.ls) <- names(HMP_95pct.ls)


# subset to only matched ASV
for(i in Biomes){
  HMP_95pct.ls[[i]] <- data.frame(HMP = subset(HMP_95pct.ls[[i]], names(HMP_95pct.ls[[i]]) %in% rownames(Shared_5pct_ASVs.ls[[i]])),
                                  ASV = names(subset(HMP_95pct.ls[[i]], names(HMP_95pct.ls[[i]]) %in% rownames(Shared_5pct_ASVs.ls[[i]]))))
}


# combine data frames
Shared_5pct_ASVs_dfs.ls <- list()
for(i in Biomes){
  Shared_5pct_ASVs_dfs.ls[[i]] <- merge(Shared_5pct_ASVs.ls[[i]], HMP_95pct.ls[[i]], by = "ASV")
}
names(Shared_5pct_ASVs_dfs.ls) <- names(HMP_95pct.ls)


# 1. show which have min wastewater abundance that's greater than max HMP abundance,
# 2. subset to ones that were true
for(i in Biomes){
  Shared_5pct_ASVs_dfs.ls[[i]] <- data.frame(Shared_5pct_ASVs_dfs.ls[[i]], minWW_gt_maxHMP =
                                        as.logical((Shared_5pct_ASVs_dfs.ls[[i]][,2] > Shared_5pct_ASVs_dfs.ls[[i]][,3])))
  Shared_5pct_ASVs_dfs.ls[[i]] <- subset(Shared_5pct_ASVs_dfs.ls[[i]], Shared_5pct_ASVs_dfs.ls[[i]][,4] == TRUE)
}


# add HMP and wastewater means, and taxonomy
for(i in Biomes){
  Shared_5pct_ASVs_dfs.ls[[i]] <- merge(merge(merge(Shared_5pct_ASVs_dfs.ls[[i]], HMP_aligned_means.ls[[i]], by = "ASV"),
                                        WWTP_means, by = "ASV"),
                                        Taxonomy_all, by = "ASV")
}


# which are they?
ReclassASVs_95pct <- unique(c(as.character(Shared_5pct_ASVs_dfs.ls[[1]]$ASV), as.character(Shared_5pct_ASVs_dfs.ls[[2]]$ASV),
                       as.character(Shared_5pct_ASVs_dfs.ls[[3]]$ASV), as.character(Shared_5pct_ASVs_dfs.ls[[4]]$ASV)))




###############################
### threshold within biomes ###
###############################

# need to make sure that just because an ASV passes the threshold for one body site, it doesn't for another.
# for example, Blautia is low in in vaginal microbiome, but shouldn't be kicked out because it has normal
# abundance in stool. so pull out 95pctiles of all of reclassified ASVs!

# subset 95pctile data to only those that were reclassified from "human" to "sewer", add ASV variable
Reclassed_95pcts.ls <- list()
for(i in Biomes){
  Reclassed_95pcts.ls[[i]] <- data.frame(HMP_95pct = subset(HMP_95pct.ls[[i]], ASV = rownames(HMP_95pct.ls[[i]]) %in% ReclassASVs_95pct))
  colnames(Reclassed_95pcts.ls[[i]]) <- c("HMP_95pct", "ASV")
}
names(Reclassed_95pcts.ls) <- Biomes


# same for their abundance in WWTP
WWTP_5pct_inReclass <- data.frame(WWTP = subset(WWTP_5pct, names(WWTP_5pct) %in% ReclassASVs_95pct),
                                  ASV = names(subset(WWTP_5pct, names(WWTP_5pct) %in% ReclassASVs_95pct)))

WWTP_means_inReclass <- data.frame(WWTP = subset(WWTP_means, rownames(WWTP_means) %in% ReclassASVs_95pct))


# combine them into one data frame
Reclassed_95pcts.df <- merge(merge(merge(Reclassed_95pcts.ls[["Stool"]], Reclassed_95pcts.ls[["Skin"]], by = "ASV", all = TRUE),
                             Reclassed_95pcts.ls[["Oral"]], by = "ASV", all = TRUE),
                             Reclassed_95pcts.ls[["Vaginal"]], by = "ASV", all = TRUE)
colnames(Reclassed_95pcts.df)[2:5] <- Biomes

Reclassed_95pcts.df <- merge(Reclassed_95pcts.df, WWTP_5pct_inReclass, by = "ASV", all = TRUE)


# which sources have highest values?
Reclassed_95pcts.df$Greatest_source <- colnames(Reclassed_95pcts.df)[apply(Reclassed_95pcts.df, 1, which.max)]


# save
Final_human_ASVs <- as.character(droplevels(
  subset(HMP_align$ASV, ! HMP_align$ASV %in% 
           subset(Reclassed_95pcts.df, Greatest_source == "WWTP")$ASV)))
#saveRDS(Final_human_ASVs, "./RData/Final_human_ASVs.RData")




#################################
### designate sources of ASVs ###
#################################

# want each ASV to have one designated source: 
# stool, oral, skin, vaginal, or sewer.
# for human associated ones, i'm comparing the 
# mean relative abundances of final human ASVs
# in the human microbiome. since some ASVs are
# found in >1 body site, whichever body site has
# the greatest abundance of a given ASV, that will
# be its source. 

# subset HMP means from before, keeping only ASVs considered human after threshold
for(i in Biomes){
  HMP_aligned_means.ls[[i]] <- subset(HMP_aligned_means.ls[[i]], rownames(HMP_aligned_means.ls[[i]]) %in% Final_human_ASVs)
}


# combine them in order to figure out which is the most
HMP_means_all <- merge(merge(merge(HMP_aligned_means.ls[["Stool"]], HMP_aligned_means.ls[["Skin"]], by = "ASV", all = TRUE),
                       HMP_aligned_means.ls[["Oral"]], by = "ASV", all = TRUE),
                       HMP_aligned_means.ls[["Vaginal"]], by = "ASV", all = TRUE)
colnames(HMP_means_all)[2:5] <- Biomes


# determine which is the most
HMP_means_all$Greatest_source <- colnames(HMP_means_all)[apply(HMP_means_all, 1, which.max)]


# save
#saveRDS(HMP_means_all, "./RData/HMP/Greatest_biome_sources.RData")
#save.image("./RData/HMP/Threshold_env.RData")
