#########################################
### Results section comparing resident
### sewer and human fecal ASVs.
### Lou LaMartina, finalized Jan 31, 2019
#########################################



setwd("~/Desktop/TimeSeries_final")
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(ggalluvial)
library(reshape2)


# load data (refer to TimeSeries_DataPrep.R)
TimeSeries_object <- readRDS("./RData/TimeSeries_phyloseq_object.RData")
Cities_object <- readRDS("./RData/Cities_phyloseq_object.RData")
Neighborhood_object <- readRDS("./RData/Neighborhood_phyloseq_object.RData")
HMP_object <- readRDS("./RData/HMP/HMP_phyloseq_object.RData")


# subset to stool
Stool_counts_object <- subset_samples(HMP_object, Biome == "Stool")


# convert to relative abundance
Stool_relabun_object <- transform_sample_counts(Stool_counts_object, function(x) x / sum(x))


# human-associated ASVs (refer to TimeSeries_threshold.R)
Final_human_ASVs <- readRDS("./RData/Final_human_ASVs.RData")
Greatest_biome_sources <- readRDS("./RData/HMP/Greatest_biome_sources.RData")


# convert to relative abundance
TimeSeries_relabun_object <- transform_sample_counts(TimeSeries_object, function(x) x / sum(x))
Cities_relabun_object <- transform_sample_counts(Cities_object, function(x) x / sum(x))
Neighborhood_relabun_object <- transform_sample_counts(Neighborhood_object, function(x) x / sum(x))


# extract data
TimeSeries_counts <- data.frame(TimeSeries_object@otu_table@.Data)
Cities_counts <- data.frame(Cities_object@otu_table@.Data)
Neighborhood_counts <- data.frame(Neighborhood_object@otu_table@.Data)
Stool_counts <- data.frame(Stool_counts_object@otu_table@.Data)
Taxonomy_all <- read.csv("./RData/Taxonomy_all.csv")
HMP_tax <- data.frame(HMP_object@tax_table@.Data)
HMP_tax$SeqID <- rownames(HMP_tax)




#############################
### proportions of genera ###
#############################

# proportion Bacteroides
sum(TimeSeries_counts[, colnames(TimeSeries_counts) %in% subset(Taxonomy_all, Genus == "Bacteroides")$ASV]) / sum(TimeSeries_counts) * 100
# [1] 3.374705

sum(Cities_counts[, colnames(Cities_counts) %in% subset(Taxonomy_all, Genus == "Bacteroides")$ASV]) / sum(Cities_counts) * 100
# [1] 6.946048

sum(Neighborhood_counts[, colnames(Neighborhood_counts) %in% subset(Taxonomy_all, Genus == "Bacteroides")$ASV]) / sum(Neighborhood_counts) * 100
# [1] 11.44454

sum(Stool_counts[rownames(Stool_counts) %in% subset(HMP_tax, Genus == "Bacteroides")$SeqID, ]) / sum(Stool_counts) * 100
# [1] 53.27712


# proportion Acinetobacter
sum(TimeSeries_counts[, colnames(TimeSeries_counts) %in% subset(Taxonomy_all, Genus == "Acinetobacter")$ASV]) / sum(TimeSeries_counts) * 100
# [1] 10.76926

sum(Cities_counts[, colnames(Cities_counts) %in% subset(Taxonomy_all, Genus == "Acinetobacter")$ASV]) / sum(Cities_counts) * 100
# [1] 8.760804

sum(Neighborhood_counts[, colnames(Neighborhood_counts) %in% subset(Taxonomy_all, Genus == "Acinetobacter")$ASV]) / sum(Neighborhood_counts) * 100
# [1] 5.280948

sum(Stool_counts[rownames(Stool_counts) %in% subset(HMP_tax, Genus == "Acinetobacter")$SeqID, ]) / sum(Stool_counts) * 100
# [1] 0.004532315




####################################
### proportions of genera in HMP ###
####################################

# combine to genus level
Stool_genus_glom <- speedyseq::tax_glom(Stool_relabun_object, "Genus")
TimeSeries_genus_glom <- speedyseq::tax_glom(TimeSeries_relabun_object, "Genus")
Cities_genus_glom <- speedyseq::tax_glom(Cities_relabun_object, "Genus")
Neighborhood_genus_glom <- speedyseq::tax_glom(Neighborhood_relabun_object, "Genus")


# extract most abundant of each
Stool_top_genus <- prune_taxa(names(sort(taxa_sums(Stool_genus_glom), decreasing = TRUE))[1:5], 
                              Stool_genus_glom)
Stool_top_genera <- as.character(tax_table(Stool_top_genus)[,6])

TimeSeries_top_genus <- prune_taxa(names(sort(taxa_sums(TimeSeries_genus_glom), decreasing = TRUE))[1:5], 
                                   TimeSeries_genus_glom)
TimeSeries_top_genera <- as.character(tax_table(TimeSeries_top_genus)[,6])

Cities_top_genus <- prune_taxa(names(sort(taxa_sums(Cities_genus_glom), decreasing = TRUE))[1:5], 
                               Cities_genus_glom)
Cities_top_genera <- as.character(tax_table(Cities_top_genus)[,6])

Neighborhood_top_genus <- prune_taxa(names(sort(taxa_sums(Neighborhood_genus_glom), decreasing = TRUE))[1:5], 
                                     Neighborhood_genus_glom)
Neighborhood_top_genera <- as.character(tax_table(Neighborhood_top_genus)[,6])


# what are they?
top_genera <- unique(c(tax_table(Stool_top_genus)[,6], tax_table(TimeSeries_top_genus)[,6], 
                       tax_table(Cities_top_genus)[,6], tax_table(Neighborhood_top_genus)[,6]))


# extract top genera from genus glom objects
Stool_top_genus <- subset_taxa(Stool_genus_glom, Genus %in% top_genera)
TimeSeries_top_genus <- subset_taxa(TimeSeries_genus_glom, Genus %in% top_genera)
Cities_top_genus <- subset_taxa(Cities_genus_glom, Genus %in% top_genera)
Neighborhood_top_genus <- subset_taxa(Neighborhood_genus_glom, Genus %in% top_genera)


# get means of relative abundances and genus names
# stool
Stool_top_genus.df <- data.frame(Genus_avg = rowSums(as.matrix(Stool_top_genus@otu_table@.Data)) / 
                                   ncol(as.matrix(Stool_top_genus@otu_table@.Data)))
Stool_top_genus.df$SeqID <- rownames(Stool_top_genus.df)
Stool_top_genus.df <- merge(Stool_top_genus.df, HMP_tax[6:7], by = "SeqID")
Stool_top_genus.df$Source <- "HMP"


# time series
TimeSeries_top_genus.df <- data.frame(Genus_avg = colSums(as.matrix(TimeSeries_top_genus@otu_table@.Data)) / 
                                        nrow(as.matrix(TimeSeries_top_genus@otu_table@.Data)))
TimeSeries_top_genus.df$ASV <- rownames(TimeSeries_top_genus.df)
TimeSeries_top_genus.df <- merge(TimeSeries_top_genus.df, Taxonomy_all[c(1,7)], by = "ASV")
TimeSeries_top_genus.df$Source <- "TimeSeries"


# cities
Cities_top_genus.df <- data.frame(Genus_avg = colSums(as.matrix(Cities_top_genus@otu_table@.Data)) / 
                                    nrow(as.matrix(Cities_top_genus@otu_table@.Data)))
Cities_top_genus.df$ASV <- rownames(Cities_top_genus.df)
Cities_top_genus.df <- merge(Cities_top_genus.df, Taxonomy_all[c(1,7)], by = "ASV")
Cities_top_genus.df$Source <- "Cities"


# neighborhoods
Neighborhood_top_genus.df <- data.frame(Genus_avg = colSums(as.matrix(Neighborhood_top_genus@otu_table@.Data)) / 
                                          nrow(as.matrix(Neighborhood_top_genus@otu_table@.Data)))
Neighborhood_top_genus.df$ASV <- rownames(Neighborhood_top_genus.df)
Neighborhood_top_genus.df <- merge(Neighborhood_top_genus.df, Taxonomy_all[c(1,7)], by = "ASV")
Neighborhood_top_genus.df$Source <- "Neighborhoods"


# convert genus average to relative abundance of averages (want them to be proportions that = 1)
Stool_top_genus.df$Rel_avg <- Stool_top_genus.df$Genus_avg / sum(Stool_top_genus.df$Genus_avg)
TimeSeries_top_genus.df$Rel_avg <- TimeSeries_top_genus.df$Genus_avg / sum(TimeSeries_top_genus.df$Genus_avg)
Cities_top_genus.df$Rel_avg <- Cities_top_genus.df$Genus_avg / sum(Cities_top_genus.df$Genus_avg)
Neighborhood_top_genus.df$Rel_avg <- Neighborhood_top_genus.df$Genus_avg / sum(Neighborhood_top_genus.df$Genus_avg)


# combine
Top_genera.df <- rbind(Stool_top_genus.df[-1], TimeSeries_top_genus.df[-1], Cities_top_genus.df[-1], Neighborhood_top_genus.df[-1])


# colors
bar_colors <- colorRampPalette(brewer.pal(11, "Spectral"))
plot_colors <- bar_colors(length(top_genera))
names(plot_colors) <- sort(top_genera)


### figure 2A ###
bars <- 
  ggplot(Top_genera.df, aes(x = Source, y = Rel_avg, stratum = Genus, alluvium = Genus, label = Genus, fill = Genus)) +
  geom_flow(stat = "alluvium", lode.guidance = "rightleft", alpha = 0.3) +
  geom_stratum(color = NA, width = 0.5) +
  theme_classic() +
  scale_x_discrete(limits = c("HMP", "Neighborhoods", "TimeSeries", "Cities"),
                   labels = c("Human\nstool", "Residential\nsewers", "Milwaukee\ninfluent", "Average\nU.S. influent")) +
  scale_fill_manual(values = plot_colors, breaks = sort(top_genera)) +
  theme(axis.text.x = element_text(size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black", face = "bold"),
        axis.title.x = element_blank(),
        axis.line = element_line(size = 0.25),
        axis.ticks = element_line(size = 0.25),
        panel.border = element_rect(size = 0.25, color = "grey80", fill = NA),
        legend.text = element_text(size = 5, color = "black", face = "italic"),
        legend.title = element_text(size = 6, color = "black", face = "bold")) +
  guides(fill = guide_legend(keyheight = 0.1, keywidth = 0.4, units = "in", ncol = 1)) +
  labs(y = "Relative proportions of top genera\nin each dataset", x = "Source", fill = "Genus")
bars

#ggsave("./Plots/bars.pdf", plot = bars, device = "pdf", width = 3.5, height = 2.0, units = "in")




##########################
### top genus dot plot ###
##########################

# combine WWTP datasets
WWTP_relabun_object <- merge_phyloseq(Cities_relabun_object, TimeSeries_relabun_object)


# glom to genus level
WWTP_genus_glom <- speedyseq::tax_glom(WWTP_relabun_object, "Genus")


# subset to most abundant genera from previous step
WWTP_top_genus <- subset_taxa(WWTP_genus_glom, Genus %in% top_genera)


# extract data
WWTP_top_genus.df <- data.frame(WWTP_top_genus@otu_table@.Data)


# change column names to genus
colnames(WWTP_top_genus.df) <- tax_table(WWTP_top_genus)[,6]


# melt for plotting
WWTP_top_genus.m <- melt(WWTP_top_genus.df, variable.name = "Genus", value.name = "Relabun")


# get stats
WWTPmeans <- aggregate(. ~ Genus, mean, data = WWTP_top_genus.m)


### figure 2B ###
genus.dots <- 
  ggplot(WWTP_top_genus.m, aes(x = Genus, y = Relabun, color = Genus)) +
  geom_jitter(size = 1, alpha = 0.5, shape = 1, width = 0.2) +
  theme_classic() +
  scale_x_discrete(limits = names(sort(colSums(WWTP_top_genus.df), decreasing = TRUE))) +
  scale_color_manual(values = plot_colors) +
  theme(axis.text.x = element_text(size = 6, color = "black", face = "italic", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black", face = "bold"),
        axis.title.x = element_blank(),
        panel.border = element_rect(size = 0.25, color = "grey80", fill = NA),
        axis.line = element_line(size = 0.25),
        axis.ticks = element_line(size = 0.25),
        legend.position = "none") +
  geom_crossbar(data = WWTPmeans, aes(ymin = Relabun, ymax = Relabun), size = 0.1, color = "black") +
  labs(y = "Proportions of top genera\nin WWTP influent", x = "Genus", fill = "Genus")
genus.dots

#ggsave("./Plots/genusdots.pdf", plot = genus.dots, device = "pdf", width = 3.5, height = 2, units = "in")




###########################################################
### proportions of human-associated reads in wastewater ###
###########################################################

# find proportions of reads in each dataset
TimeSeries_human_props <- data.frame(Sample_name = rownames(TimeSeries_counts),
                                     Prop = rowSums(TimeSeries_counts[, colnames(TimeSeries_counts) %in% Final_human_ASVs]) / rowSums(TimeSeries_counts),
                                     Source = "TimeSeries")

Cities_human_props <- data.frame(Sample_name = rownames(Cities_counts),
                                 Prop = rowSums(Cities_counts[, colnames(Cities_counts) %in% Final_human_ASVs]) / rowSums(Cities_counts),
                                 Source = "Cities")

Neighborhood_human_props <- data.frame(Sample_name = rownames(Neighborhood_counts),
                                       Prop = rowSums(Neighborhood_counts[, colnames(Neighborhood_counts) %in% Final_human_ASVs]) / rowSums(Neighborhood_counts),
                                       Source = "Neighborhoods")


# combine
Human_props <- rbind(Neighborhood_human_props, TimeSeries_human_props, Cities_human_props)


# get stats
propmeans <- aggregate(. ~ Source, mean, data = Human_props)


### figure 2C ###
props <- 
  ggplot(Human_props, aes(x = Source, y = Prop, color = Source)) +
  geom_jitter(size = 1, alpha = 0.5, shape = 1, width = 0.2) +
  theme_classic() +
  scale_color_manual(values = brewer.pal(9, "Blues")[c(9, 5, 7)]) +
  scale_x_discrete(labels = c("Residential\nsewers", "Average\nU.S. influent", "Milwaukee\ninfluent"),
                   limits = c("Neighborhoods", "Cities", "TimeSeries")) +
  theme(axis.text.x = element_text(size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6, color = "black", face = "bold"),
        legend.position = "none",
        panel.border = element_rect(size = 0.25, color = "grey80", fill = NA),
        axis.line = element_line(size = 0.25),
        axis.ticks = element_line(size = 0.25)) +
  geom_crossbar(data = propmeans, aes(ymin = Prop, ymax = Prop), size = 0.1, color = "black") +
  labs(y = "Proportions of human-associated ASVs\nin raw wastewater", x = "Source",
       color = "Mean\npercents")
props

#ggsave("./Plots/props.pdf", plot = props, device = "pdf", width = 2.2, height = 2.0, units = "in")




######################################################
### proportions of source body sites in wastewater ###
######################################################

# define vector you'll be using
Biomes <- c("Stool", "Skin", "Oral", "Vaginal")


# transpose so ASVs are rows
WWTP_counts <- data.frame(merge_phyloseq(TimeSeries_object, Cities_object)@otu_table@.Data)
Proportion_WWTP <- data.frame(t(WWTP_counts))


# add greatest source for each ASV
Proportion_WWTP$Source <- NA
for(i in Biomes){
  Proportion_WWTP$Source[rownames(Proportion_WWTP) %in% subset(Greatest_biome_sources, Greatest_source == i)$ASV] <- i
  Proportion_WWTP$Source[is.na(Proportion_WWTP$Source)] <- "Sewer"
}


# split them up by source, get rid of source column, convert to data matrices
Proportions_in_WWTP_dms.ls <- list()
for(i in Biomes){
  Proportions_in_WWTP_dms.ls[[i]] <- as.matrix(subset(Proportion_WWTP, Source == i)[,-ncol(Proportion_WWTP)])
}
names(Proportions_in_WWTP_dms.ls) <- Biomes

Proportions_in_WWTP_dms.ls[["Sewer"]] <- as.matrix(subset(Proportion_WWTP, Source == "Sewer")[,-ncol(Proportion_WWTP)])


# get total number of human-associated reads in each sample, then divide by the total number of reads
WWTP_meanreads.ls <- list()
for(i in Biomes){
  WWTP_meanreads.ls[[i]] <- colSums(Proportions_in_WWTP_dms.ls[[i]]) / rowSums(data.matrix(WWTP_counts))
}

WWTP_meanreads.ls[["Sewer"]] <- colSums(Proportions_in_WWTP_dms.ls[["Sewer"]]) / rowSums(data.matrix(WWTP_counts))


# combine them
Mean_sources_inWWTP <- data.frame(Stool = WWTP_meanreads.ls[["Stool"]], Skin = WWTP_meanreads.ls[["Skin"]],
                                  Oral = WWTP_meanreads.ls[["Oral"]], Vaginal = WWTP_meanreads.ls[["Vaginal"]],
                                  Sewer = WWTP_meanreads.ls[["Sewer"]])


# R adds X in front of cities names because they start with a number. change them back
rownames(Mean_sources_inWWTP) <- rownames(WWTP_counts)


# melt
Mean_sources_inWWTP.m <- melt(Mean_sources_inWWTP, variable.name = "Source", value.name = "Proportion")

means <- aggregate(. ~ Source, mean, data = Mean_sources_inWWTP.m)
sd <- aggregate(. ~ Source, sd, data = Mean_sources_inWWTP.m)


### figure 2D ###
dots <- 
  ggplot(Mean_sources_inWWTP.m, aes(x = Source, y = Proportion, color = Source)) +
  geom_jitter(size = 1.5, alpha = 0.5, shape = 1) +
  theme_classic() +
  scale_x_discrete(limits = c("Sewer", "Stool", "Skin", "Oral", "Vaginal")) +
  scale_color_manual(values = brewer.pal(9, "Blues")[c(8:5,9)]) +
  theme(axis.text.x = element_text(size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black", face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(size = 0.25, color = "grey80", fill = NA),
        axis.line = element_line(size = 0.25),
        axis.ticks = element_line(size = 0.25)) +
  geom_crossbar(data = means, aes(ymin = Proportion, ymax = Proportion), size = 0.2, color = "black") +
  labs(y = "Proportions of body site associations\nof WWTP influent ASVs")
dots

#ggsave("./Plots/dots.pdf", plot = dots, device = "pdf", width = 2.25, height = 1.9, units = "in")


# save
#save.image("./RData/Results2_env.RData")


