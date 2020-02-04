#########################################
### Results section analyzing seasonality
### Lou LaMartina, finalized Feb 4, 2019
#########################################


setwd("~/Desktop/TimeSeries_final")
library(phyloseq)
library(vegan)
library(ggplot2)
library(OTUtable)
library(RColorBrewer)
library(indicspecies)


# load data (refer to TimeSeries1_DataPrep.R)
TimeSeries_object <- readRDS("./RData/TimeSeries_phyloseq_object.RData")
TimeSeries_info <- read.csv("./RData/TimeSeries_sample_info.csv")
TimeSeries_info$Collection_date <- as.Date(TimeSeries_info$Collection_date, format = "%m/%d/%y")
rownames(TimeSeries_info) <- TimeSeries_info$Sample_name
Taxonomy_all <- read.csv("./RData/Taxonomy_all.csv")


# human ASVs (refer to TimeSeries2_Threshold.R)
Final_human_ASVs <- readRDS("./RData/Final_human_ASVs.RData")

# subset the two treatment plants
JI_object <- subset_samples(TimeSeries_object, Treatment_plant == "Jones_Island")
SS_object <- subset_samples(TimeSeries_object, Treatment_plant == "South_Shore")
JI_info <- subset(TimeSeries_info, Treatment_plant == "Jones_Island")
SS_info <- subset(TimeSeries_info, Treatment_plant == "South_Shore")


# convert to relative abundance
TimeSeries_relabun_object <- transform_sample_counts(TimeSeries_object, function(x) x / sum(x))
JI_relabun_object <- transform_sample_counts(JI_object, function(x) x / sum(x))
SS_relabun_object <- transform_sample_counts(SS_object, function(x) x / sum(x))


# extract data
TimeSeries_relabun <- data.frame(TimeSeries_relabun_object@otu_table@.Data)
JI_counts <- data.frame(JI_object@otu_table@.Data)
JI_relabun <- data.frame(JI_relabun_object@otu_table@.Data)
SS_relabun <- data.frame(SS_relabun_object@otu_table@.Data)


# split up human- and sewer-associated ASVs, remove empty ASVs, convert to relative abundances
JI_human_counts <- JI_counts[, colnames(JI_counts) %in% Final_human_ASVs]
JI_human_counts <- JI_human_counts[, colSums(JI_human_counts != 0) > 0]
JI_human_relabun <- JI_human_counts / rowSums(JI_human_counts)

JI_sewer_counts <- JI_counts[, colnames(JI_counts) %in% Final_human_ASVs]
JI_sewer_counts <- JI_sewer_counts[, colSums(JI_sewer_counts != 0) > 0]
JI_sewer_relabun <- JI_sewer_counts / rowSums(JI_sewer_counts)


# remove ASVs with less than 1% as their maximum relative abundance
maxtax <- apply(JI_sewer_relabun, 2, max)
mintax <- names(which(maxtax < 0.01))
JI_sewer_filt <- JI_sewer_relabun[, -which(colnames(JI_sewer_relabun) %in% mintax)]
dim(JI_sewer_filt)
# [1] 60 69

maxtax <- apply(JI_human_relabun, 2, max)
mintax <- names(which(maxtax < 0.01))
JI_human_filt <- JI_human_relabun[, -which(colnames(JI_human_relabun) %in% mintax)]
dim(JI_human_filt)
# [1] 60 67




##################
### ordination ###
##################

# stat
TimeSeries_ord <- ordinate(TimeSeries_relabun_object, method = "PCoA")


# extract values
TimeSeries_PCoA.df <- data.frame(TimeSeries_ord$vectors[,1:2])


# add info
TimeSeries_PCoA.df$Sample_name <- rownames(TimeSeries_PCoA.df)
TimeSeries_PCoA.df <- merge(TimeSeries_PCoA.df, TimeSeries_info[1:9], by = "Sample_name")


### figure 3A ###
months <- 
  ggplot(TimeSeries_PCoA.df, aes(x = Axis.2, y = Axis.1, color = Month, shape = Treatment_plant)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80", size = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80", size = 0.2) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_classic() +
  scale_color_manual(values = c("#ABDDA4", "#FDAE61", "#4D004B", "#3288BD", 
                                "#5E4FA2", "#FEE08B", "#FFFFBF", "#66C2A5", 
                                "#E6F598", "#9E0142", "#D53E4F", "#F46D43"),
                     breaks = c("January","February", "March", "April",
                                "May", "June", "July", "August", "September", 
                                "October", "November", "December")) +
  scale_shape_discrete(labels = c("Jones Island", "South Shore")) +
  theme(axis.text.x = element_text(size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black"),
        axis.title.x = element_text(size = 6, color = "black"),
        legend.text = element_text(size = 5, color = "black"),
        legend.title = element_text(size = 6, color = "black", face = "bold"),
        legend.position = "left",
        axis.line = element_line(size = 0.25),
        axis.ticks = element_line(size = 0.25),
        panel.border = element_rect(color = "grey80", fill = NA, size = 0.25)) +
  guides(color = guide_legend(keyheight = 0.2, keywidth = 0.2, units = "in", ncol = 1),
         shape = guide_legend(keyheight = 0.2, keywidth = 0.2, units = "in", ncol = 1)) +
  labs(y = "Axis 1\n(27.4%)", x = "Axis 2\n(12.5%)", color = "Sample\nmonth",
       shape = "Treatment\nplant")
months

#ggsave("./Plots/months.pdf", plot = months, device = "pdf", width = 3, height = 2.8, units = "in")



### figure 3B ###
line <- 
  ggplot(TimeSeries_PCoA.df, aes(x = Collection_date, y = Axis.1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80", size = 0.2) +
  geom_line(aes(color = Treatment_plant), size = 0.3) +
  geom_point(aes(shape = Treatment_plant, color = Treatment_plant), size = 1.5, alpha = 0.8) +
  theme_classic() +
  scale_color_manual(values = c("grey20", "grey60")) +
  ylim(min(TimeSeries_PCoA.df$Axis.1), max(TimeSeries_PCoA.df$Axis.1)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme(axis.text.x = element_text(size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 0.25),
        axis.ticks = element_line(size = 0.25),
        panel.border = element_rect(color = "grey80", fill = NA, size = 0.25))
line

#ggsave("./Plots/line.pdf", plot = line, device = "pdf", width = 3, height = 2.6, units = "in")


# add season variable
TimeSeries_PCoA.df$Season[TimeSeries_PCoA.df$Month == "February" | TimeSeries_PCoA.df$Month == "March" |
                            TimeSeries_PCoA.df$Month == "April" | TimeSeries_PCoA.df$Month == "May" | 
                            TimeSeries_PCoA.df$Month == "June"] <- "Spring"

TimeSeries_PCoA.df$Season[TimeSeries_PCoA.df$Month == "August" | TimeSeries_PCoA.df$Month == "September" |
                            TimeSeries_PCoA.df$Month == "October" | TimeSeries_PCoA.df$Month == "November" | 
                            TimeSeries_PCoA.df$Month == "December"] <- "Fall"



### figure 3C ###
seasons <- 
  ggplot(na.omit(TimeSeries_PCoA.df), aes(x = Treatment_plant, y = Axis.1, color = Season, fill = Season)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80", size = 0.2) +
  geom_boxplot(width = 0.4, alpha = 0.5, outlier.size = 0.5, size = 0.3, outlier.shape = 1) +
  theme_classic() +
  scale_x_discrete(labels = c("Jones\nIsland", "South\nShore")) +
  ylim(min(TimeSeries_PCoA.df$Axis.1), max(TimeSeries_PCoA.df$Axis.1)) +
  scale_color_manual(values = c("#F46D43", "#66C2A5")) + 
  scale_fill_manual(values = c("#F46D43", "#66C2A5")) +
  theme(axis.text.x = element_text(size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 5, color = "black"),
        legend.title = element_text(size = 6, color = "black", face = "bold"),
        legend.position = c(0.29,0.93),
        legend.background = element_rect(fill = NA),
        axis.line = element_line(size = 0.25),
        axis.ticks = element_line(size = 0.25),
        panel.border = element_rect(color = "grey80", fill = NA, size = 0.25)) +
  guides(color = guide_legend(keyheight = 0.3, keywidth = 0.3, units = "in", ncol = 1))
seasons

#ggsave("./Plots/seasons.pdf", plot = seasons, device = "pdf", width = 1.1, height = 2.7, units = "in")




##########################
### analyze ordination ###
##########################

# which months are in + of Axis 1?
cbind(Pos = data.frame(table(TimeSeries_PCoA.df[TimeSeries_PCoA.df$Axis.1 > 0, "Month"])), 
      Neg = data.frame(table(TimeSeries_PCoA.df[TimeSeries_PCoA.df$Axis.1 <= 0, "Month"])))[-3]
#     Pos.Var1 Pos.Freq Neg.Freq
# 1      April        0        8
# 2     August        8        0
# 3   December        4        5
# 4   February        0        9
# 5    January        0        9
# 6       July        6        2
# 7       June        2        5
# 8      March        0        6
# 9        May        0        7
# 10  November        6        0
# 11   October        9        0
# 12 September        8        0




################################
### determine indicator ASVs ###
################################

# add order, family, genus, species to ASV ID:
# transpose to ASVs are row names
JI_sewer_filt.df.t <- data.frame(t(JI_sewer_filt))
JI_sewer_filt.df.t$ASV <- rownames(JI_sewer_filt.df.t)

JI_human_filt.df.t <- data.frame(t(JI_human_filt))
JI_human_filt.df.t$ASV <- rownames(JI_human_filt.df.t)


# merge with tax info
JI_sewer_filt.df.t <- merge(Taxonomy_all[-8], JI_sewer_filt.df.t, by = "ASV")
JI_human_filt.df.t <- merge(Taxonomy_all[-8], JI_human_filt.df.t, by = "ASV")


# combine into one variable
JI_sewer_filt.df.t$Names <- paste0(JI_sewer_filt.df.t$Order, "__", JI_sewer_filt.df.t$Family, "__", 
                                   JI_sewer_filt.df.t$Genus, "__", JI_sewer_filt.df.t$ASV)

JI_human_filt.df.t$Names <- paste0(JI_human_filt.df.t$Order, "__", JI_human_filt.df.t$Family, "__", 
                                   JI_human_filt.df.t$Genus, "__", JI_human_filt.df.t$ASV)


# make that row names
rownames(JI_sewer_filt.df.t) <- JI_sewer_filt.df.t$Names
rownames(JI_human_filt.df.t) <- JI_human_filt.df.t$Names


# transpose again, remove tax data
JI_sewer_indic <- data.frame(t(JI_sewer_filt.df.t[,-c(1:7, ncol(JI_sewer_filt.df.t))]))
JI_human_indic <- data.frame(t(JI_human_filt.df.t[,-c(1:7, ncol(JI_human_filt.df.t))]))


# indicpsecies analysis
JI_sewer_3mo <- multipatt(JI_sewer_indic, as.vector(JI_info$Month), min.order = 3, max.order = 3, control = how(nperm = 999))
summary(JI_sewer_3mo, indvalcomp = TRUE)

JI_human_3mo <- multipatt(JI_human_indic, as.vector(JI_info$Month), min.order = 3, max.order = 3, control = how(nperm = 999))
summary(JI_human_3mo, indvalcomp = TRUE)
# output was copy & pasted into excel and organized




##############################
### environmental metadata ###
##############################


################
### Jones Island

# stats to test: only variables with enough data
variables <- colnames(JI_info)[c(2:5, 7, 10:11, 24:25)]


# CCA and environmental fit
identical(rownames(JI_relabun), rownames(JI_info))
JI.cca <- cca(JI_relabun ~ ., JI_info[, variables], na.action = na.exclude)
JI.fit <- envfit(JI.cca, JI_info[, variables], permu = 999, na.rm = TRUE)
JI.fit
#                       CCA1     CCA2     r2 Pr(>r)    
#   Year             0.84769 -0.53049 0.1183  0.038 *  
#   Air_temp_F       0.97133  0.23772 0.2439  0.001 ***
#   Precip_48hrs_in -0.93582 -0.35248 0.0036  0.911    
#   Flow_MGD        -0.80354 -0.59525 0.1577  0.031 *  
#   Ammonia_mgL      0.65490  0.75571 0.2720  0.002 ** 
#   BOD5_mgL         0.95243  0.30477 0.2798  0.001 ***
#   Phosphorus_mgL   0.97352 -0.22861 0.1672  0.020 *  
#   TSS_mgL          0.99846  0.05546 0.0536  0.120    
#
#           r2 Pr(>r)    
# Month 0.5633  0.001 ***


# is the metadata itself seasonal?
# omit NA's because zscore can't handle them
JI_metadata <- na.omit(as.matrix(JI_info[, c(4, 5, 7, 10:11, 24:25)]))


# normalize to z score
JI_metadata.z <- data.frame(t(zscore(t(JI_metadata))))


# ANOVA to check if environmental metadata is predicted by months
JI_aov.ls <- list()
for(i in colnames(JI_metadata.z)){
  JI_aov.ls[[i]] <- summary(aov(JI_metadata.z[[i]] ~ 
                                  subset(JI_info, Sample_name %in% rownames(JI_metadata.z))$Month))[[1]][["Pr(>F)"]]
}

JI_aov.df <- data.frame(p = unlist(JI_aov.ls))
JI_aov.df$variable <- "na"
JI_aov.df <- JI_aov.df[complete.cases(JI_aov.df),]
JI_aov.df$variable <- colnames(JI_metadata.z)
JI_aov.df$variable[JI_aov.df$p <= 0.05]
# [1] "Air_temp_F"  "Flow_MGD"    "Ammonia_mgL"



###############
### South Shore

identical(rownames(SS_relabun), rownames(SS_info))
SS.cca <- cca(SS_relabun ~ ., SS_info[, variables], na.action = na.exclude)
SS.fit <- envfit(SS.cca, SS_info[, variables], permu = 999, na.rm = TRUE)
SS.fit
#                       CCA1     CCA2     r2 Pr(>r)    
#   Year            -0.92048 -0.39079 0.3499  0.011 *  
#   Air_temp_F       0.34349 -0.93916 0.4875  0.003 ** 
#   Precip_48hrs_in -0.73783 -0.67499 0.3634  0.020 *  
#   Flow_MGD        -0.99183 -0.12755 0.8533  0.001 ***
#   Ammonia_mgL      0.97995  0.19926 0.5816  0.002 ** 
#   BOD5_mgL         0.96292 -0.26979 0.3809  0.011 *  
#   Phosphorus_mgL   0.99902  0.04420 0.4521  0.004 ** 
#   TSS_mgL          0.72119 -0.69274 0.5103  0.001 ***
#
#           r2 Pr(>r)  
# Month 0.6343  0.043 *


# omit NA's because zscore can't handle them
SS_metadata <- na.omit(as.matrix(SS_info[, c(4, 5, 7, 10:11, 24:25)]))


# normalize to z score
SS_metadata.z <- data.frame(t(zscore(t(SS_metadata))))


# ANOVA to check if environmental metadata is predicted by months
SS_aov.ls <- list()
for(i in colnames(SS_metadata.z)){
  SS_aov.ls[[i]] <- summary(aov(SS_metadata.z[[i]] ~ 
                                  subset(SS_info, Sample_name %in% rownames(SS_metadata.z))$Month))[[1]][["Pr(>F)"]]
}

SS_aov.df <- data.frame(p = unlist(SS_aov.ls))
SS_aov.df$variable <- "na"
SS_aov.df <- SS_aov.df[complete.cases(SS_aov.df),]
SS_aov.df$variable <- colnames(SS_metadata.z)
SS_aov.df$variable[SS_aov.df$p <= 0.05]
# [1] "Air_temp_F"

