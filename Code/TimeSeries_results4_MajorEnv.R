##########################################
### Results section analyzing seasonality
### of human vs environmental communities.
### Lou LaMartina, finalized Jan 11, 2019
##########################################


setwd("~/Desktop/TimeSeries_final")
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(OTUtable)
library(reshape2)



#######################
### getting started ###
#######################

# prep counts
TimeSeries_counts <- read.csv("./RData/TimeSeries_counts.csv")
rownames(TimeSeries_counts) <- TimeSeries_counts$Sample_name
TimeSeries_counts <- TimeSeries_counts[-1]

# sample info
TimeSeries_info <- read.csv("./RData/TimeSeries_sample_info.csv")
rownames(TimeSeries_info) <- TimeSeries_info$Sample_name
TimeSeries_info$Collection_date <- as.Date(TimeSeries_info$Collection_date, format = "%m/%d/%y")

# tax
Taxonomy_all <- read.csv("./RData/Taxonomy_all.csv")
rownames(Taxonomy_all) <- Taxonomy_all$ASV

# load ddPCR data
ddPCR_data <- read.csv("./RData/Flavo_Cloaci_HF183_ddPCR.csv")
rownames(ddPCR_data) <- ddPCR_data$Sample_name

# human ASVs: how code on how this was made, go to Threshold_analysis.R
Final_human_ASVs <- readRDS("./RData/Final_human_ASVs.RData")

# contaminating sequences
Contaminants <- readRDS("./RData/Contaminants.RData")

# indicspecies results (see TimeSeries_indicspecies.R for code)
indic_ASVs <- read.csv("./RData/indic_ASVs.csv")

# function that allows you to extract p values from linear regression
# https://www.gettinggeneticsdone.com/2011/01/rstats-function-for-extracting-f-test-p.html
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1], f[2], f[3],lower.tail = F)
  attributes(p) <- NULL
  return(p)
}



# remove contaminating sequences
TimeSeries_counts <- TimeSeries_counts[, ! colnames(TimeSeries_counts) %in% Contaminants$ASV]

# remove empty ASVs and convert to relative abundance
TimeSeries_counts <- TimeSeries_counts[, colSums(TimeSeries_counts != 0) > 0]
TimeSeries_relabun <- TimeSeries_counts / rowSums(TimeSeries_counts)

# remove contaminants from taxonomy
Taxonomy_filt <- subset(Taxonomy_all, ! rownames(Taxonomy_all) %in% Contaminants$ASV)




###########################
###  subset & filter JI ###
###########################

JI_counts <- TimeSeries_counts[grepl("JI", rownames(TimeSeries_counts)), ]
JI_counts <- JI_counts[rowSums(JI_counts != 0) > 0, colSums(JI_counts != 0) > 0]
JI_relabun <- JI_counts / rowSums(JI_counts)

JI_info <- subset(TimeSeries_info, Treatment_plant == "Jones_Island")
rownames(JI_info) <- JI_info$Sample_name


# split up human and sewer
JI_human_counts <- JI_counts[, colnames(JI_counts) %in% Final_human_ASVs]
JI_human_relabun <- JI_human_counts / rowSums(JI_human_counts)

JI_sewer_counts <- JI_counts[, ! colnames(JI_counts) %in% Final_human_ASVs]
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
# [1] 60 63




###########################
###  subset & filter SS ###
###########################

SS_counts <- TimeSeries_counts[grepl("SS", rownames(TimeSeries_counts)), ]
SS_counts <- SS_counts[rowSums(SS_counts != 0) > 0, colSums(SS_counts != 0) > 0]
SS_relabun <- SS_counts / rowSums(SS_counts)

SS_info <- subset(TimeSeries_info, Treatment_plant == "South_Shore")
rownames(SS_info) <- SS_info$Sample_name


# split up human and sewer
SS_human_counts <- SS_counts[, colnames(SS_counts) %in% Final_human_ASVs]
SS_human_relabun <- SS_human_counts / rowSums(SS_human_counts)

SS_sewer_counts <- SS_counts[, ! colnames(SS_counts) %in% Final_human_ASVs]
SS_sewer_relabun <- SS_sewer_counts / rowSums(SS_sewer_counts)


# remove ASVs with less than 1% as their maximum relative abundance
maxtax <- apply(SS_sewer_relabun, 2, max)
mintax <- names(which(maxtax < 0.01))
SS_sewer_filt <- SS_sewer_relabun[, -which(colnames(SS_sewer_relabun) %in% mintax)]
dim(SS_sewer_filt)
# [1] 34 66

maxtax <- apply(SS_human_relabun, 2, max)
mintax <- names(which(maxtax < 0.01))
SS_human_filt <- SS_human_relabun[, -which(colnames(SS_human_relabun) %in% mintax)]
dim(SS_human_filt)
# [1] 34 57




####################################
### human vs sewer seasonal ASVs ###
####################################

################
### JONES ISLAND


### test whole dataset against months in MANOVA ###

JI_human.adonis <- adonis(JI_human_filt ~ JI_info$Month, method = "bray", perm = 999)
JI_human.adonis
# R2 0.41515
# p 0.001

JI_sewer.adonis <- adonis(JI_sewer_filt ~ JI_info$Month, method = "bray", perm = 999)
JI_sewer.adonis
# R2 0.55037
# p 0.001


### test relative abundances of individual ASVs against months ###

### human
temp <- cbind(JI_human_filt, Month = JI_info$Month)

JI_human_R2 <- list()
JI_human_p <- list()
for(i in names(temp)[-ncol(temp)]){
  JI_human_R2[[i]] <- summary(lm(paste(i , "~ Month"), temp))$adj.r.squared
  JI_human_p[[i]] <- lmp(lm(paste(i , "~ Month"), temp))
}

JI_human_R2.df <- data.frame(R2 = unlist(JI_human_R2))
JI_human_p.df <- data.frame(p = unlist(JI_human_p))
JI_human_lm.df <- data.frame(JI_human_R2.df, JI_human_p.df)

# proportion whose p is < 0.05?
length(which(JI_human_lm.df$p < 0.05)) / nrow(JI_human_R2.df) * 100
# [1] 38.46154 % significantly correlated to month

# add tax info
JI_human_lm.df$ASV <- rownames(JI_human_lm.df)
JI_human_lm.df <- merge(JI_human_lm.df, Taxonomy_all[c(1, 6:7)], by = "ASV")
JI_human_lm.df <- subset(JI_human_lm.df, p < 0.05)


### sewer
temp <- cbind(JI_sewer_filt, Month = JI_info$Month)

JI_sewer_R2 <- list()
JI_sewer_p <- list()
for(i in names(temp)[-ncol(temp)]){
  JI_sewer_R2[[i]] <- summary(lm(paste(i , "~ Month"), temp))$adj.r.squared
  JI_sewer_p[[i]] <- lmp(lm(paste(i , "~ Month"), temp))
}

JI_sewer_R2.df <- data.frame(R2 = unlist(JI_sewer_R2))
JI_sewer_p.df <- data.frame(p = unlist(JI_sewer_p))
JI_sewer_lm.df <- data.frame(JI_sewer_R2.df, JI_sewer_p.df)

# proportion whose p is < 0.05?
length(which(JI_sewer_lm.df$p < 0.05)) / nrow(JI_sewer_R2.df) * 100
# [1] 65.21739 % significantly correlated to month

# add tax info
JI_sewer_lm.df$ASV <- rownames(JI_sewer_lm.df)
JI_sewer_lm.df <- merge(JI_sewer_lm.df, Taxonomy_all[c(1, 6:7)], by = "ASV")
JI_sewer_lm.df <- subset(JI_sewer_lm.df, p < 0.05)





###############
### SOUTH SHORE 

### test whole dataset against months in MANOVA
SS_human.adonis <- adonis(SS_human_filt ~ SS_info$Month, method = "bray", perm = 999)
SS_human.adonis
# p 0.084

SS_sewer.adonis <- adonis(SS_sewer_filt ~ SS_info$Month, method = "bray", perm = 999)
SS_sewer.adonis
# p 0.005


### test individual ASVs ###

### human
temp <- cbind(SS_human_filt, Month = SS_info$Month)

SS_human_R2 <- list()
SS_human_p <- list()
for(i in names(temp)[-ncol(temp)]){
  SS_human_R2[[i]] <- summary(lm(paste(i , "~ Month"), temp))$adj.r.squared
  SS_human_p[[i]] <- lmp(lm(paste(i , "~ Month"), temp))
}

SS_human_R2.df <- data.frame(R2 = unlist(SS_human_R2))
SS_human_p.df <- data.frame(p = unlist(SS_human_p))
SS_human_lm.df <- data.frame(SS_human_R2.df, SS_human_p.df)

# proportion whose p is < 0.05?
length(which(SS_human_lm.df$p < 0.05)) / nrow(SS_human_R2.df) * 100
# [1] 10.52632 % significantly correlated to month

# add tax info
SS_human_lm.df$ASV <- rownames(SS_human_lm.df)
SS_human_lm.df <- merge(SS_human_lm.df, Taxonomy_all[c(1, 6:7)], by = "ASV")
SS_human_lm.df <- subset(SS_human_lm.df, p < 0.05)



### sewer
temp <- cbind(SS_sewer_filt, Month = SS_info$Month)

SS_sewer_R2 <- list()
SS_sewer_p <- list()
for(i in names(temp)[-ncol(temp)]){
  SS_sewer_R2[[i]] <- summary(lm(paste(i , "~ Month"), temp))$adj.r.squared
  SS_sewer_p[[i]] <- lmp(lm(paste(i , "~ Month"), temp))
}

SS_sewer_R2.df <- data.frame(R2 = unlist(SS_sewer_R2))
SS_sewer_p.df <- data.frame(p = unlist(SS_sewer_p))
SS_sewer_lm.df <- data.frame(SS_sewer_R2.df, SS_sewer_p.df)

# proportion whose p is < 0.05?
length(which(SS_sewer_lm.df$p < 0.05)) / nrow(SS_sewer_R2.df) * 100
# [1] 36.36364 % significantly correlated to month

# add tax info
SS_sewer_lm.df$ASV <- rownames(SS_sewer_lm.df)
SS_sewer_lm.df <- merge(SS_sewer_lm.df, Taxonomy_all[c(1, 6:7)], by = "ASV")
SS_sewer_lm.df <- subset(SS_sewer_lm.df, p < 0.05)


### plot
JI_sewer_lm.df <- droplevels(JI_sewer_lm.df) # drop levels leftover from subsetting
JI_sewer_lm_genera <- data.frame(table(JI_sewer_lm.df$Genus))
colnames(JI_sewer_lm_genera)[1] <- "Genus"
JI_sewer_lm_genera$WWTP <- "Jones_Island"
JI_sewer_lm_genera$Source <- "Sewer"
JI_sewer_lm_genera <- unique(merge(JI_sewer_lm_genera, JI_sewer_lm.df[4:5], by = "Genus"))

JI_human_lm.df <- droplevels(JI_human_lm.df) 
JI_human_lm_genera <- data.frame(table(JI_human_lm.df$Genus))
colnames(JI_human_lm_genera)[1] <- "Genus"
JI_human_lm_genera$WWTP <- "Jones_Island"
JI_human_lm_genera$Source <- "Human"
JI_human_lm_genera <- unique(merge(JI_human_lm_genera, JI_human_lm.df[4:5], by = "Genus"))

SS_sewer_lm.df <- droplevels(SS_sewer_lm.df) 
SS_sewer_lm_genera <- data.frame(table(SS_sewer_lm.df$Genus))
colnames(SS_sewer_lm_genera)[1] <- "Genus"
SS_sewer_lm_genera$WWTP <- "South_Shore"
SS_sewer_lm_genera$Source <- "Sewer"
SS_sewer_lm_genera <- unique(merge(SS_sewer_lm_genera, SS_sewer_lm.df[4:5], by = "Genus"))

SS_human_lm.df <- droplevels(SS_human_lm.df)
SS_human_lm_genera <- data.frame(table(SS_human_lm.df$Genus))
colnames(SS_human_lm_genera)[1] <- "Genus"
SS_human_lm_genera$WWTP <- "South_Shore"
SS_human_lm_genera$Source <- "Human"
SS_human_lm_genera <- unique(merge(SS_human_lm_genera, SS_human_lm.df[4:5], by = "Genus"))

# combine
All_lm_genera <- rbind(JI_sewer_lm_genera, JI_human_lm_genera, SS_sewer_lm_genera, SS_human_lm_genera)
length(unique(All_lm_genera$Genus))
length(unique(All_lm_genera$Family))

# most common seasonal sewer family?
temp <- subset(All_lm_genera, Source == "Sewer")
temp <- ddply(temp, .(Family), summarize, Freq.sum = sum(Freq))
temp[order(temp$Freq.sum, decreasing = TRUE),]$Family[1:3]
# [1] Acinetobacter  Arcobacter     Flavobacterium

# most common seasonal human family?
temp <- subset(All_lm_genera, Source == "Human")
temp <- ddply(temp, .(Family), summarize, Freq.sum = sum(Freq))
temp[order(temp$Freq.sum, decreasing = TRUE),]$Family[1:3]
# [1] Bacteroides   Klebsiella    Acinetobacter

# plot
bar_colors <- colorRampPalette(brewer.pal(11, "Spectral"))

human_bars <- 
  ggplot(subset(All_lm_genera, Source == "Human"), aes(x = WWTP, y = Freq, fill = Family, color = Family)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = bar_colors(length(unique(subset(All_lm_genera, Source == "Human")$Family)))) +
  scale_color_manual(values = bar_colors(length(unique(subset(All_lm_genera, Source == "Human")$Family)))) +
  theme_classic() +
  ylim(0, 45) +
  scale_x_discrete(labels = c("Jones\nIsland", "South\nShore")) +
  theme(axis.text.x = element_text(size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black", face = "bold"),
        axis.title.x = element_text(size = 6, color = "black", face = "bold"),
        legend.text = element_text(size = 6, color = "black", face = "italic"),
        legend.title = element_text(size = 6, color = "black", face = "bold")) +
  guides(fill = guide_legend(keyheight = 0.4, keywidth = 0.4, units = "in", ncol = 1)) +
  labs(y = "Families of seasonal human-associated ASVs", x = "Wastewater treatment plant")
human_bars

#ggsave("./Plots/human_bars.pdf", plot = human_bars, device = "pdf", width = 2.5, height = 3, units = "in")


sewer_bars <- 
  ggplot(subset(All_lm_genera, Source == "Sewer"), aes(x = WWTP, y = Freq, fill = Family, color = Family)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = bar_colors(length(unique(subset(All_lm_genera, Source == "Sewer")$Family)))) +
  scale_color_manual(values = bar_colors(length(unique(subset(All_lm_genera, Source == "Sewer")$Family)))) +
  theme_classic() +
  scale_x_discrete(labels = c("Jones\nIsland", "South\nShore")) +
  theme(axis.text.x = element_text(size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black", face = "bold"),
        axis.title.x = element_text(size = 6, color = "black", face = "bold"),
        legend.text = element_text(size = 6, color = "black", face = "italic"),
        legend.title = element_text(size = 6, color = "black", face = "bold")) +
  guides(fill = guide_legend(keyheight = 0.4, keywidth = 0.4, units = "in", ncol = 1)) +
  labs(y = "Families of seasonal sewer-associated ASVs", x = "Wastewater treatment plant")
sewer_bars

#ggsave("./Plots/sewer_bars.pdf", plot = sewer_bars, device = "pdf", width = 2.4, height = 3, units = "in")




################################
### Flavo vs Cloaci vs HF183 ###
################################

# add relative abundances to ddPCR data.
# ASV44 matched HF183 sequence determined with local BLAST.
Sequence_data <- JI_relabun[c("ASV8", "ASV42", "ASV44")]
Sequence_data$Sample_name <- rownames(Sequence_data)

# add month point to abundance data
ddPCR_data <- merge(ddPCR_data, TimeSeries_info[c(1,6)], by = "Sample_name")
Sequence_data <- merge(Sequence_data, TimeSeries_info[c(1,2,6)], by = "Sample_name")


####################
### autocorrelations

significance_level <- qnorm((1 + 0.95)/2)/sqrt(sum(JI_relabun[-1]))

# stat
ASV8_acf <- acf(Sequence_data$ASV8, lag.max = 60, plot = FALSE)
ASV8_acf.df <- with(ASV8_acf, data.frame(lag, acf))
ASV8_acf.df$Genus <- "Cloacibacterium"

ASV42_acf <- acf(Sequence_data$ASV42, lag.max = 60, plot = FALSE)
ASV42_acf.df <- with(ASV42_acf, data.frame(lag, acf))
ASV42_acf.df$Genus <- "Flavobacterium"

ASV44_acf <- acf(Sequence_data$ASV44, lag.max = 60, plot = FALSE)
ASV44_acf.df <- with(ASV44_acf, data.frame(lag, acf))
ASV44_acf.df$Genus <- "Bacteroides"

# combine (remove lag 0)
ACFs <- rbind(ASV8_acf.df[-1,], ASV42_acf.df[-1,], ASV44_acf.df[-1,])

# plot
autocorr <-
  ggplot(ACFs, aes(x = lag, y = acf)) +
  geom_hline(yintercept = significance_level, linetype = "dashed", color = "grey90", size = 0.25) +
  geom_hline(yintercept = -significance_level, linetype = "dashed", color = "grey90", size = 0.25) +
  geom_segment(aes(xend = lag, yend = 0), size = 0.6, color = "black") +
  facet_wrap(. ~ Genus, ncol = 3) +
  theme_classic() +
  scale_y_continuous(breaks = c(-0.2, 0, 0.2, 0.4, 0.6)) +
  scale_x_continuous(breaks = c(1, 12, 24, 36, 48, 60), limits = c(1, 60)) +
  theme(axis.text.x = element_text(size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black", face = "bold"),
        axis.title.x = element_text(size = 6, color = "black", face = "bold"),
        legend.position = "none",
        panel.spacing = unit(1.3, "lines"),
        panel.border = element_rect(color = "grey80", fill = NA, size = 0.5), 
        strip.text = element_text(size = 6, color = "black", face = "bold.italic"),
        strip.background = element_rect(colour = "white", fill = "white"),
        axis.line = element_line(size = 0.25),
        axis.ticks = element_line(size = 0.25)) +
  geom_hline(aes(yintercept = 0), size = 0.25) +
  labs(y = "Autocorrelation", x = "Time lag (months)")
autocorr

ggsave("./Plots/autocorr.pdf", plot = autocorr, device = "pdf", width = 5.2, height = 1.5, units = "in")



#########
### ddPCR

# linear regression of month and quantitative data
temp <- ddPCR_data[c(2:4,6)]

ddPCR_p <- list()
ddPCR_R2 <- list()
for(i in names(temp)[-ncol(temp)]){
  ddPCR_p[[i]] <- lmp(lm(paste(i , "~ Month"), temp))
  ddPCR_R2[[i]] <- summary(lm(paste(i , "~ Month"), temp))$adj.r.squared
}
ddPCR_p
# $flavo42_permL
# [1] 0.02305633
# 
# $HF183_permL
# [1] 0.02630589
# 
# $cloaci08_permL
# [1] 4.884896e-10

ddPCR_R2
# $flavo42_permL
# [1] 0.1959303
# 
# $HF183_permL
# [1] 0.2516699
# 
# $cloaci08_permL
# [1] 0.6605403

# can't have NA when getting z scores
ddPCR.na <- na.omit(ddPCR_data)

# calculate z score
# ddPCR.z <- data.frame(t(zscore(t(ddPCR.na[2:4]))))

# add z scores to data frame
# ddPCR.z <- data.frame(ddPCR.z[2], ddPCR.z[3], ddPCR.z[1]) # change order, for facet later
# ddPCR.na <- data.frame(ddPCR.na[-c(2:4)], ddPCR.z)
# colnames(ddPCR.na)[6:8] <- c("Bacteroides", "Cloacibacterium", "Flavobacterium")

# divide concentrations by 1 million
ddPCR.na[2:4] <- ddPCR.na[2:4] / 1000000

# melt
colnames(ddPCR.na)
ddPCR.na <- ddPCR.na[c(1, 3, 4, 2, 5:8)]
colnames(ddPCR.na)[2:4] <- c("Bacteroides", "Cloacibacterium", "Flavobacterium") # change order, for facet later
ddPCR.m <- melt(ddPCR.na, variable.name = "Genus", value.name = "Concentration", 
                id.vars = c("Collection_date", "Month", "Month_point", "Year", "Sample_name"))

# boxplot
scaleFUN <- function(x) sprintf("%.1f", x)
ddPCR <-
  ggplot(ddPCR.m, aes(x = Month, y = Concentration)) +
  geom_boxplot(width = 0.4, outlier.size = 0.5, size = 0.3, color = "black", fill = "grey90") +
  facet_wrap(. ~ Genus, ncol = 3, scales = "free") +
  theme_classic() +
  scale_x_discrete(limits = c("January", "February", "March", "April",
                              "May", "June", "July", "August",
                              "September", "October", "November", "December"),
                   labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) +
  scale_y_continuous(labels = scaleFUN) +
    theme(axis.text.x = element_text(size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black", face = "bold"),
        axis.title.x = element_text(size = 6, color = "black", face = "bold"),
        legend.position = "none",
        panel.spacing = unit(0.5, "lines"),
        panel.border = element_rect(color = "grey80", fill = NA, size = 0.25), 
        strip.text = element_text(size = 6, color = "black", face = "bold.italic"),
        strip.background = element_rect(colour = "white", fill = "white"),
        axis.line = element_line(size = 0.25),
        axis.ticks = element_line(size = 0.25)) +
  labs(y = "Counts per mL sewage\nx 1 million", x = "Month")
ddPCR

ggsave("./Plots/ddPCR.pdf", plot = ddPCR, device = "pdf", width = 5.3, height = 1.5, units = "in")


### means
# HF183
mean(Sequence_data$ASV44) # 0.002102685
sd(Sequence_data$ASV44)   # 0.0006555036

mean(ddPCR_data$HF183_permL, na.rm = TRUE) # 428324
sd(ddPCR_data$HF183_permL, na.rm = TRUE)   # 184483

# Cloaci
mean(Sequence_data$ASV8) # 0.01124779
sd(Sequence_data$ASV8)   # 0.006884401

mean(ddPCR_data$cloaci08_permL, na.rm = TRUE) # 1546415
sd(ddPCR_data$cloaci08_permL, na.rm = TRUE)   # 1111887

# Flavo
mean(Sequence_data$ASV42) # 0.004066128
sd(Sequence_data$ASV42)   # 0.002467459

mean(ddPCR_data$flavo42_permL, na.rm = TRUE) # 488665
sd(ddPCR_data$flavo42_permL, na.rm = TRUE)   # 461554



##############
### Sequencing

### linear regression values from above

# Cloaci
subset(JI_sewer_lm.df, ASV == "ASV8")[,2:3]
#        R2           p
# 0.7103748 1.36339e-11

# Flavo
subset(JI_sewer_lm.df, ASV == "ASV11")[,2:3]
#        R2           p
# 0.4021164 9.223368e-05

# Bacteroides
JI_human_lm.df2 <- data.frame(JI_human_R2.df, JI_human_p.df)
JI_human_lm.df2$ASV <- rownames(JI_human_lm.df2)
JI_human_lm.df2 <- merge(JI_human_lm.df2, Taxonomy_all[c(1, 6:7)], by = "ASV")
subset(JI_human_lm.df2, ASV == "ASV44")[,2:3]
#        R2           p
# 0.1547701 0.05155181




#############################
### sewer dendro heat map ###
#############################

# convert to z score  = (x-μ)/σ
JI_sewer.z <- t(zscore(t(JI_sewer_filt)))
JI_human.z <- t(zscore(t(JI_human_filt)))

# change to chronological order
JI_info <- JI_info[order(JI_info$Month_point), ]
JI_sewer.z <- JI_sewer.z[match(rownames(JI_info), rownames(JI_sewer.z)),]
identical(rownames(JI_sewer.z), rownames(JI_info))

JI_human.z <- JI_human.z[match(rownames(JI_info), rownames(JI_human.z)),]
identical(rownames(JI_human.z), rownames(JI_info))

# add all tax
temp <- data.frame(t(JI_sewer.z))
temp$ASV <- rownames(temp)
temp <- merge(Taxonomy_all, temp, by = "ASV")
rownames(temp) <- paste0(temp$Phylum, "__", temp$Genus, "__", temp$ASV)
genera <- as.character(temp$Genus)
temp <- as.matrix(temp[,9:ncol(temp)])
JI_sewer.z <- t(temp)

# change names
rownames(JI_sewer.z) <- paste(JI_info$Month, JI_info$Year)


### stats

# euclidian distances on z score matrices
JI_ASV_dist <- vegdist(t(JI_sewer.z), method = "euclidian")

# cluster ASVs based on euclidian distances
JI_ASV_clus <- hclust(JI_ASV_dist, "aver")

# add taxa info
ASVs.df <- data.frame(dendro_order = 1:length(colnames(JI_sewer.z)), ASV = colnames(JI_sewer.z))
ASVs.df$Phylum <- sapply(strsplit(as.character(ASVs.df$ASV), "__"), `[`, 1)
ASVs.df$Genus <- sapply(strsplit(as.character(ASVs.df$ASV), "__"), `[`, 2)

# define colors: blue is low, red is high
heat_colors <- rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(299))

# extract order of dendro tree leaves
JI_ASV_clus_order <- JI_ASV_clus$labels[c(JI_ASV_clus$order)]
JI_ASV_clus_order.df <- data.frame(t(data.frame(strsplit(JI_ASV_clus_order, "__"))))
colnames(JI_ASV_clus_order.df) <- c("Phylum", "Genus", "ASV_ID")
JI_ASV_clus_order.df$ASV <- JI_ASV_clus_order
rownames(JI_ASV_clus_order.df) <- JI_ASV_clus_order.df$ASV_ID

# indicator ASV colors
indic_ASVs$Season <- "Spring"
indic_ASVs$Season[indic_ASVs$Window > 6] <- "Fall"
ASVs.df$ASV <- sapply(strsplit(as.character(ASVs.df$ASV), "__", 1), '[', 3)
Seasonal_colors <- merge(indic_ASVs, ASVs.df, by = "ASV", all.y = TRUE)
Seasonal_colors$Seasonal[is.na(Seasonal_colors$Season) == TRUE] <- "FALSE"
Seasonal_colors <- Seasonal_colors[order(Seasonal_colors$dendro_order), ]
Seasonal_colors$Seas_colors[Seasonal_colors$Season == "Spring"] <- "#66C2A5"
Seasonal_colors$Seas_colors[Seasonal_colors$Season == "Fall"] <- "#F46D43"

# heatmap with seasonal ASV colors
heatmap(as.matrix(t(JI_sewer.z)),
        Rowv = as.dendrogram(JI_ASV_clus),
        Colv = NA, margins = c(2, 12),
        col = heat_colors, 
        RowSideColors = as.character(Seasonal_colors$Seas_colors))


### make same heatmap in ggplot2, so i can actually modify it
JI_sewer.z.df <- data.frame(JI_sewer.z)
JI_sewer.z.df$Sample_name <- rownames(JI_sewer.z.df)
JI_sewer.z.df <- JI_sewer.z.df[match(paste(JI_info$Month, JI_info$Year), JI_sewer.z.df$Sample_name), ]
JI_sewer.z.df$month_order <- 1:nrow(JI_sewer.z.df)
JI_sewer.z.df.m <- melt(JI_sewer.z.df, id.vars = c("month_order", "Sample_name"), variable.name = "ASV", value.name = "Z")
JI_sewer.z.df.m$ASV_ID <- sapply(strsplit(as.character(JI_sewer.z.df.m$ASV), "__"), '[', 3)

# need to keep order of ASVs in dendrogram
JI_ASV_clus_order.df$dendro_order <- 1:nrow(JI_ASV_clus_order.df)
JI_sewer.z.df.m <- merge(JI_sewer.z.df.m, JI_ASV_clus_order.df[c(3,5)], by = "ASV_ID")
JI_sewer.z.df.m <- JI_sewer.z.df.m[order(JI_sewer.z.df.m$dendro_order), ]

# make ASVs an ordered factor
JI_sewer.z.df.m$ASV_ID <- factor(JI_sewer.z.df.m$ASV_ID, levels = unique(JI_sewer.z.df.m$ASV_ID))

# specify color breaks (roughly using quantiles)
heat_colors <- c(rev(brewer.pal(9, "Blues")), "white", brewer.pal(9, "YlOrRd")[1:4])
heat_values <- c(seq(0, 0.2, length.out = 9), seq(0.21, 0.4, length.out = 2), seq(0.41, 1, length.out = 2))

# custom y axis labels
sewer_labels <- as.character(JI_ASV_clus_order.df$Genus)
sewer_labels[sewer_labels == ""] <- "sp."
sewer_labels <- paste0(JI_ASV_clus_order.df$Phylum, " (", sewer_labels, ")")


# plot
sewer_dendro <- 
  ggplot(JI_sewer.z.df.m, aes(x = month_order, y = ASV_ID, fill = Z)) +
  geom_tile() +
  scale_fill_gradientn(colors = heat_colors, values = heat_values, breaks = c(-2, 0, 2, 4, 6)) +
  theme_classic() +
  scale_y_discrete(labels = sewer_labels, position = "right") +
  scale_x_continuous(breaks = c(1, 12, 24, 36, 48, 60)) +
  theme(axis.text.x = element_text(size = 6, color = "black"),
        axis.text.y = element_text(size = 2.5, color = "black", face = "italic"),
        axis.title.y = element_text(size = 6, color = "black", face = "bold"),
        axis.title.x = element_text(size = 6, color = "black", face = "bold"),
        legend.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 6, color = "black", face = "bold"),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_blank(),
        axis.ticks = element_line(size = 0.25),
        axis.ticks.y = element_blank(),
        panel.border = element_blank()) +
  guides(fill = guide_colorbar(barheight = 2, barwidth = 0.3, units = "in")) +
  labs(x = "Month", y = "Sewer-associated ASVs\nPhylum (genus)", fill = "Z score")
sewer_dendro

#ggsave("./Plots/sewer heatmap.pdf", plot = sewer_dendro, device = "pdf", width = 4.5, height = 3.5, units = "in")




#############################
### human dendro heat map ###
#############################

# add all tax
temp <- data.frame(t(JI_human.z))
temp$ASV <- rownames(temp)
temp <- merge(Taxonomy_all, temp, by = "ASV")
rownames(temp) <- paste0(temp$Phylum, "__", temp$Genus, "__", temp$ASV)
genera <- as.character(temp$Genus)
temp <- as.matrix(temp[,9:ncol(temp)])
JI_human.z <- t(temp)

# change names
rownames(JI_human.z) <- paste(JI_info$Month, JI_info$Year)


### stats

# euclidian distances on z score matrices
JI_ASV_dist <- vegdist(t(JI_human.z), method = "euclidian")

# cluster ASVs based on euclidian distances
JI_ASV_clus <- hclust(JI_ASV_dist, "aver")

# add taxa info
ASVs.df <- data.frame(dendro_order = 1:length(colnames(JI_human.z)), ASV = colnames(JI_human.z))
ASVs.df$Phylum <- sapply(strsplit(as.character(ASVs.df$ASV), "__"), `[`, 1)
ASVs.df$Genus <- sapply(strsplit(as.character(ASVs.df$ASV), "__"), `[`, 2)

# define colors: blue is low, red is high
heat_colors <- rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(299))

# extract order of dendro tree leaves
JI_ASV_clus_order <- JI_ASV_clus$labels[c(JI_ASV_clus$order)]
JI_ASV_clus_order.df <- data.frame(t(data.frame(strsplit(JI_ASV_clus_order, "__"))))
colnames(JI_ASV_clus_order.df) <- c("Phylum", "Genus", "ASV_ID")
JI_ASV_clus_order.df$ASV <- JI_ASV_clus_order
rownames(JI_ASV_clus_order.df) <- JI_ASV_clus_order.df$ASV_ID

# indicator ASV colors
indic_ASVs$Season <- "Spring"
indic_ASVs$Season[indic_ASVs$Window > 6] <- "Fall"
ASVs.df$ASV <- sapply(strsplit(as.character(ASVs.df$ASV), "__", 1), '[', 3)
Seasonal_colors <- merge(indic_ASVs, ASVs.df, by = "ASV", all.y = TRUE)
Seasonal_colors$Seasonal[is.na(Seasonal_colors$Season) == TRUE] <- "FALSE"
Seasonal_colors <- Seasonal_colors[order(Seasonal_colors$dendro_order), ]
Seasonal_colors$Seas_colors[Seasonal_colors$Season == "Spring"] <- "#66C2A5"
Seasonal_colors$Seas_colors[Seasonal_colors$Season == "Fall"] <- "#F46D43"

# heatmap with phylum colors
heatmap(as.matrix(t(JI_human.z)),
        Rowv = as.dendrogram(JI_ASV_clus),
        Colv = NA, margins = c(2, 12),
        col = heat_colors, 
        RowSideColors = as.character(Seasonal_colors$Seas_colors))
# save as 10 in x 14 in



### make same heatmap in ggplot2, so i can actually modify it
JI_human.z.df <- data.frame(JI_human.z)
JI_human.z.df$Sample_name <- rownames(JI_human.z.df)
JI_human.z.df <- JI_human.z.df[match(paste(JI_info$Month, JI_info$Year), JI_human.z.df$Sample_name), ]
JI_human.z.df$month_order <- 1:nrow(JI_human.z.df)
JI_human.z.df.m <- melt(JI_human.z.df, id.vars = c("month_order", "Sample_name"), variable.name = "ASV", value.name = "Z")
JI_human.z.df.m$ASV_ID <- sapply(strsplit(as.character(JI_human.z.df.m$ASV), "__"), '[', 3)

# need to keep order of ASVs in dendrogram
JI_ASV_clus_order.df$dendro_order <- 1:nrow(JI_ASV_clus_order.df)
JI_human.z.df.m <- merge(JI_human.z.df.m, JI_ASV_clus_order.df[c(3,5)], by = "ASV_ID")
JI_human.z.df.m <- JI_human.z.df.m[order(JI_human.z.df.m$dendro_order), ]

# make ASVs an ordered factor
JI_human.z.df.m$ASV_ID <- factor(JI_human.z.df.m$ASV_ID, levels = unique(JI_human.z.df.m$ASV_ID))

# specify color breaks (roughly using quantiles)
heat_colors <- c(rev(brewer.pal(9, "Blues")), "white", brewer.pal(9, "YlOrRd")[1:4])
heat_values <- c(seq(0, 0.3, length.out = 9), seq(0.31, 0.4, length.out = 2), seq(0.41, 1, length.out = 2))

# custom y axis labels
human_labels <- as.character(JI_ASV_clus_order.df$Genus)
human_labels[human_labels == ""] <- "sp."
human_labels <- paste0(JI_ASV_clus_order.df$Phylum, " (", human_labels, ")")

# plot
human_dendro <- 
  ggplot(JI_human.z.df.m, aes(x = month_order, y = ASV_ID, fill = Z)) +
  geom_tile() +
  scale_fill_gradientn(colors = heat_colors, values = heat_values, breaks = c(-2, 0, 2, 4, 6)) +
  theme_classic() +
  scale_y_discrete(labels = human_labels, position = "right") +
  scale_x_continuous(breaks = c(1, 12, 24, 36, 48, 60)) +
  theme(axis.text.x = element_text(size = 6, color = "black"),
        axis.text.y = element_text(size = 2.5, color = "black", face = "italic"),
        axis.title.y = element_text(size = 6, color = "black", face = "bold"),
        axis.title.x = element_text(size = 6, color = "black", face = "bold"),
        legend.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 6, color = "black", face = "bold"),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_blank(),
        axis.ticks = element_line(size = 0.25),
        axis.ticks.y = element_blank(),
        panel.border = element_blank()) +
  guides(fill = guide_colorbar(barheight = 2, barwidth = 0.3, units = "in")) +
  labs(x = "Month", y = "Human-associated ASVs\nPhylum (genus)", fill = "Z score")
human_dendro

#ggsave("./Plots/human heatmap.pdf", plot = human_dendro, device = "pdf", width = 4.6, height = 3.5, units = "in")


### plots were overlaid and converted to PDF in Keynote ###

