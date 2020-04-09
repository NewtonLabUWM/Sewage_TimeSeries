##########################################
### Results section "Sewer bacteria drive 
### temporal trends"
### Lou LaMartina, finalized Mar 19, 2020
##########################################


setwd("~/Desktop/TimeSeries_final")
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(OTUtable)
library(reshape2)
library(phyloseq)


# load data (refer to TimeSeries1_DataPrep.R)
TimeSeries_object <- readRDS("./RData/TimeSeries_phyloseq_object.RData")
TimeSeries_info <- read.csv("./RData/TimeSeries_sample_info.csv")
rownames(TimeSeries_info) <- TimeSeries_info$Sample_name
TimeSeries_info$Collection_date <- as.Date(TimeSeries_info$Collection_date, format = "%m/%d/%y")
Neighborhood_info <- read.csv("./RData/Neighborhood_sample_info.csv")
Final_human_ASVs <- readRDS("./RData/Final_human_ASVs.RData")
Taxonomy_all <- read.csv("./RData/Taxonomy_all.csv")


# load ddPCR data
TimeSeries_ddPCR_data <- read.csv("./RData/TimeSeries_Flavo_Cloaci_HF183_ddPCR.csv")
Neighborhood_ddPCR_data <- read.csv("./RData/Neighborhood_Flavo_Cloaci_HF183_ddPCR.csv")


# human ASVs (refer to TimeSeries2_Threshold.R)
Final_human_ASVs <- readRDS("./RData/Final_human_ASVs.RData")


# indicspecies results (refer to TimeSeries5_Seasonal.R)
indic_ASVs <- read.csv("./RData/Indicator_ASVs.csv")


# subset the two treatment plants
JI_object <- subset_samples(TimeSeries_object, Treatment_plant == "Jones_Island")
JI_info <- subset(TimeSeries_info, Treatment_plant == "Jones_Island")


# convert to relative abundance
TimeSeries_relabun_object <- transform_sample_counts(TimeSeries_object, function(x) x / sum(x))
JI_relabun_object <- transform_sample_counts(JI_object, function(x) x / sum(x))


# extract data
TimeSeries_relabun <- data.frame(TimeSeries_relabun_object@otu_table@.Data)
JI_counts <- data.frame(JI_object@otu_table@.Data)
JI_relabun <- data.frame(JI_relabun_object@otu_table@.Data)


# split up human- and sewer-associated ASVs, remove empty ASVs, convert to relative abundances
JI_human_counts <- JI_counts[, colnames(JI_counts) %in% Final_human_ASVs]
JI_human_counts <- JI_human_counts[, colSums(JI_human_counts != 0) > 0]
JI_human_relabun <- JI_human_counts / rowSums(JI_human_counts)

JI_sewer_counts <- JI_counts[, ! colnames(JI_counts) %in% Final_human_ASVs]
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


# Function to calculate the mean and the standard deviation
# for each group
# http://www.sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


# decimal points in ggplot
# https://stackoverflow.com/questions/38722202/how-do-i-change-the-number-of-decimal-places-on-axis-labels-in-ggplot2/38722547
scaleFUN <- function(x) sprintf("%.1f", x)



#############################
### sewer dendro heat map ###
#############################

# convert to z score  = (x-μ)/σ
JI_sewer.z <- t(zscore(t(JI_sewer_filt)))


# change to chronological order
JI_info <- JI_info[order(JI_info$Month_point), ]
JI_sewer.z <- JI_sewer.z[match(rownames(JI_info), rownames(JI_sewer.z)),]
identical(rownames(JI_sewer.z), rownames(JI_info))


# add all tax
temp_sewer <- data.frame(t(JI_sewer.z))
temp_sewer$ASV <- rownames(temp_sewer)
temp_sewer <- merge(Taxonomy_all, temp_sewer, by = "ASV")
rownames(temp_sewer) <- paste0(temp_sewer$Phylum, "__", temp_sewer$Genus, "__", temp_sewer$ASV)
genera <- as.character(temp_sewer$Genus)
temp_sewer <- as.matrix(temp_sewer[,9:ncol(temp_sewer)])
JI_sewer.z <- t(temp_sewer)


# change names
rownames(JI_sewer.z) <- paste(JI_info$Month, JI_info$Year)


# euclidian distances on z score matrices
JI_ASV_dist <- vegdist(t(JI_sewer.z), method = "euclidian")


# cluster ASVs based on euclidian distances
JI_ASV_clus <- hclust(JI_ASV_dist, method = "average")


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


# heatmap with seasonal ASV colors
heatmap(as.matrix(t(JI_sewer.z)),
        Rowv = as.dendrogram(JI_ASV_clus),
        Colv = NA, margins = c(2, 12),
        col = heat_colors)


# make same heatmap in ggplot2, so i can actually modify it
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


### figure 4 sewer ###
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

#ggsave("./Plots/sewer heatmap.pdf", plot = sewer_dendro, device = "pdf", width = 3.5, height = 3.5, units = "in")




#############################
### human dendro heat map ###
#############################

# convert to z score  = (x-μ)/σ
JI_human.z <- t(zscore(t(JI_human_filt)))


# change to chronological order
JI_info <- JI_info[order(JI_info$Month_point), ]
JI_human.z <- JI_human.z[match(rownames(JI_info), rownames(JI_human.z)),]
identical(rownames(JI_human.z), rownames(JI_info))


# add all tax
temp_human <- data.frame(t(JI_human.z))
temp_human$ASV <- rownames(temp_human)
temp_human <- merge(Taxonomy_all, temp_human, by = "ASV")
rownames(temp_human) <- paste0(temp_human$Phylum, "__", temp_human$Genus, "__", temp_human$ASV)
genera <- as.character(temp_human$Genus)
temp_human <- as.matrix(temp_human[,9:ncol(temp_human)])
JI_human.z <- t(temp_human)


# change names
rownames(JI_human.z) <- paste(JI_info$Month, JI_info$Year)


# euclidian distances on z score matrices
JI_ASV_dist <- vegdist(t(JI_human.z), method = "euclidian")


# cluster ASVs based on euclidian distances
JI_ASV_clus <- hclust(JI_ASV_dist, method = "average")


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


# heatmap with phylum colors
heatmap(as.matrix(t(JI_human.z)),
        Rowv = as.dendrogram(JI_ASV_clus),
        Colv = NA, margins = c(2, 12),
        col = heat_colors)


# make same heatmap in ggplot2, so i can actually modify it
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


### figure 4 human ###
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

#ggsave("./Plots/human heatmap.pdf", plot = human_dendro, device = "pdf", width = 3.6, height = 3.5, units = "in")




########################
### autocorrelations ###
########################

# add relative abundances to ddPCR data.
# ASV44 matched HF183 sequence determined with local BLAST.
Sequence_data <- JI_relabun[c("ASV8", "ASV42", "ASV44")]
Sequence_data$Sample_name <- rownames(Sequence_data)


# add month point to abundance data
Sequence_data <- merge(Sequence_data, TimeSeries_info[c(1,3,8)], by = "Sample_name")


# get significance level - this is the default line that appears
# when you do the acf function with plot = TRUE
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


### figure 5A ###
autocorr <-
  ggplot(ACFs, aes(x = lag, y = acf)) +
  geom_hline(yintercept = significance_level, linetype = "dashed", color = "grey70", size = 0.25) +
  geom_hline(yintercept = -significance_level, linetype = "dashed", color = "grey70", size = 0.25) +
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
        panel.spacing = unit(1.2, "lines"),
        panel.border = element_rect(color = "grey80", fill = NA, size = 0.5), 
        strip.text = element_text(size = 6, color = "black", face = "bold.italic"),
        strip.background = element_rect(colour = "grey80", fill = "grey80"),
        axis.ticks = element_line(size = 0.25),
        axis.line = element_line(size = 0.25)) +
  geom_hline(aes(yintercept = 0), size = 0.25) +
  labs(y = "Autocorrelation", x = "Time lag (months)")
autocorr

#ggsave("./Plots/autocorr.pdf", plot = autocorr, device = "pdf", width = 5.2, height = 1.5, units = "in")




#########################
### time series ddPCR ###
#########################

# change order (to match facets from acf plot), change column names
TimeSeries_ddPCR_data <- TimeSeries_ddPCR_data[c(1,3,4,2,5:7)]
colnames(TimeSeries_ddPCR_data)[2:4]
# [1] "HF183_permL"    "cloaci08_permL" "flavo42_permL" 
colnames(TimeSeries_ddPCR_data)[2:4] <- c("Bacteroides", "Cloacibacterium", "Flavobacterium")


# melt
TS_ddPCR.m <- melt(TimeSeries_ddPCR_data, variable.name = "Genus", value.name = "Concentration", 
                id.vars = c("Collection_date", "Month", "Year", "Sample_name"))


# divide concentrations by 1 million, for simpler y axis in plot
TS_ddPCR.m$Concentration <- TS_ddPCR.m$Concentration / 1000000


# summarize
TS_summary <- data_summary(TS_ddPCR.m, groupnames = c("Genus", "Month"), varname = "Concentration")



### figure 5B ###
TSddPCR <-
  ggplot(TS_summary, aes(x = Month, y = Concentration, group = 1)) +
  geom_point(size = 0.8) +
  stat_summary(geom = "line", fun.y = "mean", size = 0.25) +
  geom_errorbar(aes(ymin = Concentration - sd, ymax = Concentration + sd), width = 0.2, size = 0.25) +
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
        panel.border = element_rect(color = "grey80", fill = NA, size = 0.5), 
        strip.text = element_text(size = 6, color = "black", face = "bold.italic"),
        strip.background = element_rect(colour = "white", fill = "white"),
        axis.ticks = element_line(size = 0.25),
        axis.line = element_line(size = 0.25)) +
  labs(y = bquote(bold('Gene copies per mL (x 10'^6*')')), x = "Month")
TSddPCR

#ggsave("./Plots/TSddPCR.pdf", plot = TSddPCR, device = "pdf", width = 5.2, height = 1.5, units = "in")


# how close are the sequence and absolute abundance data?
identical(as.character(TimeSeries_ddPCR_data$Sample_name), as.character(Sequence_data$Sample_name))

cor.test(Sequence_data$ASV8, TimeSeries_ddPCR_data$Cloacibacterium, method = "spearman") 
# rho = 0.8548167, p-value < 2.2e-16

cor.test(Sequence_data$ASV42, TimeSeries_ddPCR_data$Flavobacterium, method = "spearman") 
# rho = 0.8310086, p-value < 2.2e-16

cor.test(Sequence_data$ASV44, TimeSeries_ddPCR_data$Bacteroides, method = "spearman") 
# rho = 0.4941721, p-value = 0.0005657


# are abudances predicted by water temperatre?
JI_info <- JI_info[order(JI_info$Sample_name),]
identical(as.character(TimeSeries_ddPCR_data$Sample_name), as.character(JI_info$Sample_name))

cor.test(Sequence_data$ASV8, JI_info$Sewer_temp_F, method = "spearman") 
# rho = 0.8973468, p-value = 2.84e-09

cor.test(Sequence_data$ASV42, JI_info$Sewer_temp_F, method = "spearman") 
# rho = -0.8890823, p-value = p-value = 6.396e-09

cor.test(Sequence_data$ASV44, JI_info$Sewer_temp_F, method = "spearman") 
# rho = -0.2218356, p-value = 0.2975




##########################
### neighborhood ddPCR ###
##########################


# add time series quants
colnames(Neighborhood_ddPCR_data)[2:4]
colnames(Neighborhood_ddPCR_data)[2:4] <- c("Bacteroides", "Cloacibacterium", "Flavobacterium")
Neighborhood_ddPCR_data <- rbind(Neighborhood_ddPCR_data, TimeSeries_ddPCR_data[1:4])
Neighborhood_ddPCR_data$Source <- "Neighborhood"
Neighborhood_ddPCR_data$Source[Neighborhood_ddPCR_data$Sample_name %in% TimeSeries_info$Sample_name] <- "TimeSeries"


# melt
Neigh_ddPCR.m <- melt(Neighborhood_ddPCR_data, variable.name = "Genus", value.name = "Concentration", 
                   id.vars = c("Sample_name", "Source"))


# divide concentrations by 1 million, for simpler y axis in plot (same as time series)
Neigh_ddPCR.m$Concentration <- Neigh_ddPCR.m$Concentration / 1000000


# get means
means <- aggregate(. ~ Source + Genus, mean, data = Neigh_ddPCR.m[-1])
means
#         Source           Genus Concentration
# 1 Neighborhood     Bacteroides     1.6007460
# 2   TimeSeries     Bacteroides     0.4283243
# 3 Neighborhood Cloacibacterium     0.2190466
# 4   TimeSeries Cloacibacterium     1.5464150
# 5 Neighborhood  Flavobacterium     0.3713728
# 6   TimeSeries  Flavobacterium     0.4886650


### figure 5C ###
neigh_ddPCR <-
  ggplot(Neigh_ddPCR.m, aes(x = Source, y = Concentration)) +
  geom_jitter(size = 1, alpha = 0.8, shape = 1, color = "grey60") +
  facet_wrap(. ~ Genus, ncol = 3, scales = "free_y") +
  scale_x_discrete(labels = c("Residential\nsewers", "Treatment\nplant influent")) +
  scale_y_continuous(labels = scaleFUN) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black", face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "grey80", fill = NA, size = 0.5), 
        strip.text = element_text(size = 6, color = "black", face = "bold.italic"),
        strip.background = element_rect(colour = "white", fill = "white"),
        axis.ticks = element_line(size = 0.25),
        axis.line = element_line(size = 0.25)) +
  geom_crossbar(data = means, aes(ymin = Concentration, ymax = Concentration), size = 0.1, width = 0.8, color = "black") +
  labs(y = bquote(bold('Gene copies per mL (x 10'^6*')')))
neigh_ddPCR

#ggsave("./Plots/neigh_ddPCR.pdf", plot = neigh_ddPCR, device = "pdf", width = 5.2, height = 1.5, units = "in")




###################
### quant stats ###
###################

# time series ddPCR
TS_ddPCR_stats <- data.frame(aggregate(. ~ Genus, mean, data = TS_ddPCR.m[5:6]),
                        aggregate(. ~ Genus, sd, data = TS_ddPCR.m[5:6])[2])
colnames(TS_ddPCR_stats)[2:3] <- c("mean", "sd")
TS_ddPCR_stats
#             Genus      mean        sd
# 1     Bacteroides 0.4283243 0.1844833
# 2 Cloacibacterium 1.5464150 1.1118871
# 3  Flavobacterium 0.4886650 0.4615540


# increase Flavo neighborhoods to WWTP
means[6,3] / means[5,3]
# [1] 1.315834


# increase Cloaci neighborhoods to WWTP
means[4,3] / means[3,3]
# [1] 7.059754


# decrease Bacteroides neighborhoods to WWTP
means[1,3] / means[2,3]
# [1] 7.059754


# sequencing
TS_seq_stats <- data.frame(Genus = TS_ddPCR_stats$Genus, 
                           mean = c(mean(Sequence_data$ASV44), mean(Sequence_data$ASV8), mean(Sequence_data$ASV42)),
                           sd = c(sd(Sequence_data$ASV44), sd(Sequence_data$ASV8), sd(Sequence_data$ASV42)))
cbind(TS_seq_stats[1], TS_seq_stats[2:3] * 100)
#             Genus      mean         sd
# 1     Bacteroides 0.2103109 0.06555714
# 2 Cloacibacterium 1.1249916 0.68851254
# 3  Flavobacterium 0.4067064 0.24681467


