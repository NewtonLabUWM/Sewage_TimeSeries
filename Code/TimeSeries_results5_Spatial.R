##########################################
### Results section comparing MKE temporal
### trends to US spatial trends.
### Lou LaMartina, finalized Dec 12, 2019
##########################################


setwd("~/Desktop/TimeSeries_final")
library(vegan)
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(usmap)
library(ggrepel)



#######################
### getting started ###
#######################

### prep counts
TimeSeries_counts <- read.csv("./RData/TimeSeries_counts.csv")
rownames(TimeSeries_counts) <- TimeSeries_counts$Sample_name
TimeSeries_counts <- TimeSeries_counts[-1]
TimeSeries_counts <- TimeSeries_counts[rowSums(TimeSeries_counts != 0) > 0, colSums(TimeSeries_counts != 0) > 0]

Cities_counts <- read.csv("./RData/Cities_counts.csv")
rownames(Cities_counts) <- Cities_counts$Sample_name
Cities_counts <- Cities_counts[-1]
Cities_counts <- Cities_counts[rowSums(Cities_counts != 0) > 0, colSums(Cities_counts != 0) > 0]
Cities_counts <- Cities_counts[- which(grepl("Reus", rownames(Cities_counts))), ] # remove Reus

### sample info
TimeSeries_info <- read.csv("./RData/TimeSeries_sample_info.csv")
rownames(TimeSeries_info) <- TimeSeries_info$Sample_name

Cities_info <- read.csv("./RData/Cities_sample_info.csv")
rownames(Cities_info) <- Cities_info$Sample_name
Cities_info <- Cities_info[Cities_info$City != "Reus", ] # remove Reus
Cities_info$Temp_profile_metric <- as.numeric(as.character(Cities_info$Temp_profile_metric))
Cities_info$Collection_date <- as.Date(Cities_info$Collection_date, format = "%m/%d/%y")

### tax
Taxonomy_all <- read.csv("./RData/Taxonomy_all.csv")
rownames(Taxonomy_all) <- Taxonomy_all$ASV
Taxonomy_all <- Taxonomy_all[-1]



############################
### organize in phyloseq ###
############################

### make objects
TimeSeries_object <- phyloseq(otu_table(as.matrix(TimeSeries_counts), taxa_are_rows = F),
                              tax_table(as.matrix(Taxonomy_all[colnames(TimeSeries_counts), ])),
                              sample_data(TimeSeries_info))

Cities_object <- phyloseq(otu_table(as.matrix(Cities_counts), taxa_are_rows = F),
                          tax_table(as.matrix(Taxonomy_all[colnames(Cities_counts), ])),
                          sample_data(Cities_info))


### keep only bacteria
TimeSeries_object <- subset_taxa(TimeSeries_object, Kingdom == "Bacteria")
Cities_object <- subset_taxa(Cities_object, Kingdom == "Bacteria")

### remove contaminant sequences from time series (see results pt1)
TimeSeries_object  <- subset_taxa(TimeSeries_object, ! taxa_names(TimeSeries_object) %in% c("ASV3695", "ASV3747", "ASV6808"))

### convert to relative abundance
TimeSeries_relabun_object <- transform_sample_counts(TimeSeries_object, function(x) x/sum(x))
TimeSeries_relabun <- data.frame(TimeSeries_relabun_object@otu_table@.Data)
Cities_relabun_object <- transform_sample_counts(Cities_object, function(x) x/sum(x))




######################################
### subset warmest & coldest WWTPs ###
######################################

# 10 warmest and cold cities
hot_cities <- unique(Cities_info[order(Cities_info$Temp_profile_metric, decreasing = TRUE),]$City)[1:10]
cold_cities <- unique(Cities_info[order(Cities_info$Temp_profile_metric, decreasing = FALSE),]$City)[1:10]

# merge phyloseq objects
HotCold_cities_object <- subset_samples(Cities_relabun_object, City %in% hot_cities | City %in% cold_cities)

# merge with time series
HotCold_cities_object <- merge_phyloseq(TimeSeries_relabun_object, HotCold_cities_object)

# add to sample info
sample_data(HotCold_cities_object)$Source[sample_data(HotCold_cities_object)$City %in% hot_cities] <- "Southern city"
sample_data(HotCold_cities_object)$Source[sample_data(HotCold_cities_object)$City %in% cold_cities] <- "Northern city"
sample_data(HotCold_cities_object)$Source[sample_data(HotCold_cities_object)$Sample_name %in% TimeSeries_info$Sample_name] <- "Milwaukee time series"

# PCoA
hotcold_bray <- ordinate(HotCold_cities_object, method = "PCoA", distance = "bray")

# extract results
Axes <- data.frame(hotcold_bray$vectors)[1:2]
Axes$Sample_name <- rownames(Axes)
Axes <- merge(Cities_info[c(1:2, 7)], Axes, by = "Sample_name", all.y = TRUE)

# add info
Axes$City[is.na(Axes$City) == TRUE] <- "Milwaukee"
Axes$Source[Axes$City %in% hot_cities] <- "Southern city"
Axes$Source[Axes$City %in% cold_cities] <- "Northern city"
Axes$Source[Axes$City == "Milwaukee"] <- "Milwaukee time series"

# customize sampling periods
Axes <- merge(Axes, TimeSeries_info[1:2], by = "Sample_name", all.x = TRUE)
Axes$Month <- c(as.character(Axes$Sample_period[1:55]), as.character(Axes$Month[56:149]))
Axes$Month[Axes$Month == "Aug"] <- "August"
Axes$Month[Axes$Month == "Jan"] <- "January"

# what months are in each city sampling period?
Cities_info$Collection_date <- as.Date(Cities_info$Collection_date, format = "%m/%d/%y")
range(subset(Cities_info, Sample_period == "Aug")$Collection_date)
# [1] "2012-08-07" "2012-09-07"

range(subset(Cities_info, Sample_period == "Jan")$Collection_date)
# [1] "2013-01-09" "2013-04-02"

range(subset(Cities_info, Sample_period == "May")$Collection_date)
# [1] "2013-04-28" "2013-06-04"

Axes$Period[Axes$Month == "January" | Axes$Month == "February" |
             Axes$Month == "March" | Axes$Month == "April"] <- "January"
Axes$Period[Axes$Month == "April" | Axes$Month == "May" |
             Axes$Month == "June"] <- "May"
Axes$Period[Axes$Month == "August" | Axes$Month == "September"] <- "August"
Axes <- Axes[which(is.na(Axes$Period) == FALSE), -3]



# extract info for stats
HotCities_relabun <- data.frame(subset_samples(Cities_relabun_object, City %in% hot_cities)@otu_table@.Data)
HotCities_relabun <- HotCities_relabun[rowSums(HotCities_relabun != 0) > 0, colSums(HotCities_relabun != 0) > 0]

ColdCities_relabun <- data.frame(subset_samples(Cities_relabun_object, City %in% cold_cities)@otu_table@.Data)
ColdCities_relabun <- ColdCities_relabun[rowSums(ColdCities_relabun != 0) > 0, colSums(ColdCities_relabun != 0) > 0]

# combine with time series
TimeSeries_relabun.t <- t(TimeSeries_relabun)
TimeSeries_relabun.t <- data.frame(ASV = rownames(TimeSeries_relabun.t), TimeSeries_relabun.t)
TimeSeries_relabun.t[1:5,1:5]

HotCities_relabun.t <- t(HotCities_relabun)
HotCities_relabun.t <- data.frame(ASV = rownames(HotCities_relabun.t), HotCities_relabun.t)
HotCities_relabun.t[1:5,1:5]

ColdCities_relabun.t <- t(ColdCities_relabun)
ColdCities_relabun.t <- data.frame(ASV = rownames(ColdCities_relabun.t), ColdCities_relabun.t)
ColdCities_relabun.t[1:5,1:5]

TS_HotCold_relabun <- merge(TimeSeries_relabun.t, HotCities_relabun.t, by = "ASV", all = TRUE)
TS_HotCold_relabun <- merge(TS_HotCold_relabun, ColdCities_relabun.t, by = "ASV", all = TRUE)
rownames(TS_HotCold_relabun) <- TS_HotCold_relabun$ASV
TS_HotCold_relabun <- data.frame(t(TS_HotCold_relabun[-1]))
TS_HotCold_relabun[1:5,1:5]

TS_HotCold_relabun$Source <- "Milwaukee time series"
TS_HotCold_relabun$Source[rownames(TS_HotCold_relabun) %in% paste0("X", rownames(HotCities_relabun))] <- "Southern city"
TS_HotCold_relabun$Source[rownames(TS_HotCold_relabun) %in% paste0("X", rownames(ColdCities_relabun))] <- "Northern city"


# MANOVA
temp <- subset(TS_HotCold_relabun, Source == "Milwaukee time series" | Source == "Southern city")
TS_hot.adonis <- adonis(temp[,-ncol(temp)] ~ temp$Source, method = "bray", perm = 999, na.rm = TRUE)
TS_hot.adonis
# p 0.001 ***
# R2 0.31244

temp <- subset(TS_HotCold_relabun, Source == "Milwaukee time series" | Source == "Northern city")
TS_cold.adonis <- adonis(temp[,-ncol(temp)] ~ temp$Source, method = "bray", perm = 999, na.rm = TRUE)
TS_cold.adonis
# p 0.001 ***
# R2 0.05062




###########
### map ###
###########

# subset latitudes of examined cities
points <- subset(Cities_info, City %in% hot_cities | City %in% cold_cities | City == "Milwaukee")
points <- points[c(5:2)]
points$source[points$City %in% hot_cities] <- "south"
points$source[points$City %in% cold_cities] <- "north"
points$source[points$City == "Milwaukee"] <- "Milwaukee"

# change order by latitude
points <- points[order(points$Latitude, decreasing = TRUE),]

# only unique cities
points <- unique(points)

# add numbers
points$label <- seq(1:nrow(points))

# transform from lat-long to the format that usmap package uses
points <- usmap_transform(points)

# plot
map <- 
  plot_usmap("states", size = 0.1, fill = "azure3", color = "white") +
  scale_color_manual(values = c("#FEE08B", "#3288BD", "#D53E4F"),
                     labels = c("Milwaukee", "Northern U.S.", "Southern U.S.")) +
  geom_label_repel(data = points, aes(x = Longitude.1, y = Latitude.1, label = label, color = source), 
                  size = 2, fontface = "bold") +
  theme(legend.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 6, color = "black", face = "bold"),
        legend.position = "none")
map

#ggsave("./Plots/map.pdf", plot = map, device = "pdf", width = 3, height = 2, units = "in")





### plot ###
# update with numbered points
Axes <- merge(Axes, points[c(4,6)], by = "City")

rbind(subset(Axes, Source == "Northern city" & Month == "August"),
      subset(Axes, Source == "Southern city" & Month == "January"))

hotcold <- 
  ggplot(Axes, aes(x = Axis.1, y = Axis.2, color = Source, shape = Period)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80", size = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80", size = 0.2) +
  geom_point(size = 2, alpha = 0.8) +
  theme_classic() +
  geom_text_repel(data = rbind(subset(Axes, Source == "Northern city" & Month == "August"),
                               subset(Axes, Source == "Southern city" & Month == "January")), 
                  aes(x = Axis.1, y = Axis.2, label = label), size = 2) +
  scale_color_manual(values = c("#FEE08B", "#3288BD", "#D53E4F"),
                     labels = c("Milwaukee", "Northern U.S.", "Southern U.S.")) +
  scale_shape_manual(values = c(20, 17, 15)) +
  theme(axis.text.x = element_text(size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black"),
        axis.title.x = element_text(size = 6, color = "black"),
        legend.text = element_text(size = 4, color = "black"),
        legend.title = element_text(size = 6, color = "black", face = "bold"),
        axis.line = element_line(size = 0.25),
        axis.ticks = element_line(size = 0.25),
        panel.border = element_rect(color = "grey80", fill = NA, size = 0.5)) +
  guides(color = guide_legend(keyheight = 0.2, keywidth = 0.2, units = "in", ncol = 1),
         shape = guide_legend(keyheight = 0.2, keywidth = 0.2, units = "in")) +
  labs(x = "Axis 1\n(32.5%)", y = "Axis 2\n(9.9%)", color = "Treatment\nplant",
       shape = "Sample\nperiod")
hotcold

#ggsave("./Plots/hotcold.pdf", plot = hotcold, device = "pdf", width = 3.8, height = 2.2, units = "in")
