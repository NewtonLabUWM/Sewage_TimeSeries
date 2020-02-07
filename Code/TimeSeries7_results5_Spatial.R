##########################################
### Results section comparing MKE temporal
### trends to US spatial trends.
### Lou LaMartina,  finalized Feb 4, 2020
##########################################


setwd("~/Desktop/TimeSeries_final")
library(phyloseq)
library(ggplot2)
library(ggrepel)
library(vegan)
library(usmap)
library(ape)


# load data (refer to TimeSeries1_DataPrep.R)
TimeSeries_object <- readRDS("./RData/TimeSeries_phyloseq_object.RData")
Cities_object <- readRDS("./RData/Cities_phyloseq_object.RData")
Cities_info <- read.csv("./RData/Cities_sample_info.csv")
TimeSeries_info <- read.csv("./RData/TimeSeries_sample_info.csv")
Neighborhood_info <- read.csv("./RData/Neighborhood_sample_info.csv")
Cities_object <- subset_samples(Cities_object, City != "Reus")


# convert to relative abundance
TimeSeries_relabun_object <- transform_sample_counts(TimeSeries_object, function(x) x / sum(x))
Cities_relabun_object <- transform_sample_counts(Cities_object, function(x) x / sum(x))


# extract data
TimeSeries_relabun <- data.frame(TimeSeries_relabun_object@otu_table@.Data)


# subset warmest & northest WWTPs
south_cities <- unique(Cities_info[order(Cities_info$Temp_profile_metric, decreasing = TRUE),]$City)[1:10]
north_cities <- unique(Cities_info[order(Cities_info$Temp_profile_metric, decreasing = FALSE),]$City)[1:10]


# merge phyloseq objects
northsouth_cities_object <- subset_samples(Cities_relabun_object, City %in% south_cities | City %in% north_cities)


# merge with time series
northsouth_cities_object <- merge_phyloseq(TimeSeries_relabun_object, northsouth_cities_object)


# add to sample info
sample_data(northsouth_cities_object)$Source[sample_data(northsouth_cities_object)$City %in% south_cities] <- "Southern city"
sample_data(northsouth_cities_object)$Source[sample_data(northsouth_cities_object)$City %in% north_cities] <- "Northern city"
sample_data(northsouth_cities_object)$Source[sample_data(northsouth_cities_object)$Sample_name %in% TimeSeries_info$Sample_name] <- "Milwaukee time series"




###########
### map ###
###########

# subset latitudes of warmest & northest cities
points <- subset(Cities_info, City %in% south_cities | City %in% north_cities | City == "Milwaukee")
points <- unique(points[c(5:2)])
points$source[points$City %in% south_cities] <- "south"
points$source[points$City %in% north_cities] <- "north"
points$source[points$City == "Milwaukee"] <- "Milwaukee"


# change order by latitude
points <- points[order(points$Latitude, decreasing = TRUE),]


# add numbers
points$label <- seq(1:nrow(points))


# transform from lat-long to the format that usmap package uses
points <- usmap_transform(points)


### figure 6A ###
map <- 
  plot_usmap("states", size = 0.1, fill = "azure3", color = "white") +
  scale_color_manual(values = c("#FEE08B", "#3288BD", "#D53E4F"),
                     labels = c("Milwaukee", "Northern U.S.", "Southern U.S.")) +
  geom_text_repel(data = points, aes(x = Longitude.1, y = Latitude.1, label = label, color = source), 
                  size = 2, fontface = "bold", segment.size = 0.3) +
  geom_point(data = points, aes(x = Longitude.1, y = Latitude.1, color = source), 
             size = 0.5) +
  theme(legend.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 6, color = "black", face = "bold"),
        legend.position = "none")
map

#ggsave("./Plots/map.pdf", plot = map, device = "pdf", width = 2.5, height = 1.5, units = "in")




############
### PCoA ###
############

# stat
hotcold_PCoA <- pcoa(vegdist(data.frame(HotCold_cities_object@otu_table@.Data), method = "bray"))


# extract results
Axes <- data.frame(hotcold_PCoA$vectors)[1:2]
Axes$Sample_name <- rownames(Axes)
Axes <- merge(Cities_info[c(1:2, 7)], Axes, by = "Sample_name", all.y = TRUE)


# add info
Axes$City[is.na(Axes$City) == TRUE] <- "Milwaukee"
Axes$Source[Axes$City %in% south_cities] <- "Southern city"
Axes$Source[Axes$City %in% north_cities] <- "Northern city"
Axes$Source[Axes$City == "Milwaukee"] <- "Milwaukee time series"


# customize sampling periods
Axes <- merge(Axes, TimeSeries_info[1:2], by = "Sample_name", all.x = TRUE)
Axes$Month <- c(as.character(Axes$Sample_period[which(is.na(Axes$Sample_period) == FALSE)]), 
                as.character(Axes$Month[which(is.na(Axes$Month) == FALSE)]))
Axes$Month[Axes$Month == "Aug"] <- "August"
Axes$Month[Axes$Month == "Jan"] <- "January"


# make months match sampling periods
Axes$Period[Axes$Month == "January" | Axes$Month == "February" |
             Axes$Month == "March" | Axes$Month == "April"] <- "January"
Axes$Period[Axes$Month == "April" | Axes$Month == "May" |
             Axes$Month == "June"] <- "May"
Axes$Period[Axes$Month == "August" | Axes$Month == "September"] <- "August"
Axes <- Axes[which(is.na(Axes$Period) == FALSE),]


# update with numbered points
Axes <- merge(Axes, points[c("City", "label")], by = "City")


### figure 6B ###
northsouth <- 
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
        panel.border = element_rect(color = "grey80", fill = NA, size = 0.25)) +
  guides(color = guide_legend(keyheight = 0.2, keywidth = 0.2, units = "in", ncol = 1),
         shape = guide_legend(keyheight = 0.2, keywidth = 0.2, units = "in")) +
  labs(x = "Axis 1\n(32.5%)", y = "Axis 2\n(9.9%)", color = "Treatment\nplant",
       shape = "Sample\nperiod")
northsouth

#ggsave("./Plots/northsouth.pdf", plot = northsouth, device = "pdf", width = 4, height = 2.2, units = "in")





#######################################
### which dataset is more variable? ###
#######################################


# extract info
southCities_relabun <- data.frame(subset_samples(Cities_relabun_object, City %in% south_cities)@otu_table@.Data)
southCities_relabun <- southCities_relabun[rowSums(southCities_relabun != 0) > 0, colSums(southCities_relabun != 0) > 0]

northCities_relabun <- data.frame(subset_samples(Cities_relabun_object, City %in% north_cities)@otu_table@.Data)
northCities_relabun <- northCities_relabun[rowSums(northCities_relabun != 0) > 0, colSums(northCities_relabun != 0) > 0]


# combine with time series
TimeSeries_relabun.t <- t(TimeSeries_relabun)
TimeSeries_relabun.t <- data.frame(ASV = rownames(TimeSeries_relabun.t), TimeSeries_relabun.t)

southCities_relabun.t <- t(southCities_relabun)
southCities_relabun.t <- data.frame(ASV = rownames(southCities_relabun.t), southCities_relabun.t)

northCities_relabun.t <- t(northCities_relabun)
northCities_relabun.t <- data.frame(ASV = rownames(northCities_relabun.t), northCities_relabun.t)

TS_northsouth_relabun <- merge(TimeSeries_relabun.t, southCities_relabun.t, by = "ASV", all = TRUE)
TS_northsouth_relabun <- merge(TS_northsouth_relabun, northCities_relabun.t, by = "ASV", all = TRUE)
rownames(TS_northsouth_relabun) <- TS_northsouth_relabun$ASV
TS_northsouth_relabun <- data.frame(t(TS_northsouth_relabun[-1]))


# convert NA to zero (takes a minute)
TS_northsouth_relabun[is.na(TS_northsouth_relabun)] <- 0


# calculate bray curtis distance, then do cluster analysis on that
northsouth.bray <- vegdist(TS_northsouth_relabun[-ncol(TS_northsouth_relabun)], distance = "bray")
northsouth.hclust <- hclust(northsouth.bray, method = "centroid")


# extract info from cluster analysis
merge_height <- data.frame(northsouth.hclust$merge, northsouth.hclust$height, stringsAsFactors = FALSE)
colnames(merge_height) <- c("merge1", "merge2", "height")
merge_height$Compare <- paste0("Compare", seq(1:nrow(merge_height)))
rownames(merge_height) <- merge_height$Compare

labels_order <- data.frame(northsouth.hclust$labels, northsouth.hclust$order, stringsAsFactors = FALSE)
colnames(labels_order) <- c("Sample_name", "order")
labels_order$Sample_name[grepl("X", labels_order$Sample_name)] <- 
  sapply(strsplit(as.character(labels_order$Sample_name[grepl("X", labels_order$Sample_name)]), "X"), '[', 2)


# combine by first comparison
merge_height$x1 <- abs(merge_height$merge1)
labels_order$x1 <- labels_order$order
combine1 <- merge(labels_order, merge_height, by = "x1")[c(2,4,7)]
colnames(combine1)[1] <- "name1"


# combine by second comparison
colnames(labels_order)[3] <- "x2"
merge_height$x2 <- abs(merge_height$merge2)
combine2 <- merge(labels_order, merge_height[-5], by = "x2")[c(2,5,6,7)]
colnames(combine2)[1] <- "name2"


# combine them both
comparisons <- merge(combine1, combine2, by = "Compare")


# subset to types of comparisons
south_comparisons <- subset(comparisons, name1  %in% subset(Axes, Source == "Southern city")$Sample_name &
                            name2 %in% subset(Axes, Source == "Southern city")$Sample_name)
south_comparisons$type <- "south"

north_comparisons <- subset(comparisons, name1  %in% subset(Axes, Source == "Northern city")$Sample_name &
                            name2 %in% subset(Axes, Source == "Northern city")$Sample_name)
north_comparisons$type <- "north"

TS_comparisons <- subset(comparisons, name1  %in% subset(Axes, Source == "Milwaukee time series")$Sample_name &
                           name2 %in% subset(Axes, Source == "Milwaukee time series")$Sample_name)
TS_comparisons$type <- "mke"


# combine for boxplot
comparisons <- rbind(south_comparisons, north_comparisons, TS_comparisons)

clust <-
  ggplot(comparisons, aes(x = type, y = height)) +
  geom_boxplot(width = 0.5, outlier.size = 0.5, size = 0.3, outlier.shape = 1, color = "black", fill = "grey90") +
  theme_classic() +
  scale_x_discrete(limits = c("mke", "north", "south"),
                   labels = c(mke = "Milwaukee\ntime\nseries", north = "Northern\nU.S.",
                              south = "Southern\nU.S.")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black", face = "bold"),
        axis.title.x = element_text(size = 6, color = "black", face = "bold"),
        axis.line = element_line(size = 0.25),
        axis.ticks = element_line(size = 0.25),
        panel.border = element_rect(color = "grey80", fill = NA, size = 0.25)) +
  labs(y = "Centroid cluster analysis height\nusing Bray-Curtis distances", x = "Sample comparisons")
clust

#ggsave("./Plots/clust.pdf", plot = clust, device = "pdf", width = 2, height = 2.2, units = "in")


# which is greater?
shapiro.test(south_comparisons$height) # p-value = 0.7472
shapiro.test(north_comparisons$height) # p-value = 0.4412
shapiro.test(TS_comparisons$height) # p-value = 0.0007935 -> not normally distributed, can't use ttest

wilcox.test(x = TS_comparisons$height,
            y = south_comparisons$height,
            paired = FALSE, alternative = "greater")
# p-value = 0.003764
# time series > south

wilcox.test(x = north_comparisons$height,
            y = south_comparisons$height,
            paired = FALSE, alternative = "greater")
# p-value = 7.648e-07
# north > south

wilcox.test(x = TS_comparisons$height,
            y = north_comparisons$height,
            paired = FALSE)
# p-value = 0.8305
# time series ~ north




# # # # # # # # # # # # #
#  supporting info maps #
# # # # # # # # # # # # #



###############
### USA map ###
###############

# which cities had wastewater collected during which sampling periods?
City_periods <- droplevels(Cities_info[c(2,7)])


# determine the sampling periods covered by the cities 
City_periods.ls <- list()
for(i in as.character(unique(City_periods$City))){
  City_periods.ls[[i]] <- subset(City_periods, City == i)$Sample_period
  City_periods.ls[[i]] <- droplevels(City_periods.ls[[i]])
}


# convert them to true/false
AugJanMay_cities.ls <- list()
for(i in names(City_periods.ls)){
  AugJanMay_cities.ls[[i]] <- c("Aug","Jan","May") %in% City_periods.ls[[i]]
}


# convert to data frame
AugJanMay_cities.df <- data.frame(City = names(AugJanMay_cities.ls), matrix(unlist(AugJanMay_cities.ls), 
                                         nrow = length(AugJanMay_cities.ls), byrow = T))
colnames(AugJanMay_cities.df)[2:4] <- c("Aug", "Jan", "May")


# add variable to which sampling periods were covered
AugJanMay_cities.df$Coverage <- NA

AugJanMay_cities.df$Coverage[AugJanMay_cities.df$Aug == TRUE & 
                               AugJanMay_cities.df$Jan == TRUE &
                               AugJanMay_cities.df$May == TRUE] <- "AugJanMay"

AugJanMay_cities.df$Coverage[AugJanMay_cities.df$Aug == TRUE & 
                               AugJanMay_cities.df$Jan == TRUE &
                               AugJanMay_cities.df$May == FALSE] <- "AugJan"

AugJanMay_cities.df$Coverage[AugJanMay_cities.df$Aug == TRUE & 
                               AugJanMay_cities.df$Jan == FALSE &
                               AugJanMay_cities.df$May == TRUE] <- "AugMay"

AugJanMay_cities.df$Coverage[AugJanMay_cities.df$Aug == FALSE & 
                               AugJanMay_cities.df$Jan == TRUE &
                               AugJanMay_cities.df$May == TRUE] <- "JanMay"

AugJanMay_cities.df$Coverage[AugJanMay_cities.df$Aug == TRUE & 
                               AugJanMay_cities.df$Jan == FALSE &
                               AugJanMay_cities.df$May == FALSE] <- "Aug"

AugJanMay_cities.df$Coverage[AugJanMay_cities.df$Aug == FALSE & 
                               AugJanMay_cities.df$Jan == TRUE &
                               AugJanMay_cities.df$May == FALSE] <- "Jan"

AugJanMay_cities.df$Coverage[AugJanMay_cities.df$Aug == FALSE & 
                               AugJanMay_cities.df$Jan == FALSE &
                               AugJanMay_cities.df$May == TRUE] <- "May"


# add lat, long, state
Cities_points <- Cities_info[2:5]
Cities_points <- Cities_points[match(unique(Cities_points$City), Cities_points$City),]
Cities_points <- merge(AugJanMay_cities.df[c(1,5)], Cities_points, by = "City")


# rearrage - next function needs long first and lat second
Cities_points <- Cities_points[c(5,4,3,1,2)]


# transform lat/long to points
Cities_points <- usmap_transform(Cities_points)


# add numbers to points
Cities_points <- Cities_points[order(Cities_points$Latitude, decreasing = TRUE),]
Cities_points$Label <- seq(1:nrow(Cities_points))
rownames(Cities_points) <- Cities_points$Label


# Brockton, MA translated weird - just make the latitude same as Gloucester
subset(Cities_points, City == "Gloucester")$Latitude.1 # [1] 161912.2
subset(Cities_points, City == "Brockton")$Latitude.1 # [1] 9148306
Cities_points[22,7] # [1] 9148306
Cities_points[22,7] = 161912.2


### figure S1 ###
usmap <- 
  plot_usmap("states", size = 0.1, fill = "white", color = "azure3") +
  geom_text_repel(data = Cities_points, size = 2, fontface = "bold", segment.size = 0.3,
                  aes(label = Label, x = Longitude.1, y = Latitude.1, color = Coverage)) +
  geom_point(data = Cities_points, size = 0.5, aes(x = Longitude.1, y = Latitude.1, color = Coverage)) +
  scale_color_discrete(breaks = c("Aug", "Jan", 
                                  "AugJan", "AugMay", "JanMay",
                                  "AugJanMay")) +
  theme(legend.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 6, color = "black", face = "bold"),
        legend.background = element_rect(fill = NA)) +
  guides(color = guide_legend(keyheight = 0.2, keywidth = 0.2, units = "in", ncol = 1))
usmap

#ggsave("./Plots/usmap.pdf", plot = usmap, device = "pdf", width = 4, height = 3, units = "in")
#write.csv(Cities_points[c(8,1:5)], "./Manuscript/Supporting info/usmap_points.csv", row.names = FALSE)




#########################
### neighborhoods map ###
#########################

# keep only neighborhood name and lat/lon
Neighborhood_points <- Neighborhood_info[6:4]
colnames(Neighborhood_points)[3] <- "Site"


# add WWTP coordinates
WWTP_info <- data.frame(Longitide = c(-87.898647, -87.851993),
                        Latitude = c(43.022104, 42.887431),
                        Site = c("Jones_Island", "South_Shore"))
colnames(WWTP_info) <- colnames(Neighborhood_points)


# combine, add source variables
Neighborhood_points <- rbind(Neighborhood_points, WWTP_info)
Neighborhood_points$Source <- "Neighborhood"
Neighborhood_points$Source[Neighborhood_points$Site == "Jones_Island" |
                             Neighborhood_points$Site == "South_Shore"] <- "WWTP"


# transform lat/long to points
Neighborhood_points <- usmap_transform(Neighborhood_points)
Neighborhood_points$Label <- seq(1:nrow(Neighborhood_points))


### figure S2 ###
mkemap <-
  ggplot(subset(us_map(regions = "counties"), county == "Milwaukee County"), aes(x = x, y = y)) +
  geom_polygon(fill = "white", color = "azure3") +
  geom_text_repel(data = Neighborhood_points, size = 2, fontface = "bold", segment.size = 0.3,
                  aes(label = Label, x = Longitude.1, y = Latitude.1, color = Site)) +
  geom_point(data = Neighborhood_points, size = 1, aes(x = Longitude.1, y = Latitude.1, 
                                                         color = Site, shape = Source)) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 6, color = "black", face = "bold"),
        legend.background = element_rect(fill = NA)) +
  guides(color = guide_legend(keyheight = 0.2, keywidth = 0.2, units = "in", ncol = 1),
         shape = guide_legend(keyheight = 0.2, keywidth = 0.2, units = "in", ncol = 1))
mkemap

#ggsave("./Plots/mkemap.pdf", plot = mkemap, device = "pdf", width = 3, height = 3, units = "in")
