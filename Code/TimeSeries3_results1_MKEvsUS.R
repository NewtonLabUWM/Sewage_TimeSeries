#########################################
### Results section comparing time series
### and US cities microbial communities.
### Lou LaMartina, finalized Feb 4, 2020
#########################################


setwd("~/Desktop/TimeSeries_final")
library(phyloseq)
library(ggplot2)
library(vegan)
library(reshape2)


# load data (refer to TimeSeries1_DataPrep.R)
TimeSeries_object <- readRDS("./RData/TimeSeries_phyloseq_object.RData")
Cities_object <- readRDS("./RData/Cities_phyloseq_object.RData")
Neighborhood_object <- readRDS("./RData/Neighborhood_phyloseq_object.RData")
TimeSeries_info <- read.csv("./RData/TimeSeries_sample_info.csv")
JI_info <- subset(TimeSeries_info, Treatment_plant == "Jones_Island")


# convert to relative abundance
TimeSeries_relabun_object <- transform_sample_counts(TimeSeries_object, function(x) x / sum(x))
Cities_relabun_object <- transform_sample_counts(Cities_object, function(x) x / sum(x))
Neighborhood_relabun_object <- transform_sample_counts(Neighborhood_object, function(x) x / sum(x))


# extract info
JI_relabun <- data.frame(subset_samples(TimeSeries_relabun_object, Treatment_plant == "Jones_Island")@otu_table@.Data)
JI_relabun <- JI_relabun[, colSums(JI_relabun) > 0]
TimeSeries_relabun <- data.frame(TimeSeries_relabun_object@otu_table@.Data)
Cities_relabun <- data.frame(Cities_relabun_object@otu_table@.Data)
Neighborhood_relabun <- data.frame(Neighborhood_relabun_object@otu_table@.Data)



#######################
### alpha diversity ###
#######################

# jones island only
JI.alpha <- data.frame(shannon = diversity(JI_relabun, "shannon"))
JI.alpha$Source <- "Jones_Island"

# jones island and south shore
JI_SS.alpha <- data.frame(shannon = diversity(TimeSeries_relabun, "shannon"))
JI_SS.alpha$Source <- "JI_SS"

# neighborhoods only
Neighborhood_relabun <- data.frame(Neighborhood_relabun_object@otu_table@.Data)
Neigh.alpha <- data.frame(shannon = diversity(Neighborhood_relabun, "shannon"))
Neigh.alpha$Source <- "Neighborhood"

# cities only
Cities.alpha <- data.frame(shannon = diversity(Cities_relabun, "shannon"))
Cities.alpha$Source <- "Cities"

# combine
All_alpha <- rbind(JI.alpha, JI_SS.alpha, Neigh.alpha, Cities.alpha)

# get means
aggregate(. ~ Source, mean, data = All_alpha[-3])

# st dev
aggregate(. ~ Source, sd, data = All_alpha[-3])

# are t-tests appropriate to use? (is the data normally distributed?)
for(i in unique(All_alpha$Source)){
  print(i)
  print(shapiro.test(subset(All_alpha[-3], Source == i)$shannon))
}

# the cities is not. so i will use mann whitney 
# (these samples are independent of each other, i will not use wilcoxon -> paired = F)


# does influent from across US have more diversity than in one city?
wilcox.test(x = subset(All_alpha, Source == "Cities")$shannon,
       y = subset(All_alpha, Source == "JI_SS")$shannon,
       paired = FALSE, alternative = "greater")
# p-value = 0.3335
# no


# do neighborhood sewers have more diversity than any WWTP influent?
wilcox.test(x = subset(All_alpha, Source == "Neighborhood")$shannon,
       y = subset(All_alpha, Source != "Neighborhood")$shannon,
       paired = FALSE, alternative = "greater")
# p-value = 4.897e-11
# yes




######################
### beta diversity ###
######################

# bray curtis scores within Jones Island
JI.bray <- scores(vegdist(JI_relabun, method = "bray"))
JI.bray <- melt(JI.bray)
JI.bray <- JI.bray[JI.bray$value > 0, ]
JI.bray$Source <- "Jones_Island"


# bray curtis scores comparing Jones Island to South Shore
JI_SS.bray <- scores(vegdist(TimeSeries_relabun, method = "bray"))
JI_SS.bray <- JI_SS.bray[grepl("JI", rownames(JI_SS.bray)), grepl("SS", colnames(JI_SS.bray))]
JI_SS.bray <- melt(JI_SS.bray)
JI_SS.bray$Source <- "JI_SS"


# bray curtis scores within neighborhoods
Neigh.bray <- scores(vegdist(Neighborhood_relabun, method = "bray"))
Neigh.bray <- melt(Neigh.bray)
Neigh.bray <- Neigh.bray[Neigh.bray$value > 0,]
Neigh.bray$Source <- "Neighborhood"


# bray curtis scores within Cities
Cities.bray <- scores(vegdist(Cities_relabun, method = "bray"))
Cities.bray <- melt(Cities.bray)
Cities.bray <- Cities.bray[Cities.bray$value > 0,]
Cities.bray$Source <- "Cities"


# bray curtis in same months at jones island
colnames(JI.bray)[1:2] <- c("Sample_name", "Sample_name2")
JI.bray <- merge(JI.bray, JI_info[c(1,3)], by = "Sample_name")
JI_info$Sample_name2 <- JI_info$Sample_name
JI.bray <- merge(JI.bray, JI_info[c(3,29)], by = "Sample_name2")


# loop through and extract values that compared samples in the same months
month.bray <- list()
for(i in JI.bray$Month.x) {
  month.bray[[i]] <- subset(JI.bray, Month.x == i & Month.y == i)$value
  }

month.bray <- data.frame(value = unlist(month.bray))
month.bray$Source <- "Month"


# bray curtis in same years at jones island
JI.bray <- scores(vegdist(JI_relabun, method = "bray"))
JI.bray <- melt(JI.bray)
JI.bray <- JI.bray[JI.bray$value > 0, ]
JI.bray$Source <- "Jones_Island"
colnames(JI.bray)[1:2] <- c("Sample_name", "Sample_name2")
JI.bray <- merge(JI.bray, JI_info[c(1,4)], by = "Sample_name")
JI_info$Sample_name2 <- JI_info$Sample_name
JI.bray <- merge(JI.bray, JI_info[c(4,29)], by = "Sample_name2")


# loop through and extract values that compared samples in the same months
year.bray <- list()
for(i in JI.bray$Year.x) {
  year.bray[[i]] <- subset(JI.bray, Year.x == i & Year.y == i)$value
}

year.bray <- data.frame(value = unlist(year.bray))
year.bray$Source <- "Year"


# combine
All_bray <- rbind(month.bray, year.bray, JI.bray[3:4], JI_SS.bray[3:4], 
                  Neigh.bray[3:4], Cities.bray[3:4])

# means and standard deviations
cbind(aggregate(. ~ Source, mean, data = All_bray), aggregate(. ~ Source, sd, data = All_bray))[-3]
#         Source     value    value.1
# 1       Cities 0.5963390 0.13786574
# 2        JI_SS 0.4615963 0.08547056
# 3 Jones_Island 0.3796046 0.08741774
# 4        Month 0.3095080 0.05822332
# 5 Neighborhood 0.4646036 0.09086621
# 6         Year 0.3688066 0.09263566


# is beta diversity btwn neighborhoods and time series different?
wilcox.test(x = subset(All_bray, Source == "Neighborhood")$Score,
            y = subset(All_bray, Source == "JI_SS")$Score,
            paired = FALSE)
# p-value = 0.116
# nope


# does the time series have greater beta diversity than neighborhoods?
wilcox.test(x = subset(All_bray, Source == "Cities")$Score,
            y = subset(All_bray, Source == "Neighborhood")$Score,
            paired = FALSE, alternative = "greater")
# p-value < 2.2e-16
# yes


                                                       

###########################
### alpha and beta plot ###
###########################

# combine
rownames(All_alpha) <- 1:nrow(All_alpha)
colnames(All_alpha)[1] <- "Score"
All_alpha$Metric <- "Alpha"

All_bray <- All_bray
rownames(All_bray) <- 1:nrow(All_bray)
colnames(All_bray)[1] <- "Score"
All_bray$Metric <- "Beta"

All_diversity <- rbind(All_alpha, All_bray)


### figure 1 ###
div <-
  ggplot(All_diversity, aes(x = Source, y = Score)) +
  geom_boxplot(width = 0.5, outlier.size = 0.5, size = 0.3, outlier.shape = 1, color = "black", fill = "grey90") +
  facet_wrap(Metric ~ ., scales = "free_y", strip.position = "top", ncol = 1,
             labeller = labeller(Metric = c(Alpha = "Within-sample diversity", Beta = "Between-sample diversity"))) +
  scale_x_discrete(limits = c("Month", "Year", "Jones_Island", "JI_SS", "Cities", "Neighborhood"),
                   labels = c(Month = "Jones\nIsland\n(single\nmonth)", Year = "Jones\nIsland\n(single\nyear)",
                              Jones_Island = "Jones\nIsland\n(all\nmonths)", 
                              JI_SS = "Jones\nIsland\nand South\nShore",
                              Cities = "71 U.S.\ncities", Neighborhood = "Milwaukee\nneighborhoods")) +
  scale_y_continuous(labels = function(x) sprintf("%.1f", x)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_line(size = 0.25),
        axis.ticks = element_line(size = 0.25),
        panel.border = element_rect(color = "grey80", fill = NA, size = 0.25),
        panel.spacing = unit(0.1, "inches"),
        strip.text = element_text(size = 6, color = "black", face = "bold"),
        strip.background = element_rect(colour = "white", fill = "white"))
div

#ggsave("./Plots/diversity.pdf", plot = div, device = "pdf", width = 3.1, height = 3, units = "in")
