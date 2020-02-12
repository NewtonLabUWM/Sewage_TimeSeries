###########################################
### Prepare data for LaMartina et al., 2020
### Lou LaMartina, finalized Feb 4, 2020
###########################################


setwd("~/Desktop/TimeSeries_final")
library(phyloseq)
library(decontam)




##########################
### abundance matrices ###
##########################

##############
### wastewater
# all fastq files were processed together in dada2, resulting 35,332 ASVs.
# these abundance matrices were subsetted from the total and saved as CSVs.


# 5yr time series of 2 WWTPs in MKE, WI
TimeSeries_counts <- read.csv("./RData/TimeSeries_counts.csv")
rownames(TimeSeries_counts) <- TimeSeries_counts$Sample_name
TimeSeries_counts <- TimeSeries_counts[-1]


# 72 WWTPs across the US
Cities_counts <- read.csv("./RData/Cities_counts.csv")
rownames(Cities_counts) <- Cities_counts$Sample_name
Cities_counts <- Cities_counts[-1]


# 5 MKE neighborhood sewers
Neighborhood_counts <- read.csv("./RData/Neighborhood_counts.csv")
rownames(Neighborhood_counts) <- Neighborhood_counts$Sample_name
Neighborhood_counts <- Neighborhood_counts[-1]




############################
### Human Microbiome Project
# these are abundances of sequence IDs in each sample pulled from HMP

# load abundance matrix of HMP sequence IDs by samples (9.496255 mins)
HMP_sample_counts <- read.csv("./RData/HMP_sample_counts.csv")
rownames(HMP_sample_counts) <- HMP_sample_counts$SeqID




##########################
### sample information ###
##########################


##############
### wastewater 
# 5yr time series of 2 WWTPs in MKE, WI
TimeSeries_info <- read.csv("./RData/TimeSeries_sample_info.csv")
rownames(TimeSeries_info) <- TimeSeries_info$Sample_name


# 73 WWTPs across the US
Cities_info <- read.csv("./RData/Cities_sample_info.csv")
rownames(Cities_info) <- Cities_info$Sample_name


# 5 MKE neighborhood sewers
Neighborhood_info <- read.csv("./RData/Neighborhood_sample_info.csv")
rownames(Neighborhood_info) <- Neighborhood_info$Sample_name



#######
### HMP

# primer information - how what 16S region was amplified
# (only want v3v5 because we used v4v5 in wastewater)
HMP_primers <- read.csv("./RData/HMP_primers.csv")
rownames(HMP_primers) <- HMP_primers$SeqID


# samples corresponding to which body sites
HMP_bodysites <- read.csv("./RData/HMP_bodysites.csv")
rownames(HMP_bodysites) <- HMP_bodysites$Sample_name




################
### taxonomy ###
################

##############
### wastewater
# taxonomy table from dada2
Taxonomy_all <- read.csv("./RData/Taxonomy_all.csv")
rownames(Taxonomy_all) <- Taxonomy_all$ASV
Taxonomy_all <- Taxonomy_all[-1]



#######
### HMP

# assigned by uploading curated fasta to SILVAngs
HMP_tax <- read.csv("./RData/HMP_taxonomy.csv")
rownames(HMP_tax) <- HMP_tax$SeqID


# subset taxonomy to only samples with counts
HMP_tax <- subset(HMP_tax, SeqID %in% rownames(HMP_sample_counts))


# make abundance and taxonomy tables have same sequence IDs
HMP_sample_counts <- HMP_sample_counts[rownames(HMP_sample_counts) %in% rownames(HMP_tax),]


# make sample counts and body site info match
HMP_sample_counts <- HMP_sample_counts[, colnames(HMP_sample_counts) %in% HMP_bodysites$Sample_name]




############################
### organize in phyloseq ###
############################
# phyloseq stores abundance matrices, taxonomy tables, and sample information into single objects


##############
### wastewater
TimeSeries_object <- phyloseq(otu_table(as.matrix(TimeSeries_counts), taxa_are_rows = F),
                              tax_table(as.matrix(Taxonomy_all[colnames(TimeSeries_counts), ])),
                              sample_data(TimeSeries_info))

Cities_object <- phyloseq(otu_table(as.matrix(Cities_counts), taxa_are_rows = F),
                          tax_table(as.matrix(Taxonomy_all[colnames(Cities_counts), ])),
                          sample_data(Cities_info))

Neighborhood_object <- phyloseq(otu_table(as.matrix(Neighborhood_counts), taxa_are_rows = F),
                                tax_table(as.matrix(Taxonomy_all[colnames(Neighborhood_counts), ])),
                                sample_data(Neighborhood_info))


# remove MINORTH3_15; very low number of reads
Neighborhood_object <- subset_samples(Neighborhood_object, sample_names(Neighborhood_object) != "MINORTH3_15")




#######
### HMP

# make sure taxa names match both abundance and taxonomy tables,
# and sample names match both sample info and abundance tables
HMP_tax <- HMP_tax[order(rownames(HMP_tax)),]

identical(taxa_names(tax_table(as.matrix(HMP_tax[-c(1:2)]))),
          taxa_names(otu_table(as.matrix(HMP_sample_counts), taxa_are_rows = TRUE)))

identical(sample_names(sample_data(HMP_bodysites)), 
          sample_names(otu_table(as.matrix(HMP_sample_counts), taxa_are_rows = TRUE)))

HMP_object <- phyloseq(otu_table(as.matrix(HMP_sample_counts), taxa_are_rows = TRUE),
                       tax_table(as.matrix(HMP_tax[-c(1:2)])),
                       sample_data(HMP_bodysites))




###########################
### remove contaminants ###
###########################
# identify potential contaminating reads that came from the 
# mock community or negative control using "decontam"
# note: only time series run included mock and no template control


# load counts that have negative and mock samples (from dada2)
TimeSeries_wMockNTC_counts <- read.csv("./RData/TimeSeries_wMockNTC.csv")
rownames(TimeSeries_wMockNTC_counts) <- TimeSeries_wMockNTC_counts$Sample_name
TimeSeries_wMockNTC_counts <- TimeSeries_wMockNTC_counts[-1]

TimeSeries_wMockNTC_taxa <- read.csv("./RData/TimeSeries_wMockNTC_tax.csv")
rownames(TimeSeries_wMockNTC_taxa) <- TimeSeries_wMockNTC_taxa$ASV
TimeSeries_wMockNTC_taxa <- TimeSeries_wMockNTC_taxa[-1]


# make sample info
TimeSeries_wMockNTC_info <- data.frame(Sample_name = rownames(TimeSeries_wMockNTC_counts))
rownames(TimeSeries_wMockNTC_info) <- TimeSeries_wMockNTC_info$Sample_name


# make phyloseq object
Decontam_object <- phyloseq(otu_table(as.matrix(TimeSeries_wMockNTC_counts), taxa_are_rows = FALSE),
                            tax_table(as.matrix(TimeSeries_wMockNTC_taxa)),
                            sample_data(TimeSeries_wMockNTC_info))


# add logical variable of if it's negative control or not
sample_data(Decontam_object)$NTC <- sample_data(Decontam_object)$Sample_name == "NTC"
sample_data(Decontam_object)$Mock <- sample_data(Decontam_object)$Sample_name == "Mock"


# run check for contaminants using prevalence method - used for sequencing data
Contam_prev_NTC <- isContaminant(Decontam_object, method = "prevalence", neg = "NTC")
Contam_prev_mock <- isContaminant(Decontam_object, method = "prevalence", neg = "Mock")


# just curious - what taxa are they?
Contam_prev_NTC <- subset(Contam_prev_NTC, contaminant == TRUE)
Contam_prev_mock <- subset(Contam_prev_mock, contaminant == TRUE)
Contam_prev <- rbind(Contam_prev_NTC, Contam_prev_mock)
Contam_prev$ASV <- rownames(Contam_prev)
TimeSeries_wMockNTC_taxa <- data.frame(TimeSeries_wMockNTC_taxa)
TimeSeries_wMockNTC_taxa$ASV <- rownames(TimeSeries_wMockNTC_taxa)
Contaminants <- merge(Contam_prev, TimeSeries_wMockNTC_taxa, by = "ASV")


# also find anything classified as chloroplast, mitochondria, or eukarya.
# V4-V5 are only designed to target bacteria
Badtax <- subset(Taxonomy_all, Family == "Mitochondria" | Order == "Chloroplast" | Kingdom != "Bacteria")


# combine with contaminants
Badtax <- rbind(Badtax, Contaminants[8:14])
Badtax$ASV <- rownames(Badtax)


# remove those from objects
TimeSeries_object <- subset_taxa(TimeSeries_object, ! taxa_names(TimeSeries_object) %in% Badtax$ASV)
Cities_object <- subset_taxa(Cities_object, ! taxa_names(Cities_object) %in% Badtax$ASV)
Neighborhood_object <- subset_taxa(Neighborhood_object, ! taxa_names(Neighborhood_object) %in% Badtax$ASV)




################
### clean up ###
################

# remove empty ASVs
TimeSeries_object <- subset_taxa(TimeSeries_object, taxa_sums(TimeSeries_object) > 0)
Cities_object <- subset_taxa(Cities_object, taxa_sums(Cities_object) > 0)
Neighborhood_object <- subset_taxa(Neighborhood_object, taxa_sums(Neighborhood_object) > 0)
HMP_object <- subset_taxa(HMP_object, taxa_sums(HMP_object) > 0)


# save phyloseq objects, to load into later scripts
saveRDS(TimeSeries_object, "./RData/TimeSeries_phyloseq_object.RData")
saveRDS(Cities_object, "./RData/Cities_phyloseq_object.RData")
saveRDS(Neighborhood_object, "./RData/Neighborhood_phyloseq_object.RData")
saveRDS(HMP_object, "./RData/HMP_phyloseq_object.RData")
