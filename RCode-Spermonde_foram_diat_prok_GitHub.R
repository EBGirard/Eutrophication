#Processing the foram, diatoms and prokaryote datasets
#Author: Elsa B. Girard

options(rstudio.help.showDataPreview = FALSE)

rm(list=ls())
setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
setSessionTimeLimit(cpu = Inf, elapsed = Inf)

#package needed
library(stringr)
library(pheatmap)
library(vegan)
library(ggplot2)
library(reshape2)
library(dplyr)
library(car)
library(patchwork)
library(ggvenn)
library(ggbreak)
library(scales)
library(cowplot)
library(ggpubr)
library(gghalves)
library(tidyverse)
library(rstatix)
library(decontam)
library(phyloseq)
library(viridis)  
library(cooccur)
library(visNetwork)
library(igraph)
library(RColorBrewer) 
library(heatmap3)
library(hrbrthemes)
library(indicspecies)
library(zetadiv)

"%not%" <- Negate("%in%") #to create a "not in" sign, the opposite of %in%
"is.not.na" <- Negate("is.na")

#________________________________

#--------prep foram data--------------------
foram <- read.csv("~downloads/Datasets/forams_esv.csv", row.names = 1)
foram_tax <- read.csv("~downloads/Datasets/forams_taxo.csv", row.names = 1)
foram_met <- read.csv("~downloads/Datasets/METADATA-allNGSsamples_Spermonde2022_Foram.csv", row.names = 1)

#create phyloseq object (make row names for metadata, taxonomy and esv table)
foram <- as.matrix(foram)
foram_tax <- as.matrix(foram_tax)

OTU = otu_table(foram, taxa_are_rows = TRUE)
TAX = tax_table(foram_tax)
samples = sample_data(foram_met)

foram_data <- phyloseq(OTU, TAX, samples)

#________________decontam package to filter from tag switching_____________________

f_neg <- read.csv("~downloads/Datasets/foram_negatives.csv")

sample_data(foram_data)$Sample_or_Control <- NA
sample_data(foram_data)$Sample_or_Control[sample_data(foram_data)$type == "neg_foram"] <- "Control Sample"
sample_data(foram_data)$Sample_or_Control[sample_data(foram_data)$type == "forams"] <- "True Sample"
sample_data(foram_data)$total_reads <- NA
sample_data(foram_data)$total_reads <- colSums(otu_table(foram_data))

#identify which samples with which negatives

#extraction negative 1
extract1 <- subset_samples(foram_data, field.nmbr. %in% f_neg$extract1)
sample_data(extract1)$is.neg <- sample_data(extract1)$Sample_or_Control == "Control Sample"
contamdf.prev1 <- isContaminant(extract1, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev1$contaminant)
contam_extr1 <- taxa_names(extract1)[which(contamdf.prev1$contaminant)]

otu_table(foram_data)[rownames(otu_table(foram_data)) %in% contam_extr1, colnames(otu_table(foram_data)) %in% rownames(sample_data(extract1))] <- 0

# ps.pa <- transform_sample_counts(extract1, function(abund) 1*(abund>0))
# ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
# ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# # Make data.frame of prevalence in positive and negative samples
# df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
#                     contaminant=contamdf.prev1$contaminant)
# 
# plot_frequency(extract1, taxa_names(extract1)[head(which(contamdf.prev1$contaminant),20)], conc="total_reads") + 
#   xlab("DNA Concentration (PicoGreen fluorescent intensity)")

#extraction negative 2
extract2 <- subset_samples(foram_data, field.nmbr. %in% f_neg$extract2)
sample_data(extract2)$is.neg <- sample_data(extract2)$Sample_or_Control == "Control Sample"
contamdf.prev2 <- isContaminant(extract2, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev2$contaminant)
contam_extr2 <- taxa_names(extract2)[which(contamdf.prev2$contaminant)]

otu_table(foram_data)[rownames(otu_table(foram_data)) %in% contam_extr2, colnames(otu_table(foram_data)) %in% rownames(sample_data(extract2))] <- 0


#extraction negative 3
extract3 <- subset_samples(foram_data, field.nmbr. %in% f_neg$extract3)
sample_data(extract3)$is.neg <- sample_data(extract3)$Sample_or_Control == "Control Sample"
contamdf.prev3 <- isContaminant(extract3, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev3$contaminant)
contam_extr3 <- taxa_names(extract3)[which(contamdf.prev3$contaminant)]

otu_table(foram_data)[rownames(otu_table(foram_data)) %in% contam_extr3, colnames(otu_table(foram_data)) %in% rownames(sample_data(extract3))] <- 0


#extraction negative 4
extract4 <- subset_samples(foram_data, field.nmbr. %in% f_neg$extract4)
sample_data(extract4)$is.neg <- sample_data(extract4)$Sample_or_Control == "Control Sample"
contamdf.prev4 <- isContaminant(extract4, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev4$contaminant)
contam_extr4 <- taxa_names(extract4)[which(contamdf.prev4$contaminant)]

otu_table(foram_data)[rownames(otu_table(foram_data)) %in% contam_extr4, colnames(otu_table(foram_data)) %in% rownames(sample_data(extract4))] <- 0


#extraction negative 5
extract5 <- subset_samples(foram_data, field.nmbr. %in% f_neg$extract5)
sample_data(extract5)$is.neg <- sample_data(extract5)$Sample_or_Control == "Control Sample"
contamdf.prev5 <- isContaminant(extract5, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev5$contaminant)
contam_extr5 <- taxa_names(extract5)[which(contamdf.prev5$contaminant)]

otu_table(foram_data)[rownames(otu_table(foram_data)) %in% contam_extr5, colnames(otu_table(foram_data)) %in% rownames(sample_data(extract5))] <- 0


#extraction negative 7
extract7 <- subset_samples(foram_data, field.nmbr. %in% f_neg$extract7)
sample_data(extract7)$is.neg <- sample_data(extract7)$Sample_or_Control == "Control Sample"
contamdf.prev7 <- isContaminant(extract7, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev7$contaminant)
contam_extr7 <- taxa_names(extract7)[which(contamdf.prev7$contaminant)]

otu_table(foram_data)[rownames(otu_table(foram_data)) %in% contam_extr7, colnames(otu_table(foram_data)) %in% rownames(sample_data(extract7))] <- 0


#extraction negative 8
extract8 <- subset_samples(foram_data, field.nmbr. %in% f_neg$extract8)
sample_data(extract8)$is.neg <- sample_data(extract8)$Sample_or_Control == "Control Sample"
contamdf.prev8 <- isContaminant(extract8, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev8$contaminant)
contam_extr8 <- taxa_names(extract8)[which(contamdf.prev8$contaminant)]

otu_table(foram_data)[rownames(otu_table(foram_data)) %in% contam_extr8, colnames(otu_table(foram_data)) %in% rownames(sample_data(extract8))] <- 0


#PCR negative EG009
pcr009 <- subset_samples(foram_data, field.nmbr. %in% f_neg$pcrEG009 & sample_data(foram_data)$PCR2 == "EG009")
sample_data(pcr009)$is.neg <- sample_data(pcr009)$Sample_or_Control == "Control Sample"
contamdf.prev009 <- isContaminant(pcr009, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev009$contaminant)
contam_EG009 <- taxa_names(pcr009)[which(contamdf.prev009$contaminant)]

otu_table(foram_data)[rownames(otu_table(foram_data)) %in% contam_EG009, colnames(otu_table(foram_data)) %in% rownames(sample_data(pcr009))] <- 0


#PCR negative EG010
pcr010 <- subset_samples(foram_data, field.nmbr. %in% f_neg$pcrEG010 & sample_data(foram_data)$PCR2 == "EG010")
sample_data(pcr010)$is.neg <- sample_data(pcr010)$Sample_or_Control == "Control Sample"
contamdf.prev010 <- isContaminant(pcr010, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev010$contaminant)
contam_EG010 <- taxa_names(pcr010)[which(contamdf.prev010$contaminant)]

otu_table(foram_data)[rownames(otu_table(foram_data)) %in% contam_EG010, colnames(otu_table(foram_data)) %in% rownames(sample_data(pcr010))] <- 0


#PCR negative EG011
pcr011 <- subset_samples(foram_data, field.nmbr. %in% f_neg$pcrEG011 & sample_data(foram_data)$PCR2 == "EG011")
sample_data(pcr011)$is.neg <- sample_data(pcr011)$Sample_or_Control == "Control Sample"
contamdf.prev011 <- isContaminant(pcr011, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev011$contaminant)
contam_EG011 <- taxa_names(pcr011)[which(contamdf.prev011$contaminant)]

otu_table(foram_data)[rownames(otu_table(foram_data)) %in% contam_EG011, colnames(otu_table(foram_data)) %in% rownames(sample_data(pcr011))] <- 0

foram <- data.frame(otu_table(foram_data))
foram_met <- data.frame(sample_data(foram_data))
foram_tax <- data.frame(tax_table(foram_data))

#________________Filtering 0.1% + assignment________________________

#df[row,column]
df1 <- foram[rowSums(foram) > 2,]

#transformation switching rows and columns
df1 <- t(df1)
df2 <- as.data.frame(df1)

otudf <- df2
sum(otudf)

#transform raw values into percentages 
df2 <- t(apply(df2, 1, function(x) x/sum(x)))
#make sure row sum equal 1
rowSums(df2)

#remove asv < 0.1% of reads of the total reads but keeping reads number
otudf[df2 < 0.001] <- 0

#keep ESVs that still have reads
otudf <- as.data.frame(otudf[,colSums(otudf) > 0]) #from 7996 esvs to 2407 esvs

sum(otudf)

#number of reads retained after 0.1% quality filter
59828661/70794921*100

#transform back data (switch rows and columns)
datafilt <- as.data.frame(t(otudf))

#create a column with ESVs
datafilt$ESVs <- row.names(datafilt) #filtered data
foram_tax$ESVs <- row.names(foram_tax)

#combine datasets
datafilt1 <- merge(datafilt,foram_tax, by = "ESVs", all.x = TRUE)

#rearrange column order
datafilt2 <- datafilt1[, c(1,356,358,2:353)]

write.csv(datafilt2, "~downloads/Datasets/filtered_forams_esv_taxo.csv")


#________________Add metadata________________________

#melt dataset to creat baseclear column
df_m <- melt(datafilt2, id.vars = c("ESVs", "species", "identity"), 
             value.name = "reads", variable.name = "table_ID")

foram_met$table_ID <- as.character(row.names(foram_met))
df_m$table_ID <- as.character(as.factor(df_m$table_ID))

df2 <- merge(df_m, foram_met, by = "table_ID", all = TRUE)

write.csv(df2, "~downloads/Datasets/filtered_forams_esv_taxo_meta.csv", row.names=FALSE)

#________________get to species level "99.4 % ID)________________________

df2 <- read.csv("~downloads/Datasets/filtered_forams_esv_taxo_meta.csv")
taxo <- read.csv("~downloads/Datasets/taxo_forams_family.csv")

df2 <- df2[!is.na(df2$species),]

df2a <- merge(df2, taxo, by = "species")

df2 <- df2a

df2$family[df2$identity < 96] <- "other_family"

df2$species[df2$identity < 99.4] <- "Z_other_foraminifera"
df2$species[is.na(df2$species)] <- "Z_not_foraminifera"
df2 <- df2[is.not.na(df2$reads),]

df2 <- df2[df2$reads > 0,]

write.csv(df2, "~downloads/Datasets/filtered_forams_esv_taxo_meta_species.csv", row.names=FALSE)

#________________filter esvs that are in at least 2 replicates_________________

foram <- read.csv("~downloads/Datasets/filtered_forams_esv_taxo_meta_species.csv")

df1 <- foram
df1$pa[df1$reads > 0] <- 1
df1$pa[df1$reads == 0] <- 0

df2 <- df1 %>% group_by(ESVs, field.nmbr., island) %>% summarize(counts = sum(pa)) 

df3 <- df2[!(df2$counts  < 2 & df2$island != "negative"), ]

df4 <- merge(df3[,1:3], df1, by = c("ESVs", "field.nmbr.", "island"))
df5 <- df4[,1:(ncol(df4)-1)]

new_df <- df5

new_df$dist_coast_km <- NA

#correct distances in km
new_df$dist_coast_km[new_df$island == "Langkai"] <- 42
new_df$dist_coast_km[new_df$island == "Samalona"] <- 7
new_df$dist_coast_km[new_df$island == "Badi"] <- 23
new_df$dist_coast_km[new_df$island == "K. Keke"] <- 15
new_df$dist_coast_km[new_df$island == "Langkadea"] <- 14
new_df$dist_coast_km[new_df$island == "Karanrang"] <- 14
new_df$dist_coast_km[new_df$island == "Barang Baringan"] <- 5
new_df$dist_coast_km[new_df$island == "Lae Lae"] <- 1
new_df$dist_coast_km[new_df$island == "Lumulumu"] <- 31
new_df$dist_coast_km[new_df$island == "Polewali"] <- 11
new_df$dist_coast_km[new_df$island == "Padjenekang"] <- 19
new_df$dist_coast_km[new_df$island == "Kapoposang Masdar point"] <- 62

write.csv(new_df, "~downloads/Datasets/filtered_forams_esv_taxo_meta_species_replifilt.csv")

#_________merge replicates___________________________

df <- read.csv("~downloads/Datasets/filtered_forams_esv_taxo_meta_species_replifilt.csv")

colnames(df)

df1 <- df %>% 
  group_by(field.nmbr.,island, ESVs, substrate, depth..below.the.local.surface..in.meters.,species, identity, collector, collecting.date..dd.mm.yyyy., country,
           locality,dist_coast_km, family, group) %>%
  summarise(reads = sum(reads))


write.csv(df1, "~downloads/Datasets/merged_forams_data.csv")



#--------prep prok data---------------------
prok <- read.csv("~downloads/Datasets/prok_esv.csv", row.names = 1)
prok_tax <- read.csv("~downloads/Datasets/prok_taxo.csv", row.names = 1)
prok_met <- read.csv("~downloads/Datasets/METADATA-allNGSsamples_Spermonde2022_Prokaryota.csv", row.names = 1)

#create phyloseq object (make row names for metadata, taxonomy and esv table)
prok <- as.matrix(prok[,-732])
prok_tax <- as.matrix(prok_tax)

OTU = otu_table(prok, taxa_are_rows = TRUE)
TAX = tax_table(prok_tax)
samples = sample_data(prok_met)

p_data <- phyloseq(OTU, TAX, samples)

#________________decontam package to filter from tag switching_____________________

p_neg <- read.csv("~downloads/Datasets/Prok_negatives.csv")

sample_data(p_data)$Sample_or_Control <- NA
sample_data(p_data)$Sample_or_Control[sample_data(p_data)$Substrate.type %in% c("PCR NEGATIVE","Extraction NEGATIVE")] <- "Control Sample"
sample_data(p_data)$Sample_or_Control[sample_data(p_data)$Substrate.type %not% c("PCR NEGATIVE","Extraction NEGATIVE")] <- "True Sample"
sample_data(p_data)$total_reads <- NA
sample_data(p_data)$total_reads <- colSums(otu_table(p_data))

#NEGATIVES FOR SUBSTRATES - EXTRACTION
#extraction negative 1
extract1 <- subset_samples(p_data, field.nmbr. %in% p_neg$extract1)
sample_data(extract1)$is.neg <- sample_data(extract1)$Sample_or_Control == "Control Sample"
contamdf.prev1 <- isContaminant(extract1, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev1$contaminant)
contam_extr1 <- taxa_names(extract1)[which(contamdf.prev1$contaminant)]

otu_table(p_data)[rownames(otu_table(p_data)) %in% contam_extr1, colnames(otu_table(p_data)) %in% rownames(sample_data(extract1))] <- 0
#extraction negative 2
extract2 <- subset_samples(p_data, field.nmbr. %in% p_neg$extract2)
sample_data(extract2)$is.neg <- sample_data(extract2)$Sample_or_Control == "Control Sample"
contamdf.prev2 <- isContaminant(extract2, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev2$contaminant)
contam_extr2 <- taxa_names(extract2)[which(contamdf.prev2$contaminant)]

otu_table(p_data)[rownames(otu_table(p_data)) %in% contam_extr2, colnames(otu_table(p_data)) %in% rownames(sample_data(extract2))] <- 0
#extraction negative 3
extract3 <- subset_samples(p_data, field.nmbr. %in% p_neg$extract3)
sample_data(extract3)$is.neg <- sample_data(extract3)$Sample_or_Control == "Control Sample"
contamdf.prev3 <- isContaminant(extract3, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev3$contaminant)
contam_extr3 <- taxa_names(extract3)[which(contamdf.prev3$contaminant)]

otu_table(p_data)[rownames(otu_table(p_data)) %in% contam_extr3, colnames(otu_table(p_data)) %in% rownames(sample_data(extract3))] <- 0
#extraction negative 4
extract4 <- subset_samples(p_data, field.nmbr. %in% p_neg$extract4)
sample_data(extract4)$is.neg <- sample_data(extract4)$Sample_or_Control == "Control Sample"
contamdf.prev4 <- isContaminant(extract4, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev4$contaminant)
contam_extr4 <- taxa_names(extract4)[which(contamdf.prev4$contaminant)]

otu_table(p_data)[rownames(otu_table(p_data)) %in% contam_extr4, colnames(otu_table(p_data)) %in% rownames(sample_data(extract4))] <- 0
#extraction negative 5
extract5 <- subset_samples(p_data, field.nmbr. %in% p_neg$extract5)
sample_data(extract5)$is.neg <- sample_data(extract5)$Sample_or_Control == "Control Sample"
contamdf.prev5 <- isContaminant(extract5, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev5$contaminant)
contam_extr5 <- taxa_names(extract5)[which(contamdf.prev5$contaminant)]

otu_table(p_data)[rownames(otu_table(p_data)) %in% contam_extr5, colnames(otu_table(p_data)) %in% rownames(sample_data(extract5))] <- 0
#extraction negative 6
extract6 <- subset_samples(p_data, field.nmbr. %in% p_neg$extract6)
sample_data(extract6)$is.neg <- sample_data(extract6)$Sample_or_Control == "Control Sample"
contamdf.prev6 <- isContaminant(extract6, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev6$contaminant)
contam_extr6 <- taxa_names(extract6)[which(contamdf.prev6$contaminant)]

otu_table(p_data)[rownames(otu_table(p_data)) %in% contam_extr6, colnames(otu_table(p_data)) %in% rownames(sample_data(extract6))] <- 0

#NEGATIVES FOR SUBSTRATES - PCR
#pcr negative LR029
pcr029 <- subset_samples(p_data, field.nmbr. %in% p_neg$pcrLR029 & PCR2_bacteria == "LR029")
sample_data(pcr029)$is.neg <- sample_data(pcr029)$Sample_or_Control == "Control Sample"
contamdf.prev029 <- isContaminant(pcr029, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev029$contaminant)
contam_pcr029 <- taxa_names(pcr029)[which(contamdf.prev029$contaminant)]

otu_table(p_data)[rownames(otu_table(p_data)) %in% contam_pcr029, colnames(otu_table(p_data)) %in% rownames(sample_data(pcr029))] <- 0
#pcr negative LR030
pcr030 <- subset_samples(p_data, field.nmbr. %in% p_neg$pcrLR030 & PCR2_bacteria == "LR030")
sample_data(pcr030)$is.neg <- sample_data(pcr030)$Sample_or_Control == "Control Sample"
contamdf.prev030 <- isContaminant(pcr030, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev030$contaminant)
contam_pcr030 <- taxa_names(pcr030)[which(contamdf.prev030$contaminant)]

otu_table(p_data)[rownames(otu_table(p_data)) %in% contam_pcr030, colnames(otu_table(p_data)) %in% rownames(sample_data(pcr030))] <- 0
#pcr negative LR031
pcr031 <- subset_samples(p_data, field.nmbr. %in% p_neg$pcrLR031 & PCR2_bacteria == "LR031")
sample_data(pcr031)$is.neg <- sample_data(pcr031)$Sample_or_Control == "Control Sample"
contamdf.prev031 <- isContaminant(pcr031, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev031$contaminant)
contam_pcr031 <- taxa_names(pcr031)[which(contamdf.prev031$contaminant)]

otu_table(p_data)[rownames(otu_table(p_data)) %in% contam_pcr031, colnames(otu_table(p_data)) %in% rownames(sample_data(pcr031))] <- 0
#pcr negative LR033
pcrLR033 <- subset_samples(p_data, field.nmbr. %in% p_neg$pcrLR033 & PCR2_bacteria == "LR033")
sample_data(pcrLR033)$is.neg <- sample_data(pcrLR033)$Sample_or_Control == "Control Sample"
contamdf.prev033 <- isContaminant(pcrLR033, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev033$contaminant)
contam_pcrLR033 <- taxa_names(pcrLR033)[which(contamdf.prev033$contaminant)]

otu_table(p_data)[rownames(otu_table(p_data)) %in% contam_pcrLR033, colnames(otu_table(p_data)) %in% rownames(sample_data(pcrLR033))] <- 0


#NEGATIVES FOR WATER - EXTRACTION&PCR
#pcr negative LR034
pcr034 <- subset_samples(p_data, field.nmbr. %in% p_neg$pcrLR034 & PCR2_bacteria == "LR034")
sample_data(pcr034)$is.neg <- sample_data(pcr034)$Sample_or_Control == "Control Sample"
contamdf.prev034 <- isContaminant(pcr034, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev034$contaminant)
contam_pcr034 <- taxa_names(pcr034)[which(contamdf.prev034$contaminant)]

otu_table(p_data)[rownames(otu_table(p_data)) %in% contam_pcr034, colnames(otu_table(p_data)) %in% rownames(sample_data(pcr034))] <- 0
#pcr negative LR027
pcr027 <- subset_samples(p_data, field.nmbr. %in% p_neg$pcrLR027 & PCR2_bacteria == "LR027")
sample_data(pcr027)$is.neg <- sample_data(pcr027)$Sample_or_Control == "Control Sample"
contamdf.prev027 <- isContaminant(pcr027, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev027$contaminant)
contam_pcr027 <- taxa_names(pcr027)[which(contamdf.prev027$contaminant)]

otu_table(p_data)[rownames(otu_table(p_data)) %in% contam_pcr027, colnames(otu_table(p_data)) %in% rownames(sample_data(pcr027))] <- 0
#pcr negative LR028
pcr028 <- subset_samples(p_data, field.nmbr. %in% p_neg$prcLR028 & PCR2_bacteria == "LR028")
sample_data(pcr028)$is.neg <- sample_data(pcr028)$Sample_or_Control == "Control Sample"
contamdf.prev028 <- isContaminant(pcr028, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev028$contaminant)
contam_pcr028 <- taxa_names(pcr028)[which(contamdf.prev028$contaminant)]

otu_table(p_data)[rownames(otu_table(p_data)) %in% contam_pcr028, colnames(otu_table(p_data)) %in% rownames(sample_data(pcr028))] <- 0
#pcr negative LR032
pcr032 <- subset_samples(p_data, field.nmbr. %in% p_neg$pcrLR032 & PCR2_bacteria == "LR032")
sample_data(pcr032)$is.neg <- sample_data(pcr032)$Sample_or_Control == "Control Sample"
contamdf.prev032 <- isContaminant(pcr032, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev032$contaminant)
contam_pcr032 <- taxa_names(pcr032)[which(contamdf.prev032$contaminant)]

otu_table(p_data)[rownames(otu_table(p_data)) %in% contam_pcr032, colnames(otu_table(p_data)) %in% rownames(sample_data(pcr032))] <- 0

prok <- data.frame(otu_table(p_data))
prok_met <- data.frame(sample_data(p_data))
prok_tax <- data.frame(tax_table(p_data))

#________________Filtering 0.1% + assignment________________________

#df[row,column]
df1 <- prok[rowSums(prok) > 2,]

#transformation switching rows and columns
df1 <- t(df1)
df2 <- as.data.frame(df1)

#remove samples with less than 1000 reads
otudf <- df2
sum(otudf)

#transform raw values into percentages 
df2 <- t(apply(df2, 1, function(x) x/sum(x)))
#make sure row sum equal 1
rowSums(df2)

#remove asv < 0.1% of reads of the total reads but keeping reads number
otudf[df2 < 0.001] <- 0

#keep ESVs that still have reads
otudf <- as.data.frame(otudf[,colSums(otudf) > 0]) #from 94858 esvs to 8468 esvs

sum(otudf)

#number of reads retained after 0.1% quality filter
57675403/87402769*100

#transform back data (switch rows and columns)
datafilt <- as.data.frame(t(otudf))

#create a column with ESVs
datafilt$ESVs <- row.names(datafilt) #filtered data
prok_tax$ESVs <- row.names(prok_tax)

#combine datasets
datafilt1 <- merge(datafilt,prok_tax, by = "ESVs", all.x = TRUE)

#rearrange column order
datafilt2 <- datafilt1[, c(1,741,736,2:731)]

write.csv(datafilt2, "~downloads/Datasets/filtered_prok_esv_taxo.csv")


#________________Add metadata________________________

#melt dataset to creat baseclear column
df_m <- melt(datafilt2, id.vars = c("ESVs", "Taxonomy", "Identity.percentage"), 
             value.name = "reads", variable.name = "real_ID")

df_m$real_ID <- as.character(df_m$real_ID)
prok_met$real_ID <- as.character(row.names(prok_met))

df2 <- merge(df_m, prok_met, by = "real_ID", all = TRUE)

write.csv(df2, "~downloads/Datasets/filtered_prok_esv_taxo_meta.csv", row.names=FALSE)


#________________taxonomy________________________

df2 <- read.csv("~downloads/Datasets/filtered_prok_esv_taxo_meta.csv")


splitspecies <- str_split_fixed(df2$Taxonomy, " / ", 7)
colnames(splitspecies) <- c("Kingdom", "Phylum", "Class","Order", "Family", "Genus", "Species")
df2 <- cbind(df2, splitspecies)

df2 <- df2[df2$reads > 0,]

write.csv(df2, "~downloads/Datasets/filtered_prok_esv_taxo_meta_species.csv", row.names=FALSE)

#________________filter esvs that are in at least 2 replicates_________________

prok <- read.csv("~downloads/Datasets/filtered_prok_esv_taxo_meta_species.csv")

df1 <- prok
df1$pa[df1$reads > 0] <- 1
df1$pa[df1$reads == 0] <- 0

df2 <- df1 %>% group_by(ESVs, Substrate.type, island, depth..below.the.local.surface..in.meters.) %>% summarize(counts = sum(pa)) 

df3 <- df2[!(df2$counts  < 2), ]

df4 <- merge(df3[,1:4], df1, by = c("ESVs", "Substrate.type", "island", "depth..below.the.local.surface..in.meters."))
df5 <- df4[,1:(ncol(df4)-1)]

new_df <- df5

new_df$dist_coast_km <- NA

#correct distances in km
new_df$dist_coast_km[new_df$island == "Langkai"] <- 42
new_df$dist_coast_km[new_df$island == "Samalona"] <- 7
new_df$dist_coast_km[new_df$island == "Badi"] <- 23
new_df$dist_coast_km[new_df$island == "K. Keke"] <- 15
new_df$dist_coast_km[new_df$island == "Langkadea"] <- 14
new_df$dist_coast_km[new_df$island == "Karanrang"] <- 14
new_df$dist_coast_km[new_df$island == "Barang Baringan"] <- 5
new_df$dist_coast_km[new_df$island == "Lae Lae"] <- 1
new_df$dist_coast_km[new_df$island == "Lumulumu"] <- 31
new_df$dist_coast_km[new_df$island == "Polewali"] <- 11
new_df$dist_coast_km[new_df$island == "Padjenekang"] <- 19
new_df$dist_coast_km[new_df$island == "Kapoposang Masdar point"] <- 62




write.csv(new_df, "~downloads/Datasets/filtered_prok_esv_taxo_meta_species_replifilt.csv")

#_________remove chloroplast + mitochondria and organelle___________________

df <- read.csv("~downloads/Datasets/filtered_prok_esv_taxo_meta_species_replifilt.csv")

new_df <- df[df$Order != "Chloroplast",]
new_df <- new_df[new_df$Family %not% "Mitochondria",]

write.csv(new_df, "~downloads/Datasets/filtered_prok_esv_taxo_meta_species_replifilt.csv")

#_________merge replicates___________________________

df <- read.csv("~downloads/Datasets/filtered_prok_esv_taxo_meta_species_replifilt.csv")

colnames(df)

df1 <- df %>% 
  group_by(field.nmbr.,island, ESVs, Substrate.type,substrate, depth..below.the.local.surface..in.meters.,Taxonomy, Identity.percentage, collector, collecting.date..dd.mm.yyyy., country, state_province,
           locality,Kingdom,Phylum, Class,
           Order,Family,Genus,Species,dist_coast_km) %>%
  summarise(reads = sum(reads))


write.csv(df1, "~downloads/Datasets/merged_prok_data.csv")



#--------prep diat data---------------------
diat <- read.csv("~downloads/Datasets/diat_esv.csv", row.names = 1)
diat_tax <- read.csv("~downloads/Datasets/diat_taxo.csv", row.names = 1)
diat_met <- read.csv("~downloads/Datasets/METADATA-allNGSsamples_Spermonde2022_Diatoms.csv", row.names = 1)

#create phyloseq object (make row names for metadata, taxonomy and esv table)
diat <- as.matrix(diat[,-731])
diat_tax <- as.matrix(diat_tax)

OTU = otu_table(diat, taxa_are_rows = TRUE)
TAX = tax_table(diat_tax)
samples = sample_data(diat_met)

d_data <- phyloseq(OTU, TAX, samples)

#________________decontam package to filter from tag switching_____________________

d_neg <- read.csv("~downloads/Datasets/Diat_negatives.csv")

sample_data(d_data)$Sample_or_Control <- NA
sample_data(d_data)$Sample_or_Control[sample_data(d_data)$Substrate.type %in% c("PCR NEGATIVE","Extraction NEGATIVE")] <- "Control Sample"
sample_data(d_data)$Sample_or_Control[sample_data(d_data)$Substrate.type %not% c("PCR NEGATIVE","Extraction NEGATIVE")] <- "True Sample"
sample_data(d_data)$total_reads <- NA
sample_data(d_data)$total_reads <- colSums(otu_table(d_data))

#NEGATIVES FOR SUBSTRATES - EXTRACTION
#extraction negative 1
extract1 <- subset_samples(d_data, field.nmbr. %in% d_neg$extract1)
sample_data(extract1)$is.neg <- sample_data(extract1)$Sample_or_Control == "Control Sample"
contamdf.prev1 <- isContaminant(extract1, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev1$contaminant)
contam_extr1 <- taxa_names(extract1)[which(contamdf.prev1$contaminant)]

otu_table(d_data)[rownames(otu_table(d_data)) %in% contam_extr1, colnames(otu_table(d_data)) %in% rownames(sample_data(extract1))] <- 0
#extraction negative 2
extract2 <- subset_samples(d_data, field.nmbr. %in% d_neg$extract2)
sample_data(extract2)$is.neg <- sample_data(extract2)$Sample_or_Control == "Control Sample"
contamdf.prev2 <- isContaminant(extract2, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev2$contaminant)
contam_extr2 <- taxa_names(extract2)[which(contamdf.prev2$contaminant)]

otu_table(d_data)[rownames(otu_table(d_data)) %in% contam_extr2, colnames(otu_table(d_data)) %in% rownames(sample_data(extract2))] <- 0
#extraction negative 3
extract3 <- subset_samples(d_data, field.nmbr. %in% d_neg$extract3)
sample_data(extract3)$is.neg <- sample_data(extract3)$Sample_or_Control == "Control Sample"
contamdf.prev3 <- isContaminant(extract3, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev3$contaminant)
contam_extr3 <- taxa_names(extract3)[which(contamdf.prev3$contaminant)]

otu_table(d_data)[rownames(otu_table(d_data)) %in% contam_extr3, colnames(otu_table(d_data)) %in% rownames(sample_data(extract3))] <- 0
#extraction negative 4
extract4 <- subset_samples(d_data, field.nmbr. %in% d_neg$extract4)
sample_data(extract4)$is.neg <- sample_data(extract4)$Sample_or_Control == "Control Sample"
contamdf.prev4 <- isContaminant(extract4, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev4$contaminant)
contam_extr4 <- taxa_names(extract4)[which(contamdf.prev4$contaminant)]

otu_table(d_data)[rownames(otu_table(d_data)) %in% contam_extr4, colnames(otu_table(d_data)) %in% rownames(sample_data(extract4))] <- 0
#extraction negative 5
extract5 <- subset_samples(d_data, field.nmbr. %in% d_neg$extract5)
sample_data(extract5)$is.neg <- sample_data(extract5)$Sample_or_Control == "Control Sample"
contamdf.prev5 <- isContaminant(extract5, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev5$contaminant)
contam_extr5 <- taxa_names(extract5)[which(contamdf.prev5$contaminant)]

otu_table(d_data)[rownames(otu_table(d_data)) %in% contam_extr5, colnames(otu_table(d_data)) %in% rownames(sample_data(extract5))] <- 0
#extraction negative 6
extract6 <- subset_samples(d_data, field.nmbr. %in% d_neg$extract6)
sample_data(extract6)$is.neg <- sample_data(extract6)$Sample_or_Control == "Control Sample"
contamdf.prev6 <- isContaminant(extract6, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev6$contaminant)
contam_extr6 <- taxa_names(extract6)[which(contamdf.prev6$contaminant)]

otu_table(d_data)[rownames(otu_table(d_data)) %in% contam_extr6, colnames(otu_table(d_data)) %in% rownames(sample_data(extract6))] <- 0

#NEGATIVES FOR SUBSTRATES - PCR
#pcr negative SV019
pcr019 <- subset_samples(d_data, field.nmbr. %in% d_neg$pcrSV019 & PCR2_diatom == "SV019")
sample_data(pcr019)$is.neg <- sample_data(pcr019)$Sample_or_Control == "Control Sample"
contamdf.prev019 <- isContaminant(pcr019, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev019$contaminant)
contam_pcr019 <- taxa_names(pcr019)[which(contamdf.prev019$contaminant)]

otu_table(d_data)[rownames(otu_table(d_data)) %in% contam_pcr019, colnames(otu_table(d_data)) %in% rownames(sample_data(pcr019))] <- 0
#pcr negative SV020
pcr020 <- subset_samples(d_data, field.nmbr. %in% d_neg$pcrSV020 & PCR2_diatom == "SV020")
sample_data(pcr020)$is.neg <- sample_data(pcr020)$Sample_or_Control == "Control Sample"
contamdf.prev020 <- isContaminant(pcr020, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev020$contaminant)
contam_pcr020 <- taxa_names(pcr020)[which(contamdf.prev020$contaminant)]

otu_table(d_data)[rownames(otu_table(d_data)) %in% contam_pcr020, colnames(otu_table(d_data)) %in% rownames(sample_data(pcr020))] <- 0
#pcr negative SV021
pcr021 <- subset_samples(d_data, field.nmbr. %in% d_neg$pcrSV021 & PCR2_diatom == "SV021")
sample_data(pcr021)$is.neg <- sample_data(pcr021)$Sample_or_Control == "Control Sample"
contamdf.prev021 <- isContaminant(pcr021, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev021$contaminant)
contam_pcr021 <- taxa_names(pcr021)[which(contamdf.prev021$contaminant)]

otu_table(d_data)[rownames(otu_table(d_data)) %in% contam_pcr021, colnames(otu_table(d_data)) %in% rownames(sample_data(pcr021))] <- 0
#pcr negative TM018
pcrTM018 <- subset_samples(d_data, field.nmbr. %in% d_neg$pcrTM018 & PCR2_diatom == "TM018")
sample_data(pcrTM018)$is.neg <- sample_data(pcrTM018)$Sample_or_Control == "Control Sample"
contamdf.prev018 <- isContaminant(pcrTM018, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev018$contaminant)
contam_pcrTM018 <- taxa_names(pcrTM018)[which(contamdf.prev018$contaminant)]

otu_table(d_data)[rownames(otu_table(d_data)) %in% contam_pcrTM018, colnames(otu_table(d_data)) %in% rownames(sample_data(pcrTM018))] <- 0


#NEGATIVES FOR WATER - EXTRACTION&PCR
#pcr negative SV023
pcr023 <- subset_samples(d_data, field.nmbr. %in% d_neg$pcrSV023 & PCR2_diatom == "SV023")
sample_data(pcr023)$is.neg <- sample_data(pcr023)$Sample_or_Control == "Control Sample"
contamdf.prev023 <- isContaminant(pcr023, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev023$contaminant)
contam_pcr023 <- taxa_names(pcr023)[which(contamdf.prev023$contaminant)]

otu_table(d_data)[rownames(otu_table(d_data)) %in% contam_pcr023, colnames(otu_table(d_data)) %in% rownames(sample_data(pcr023))] <- 0
#pcr negative SV016
pcr016 <- subset_samples(d_data, field.nmbr. %in% d_neg$pcrSV016 & PCR2_diatom == "SV016")
sample_data(pcr016)$is.neg <- sample_data(pcr016)$Sample_or_Control == "Control Sample"
contamdf.prev016 <- isContaminant(pcr016, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev016$contaminant)
contam_pcr016 <- taxa_names(pcr016)[which(contamdf.prev016$contaminant)]

otu_table(d_data)[rownames(otu_table(d_data)) %in% contam_pcr016, colnames(otu_table(d_data)) %in% rownames(sample_data(pcr016))] <- 0
#pcr negative SV017
pcr017 <- subset_samples(d_data, field.nmbr. %in% d_neg$prcSV017 & PCR2_diatom == "SV017")
sample_data(pcr017)$is.neg <- sample_data(pcr017)$Sample_or_Control == "Control Sample"
contamdf.prev017 <- isContaminant(pcr017, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev017$contaminant)
contam_pcr017 <- taxa_names(pcr017)[which(contamdf.prev017$contaminant)]

otu_table(d_data)[rownames(otu_table(d_data)) %in% contam_pcr017, colnames(otu_table(d_data)) %in% rownames(sample_data(pcr017))] <- 0
#pcr negative SV018
pcr018 <- subset_samples(d_data, field.nmbr. %in% d_neg$pcrSV018 & PCR2_diatom == "SV018")
sample_data(pcr018)$is.neg <- sample_data(pcr018)$Sample_or_Control == "Control Sample"
contamdf.prev018 <- isContaminant(pcr018, method="prevalence", neg="is.neg",threshold=0.2)

table(contamdf.prev018$contaminant)
contam_pcr018 <- taxa_names(pcr018)[which(contamdf.prev018$contaminant)]

otu_table(d_data)[rownames(otu_table(d_data)) %in% contam_pcr018, colnames(otu_table(d_data)) %in% rownames(sample_data(pcr018))] <- 0

diat <- data.frame(otu_table(d_data))
diat_met <- data.frame(sample_data(d_data))
diat_tax <- data.frame(tax_table(d_data))


#________________Filtering 0.1% + assignment________________________

#df[row,column]
df1 <- diat[rowSums(diat) > 2,]

#transformation switching rows and columns
df1 <- t(df1)
df2 <- as.data.frame(df1)

otudf <- df2
sum(otudf)

#transform raw values into percentages 
df2 <- t(apply(df2, 1, function(x) x/sum(x)))
#make sure row sum equal 1
rowSums(df2)

#remove asv < 0.1% of reads of the total reads but keeping reads number
otudf[df2 < 0.001] <- 0

#keep ESVs that still have reads
otudf <- as.data.frame(otudf[,colSums(otudf) > 0]) #from 22193 esvs to 10046 esvs

sum(otudf)

#number of reads retained after 0.1% quality filter
47403615/54193779*100

#transform back data (switch rows and columns)
datafilt <- as.data.frame(t(otudf))

#create a column with ESVs
datafilt$ESVs <- row.names(datafilt) #filtered data
diat_tax$ESVs <- row.names(diat_tax)

#combine datasets
datafilt1 <- merge(datafilt,diat_tax, by = "ESVs", all.x = TRUE)

#rearrange column order
datafilt2 <- datafilt1[, c(1,740,735,2:731)]

write.csv(datafilt2, "~downloads/Datasets/filtered_diat_esv_taxo.csv")


#________________Add metadata________________________

#melt dataset to creat baseclear column
df_m <- melt(datafilt2, id.vars = c("ESVs", "Taxonomy", "Identity.percentage"), 
             value.name = "reads", variable.name = "real_ID")

diat_met$real_ID <- as.character(row.names(diat_met))
df_m$real_ID <- as.character(as.factor(df_m$real_ID))

df2 <- merge(df_m, diat_met, by = "real_ID", all = TRUE)

write.csv(df2, "~downloads/Datasets/filtered_diat_esv_taxo_meta.csv", row.names=FALSE)
df2<- read.csv("~downloads/Datasets/filtered_diat_esv_taxo_meta.csv")

#________________taxonomy________________________from R-syst::diatom website
#cite diatom r-syst:
#http://dx.doi.org/10.5281/zenodo.31137
#http://database.oxfordjournals.org/content/2016/baw016.full?keytype=ref&ijkey=H324uA95JzzEomz

df2 <- df2[df2$reads > 0,]

df2$Taxonomy[is.na(df2$Taxonomy)] <- "unknown kingdom/unknown phylum/unknown class/unknown order/unknown family/unknown genus/unknown species"
splitspecies <- str_split_fixed(df2$Taxonomy, "/", 7)
colnames(splitspecies) <- c("Kindgom", "Phylum", "Class","Order", "Family", "Genus", "Species")
df2 <- cbind(df2, splitspecies)

write.csv(df2, "~downloads/Datasets/filtered_diat_esv_taxo_meta_species.csv", row.names=FALSE)

#________________filter esvs that are in at least 2 replicates_________________

diat <- read.csv("~downloads/Datasets/filtered_diat_esv_taxo_meta_species.csv")

df1 <- diat
df1$pa[df1$reads > 0] <- 1
df1$pa[df1$reads == 0] <- 0

df2 <- df1 %>% group_by(ESVs, Substrate.type, island, depth..below.the.local.surface..in.meters.) %>% summarize(counts = sum(pa)) 

df3 <- df2[!(df2$counts  < 2), ]

df4 <- merge(df3[,1:4], df1, by = c("ESVs", "Substrate.type", "island", "depth..below.the.local.surface..in.meters."))
df5 <- df4[,1:(ncol(df4)-1)]

new_df <- df5

new_df$dist_coast_km <- NA

#correct distances in km
new_df$dist_coast_km[new_df$island == "Langkai"] <- 42
new_df$dist_coast_km[new_df$island == "Samalona"] <- 7
new_df$dist_coast_km[new_df$island == "Badi"] <- 23
new_df$dist_coast_km[new_df$island == "K. Keke"] <- 15
new_df$dist_coast_km[new_df$island == "Langkadea"] <- 14
new_df$dist_coast_km[new_df$island == "Karanrang"] <- 14
new_df$dist_coast_km[new_df$island == "Barang Baringan"] <- 5
new_df$dist_coast_km[new_df$island == "Lae Lae"] <- 1
new_df$dist_coast_km[new_df$island == "Lumulumu"] <- 31
new_df$dist_coast_km[new_df$island == "Polewali"] <- 11
new_df$dist_coast_km[new_df$island == "Padjenekang"] <- 19
new_df$dist_coast_km[new_df$island == "Kapoposang Masdar point"] <- 62

write.csv(new_df, "~downloads/Datasets/filtered_diat_esv_taxo_meta_species_replifilt.csv")

#_________merge replicates___________________________

df <- read.csv("~downloads/Datasets/filtered_diat_esv_taxo_meta_species_replifilt.csv")

df1 <- df %>% 
  group_by(field.nmbr.,island, ESVs, Substrate.type,substrate, depth..below.the.local.surface..in.meters.,Taxonomy, Identity.percentage, collector, collecting.date..dd.mm.yyyy., country, state_province,
                       locality,Kindgom,Phylum, Class,
                       Order,Family,Genus,Species,dist_coast_km) %>%
  summarise(reads = sum(reads))


write.csv(df1, "~downloads/Datasets/merged_diat_data.csv")



#--------Fig. 2b,c --- Foram NMDS - lbf ESVs > 99.4------------------------
df <- read.csv("~downloads/Datasets/merged_forams_data.csv")
tsm_chl <- read.csv("~downloads/Datasets/filtered_tsm_chl_mean_islands.csv")

df$depth_water <- df$depth..below.the.local.surface..in.meters.

df <- merge(df, tsm_chl[,c(1,4,5)], by = "island", all = TRUE)

#rarefy to the sample (none negative) with the lowest read number: 1924)

dfx <- df %>% 
  group_by(field.nmbr.) %>% 
  mutate(N=reads/sum(reads)*1924)

dfx <- dfx %>% 
  group_by(field.nmbr.) %>% 
  mutate(percent_N=reads/sum(reads)*100)

dfx1 <- dfx[dfx$group == "LBF" & dfx$identity > 99.4, ]

df2 <- dfx1 %>% 
  group_by(field.nmbr., ESVs, species,dist_coast_km, island, mean_chl_filtered,mean_tsm_filtered, depth_water, substrate, group, family) %>% 
  summarise(N_mean=mean(N), percent_N_mean = mean(percent_N))

df2$reef <- NA
df2$reef[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- "flat"
df2$reef[df2$depth_water %not% c("RF1", "RF2", "RF3")] <- "slope"
df2$reef[df2$depth_water %in% c("negative")] <- NA
df2$depth_water[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- 1

df2$depth_water <- as.numeric(paste(df2$depth_water))

df2$turbidity <- NA
df2$turbidity[df2$mean_tsm_filtered > 1.5] <- "high turbidity"
df2$turbidity[df2$mean_tsm_filtered < 0.5] <- "low turbidity"
df2$turbidity[is.na(df2$turbidity)] <- "moderate turbidity"

df2$substrate[df2$substrate == "rubble + seagrass + algae"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "rubble + mud"] <- "rubble + sand"
df2$substrate[df2$substrate == "rubble + algae + seagrass"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "thick layer algae on rubble"] <- "rubble + algae"
df2$substrate[df2$substrate == "seagrass + algae + rubble + sand"] <- "rubble + seagrass"

df2$morpho_species <- NA

df2$morpho_species[df2$species == "Alveolinella_quoyi_Spermonde"] <- "Alveolinella quoyi"
df2$morpho_species[df2$species == "Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS_Spermonde"] <- "Amphisorus sp."
df2$morpho_species[df2$species == "Amphistegina_lessonii"] <- "Amphistegina lessonii"
df2$morpho_species[df2$species == "Amphistegina_lessonii_5174"] <- "Amphistegina lessonii"
df2$morpho_species[df2$species == "Amphistegina_lobifera"] <- "Amphistegina lobifera"
df2$morpho_species[df2$species == "Amphistegina_papillosa_Rik_9440_consensus0_1051_SaudiArabia"] <- "Amphistegina papillosa"
df2$morpho_species[df2$species == "Amphistegina_papillosa_Spermonde_3741"] <- "Amphistegina papillosa"
df2$morpho_species[df2$species == "Amphistegina_radiata_Spermonde"] <- "Amphistegina radiata"
df2$morpho_species[df2$species == "Baculogypsinoides_spinosus_8078_8079_8080"] <- "Baculogypsinoides spinosus"
df2$morpho_species[df2$species == "Borelis_schlumbergeri_3737_Spermonde"] <- "Borelis schlumbergeri"
df2$morpho_species[df2$species == "Calcarina_hispida_&_Calcarina_sp_5163_Spermonde"] <- "Calcarina hispida"
df2$morpho_species[df2$species == "Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._5247_Spermonde"] <- "Calcarina spengleri"
df2$morpho_species[df2$species == "Calcarina_sp._5164_Spermonde"] <- "Calcarina sp."
df2$morpho_species[df2$species == "Calcarina_sp._5165_Spermonde"] <- "Calcarina sp."                           
df2$morpho_species[df2$species == "Calcarina_sp._Spermonde"] <- "Calcarina sp."                  
df2$morpho_species[df2$species == "Calcarina_spengleri_Rik_9352_consensus0_1732_Spermonde"] <- "Calcarina spengleri"   
df2$morpho_species[df2$species == "Cyrea_szymborska_17247"] <- "Z_other_foraminifera"   
df2$morpho_species[df2$species == "Eggerelloides_scaber_12301"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Elphidium_Rik_9327_9337_9350_Spermonde"] <- "Elphidium sp."
df2$morpho_species[df2$species == "Encrusting_Rik_9441_consensus0_977"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Glabratellina_sp._17921_&_Planoglabratella_opercularis_18053_&_Planoglabratella_sp"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Globocassidulina_biora_17225"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Haynesina_germanica_18209"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Heterostegina_depressa_sp2_Spermonde"] <- "Heterostegina depressa"
df2$morpho_species[df2$species == "Heterostegina_depressa_sp1_Spermonde"] <- "Heterostegina depressa"
df2$morpho_species[df2$species == "Miliolid_Rik_9400_9405_Israel"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Milliolid_red_Rik_9370_9371_9372_Israel"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Murrayinella_globosa_18043"] <- "Z_other_foraminifera"                                                         
df2$morpho_species[df2$species == "Neorotalia_calcar_Spermonde"] <- "Neorotalia calcar"
df2$morpho_species[df2$species == "Neorotalia_gaimardi_&_Baculogypsina_sphaerulata_Spermonde"] <- "Neorotalia gaimardi"
df2$morpho_species[df2$species == "Nummulites_venosus"] <- "Nummulites venosus"
df2$morpho_species[df2$species == "Operculina_ammonoides"] <- "Operculina ammonoides"
df2$morpho_species[df2$species == "Parasorites_sp"] <- "Parasorites sp."
df2$morpho_species[df2$species == "Peneroplis_sp1_Spermonde"] <- "Peneroplis sp."
df2$morpho_species[df2$species == "Peneroplis_sp2_&_Peneroplis_pertusus_5117_&_Dendritina_ambigua_Spermonde"] <- "Peneroplis sp."
df2$morpho_species[df2$species == "Planorbulinella_sp._17893"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Psammophaga_sp._19260_&_Psammophaga_sp._19296_&_Psammophaga_sp2"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Pyrgo_Rik_9347_consensus0_1557"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Sorites_sp1_Spermonde"] <- "Sorites sp."
df2$morpho_species[df2$species == "Sorites_sp2_Spermonde"] <- "Sorites sp."
df2$morpho_species[df2$species == "Stainforthia_sp._17321"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Textularid_Rik_9397_consensus4_105"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Textularid_Rik_9398_consensus128_73"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Trifarina_earlandi_&_Uvigerina_bifurcata"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Vertebralina_Rik_9325_consensus0_1751"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Z_not_foraminifera"] <- "Z_not_foraminifera" 
df2$morpho_species[df2$species == "Z_other_foraminifera"] <- "Z_other_foraminifera" 

unique(df2$morpho_species)

df2a <- df2[df2$island != "negative",]
df2a <- df2a[!is.na(df2a$island),]

foram_tot <- unique(df$ESVs)
foram_order <- unique(df2a$ESVs)
#foram_group <- unique(df2b$ESVs)

sum(df2a$N_mean, na.rm = TRUE)/sum(dfx$N,na.rm = TRUE)

write.csv(df2a,"~downloads/Datasets/foram_working_dataset_turnover.csv", row.names=FALSE)


#____reef flat_______
df2b <- df2a[df2a$reef == "flat",]

df3 <- dcast(df2b, field.nmbr. + depth_water + island + dist_coast_km + mean_chl_filtered + mean_tsm_filtered + turbidity + reef + substrate ~ ESVs, 
             fun.aggregate = sum, value.var = "N_mean")

df3[df3 == "NaN"] <- 0

df_NMDS_vegan <- dcast(df2b, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
ord <- metaMDS(df_NMDS_vegan, distance = "bray")
env <- df3[,10:ncol(df3)]

#to calculate the species vectors going over the NMDS
en <- envfit(ord, env, permutations = 999, na.rm = TRUE)

#groups for plots
island <- droplevels(as.factor(df3$island))
depth <- as.factor(round(as.numeric(df3$depth_water, digits = 0)))
field <- droplevels(as.factor(df3$field.nmbr.))
dist <- as.numeric(df3$dist_coast_km)
tsm <- as.numeric(df3$mean_tsm_filtered)
chl <- as.numeric(df3$mean_chl_filtered)
subst <- droplevels(as.factor(df3$substrate))
reef <- droplevels(as.factor(df3$reef))
turb <- droplevels(as.factor(df3$turbidity))

colfield <- hcl.colors(12, "Roma")

islandord <- c("Lae Lae","Barang Baringan","Samalona","Polewali","Karanrang",
               "Langkadea","K. Keke","Padjenekang","Badi","Lumulumu","Langkai",
               "Kapoposang Masdar point")

turbidity <- c("high turbidity", "moderate turbidity", "low turbidity")
colturb <- c("#9E4E09" ,"#C7AE4C" ,"#008FA7" )
pchils <- c(0:6,8,15,17:19)

#Significance is the p value
ano <- anosim(df_NMDS_vegan, turb)

#for plot _ distance to shore
plot(ord, disp="sites", type = "n",main = paste ("Reef flat: Similarity between turbidity levels - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3)))
ordisurf(ord, dist, col = "black",linetype = "dashed", add = TRUE)
ordihull(ord, groups=turb,  draw = "polygon", col = colturb[c(1,3,2)])
points(ord, disp="sites", pch = pchils[factor(island, islandord)] , cex=1, col = "black")
with(ord, legend(x = "topleft", legend = levels(factor(island, islandord)), col = "black", pch = pchils))
with(ord, legend(x = "bottomleft", legend = levels(factor(turb, turbidity)), col = colturb, pch = 15))
#orditorp(ord,display="sites",cex=0.7,air=0.25)


#____reef slope_______
df2b <- df2a[df2a$reef == "slope",]

df3 <- dcast(df2b, field.nmbr. + depth_water + island + turbidity + dist_coast_km + reef + substrate ~ ESVs, 
             fun.aggregate = sum, value.var = "N_mean")

df3[df3 == "NaN"] <- 0

df_NMDS_vegan <- dcast(df2b, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
ord <- metaMDS(df_NMDS_vegan, distance = "bray")
env <- df3[,8:ncol(df3)]

#to calculate the species vectors going over the NMDS
en <- envfit(ord, env, permutations = 999, na.rm = TRUE)

#groups for plots
island <- droplevels(as.factor(df3$island))
depth <- as.numeric(df3$depth_water)
field <- droplevels(as.factor(df3$field.nmbr.))
dist <- as.numeric(df3$dist_coast_km)
subst <- droplevels(as.factor(df3$substrate))
reef <- droplevels(as.factor(df3$reef))
turb <- droplevels(as.factor(df3$turbidity))

colfield <- hcl.colors(3, "Roma")
shapes <- c(1,4,13,17,16,6)
substord <- c("sand","rubble", "rubble + sand", "rubble + algae",
              "rubble + sand + algae", "rubble + seagrass")
islandord <- c("Lae Lae","Barang Baringan","Samalona","Polewali","Karanrang",
               "Langkadea","K. Keke","Padjenekang","Badi","Lumulumu","Langkai",
               "Kapoposang Masdar point")

turbidity <- c("high turbidity", "moderate turbidity", "low turbidity")
colturb <- c("#9E4E09" ,"#C7AE4C" ,"#008FA7" )

#for plot _ distance to shore
pchils <- c(0:6,8,15,17:19)

#Significance is the p value
ano <- anosim(df_NMDS_vegan, turb)

#for plot _ distance to shore
plot(ord, disp="sites", type = "n",main = paste ("Reef slope: Similarity between turbidity levels - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3)))
ordisurf(ord, dist, col = "black",linetype = "dashed", add = TRUE)
ordihull(ord, groups=turb,  draw = "polygon", col = colturb[c(1,3,2)])
points(ord, disp="sites", pch = pchils[factor(island, islandord)] , cex=1, col = "black")
with(ord, legend(x = "topleft", legend = levels(factor(island, islandord)), col = "black", pch = pchils))
with(ord, legend(x = "bottomleft", legend = levels(factor(turb, turbidity)), col = colturb, pch = 15))




#--------Fig. 2d,e --- Foram NMDS - all foram ESVs------------------------
df <- read.csv("~downloads/Datasets/merged_forams_data.csv")
tsm_chl <- read.csv("~downloads/Datasets/filtered_tsm_chl_mean_islands.csv")

df$depth_water <- df$depth..below.the.local.surface..in.meters.

df <- merge(df, tsm_chl[,c(1,4,5)], by = "island", all = TRUE)

#rarefy to the sample (none negative) with the lowest read number: 1924)

dfx <- df %>% 
  group_by(field.nmbr.) %>% 
  mutate(N=reads/sum(reads)*1924)

dfx <- dfx %>% 
  group_by(field.nmbr.) %>% 
  mutate(percent_N=reads/sum(reads)*100)

dfx1 <- dfx

df2 <- dfx1 %>% 
  group_by(field.nmbr., ESVs, species,dist_coast_km, island, mean_chl_filtered,mean_tsm_filtered, depth_water, substrate, group, family) %>% 
  summarise(N_mean=mean(N), percent_N_mean = mean(percent_N))

df2$reef <- NA
df2$reef[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- "flat"
df2$reef[df2$depth_water %not% c("RF1", "RF2", "RF3")] <- "slope"
df2$reef[df2$depth_water %in% c("negative")] <- NA
df2$depth_water[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- 1

df2$depth_water <- as.numeric(paste(df2$depth_water))

df2$turbidity <- NA
df2$turbidity[df2$mean_tsm_filtered > 1.5] <- "high turbidity"
df2$turbidity[df2$mean_tsm_filtered < 0.5] <- "low turbidity"
df2$turbidity[is.na(df2$turbidity)] <- "moderate turbidity"

df2$substrate[df2$substrate == "rubble + seagrass + algae"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "rubble + mud"] <- "rubble + sand"
df2$substrate[df2$substrate == "rubble + algae + seagrass"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "thick layer algae on rubble"] <- "rubble + algae"
df2$substrate[df2$substrate == "seagrass + algae + rubble + sand"] <- "rubble + seagrass"

df2$morpho_species <- NA

df2$morpho_species[df2$species == "Alveolinella_quoyi_Spermonde"] <- "Alveolinella quoyi"
df2$morpho_species[df2$species == "Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS_Spermonde"] <- "Amphisorus sp."
df2$morpho_species[df2$species == "Amphistegina_lessonii"] <- "Amphistegina lessonii"
df2$morpho_species[df2$species == "Amphistegina_lessonii_5174"] <- "Amphistegina lessonii"
df2$morpho_species[df2$species == "Amphistegina_lobifera"] <- "Amphistegina lobifera"
df2$morpho_species[df2$species == "Amphistegina_papillosa_Rik_9440_consensus0_1051_SaudiArabia"] <- "Amphistegina papillosa"
df2$morpho_species[df2$species == "Amphistegina_papillosa_Spermonde_3741"] <- "Amphistegina papillosa"
df2$morpho_species[df2$species == "Amphistegina_radiata_Spermonde"] <- "Amphistegina radiata"
df2$morpho_species[df2$species == "Baculogypsinoides_spinosus_8078_8079_8080"] <- "Baculogypsinoides spinosus"
df2$morpho_species[df2$species == "Borelis_schlumbergeri_3737_Spermonde"] <- "Borelis schlumbergeri"
df2$morpho_species[df2$species == "Calcarina_hispida_&_Calcarina_sp_5163_Spermonde"] <- "Calcarina hispida"
df2$morpho_species[df2$species == "Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._5247_Spermonde"] <- "Calcarina spengleri"
df2$morpho_species[df2$species == "Calcarina_sp._5164_Spermonde"] <- "Calcarina sp."
df2$morpho_species[df2$species == "Calcarina_sp._5165_Spermonde"] <- "Calcarina sp."                           
df2$morpho_species[df2$species == "Calcarina_sp._Spermonde"] <- "Calcarina sp."                  
df2$morpho_species[df2$species == "Calcarina_spengleri_Rik_9352_consensus0_1732_Spermonde"] <- "Calcarina spengleri"   
df2$morpho_species[df2$species == "Cyrea_szymborska_17247"] <- "Z_other_foraminifera"   
df2$morpho_species[df2$species == "Eggerelloides_scaber_12301"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Elphidium_Rik_9327_9337_9350_Spermonde"] <- "Elphidium sp."
df2$morpho_species[df2$species == "Encrusting_Rik_9441_consensus0_977"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Glabratellina_sp._17921_&_Planoglabratella_opercularis_18053_&_Planoglabratella_sp"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Globocassidulina_biora_17225"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Haynesina_germanica_18209"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Heterostegina_depressa_sp2_Spermonde"] <- "Heterostegina depressa"
df2$morpho_species[df2$species == "Heterostegina_depressa_sp1_Spermonde"] <- "Heterostegina depressa"
df2$morpho_species[df2$species == "Miliolid_Rik_9400_9405_Israel"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Milliolid_red_Rik_9370_9371_9372_Israel"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Murrayinella_globosa_18043"] <- "Z_other_foraminifera"                                                         
df2$morpho_species[df2$species == "Neorotalia_calcar_Spermonde"] <- "Neorotalia calcar"
df2$morpho_species[df2$species == "Neorotalia_gaimardi_&_Baculogypsina_sphaerulata_Spermonde"] <- "Neorotalia gaimardi"
df2$morpho_species[df2$species == "Nummulites_venosus"] <- "Nummulites venosus"
df2$morpho_species[df2$species == "Operculina_ammonoides"] <- "Operculina ammonoides"
df2$morpho_species[df2$species == "Parasorites_sp"] <- "Parasorites sp."
df2$morpho_species[df2$species == "Peneroplis_sp1_Spermonde"] <- "Peneroplis sp."
df2$morpho_species[df2$species == "Peneroplis_sp2_&_Peneroplis_pertusus_5117_&_Dendritina_ambigua_Spermonde"] <- "Peneroplis sp."
df2$morpho_species[df2$species == "Planorbulinella_sp._17893"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Psammophaga_sp._19260_&_Psammophaga_sp._19296_&_Psammophaga_sp2"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Pyrgo_Rik_9347_consensus0_1557"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Sorites_sp1_Spermonde"] <- "Sorites sp."
df2$morpho_species[df2$species == "Sorites_sp2_Spermonde"] <- "Sorites sp."
df2$morpho_species[df2$species == "Stainforthia_sp._17321"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Textularid_Rik_9397_consensus4_105"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Textularid_Rik_9398_consensus128_73"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Trifarina_earlandi_&_Uvigerina_bifurcata"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Vertebralina_Rik_9325_consensus0_1751"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Z_not_foraminifera"] <- "Z_not_foraminifera" 
df2$morpho_species[df2$species == "Z_other_foraminifera"] <- "Z_other_foraminifera" 

unique(df2$morpho_species)

df2a <- df2[df2$island != "negative",]
df2a <- df2a[!is.na(df2a$island),]

foram_tot <- unique(df$ESVs)
foram_order <- unique(df2a$ESVs)
#foram_group <- unique(df2b$ESVs)

sum(df2a$N_mean, na.rm = TRUE)/sum(dfx$N,na.rm = TRUE)

write.csv(df2a,"~downloads/Datasets/Allforam_working_dataset_turnover.csv", row.names=FALSE)


#____reef flat_______
df2b <- df2a[df2a$reef == "flat",]

df3 <- dcast(df2b, field.nmbr. + depth_water + island + dist_coast_km + mean_chl_filtered + mean_tsm_filtered + turbidity + reef + substrate ~ ESVs, 
             fun.aggregate = sum, value.var = "N_mean")

df3[df3 == "NaN"] <- 0

df_NMDS_vegan <- dcast(df2b, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
ord <- metaMDS(df_NMDS_vegan, distance = "bray")
env <- df3[,10:ncol(df3)]

#to calculate the species vectors going over the NMDS
en <- envfit(ord, env, permutations = 999, na.rm = TRUE)

#groups for plots
island <- droplevels(as.factor(df3$island))
depth <- as.factor(round(as.numeric(df3$depth_water, digits = 0)))
field <- droplevels(as.factor(df3$field.nmbr.))
dist <- as.numeric(df3$dist_coast_km)
tsm <- as.numeric(df3$mean_tsm_filtered)
chl <- as.numeric(df3$mean_chl_filtered)
subst <- droplevels(as.factor(df3$substrate))
reef <- droplevels(as.factor(df3$reef))
turb <- droplevels(as.factor(df3$turbidity))

colfield <- hcl.colors(12, "Roma")

islandord <- c("Lae Lae","Barang Baringan","Samalona","Polewali","Karanrang",
               "Langkadea","K. Keke","Padjenekang","Badi","Lumulumu","Langkai",
               "Kapoposang Masdar point")

turbidity <- c("high turbidity", "moderate turbidity", "low turbidity")
colturb <- c("#9E4E09" ,"#C7AE4C" ,"#008FA7" )
pchils <- c(0:6,8,15,17:19)

#Significance is the p value
ano <- anosim(df_NMDS_vegan, turb)

#for plot _ distance to shore
plot(ord, disp="sites", type = "n",main = paste ("Reef flat: Similarity between turbidity levels - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3)))
ordisurf(ord, dist, col = "black",linetype = "dashed", add = TRUE)
ordihull(ord, groups=turb,  draw = "polygon", col = colturb[c(1,3,2)])
points(ord, disp="sites", pch = pchils[factor(island, islandord)] , cex=1, col = "black")
with(ord, legend(x = "topleft", legend = levels(factor(island, islandord)), col = "black", pch = pchils))
with(ord, legend(x = "bottomleft", legend = levels(factor(turb, turbidity)), col = colturb, pch = 15))
#orditorp(ord,display="sites",cex=0.7,air=0.25)


#____reef slope_______
df2b <- df2a[df2a$reef == "slope",]

df3 <- dcast(df2b, field.nmbr. + depth_water + island + turbidity + dist_coast_km + reef + substrate ~ ESVs, 
             fun.aggregate = sum, value.var = "N_mean")

df3[df3 == "NaN"] <- 0

df_NMDS_vegan <- dcast(df2b, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
ord <- metaMDS(df_NMDS_vegan, distance = "bray")
env <- df3[,8:ncol(df3)]

#to calculate the species vectors going over the NMDS
en <- envfit(ord, env, permutations = 999, na.rm = TRUE)

#groups for plots
island <- droplevels(as.factor(df3$island))
depth <- as.numeric(df3$depth_water)
field <- droplevels(as.factor(df3$field.nmbr.))
dist <- as.numeric(df3$dist_coast_km)
subst <- droplevels(as.factor(df3$substrate))
reef <- droplevels(as.factor(df3$reef))
turb <- droplevels(as.factor(df3$turbidity))

colfield <- hcl.colors(3, "Roma")
shapes <- c(1,4,13,17,16,6)
substord <- c("sand","rubble", "rubble + sand", "rubble + algae",
              "rubble + sand + algae", "rubble + seagrass")
islandord <- c("Lae Lae","Barang Baringan","Samalona","Polewali","Karanrang",
               "Langkadea","K. Keke","Padjenekang","Badi","Lumulumu","Langkai",
               "Kapoposang Masdar point")

turbidity <- c("high turbidity", "moderate turbidity", "low turbidity")
colturb <- c("#9E4E09" ,"#C7AE4C" ,"#008FA7" )

#for plot _ distance to shore
pchils <- c(0:6,8,15,17:19)

#Significance is the p value
ano <- anosim(df_NMDS_vegan, turb)

#for plot _ distance to shore
plot(ord, disp="sites", type = "n",main = paste ("Reef slope: Similarity between turbidity levels - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3)))
ordisurf(ord, dist, col = "black",linetype = "dashed", add = TRUE)
ordihull(ord, groups=turb,  draw = "polygon", col = colturb[c(1,3,2)])
points(ord, disp="sites", pch = pchils[factor(island, islandord)] , cex=1, col = "black")
with(ord, legend(x = "topleft", legend = levels(factor(island, islandord)), col = "black", pch = pchils))
with(ord, legend(x = "bottomleft", legend = levels(factor(turb, turbidity)), col = colturb, pch = 15))

#--------FORAMS - community composition------------------------
df <- read.csv("~downloads/Datasets/merged_forams_data.csv")

df$depth_water <- df$depth..below.the.local.surface..in.meters.

#rarefy to the sample (none negative) with the lowest read number: 1924)
dfx <- df %>% 
  group_by(field.nmbr.) %>% 
  mutate(N=reads/sum(reads)*1924)

dfx <- dfx %>% 
  group_by(field.nmbr.,) %>% 
  mutate(percent_N=reads/sum(reads)*100)

df2 <- dfx %>% 
  group_by(field.nmbr., ESVs, species,dist_coast_km, island, depth_water, substrate, group, family) %>% 
  summarise(N_mean=mean(N), percent_N_mean = mean(percent_N))

df2$counts <- 1

df2$morpho_species <- NA

df2$morpho_species[df2$species == "Alveolinella_quoyi_Spermonde"] <- "Alveolinella quoyi"
df2$morpho_species[df2$species == "Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS_Spermonde"] <- "Amphisorus sp."
df2$morpho_species[df2$species == "Amphistegina_lessonii"] <- "Amphistegina lessonii"
df2$morpho_species[df2$species == "Amphistegina_lessonii_5174"] <- "Amphistegina lessonii"
df2$morpho_species[df2$species == "Amphistegina_lobifera"] <- "Amphistegina lobifera"
df2$morpho_species[df2$species == "Amphistegina_papillosa_Rik_9440_consensus0_1051_SaudiArabia"] <- "Amphistegina papillosa"
df2$morpho_species[df2$species == "Amphistegina_papillosa_Spermonde_3741"] <- "Amphistegina papillosa"
df2$morpho_species[df2$species == "Amphistegina_radiata_Spermonde"] <- "Amphistegina radiata"
df2$morpho_species[df2$species == "Baculogypsinoides_spinosus_8078_8079_8080"] <- "Baculogypsinoides spinosus"
df2$morpho_species[df2$species == "Borelis_schlumbergeri_3737_Spermonde"] <- "Borelis schlumbergeri"
df2$morpho_species[df2$species == "Calcarina_hispida_&_Calcarina_sp_5163_Spermonde"] <- "Calcarina hispida"
df2$morpho_species[df2$species == "Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._5247_Spermonde"] <- "Calcarina spengleri"
df2$morpho_species[df2$species == "Calcarina_sp._5164_Spermonde"] <- "Calcarina sp."
df2$morpho_species[df2$species == "Calcarina_sp._5165_Spermonde"] <- "Calcarina sp."                           
df2$morpho_species[df2$species == "Calcarina_sp._Spermonde"] <- "Calcarina sp."                  
df2$morpho_species[df2$species == "Calcarina_spengleri_Rik_9352_consensus0_1732_Spermonde"] <- "Calcarina spengleri"   
df2$morpho_species[df2$species == "Cyrea_szymborska_17247"] <- "Z_other_foraminifera"   
df2$morpho_species[df2$species == "Eggerelloides_scaber_12301"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Elphidium_Rik_9327_9337_9350_Spermonde"] <- "Elphidium sp."
df2$morpho_species[df2$species == "Encrusting_Rik_9441_consensus0_977"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Glabratellina_sp._17921_&_Planoglabratella_opercularis_18053_&_Planoglabratella_sp"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Globocassidulina_biora_17225"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Haynesina_germanica_18209"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Heterostegina_depressa_sp2_Spermonde"] <- "Heterostegina depressa"
df2$morpho_species[df2$species == "Heterostegina_depressa_sp1_Spermonde"] <- "Heterostegina depressa"
df2$morpho_species[df2$species == "Miliolid_Rik_9400_9405_Israel"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Milliolid_red_Rik_9370_9371_9372_Israel"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Murrayinella_globosa_18043"] <- "Z_other_foraminifera"                                                         
df2$morpho_species[df2$species == "Neorotalia_calcar_Spermonde"] <- "Neorotalia calcar"
df2$morpho_species[df2$species == "Neorotalia_gaimardi_&_Baculogypsina_sphaerulata_Spermonde"] <- "Neorotalia gaimardi"
df2$morpho_species[df2$species == "Nummulites_venosus"] <- "Nummulites venosus"
df2$morpho_species[df2$species == "Operculina_ammonoides"] <- "Operculina ammonoides"
df2$morpho_species[df2$species == "Parasorites_sp"] <- "Parasorites sp."
df2$morpho_species[df2$species == "Peneroplis_sp1_Spermonde"] <- "Peneroplis sp."
df2$morpho_species[df2$species == "Peneroplis_sp2_&_Peneroplis_pertusus_5117_&_Dendritina_ambigua_Spermonde"] <- "Peneroplis sp."
df2$morpho_species[df2$species == "Planorbulinella_sp._17893"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Psammophaga_sp._19260_&_Psammophaga_sp._19296_&_Psammophaga_sp2"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Pyrgo_Rik_9347_consensus0_1557"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Sorites_sp1_Spermonde"] <- "Sorites sp."
df2$morpho_species[df2$species == "Sorites_sp2_Spermonde"] <- "Sorites sp."
df2$morpho_species[df2$species == "Stainforthia_sp._17321"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Textularid_Rik_9397_consensus4_105"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Textularid_Rik_9398_consensus128_73"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Trifarina_earlandi_&_Uvigerina_bifurcata"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Vertebralina_Rik_9325_consensus0_1751"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Z_not_foraminifera"] <- "Z_not_foraminifera" 
df2$morpho_species[df2$species == "Z_other_foraminifera"] <- "Z_other_foraminifera" 

#number of molecular species
unique(df2$species[df2$group == "LBF"])

#number of morpho species
unique(df2$morpho_species[df2$group == "LBF"])

#number of families
unique(df2$family[df2$group == "LBF"])


#--------Fig. 3b,c,d,e --- Diatom NMDS - all diatom ESVs------------------------
df <- read.csv("~downloads/Datasets/merged_diat_data.csv")

df$depth_water <- df$depth..below.the.local.surface..in.meters.
  
tsm_chl <- read.csv("~downloads/Datasets/filtered_tsm_chl_mean_islands.csv")

df <- merge(df, tsm_chl[,c(1,4,5)], by = "island", all = TRUE)


#rarefy to the sample (none negative) with the lowest read number)
test <- df %>% 
  group_by(field.nmbr.,Substrate.type) %>% 
  summarise(N_sum=sum(reads))
test1 <- test[test$N_sum > 1000,] #lowest read number for sample > 1000 = 1654
df1 <- merge(df, test1[,1:3], by = c("field.nmbr.", "Substrate.type"))

dfx <- df1 %>% 
  group_by(field.nmbr.) %>% 
  mutate(N=reads/sum(reads)*1654)

dfx <- dfx[dfx$Phylum == "Bacillariophyta",] #assigned only to diatoms

dfx1 <- dfx


df2 <- dfx1 %>% 
  group_by(field.nmbr., ESVs, Species,dist_coast_km, island, mean_chl_filtered,mean_tsm_filtered, depth_water, substrate, Substrate.type) %>% 
  summarise(N_mean=mean(N))


df2$reef <- NA
df2$reef[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- "flat"
df2$reef[df2$depth_water %not% c("RF1", "RF2", "RF3")] <- "slope"
df2$reef[is.na(df2$depth_water)] <- NA
df2$depth_water[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- 1
df2$depth_water <- as.numeric(df2$depth_water)
df2$depth_group <- NA
df2$depth_group[df2$depth_water < 10 & df2$reef == "slope"] <- "<10m"
df2$depth_group[df2$depth_water > 9 & df2$depth_water < 21] <- "10-20m"
df2$depth_group[df2$depth_water > 20 & df2$reef == "slope"] <- ">20m"

df2$turbidity <- NA
df2$turbidity[df2$mean_tsm_filtered > 1.5] <- "high turbidity"
df2$turbidity[df2$mean_tsm_filtered < 0.5] <- "low turbidity"
df2$turbidity[is.na(df2$turbidity)] <- "moderate turbidity"

df2$substrate[df2$substrate == "rubble + seagrass + algae"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "rubble + mud"] <- "rubble + sand"
df2$substrate[df2$substrate == "rubble + algae + seagrass"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "thick layer algae on rubble"] <- "rubble + algae"
df2$substrate[df2$substrate == "seagrass + algae + rubble + sand"] <- "rubble + seagrass"

df2a <- df2[df2$island != "negative",]
df2a <- df2a[!is.na(df2a$island),]

diat_tot <- unique(dfx$ESVs)
diat_genus <- unique(df2a$ESVs)
#diat_group <- unique(df2b$ESVs)

sum(df2a$N_mean, na.rm = TRUE)/sum(dfx$N,na.rm = TRUE)

write.csv(df2a,"~downloads/Datasets/diat_working_dataset_turnover.csv", row.names=FALSE)

#____reef flat - seawater_______
df2b <- df2a[df2a$reef == "flat" & df2a$Substrate.type == "seawater",]

#remove 1 outlier: W030 (Langkai)
df2b <- df2b[df2b$field.nmbr. != "UPG-EG22-W030",]

df3 <- dcast(df2b, field.nmbr. + depth_group + island + dist_coast_km + reef + substrate + Substrate.type + turbidity ~ ESVs, 
             fun.aggregate = sum, value.var = "N_mean")

df3[df3 == "NaN"] <- 0

df_NMDS_vegan <- dcast(df2b, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
ord <- metaMDS(df_NMDS_vegan, distance = "bray")
env <- df3[,9:ncol(df3)]

#to calculate the species vectors going over the NMDS
en <- envfit(ord, env, permutations = 999, na.rm = TRUE)

#groups for plots
island <- droplevels(as.factor(df3$island))
depth <- as.factor(df3$depth_group)
field <- droplevels(as.factor(df3$field.nmbr.))
dist <- as.numeric(df3$dist_coast_km)
subst <- droplevels(as.factor(df3$substrate))
reef <- droplevels(as.factor(df3$reef))
type <- droplevels(as.factor(df3$Substrate.type))
turb <- droplevels(as.factor(df3$turbidity))

colfield <- hcl.colors(12, "Roma")
shapes <- c(19,17)

substord <- c("sand","rubble", "rubble + sand", "rubble + algae",
              "rubble + sand + algae", "rubble + seagrass")
islandord <- c("Lae Lae","Barang Baringan","Samalona","Polewali","Karanrang",
               "Langkadea","K. Keke","Padjenekang","Badi","Lumulumu","Langkai",
               "Kapoposang Masdar point")

turbidity <- c("high turbidity", "moderate turbidity", "low turbidity")
colturb <- c("#9E4E09" ,"#C7AE4C" ,"#008FA7" )
pchils <- c(0:6,8,15,17:19)

#Significance is the p value
ano <- anosim(df_NMDS_vegan, turb)

#for plot _ distance to shore
plot(ord, disp="sites", type = "n",main = paste ("diat Reef flat SW: turbidity levels - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3)))
ordisurf(ord, dist, col = "black",linetype = "dashed", add = TRUE)
ordihull(ord, groups=turb,  draw = "polygon", col = colturb[c(1,3,2)])
points(ord, disp="sites", pch = pchils[factor(island, islandord)] , cex=1, col = "black")
with(ord, legend(x = "topleft", legend = levels(factor(island, islandord)), col = "black", pch = pchils))
with(ord, legend(x = "bottomleft", legend = levels(factor(turb, turbidity)), col = colturb, pch = 15))


#____reef slope - seawater_______
df2b <- df2a[df2a$reef == "slope" & df2a$Substrate.type == "seawater",]

#remove 1 outlier: W027 (Langkai)
df2b <- df2b[df2b$field.nmbr. != "UPG-EG22-W027",]

df3 <- dcast(df2b, field.nmbr. + depth_group + island + dist_coast_km + reef + substrate + Substrate.type + turbidity ~ ESVs, 
             fun.aggregate = sum, value.var = "N_mean")

df3[df3 == "NaN"] <- 0

df_NMDS_vegan <- dcast(df2b, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
ord <- metaMDS(df_NMDS_vegan, distance = "bray")
env <- df3[,9:ncol(df3)]

#to calculate the species vectors going over the NMDS
en <- envfit(ord, env, permutations = 999, na.rm = TRUE)

#groups for plots
island <- droplevels(as.factor(df3$island))
depth <- as.factor(df3$depth_group)
field <- droplevels(as.factor(df3$field.nmbr.))
dist <- as.numeric(df3$dist_coast_km)
subst <- droplevels(as.factor(df3$substrate))
reef <- droplevels(as.factor(df3$reef))
type <- droplevels(as.factor(df3$Substrate.type))
turb <- droplevels(as.factor(df3$turbidity))

colfield <- hcl.colors(12, "Roma")
shapes <- c(19,17)

substord <- c("sand","rubble", "rubble + sand", "rubble + algae",
              "rubble + sand + algae", "rubble + seagrass")
islandord <- c("Lae Lae","Barang Baringan","Samalona","Polewali","Karanrang",
               "Langkadea","K. Keke","Padjenekang","Badi","Lumulumu","Langkai",
               "Kapoposang Masdar point")

turbidity <- c("high turbidity", "moderate turbidity", "low turbidity")
colturb <- c("#9E4E09" ,"#C7AE4C" ,"#008FA7" )
pchils <- c(0:6,8,15,17:19)

#Significance is the p value
ano <- anosim(df_NMDS_vegan, turb)

#for plot _ distance to shore
plot(ord, disp="sites", type = "n",main = paste ("diat Reef slope sw: turbidity levels - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3)))
ordisurf(ord, dist, col = "black",linetype = "dashed", add = TRUE)
ordihull(ord, groups=turb,  draw = "polygon", col = colturb[c(1,3,2)])
points(ord, disp="sites", pch = pchils[factor(island, islandord)] , cex=1, col = "black")
with(ord, legend(x = "topleft", legend = levels(factor(island, islandord)), col = "black", pch = pchils))
with(ord, legend(x = "bottomleft", legend = levels(factor(turb, turbidity)), col = colturb, pch = 15))


#____reef flat - substrate_______
df2b <- df2a[df2a$reef == "flat" & df2a$Substrate.type == "substrate",]

df3 <- dcast(df2b, field.nmbr. + depth_group + island + dist_coast_km + reef + turbidity + substrate + Substrate.type ~ ESVs, 
             fun.aggregate = sum, value.var = "N_mean")

df3[df3 == "NaN"] <- 0

df_NMDS_vegan <- dcast(df2b, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
ord <- metaMDS(df_NMDS_vegan, distance = "bray")
env <- df3[,9:ncol(df3)]

#to calculate the species vectors going over the NMDS
en <- envfit(ord, env, permutations = 999, na.rm = TRUE)

#groups for plots
island <- droplevels(as.factor(df3$island))
depth <- as.factor(df3$depth_group)
field <- droplevels(as.factor(df3$field.nmbr.))
dist <- as.numeric(df3$dist_coast_km)
subst <- droplevels(as.factor(df3$substrate))
reef <- droplevels(as.factor(df3$reef))
type <- droplevels(as.factor(df3$Substrate.type))
turb <- droplevels(as.factor(df3$turbidity))

colfield <- hcl.colors(12, "Roma")
shapes <- c(19,17)

substord <- c("sand","rubble", "rubble + sand", "rubble + algae",
              "rubble + sand + algae", "rubble + seagrass")
islandord <- c("Lae Lae","Barang Baringan","Samalona","Polewali","Karanrang",
               "Langkadea","K. Keke","Padjenekang","Badi","Lumulumu","Langkai",
               "Kapoposang Masdar point")

turbidity <- c("high turbidity", "moderate turbidity", "low turbidity")
colturb <- c("#9E4E09" ,"#C7AE4C" ,"#008FA7" )
pchils <- c(0:6,8,15,17:19)

#Significance is the p value
ano <- anosim(df_NMDS_vegan, turb)

#for plot _ distance to shore
plot(ord, disp="sites", type = "n",main = paste ("diat Reef flat subs:turbidity levels - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3)))
ordisurf(ord, dist, col = "black",linetype = "dashed", add = TRUE)
ordihull(ord, groups=turb,  draw = "polygon", col = colturb[c(1,3,2)])
points(ord, disp="sites", pch = pchils[factor(island, islandord)] , cex=1, col = "black")
with(ord, legend(x = "topleft", legend = levels(factor(island, islandord)), col = "black", pch = pchils))
with(ord, legend(x = "bottomleft", legend = levels(factor(turb, turbidity)), col = colturb, pch = 15))


#____reef slope - substrate_______
df2b <- df2a[df2a$reef == "slope" & df2a$Substrate.type == "substrate",]

df3 <- dcast(df2b, field.nmbr. + depth_group + island + dist_coast_km + reef + turbidity + substrate + Substrate.type ~ ESVs, 
             fun.aggregate = sum, value.var = "N_mean")

df3[df3 == "NaN"] <- 0

df_NMDS_vegan <- dcast(df2b, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
ord <- metaMDS(df_NMDS_vegan, distance = "bray")
env <- df3[,9:ncol(df3)]

#to calculate the species vectors going over the NMDS
en <- envfit(ord, env, permutations = 999, na.rm = TRUE)

#groups for plots
island <- droplevels(as.factor(df3$island))
depth <- as.factor(df3$depth_group)
field <- droplevels(as.factor(df3$field.nmbr.))
dist <- as.numeric(df3$dist_coast_km)
subst <- droplevels(as.factor(df3$substrate))
reef <- droplevels(as.factor(df3$reef))
type <- droplevels(as.factor(df3$Substrate.type))
turb <- droplevels(as.factor(df3$turbidity))

colfield <- hcl.colors(12, "Roma")
shapes <- c(19,17)

substord <- c("sand","rubble", "rubble + sand", "rubble + algae",
              "rubble + sand + algae", "rubble + seagrass")
islandord <- c("Lae Lae","Barang Baringan","Samalona","Polewali","Karanrang",
               "Langkadea","K. Keke","Padjenekang","Badi","Lumulumu","Langkai",
               "Kapoposang Masdar point")

turbidity <- c("high turbidity", "moderate turbidity", "low turbidity")
colturb <- c("#9E4E09" ,"#C7AE4C" ,"#008FA7" )
pchils <- c(0:6,8,15,17:19)

#Significance is the p value
ano <- anosim(df_NMDS_vegan, turb)

#for plot _ distance to shore
plot(ord, disp="sites", type = "n",main = paste ("diat Reef slope subs:turbidity levels - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3)))
ordisurf(ord, dist, col = "black",linetype = "dashed", add = TRUE)
ordihull(ord, groups=turb,  draw = "polygon", col = colturb[c(1,3,2)])
points(ord, disp="sites", pch = pchils[factor(island, islandord)] , cex=1, col = "black")
with(ord, legend(x = "topleft", legend = levels(factor(island, islandord)), col = "black", pch = pchils))
with(ord, legend(x = "bottomleft", legend = levels(factor(turb, turbidity)), col = colturb, pch = 15))



#--------DIATOMS - community composition---------------
df <- read.csv("~downloads/Datasets/merged_diat_data.csv")

df$depth_water <- df$depth..below.the.local.surface..in.meters.
df <- df[df$Phylum == "Bacillariophyta",]

#rarefy to the sample (none negative) with the lowest read number: 1924)
dfx <- df %>% 
  group_by(field.nmbr.) %>% 
  mutate(N=reads/sum(reads)*1654)

dfx <- dfx %>% 
  group_by(field.nmbr.) %>% 
  mutate(percent_N=reads/sum(reads)*100)

df2 <- dfx %>% 
  group_by(field.nmbr., ESVs,dist_coast_km, island, depth_water, substrate, Substrate.type,Identity.percentage, Class, Order, Family, Genus, Species) %>% 
  summarise(N_mean=mean(N), percent_N_mean = mean(percent_N))

df2$counts <- 1

df2$Species[df2$Identity.percentage < 95] <- "unknown species"
df2$Genus[df2$Identity.percentage < 90] <- "unknown genus"
df2$Family[df2$Identity.percentage < 85] <- "unknown family"


sort(unique(df2$Class[df2$Identity.percentage > 90]))
sort(unique(df2$Order[df2$Identity.percentage > 90]))
sort(unique(df2$Family[df2$Identity.percentage > 90]))
sort(unique(df2$Genus[df2$Identity.percentage > 90]))
sort(unique(df2$Species[df2$Identity.percentage > 95]))

a <- unique(df$ESVs)


list_unknown <- c("unknown class","uncultured","unknown order","Unknown Family","unknown family",
                  "unknown genus")

f <- unique(df$ESVs[df$Identity.percentage > 90 & df$Family %not% list_unknown])
g <- unique(df$ESVs[df$Identity.percentage > 95 & df$Genus %not% list_unknown])

h <- sum(df$reads[df$Identity.percentage > 90 & df$Family %not% list_unknown], na.rm = TRUE)/sum(df$reads, na.rm = TRUE)

#--------Fig. 4b,c,d,e --- Prokaryotes NMDS - esv prok esvs------------------------
df <- read.csv("~downloads/Datasets/merged_prok_data.csv")

df$depth_water <- df$depth..below.the.local.surface..in.meters.

tsm_chl <- read.csv("~downloads/Datasets/filtered_tsm_chl_mean_islands.csv")

df <- merge(df, tsm_chl[,c(1,4,5)], by = "island", all = TRUE)


#rarefy to the sample (none negative) with the lowest read number)
test <- df %>% 
  group_by(field.nmbr.,Substrate.type) %>% 
  summarise(N_sum=sum(reads))
test1 <- test[test$N_sum > 1000,] #lowest read number for sample > 1000 = 1905
df1 <- merge(df, test1[,1:3], by = c("field.nmbr.", "Substrate.type"))

dfx <- df1 %>% 
  group_by(field.nmbr.) %>% 
  mutate(N=reads/sum(reads)*1905)

dfx1 <- dfx

df2 <- dfx1 %>% 
  group_by(field.nmbr., ESVs, Kingdom, Species,dist_coast_km, island, depth_water, substrate, Substrate.type, mean_chl_filtered, mean_tsm_filtered) %>% 
  summarise(N_mean=mean(N))


df2$reef <- NA
df2$reef[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- "flat"
df2$reef[df2$depth_water %not% c("RF1", "RF2", "RF3")] <- "slope"
df2$reef[is.na(df2$depth_water)] <- NA
df2$depth_water[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- 1
df2$depth_water <- as.numeric(df2$depth_water)
df2$depth_group <- NA
df2$depth_group[df2$depth_water < 10 & df2$reef == "slope"] <- "<10m"
df2$depth_group[df2$depth_water > 9 & df2$depth_water < 21] <- "10-20m"
df2$depth_group[df2$depth_water > 20 & df2$reef == "slope"] <- ">20m"

df2$turbidity <- NA
df2$turbidity[df2$mean_tsm_filtered > 1.5] <- "high turbidity"
df2$turbidity[df2$mean_tsm_filtered < 0.5] <- "low turbidity"
df2$turbidity[is.na(df2$turbidity)] <- "moderate turbidity"

df2$substrate[df2$substrate == "rubble + seagrass + algae"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "rubble + mud"] <- "rubble + sand"
df2$substrate[df2$substrate == "rubble + algae + seagrass"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "thick layer algae on rubble"] <- "rubble + algae"
df2$substrate[df2$substrate == "seagrass + algae + rubble + sand"] <- "rubble + seagrass"

df2a <- df2[df2$island != "negative",]
df2a <- df2a[!is.na(df2a$island),]

sum(dfx$N[dfx$Kingdom == "Bacteria "], na.rm = TRUE)/sum(dfx$N,na.rm = TRUE)
sum(dfx$N[dfx$Kingdom == "Archaea "], na.rm = TRUE)/sum(dfx$N,na.rm = TRUE)

sum(df2a$N_mean[df2a$Kingdom == "Bacteria "], na.rm = TRUE)/sum(dfx$N,na.rm = TRUE)
sum(df2a$N_mean[df2a$Kingdom == "Archaea "], na.rm = TRUE)/sum(dfx$N,na.rm = TRUE)

prok_tot <- unique(df1$ESVs)
prok_family <- unique(df2a$ESVs)
prok_group <- unique(df2b$ESVs)

prok_tot <- unique(df1$ESVs[df1$Kingdom == "Bacteria "])
prok_family <- unique(df2a$ESVs[df2a$Kingdom == "Bacteria "])
prok_group <- unique(df2b$ESVs[df2b$Kingdom == "Bacteria "])

prok_tot <- unique(df1$ESVs[df1$Kingdom == "Archaea "])
prok_family <- unique(df2a$ESVs[df2a$Kingdom == "Archaea "])
#prok_group <- unique(df2b$ESVs[df2b$Kingdom == "Archaea "])

write.csv(df2a,"~downloads/Datasets/prok_working_dataset_turnover.csv", row.names=FALSE)

#____reef flat - seawater_______
df2b <- df2a[df2a$reef == "flat" & df2a$Substrate.type == "seawater",]

df3 <- dcast(df2b, field.nmbr. + depth_group + island + dist_coast_km + reef + substrate + Substrate.type + turbidity ~ ESVs, 
             fun.aggregate = sum, value.var = "N_mean")

df3[df3 == "NaN"] <- 0

df_NMDS_vegan <- dcast(df2b, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
ord <- metaMDS(df_NMDS_vegan, distance = "bray")
env <- df3[,9:ncol(df3)]

#to calculate the species vectors going over the NMDS
en <- envfit(ord, env, permutations = 999, na.rm = TRUE)

#groups for plots
island <- droplevels(as.factor(df3$island))
depth <- as.factor(df3$depth_group)
field <- droplevels(as.factor(df3$field.nmbr.))
dist <- as.numeric(df3$dist_coast_km)
subst <- droplevels(as.factor(df3$substrate))
reef <- droplevels(as.factor(df3$reef))
type <- droplevels(as.factor(df3$Substrate.type))
turb <- droplevels(as.factor(df3$turbidity))

colfield <- hcl.colors(12, "Roma")
shapes <- c(19,17)

substord <- c("sand","rubble", "rubble + sand", "rubble + algae",
              "rubble + sand + algae", "rubble + seagrass")
islandord <- c("Lae Lae","Barang Baringan","Samalona","Polewali","Karanrang",
               "Langkadea","K. Keke","Padjenekang","Badi","Lumulumu","Langkai",
               "Kapoposang Masdar point")

turbidity <- c("high turbidity", "moderate turbidity", "low turbidity")
colturb <- c("#9E4E09" ,"#C7AE4C" ,"#008FA7" )
pchils <- c(0:6,8,15,17:19)

#Significance is the p value
ano <- anosim(df_NMDS_vegan, turb)

#for plot _ distance to shore
plot(ord, disp="sites", type = "n",main = paste ("prok Reef flat SW: turbidity levels - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3)))
ordisurf(ord, dist, col = "black",linetype = "dashed", add = TRUE)
ordihull(ord, groups=turb,  draw = "polygon", col = colturb[c(1,3,2)])
points(ord, disp="sites", pch = pchils[factor(island, islandord)] , cex=1, col = "black")
with(ord, legend(x = "topleft", legend = levels(factor(island, islandord)), col = "black", pch = pchils))
with(ord, legend(x = "bottomleft", legend = levels(factor(turb, turbidity)), col = colturb, pch = 15))


#____reef slope - seawater_______
df2b <- df2a[df2a$reef == "slope" & df2a$Substrate.type == "seawater",]

df3 <- dcast(df2b, field.nmbr. + depth_group + island + dist_coast_km + reef + substrate + Substrate.type + turbidity ~ ESVs, 
             fun.aggregate = sum, value.var = "N_mean")

df3[df3 == "NaN"] <- 0

df_NMDS_vegan <- dcast(df2b, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
ord <- metaMDS(df_NMDS_vegan, distance = "bray")
env <- df3[,9:ncol(df3)]

#to calculate the species vectors going over the NMDS
en <- envfit(ord, env, permutations = 999, na.rm = TRUE)

#groups for plots
island <- droplevels(as.factor(df3$island))
depth <- as.factor(df3$depth_group)
field <- droplevels(as.factor(df3$field.nmbr.))
dist <- as.numeric(df3$dist_coast_km)
subst <- droplevels(as.factor(df3$substrate))
reef <- droplevels(as.factor(df3$reef))
type <- droplevels(as.factor(df3$Substrate.type))
turb <- droplevels(as.factor(df3$turbidity))

colfield <- hcl.colors(12, "Roma")
shapes <- c(19,17)

substord <- c("sand","rubble", "rubble + sand", "rubble + algae",
              "rubble + sand + algae", "rubble + seagrass")
islandord <- c("Lae Lae","Barang Baringan","Samalona","Polewali","Karanrang",
               "Langkadea","K. Keke","Padjenekang","Badi","Lumulumu","Langkai",
               "Kapoposang Masdar point")

turbidity <- c("high turbidity", "moderate turbidity", "low turbidity")
colturb <- c("#9E4E09" ,"#C7AE4C" ,"#008FA7" )
pchils <- c(0:6,8,15,17:19)

#Significance is the p value
ano <- anosim(df_NMDS_vegan, turb)

#for plot _ distance to shore
plot(ord, disp="sites", type = "n",main = paste ("prok Reef slope SW: turbidity levels - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3)))
ordisurf(ord, dist, col = "black",linetype = "dashed", add = TRUE)
ordihull(ord, groups=turb,  draw = "polygon", col = colturb[c(1,3,2)])
points(ord, disp="sites", pch = pchils[factor(island, islandord)] , cex=1, col = "black")
with(ord, legend(x = "topleft", legend = levels(factor(island, islandord)), col = "black", pch = pchils))
with(ord, legend(x = "bottomleft", legend = levels(factor(turb, turbidity)), col = colturb, pch = 15))

#____reef flat - substrate_______
df2b <- df2a[df2a$reef == "flat" & df2a$Substrate.type == "substrate",]

df3 <- dcast(df2b, field.nmbr. + depth_group + island + dist_coast_km + reef + substrate + Substrate.type + turbidity ~ ESVs, 
             fun.aggregate = sum, value.var = "N_mean")

df3[df3 == "NaN"] <- 0

df_NMDS_vegan <- dcast(df2b, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
ord <- metaMDS(df_NMDS_vegan, distance = "bray")
env <- df3[,9:ncol(df3)]

#to calculate the species vectors going over the NMDS
en <- envfit(ord, env, permutations = 999, na.rm = TRUE)

#groups for plots
island <- droplevels(as.factor(df3$island))
depth <- as.factor(df3$depth_group)
field <- droplevels(as.factor(df3$field.nmbr.))
dist <- as.numeric(df3$dist_coast_km)
subst <- droplevels(as.factor(df3$substrate))
reef <- droplevels(as.factor(df3$reef))
type <- droplevels(as.factor(df3$Substrate.type))
turb <- droplevels(as.factor(df3$turbidity))

colfield <- hcl.colors(12, "Roma")
shapes <- c(19,17)

substord <- c("sand","rubble", "rubble + sand", "rubble + algae",
              "rubble + sand + algae", "rubble + seagrass")
islandord <- c("Lae Lae","Barang Baringan","Samalona","Polewali","Karanrang",
               "Langkadea","K. Keke","Padjenekang","Badi","Lumulumu","Langkai",
               "Kapoposang Masdar point")

turbidity <- c("high turbidity", "moderate turbidity", "low turbidity")
colturb <- c("#9E4E09" ,"#C7AE4C" ,"#008FA7" )
pchils <- c(0:6,8,15,17:19)

#Significance is the p value
ano <- anosim(df_NMDS_vegan, turb)

#for plot _ distance to shore
plot(ord, disp="sites", type = "n",main = paste ("prok Reef flat subs: turbidity levels - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3)))
ordisurf(ord, dist, col = "black",linetype = "dashed", add = TRUE)
ordihull(ord, groups=turb,  draw = "polygon", col = colturb[c(1,3,2)])
points(ord, disp="sites", pch = pchils[factor(island, islandord)] , cex=1, col = "black")
with(ord, legend(x = "topleft", legend = levels(factor(island, islandord)), col = "black", pch = pchils))
with(ord, legend(x = "bottomleft", legend = levels(factor(turb, turbidity)), col = colturb, pch = 15))


#____reef slope - substrate_______
df2b <- df2a[df2a$reef == "slope" & df2a$Substrate.type == "substrate",]

#impressive outlier, removed
#df2b <- df2b[df2b$real_ID != "e1100022118_NBCLAB10696_UPG.EG22.S107_LR033_D8_D07_16SPool_102190.1299_AGTGTCGCGAGTTGCATCGT_L001",]

df3 <- dcast(df2b, field.nmbr. + depth_group + island + dist_coast_km + reef + substrate + Substrate.type + turbidity ~ ESVs, 
             fun.aggregate = sum, value.var = "N_mean")

df3[df3 == "NaN"] <- 0

df_NMDS_vegan <- dcast(df2b, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
ord <- metaMDS(df_NMDS_vegan, distance = "bray")
env <- df3[,9:ncol(df3)]

#to calculate the species vectors going over the NMDS
en <- envfit(ord, env, permutations = 999, na.rm = TRUE)

#groups for plots
island <- droplevels(as.factor(df3$island))
depth <- as.factor(df3$depth_group)
field <- droplevels(as.factor(df3$field.nmbr.))
dist <- as.numeric(df3$dist_coast_km)
subst <- droplevels(as.factor(df3$substrate))
reef <- droplevels(as.factor(df3$reef))
type <- droplevels(as.factor(df3$Substrate.type))
turb <- droplevels(as.factor(df3$turbidity))

colfield <- hcl.colors(12, "Roma")
shapes <- c(19,17)

substord <- c("sand","rubble", "rubble + sand", "rubble + algae",
              "rubble + sand + algae", "rubble + seagrass")
islandord <- c("Lae Lae","Barang Baringan","Samalona","Polewali","Karanrang",
               "Langkadea","K. Keke","Padjenekang","Badi","Lumulumu","Langkai",
               "Kapoposang Masdar point")

turbidity <- c("high turbidity", "moderate turbidity", "low turbidity")
colturb <- c("#9E4E09" ,"#C7AE4C" ,"#008FA7" )
pchils <- c(0:6,8,15,17:19)

#Significance is the p value
ano <- anosim(df_NMDS_vegan, turb)

#for plot _ distance to shore
plot(ord, disp="sites", type = "n",main = paste ("prok Reef slope subs: turbidity levels - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3)))
ordisurf(ord, dist, col = "black",linetype = "dashed", add = TRUE)
ordihull(ord, groups=turb,  draw = "polygon", col = colturb[c(1,3,2)])
points(ord, disp="sites", pch = pchils[factor(island, islandord)] , cex=1, col = "black")
with(ord, legend(x = "topleft", legend = levels(factor(island, islandord)), col = "black", pch = pchils))
with(ord, legend(x = "bottomleft", legend = levels(factor(turb, turbidity)), col = colturb, pch = 15))



#____reef flat - water vs substrate_______
df2b <- df2a[df2a$reef %in% c("flat","slope"),]

df3 <- dcast(df2b, real_ID + field.nmbr. + depth_group + island + dist_coast_km + reef + substrate + Substrate.type ~ Species, 
             fun.aggregate = sum, value.var = "N_mean")

df3[df3 == "NaN"] <- 0

df_NMDS_vegan <- dcast(df2b, real_ID ~ Species, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$real_ID
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
ord <- metaMDS(df_NMDS_vegan, distance = "bray")
env <- df3[,8:ncol(df3)]

#to calculate the species vectors going over the NMDS
en <- envfit(ord, env, permutations = 999, na.rm = TRUE)

#groups for plots
island <- droplevels(as.factor(df3$island))
depth <- as.factor(df3$depth_group)
field <- droplevels(as.factor(df3$field.nmbr.))
dist <- as.numeric(df3$dist_coast_km)
subst <- droplevels(as.factor(df3$substrate))
reef <- droplevels(as.factor(df3$reef))
type <- droplevels(as.factor(df3$Substrate.type))


colfield <- hcl.colors(12, "Roma")
coltype <- c("#008FA7","#9E4E09")
shapes <- c(17,19)

substord <- c("sand","rubble", "rubble + sand", "rubble + algae",
              "rubble + sand + algae", "rubble + seagrass")
islandord <- c("Lae Lae","Barang Baringan","Samalona","Polewali","Karanrang",
               "Langkadea","K. Keke","Padjenekang","Badi","Lumulumu","Langkai",
               "Kapoposang Masdar point")

#Significance is the p value
ano <- anosim(df_NMDS_vegan, type)

#for plot
plot(ord, disp="sites", type = "n",main = paste ("Similarity between type - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3)))
#ordisurf(ord, dist, col = "black",linetype = "dashed", add = TRUE)
points(ord, disp="sites", pch = shapes[type], cex=1, col = coltype[type])
with(ord, legend(x = "topleft", legend = levels(type), col = coltype, pch = shapes))


#--------PROKARYOTES - community composition-------------------------
df <- read.csv("~downloads/Datasets/merged_prok_data.csv")

df$depth_water <- df$depth..below.the.local.surface..in.meters.

#rarefy to the sample (none negative) with the lowest read number: 1924)
dfx <- df %>% 
  group_by(field.nmbr.) %>% 
  mutate(N=reads/sum(reads)*1654)

dfx <- dfx %>% 
  group_by(field.nmbr.) %>% 
  mutate(percent_N=reads/sum(reads)*100)

df2 <- dfx %>% 
  group_by(field.nmbr., ESVs,dist_coast_km, island, depth_water, substrate, Substrate.type,Identity.percentage, Phylum, Class, Order, Family, Genus, Species) %>% 
  summarise(N_mean=mean(N), percent_N_mean = mean(percent_N))

df2$counts <- 1

sort(unique(df2$Phylum[df2$Identity.percentage > 86.5]))
sort(unique(df2$Class[df2$Identity.percentage > 86.5]))
sort(unique(df2$Order[df2$Identity.percentage > 86.5]))
sort(unique(df2$Family[df2$Identity.percentage > 86.5]))
sort(unique(df2$Genus[df2$Identity.percentage > 94.5]))

a <- unique(df$ESVs)
b <- unique(df$ESVs[df$Kingdom == "Bacteria "])
c <- unique(df$ESVs[df$Kingdom == "Archaea "])

d <- sum(df$reads[df$Kingdom == "Bacteria "], na.rm = TRUE)/sum(df$reads, na.rm = TRUE)
e <- sum(df$reads[df$Kingdom == "Archaea "], na.rm = TRUE)/sum(df$reads, na.rm = TRUE)

list_unknown <- c(" unknown class "," uncultured "," unknown order "," Unknown Family "," unknown family ",
                  " unknown genus ")

f <- unique(df$ESVs[df$Identity.percentage > 86.5 & df$Family %not% list_unknown])
g <- unique(df$ESVs[df$Identity.percentage > 94.5 & df$Genus %not% list_unknown])

h <- sum(df$reads[df$Identity.percentage > 86.5 & df$Family %not% list_unknown], na.rm = TRUE)/sum(df$reads, na.rm = TRUE)


#--------Table 2 --- db-RDA biplot LBF----------------------

df <- read.csv("~downloads/Datasets/merged_forams_data.csv")
tsm_chl <- read.csv("~downloads/Datasets/filtered_tsm_chl_mean_islands.csv")

df$depth_water <- df$depth..below.the.local.surface..in.meters.

df <- merge(df, tsm_chl[,c(1,4,5)], by = "island", all = TRUE)

#rarefy to the sample (none negative) with the lowest read number: 1924)
dfx <- df %>% 
  group_by(field.nmbr.) %>% 
  mutate(N=reads/sum(reads)*1924)

dfx1 <- dfx[dfx$group == "LBF" & dfx$identity > 99.4, ]

df2 <- dfx1 %>% 
  group_by(field.nmbr.,ESVs, species,dist_coast_km, island, mean_chl_filtered,mean_tsm_filtered, depth_water, substrate) %>% 
  summarise(N_mean=mean(N))

df2$reef <- NA
df2$reef[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- "flat"
df2$reef[df2$depth_water %not% c("RF1", "RF2", "RF3")] <- "slope"
df2$reef[df2$depth_water %in% c("negative")] <- NA
df2$depth_water[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- 1

df2$turbidity <- NA
df2$turbidity[df2$mean_tsm_filtered > 1.5] <- "high turbidity"
df2$turbidity[df2$mean_tsm_filtered < 0.5] <- "low turbidity"
df2$turbidity[is.na(df2$turbidity)] <- "moderate turbidity"

df2$substrate[df2$substrate == "rubble + seagrass + algae"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "rubble + mud"] <- "rubble + sand"
df2$substrate[df2$substrate == "rubble + algae + seagrass"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "thick layer algae on rubble"] <- "rubble + algae"
df2$substrate[df2$substrate == "seagrass + algae + rubble + sand"] <- "rubble + seagrass"

df2$morpho_species <- NA

df2$morpho_species[df2$species == "Alveolinella_quoyi_Spermonde"] <- "Alveolinella quoyi"
df2$morpho_species[df2$species == "Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS_Spermonde"] <- "Amphisorus sp."
df2$morpho_species[df2$species == "Amphistegina_lessonii"] <- "Amphistegina lessonii"
df2$morpho_species[df2$species == "Amphistegina_lessonii_5174"] <- "Amphistegina lessonii"
df2$morpho_species[df2$species == "Amphistegina_lobifera"] <- "Amphistegina lobifera"
df2$morpho_species[df2$species == "Amphistegina_papillosa_Rik_9440_consensus0_1051_SaudiArabia"] <- "Amphistegina papillosa"
df2$morpho_species[df2$species == "Amphistegina_papillosa_Spermonde_3741"] <- "Amphistegina papillosa"
df2$morpho_species[df2$species == "Amphistegina_radiata_Spermonde"] <- "Amphistegina radiata"
df2$morpho_species[df2$species == "Baculogypsinoides_spinosus_8078_8079_8080"] <- "Baculogypsinoides spinosus"
df2$morpho_species[df2$species == "Borelis_schlumbergeri_3737_Spermonde"] <- "Borelis schlumbergeri"
df2$morpho_species[df2$species == "Calcarina_hispida_&_Calcarina_sp_5163_Spermonde"] <- "Calcarina hispida"
df2$morpho_species[df2$species == "Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._5247_Spermonde"] <- "Calcarina spengleri"
df2$morpho_species[df2$species == "Calcarina_sp._5164_Spermonde"] <- "Calcarina sp."
df2$morpho_species[df2$species == "Calcarina_sp._5165_Spermonde"] <- "Calcarina sp."                           
df2$morpho_species[df2$species == "Calcarina_sp._Spermonde"] <- "Calcarina sp."                  
df2$morpho_species[df2$species == "Calcarina_spengleri_Rik_9352_consensus0_1732_Spermonde"] <- "Calcarina spengleri"   
df2$morpho_species[df2$species == "Cyrea_szymborska_17247"] <- "Z_other_foraminifera"   
df2$morpho_species[df2$species == "Eggerelloides_scaber_12301"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Elphidium_Rik_9327_9337_9350_Spermonde"] <- "Elphidium sp."
df2$morpho_species[df2$species == "Encrusting_Rik_9441_consensus0_977"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Glabratellina_sp._17921_&_Planoglabratella_opercularis_18053_&_Planoglabratella_sp"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Globocassidulina_biora_17225"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Haynesina_germanica_18209"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Heterostegina_depressa_sp2_Spermonde"] <- "Heterostegina depressa"
df2$morpho_species[df2$species == "Heterostegina_depressa_sp1_Spermonde"] <- "Heterostegina depressa"
df2$morpho_species[df2$species == "Miliolid_Rik_9400_9405_Israel"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Milliolid_red_Rik_9370_9371_9372_Israel"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Murrayinella_globosa_18043"] <- "Z_other_foraminifera"                                                         
df2$morpho_species[df2$species == "Neorotalia_calcar_Spermonde"] <- "Neorotalia calcar"
df2$morpho_species[df2$species == "Neorotalia_gaimardi_&_Baculogypsina_sphaerulata_Spermonde"] <- "Neorotalia gaimardi"
df2$morpho_species[df2$species == "Nummulites_venosus"] <- "Nummulites venosus"
df2$morpho_species[df2$species == "Operculina_ammonoides"] <- "Operculina ammonoides"
df2$morpho_species[df2$species == "Parasorites_sp"] <- "Parasorites sp."
df2$morpho_species[df2$species == "Peneroplis_sp1_Spermonde"] <- "Peneroplis sp."
df2$morpho_species[df2$species == "Peneroplis_sp2_&_Peneroplis_pertusus_5117_&_Dendritina_ambigua_Spermonde"] <- "Peneroplis sp."
df2$morpho_species[df2$species == "Planorbulinella_sp._17893"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Psammophaga_sp._19260_&_Psammophaga_sp._19296_&_Psammophaga_sp2"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Pyrgo_Rik_9347_consensus0_1557"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Sorites_sp1_Spermonde"] <- "Sorites sp."
df2$morpho_species[df2$species == "Sorites_sp2_Spermonde"] <- "Sorites sp."
df2$morpho_species[df2$species == "Stainforthia_sp._17321"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Textularid_Rik_9397_consensus4_105"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Textularid_Rik_9398_consensus128_73"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Trifarina_earlandi_&_Uvigerina_bifurcata"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Vertebralina_Rik_9325_consensus0_1751"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Z_not_foraminifera"] <- "Z_not_foraminifera" 
df2$morpho_species[df2$species == "Z_other_foraminifera"] <- "Z_other_foraminifera" 

unique(df2$morpho_species)

df2a <- df2[df2$island != "negative",]
df2a <- df2a[!is.na(df2a$island),]
df2a$depth_water <- as.numeric(df2a$depth_water)


#_______________________________prep dataset with only 5% of the samples______________

df_NMDS_vegan <- dcast(df2a, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

percent_df <- decostand(df_NMDS_vegan, "total")
pa_df <- df_NMDS_vegan
pa_df[pa_df > 0] <- 1

colSums(percent_df)
esv_5_sample <- colnames(pa_df[colSums(pa_df)/109 * 100 >= 5]) #in at least 5% of the samples

#resulting ESVs fitting the conditions

df2d <- df2a[df2a$ESVs %in% esv_5_sample,]




#____reef flat_______
df2b <- df2a[df2a$reef == "flat",]

df3 <- dcast(df2b, field.nmbr. + depth_water + island + dist_coast_km + mean_chl_filtered + mean_tsm_filtered + reef + substrate ~ ESVs, 
             fun.aggregate = sum, value.var = "N_mean")

df3[df3 == "NaN"] <- 0

df_NMDS_vegan <- dcast(df2b, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

env <- df3[,1:8]
cca_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns

#CCA
#ord_dist_sub <- cca(cca_vegan~dist_coast_km + substrate + mean_chl_filtered + mean_tsm_filtered, data = env)
#ord_dist_sub <- cca(cca_vegan~dist_coast_km + mean_chl_filtered + mean_tsm_filtered, data = env)
#ord_dist_sub <- cca(cca_vegan~substrate, data = env)
#summary(eigenvals(ord_dist_sub))
#summary((ord_dist_sub))
#plot(ord_dist_sub)
#autoplot(ord_dist_sub) + theme_bw() + 
#  xlab("CCA1: 12.77%") + ylab("CCA2: 11.34%")+
#  ggtitle("Reef flat forams")


#db-RDA analysis
#all
dbRDA <- capscale(cca_vegan~dist_coast_km + substrate + mean_chl_filtered + mean_tsm_filtered, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables

#eutrophication only
dbRDA <- capscale(cca_vegan~dist_coast_km +  mean_chl_filtered + mean_tsm_filtered, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables

#substrate only
dbRDA <- capscale(cca_vegan~substrate, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables


#____reef slope_______
df2b <- df2a[df2a$reef == "slope",]

df3 <- dcast(df2b, field.nmbr. + depth_water + island + dist_coast_km + mean_chl_filtered + mean_tsm_filtered + reef + substrate ~ ESVs, 
             fun.aggregate = sum, value.var = "N_mean")

df3[df3 == "NaN"] <- 0

df_NMDS_vegan <- dcast(df2b, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

env <- df3[,1:8]
cca_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns

#db-RDA analysis
#all
dbRDA <- capscale(cca_vegan~dist_coast_km + substrate + mean_chl_filtered + mean_tsm_filtered + depth_water, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables

#eutrophication only
dbRDA <- capscale(cca_vegan~dist_coast_km +  mean_chl_filtered + mean_tsm_filtered, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables

#substrate only
dbRDA <- capscale(cca_vegan~substrate, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables

#water depth only
dbRDA <- capscale(cca_vegan~depth_water, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables


#--------Table 2 --- NEW_db-RDA biplot all forams----------------------

df <- read.csv("~downloads/Datasets/merged_forams_data.csv")
tsm_chl <- read.csv("~downloads/Datasets/filtered_tsm_chl_mean_islands.csv")

df$depth_water <- df$depth..below.the.local.surface..in.meters.

df <- merge(df, tsm_chl[,c(1,4,5)], by = "island", all = TRUE)

#rarefy to the sample (none negative) with the lowest read number: 1924)
dfx <- df %>% 
  group_by(field.nmbr.) %>% 
  mutate(N=reads/sum(reads)*1924)

dfx1 <- dfx

df2 <- dfx1 %>% 
  group_by(field.nmbr.,ESVs, species,dist_coast_km, island, mean_chl_filtered,mean_tsm_filtered, depth_water, substrate) %>% 
  summarise(N_mean=mean(N))

df2$reef <- NA
df2$reef[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- "flat"
df2$reef[df2$depth_water %not% c("RF1", "RF2", "RF3")] <- "slope"
df2$reef[df2$depth_water %in% c("negative")] <- NA
df2$depth_water[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- 1

df2$turbidity <- NA
df2$turbidity[df2$mean_tsm_filtered > 1.5] <- "high turbidity"
df2$turbidity[df2$mean_tsm_filtered < 0.5] <- "low turbidity"
df2$turbidity[is.na(df2$turbidity)] <- "moderate turbidity"

df2$substrate[df2$substrate == "rubble + seagrass + algae"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "rubble + mud"] <- "rubble + sand"
df2$substrate[df2$substrate == "rubble + algae + seagrass"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "thick layer algae on rubble"] <- "rubble + algae"
df2$substrate[df2$substrate == "seagrass + algae + rubble + sand"] <- "rubble + seagrass"

df2$morpho_species <- NA

df2$morpho_species[df2$species == "Alveolinella_quoyi_Spermonde"] <- "Alveolinella quoyi"
df2$morpho_species[df2$species == "Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS_Spermonde"] <- "Amphisorus sp."
df2$morpho_species[df2$species == "Amphistegina_lessonii"] <- "Amphistegina lessonii"
df2$morpho_species[df2$species == "Amphistegina_lessonii_5174"] <- "Amphistegina lessonii"
df2$morpho_species[df2$species == "Amphistegina_lobifera"] <- "Amphistegina lobifera"
df2$morpho_species[df2$species == "Amphistegina_papillosa_Rik_9440_consensus0_1051_SaudiArabia"] <- "Amphistegina papillosa"
df2$morpho_species[df2$species == "Amphistegina_papillosa_Spermonde_3741"] <- "Amphistegina papillosa"
df2$morpho_species[df2$species == "Amphistegina_radiata_Spermonde"] <- "Amphistegina radiata"
df2$morpho_species[df2$species == "Baculogypsinoides_spinosus_8078_8079_8080"] <- "Baculogypsinoides spinosus"
df2$morpho_species[df2$species == "Borelis_schlumbergeri_3737_Spermonde"] <- "Borelis schlumbergeri"
df2$morpho_species[df2$species == "Calcarina_hispida_&_Calcarina_sp_5163_Spermonde"] <- "Calcarina hispida"
df2$morpho_species[df2$species == "Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._5247_Spermonde"] <- "Calcarina spengleri"
df2$morpho_species[df2$species == "Calcarina_sp._5164_Spermonde"] <- "Calcarina sp."
df2$morpho_species[df2$species == "Calcarina_sp._5165_Spermonde"] <- "Calcarina sp."                           
df2$morpho_species[df2$species == "Calcarina_sp._Spermonde"] <- "Calcarina sp."                  
df2$morpho_species[df2$species == "Calcarina_spengleri_Rik_9352_consensus0_1732_Spermonde"] <- "Calcarina spengleri"   
df2$morpho_species[df2$species == "Cyrea_szymborska_17247"] <- "Z_other_foraminifera"   
df2$morpho_species[df2$species == "Eggerelloides_scaber_12301"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Elphidium_Rik_9327_9337_9350_Spermonde"] <- "Elphidium sp."
df2$morpho_species[df2$species == "Encrusting_Rik_9441_consensus0_977"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Glabratellina_sp._17921_&_Planoglabratella_opercularis_18053_&_Planoglabratella_sp"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Globocassidulina_biora_17225"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Haynesina_germanica_18209"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Heterostegina_depressa_sp2_Spermonde"] <- "Heterostegina depressa"
df2$morpho_species[df2$species == "Heterostegina_depressa_sp1_Spermonde"] <- "Heterostegina depressa"
df2$morpho_species[df2$species == "Miliolid_Rik_9400_9405_Israel"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Milliolid_red_Rik_9370_9371_9372_Israel"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Murrayinella_globosa_18043"] <- "Z_other_foraminifera"                                                         
df2$morpho_species[df2$species == "Neorotalia_calcar_Spermonde"] <- "Neorotalia calcar"
df2$morpho_species[df2$species == "Neorotalia_gaimardi_&_Baculogypsina_sphaerulata_Spermonde"] <- "Neorotalia gaimardi"
df2$morpho_species[df2$species == "Nummulites_venosus"] <- "Nummulites venosus"
df2$morpho_species[df2$species == "Operculina_ammonoides"] <- "Operculina ammonoides"
df2$morpho_species[df2$species == "Parasorites_sp"] <- "Parasorites sp."
df2$morpho_species[df2$species == "Peneroplis_sp1_Spermonde"] <- "Peneroplis sp."
df2$morpho_species[df2$species == "Peneroplis_sp2_&_Peneroplis_pertusus_5117_&_Dendritina_ambigua_Spermonde"] <- "Peneroplis sp."
df2$morpho_species[df2$species == "Planorbulinella_sp._17893"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Psammophaga_sp._19260_&_Psammophaga_sp._19296_&_Psammophaga_sp2"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Pyrgo_Rik_9347_consensus0_1557"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Sorites_sp1_Spermonde"] <- "Sorites sp."
df2$morpho_species[df2$species == "Sorites_sp2_Spermonde"] <- "Sorites sp."
df2$morpho_species[df2$species == "Stainforthia_sp._17321"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Textularid_Rik_9397_consensus4_105"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Textularid_Rik_9398_consensus128_73"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Trifarina_earlandi_&_Uvigerina_bifurcata"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Vertebralina_Rik_9325_consensus0_1751"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Z_not_foraminifera"] <- "Z_not_foraminifera" 
df2$morpho_species[df2$species == "Z_other_foraminifera"] <- "Z_other_foraminifera" 

unique(df2$morpho_species)

df2a <- df2[df2$island != "negative",]
df2a <- df2a[!is.na(df2a$island),]
df2a$depth_water <- as.numeric(df2a$depth_water)


#_______________________________prep dataset with only 5% of the samples______________

df_NMDS_vegan <- dcast(df2a, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

percent_df <- decostand(df_NMDS_vegan, "total")
pa_df <- df_NMDS_vegan
pa_df[pa_df > 0] <- 1

colSums(percent_df)
esv_5_sample <- colnames(pa_df[colSums(pa_df)/109 * 100 >= 5]) #in at least 5% of the samples

#resulting ESVs fitting the conditions

df2d <- df2a[df2a$ESVs %in% esv_5_sample,]




#____reef flat_______
df2b <- df2a[df2a$reef == "flat",]

df3 <- dcast(df2b, field.nmbr. + depth_water + island + dist_coast_km + mean_chl_filtered + mean_tsm_filtered + reef + substrate ~ ESVs, 
             fun.aggregate = sum, value.var = "N_mean")

df3[df3 == "NaN"] <- 0

df_NMDS_vegan <- dcast(df2b, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

env <- df3[,1:8]
cca_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns

#CCA
#ord_dist_sub <- cca(cca_vegan~dist_coast_km + substrate + mean_chl_filtered + mean_tsm_filtered, data = env)
#ord_dist_sub <- cca(cca_vegan~dist_coast_km + mean_chl_filtered + mean_tsm_filtered, data = env)
#ord_dist_sub <- cca(cca_vegan~substrate, data = env)
#summary(eigenvals(ord_dist_sub))
#summary((ord_dist_sub))
#plot(ord_dist_sub)
#autoplot(ord_dist_sub) + theme_bw() + 
#  xlab("CCA1: 12.77%") + ylab("CCA2: 11.34%")+
#  ggtitle("Reef flat forams")


#db-RDA analysis
#all
dbRDA <- capscale(cca_vegan~dist_coast_km + substrate + mean_chl_filtered + mean_tsm_filtered, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables

#eutrophication only
dbRDA <- capscale(cca_vegan~dist_coast_km +  mean_chl_filtered + mean_tsm_filtered, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables

#substrate only
dbRDA <- capscale(cca_vegan~substrate, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables


#____reef slope_______
df2b <- df2a[df2a$reef == "slope",]

df3 <- dcast(df2b, field.nmbr. + depth_water + island + dist_coast_km + mean_chl_filtered + mean_tsm_filtered + reef + substrate ~ ESVs, 
             fun.aggregate = sum, value.var = "N_mean")

df3[df3 == "NaN"] <- 0

df_NMDS_vegan <- dcast(df2b, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

env <- df3[,1:8]
cca_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns

#db-RDA analysis
#all
dbRDA <- capscale(cca_vegan~dist_coast_km + substrate + mean_chl_filtered + mean_tsm_filtered + depth_water, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables

#eutrophication only
dbRDA <- capscale(cca_vegan~dist_coast_km +  mean_chl_filtered + mean_tsm_filtered, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables

#substrate only
dbRDA <- capscale(cca_vegan~substrate, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables

#water depth only
dbRDA <- capscale(cca_vegan~depth_water, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables
#--------Table 2 --- NEW_db-RDA biplot all diatoms----------------------

df <- read.csv("~downloads/Datasets/merged_diat_data.csv")

df$depth_water <- df$depth..below.the.local.surface..in.meters.

tsm_chl <- read.csv("~downloads/Datasets/filtered_tsm_chl_mean_islands.csv")

df <- merge(df, tsm_chl[,c(1,4,5)], by = "island", all = TRUE)


#rarefy to the sample (none negative) with the lowest read number)
test <- df %>% 
  group_by(field.nmbr.,Substrate.type) %>% 
  summarise(N_sum=sum(reads))
test1 <- test[test$N_sum > 1000,] #lowest read number for sample > 1000 = 1654
df1 <- merge(df, test1[,1:3], by = c("field.nmbr.", "Substrate.type"))

dfx <- df1 %>% 
  group_by(field.nmbr.) %>% 
  mutate(N=reads/sum(reads)*1654)

dfx <- dfx[dfx$Phylum == "Bacillariophyta",] #assigned only to diatoms

dfx1 <- dfx


df2 <- dfx1 %>% 
  group_by(field.nmbr., ESVs, Species,dist_coast_km, island, mean_chl_filtered,mean_tsm_filtered, depth_water, substrate, Substrate.type) %>% 
  summarise(N_mean=mean(N))


df2$reef <- NA
df2$reef[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- "flat"
df2$reef[df2$depth_water %not% c("RF1", "RF2", "RF3")] <- "slope"
df2$reef[is.na(df2$depth_water)] <- NA
df2$depth_water[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- 1
df2$depth_water <- as.numeric(df2$depth_water)
df2$depth_group <- NA
df2$depth_group[df2$depth_water < 10 & df2$reef == "slope"] <- "<10m"
df2$depth_group[df2$depth_water > 9 & df2$depth_water < 21] <- "10-20m"
df2$depth_group[df2$depth_water > 20 & df2$reef == "slope"] <- ">20m"

df2$turbidity <- NA
df2$turbidity[df2$mean_tsm_filtered > 1.5] <- "high turbidity"
df2$turbidity[df2$mean_tsm_filtered < 0.5] <- "low turbidity"
df2$turbidity[is.na(df2$turbidity)] <- "moderate turbidity"

df2$substrate[df2$substrate == "rubble + seagrass + algae"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "rubble + mud"] <- "rubble + sand"
df2$substrate[df2$substrate == "rubble + algae + seagrass"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "thick layer algae on rubble"] <- "rubble + algae"
df2$substrate[df2$substrate == "seagrass + algae + rubble + sand"] <- "rubble + seagrass"

df2a <- df2[df2$island != "negative",]
df2a <- df2a[!is.na(df2a$island),]
df2a <- df2a[ df2a$Substrate.type == "substrate",]
df2a$depth_water <- as.numeric(df2a$depth_water)


#_______________________________prep dataset 5% of the samples

df_NMDS_vegan <- dcast(df2a, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

percent_df <- decostand(df_NMDS_vegan, "total")
pa_df <- df_NMDS_vegan
pa_df[pa_df > 0] <- 1

esv_5_sample <- colnames(pa_df[colSums(pa_df)/109 * 100 >= 5]) #in at least 5% of the samples

df2a <- df2a[df2a$ESVs %in% esv_5_sample,]



#____reef flat_______
df2b <- df2a[df2a$reef == "flat",]

df3 <- dcast(df2b, field.nmbr. + depth_water + island + dist_coast_km + mean_chl_filtered + mean_tsm_filtered + reef + substrate ~ ESVs, 
             fun.aggregate = sum, value.var = "N_mean")

df3[df3 == "NaN"] <- 0

df_NMDS_vegan <- dcast(df2b, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

env <- df3[,1:8]
cca_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns

#db-RDA analysis
#all
dbRDA <- capscale(cca_vegan~dist_coast_km + substrate + mean_chl_filtered + mean_tsm_filtered, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables

#eutrophication only
dbRDA <- capscale(cca_vegan~dist_coast_km +  mean_chl_filtered + mean_tsm_filtered, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables

#substrate only
dbRDA <- capscale(cca_vegan~substrate, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables

#____reef slope_______
df2b <- df2a[df2a$reef == "slope",]

df3 <- dcast(df2b, field.nmbr. + depth_water + island + dist_coast_km + mean_chl_filtered + mean_tsm_filtered + reef + substrate ~ ESVs, 
             fun.aggregate = sum, value.var = "N_mean")

df3[df3 == "NaN"] <- 0

df_NMDS_vegan <- dcast(df2b, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

env <- df3[,1:8]
cca_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns

#db-RDA analysis
#all
dbRDA <- capscale(cca_vegan~dist_coast_km + substrate + mean_chl_filtered + mean_tsm_filtered + depth_water, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables

#eutrophication only
dbRDA <- capscale(cca_vegan~dist_coast_km +  mean_chl_filtered + mean_tsm_filtered, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables

#substrate only
dbRDA <- capscale(cca_vegan~substrate, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables

#water depth only
dbRDA <- capscale(cca_vegan~depth_water, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables

#--------Table 2 --- NEW_db-RDA biplot all prok----------------------

df <- read.csv("~downloads/Datasets/merged_prok_data.csv")

df$depth_water <- df$depth..below.the.local.surface..in.meters.

tsm_chl <- read.csv("~downloads/Datasets/filtered_tsm_chl_mean_islands.csv")

df <- merge(df, tsm_chl[,c(1,4,5)], by = "island", all = TRUE)

#rarefy to the sample (none negative) with the lowest read number)
test <- df %>% 
  group_by(field.nmbr.,Substrate.type) %>% 
  summarise(N_sum=sum(reads))
test1 <- test[test$N_sum > 1000,] #lowest read number for sample > 1000 = 1905
df1 <- merge(df, test1[,1:3], by = c("field.nmbr.", "Substrate.type"))

dfx <- df1 %>% 
  group_by(field.nmbr.) %>% 
  mutate(N=reads/sum(reads)*1905)

dfx1 <- dfx

df2 <- dfx1 %>% 
  group_by(field.nmbr., ESVs, Species,dist_coast_km, island, mean_chl_filtered,mean_tsm_filtered, depth_water, substrate, Substrate.type) %>% 
  summarise(N_mean=mean(N))


df2$reef <- NA
df2$reef[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- "flat"
df2$reef[df2$depth_water %not% c("RF1", "RF2", "RF3")] <- "slope"
df2$reef[is.na(df2$depth_water)] <- NA
df2$depth_water[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- 1
df2$depth_water <- as.numeric(df2$depth_water)
df2$depth_group <- NA
df2$depth_group[df2$depth_water < 10 & df2$reef == "slope"] <- "<10m"
df2$depth_group[df2$depth_water > 9 & df2$depth_water < 21] <- "10-20m"
df2$depth_group[df2$depth_water > 20 & df2$reef == "slope"] <- ">20m"

df2$turbidity <- NA
df2$turbidity[df2$mean_tsm_filtered > 1.5] <- "high turbidity"
df2$turbidity[df2$mean_tsm_filtered < 0.5] <- "low turbidity"
df2$turbidity[is.na(df2$turbidity)] <- "moderate turbidity"

df2$substrate[df2$substrate == "rubble + seagrass + algae"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "rubble + mud"] <- "rubble + sand"
df2$substrate[df2$substrate == "rubble + algae + seagrass"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "thick layer algae on rubble"] <- "rubble + algae"
df2$substrate[df2$substrate == "seagrass + algae + rubble + sand"] <- "rubble + seagrass"

df2a <- df2[df2$island != "negative",]
df2a <- df2a[!is.na(df2a$island),]
df2a <- df2a[ df2a$Substrate.type == "substrate",]
df2a$depth_water <- as.numeric(df2a$depth_water)


#_______________________________prep dataset 5% of the samples

df_NMDS_vegan <- dcast(df2a, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

percent_df <- decostand(df_NMDS_vegan, "total")
pa_df <- df_NMDS_vegan
pa_df[pa_df > 0] <- 1

esv_5_sample <- colnames(pa_df[colSums(pa_df)/109 * 100 >= 5]) #in at least 5% of the samples

#resulting ESVs fitting the conditions
df2a <- df2a[df2a$ESVs %in% esv_5_sample,]



#____reef flat_______
df2b <- df2a[df2a$reef == "flat",]

df3 <- dcast(df2b, field.nmbr. + depth_water + island + dist_coast_km + mean_chl_filtered + mean_tsm_filtered + reef + substrate ~ ESVs, 
             fun.aggregate = sum, value.var = "N_mean")

df3[df3 == "NaN"] <- 0

df_NMDS_vegan <- dcast(df2b, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

env <- df3[,1:8]
cca_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns

#db-RDA analysis
#all
dbRDA <- capscale(cca_vegan~dist_coast_km + substrate + mean_chl_filtered + mean_tsm_filtered, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables

#eutrophication only
dbRDA <- capscale(cca_vegan~dist_coast_km +  mean_chl_filtered + mean_tsm_filtered, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables

#substrate only
dbRDA <- capscale(cca_vegan~substrate, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables


#____reef slope_______
df2b <- df2a[df2a$reef == "slope",]

df3 <- dcast(df2b, field.nmbr. + depth_water + island + dist_coast_km + mean_chl_filtered + mean_tsm_filtered + reef + substrate ~ ESVs, 
             fun.aggregate = sum, value.var = "N_mean")

df3[df3 == "NaN"] <- 0

df_NMDS_vegan <- dcast(df2b, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[df_NMDS_vegan == "NaN"] <- 0

env <- df3[,1:8]
cca_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns

#db-RDA analysis
#all
dbRDA <- capscale(cca_vegan~dist_coast_km + substrate + mean_chl_filtered + mean_tsm_filtered + depth_water, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables

#eutrophication only
dbRDA <- capscale(cca_vegan~dist_coast_km +  mean_chl_filtered + mean_tsm_filtered, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables

#substrate only
dbRDA <- capscale(cca_vegan~substrate, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables

#water depth only
dbRDA <- capscale(cca_vegan~depth_water, data = env, dist="bray")
plot(dbRDA)
summary(dbRDA)
anova(dbRDA) # overall test of the significant of the analysis
anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=200) # test for sign. environ. variables


#--------check normality of the data------------------

#The R function shapiro.test() can be used 
#to perform the Shapiro-Wilk test of normality for one variable (univariate):
#the p-value > 0.05 implying that the distribution of the data are not 
#significantly different from normal distribution. In other words, we can assume the normality.

shapiro.test(df2a$dist_coast_km) #p < 0.001 NOT normal
shapiro.test(df2a$mean_chl_filtered) #p < 0.001 NOT normal
shapiro.test(df2a$mean_tsm_filtered) #p < 0.001 NOT normal
shapiro.test(df2a$depth_water) #p < 0.001 NOT normal


#--------Indicator species analysis lbf-------------------------
df <- read.csv("~downloads/Datasets/merged_forams_data.csv")
tsm_chl <- read.csv("~downloads/Datasets/filtered_tsm_chl_mean_islands.csv")

df$depth_water <- df$depth..below.the.local.surface..in.meters.

df <- merge(df, tsm_chl[,c(1,4,5)], by = "island", all = TRUE)

#rarefy to the sample (none negative) with the lowest read number: 1924)
dfx <- df %>% 
  group_by(field.nmbr.) %>% 
  mutate(N=reads/sum(reads)*1924)

dfx1 <- dfx[dfx$group == "LBF" & dfx$identity > 99.4, ]

df2 <- dfx1 %>% 
  group_by(field.nmbr.,ESVs, species,dist_coast_km, island, mean_chl_filtered,mean_tsm_filtered, depth_water, substrate) %>% 
  summarise(N_mean=mean(N))

df2$reef <- NA
df2$reef[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- "flat"
df2$reef[df2$depth_water %not% c("RF1", "RF2", "RF3")] <- "slope"
df2$reef[df2$depth_water %in% c("negative")] <- NA
df2$depth_water[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- 1

df2$turbidity <- NA
df2$turbidity[df2$mean_tsm_filtered > 1.5] <- "high turbidity"
df2$turbidity[df2$mean_tsm_filtered < 0.5] <- "low turbidity"
df2$turbidity[is.na(df2$turbidity)] <- "moderate turbidity"

df2$substrate[df2$substrate == "rubble + seagrass + algae"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "rubble + mud"] <- "rubble + sand"
df2$substrate[df2$substrate == "rubble + algae + seagrass"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "thick layer algae on rubble"] <- "rubble + algae"
df2$substrate[df2$substrate == "seagrass + algae + rubble + sand"] <- "rubble + seagrass"

df2$morpho_species <- NA

df2$morpho_species[df2$species == "Alveolinella_quoyi_Spermonde"] <- "Alveolinella quoyi"
df2$morpho_species[df2$species == "Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS_Spermonde"] <- "Amphisorus sp."
df2$morpho_species[df2$species == "Amphistegina_lessonii"] <- "Amphistegina lessonii"
df2$morpho_species[df2$species == "Amphistegina_lessonii_5174"] <- "Amphistegina lessonii"
df2$morpho_species[df2$species == "Amphistegina_lobifera"] <- "Amphistegina lobifera"
df2$morpho_species[df2$species == "Amphistegina_papillosa_Rik_9440_consensus0_1051_SaudiArabia"] <- "Amphistegina papillosa"
df2$morpho_species[df2$species == "Amphistegina_papillosa_Spermonde_3741"] <- "Amphistegina papillosa"
df2$morpho_species[df2$species == "Amphistegina_radiata_Spermonde"] <- "Amphistegina radiata"
df2$morpho_species[df2$species == "Baculogypsinoides_spinosus_8078_8079_8080"] <- "Baculogypsinoides spinosus"
df2$morpho_species[df2$species == "Borelis_schlumbergeri_3737_Spermonde"] <- "Borelis schlumbergeri"
df2$morpho_species[df2$species == "Calcarina_hispida_&_Calcarina_sp_5163_Spermonde"] <- "Calcarina hispida"
df2$morpho_species[df2$species == "Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._5247_Spermonde"] <- "Calcarina spengleri"
df2$morpho_species[df2$species == "Calcarina_sp._5164_Spermonde"] <- "Calcarina sp."
df2$morpho_species[df2$species == "Calcarina_sp._5165_Spermonde"] <- "Calcarina sp."                           
df2$morpho_species[df2$species == "Calcarina_sp._Spermonde"] <- "Calcarina sp."                  
df2$morpho_species[df2$species == "Calcarina_spengleri_Rik_9352_consensus0_1732_Spermonde"] <- "Calcarina spengleri"   
df2$morpho_species[df2$species == "Cyrea_szymborska_17247"] <- "Z_other_foraminifera"   
df2$morpho_species[df2$species == "Eggerelloides_scaber_12301"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Elphidium_Rik_9327_9337_9350_Spermonde"] <- "Elphidium sp."
df2$morpho_species[df2$species == "Encrusting_Rik_9441_consensus0_977"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Glabratellina_sp._17921_&_Planoglabratella_opercularis_18053_&_Planoglabratella_sp"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Globocassidulina_biora_17225"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Haynesina_germanica_18209"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Heterostegina_depressa_sp2_Spermonde"] <- "Heterostegina depressa"
df2$morpho_species[df2$species == "Heterostegina_depressa_sp1_Spermonde"] <- "Heterostegina depressa"
df2$morpho_species[df2$species == "Miliolid_Rik_9400_9405_Israel"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Milliolid_red_Rik_9370_9371_9372_Israel"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Murrayinella_globosa_18043"] <- "Z_other_foraminifera"                                                         
df2$morpho_species[df2$species == "Neorotalia_calcar_Spermonde"] <- "Neorotalia calcar"
df2$morpho_species[df2$species == "Neorotalia_gaimardi_&_Baculogypsina_sphaerulata_Spermonde"] <- "Neorotalia gaimardi"
df2$morpho_species[df2$species == "Nummulites_venosus"] <- "Nummulites venosus"
df2$morpho_species[df2$species == "Operculina_ammonoides"] <- "Operculina ammonoides"
df2$morpho_species[df2$species == "Parasorites_sp"] <- "Parasorites sp."
df2$morpho_species[df2$species == "Peneroplis_sp1_Spermonde"] <- "Peneroplis sp."
df2$morpho_species[df2$species == "Peneroplis_sp2_&_Peneroplis_pertusus_5117_&_Dendritina_ambigua_Spermonde"] <- "Peneroplis sp."
df2$morpho_species[df2$species == "Planorbulinella_sp._17893"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Psammophaga_sp._19260_&_Psammophaga_sp._19296_&_Psammophaga_sp2"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Pyrgo_Rik_9347_consensus0_1557"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Sorites_sp1_Spermonde"] <- "Sorites sp."
df2$morpho_species[df2$species == "Sorites_sp2_Spermonde"] <- "Sorites sp."
df2$morpho_species[df2$species == "Stainforthia_sp._17321"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Textularid_Rik_9397_consensus4_105"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Textularid_Rik_9398_consensus128_73"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Trifarina_earlandi_&_Uvigerina_bifurcata"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Vertebralina_Rik_9325_consensus0_1751"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Z_not_foraminifera"] <- "Z_not_foraminifera" 
df2$morpho_species[df2$species == "Z_other_foraminifera"] <- "Z_other_foraminifera" 

unique(df2$morpho_species)

df2a <- df2[df2$island != "negative",]
df2a <- df2a[!is.na(df2a$island),]
df2a$depth_water <- as.numeric(df2a$depth_water)


#_______________________________prep dataset 5% of the samples

df_NMDS_vegan <- dcast(df2a, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[is.na(df_NMDS_vegan)] <- 0

percent_df <- decostand(df_NMDS_vegan, "total")
pa_df <- df_NMDS_vegan
pa_df[pa_df > 0] <- 1

colSums(percent_df, na.rm = TRUE)
esv_5_sample <- colnames(pa_df[colSums(pa_df)/109 * 100 >= 5]) #in at least 5% of the samples

df2a <- df2a[df2a$ESVs %in% esv_5_sample,]

slope <- subset(df2a, reef == "slope")
flat <- subset(df2a, reef == "flat")
flat <- flat[flat$field.nmbr. != "UPG-EG22-053",]

#"unmelt" the dataframe with reshape and
#to isolate values (N), Species code (x) and Sample codes (y) with col and row names
slope_indicator <- dcast(slope, field.nmbr. + island + substrate + turbidity ~ ESVs, 
                         fun.aggregate = sum, value.var = "N_mean")
slope_indicator[is.na(slope_indicator)] <- 0
slope_species <- slope_indicator[,-c(2:4)]
rownames(slope_species) <- slope_species$Sample_Code
slope_species <- slope_species[,-1]
turbidity_slope <- slope_indicator$turbidity

flat_indicator <- dcast(flat, field.nmbr. + island + substrate + turbidity ~ ESVs, 
                         fun.aggregate = sum, value.var = "N_mean")
flat_indicator[is.na(flat_indicator)] <- 0
flat_species <- flat_indicator[,-c(2:4)]
rownames(flat_species) <- flat_species$Sample_Code
flat_species <- flat_species[,-1]
turbidity_flat <- flat_indicator$turbidity

#________________________________________Indicator species analysis using IndVal.g - ABUNDANCE

#Indicator species analysis for turbidity levels
indval_s<- multipatt(slope_species, turbidity_slope, control = how(nperm=2000))
summary(indval_s)
ESVs_slope <- indval_s$sign[indval_s$sign$p.value < 0.01 & !is.na(indval_s$sign$p.value),]
ESVs_slope$reef <- "slope"
ESVs_slope$ESVs <- rownames(ESVs_slope)

indval_f <- multipatt(flat_species, turbidity_flat, control = how(nperm=2000))
summary(indval_f)
ESVs_flat <- indval_f$sign[indval_f$sign$p.value < 0.01 & !is.na(indval_f$sign$p.value),]
ESVs_flat$reef <- "flat"
ESVs_flat$ESVs <- rownames(ESVs_flat)


foram_res <- rbind(ESVs_slope, ESVs_flat)

foram_res_indicspecies <- merge(unique(df[c("ESVs", "species")]), foram_res, by = "ESVs")

write.csv(foram_res_indicspecies, "~downloads/Datasets/foram_res_indicspecies.csv", row.names=FALSE)


#--------Indicator species analysis all forams-------------------------
df <- read.csv("~downloads/Datasets/merged_forams_data.csv")
tsm_chl <- read.csv("~downloads/Datasets/filtered_tsm_chl_mean_islands.csv")

df$depth_water <- df$depth..below.the.local.surface..in.meters.

df <- merge(df, tsm_chl[,c(1,4,5)], by = "island", all = TRUE)

#rarefy to the sample (none negative) with the lowest read number: 1924)
dfx <- df %>% 
  group_by(field.nmbr.) %>% 
  mutate(N=reads/sum(reads)*1924)

dfx1 <- dfx

df2 <- dfx1 %>% 
  group_by(field.nmbr.,ESVs, species,dist_coast_km, island, mean_chl_filtered,mean_tsm_filtered, depth_water, substrate) %>% 
  summarise(N_mean=mean(N))

df2$reef <- NA
df2$reef[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- "flat"
df2$reef[df2$depth_water %not% c("RF1", "RF2", "RF3")] <- "slope"
df2$reef[df2$depth_water %in% c("negative")] <- NA
df2$depth_water[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- 1

df2$turbidity <- NA
df2$turbidity[df2$mean_tsm_filtered > 1.5] <- "high turbidity"
df2$turbidity[df2$mean_tsm_filtered < 0.5] <- "low turbidity"
df2$turbidity[is.na(df2$turbidity)] <- "moderate turbidity"

df2$substrate[df2$substrate == "rubble + seagrass + algae"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "rubble + mud"] <- "rubble + sand"
df2$substrate[df2$substrate == "rubble + algae + seagrass"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "thick layer algae on rubble"] <- "rubble + algae"
df2$substrate[df2$substrate == "seagrass + algae + rubble + sand"] <- "rubble + seagrass"

df2$morpho_species <- NA

df2$morpho_species[df2$species == "Alveolinella_quoyi_Spermonde"] <- "Alveolinella quoyi"
df2$morpho_species[df2$species == "Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS_Spermonde"] <- "Amphisorus sp."
df2$morpho_species[df2$species == "Amphistegina_lessonii"] <- "Amphistegina lessonii"
df2$morpho_species[df2$species == "Amphistegina_lessonii_5174"] <- "Amphistegina lessonii"
df2$morpho_species[df2$species == "Amphistegina_lobifera"] <- "Amphistegina lobifera"
df2$morpho_species[df2$species == "Amphistegina_papillosa_Rik_9440_consensus0_1051_SaudiArabia"] <- "Amphistegina papillosa"
df2$morpho_species[df2$species == "Amphistegina_papillosa_Spermonde_3741"] <- "Amphistegina papillosa"
df2$morpho_species[df2$species == "Amphistegina_radiata_Spermonde"] <- "Amphistegina radiata"
df2$morpho_species[df2$species == "Baculogypsinoides_spinosus_8078_8079_8080"] <- "Baculogypsinoides spinosus"
df2$morpho_species[df2$species == "Borelis_schlumbergeri_3737_Spermonde"] <- "Borelis schlumbergeri"
df2$morpho_species[df2$species == "Calcarina_hispida_&_Calcarina_sp_5163_Spermonde"] <- "Calcarina hispida"
df2$morpho_species[df2$species == "Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._5247_Spermonde"] <- "Calcarina spengleri"
df2$morpho_species[df2$species == "Calcarina_sp._5164_Spermonde"] <- "Calcarina sp."
df2$morpho_species[df2$species == "Calcarina_sp._5165_Spermonde"] <- "Calcarina sp."                           
df2$morpho_species[df2$species == "Calcarina_sp._Spermonde"] <- "Calcarina sp."                  
df2$morpho_species[df2$species == "Calcarina_spengleri_Rik_9352_consensus0_1732_Spermonde"] <- "Calcarina spengleri"   
df2$morpho_species[df2$species == "Cyrea_szymborska_17247"] <- "Z_other_foraminifera"   
df2$morpho_species[df2$species == "Eggerelloides_scaber_12301"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Elphidium_Rik_9327_9337_9350_Spermonde"] <- "Elphidium sp."
df2$morpho_species[df2$species == "Encrusting_Rik_9441_consensus0_977"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Glabratellina_sp._17921_&_Planoglabratella_opercularis_18053_&_Planoglabratella_sp"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Globocassidulina_biora_17225"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Haynesina_germanica_18209"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Heterostegina_depressa_sp2_Spermonde"] <- "Heterostegina depressa"
df2$morpho_species[df2$species == "Heterostegina_depressa_sp1_Spermonde"] <- "Heterostegina depressa"
df2$morpho_species[df2$species == "Miliolid_Rik_9400_9405_Israel"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Milliolid_red_Rik_9370_9371_9372_Israel"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Murrayinella_globosa_18043"] <- "Z_other_foraminifera"                                                         
df2$morpho_species[df2$species == "Neorotalia_calcar_Spermonde"] <- "Neorotalia calcar"
df2$morpho_species[df2$species == "Neorotalia_gaimardi_&_Baculogypsina_sphaerulata_Spermonde"] <- "Neorotalia gaimardi"
df2$morpho_species[df2$species == "Nummulites_venosus"] <- "Nummulites venosus"
df2$morpho_species[df2$species == "Operculina_ammonoides"] <- "Operculina ammonoides"
df2$morpho_species[df2$species == "Parasorites_sp"] <- "Parasorites sp."
df2$morpho_species[df2$species == "Peneroplis_sp1_Spermonde"] <- "Peneroplis sp."
df2$morpho_species[df2$species == "Peneroplis_sp2_&_Peneroplis_pertusus_5117_&_Dendritina_ambigua_Spermonde"] <- "Peneroplis sp."
df2$morpho_species[df2$species == "Planorbulinella_sp._17893"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Psammophaga_sp._19260_&_Psammophaga_sp._19296_&_Psammophaga_sp2"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Pyrgo_Rik_9347_consensus0_1557"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Sorites_sp1_Spermonde"] <- "Sorites sp."
df2$morpho_species[df2$species == "Sorites_sp2_Spermonde"] <- "Sorites sp."
df2$morpho_species[df2$species == "Stainforthia_sp._17321"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Textularid_Rik_9397_consensus4_105"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Textularid_Rik_9398_consensus128_73"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Trifarina_earlandi_&_Uvigerina_bifurcata"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Vertebralina_Rik_9325_consensus0_1751"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Z_not_foraminifera"] <- "Z_not_foraminifera" 
df2$morpho_species[df2$species == "Z_other_foraminifera"] <- "Z_other_foraminifera" 

unique(df2$morpho_species)

df2a <- df2[df2$island != "negative",]
df2a <- df2a[!is.na(df2a$island),]
df2a$depth_water <- as.numeric(df2a$depth_water)


#_______________________________prep dataset 5% of the samples

df_NMDS_vegan <- dcast(df2a, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[is.na(df_NMDS_vegan)] <- 0

percent_df <- decostand(df_NMDS_vegan, "total")
pa_df <- df_NMDS_vegan
pa_df[pa_df > 0] <- 1

colSums(percent_df, na.rm = TRUE)
esv_5_sample <- colnames(pa_df[colSums(pa_df)/109 * 100 >= 5]) #in at least 5% of the samples

df2a <- df2a[df2a$ESVs %in% esv_5_sample,]

slope <- subset(df2a, reef == "slope")
flat <- subset(df2a, reef == "flat")
flat <- flat[flat$field.nmbr. != "UPG-EG22-053",]

#"unmelt" the dataframe with reshape and
#to isolate values (N), Species code (x) and Sample codes (y) with col and row names
slope_indicator <- dcast(slope, field.nmbr. + island + substrate + turbidity ~ ESVs, 
                         fun.aggregate = sum, value.var = "N_mean")
slope_indicator[is.na(slope_indicator)] <- 0
slope_species <- slope_indicator[,-c(2:4)]
rownames(slope_species) <- slope_species$Sample_Code
slope_species <- slope_species[,-1]
turbidity_slope <- slope_indicator$turbidity

flat_indicator <- dcast(flat, field.nmbr. + island + substrate + turbidity ~ ESVs, 
                        fun.aggregate = sum, value.var = "N_mean")
flat_indicator[is.na(flat_indicator)] <- 0
flat_species <- flat_indicator[,-c(2:4)]
rownames(flat_species) <- flat_species$Sample_Code
flat_species <- flat_species[,-1]
turbidity_flat <- flat_indicator$turbidity

#________________________________________Indicator species analysis using IndVal.g - ABUNDANCE

#Indicator species analysis for turbidity levels
indval_s<- multipatt(slope_species, turbidity_slope, control = how(nperm=2000))
summary(indval_s)
ESVs_slope <- indval_s$sign[indval_s$sign$p.value < 0.01 & !is.na(indval_s$sign$p.value),]
ESVs_slope$reef <- "slope"
ESVs_slope$ESVs <- rownames(ESVs_slope)

indval_f <- multipatt(flat_species, turbidity_flat, control = how(nperm=2000))
summary(indval_f)
ESVs_flat <- indval_f$sign[indval_f$sign$p.value < 0.01 & !is.na(indval_f$sign$p.value),]
ESVs_flat$reef <- "flat"
ESVs_flat$ESVs <- rownames(ESVs_flat)


allforam_res <- rbind(ESVs_slope, ESVs_flat)

allforam_res_indicspecies <- merge(unique(df[c("ESVs", "family", "species")]), allforam_res, by = "ESVs")

write.csv(allforam_res_indicspecies, "~downloads/Datasets/allforam_res_indicspecies.csv", row.names=FALSE)


#--------Indicator species analysis DIATOM-------------------------
df <- read.csv("~downloads/Datasets/merged_diat_data.csv")

df$depth_water <- df$depth..below.the.local.surface..in.meters.

tsm_chl <- read.csv("~downloads/Datasets/filtered_tsm_chl_mean_islands.csv")

df <- merge(df, tsm_chl[,c(1,4,5)], by = "island", all = TRUE)


#rarefy to the sample (none negative) with the lowest read number)
test <- df %>% 
  group_by(field.nmbr.,Substrate.type) %>% 
  summarise(N_sum=sum(reads))
test1 <- test[test$N_sum > 1000,] #lowest read number for sample > 1000 = 1654
df1 <- merge(df, test1[,1:3], by = c("field.nmbr.", "Substrate.type"))

dfx <- df1 %>% 
  group_by(field.nmbr.) %>% 
  mutate(N=reads/sum(reads)*1654)

dfx1 <- dfx[dfx$Identity.percentage > 90, ] #assigned to species level


df2 <- dfx1 %>% 
  group_by(field.nmbr., ESVs, Species,dist_coast_km, island, mean_chl_filtered,mean_tsm_filtered, depth_water, substrate, Substrate.type) %>% 
  summarise(N_mean=mean(N))


df2$reef <- NA
df2$reef[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- "flat"
df2$reef[df2$depth_water %not% c("RF1", "RF2", "RF3")] <- "slope"
df2$reef[is.na(df2$depth_water)] <- NA
df2$depth_water[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- 1
df2$depth_water <- as.numeric(df2$depth_water)
df2$depth_group <- NA
df2$depth_group[df2$depth_water < 10 & df2$reef == "slope"] <- "<10m"
df2$depth_group[df2$depth_water > 9 & df2$depth_water < 21] <- "10-20m"
df2$depth_group[df2$depth_water > 20 & df2$reef == "slope"] <- ">20m"

df2$turbidity <- NA
df2$turbidity[df2$mean_tsm_filtered > 1.5] <- "high turbidity"
df2$turbidity[df2$mean_tsm_filtered < 0.5] <- "low turbidity"
df2$turbidity[is.na(df2$turbidity)] <- "moderate turbidity"

df2$substrate[df2$substrate == "rubble + seagrass + algae"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "rubble + mud"] <- "rubble + sand"
df2$substrate[df2$substrate == "rubble + algae + seagrass"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "thick layer algae on rubble"] <- "rubble + algae"
df2$substrate[df2$substrate == "seagrass + algae + rubble + sand"] <- "rubble + seagrass"

df2a <- df2[df2$island != "negative",]
df2a <- df2a[!is.na(df2a$island),]
df2a <- df2a[ df2a$Substrate.type == "substrate",]
df2a$depth_water <- as.numeric(df2a$depth_water)


#_______________________________prep dataset 5% of the samples

df_NMDS_vegan <- dcast(df2a, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[is.na(df_NMDS_vegan)] <- 0

percent_df <- decostand(df_NMDS_vegan, "total")
pa_df <- df_NMDS_vegan
pa_df[pa_df > 0] <- 1

colSums(percent_df, na.rm = TRUE)
esv_5_sample <- colnames(pa_df[colSums(pa_df)/109 * 100 >= 5]) #in at least 5% of the samples

df2a <- df2a[df2a$ESVs %in% esv_5_sample,]

slope <- subset(df2a, reef == "slope")
flat <- subset(df2a, reef == "flat")

#"unmelt" the dataframe with reshape and
#to isolate values (N), Species code (x) and Sample codes (y) with col and row names
slope_indicator <- dcast(slope, field.nmbr. + island + substrate + turbidity ~ ESVs, 
                         fun.aggregate = sum, value.var = "N_mean")
slope_indicator[is.na(slope_indicator)] <- 0
slope_species <- slope_indicator[,-c(2:4)]
rownames(slope_species) <- slope_species$Sample_Code
slope_species <- slope_species[,-1]
turbidity_slope <- slope_indicator$turbidity

flat_indicator <- dcast(flat, field.nmbr. + island + substrate + turbidity ~ ESVs, 
                        fun.aggregate = sum, value.var = "N_mean")
flat_indicator[is.na(flat_indicator)] <- 0
flat_species <- flat_indicator[,-c(2:4)]
rownames(flat_species) <- flat_species$Sample_Code
flat_species <- flat_species[,-1]
turbidity_flat <- flat_indicator$turbidity

#________________________________________Indicator species analysis using IndVal.g - ABUNDANCE

#Indicator species analysis for turbidity levels
indval_s<- multipatt(slope_species, turbidity_slope, control = how(nperm=2000))
summary(indval_s)
ESVs_slope <- indval_s$sign[indval_s$sign$p.value < 0.01 & !is.na(indval_s$sign$p.value),]
ESVs_slope$reef <- "slope"
ESVs_slope$ESVs <- rownames(ESVs_slope)

indval_f <- multipatt(flat_species, turbidity_flat, control = how(nperm=2000))
summary(indval_f)
ESVs_flat <- indval_f$sign[indval_f$sign$p.value < 0.01 & !is.na(indval_f$sign$p.value),]
ESVs_flat$reef <- "flat"
ESVs_flat$ESVs <- rownames(ESVs_flat)


diat_res <- rbind(ESVs_slope, ESVs_flat)

diat_res_indicspecies <- merge(unique(df[c("ESVs", "Genus", "Species", "Identity.percentage")]),diat_res, by = "ESVs")

write.csv(diat_res_indicspecies,"~downloads/Datasets/diat_res_indicspecies.csv", row.names=FALSE)

#--------Indicator species analysis PROK-------------------------
df <- read.csv("~downloads/Datasets/merged_prok_data.csv")

df$depth_water <- df$depth..below.the.local.surface..in.meters.

tsm_chl <- read.csv("~downloads/Datasets/filtered_tsm_chl_mean_islands.csv")

df <- merge(df, tsm_chl[,c(1,4,5)], by = "island", all = TRUE)


#rarefy to the sample (none negative) with the lowest read number)
test <- df %>% 
  group_by(field.nmbr.,Substrate.type) %>% 
  summarise(N_sum=sum(reads))
test1 <- test[test$N_sum > 1000,] #lowest read number for sample > 1000 = 1905
df1 <- merge(df, test1[,1:3], by = c("field.nmbr.", "Substrate.type"))

dfx <- df1 %>% 
  group_by(field.nmbr.) %>% 
  mutate(N=reads/sum(reads)*1905)

dfx1 <- dfx[dfx$Identity.percentage > 90, ] #assigned to species level


df2 <- dfx1 %>% 
  group_by(field.nmbr., ESVs, Species,dist_coast_km, island, mean_chl_filtered,mean_tsm_filtered, depth_water, substrate, Substrate.type) %>% 
  summarise(N_mean=mean(N))


df2$reef <- NA
df2$reef[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- "flat"
df2$reef[df2$depth_water %not% c("RF1", "RF2", "RF3")] <- "slope"
df2$reef[is.na(df2$depth_water)] <- NA
df2$depth_water[df2$depth_water %in% c("RF1", "RF2", "RF3")] <- 1
df2$depth_water <- as.numeric(df2$depth_water)
df2$depth_group <- NA
df2$depth_group[df2$depth_water < 10 & df2$reef == "slope"] <- "<10m"
df2$depth_group[df2$depth_water > 9 & df2$depth_water < 21] <- "10-20m"
df2$depth_group[df2$depth_water > 20 & df2$reef == "slope"] <- ">20m"

df2$turbidity <- NA
df2$turbidity[df2$mean_tsm_filtered > 1.5] <- "high turbidity"
df2$turbidity[df2$mean_tsm_filtered < 0.5] <- "low turbidity"
df2$turbidity[is.na(df2$turbidity)] <- "moderate turbidity"

df2$substrate[df2$substrate == "rubble + seagrass + algae"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "rubble + mud"] <- "rubble + sand"
df2$substrate[df2$substrate == "rubble + algae + seagrass"] <- "rubble + seagrass"
df2$substrate[df2$substrate == "thick layer algae on rubble"] <- "rubble + algae"
df2$substrate[df2$substrate == "seagrass + algae + rubble + sand"] <- "rubble + seagrass"

df2a <- df2[df2$island != "negative",]
df2a <- df2a[!is.na(df2a$island),]
df2a <- df2a[ df2a$Substrate.type == "substrate",]
df2a$depth_water <- as.numeric(df2a$depth_water)


#_______________________________prep dataset 5% of the samples

df_NMDS_vegan <- dcast(df2a, field.nmbr. ~ ESVs, fun.aggregate = sum, value.var = "N_mean")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$field.nmbr.
df_NMDS_vegan <- df_NMDS_vegan[,-1]
df_NMDS_vegan[is.na(df_NMDS_vegan)] <- 0

percent_df <- decostand(df_NMDS_vegan, "total")
pa_df <- df_NMDS_vegan
pa_df[pa_df > 0] <- 1

colSums(percent_df, na.rm = TRUE)
esv_5_sample <- colnames(pa_df[colSums(pa_df)/109 * 100 >= 5]) #in at least 5% of the samples

df2a <- df2a[df2a$ESVs %in% esv_5_sample,]

slope <- subset(df2a, reef == "slope")
flat <- subset(df2a, reef == "flat")

#"unmelt" the dataframe with reshape and
#to isolate values (N), Species code (x) and Sample codes (y) with col and row names
slope_indicator <- dcast(slope, field.nmbr. + island + substrate + turbidity ~ ESVs, 
                         fun.aggregate = sum, value.var = "N_mean")
slope_indicator[is.na(slope_indicator)] <- 0
slope_species <- slope_indicator[,-c(2:4)]
rownames(slope_species) <- slope_species$Sample_Code
slope_species <- slope_species[,-1]
turbidity_slope <- slope_indicator$turbidity

flat_indicator <- dcast(flat, field.nmbr. + island + substrate + turbidity ~ ESVs, 
                        fun.aggregate = sum, value.var = "N_mean")
flat_indicator[is.na(flat_indicator)] <- 0
flat_species <- flat_indicator[,-c(2:4)]
rownames(flat_species) <- flat_species$Sample_Code
flat_species <- flat_species[,-1]
turbidity_flat <- flat_indicator$turbidity

#________________________________________Indicator species analysis using IndVal.g - ABUNDANCE

#Indicator species analysis for turbidity levels
indval_s<- multipatt(slope_species, turbidity_slope, control = how(nperm=2000))
summary(indval_s)
ESVs_slope <- indval_s$sign[indval_s$sign$p.value < 0.01 & !is.na(indval_s$sign$p.value),]
ESVs_slope$reef <- "slope"
ESVs_slope$ESVs <- rownames(ESVs_slope)

indval_f <- multipatt(flat_species, turbidity_flat, control = how(nperm=2000))
summary(indval_f)
ESVs_flat <- indval_f$sign[indval_f$sign$p.value < 0.01 & !is.na(indval_f$sign$p.value),]
ESVs_flat$reef <- "flat"
ESVs_flat$ESVs <- rownames(ESVs_flat)


prok_res <- rbind(ESVs_slope, ESVs_flat)

prok_res_indicspecies <- merge(unique(df[c("ESVs", "Family","Genus", "Species", "Identity.percentage")]),prok_res, by = "ESVs")

write.csv(prok_res_indicspecies,"~downloads/Datasets/prok_res_indicspecies.csv", row.names=FALSE)

#--------Supp. Fig S3 --- FORAMS species LBF - compositionnal turnover (zeta diversity)-----------------------

df <- read.csv("~downloads/Datasets/foram_working_dataset_turnover.csv")



#_________________reef flat__________________________
df2b <- df[df$reef == "flat",] #select only reef flat samples

#put ESVs into columns while keeping metadata for samples
df3 <- reshape2::dcast(df2b, field.nmbr. + depth_water + island + dist_coast_km + mean_chl_filtered + mean_tsm_filtered + turbidity + reef + substrate ~ ESVs, 
                       fun.aggregate = sum, value.var = "N_mean")

#change NAs for 0
df3[df3 == "NaN"] <- 0

#make sure all variables are factor or numeric
df3[,1] <- as.factor(df3[,1])
df3[,3] <- as.factor(df3[,3])
df3[,4] <- as.factor(df3[,4])
df3[,7] <- as.factor(df3[,7])
df3[,8] <- as.factor(df3[,8])
df3[,9] <- as.factor(df3[,9])

#create spec data site-ESVs (row-column) with abundance
data_spec <- df3[,10:ncol(df3)]

#transform spec data to presence-absence
data_spec_pa_flat <- decostand(data_spec, "pa")


### REEF FLAT PLOT

# Update num_sites and max_order
num_sites <- nrow(data_spec_pa_flat)
max_order <- ifelse(num_sites > 1, num_sites - 1, 1)


# Calculate zeta diversity decay parameters and plots
zetaDecay.ex_asvs_HW1_site <- Zeta.decline.ex(data_spec_pa_flat, orders = 1:max_order, plot = TRUE)

zetaDecay.ex_asvs_HW1_site[5]$ratio <- c(zetaDecay.ex_asvs_HW1_site[5]$ratio, NA)
zeta_ratios <- zetaDecay.ex_asvs_HW1_site[5]$ratio

# Create a tibble from the zeta_ratios vector
zeta_HW1_site <- tibble(
  zeta.order = 1:length(zeta_ratios),
  ratio = zeta_ratios,
  Region = "HW1",
  Level = "Sites"
)

zeta_HW1_site



#____________________reef slope______________
df2b <- df[df$reef == "slope",]

df3 <- reshape2::dcast(df2b, field.nmbr. + depth_water + island + dist_coast_km + mean_chl_filtered + mean_tsm_filtered + turbidity + reef + substrate ~ ESVs, 
                       fun.aggregate = sum, value.var = "N_mean")

df3[df3 == "NaN"] <- 0

df3[,1] <- as.factor(df3[,1])
df3[,3] <- as.factor(df3[,3])
df3[,4] <- as.factor(df3[,4])
df3[,7] <- as.factor(df3[,7])
df3[,8] <- as.factor(df3[,8])
df3[,9] <- as.factor(df3[,9])

data_spec <- df3[,10:ncol(df3)]
data_spec_pa_slope <- decostand(data_spec, "pa")

### REEF SLOPE PLOT

# Update num_sites and max_order
num_sites <- nrow(data_spec_pa_slope)
max_order <- ifelse(num_sites > 1, num_sites - 1, 1)


# Calculate zeta diversity decay parameters and plots
zetaDecay.ex_asvs_HW2_site <- Zeta.decline.ex(data_spec_pa_slope, orders = 1:max_order, plot = TRUE)

zetaDecay.ex_asvs_HW2_site[5]$ratio <- c(zetaDecay.ex_asvs_HW2_site[5]$ratio, NA)
zeta_ratios <- zetaDecay.ex_asvs_HW2_site[5]$ratio

# Create a tibble from the zeta_ratios vector
zeta_HW2_site <- tibble(
  zeta.order = 1:length(zeta_ratios),
  ratio = zeta_ratios,
  Region = "HW1",
  Level = "Sites"
)

zeta_HW2_site


### Zeta Decay Figure

zetaDecay.ex_asvs_HW1_site[1:4] %>% as_tibble() %>% mutate(Region = "Flat")-> zdx_HW1_tib
zetaDecay.ex_asvs_HW2_site[1:4] %>% as_tibble() %>% mutate(Region = "Slope")-> zdx_HW2_tib

bind_rows(zdx_HW1_tib,zdx_HW2_tib) -> zdx_regions


#________________PLOTs______________________

# Define the colors for each region
col.8 <- c("#C7AE4C" ,"#008FA7" )
names(col.8) <- c("Flat", "Slope")

##PLOT##
####
zdv_decay_sites <- zdx_regions %>%
  ggplot(aes(x = zeta.order, y = zeta.val, group = Region, color = Region, fill = Region)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = zeta.val - zeta.val.sd, ymax = zeta.val + zeta.val.sd),
              alpha = 0.1, colour = NA) +
  scale_fill_manual(values = col.8, limits = c("Flat", "Slope")) +
  scale_colour_manual(values = col.8, limits = c("Flat", "Slope")) +
  labs(x = "Zeta Order - Sites", y = "Zeta Diversity") +
  theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    legend.position = "right"
  ) +
  scale_x_continuous(breaks = seq(0, 75, 5), limits = c(0, 75)) +
  scale_y_continuous(breaks = seq(0, 15, 2), limits = c(0, 15))

zdv_decay_sites

### Species Retention Figure

zetaDecay.ex_asvs_HW1_site[5] %>% as_tibble() %>% mutate(Region = "Flat",zeta.order = 1:n()) -> ret_HW1_tib
zetaDecay.ex_asvs_HW2_site[5] %>% as_tibble() %>% mutate(Region = "Slope",zeta.order = 1:n()) -> ret_HW2_tib

bind_rows(ret_HW1_tib,ret_HW2_tib) -> zdx_regions_ret_sites

spec_ret_sites <- zdx_regions_ret_sites %>% 
  ggplot(aes(x = zeta.order, y = ratio, group = Region, color = Region)) +
  geom_point() +
  geom_line() +
  scale_fill_manual(values = col.8) +
  scale_colour_manual(values = col.8) +
  ylab("Zeta Ratio") +
  xlab("Zeta Order - Sites") +
  ylim(0, 1) +
  scale_x_continuous(breaks = seq(0, 75, 5), limits = c(0, 75)) +
  theme(
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    legend.position = "right" # Adjust legend position as needed
  )

spec_ret_sites


#--------Supp. Fig S3 --- FORAMS all - compositionnal turnover (zeta diversity)-----------------------

df <- read.csv("~downloads/Datasets/Allforam_working_dataset_turnover.csv")



#_________________reef flat__________________________
df2b <- df[df$reef == "flat",] #select only reef flat samples

#put ESVs into columns while keeping metadata for samples
df3 <- reshape2::dcast(df2b, field.nmbr. + depth_water + island + dist_coast_km + mean_chl_filtered + mean_tsm_filtered + turbidity + reef + substrate ~ ESVs, 
                       fun.aggregate = sum, value.var = "N_mean")

#change NAs for 0
df3[df3 == "NaN"] <- 0

#make sure all variables are factor or numeric
df3[,1] <- as.factor(df3[,1])
df3[,3] <- as.factor(df3[,3])
df3[,4] <- as.factor(df3[,4])
df3[,7] <- as.factor(df3[,7])
df3[,8] <- as.factor(df3[,8])
df3[,9] <- as.factor(df3[,9])

#create spec data site-ESVs (row-column) with abundance
data_spec <- df3[,10:ncol(df3)]

#transform spec data to presence-absence
data_spec_pa_flat <- decostand(data_spec, "pa")


### REEF FLAT PLOT

# Update num_sites and max_order
num_sites <- nrow(data_spec_pa_flat)
max_order <- ifelse(num_sites > 1, num_sites - 1, 1)


# Calculate zeta diversity decay parameters and plots
zetaDecay.ex_asvs_HW1_site <- Zeta.decline.ex(data_spec_pa_flat, orders = 1:max_order, plot = TRUE)

zetaDecay.ex_asvs_HW1_site[5]$ratio <- c(zetaDecay.ex_asvs_HW1_site[5]$ratio, NA)
zeta_ratios <- zetaDecay.ex_asvs_HW1_site[5]$ratio

# Create a tibble from the zeta_ratios vector
zeta_HW1_site <- tibble(
  zeta.order = 1:length(zeta_ratios),
  ratio = zeta_ratios,
  Region = "HW1",
  Level = "Sites"
)

zeta_HW1_site



#____________________reef slope______________
df2b <- df[df$reef == "slope",]

df3 <- reshape2::dcast(df2b, field.nmbr. + depth_water + island + dist_coast_km + mean_chl_filtered + mean_tsm_filtered + turbidity + reef + substrate ~ ESVs, 
                       fun.aggregate = sum, value.var = "N_mean")

df3[df3 == "NaN"] <- 0

df3[,1] <- as.factor(df3[,1])
df3[,3] <- as.factor(df3[,3])
df3[,4] <- as.factor(df3[,4])
df3[,7] <- as.factor(df3[,7])
df3[,8] <- as.factor(df3[,8])
df3[,9] <- as.factor(df3[,9])

data_spec <- df3[,10:ncol(df3)]
data_spec_pa_slope <- decostand(data_spec, "pa")

### REEF SLOPE PLOT

# Update num_sites and max_order
num_sites <- nrow(data_spec_pa_slope)
max_order <- ifelse(num_sites > 1, num_sites - 1, 1)


# Calculate zeta diversity decay parameters and plots
zetaDecay.ex_asvs_HW2_site <- Zeta.decline.ex(data_spec_pa_slope, orders = 1:max_order, plot = TRUE)

zetaDecay.ex_asvs_HW2_site[5]$ratio <- c(zetaDecay.ex_asvs_HW2_site[5]$ratio, NA)
zeta_ratios <- zetaDecay.ex_asvs_HW2_site[5]$ratio

# Create a tibble from the zeta_ratios vector
zeta_HW2_site <- tibble(
  zeta.order = 1:length(zeta_ratios),
  ratio = zeta_ratios,
  Region = "HW1",
  Level = "Sites"
)

zeta_HW2_site


### Zeta Decay Figure

zetaDecay.ex_asvs_HW1_site[1:4] %>% as_tibble() %>% mutate(Region = "Flat")-> zdx_HW1_tib
zetaDecay.ex_asvs_HW2_site[1:4] %>% as_tibble() %>% mutate(Region = "Slope")-> zdx_HW2_tib

bind_rows(zdx_HW1_tib,zdx_HW2_tib) -> zdx_regions


#________________PLOTs______________________

# Define the colors for each region
col.8 <- c("#C7AE4C" ,"#008FA7" )
names(col.8) <- c("Flat", "Slope")

##PLOT##
####
zdv_decay_sites <- zdx_regions %>%
  ggplot(aes(x = zeta.order, y = zeta.val, group = Region, color = Region, fill = Region)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = zeta.val - zeta.val.sd, ymax = zeta.val + zeta.val.sd),
              alpha = 0.1, colour = NA) +
  scale_fill_manual(values = col.8, limits = c("Flat", "Slope")) +
  scale_colour_manual(values = col.8, limits = c("Flat", "Slope")) +
  labs(x = "Zeta Order - Sites", y = "Zeta Diversity") +
  theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    legend.position = "right"
  ) +
  scale_x_continuous(breaks = seq(0, 75, 5), limits = c(0, 75)) +
  scale_y_continuous(breaks = seq(0, 15, 2), limits = c(0, 15))

zdv_decay_sites

### Species Retention Figure

zetaDecay.ex_asvs_HW1_site[5] %>% as_tibble() %>% mutate(Region = "Flat",zeta.order = 1:n()) -> ret_HW1_tib
zetaDecay.ex_asvs_HW2_site[5] %>% as_tibble() %>% mutate(Region = "Slope",zeta.order = 1:n()) -> ret_HW2_tib

bind_rows(ret_HW1_tib,ret_HW2_tib) -> zdx_regions_ret_sites

spec_ret_sites <- zdx_regions_ret_sites %>% 
  ggplot(aes(x = zeta.order, y = ratio, group = Region, color = Region)) +
  geom_point() +
  geom_line() +
  scale_fill_manual(values = col.8) +
  scale_colour_manual(values = col.8) +
  ylab("Zeta Ratio") +
  xlab("Zeta Order - Sites") +
  ylim(0, 1) +
  scale_x_continuous(breaks = seq(0, 75, 5), limits = c(0, 75)) +
  theme(
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    legend.position = "right" # Adjust legend position as needed
  )

spec_ret_sites


#--------Supp. Fig S3 --- DIATOMS - compositionnal turnover (zeta diversity)-----------------------

df <- read.csv("~downloads/Datasets/diat_working_dataset_turnover.csv")


#____reef flat - substrate_______
df2b <- df[df$reef == "flat" & df$Substrate.type == "substrate",]

df3 <- reshape2::dcast(df2b, field.nmbr. ~ ESVs, 
             fun.aggregate = sum, value.var = "N_mean")


#change NAs for 0
df3[df3 == "NaN"] <- 0

#make sure all variables are factor or numeric
df3[,1] <- as.factor(df3[,1])

#create spec data site-ESVs (row-column) with abundance
data_spec <- df3[,2:ncol(df3)]

#transform spec data to presence-absence
data_spec_pa_flat <- decostand(data_spec, "pa")


### REEF FLAT PLOT

# Update num_sites and max_order
num_sites <- nrow(data_spec_pa_flat)
max_order <- ifelse(num_sites > 1, num_sites - 1, 1)


# Calculate zeta diversity decay parameters and plots
zetaDecay.ex_asvs_HW1_site <- Zeta.decline.ex(data_spec_pa_flat, orders = 1:max_order, plot = TRUE)

zetaDecay.ex_asvs_HW1_site[5]$ratio <- c(zetaDecay.ex_asvs_HW1_site[5]$ratio, NA)
zeta_ratios <- zetaDecay.ex_asvs_HW1_site[5]$ratio

# Create a tibble from the zeta_ratios vector
zeta_HW1_site <- tibble(
  zeta.order = 1:length(zeta_ratios),
  ratio = zeta_ratios,
  Region = "HW1",
  Level = "Sites"
)

zeta_HW1_site



##____reef slope - substrate_______
df2b <- df[df$reef == "slope" & df$Substrate.type == "substrate",]

df3 <- reshape2::dcast(df2b, field.nmbr. ~ ESVs, 
                       fun.aggregate = sum, value.var = "N_mean")


#change NAs for 0
df3[df3 == "NaN"] <- 0

#make sure all variables are factor or numeric
df3[,1] <- as.factor(df3[,1])


data_spec <- df3[,2:ncol(df3)]
data_spec_pa_slope <- decostand(data_spec, "pa")

### REEF SLOPE PLOT

# Update num_sites and max_order
num_sites <- nrow(data_spec_pa_slope)
max_order <- ifelse(num_sites > 1, num_sites - 1, 1)


# Calculate zeta diversity decay parameters and plots
zetaDecay.ex_asvs_HW2_site <- Zeta.decline.ex(data_spec_pa_slope, orders = 1:max_order, plot = TRUE)

zetaDecay.ex_asvs_HW2_site[5]$ratio <- c(zetaDecay.ex_asvs_HW2_site[5]$ratio, NA)
zeta_ratios <- zetaDecay.ex_asvs_HW2_site[5]$ratio

# Create a tibble from the zeta_ratios vector
zeta_HW2_site <- tibble(
  zeta.order = 1:length(zeta_ratios),
  ratio = zeta_ratios,
  Region = "HW1",
  Level = "Sites"
)

zeta_HW2_site



#________________PLOTS______________________

### Zeta Decay Figure

zetaDecay.ex_asvs_HW1_site[1:4] %>% as_tibble() %>% mutate(Region = "Flat")-> zdx_HW1_tib
zetaDecay.ex_asvs_HW2_site[1:4] %>% as_tibble() %>% mutate(Region = "Slope")-> zdx_HW2_tib

bind_rows(zdx_HW1_tib,zdx_HW2_tib) -> zdx_regions

# Define the colors for each region
col.8 <- c("#C7AE4C" ,"#008FA7" )
names(col.8) <- c("Flat", "Slope")

##PLOT##
####
zdv_decay_sites <- zdx_regions %>%
  ggplot(aes(x = zeta.order, y = zeta.val, group = Region, color = Region, fill = Region)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = zeta.val - zeta.val.sd, ymax = zeta.val + zeta.val.sd),
              alpha = 0.1, colour = NA) +
  scale_fill_manual(values = col.8, limits = c("Flat", "Slope")) +
  scale_colour_manual(values = col.8, limits = c("Flat", "Slope")) +
  labs(x = "Zeta Order - Sites", y = "Zeta Diversity") +
  theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    legend.position = "right"
  ) +
  scale_x_continuous(breaks = seq(0, 75, 5), limits = c(0, 75))

zdv_decay_sites

### Species Retention Figure

zetaDecay.ex_asvs_HW1_site[5] %>% as_tibble() %>% mutate(Region = "Flat",zeta.order = 1:n()) -> ret_HW1_tib
zetaDecay.ex_asvs_HW2_site[5] %>% as_tibble() %>% mutate(Region = "Slope",zeta.order = 1:n()) -> ret_HW2_tib

bind_rows(ret_HW1_tib,ret_HW2_tib) -> zdx_regions_ret_sites

spec_ret_sites <- zdx_regions_ret_sites %>% 
  ggplot(aes(x = zeta.order, y = ratio, group = Region, color = Region)) +
  geom_point() +
  geom_line() +
  scale_fill_manual(values = col.8) +
  scale_colour_manual(values = col.8) +
  ylab("Zeta Ratio") +
  xlab("Zeta Order - Sites") +
  ylim(0, 1) +
  scale_x_continuous(breaks = seq(0, 75, 5), limits = c(0, 75)) +
  theme(
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    legend.position = "right" # Adjust legend position as needed
  )

spec_ret_sites

#--------Supp. Fig S3 --- PROKS - compositionnal turnover (zeta diversity)-----------------------

df <- read.csv("~downloads/Datasets/prok_working_dataset_turnover.csv")


#____reef flat - substrate_______
df2b <- df[df$reef == "flat" & df$Substrate.type == "substrate",]

df3 <- reshape2::dcast(df2b, field.nmbr. ~ ESVs, 
                       fun.aggregate = sum, value.var = "N_mean")


#change NAs for 0
df3[df3 == "NaN"] <- 0

#make sure all variables are factor or numeric
df3[,1] <- as.factor(df3[,1])

#create spec data site-ESVs (row-column) with abundance
data_spec <- df3[,2:ncol(df3)]

#transform spec data to presence-absence
data_spec_pa_flat <- decostand(data_spec, "pa")


### REEF FLAT PLOT

# Update num_sites and max_order
num_sites <- nrow(data_spec_pa_flat)
max_order <- ifelse(num_sites > 1, num_sites - 1, 1)


# Calculate zeta diversity decay parameters and plots
zetaDecay.ex_asvs_HW1_site <- Zeta.decline.ex(data_spec_pa_flat, orders = 1:max_order, plot = TRUE)

zetaDecay.ex_asvs_HW1_site[5]$ratio <- c(zetaDecay.ex_asvs_HW1_site[5]$ratio, NA)
zeta_ratios <- zetaDecay.ex_asvs_HW1_site[5]$ratio

# Create a tibble from the zeta_ratios vector
zeta_HW1_site <- tibble(
  zeta.order = 1:length(zeta_ratios),
  ratio = zeta_ratios,
  Region = "HW1",
  Level = "Sites"
)

zeta_HW1_site



##____reef slope - substrate_______
df2b <- df[df$reef == "slope" & df$Substrate.type == "substrate",]

df3 <- reshape2::dcast(df2b, field.nmbr. ~ ESVs, 
                       fun.aggregate = sum, value.var = "N_mean")


#change NAs for 0
df3[df3 == "NaN"] <- 0

#make sure all variables are factor or numeric
df3[,1] <- as.factor(df3[,1])


data_spec <- df3[,2:ncol(df3)]
data_spec_pa_slope <- decostand(data_spec, "pa")

### REEF SLOPE PLOT

# Update num_sites and max_order
num_sites <- nrow(data_spec_pa_slope)
max_order <- ifelse(num_sites > 1, num_sites - 1, 1)


# Calculate zeta diversity decay parameters and plots
zetaDecay.ex_asvs_HW2_site <- Zeta.decline.ex(data_spec_pa_slope, orders = 1:max_order, plot = TRUE)

zetaDecay.ex_asvs_HW2_site[5]$ratio <- c(zetaDecay.ex_asvs_HW2_site[5]$ratio, NA)
zeta_ratios <- zetaDecay.ex_asvs_HW2_site[5]$ratio

# Create a tibble from the zeta_ratios vector
zeta_HW2_site <- tibble(
  zeta.order = 1:length(zeta_ratios),
  ratio = zeta_ratios,
  Region = "HW1",
  Level = "Sites"
)

zeta_HW2_site



#________________PLOTS______________________

### Zeta Decay Figure

zetaDecay.ex_asvs_HW1_site[1:4] %>% as_tibble() %>% mutate(Region = "Flat")-> zdx_HW1_tib
zetaDecay.ex_asvs_HW2_site[1:4] %>% as_tibble() %>% mutate(Region = "Slope")-> zdx_HW2_tib

bind_rows(zdx_HW1_tib,zdx_HW2_tib) -> zdx_regions

# Define the colors for each region
col.8 <- c("#C7AE4C" ,"#008FA7" )
names(col.8) <- c("Flat", "Slope")

##PLOT##
####
zdv_decay_sites <- zdx_regions %>%
  ggplot(aes(x = zeta.order, y = zeta.val, group = Region, color = Region, fill = Region)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = zeta.val - zeta.val.sd, ymax = zeta.val + zeta.val.sd),
              alpha = 0.1, colour = NA) +
  scale_fill_manual(values = col.8, limits = c("Flat", "Slope")) +
  scale_colour_manual(values = col.8, limits = c("Flat", "Slope")) +
  labs(x = "Zeta Order - Sites", y = "Zeta Diversity") +
  theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    legend.position = "right"
  ) +
  scale_x_continuous(breaks = seq(0, 75, 5), limits = c(0, 75))

zdv_decay_sites

### Species Retention Figure

zetaDecay.ex_asvs_HW1_site[5] %>% as_tibble() %>% mutate(Region = "Flat",zeta.order = 1:n()) -> ret_HW1_tib
zetaDecay.ex_asvs_HW2_site[5] %>% as_tibble() %>% mutate(Region = "Slope",zeta.order = 1:n()) -> ret_HW2_tib

bind_rows(ret_HW1_tib,ret_HW2_tib) -> zdx_regions_ret_sites

spec_ret_sites <- zdx_regions_ret_sites %>% 
  ggplot(aes(x = zeta.order, y = ratio, group = Region, color = Region)) +
  geom_point() +
  geom_line() +
  scale_fill_manual(values = col.8) +
  scale_colour_manual(values = col.8) +
  ylab("Zeta Ratio") +
  xlab("Zeta Order - Sites") +
  ylim(0, 1) +
  scale_x_continuous(breaks = seq(0, 75, 5), limits = c(0, 75)) +
  theme(
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    legend.position = "right" # Adjust legend position as needed
  )

spec_ret_sites



#--------Fig 2a, 3a 4a, Supp. Fig. S2_S4 ---- ESV counts and alpha diversity-----------------------

#_______________FORAMS______________________
df <- read.csv("~downloads/Datasets/filtered_forams_esv_taxo_meta_species_replifilt.csv")

df$depth_water <- df$depth..below.the.local.surface..in.meters.

#rarefy to the sample (none negative) with the lowest read number: 1924)
dfx <- df %>% 
  group_by(field.nmbr.,replicate) %>% 
  mutate(N=reads/sum(reads)*1924)

dfx <- dfx %>% 
  group_by(field.nmbr.,replicate) %>% 
  mutate(percent_N=reads/sum(reads)*100)

df2 <- dfx %>% 
  group_by(field.nmbr., ESVs, species,dist_coast_km, island, depth_water, substrate, group, family) %>% 
  summarise(N_mean=mean(N), percent_N_mean = mean(percent_N))

df2$counts <- 1

df2$morpho_species <- NA

df2$morpho_species[df2$species == "Alveolinella_quoyi_Spermonde"] <- "Alveolinella quoyi"
df2$morpho_species[df2$species == "Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS_Spermonde"] <- "Amphisorus sp."
df2$morpho_species[df2$species == "Amphistegina_lessonii"] <- "Amphistegina lessonii"
df2$morpho_species[df2$species == "Amphistegina_lessonii_5174"] <- "Amphistegina lessonii"
df2$morpho_species[df2$species == "Amphistegina_lobifera"] <- "Amphistegina lobifera"
df2$morpho_species[df2$species == "Amphistegina_papillosa_Rik_9440_consensus0_1051_SaudiArabia"] <- "Amphistegina papillosa"
df2$morpho_species[df2$species == "Amphistegina_papillosa_Spermonde_3741"] <- "Amphistegina papillosa"
df2$morpho_species[df2$species == "Amphistegina_radiata_Spermonde"] <- "Amphistegina radiata"
df2$morpho_species[df2$species == "Baculogypsinoides_spinosus_8078_8079_8080"] <- "Baculogypsinoides spinosus"
df2$morpho_species[df2$species == "Borelis_schlumbergeri_3737_Spermonde"] <- "Borelis schlumbergeri"
df2$morpho_species[df2$species == "Calcarina_hispida_&_Calcarina_sp_5163_Spermonde"] <- "Calcarina hispida"
df2$morpho_species[df2$species == "Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._5247_Spermonde"] <- "Calcarina spengleri"
df2$morpho_species[df2$species == "Calcarina_sp._5164_Spermonde"] <- "Calcarina sp."
df2$morpho_species[df2$species == "Calcarina_sp._5165_Spermonde"] <- "Calcarina sp."                           
df2$morpho_species[df2$species == "Calcarina_sp._Spermonde"] <- "Calcarina sp."                  
df2$morpho_species[df2$species == "Calcarina_spengleri_Rik_9352_consensus0_1732_Spermonde"] <- "Calcarina spengleri"   
df2$morpho_species[df2$species == "Cyrea_szymborska_17247"] <- "Z_other_foraminifera"   
df2$morpho_species[df2$species == "Eggerelloides_scaber_12301"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Elphidium_Rik_9327_9337_9350_Spermonde"] <- "Elphidium sp."
df2$morpho_species[df2$species == "Encrusting_Rik_9441_consensus0_977"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Glabratellina_sp._17921_&_Planoglabratella_opercularis_18053_&_Planoglabratella_sp"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Globocassidulina_biora_17225"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Haynesina_germanica_18209"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Heterostegina_depressa_sp2_Spermonde"] <- "Heterostegina depressa"
df2$morpho_species[df2$species == "Heterostegina_depressa_sp1_Spermonde"] <- "Heterostegina depressa"
df2$morpho_species[df2$species == "Miliolid_Rik_9400_9405_Israel"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Milliolid_red_Rik_9370_9371_9372_Israel"] <- "Z_other_foraminifera"
df2$morpho_species[df2$species == "Murrayinella_globosa_18043"] <- "Z_other_foraminifera"                                                         
df2$morpho_species[df2$species == "Neorotalia_calcar_Spermonde"] <- "Neorotalia calcar"
df2$morpho_species[df2$species == "Neorotalia_gaimardi_&_Baculogypsina_sphaerulata_Spermonde"] <- "Neorotalia gaimardi"
df2$morpho_species[df2$species == "Nummulites_venosus"] <- "Nummulites venosus"
df2$morpho_species[df2$species == "Operculina_ammonoides"] <- "Operculina ammonoides"
df2$morpho_species[df2$species == "Parasorites_sp"] <- "Parasorites sp."
df2$morpho_species[df2$species == "Peneroplis_sp1_Spermonde"] <- "Peneroplis sp."
df2$morpho_species[df2$species == "Peneroplis_sp2_&_Peneroplis_pertusus_5117_&_Dendritina_ambigua_Spermonde"] <- "Peneroplis sp."
df2$morpho_species[df2$species == "Planorbulinella_sp._17893"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Psammophaga_sp._19260_&_Psammophaga_sp._19296_&_Psammophaga_sp2"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Pyrgo_Rik_9347_consensus0_1557"] <- "Z_other_foraminifera"  
df2$morpho_species[df2$species == "Sorites_sp1_Spermonde"] <- "Sorites sp."
df2$morpho_species[df2$species == "Sorites_sp2_Spermonde"] <- "Sorites sp."
df2$morpho_species[df2$species == "Stainforthia_sp._17321"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Textularid_Rik_9397_consensus4_105"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Textularid_Rik_9398_consensus128_73"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Trifarina_earlandi_&_Uvigerina_bifurcata"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Vertebralina_Rik_9325_consensus0_1751"] <- "Z_other_foraminifera" 
df2$morpho_species[df2$species == "Z_not_foraminifera"] <- "Z_not_foraminifera" 
df2$morpho_species[df2$species == "Z_other_foraminifera"] <- "Z_other_foraminifera" 

df3 <- df2 %>% 
  group_by(field.nmbr., dist_coast_km, island, depth_water, substrate, group) %>% 
  summarise(sum_count = sum(counts))

depth <- c("RF1", "RF2" ,   "RF3","2.5","2.7","3.3","3.5","3.9","4.7","4.9","5.2","5.5", 
           "5.6", "6.2","6.5","6.6","6.8","7.1","7.5","7.7", "7.8","8.1","8.3","8.4","8.7", 
           "8.9","9.7","9.8", "10","10.1","10.7","11.3","11.4","11.8","11.9","12" ,
           "12.8","12.9","13", "13.6","13.9","14.6","14.7","14.8","15.3","15.9", "16.5","16.9","17",
           "17.3","17.4","17.5","17.6", "19.1","19.5","20.1","20.2","21.3","21.7", 
           "22.2","22.6","23.2","24.5","25.2","25.6","25.9","27.4","28","30.5")

islandord <- c("Lae Lae","Barang Baringan","Samalona","Polewali","Karanrang",
               "Langkadea","K. Keke","Padjenekang","Badi","Lumulumu","Langkai",
               "Kapoposang Masdar point")

#find list of top 10 ESVs per sample
df4 <- df2[df2$group == "LBF",] %>% arrange(desc(N_mean)) %>% group_by(field.nmbr.) %>% slice(1:10)
all_esvs_top5 <- unique(df4$ESVs)

df5 <- df2[df2$ESVs %in% all_esvs_top5, ]


p1 <- ggplot(df3[df3$island != "negative",], aes(x = factor(island, islandord), y = sum_count, fill = group)) +
  geom_boxplot() +
  facet_grid(~group, scales = "free")+ 
  scale_fill_manual(values = c("darkred", "orange")) +
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank())

p2 <- ggplot(df3[df3$island != "negative",], aes(x = group, y = sum_count, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("darkred", "orange")) +
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p7 <- ggplot(df2[df2$island != "negative" & df2$group == "LBF",], aes(x = factor(depth_water, depth), y = N_mean, fill = morpho_species)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~factor(island, islandord), scales = "free")+ 
  theme_bw() + 
  scale_fill_manual(values = c(brewer.pal(10, "RdYlBu"), brewer.pal(11, "PRGn")))+
  theme(axis.text.x = element_text(angle = 90), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom")

species_all <- sort(unique(df2$morpho_species[df2$island != "negative" & df2$group == "LBF"]))

#with negatives
p10 <- ggplot(df2[df2$island == "negative" & df2$group == "LBF",], aes(x = field.nmbr., y = N_mean, fill = morpho_species)) +
  geom_bar(stat = "identity", position = "fill") +
  #facet_grid(~Substrate.type, scales = "free")+ 
  theme_bw() + 
  scale_fill_manual(limits = species_all, values = c(brewer.pal(10, "RdYlBu"), brewer.pal(11, "PRGn")))+
  theme(axis.text.x = element_text(angle = 90), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom")

unique(df2$species[df2$group == "LBF"])
unique(df2$morpho_species[df2$group == "LBF"])
unique(df2$family[df2$group == "LBF"])


#_______________DIATOMS______________________
df <- read.csv("~downloads/Datasets/merged_diat_data.csv")

df$depth_water <- df$depth..below.the.local.surface..in.meters.
df <- df[df$Phylum == "Bacillariophyta",]

#rarefy to the sample (none negative) with the lowest read number: 1924)
dfx <- df %>% 
  group_by(field.nmbr.) %>% 
  mutate(N=reads/sum(reads)*1654)

dfx <- dfx %>% 
  group_by(field.nmbr.) %>% 
  mutate(percent_N=reads/sum(reads)*100)

df2 <- dfx %>% 
  group_by(field.nmbr., ESVs,dist_coast_km, island, depth_water, substrate, Substrate.type,Identity.percentage, Class, Order, Family, Genus, Species) %>% 
  summarise(N_mean=mean(N), percent_N_mean = mean(percent_N))

df2$counts <- 1

df2$Species[df2$Identity.percentage < 95] <- "unknown species"
df2$Genus[df2$Identity.percentage < 90] <- "unknown genus"
df2$Family[df2$Identity.percentage < 85] <- "unknown family"


df3 <- df2 %>% 
  group_by(field.nmbr., dist_coast_km, island, depth_water, substrate, Substrate.type) %>% 
  summarise(sum_count = sum(counts))

islandord <- c("Lae Lae","Barang Baringan","Samalona","Polewali","Karanrang",
               "Langkadea","K. Keke","Padjenekang","Badi","Lumulumu","Langkai",
               "Kapoposang Masdar point")

#find list of top 10 ESVs per sample
df4 <- df2 %>% arrange(desc(percent_N_mean)) %>% group_by(island,depth_water, Substrate.type) %>% slice(1:10)


p3 <- ggplot(df3[df3$island != "negative" & df3$Substrate.type %in% c("seawater", "substrate"),], 
             aes(x = factor(island, islandord),y = sum_count, fill = Substrate.type)) +
  geom_boxplot() +
  facet_grid(~Substrate.type, scales = "free") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = "none")+
  ylab("Number of ESVs")

p4 <- ggplot(df3[df3$island != "negative" & df3$Substrate.type %in% c("seawater", "substrate"),], 
             aes(x = Substrate.type,y = sum_count, fill = Substrate.type)) +
  geom_boxplot() +
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank()) + 
  ylab("Number of ESVs")

p8 <- ggplot(df2[df2$substrate != "negative" & df2$Identity.percentage > 90,], aes(x = factor(depth_water, depth), y = N_mean, fill = Order)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(Substrate.type~factor(island, islandord), scales = "free")+ 
  theme_bw() + 
  scale_fill_manual(values = c(brewer.pal(10, "RdYlBu"), brewer.pal(10, "PRGn"),brewer.pal(10, "BrBG"),brewer.pal(4, "PiYG")))+
  theme(axis.text.x = element_text(angle = 90), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom")

order_all <- sort(unique(df2$Order[df2$substrate != "negative" & df2$Identity.percentage > 90]))

#with negatives
p11 <- ggplot(df2[df2$substrate == "negative" & df2$Identity.percentage > 90,], aes(x = field.nmbr., y = N_mean, fill = Order)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~Substrate.type, scales = "free")+ 
  theme_bw() + 
  scale_fill_manual(limits = order_all, values = c(brewer.pal(10, "RdYlBu"), brewer.pal(10, "PRGn"),brewer.pal(10, "BrBG"),brewer.pal(4, "PiYG")))+
  theme(axis.text.x = element_text(angle = 90), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom")

unique(df2$Class[df2$Identity.percentage > 90])
unique(df2$Order[df2$Identity.percentage > 90])
unique(df2$Family[df2$Identity.percentage > 90])
unique(df2$Genus[df2$Identity.percentage > 90])
unique(df2$Species[df2$Identity.percentage > 95])

#_______________PROKARYOTES______________________
df <- read.csv("~downloads/Datasets/merged_prok_data.csv")

df$depth_water <- df$depth..below.the.local.surface..in.meters.

#rarefy to the sample (none negative) with the lowest read number: 1924)
dfx <- df %>% 
  group_by(field.nmbr.) %>% 
  mutate(N=reads/sum(reads)*1654)

dfx <- dfx %>% 
  group_by(field.nmbr.) %>% 
  mutate(percent_N=reads/sum(reads)*100)

df2 <- dfx %>% 
  group_by(field.nmbr., ESVs,dist_coast_km, island, depth_water, substrate, Substrate.type,Identity.percentage, Phylum, Class, Order, Family, Genus, Species) %>% 
  summarise(N_mean=mean(N), percent_N_mean = mean(percent_N))

df2$counts <- 1

df3 <- df2 %>% 
  group_by(field.nmbr., dist_coast_km, island, depth_water, substrate, Substrate.type) %>% 
  summarise(sum_count = sum(counts))

islandord <- c("Lae Lae","Barang Baringan","Samalona","Polewali","Karanrang",
               "Langkadea","K. Keke","Padjenekang","Badi","Lumulumu","Langkai",
               "Kapoposang Masdar point")

p5 <- ggplot(df3[df3$island != "negative" & df3$Substrate.type %in% c("seawater", "substrate"),], 
             aes(x = factor(island, islandord),y = sum_count, fill = Substrate.type)) +
  geom_boxplot() +
  facet_grid(~Substrate.type, scales = "free") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90), legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p6 <- ggplot(df3[df3$island != "negative" & df3$Substrate.type %in% c("seawater", "substrate"),], 
             aes(x = Substrate.type,y = sum_count, fill = Substrate.type)) +
  geom_boxplot() +
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


df5 <- df2[df2$Identity.percentage > 86.5 & !is.na(df2$percent_N_mean),]
df5$legend <- NA
df5$legend[df5$percent_N_mean > 1] <- df5$Family[df5$percent_N_mean > 1]
df5$legend[df5$percent_N_mean <= 1] <- "others < 1%"
df5$legend[df5$legend == "unknown family"] <- "others > 1%"
df5$legend[df5$legend == "Unknown Family"] <- "others > 1%"
df5$legend[df5$legend == "uncultured"] <- "others > 1%"

df6 <- df5 %>% group_by(Substrate.type, legend) %>% summarize(sum_N = sum(N_mean))

df7 <- df6 %>% group_by(Substrate.type) %>% top_n(22)

top_20_fam <- unique(df7$legend)

df8 <- df5

df8$new_name <- NA
df8$new_name[df8$legend %not% top_20_fam] <- "others > 1%"
df8$new_name[df8$legend == "others > 1%"] <- "others > 1%"
df8$new_name[df8$legend == "others < 1%"] <- "others < 1%"
df8$new_name[df8$new_name %not% c("others > 1%", "others < 1%")] <- paste(substr(df8$Phylum[df8$new_name %not% c("others > 1%", "others < 1%")], 1,3), "-", df8$legend[df8$new_name %not% c("others > 1%", "others < 1%")], sep = "")

order_name <- sort(unique(df8$new_name))
index <- c(1:15, 44,43, 16:42)
order_name <- order_name[order(index)]

p9 <- ggplot(df8[df8$island != "negative" & df8$Substrate.type %in% c("seawater", "substrate"),], aes(x = factor(depth_water, depth), y = percent_N_mean, fill = factor(new_name, order_name))) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(Substrate.type~factor(island, islandord), scales = "free")+ 
  theme_bw() + 
  scale_fill_manual(values = c(brewer.pal(11, "RdYlBu"), brewer.pal(11, "PRGn")[c(1:5,7:11)],brewer.pal(11, "BrBG")[c(1:5,7:11)],brewer.pal(10, "PiYG")[c(1:3,5:10)], "gray", "black"))+
  theme(axis.text.x = element_text(angle = 90), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")

#with negatives
p12 <- ggplot(df8[df8$substrate == "negative",], aes(x = field.nmbr., y = percent_N_mean, fill = factor(new_name, order_name))) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~Substrate.type, scales = "free")+ 
  theme_bw() + 
  scale_fill_manual(limits = order_name, values = c(brewer.pal(11, "RdYlBu"), brewer.pal(11, "PRGn"),brewer.pal(10, "BrBG"),brewer.pal(10, "PiYG"), "gray", "black"))+
  theme(axis.text.x = element_text(angle = 90), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom")

unique(df2$Phylum[df2$Identity.percentage > 86.5])
unique(df2$Class[df2$Identity.percentage > 86.5])
unique(df2$Order[df2$Identity.percentage > 86.5])
unique(df2$Family[df2$Identity.percentage > 86.5])
unique(df2$Genus[df2$Identity.percentage > 94.5])

#for ESV counts
layout <- "
AAAAAB
CCCCCD
EEEEEF
"
p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(design = layout)

#for alpha diversity
layout <- "
A
B
B
C
C
"
p7 + p8 + p9 + plot_layout(design = layout)

layout <- "
A
B
C
"
p10 + p11 + p12 + plot_layout(design = layout)




