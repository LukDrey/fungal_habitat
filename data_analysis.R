readRDS("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/dada2/taxonomy_assignment/tax_table_fungi.rds")
taxa_table <- readRDS("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/dada2/taxonomy_assignment/tax_table_fungi.rds")
write.csv(taxa_table, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/dada2/taxonomy_assignment/tax_table_fungi.csv" )


#phyloseq
library("phyloseq")

#ggplot
library("ggplot2")
theme_set(theme_bw())

##create venn diagrams
#ggvenn
#install.packages("ggvenn") 
library("ggvenn")
#ps_venn
#install_github("Russel88/MicEco")
library(MicEco)

# tidyverse (data structuring)
library(tidyverse)

# ape for phylogenetic tree
library("ape")

# rstatix - statistical analysis enhancement
#install.packages("rstatix")
library(rstatix)

# ggpubr
# functions for creating and customizing ‘ggplot2’- based publication ready plots
#install.packages("ggpubr")
library(ggpubr)

# ggdist, e.g. for violin plots
#install.packages("ggdist")
library("ggdist")

# vegan - analysis of ecological data
# e.g. rarefaction curves
library("vegan")

# paletteer for different colour palettes
library("paletteer")

# ranacapa - contains the ggrare package for rarefaction in ggplot
#devtools::install_github("gauravsk/ranacapa")
library("ranacapa")

#devtools - for installation of customized packages
library(devtools)

# plyr
#library("plyr")

# Fantaxtic
#if(!"devtools" %in% installed.packages()){
#  install.packages("devtools", dependencies = T, )
#}
#devtools::install_github("gmteunisse/Fantaxtic")
library(fantaxtic)

##SpiecEasi - for network analysis
#install.packages("stringi") #was required for installation
library(stringi)
#install_github("zdk123/SpiecEasi")
library(SpiecEasi)

#rgexf - work with gexf objects
#install.packages("rgexf")
library(rgexf)

#igraph
library(igraph)

#here - define file paths
library(here)

library(xlsx)

#install.packages("dplyr")
library(dplyr)

#remove.packages("rlang")
#install.packages("rlang")
library(rlang)

#install_github("nhanhocu/metamicrobiomeR")
library(metamicrobiomeR) 

#install.packages("remotes")
library(remotes)

#microbiome
library(microbiome)

#install.packages("gghalves")
library(gghalves)

#install.packages("viridis")
library(viridis)

#install.packages("pals")
library(pals)

#install.packages("influential")
library(influential)

#install.packages("ade4")
library(ade4)

######################################################
##########create & work with phyloseq object###########
#######################################################

####phyloseq object with ASV table and taxonomy table###
# -> physeq1
otumat <- readRDS("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/decontam_lulu/curated_ASV_table_fungi.rds")
ASV = otu_table(otumat, taxa_are_rows = TRUE)
taxmat <- readRDS("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/dada2/taxonomy_assignment/tax_table_fungi.rds")
taxmat = matrix(taxmat, nrow = nrow(taxmat), ncol = 7)
rownames(taxmat) <- paste0("ASV_", 1:nrow(taxmat))
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
class(taxmat)
TAX = tax_table(taxmat)
ASV
TAX
head(TAX)

physeq1 = phyloseq(ASV, TAX)

#bar plot by taxonomy rank (too many samples, too many ranks)
plot_bar(physeq1, fill = "Class")
plot_bar(physeq4, fill = "Phylum")
plot_bar(physeq5, fill = "Phylum")
#plot area is too small to display anything..
#plot area adjustment to DIN A0 3370 x 2384
#plot again even larger

#heatmap (too many samples, too many ranks)
plot_heatmap(physeq1, taxa.label="Class")
#plot area is too small to see anything...

#import file for sample_data 
install.packages("tidyverse")
library(tidyverse)
sdata <- read.csv2("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/sample_data_1.csv")
colnames(sdata) <- c("sample","exploratory","dominant_tree","leaf_shape","substrate","s_or_c")
sdata2 <- sdata %>% remove_rownames %>% column_to_rownames(var="sample")
head(sdata2)
sampledata = sample_data(sdata2)
head(sampledata)

#phylogenetic tree
library("ape")
random_tree = rtree(ntaxa(physeq1), rooted=TRUE, tip.label=taxa_names(physeq1))
plot(random_tree)
#nonsense

#add sample data to phyloseq object
# -> physeq3
#physeq3 = phyloseq(ASV, TAX, sampledata, random_tree)
#physeq3

### phyloseq object with all sample data fungi ###
#20220107 updated sample_data file without control and 
# information for additional plots
#20220224 included plot as sample data information:
#20220226 included ForMI
sdata3 <- read.csv2("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/sample_data_3_AHS.csv")
colnames(sdata3) <- c("sample","exploratory","dominant_tree","tree_type","substrate","library_size","plot","ForMI","ForMIclass","Iharv","Inonat","Idwcut")
sdata4 <- sdata3 %>% remove_rownames %>% column_to_rownames(var="sample")
head(sdata4)
sampledata2 = sample_data(sdata4)
head(sampledata2)
# phyloseq object with updated sample data -> physeq4
physeq4 = phyloseq(ASV, TAX, sampledata2, random_tree)
physeq4

## import sample data and generate phyloseq object ##
#20220129 updated sample_data file with library sizes and
# reduced to Alb + Schorfheide --> object physeq5
#20220226 included plot as sample data information:
sdata5 <- read.csv2("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/sample_data_3_AS.csv")
colnames(sdata5) <- c("sample","exploratory","dominant_tree","tree_type","substrate","library_size","plot")
sdata6 <- sdata5 %>% remove_rownames %>% column_to_rownames(var="sample")
head(sdata6)
sampledata3 = sample_data(sdata6)
head(sampledata3)
# phyloseq object Alb + Schorfheide + library size -> physeq5
physeq5 = phyloseq(ASV, TAX, sampledata3, random_tree)
physeq5
physeq5_s <- subset_samples(physeq5, substrate=="soil")
physeq5_s
physeq5_s1 <- prune_taxa(taxa_sums(physeq5_s) > 0, physeq5_s)
physeq5_s1
physeq5_b <- subset_samples(physeq5, substrate=="bark")
physeq5_b
physeq5_b1 <- prune_taxa(taxa_sums(physeq5_b) > 0, physeq5_b)
physeq5_b1

## create sub phyloseq objects for Alb and Schorfheide ##
#20220226 included plot as sample data information:
#20220226 included ForMI
#20220308 added venn class to sample data

#phyloseq object Alb#
sdata7 <- read.csv2("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/sample_data_3_A.csv")
colnames(sdata7) <- c("sample","exploratory","dominant_tree","tree_type","substrate","library_size", "plot","ForMI","ForMIclass","Iharv","Inonat","Idwcut","venn_class")
sdata8 <- sdata7 %>% remove_rownames %>% column_to_rownames(var="sample")
head(sdata8)
sampledataAlb = sample_data(sdata8)
head(sampledataAlb)
physeqAlb = phyloseq(ASV, TAX, sampledataAlb, random_tree)
physeqAlb
physeqAlb_1 <- prune_taxa(taxa_sums(physeqAlb) > 0, physeqAlb)
physeqAlb_1

#phyloseq object Schorfheide#
sdata9 <- read.csv2("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/sample_data_3_S.csv")
colnames(sdata9) <- c("sample","exploratory","dominant_tree","tree_type","substrate","library_size", "plot","ForMI","ForMIclass","Iharv","Inonat","Idwcut","venn_class")
sdata10 <- sdata9 %>% remove_rownames %>% column_to_rownames(var="sample")
head(sdata10)
sampledataSchorf = sample_data(sdata10)
head(sampledataSchorf)
physeqSchorf = phyloseq(ASV, TAX, sampledataSchorf, random_tree)
physeqSchorf
physeqSchorf_1 <- prune_taxa(taxa_sums(physeqSchorf) > 0, physeqSchorf)
physeqSchorf_1

########20220304: subsets soil / bark per exploratory Alb + Schorfheide

#phyloseq object Alb soil
physeqAlb_s <- subset_samples(physeqAlb, substrate=="soil")
physeqAlb_s
physeqAlb_s1 <- prune_taxa(taxa_sums(physeqAlb_s) > 0, physeqAlb_s)
physeqAlb_s1


#phyloseq object Alb bark
physeqAlb_b <- subset_samples(physeqAlb, substrate=="bark")
physeqAlb_b
physeqAlb_b1 <- prune_taxa(taxa_sums(physeqAlb_b) > 0, physeqAlb_b)
physeqAlb_b1
table(tax_table(physeqAlb_b1)[, 2])
GP.chl = subset_taxa(physeqAlb_b1, is.na(Phylum))
ps_Alb_b_na_rel_abund = phyloseq::transform_sample_counts(GP.chl, function(x){x / sum(x)})
taxa_sums(ps_Alb_b_na_rel_abund)
most_abundant_taxa = sort(taxa_sums(ps_Alb_b_na_rel_abund), TRUE)[1:173]


#phyloseq object Schorfheide soil
physeqSchorf_s <- subset_samples(physeqSchorf, substrate=="soil")
physeqSchorf_s
physeqSchorf_s1 <- prune_taxa(taxa_sums(physeqSchorf_s) > 0, physeqSchorf_s)
physeqSchorf_s1

#phyloseq object Schorfheide bark
physeqSchorf_b <- subset_samples(physeqSchorf, substrate=="bark")
physeqSchorf_b
physeqSchorf_b1 <- prune_taxa(taxa_sums(physeqSchorf_b) > 0, physeqSchorf_b)
physeqSchorf_b1
table(tax_table(physeqSchorf_b1)[, 2])
GP.chl2 = subset_taxa(physeqSchorf_b1, is.na(Phylum))
ps_Schorf_b_na_rel_abund = phyloseq::transform_sample_counts(GP.chl2, function(x){x / sum(x)})
taxa_sums(ps_Schorf_b_na_rel_abund)
most_abundant_taxa_S = sort(taxa_sums(ps_Schorf_b_na_rel_abund), TRUE)[1:137]

########20220311: subsets soil / bark per exploratory Alb + Schorfheide & per tree
##for more detailed networks to define hub taxa

#phyloseq object Alb soil Fagus
physeqAlb_sF <- subset_samples(physeqAlb, venn_class=="F.sylvatica_soil")
physeqAlb_sF
OTU_AsF = as(otu_table(physeqAlb_sF), "matrix")
OTU_AsF
tax_table(physeqAlb_sF)
otu_table(physeqAlb_sF)
physeqAlb_sF1 <- prune_taxa(taxa_sums(physeqAlb_sF) > 0, physeqAlb_sF)
physeqAlb_sF1

#phyloseq object Alb bark Fagus
physeqAlb_bF <- subset_samples(physeqAlb, venn_class=="F.sylvatica_bark")
physeqAlb_bF
physeqAlb_bF1 <- prune_taxa(taxa_sums(physeqAlb_bF) > 0, physeqAlb_bF)
physeqAlb_bF1

#phyloseq object Alb soil Picea
physeqAlb_sP <- subset_samples(physeqAlb, venn_class=="P.abies_soil")
physeqAlb_sP
physeqAlb_sP1 <- prune_taxa(taxa_sums(physeqAlb_sP) > 0, physeqAlb_sP)
physeqAlb_sP1

#phyloseq object Alb bark Picea
physeqAlb_bP <- subset_samples(physeqAlb, venn_class=="P.abies_bark")
physeqAlb_bP
physeqAlb_bP1 <- prune_taxa(taxa_sums(physeqAlb_bP) > 0, physeqAlb_bP)
physeqAlb_bP1

#phyloseq object Schorfheide soil Fagus
physeqSchorf_sF <- subset_samples(physeqSchorf, venn_class=="F.sylvatica_soil")
physeqSchorf_sF
physeqSchorf_sF1 <- prune_taxa(taxa_sums(physeqSchorf_sF) > 0, physeqSchorf_sF)
physeqSchorf_sF1

#phyloseq object Schorfheide bark Fagus
physeqSchorf_bF <- subset_samples(physeqSchorf, venn_class=="F.sylvatica_bark")
physeqSchorf_bF
physeqSchorf_bF1 <- prune_taxa(taxa_sums(physeqSchorf_bF) > 0, physeqSchorf_bF)
physeqSchorf_bF1

#phyloseq object Schorfheide soil Pinus
physeqSchorf_sP <- subset_samples(physeqSchorf, venn_class=="P.sylvestris_soil")
physeqSchorf_sP
physeqSchorf_sP1 <- prune_taxa(taxa_sums(physeqSchorf_sP) > 0, physeqSchorf_sP)
physeqSchorf_sP1

#phyloseq object Schorfheide bark Pinus
physeqSchorf_bP <- subset_samples(physeqSchorf, venn_class=="P.sylvestris_bark")
physeqSchorf_bP
physeqSchorf_bP1 <- prune_taxa(taxa_sums(physeqSchorf_bP) > 0, physeqSchorf_bP)
physeqSchorf_bP1

########20220313: subsets Fagus / Picea / Pinus 
##for stacked bar plots

#phyloseq object Alb Fagus
physeqAlb_F <- subset_samples(physeqAlb, dominant_tree=="Fagus_sylvatica")
physeqAlb_F

#phyloseq object Alb Picea
physeqAlb_P <- subset_samples(physeqAlb, dominant_tree=="Picea_abies")
physeqAlb_P

#phyloseq object Schorfheide Fagus
physeqSchorf_F <- subset_samples(physeqSchorf, dominant_tree=="Fagus_sylvatica")
physeqSchorf_F

#phyloseq object Schorfheide Pinus
physeqSchorf_P <- subset_samples(physeqSchorf, dominant_tree=="Pinus_sylvestris")
physeqSchorf_P


###########analyze library sizes coniferous - deciduous (unbalanced sample design)################
#boxplot library size Alb + Schorfheide
libsize <- read.csv2("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/sample_data_3_AS.csv")
boxplot(libsize$library_size[libsize$tree_type=="coniferous"],
        libsize$library_size[libsize$tree_type=="deciduous"],
        main="boxplot library size coniferous vs. deciduous (A + S)", names = c("coniferous","deciduous"), 
        col = c("darkgreen","green2")) 
text(1.5,130000,"n.s.",col = "firebrick4", cex = 2)
##test statistics
#Welch t-test
libconif <- libsize$library_size[libsize$tree_type=="coniferous"]
libdecid <- libsize$library_size[libsize$tree_type=="deciduous"]
t.test(libconif, libdecid)
#Wilcoxon test
wilcox.test(library_size~tree_type, data = libsize, exact = FALSE, correct = FALSE, conf.int = FALSE)

##Alb
#boxplot library size Alb
libsizeA <- read.csv2("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/sample_data_3_A.csv")
boxplot(libsizeA$library_size[libsizeA$tree_type=="coniferous"],
        libsizeA$library_size[libsizeA$tree_type=="deciduous"],
        main="boxplot library size coniferous vs. deciduous (Alb)", names = c("coniferous","deciduous"), 
        col = c("darkgreen","green2"))
text(1.5,130000,"n.s.",col = "firebrick4", cex = 2)
#test statistics Alb
#Welch t-test
libconifA <- libsizeA$library_size[libsizeA$tree_type=="coniferous"]
libdecidA <- libsizeA$library_size[libsizeA$tree_type=="deciduous"]
t.test(libconifA, libdecidA)
#Wilcoxon test
wilcox.test(library_size~tree_type, data = libsizeA, exact = FALSE, correct = FALSE, conf.int = FALSE)

##Schorfheide
#boxplot library size Schorfheide
libsizeS <- read.csv2("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/sample_data_3_S.csv")
boxplot(libsizeS$library_size[libsizeS$tree_type=="coniferous"],
        libsizeS$library_size[libsizeS$tree_type=="deciduous"],
        main="boxplot library size coniferous vs. deciduous (Schorfheide)", names = c("coniferous","deciduous"), 
        col = c("darkgreen","green2"))
text(1.5,130000,"n.s.",col = "firebrick4", cex = 2)
#test statistics Schorfheide
#Welch t-test
libconifS <- libsizeS$library_size[libsizeS$tree_type=="coniferous"]
libdecidS <- libsizeS$library_size[libsizeS$tree_type=="deciduous"]
t.test(libconifS, libdecidS)
#Wilcoxon test
wilcox.test(library_size~tree_type, data = libsizeS, exact = FALSE, correct = FALSE, conf.int = FALSE)

#boxplot Schorfheide + Alb separated
boxplot(libsizeA$library_size[libsizeA$tree_type=="coniferous"],
        libsizeA$library_size[libsizeA$tree_type=="deciduous"],
  libsizeS$library_size[libsizeS$tree_type=="coniferous"],
        libsizeS$library_size[libsizeS$tree_type=="deciduous"],
        main="boxplot library size coniferous vs. deciduous", 
  names = c("coniferous Alb","deciduous Alb", "coniferous Schorfh.", "deciduous Schorfh."), 
        col = c("darkgreen","green2","darkgreen","green2"))

####two way ANOVA for the boxplot library sizes 
#not applicable for the dataframe!
install.packages("rstatix")
install.packages("ggpubr")
library(rstatix)
library(ggpubr)
set.seed(123)
libsize %>% sample_n_by(exploratory, tree_type, size = 1)
libsize %>%
  group_by(exploratory, tree_type) %>%
  get_summary_stats(library_size, type = "mean_sd")
#identify outliers
libsize %>%
  group_by(exploratory, tree_type) %>%
  identify_outliers(library_size)
#check normality assumption
# Build the linear model
model  <- lm(library_size ~ exploratory*tree_type,
             data = libsize)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
#-> ANoVA nur zulässig bei normalverteilten Residuen, hier nicht der Fall
# Alternative: Kruskal-Wallis test 
describeBy(libsize$library_size,libsize$exploratory)
describeBy(libsize$library_size,libsize$tree_type)
kruskal.test(libsize$library_size~libsize$exploratory)
kruskal.test(libsize$library_size~libsize$tree_type)

###rarefaction curves Alb + Schorfheide separated###
library("vegan")
library("ranacapa")
library(devtools)
devtools::install_github("gauravsk/ranacapa")
###rarefaction curve samples Alb###
rarecurve <- vegan::rarecurve(t(otu_table(physeqAlb)), step = 50, cex = 0.5,  xlim=c(0, 90000), label = F, col = "dodgerblue4")
rarecurveAlb <- ggrare(otu_table(physeqAlb), step = 50, color = "substrate", label = "rarefaction Alb by substrate") 
rlang::last_error()
rlang::last_trace()
plot(rarecurveAlb)

rareAlb <- ggrare(physeqAlb, step = 50, color = "tree_type", se = FALSE)
ggsave(
  "rarefaction_Alb.png",
  width = 8.3,
  height = 5.7,
  dpi = 800
)
rareAlb <- ggrare(physeqAlb, step = 50, color = "substrate", se = FALSE)

#20220728 - color by venn_class
rareAlb <- ggrare(physeqAlb, step = 50, color = "venn_class", se = FALSE)
ggsave(
  "rarefaction_Alb_venn.png",
  width = 8.3,
  height = 5.7,
  dpi = 800
)


###rarefaction curve samples Schorfheide###
rarecurveSchorf <- vegan::rarecurve(t(otu_table(physeqSchorf)), step = 50, cex = 0.5,  xlim=c(0, 90000), label = F, col = "dodgerblue4")

rareSchorf <- ggrare(physeqSchorf, step = 50, color = "tree_type", se = FALSE)
ggsave(
  "rarefaction_Schorf.png",
  width = 8.3,
  height = 5.7,
  dpi = 800
)
rareSchorf <- ggrare(physeqSchorf, step = 50, color = "substrate", se = FALSE)

#20220728 - color by venn_class
rareSchorf <- ggrare(physeqSchorf, step = 50, color = "venn_class", se = FALSE)
ggsave(
  "rarefaction_Schorf_venn.png",
  width = 8.3,
  height = 5.7,
  dpi = 800
)

#######examples phyloseq3################

#trees make no sense with 11952 ASVs
#plot_tree(physeq3, color="exploratory", label.tips="substrate", ladderize="left", plot.margin=0.3)
#plot_tree(physeq3, color="substrate", label.tips="substrate", ladderize="left", plot.margin=0.3)

########################################
########alpha diversity#################
########################################


#older versions / physeq object
#plot_richness(physeq3)
#plot_richness(physeq3, measures = c("Observed"))
#plot_richness(physeq3, x="leaf_shape", measures = c("Observed","Shannon"))
#plot_richness(physeq3, x="dominant_tree", measures = c("Observed","Shannon"))
#-> control samples sollten aus den Daten rausgeschmissen werden
#plot_richness(physeq3, x="exploratory", measures = c("Observed","Shannon"))
#plot_richness(physeq3, x="exploratory", color = "substrate", measures = c("Observed","Shannon"))
#plot_richness(physeq3, x="exploratory", color = "substrate", measures = c("Chao1"))
#plot_richness(physeq3, x="leaf_shape", color = "substrate", measures = c("Observed","Shannon"))
#plot_richness(physeq3, x="substrate", color = "leaf_shape", measures = c("Observed","Shannon"))

install.packages("ggdist")
library("ggdist")

###by substrate
plot_richness(physeq4, x="substrate", measures = c("Chao1","Shannon"), color = "substrate") + scale_colour_manual(values = c("burlywood4","chocolate1")) + 
  geom_boxplot(outlier.shape  = NA) + geom_jitter(aes(), height = 0, width = .2) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),axis.ticks.x = element_blank())
#Shannon diversity all exploratories by subsrate
plot_richness(physeq4, x="substrate", measures = c("Shannon"), color = "substrate") + 
  scale_colour_manual(values = c("burlywood4","chocolate1")) + 
  geom_boxplot(outlier.shape  = NA, fill = c("burlywood", "chocolate3")) + geom_jitter(aes(), height = 0, width = .2) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),axis.ticks.x = element_blank()) +
  geom_text(x=1.5, y=5.1, label="***", color="firebrick4", size=13) +
  ## add half-violin from {ggdist} package
  ggdist::stat_halfeye(
    ## custom bandwidth
    adjust = .5, 
    ## adjust height
    width = .7, 
    ## move geom to the right
    justification = -.01,
    ## remove slab interval(transparency)
    .width = 0, 
    point_colour = NA,
    alpha = 0.5
  ) 
#annotate("text", x=1, y=5, label= "***", color = "deeppink", size = 15)
#by tree_type
plot_richness(physeq4, x="tree_type", measures = c("Chao1","Shannon"), color = "tree_type") + scale_colour_manual(values = c("darkgreen","green2")) +
  geom_boxplot(outlier.shape  = NA) + geom_jitter(aes(), height = 0, width = .2) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),axis.ticks.x = element_blank())
#by dominant_tree
plot_richness(physeq4, x="dominant_tree", measures = c("Shannon"), color = "dominant_tree") +
  geom_boxplot(outlier.shape  = NA) + geom_jitter(aes(), height = 0, width = .2)
#by dominant_tree and substrate
plot_richness(physeq4, x="dominant_tree", measures = c("Shannon"), color = "substrate") +
  scale_colour_manual(values = c("burlywood4","chocolate1")) +
  geom_boxplot(outlier.shape  = NA) + geom_jitter(aes(), height = 0, width = .2)
#by exploratory
plot_richness(physeq4, x="exploratory", color = "substrate", measures = c("Shannon")) +
  scale_colour_manual(values = c("burlywood4","chocolate1")) +
  geom_boxplot(outlier.shape  = NA) + geom_jitter(aes(), height = 0, width = .2)
plot_richness(physeq4, x="exploratory", color = "exploratory", measures = c("Shannon")) +
  scale_colour_manual(values = c("hotpink","dodgerblue","seagreen4")) +
  geom_boxplot(outlier.shape  = NA) + geom_jitter(aes(), height = 0, width = .2) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),axis.ticks.x = element_blank()) +
  geom_text(x=3, y=5.1, label="***", color="firebrick4", size=13)
#by substrate and tree_type
plot_richness(physeq4, x="substrate", color = "tree_type", measures = c("Shannon"))+ scale_colour_manual(values = c("darkgreen","green2")) +
  geom_boxplot(outlier.shape  = NA) + geom_jitter(aes(), height = 0, width = .2) 
#by tree type and substrate
plot_richness(physeq4, x="tree_type", color = "substrate", measures = c("Shannon"))+ scale_colour_manual(values = c("burlywood4","chocolate1")) +
  geom_boxplot(outlier.shape  = NA) + geom_jitter(aes(), height = 0, width = .2) +
  ## add half-violin from {ggdist} package
  ggdist::stat_halfeye(
    ## custom bandwidth
    adjust = .5, 
    ## adjust height
    width = .7, 
    ## move geom to the right
    justification = -.15, 
    ## remove slab interval
    .width = 0, 
    point_colour = NA,
    alpha = 0.5
  ) 

###by plot##
plot_richness(physeq4, x="plot", color = "substrate", measures = c("Shannon")) +
  scale_colour_manual(values = c("burlywood4","chocolate1")) +
  geom_boxplot(outlier.shape  = NA) + geom_jitter(aes(), height = 0, width = .2)


###by management intensity##
plot_richness(physeq4, x="ForMIclass", measures = c("Shannon"), color = "ForMIclass") + 
  scale_colour_manual(values = c("firebrick2","seagreen2","palegreen","seashell4")) + 
  geom_boxplot(outlier.shape  = NA, fill = c("firebrick2","seagreen2","palegreen","seashell4")) + geom_jitter(aes(), height = 0, width = .2) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),axis.ticks.x = element_blank()) +
  ## add half-violin from {ggdist} package
  ggdist::stat_halfeye(
    ## custom bandwidth
    adjust = .5, 
    ## adjust height
    width = .7, 
    ## move geom to the right
    justification = -.01,
    ## remove slab interval(transparency)
    .width = 0, 
    point_colour = NA,
    alpha = 0.5
  ) 

################Alb + Schorfheide###########################

################Swabian Alb###########################

# Calculate the Shannon Diversity.
shannon_Alb <- estimate_richness(physeqAlb, split = T, measures = 'Shannon')
# Make a new sample_data object and merge with existing phyloseq object.
new_sampledata_Alb <- sample_data(shannon_Alb)
physeqAlb_div <- merge_phyloseq(physeqAlb, new_sampledata_Alb)
fungi_sampledata_Alb <- sample_data(physeqAlb_div)

#basic boxplot for statistical data extraction
shannon_plot_Alb0 <- ggplot(fungi_sampledata_Alb, aes(y= Shannon)) + geom_boxplot(
  ## remove outliers
  outlier.color = NA## `outlier.shape = NA` works as well 
)
shannon_plot_Alb0
divdataA <- ggplot_build(shannon_plot_Alb0)$data
divdataA

####by tree type (in this case species) and substrate for Swabian Alb

shannon_plot_Alb1 <- ggplot(fungi_sampledata_Alb, aes(x = substrate, y = Shannon, fill = dominant_tree, colour = dominant_tree)) +
  scale_fill_manual(values = c("green","darkgreen"),name = "dominant tree species", labels = c("Fagus sylvatica", "Picea abies")) +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(face = "italic", size = 12))+
  theme(legend.position = c(0.15, 0.9))+
  scale_color_manual(values = c("green3","mediumseagreen"), guide = "none") +
  geom_boxplot(
    ## remove outliers
    outlier.color = NA## `outlier.shape = NA` works as well 
  ) +
  theme(axis.text = element_text(size =12))  +
  theme(axis.title = element_text(size = 12)) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point() +
  ylab('alpha diversity (Shannon)') + 
  ylim (2, 5.5) 
shannon_plot_Alb1
ggsave(
  "alpha_diversity_Alb1.png",
  width = 8.3,
  height = 5.8,
  dpi = 800
)

#statistics
ggplot_build(shannon_plot_Alb_1)$data


plot_richness(physeqAlb, x="substrate", color = "dominant_tree", measures = c("Shannon"),
              title = "fungal Shannon diversity by substrate & dominant tree Swabian Alb exploratory") +
  geom_boxplot(outlier.shape  = NA,fill = c("green3", "mediumseagreen","green3", "mediumseagreen")) + 
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .4, 
    ## add some transparency
    alpha = .5,
    size = .5) +
  scale_colour_manual(values = c("green","darkgreen")) +
  theme(legend.position = c(0.15, 0.9)) +
  ylab("alpha diversity (Shannon)") +
  ylim(1.5, 5.5) 
  ## add half-violin from {ggdist} package
  ggdist::stat_halfeye(
    ## custom bandwidth
    adjust = .5, 
    ## adjust height
    width = .7, 
    ## move geom to the right
    justification = -.15, 
    ## remove slab interval
    .width = 0, 
    point_colour = NA,
    alpha = 0.5
  ) 

  ###just by tree type (species) ##
  
  shannon_plot_Alb2 <- ggplot(fungi_sampledata_Alb, aes(x = dominant_tree, y = Shannon, fill = dominant_tree, colour = dominant_tree)) +
    scale_fill_manual(values = c("green","darkgreen")) +
    scale_color_manual(values = c("green3","mediumseagreen")) +
    geom_boxplot(
      ## remove outliers
      outlier.color = NA, show.legend = FALSE## `outlier.shape = NA` works as well 
    ) +
    ## add justified jitter from the {gghalves} package
    gghalves::geom_half_point() +
    ylab('alpha diversity (Shannon)') + 
    xlab("dominant tree species") +
    ylim (2, 5.5) +
    theme(legend.position = "none")+
    scale_x_discrete(labels=c("Fagus_sylvatica" = "Fagus sylvatica", "Picea_abies" = "Picea abies"))+
    theme(axis.text = element_text(face = "italic", size =12))  +
    theme(axis.title = element_text(size = 12)) +
    #theme(legend.position = c(0.15, 0.1)) +
    #geom_text(x=1.5, y=5.3, label="n.s.", color="firebrick4", size=6)+
    ggdist::stat_halfeye(
      ## custom bandwidth
      adjust = .5, 
      ## adjust height
      width = .7, 
      ## move geom to the right
      justification = -.15, 
      ## remove slab interval
      .width = 0, 
      point_colour = NA,
      alpha = 0.2
    )
  shannon_plot_Alb2
  ggsave(
    "alpha_diversity_Alb2.png",
    width = 8.3,
    height = 5.8,
    dpi = 800
  )
  
  
  plot_richness(physeqAlb, x="dominant_tree", measures = c("Shannon"),
                title = "Fungal Shannon diversity by dominant tree Swabian Alb exploratory") +
    geom_boxplot(outlier.shape  = NA,fill = c("green3", "mediumseagreen")) + 
    geom_jitter(aes(), height = 0, width = .2) +
    scale_colour_manual(values = c("green","darkgreen")) +
    geom_text(x=1.5, y=5.3, label="n.s.", color="firebrick4", size=13) +
    ## add half-violin from {ggdist} package
    ggdist::stat_halfeye(
      ## custom bandwidth
      adjust = .5, 
      ## adjust height
      width = .7, 
      ## move geom to the right
      justification = -.15, 
      ## remove slab interval
      .width = 0, 
      point_colour = NA,
      alpha = 0.5
    ) 

##diversity by plot Alb##
plot_richness(physeqAlb, x="plot", color = "substrate", measures = c("Shannon")) +
    scale_colour_manual(values = c("burlywood4","chocolate1")) +
    geom_boxplot(outlier.shape  = NA) + geom_jitter(aes(), height = 0, width = .2)  

  
##statistics for bark and soil Swabian Alb
  #unwrap data from phyloseq object
  testAlb = estimate_richness(physeqAlb, measures = "Shannon")
  A = sample_data(physeqAlb)
  #statistic diversity bark - soil
  bark = testAlb[A[,"substrate"] == "bark"]
  soil = testAlb[A[,"substrate"] == "soil"]
  wilcox.test(bark, soil)
  ## test within groups bark & soil
  #bark Fagus vs. Picea
  barkFagus = bark[A[,"dominant_tree"] == "Fagus_sylvatica"]
  barkPicea = bark[A[,"dominant_tree"] == "Picea_abies"]
  wilcox.test(barkFagus, barkPicea)
  #soil Fagus vs. Picea
  soilFagus = soil[A[,"dominant_tree"] == "Fagus_sylvatica"]
  soilPicea = soil[A[,"dominant_tree"] == "Picea_abies"]
  wilcox.test(soilFagus, soilPicea)
##statistics for Fagus vs. Picea Swabian Alb
  fagusA = testAlb[A[,"dominant_tree"] == "Fagus_sylvatica"]
  piceaA = testAlb[A[,"dominant_tree"] == "Picea_abies"]
  wilcox.test(fagusA, piceaA)


#example violin plot 
plot_richness(physeqAlb, x="substrate", color = "dominant_tree", measures = c("Shannon"),
              title = "fungal Shannon diversity by substrate & dominant tree Swabian Alb exploratory") +
  scale_colour_manual(values = c("green2","darkgreen")) +
  geom_violin() + scale_fill_manual(values=c("green2", "darkgreen")) +
  geom_jitter(aes(), height = 0, width = .2, shape=16) 

################Schorfheide###########################

# Calculate the Shannon Diversity.
shannon_Schorf <- estimate_richness(physeqSchorf, split = T, measures = 'Shannon')
# Make a new sample_data object and merge with existing phyloseq object.
new_sampledata_Schorf <- sample_data(shannon_Schorf)
physeqSchorf_div <- merge_phyloseq(physeqSchorf, new_sampledata_Schorf)
fungi_sampledata_Schorf <- sample_data(physeqSchorf_div)

#basic boxplot for statistical data extraction
shannon_plot_Schorf0 <- ggplot(fungi_sampledata_Schorf, aes(y= Shannon)) + geom_boxplot(
  ## remove outliers
  outlier.color = NA## `outlier.shape = NA` works as well 
)
shannon_plot_Schorf0
divdataS <- ggplot_build(shannon_plot_Schorf0)$data
divdataS


####by tree type (in this case species) and substrate for Schorfheide##

shannon_plot_Schorf1 <- ggplot(fungi_sampledata_Schorf, aes(x = substrate, y = Shannon, fill = dominant_tree, colour = dominant_tree)) +
  scale_fill_manual(values = c("green","darkgreen"),name = "dominant tree species", labels = c("Fagus sylvatica", "Pinus sylvestris")) +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(face = "italic", size = 12))+
  theme(legend.position = c(0.15, 0.9))+
  scale_color_manual(values = c("green3","mediumseagreen"), guide = "none") +
  geom_boxplot(
    ## remove outliers
    outlier.color = NA## `outlier.shape = NA` works as well 
  ) +
  theme(axis.text = element_text(size =12))  +
  theme(axis.title = element_text(size = 12)) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point() +
  ylab('alpha diversity (Shannon)') + 
  ylim (2, 5.5) 
shannon_plot_Schorf1
ggsave(
  "alpha_diversity_Schorf1.png",
  width = 8.3,
  height = 5.8,
  dpi = 800
)

#statistics
ggplot_build(shannon_plot_Schorf1)$data

plot_richness(physeqSchorf, x="substrate", color = "dominant_tree", measures = c("Shannon"),
              title = "fungal Shannon diversity by substrate & dominant tree Schorfheide exploratory") +
  geom_boxplot(outlier.shape  = NA,fill = c("green3", "mediumseagreen","green3", "mediumseagreen")) + 
  geom_jitter(aes(), height = 0, width = .2) +
  scale_colour_manual(values = c("green","darkgreen")) 
## add half-violin from {ggdist} package
ggdist::stat_halfeye(
  ## custom bandwidth
  adjust = .5, 
  ## adjust height
  width = .7, 
  ## move geom to the right
  justification = -.15, 
  ## remove slab interval
  .width = 0, 
  point_colour = NA,
  alpha = 0.5
) 

###just by tree type (species) --> shown significance##

shannon_plot_Schorf2 <- ggplot(fungi_sampledata_Schorf, aes(x = dominant_tree, y = Shannon, fill = dominant_tree, colour = dominant_tree)) +
  scale_fill_manual(values = c("green","darkgreen")) +
  scale_color_manual(values = c("green3","mediumseagreen")) +
  geom_boxplot(
    ## remove outliers
    outlier.color = NA, show.legend = FALSE## `outlier.shape = NA` works as well 
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point() +
  ylab('alpha diversity (Shannon)') + 
  xlab("dominant tree species") +
  ylim (2, 5.5) +
  theme(legend.position = "none")+
  scale_x_discrete(labels=c("Fagus_sylvatica" = "Fagus sylvatica", "Pinus_sylvestris" = "Pinus sylvestris"))+
  theme(axis.text = element_text(face = "italic", size =12))  +
  theme(axis.title = element_text(size = 12)) +
  #theme(legend.position = c(0.15, 0.1)) +
  #geom_text(x=1.5, y=5.3, label="n.s.", color="firebrick4", size=6)+
  ggdist::stat_halfeye(
    ## custom bandwidth
    adjust = .5, 
    ## adjust height
    width = .7, 
    ## move geom to the right
    justification = -.15, 
    ## remove slab interval
    .width = 0, 
    point_colour = NA,
    alpha = 0.2
  )
shannon_plot_Schorf2
ggsave(
  "alpha_diversity_Schorf2.png",
  width = 8.3,
  height = 5.8,
  dpi = 800
)



plot_richness(physeqSchorf, x="dominant_tree", measures = c("Shannon"),
              title = "Fungal Shannon diversity by dominant tree Schorfheide exploratory") +
  geom_boxplot(outlier.shape  = NA,fill = c("green3", "mediumseagreen")) + 
  geom_jitter(aes(), height = 0, width = .2) +
  scale_colour_manual(values = c("green","darkgreen")) +
  geom_text(x=1.5, y=4.5, label="*", color="firebrick4", size=13) +
  ## add half-violin from {ggdist} package
  ggdist::stat_halfeye(
    ## custom bandwidth
    adjust = .5, 
    ## adjust height
    width = .7, 
    ## move geom to the right
    justification = -.15, 
    ## remove slab interval
    .width = 0, 
    point_colour = NA,
    alpha = 0.5
  ) 

###by plot (image saved - diversity_plots_substrate_Schorf)##
#dimension 700 x 671
plot_richness(physeqSchorf, x="plot", color = "substrate", measures = c("Shannon")) +
  scale_colour_manual(values = c("burlywood4","chocolate1")) +
  geom_boxplot(outlier.shape  = NA) + geom_jitter(aes(), height = 0, width = .2)

###statistics for bark and soil Schorfheide
#unwrap data from phyloseq object
testSchorf = estimate_richness(physeqSchorf, measures = "Shannon")
S = sample_data(physeqSchorf)
#statistic diversity bark - soil
bark = testSchorf[S[,"substrate"] == "bark"]
soil = testSchorf[S[,"substrate"] == "soil"]
wilcox.test(bark, soil)
## test within groups bark & soil
#bark Fagus vs. Pinus
barkFagus = bark[S[,"dominant_tree"] == "Fagus_sylvatica"]
barkPinus = bark[S[,"dominant_tree"] == "Pinus_sylvestris"]
wilcox.test(barkFagus, barkPinus)
#soil Fagus vs. Pinus
soilFagus = soil[S[,"dominant_tree"] == "Fagus_sylvatica"]
soilPinus = soil[S[,"dominant_tree"] == "Pinus_sylvestris"]
wilcox.test(soilFagus, soilPinus)
##statistics for Fagus vs. Pinus Schorfheide
fagusS = testSchorf[S[,"dominant_tree"] == "Fagus_sylvatica"]
pinusS = testSchorf[S[,"dominant_tree"] == "Pinus_sylvestris"]
wilcox.test(fagusS, pinusS)

###example violin plot##
plot_richness(physeqSchorf, x="substrate", color = "dominant_tree", measures = c("Shannon"),
              title = "fungal Shannon diversity by substrate & dominant tree Schorfheide exploratory") +
  scale_colour_manual(values = c("green2","darkgreen")) +
  geom_violin() + scale_fill_manual(values=c("green2", "darkgreen")) +
  geom_jitter(aes(), height = 0, width = .2, shape=16) 

###20220802 - lists of alpha diversity
shannon_Alb_b_F <- estimate_richness(physeqAlb_bF1, split = T, measures = "Shannon")
shannon_Alb_b_F
write.csv(shannon_Alb_b_F, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Paper_Masterthesis/shannon_Alb_b_F.csv")

shannon_Alb_b_P <- estimate_richness(physeqAlb_bP1, split = T, measures = "Shannon")
write.csv(shannon_Alb_b_P, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Paper_Masterthesis/shannon_Alb_b_P.csv")

shannon_Alb_s_F <- estimate_richness(physeqAlb_sF1, split = T, measures = "Shannon")
write.csv(shannon_Alb_s_F, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Paper_Masterthesis/shannon_Alb_s_F.csv")

shannon_Alb_s_P <- estimate_richness(physeqAlb_sP1, split = T, measures = "Shannon")
write.csv(shannon_Alb_s_P, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Paper_Masterthesis/shannon_Alb_s_P.csv")

shannon_Schorf_b_F <- estimate_richness(physeqSchorf_bF1, split = T, measures = "Shannon")
write.csv(shannon_Schorf_b_F, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Paper_Masterthesis/shannon_Schorf_b_F.csv")

shannon_Schorf_b_P <- estimate_richness(physeqSchorf_bP1, split = T, measures = "Shannon")
write.csv(shannon_Schorf_b_P, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Paper_Masterthesis/shannon_Schorf_b_P.csv")

shannon_Schorf_s_F <- estimate_richness(physeqSchorf_sF1, split = T, measures = "Shannon")
write.csv(shannon_Schorf_s_F, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Paper_Masterthesis/shannon_Schorf_s_F.csv")

shannon_Schorf_s_P <- estimate_richness(physeqSchorf_sP1, split = T, measures = "Shannon")
write.csv(shannon_Schorf_s_P, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Paper_Masterthesis/shannon_Schorf_s_P.csv")


################taxonomy exploration#################

############Alb + Schorfheide################
##heatmaps
#by phylum
plot_heatmap(physeq5, taxa.label="Phylum")
plot_heatmap(physeq5, taxa.label="Class")


######for all EP plots#######

###test statistics for boxplots
#unwrap data from phyloseq object
library("phyloseq")
results = estimate_richness(physeq4, measures = "Shannon")
d = sample_data(physeq4)
#subsets for comparison factors
#statistic diversity bark - soil
bark = results[d[,"substrate"] == "bark"]
soil = results[d[,"substrate"] == "soil"]
wilcox.test(bark, soil)
#p < 2.2e-16
#effect size of Wilcoxon test (Björn Walter)
  # N = 305 correct number of observations?
z <- qnorm(2.2e-16)
r <- z/sqrt(305)
r
  # r= -0.4653514 absolute value < 0.5 -> medium effect size

t.test(bark, soil)
#statistic diversity tree_type
coniferous = results[d[,"tree_type"] == "coniferous"]
deciduous = results[d[,"tree_type"] == "deciduous"]
wilcox.test(coniferous, deciduous)
t.test(coniferous, deciduous)
#statistics for dominant_tree
install.packages("psych")
library(psych)
describeBy(results$Shannon,d$dominant_tree)
anova_dominant_tree <- aov(results$Shannon~d$dominant_tree)
summary(anova_dominant_tree)
#post-hoc analysis to determine between which groups difference is significant
# - pairwise t-test
pairwise.t.test(results$Shannon, d$dominant_tree, p.adjust="bonferroni")
#statistics for exploratory
describeBy(results$Shannon,d$exploratory)
anova_exploratory <- aov(results$Shannon~d$exploratory)
summary(anova_exploratory)
#post-hoc analysis to determine between which groups difference is significant
# - pairwise t-test
pairwise.t.test(results$Shannon, d$exploratory, p.adjust="bonferroni")


##subsample data to calculate test statistics
#subsample reads
(ps_rare <- phyloseq::rarefy_even_depth(physeq4, rngseed = 123, replace = FALSE))
#generate data.frame with adiv measures
adiv <- data.frame(
  "Shannon" = phyloseq::estimate_richness(ps_rare, measures = Shannon),
  "PD" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(ps_rare)))), tree = phyloseq::phy_tree(ps_rare))[, 1],
  "substrate" = phyloseq::sample_data(ps_rare)$substrate,
  "tree_type" = phyloseq::sample_data(ps_rare)$tree_type,
  "dominant_tree" = phyloseq::sample_data(ps_rare)$dominant_tree,
  "exploratory" = phyloseq::sample_data(ps_rare)$exploratory,)
head(adiv)


#das klappt noch nicht mit wrapper function, aber oben i.O.
#my_richness_boxplot <- function(myphyseq, myxcol, mytitle, myxlab)  {
  #phyloseq::plot_richness(physeq4 = myphyseq, x = {{myxcol}}, measures = c('Shannon', 'Observed')) +
   # geom_boxplot(aes(fill= intensity)) +
    #theme_bw() +
    #xlab(myxlab) +
    #ylab('Alpha Diversity') +
    #ggtitle(mytitle) +
    #theme(axis.text=element_text(size=8),axis.title = element_text(size=8),
     #     legend.position = 'none',axis.text.x = element_text(angle = 90))
#}


##phylogenetic tree#
plot_tree(physeq4, color = "tree_type", label.tips = "Phylum", ladderize = "left", justify = "left")
#man sieht leider nix 

#######barplots####
plot_bar(physeq3, x="substrate", fill="Class")
plot_bar(physeq3, x="substrate", fill="Phylum")
plot_bar(physeq3, x="leaf_shape", fill="Phylum")

######################################################
############beta diversity (ordination)#################
######################################################

############all exploratories#################

library("phyloseq")
library("ggplot2")
#ggplot package theme set
theme_set(theme_bw())
library("plyr")

#remove ASVs that do not appear more than 5 times in more than half the samples
#wh0 = genefilter_sample(physeq3, filterfun_sample(function(x) x>5), A=0.5*nsamples(physeq3))
#GP1 = prune_taxa(wh0, physeq3) #GP1 naming from tutorial
#transform to even sampling depth
#GP1 = transform_sample_counts(GP1, function(x) 1E6 * x/sum(x))
#keep only most abundant five phyla
#phylum.sum = tapply(taxa_sums(GP1), tax_table(GP1)[, "Phylum"], sum, na.rm=TRUE)
#top5phyla = names(sort(phylum.sum, TRUE))[1:5]
#GP1 = prune_taxa((tax_table(GP1)[, "Phylum"] %in% top5phyla), GP1)

#plot ordination of ASVs by phylum (NMDS)
ASVord <- ordinate(physeq3, "NMDS", "bray")
p1 = plot_ordination(physeq3, ASVord, type = "taxa", color = "Phylum", title = "taxa")
print(p1)
#facetting of the plot
p1 + facet_wrap(~Phylum, 3)

#ordination by sample types
#by leaf_shape
p2 = plot_ordination(physeq3, ASVord, type="sample", color="leaf_shape")
print(p2)
#p2 + geom_polygon(aes(fill="leaf_shape")) + geom_point(size=5) + ggtitle("samples by leaf_shape")
# sieht mit polygon kacke aus
#by exploratory
p3 = plot_ordination(physeq3, ASVord, type="sample", color="exploratory")
print(p3)
#p3 + geom_polygon(aes(fill="exploratory")) + geom_point(size=5) + ggtitle("samples by exploratory")
#by substrate
p4 = plot_ordination(physeq3, ASVord, type="sample", color="substrate")
print(p4)
#p4 + geom_polygon(aes(fill="substrate")) + geom_point(size=5) + ggtitle("samples by substrate")
p4 + ggtitle("samples by substrate")
#by dominant_tree
p5 = plot_ordination(physeq3, ASVord, type="sample", color="dominant_tree")
p5 + ggtitle("samples by dominant tree")

#biplot ordination graphic with sample and ASV information
p6 = plot_ordination(physeq3, ASVord, type="biplot", color="substrate", shape="Phylum", 
                     title="biplot ordination substrate + Phylum")
print(p6)
# Some stuff to modify the automatic shape scale
physeq3.shape.names = get_taxa_unique(physeq3, "Phylum")
physeq3.shape <- 15:(15 + length(physeq3.shape.names) - 1)
names(physeq3.shape) <- physeq3.shape.names
physeq3.shape["samples"] <- 16
p6 + scale_shape_manual(values=physeq3.shape)
#macht so keinen Sinn
p7 = plot_ordination(physeq3, ASVord, type="split", color="Phylum", shape="substrate", label="exploratory", title="split") 
p7

##ordination PCA/RDA##
ASVord2 <- ordinate(physeq4, "PCoA", "bray")
# ordination of all ASVs by Phylum
p8 = plot_ordination(physeq4, ASVord2, type = "taxa", color = "Phylum", title = "ordination by Phylum - PCoA")
print(p8) # difference to p1 ???
p8 + facet_wrap(~"substrate")

p9 = plot_ordination(physeq4, ASVord2, type="taxa", color="substrate", title="ordination by tree type - PCoA") 
p9 + facet_wrap(~Phylum, 3)

#comparison or ordination methods
dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, physeq, dist){
  ordi = ordinate(physeq4, method=i, distance=dist)
  plot_ordination(physeq4, ordi, "samples", color="SampleType")
}, ASVord2, dist)

names(plist) <- ord_meths

#########################Alb############################

#PCA (RDA) ordination
ordination_A <- ordinate(physeqAlb, "RDA")
ordination_A_rel <- ordinate(ps_rel_abund_A, "RDA")

#ordinations look different depending on rel. abundance or counts

#ordination of ASVs by Phylum -> no clustering
p1 = plot_ordination(physeqAlb, ordination_A, type = "taxa", color = "Phylum", title = "taxa")
p1 = plot_ordination(ps_rel_abund_A, ordination_A_rel, type = "taxa", color = "Phylum", title = "taxa")
print(p1)

#ordination by dominant_tree -> no clustering
p2 = plot_ordination(physeqAlb, ordination_A, type="sample", color="dominant_tree")
p2 = plot_ordination(ps_rel_abund_A, ordination_A_rel, type="sample", color="dominant_tree")
print(p2)

#ordination by substrate -> clustering
p3 = plot_ordination(physeqAlb, ordination_A, type="sample", color="substrate")
p3 = plot_ordination(ps_rel_abund_A, ordination_A_rel, type="sample", color="substrate")
print(p3)

#ordination by ForMIclass -> no clustering
p4 = plot_ordination(physeqAlb, ordination_A, type="sample", color="ForMIclass")
p4 = plot_ordination(ps_rel_abund_A, ordination_A_rel, type="sample", color="ForMIclass")
print(p4)

###approach Ollberding (preferred)#

##ordination with PCA -> CLR transformation -> values are not counts but dominance
ps_clr_A <- microbiome::transform(physeqAlb, "clr")  
phyloseq::otu_table(ps_clr_A)[1:5,1:5]
#PCA via phyloseq
ordination_clr_A <- phyloseq::ordinate(ps_clr_A, "RDA")

################################
#20220517 Ordination as NMDS
ordination_A_NMDS <- ordinate(physeqAlb, "NMDS")
ordination_clr_NMDS <- phyloseq::ordinate(ps_clr_A, "NMDS")
phyloseq::plot_ordination(physeqAlb, ordination_A_NMDS, type="samples", color="dominant_tree", shape="substrate") + 
  geom_jitter(size = 4) +
  stat_ellipse(aes(group = dominant_tree), linetype = 2) +
  scale_colour_manual(values = c("green","darkgreen"), name = "dominant tree species", 
                      labels = c("Fagus sylvatica", "Picea abies"))+
  theme(legend.text = element_text(face = "italic", size = 15), legend.title = element_text(size = 15, face = "bold")) +
  theme(legend.position = "bottom",legend.direction = "vertical")+
  theme(legend.spacing = unit(2,"cm")) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 15))


ordination_S_NMDS <- ordinate(physeqSchorf, "NMDS")
phyloseq::plot_ordination(physeqSchorf, ordination_S_NMDS, type="samples", color="dominant_tree", shape="substrate") + 
  geom_point(size = 4) +
  stat_ellipse(aes(group = dominant_tree), linetype = 2) +
  scale_colour_manual(values = c("green","darkgreen"), name = "dominant tree species", 
                      labels = c("Fagus sylvatica", "Pinus sylvestris"))+
  theme(legend.text = element_text(face = "italic", size = 15), legend.title = element_text(size = 15, face = "bold")) +
  theme(legend.position = "bottom",legend.direction = "vertical")+
  theme(legend.spacing = unit(2,"cm")) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 15))

#20220726 - change ordination color
phyloseq::plot_ordination(physeqAlb, ordination_A_NMDS, type="samples", color="dominant_tree", shape="substrate") + 
  geom_jitter(size = 4) +
  stat_ellipse(aes(group = dominant_tree), linetype = 2) +
  scale_colour_manual(values = c("black","azure4"), name = "dominant tree species", 
                      labels = c("Fagus sylvatica", "Picea abies"))+
  theme(legend.text = element_text(face = "italic", size = 20), legend.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom",legend.direction = "vertical")+
  theme(legend.spacing = unit(2,"cm")) +
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 20))

phyloseq::plot_ordination(physeqSchorf, ordination_S_NMDS, type="samples", color="dominant_tree", shape="substrate") + 
  geom_point(size = 4) +
  stat_ellipse(aes(group = dominant_tree), linetype = 2) +
  scale_colour_manual(values = c("black","azure4"), name = "dominant tree species", 
                      labels = c("Fagus sylvatica", "Pinus sylvestris"))+
  theme(legend.text = element_text(face = "italic", size = 20), legend.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom",legend.direction = "vertical")+
  theme(legend.spacing = unit(2,"cm")) +
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 20))

#######################################
###20220808 - variance partitioning#####
########################################

#define extract OTU table function
veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}

###Swabian Alb###

#extract OTU table Alb
Y <- veganotu(physeqAlb)
#extract sample data and remove sample ID (sample)
dat <- data.frame(sample_data(physeqAlb))
dat <- dat[,-1]
dat

######approach davidzeleny#####

##varpar Alb
#fractions dominant_tree + substrate Alb
rda.all <- rda(Y ~ substrate + dominant_tree, data = dat)
#fraction dominant_tree
rda.tree <- rda(Y ~ dominant_tree, data = dat)
#fraction substrate Alb
rda.substrate <- rda(Y ~ substrate, data = dat)
#calculate variations Alb
# variation substrate + dominant_tree Alb
RsquareAdj(rda.all)
#$r.squared
#[1] 0.100239
#$adj.r.squared
#[1] 0.0802443
tree_substrate_Alb <- RsquareAdj(rda.all)$adj.r.squared
#variation dominant_tree Alb
RsquareAdj(rda.tree)
#$r.squared
#[1] 0.03106127
#$adj.r.squared
#[1] 0.0204136
tree_Alb <- RsquareAdj(rda.tree)$adj.r.squared
# variation substrate Alb
RsquareAdj(rda.substrate)
#$r.squared
#[1] 0.06912488
#$adj.r.squared
#[1] 0.05889548
substrate_Alb <- RsquareAdj(rda.substrate)$adj.r.squared

#shared variation Alb
shared_var <- substrate_Alb + tree_Alb - tree_substrate_Alb #-0.0009352148 - negative value possible? Betrag?
shared_var_2 <- tree_substrate_Alb - substrate_Alb - tree_Alb
# variation substrate Alb
substrate_var <- substrate_Alb - shared_var_2 # 0.05796027
# variation tree Alb
tree_var <- tree_Alb - shared_var_2 # 0.01947838

# variation with varpart function Alb
varp <- varpart(Y, ~ substrate, ~ dominant_tree, data = dat)
varp

#plot venn diagramm variation Alb
plot(varp, digits = 2, Xnames = c("substrate", "tree species"), bg = c("navy", "hotpink"))

#anova testing Alb
#fraction substrate Alb
rda.substrate.tree <- rda(Y ~ substrate + Condition (dominant_tree), data = dat)
#fraction tree Alb
rda.tree.substrate <- rda(Y ~ dominant_tree + Condition (substrate), data = dat)
#test individual fractions Alb
anova(rda.all) #***
anova(rda.substrate)#***
anova(rda.tree) #***
anova(rda.substrate.tree) #***
anova(rda.tree.substrate) #***

### Schorfheide ###

#extract OTU table Schorf
Z <- veganotu(physeqSchorf)
#extract sample data and remove sample ID (sample)
dat2 <- data.frame(sample_data(physeqSchorf))
dat2 <- dat2[,-1]
dat2

######approach davidzeleny#####

##varpar Schorf
#fractions dominant_tree + substrate Schorf
rda.all_Schorf <- rda(Z ~ substrate + dominant_tree, data = dat2)
#fraction dominant_tree Schorf
rda.tree_Schorf <- rda(Z ~ dominant_tree, data = dat2)
#fraction substrate Schorf
rda.substrate_Schorf <- rda(Z ~ substrate, data = dat2)
#calculate variations Schorf
# variation substrate + dominant_tree Schorf
RsquareAdj(rda.all_Schorf)
#$r.squared
#[1] 0.1790991
#$adj.r.squared
#[1] 0.1595538
tree_substrate_Schorf <- RsquareAdj(rda.all_Schorf)$adj.r.squared
#variation dominant_tree Schorf
RsquareAdj(rda.tree_Schorf)
#$r.squared
#[1] 0.02713039
#$adj.r.squared
#[1] 0.01568487
tree_Schorf <- RsquareAdj(rda.tree_Schorf)$adj.r.squared
# variation substrate Schorf
RsquareAdj(rda.substrate_Schorf)
#$r.squared
#[1] 0.1519012
#$adj.r.squared
#[1] 0.1419236
substrate_Schorf <- RsquareAdj(rda.substrate_Schorf)$adj.r.squared

#shared variation Schorf
shared_var_Schorf <- substrate_Schorf + tree_Schorf - tree_substrate_Schorf #-0.001945349
# variation substrate Schorf
substrate_var_Schorf <- substrate_Schorf - shared_var_Schorf # 0.143869
# variation tree Schorf
tree_var_Schorf <- tree_Schorf - shared_var_Schorf # 0.01763022

# variation with varpart function Schorf
varp_Schorf <- varpart(Z, ~ substrate, ~ dominant_tree, data = dat2)
varp_Schorf

#plot venn diagramm variation Schorf
plot(varp_Schorf, digits = 2, Xnames = c("substrate", "tree species"), bg = c("navy", "hotpink"))

#anova testing Schorf
#fraction substrate Schorf
rda.substrate.tree_Schorf <- rda(Z ~ substrate + Condition (dominant_tree), data = dat2)
#fraction tree Schorf
rda.tree.substrate_Schorf <- rda(Z ~ dominant_tree + Condition (substrate), data = dat2)
#test individual fractions
anova(rda.all_Schorf) #***
anova(rda.substrate_Schorf)#***
anova(rda.tree_Schorf) #** - p 0.008
anova(rda.substrate.tree_Schorf) #***
anova(rda.tree.substrate_Schorf) #***


# use standardised categorical and continuous predictors in the constraint (X)
X <- ade4::dudi.mix(dat, scannf=F, nf=2)$tab

library(devtools)
install_github("umerijaz/microbiomeSeq")
library(microbiomeSeq)
p <- plot_anova_env(physeqAlb, grouping_column = "substrate", pValueCutoff = 0.05, 
                    select.variables = "dominant_tree")
print(p)


##determine the influence of PCs
#Plot scree plot -> shows proportion of variance per PC
phyloseq::plot_scree(ordination_clr_A) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
#Examine eigenvalues and % prop. variance explained
head(ordination_clr_A$CA$eig)
sapply(ordination_clr_A$CA$eig[1:5], function(x) x / sum(ordination_clr_A$CA$eig))

#Scale axes
clr1A <- ordination_clr_A$CA$eig[1] / sum(ordination_clr_A$CA$eig)
clr2A <- ordination_clr_A$CA$eig[2] / sum(ordination_clr_A$CA$eig)

#Generate distance matrix
clr_dist_matrix_A <- phyloseq::distance(ps_clr_A, method = "euclidean") 

#ordination by substrate
phyloseq::plot_ordination(physeqAlb, ordination_clr_A, type="samples", color="substrate", shape="dominant_tree") + 
  geom_point(size = 2) +
  coord_fixed(clr2A / clr1A) +
  stat_ellipse(aes(group = substrate), linetype = 2)

#ADONIS test
vegan::adonis(clr_dist_matrix_A ~ phyloseq::sample_data(ps_clr_A)$substrate)
vegan::adonis(clr_dist_matrix_A ~ phyloseq::sample_data(ps_clr_A)$substrate * dominant_tree)
#Dispersion test and plot -> plot with distances to centroid
dispr_A <- vegan::betadisper(clr_dist_matrix_A, phyloseq::sample_data(ps_clr_A)$substrate)
dispr_A
plot(dispr_A, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
#boxplot distance to centroid
boxplot(dispr_A, main ="", xlab = "")
#permutation test
permutest(dispr_A)

#ordination by dominant_tree
phyloseq::plot_ordination(physeqAlb, ordination_clr_A, type="samples", color="dominant_tree", shape="substrate") + 
  geom_point(size = 2) +
  coord_fixed(clr2A / clr1A) +
  stat_ellipse(aes(group = dominant_tree), linetype = 2) +
  scale_colour_manual(values = c("green","darkgreen"), name = "dominant tree species", 
                      labels = c("Fagus sylvatica", "Picea abies"))+
  theme(legend.text = element_text(face = "italic", size = 12), legend.title = element_text(size = 12, face = "bold")) +
  theme(legend.position = "bottom",legend.direction = "vertical")+
  theme(legend.spacing = unit(2,"cm"))
ggsave(
  "PCA_Alb.png",
  width = 8.3,
  height = 5.7,
  dpi = 800
)  


#ADONIS test
vegan::adonis(clr_dist_matrix_A ~ phyloseq::sample_data(ps_clr_A)$dominant_tree)
#Dispersion test and plot -> plot with distances to centroid
dispr_A <- vegan::betadisper(clr_dist_matrix_A, phyloseq::sample_data(ps_clr_A)$dominant_tree)
dispr_A
plot(dispr_A, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
#boxplot distance to centroid
boxplot(dispr_A, main ="", xlab = "")
#permutation test
permutest(dispr_A)

#ordination by ForMIclass -> no clustering
phyloseq::plot_ordination(physeqAlb, ordination_clr_A, type="samples", color="ForMIclass") + 
  geom_point(size = 2) +
  coord_fixed(clr2A / clr1A) +
  stat_ellipse(aes(group = ForMIclass), linetype = 2)

##phylogenetic information in beta diversity - PCoA analysis
#generate distances
ord_unifrac <- ordinate(physeqAlb, method = "PCoA", distance = "wunifrac") 
ord_unifrac_un <- ordinate(physeqAlb, method = "PCoA", distance = "unifrac")
#Plot ordinations
a <- plot_ordination(physeqAlb, ord_unifrac, color = "substrate") + geom_point(size = 2)
b <- plot_ordination(physeqAlb, ord_unifrac_un, color = "substrate") + geom_point(size = 2)
cowplot::plot_grid(a, b, nrow = 1, ncol = 2, scale = .9, labels = c("Weighted", "Unweighted"))
# idea: create cowplot for substrate + dominant_tree


######################Schorfheide#######################

#PCA (RDA) ordination
ordination_S <- ordinate(physeqSchorf, "RDA")
ordination_S_rel <- ordinate(ps_rel_abund_S, "RDA")

#ordinations look different depending on rel. abundance or counts

#ordination of ASVs by Phylum -> no clustering
#p5 = plot_ordination(physeqSchorf, ordination_S, type = "taxa", color = "Phylum", title = "taxa")
p5 = plot_ordination(ps_rel_abund_S, ordination_S_rel, type = "taxa", color = "Phylum", title = "taxa")
print(p5)

#ordination by dominant_tree -> no clustering
#p6 = plot_ordination(physeqSchorf, ordination_S, type="sample", color="dominant_tree")
p6 = plot_ordination(ps_rel_abund_S, ordination_S_rel, type="sample", color="dominant_tree")
print(p6)

#ordination by substrate -> clustering
#p7 = plot_ordination(physeqSchorf, ordination_S, type="sample", color="substrate")
p7 = plot_ordination(ps_rel_abund_S, ordination_S_rel, type="sample", color="substrate")
print(p7)

#ordination by ForMIclass -> no clustering
#p8 = plot_ordination(physeqSchorf, ordination_S, type="sample", color="ForMIclass")
p8 = plot_ordination(ps_rel_abund_S, ordination_S_rel, type="sample", color="ForMIclass")
print(p8)

###approach Ollberding (preferred)#

##ordination with PCA -> CLR transformation -> values are not counts but dominance
ps_clr_S <- microbiome::transform(physeqSchorf, "clr")  
phyloseq::otu_table(ps_clr_S)[1:5,1:5]
#PCA via phyloseq
ordination_clr_S <- phyloseq::ordinate(ps_clr_S, "RDA")

##determine the influence of PCs
#Plot scree plot -> shows proportion of variance per PC
phyloseq::plot_scree(ordination_clr_S) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
#Examine eigenvalues and % prop. variance explained
head(ordination_clr_S$CA$eig)
sapply(ordination_clr_S$CA$eig[1:5], function(x) x / sum(ordination_clr_S$CA$eig))

#Scale axes
clr1S <- ordination_clr_S$CA$eig[1] / sum(ordination_clr_S$CA$eig)
clr2S <- ordination_clr_S$CA$eig[2] / sum(ordination_clr_S$CA$eig)

#Generate distance matrix
clr_dist_matrix_S <- phyloseq::distance(ps_clr_S, method = "euclidean")

#ordination by substrate
phyloseq::plot_ordination(physeqSchorf, ordination_clr_S, type="samples", color="substrate", shape="dominant_tree") + 
  geom_point(size = 2) +
  coord_fixed(clr2S / clr1S) +
  stat_ellipse(aes(group = substrate), linetype = 2)

#ADONIS test
vegan::adonis(clr_dist_matrix_S ~ phyloseq::sample_data(ps_clr_S)$substrate)
#Dispersion test and plot -> plot with distances to centroid
dispr_S <- vegan::betadisper(clr_dist_matrix_S, phyloseq::sample_data(ps_clr_S)$substrate)
dispr_S
plot(dispr_S, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
#boxplot distance to centroid
boxplot(dispr_S, main ="", xlab = "")
#permutation test
permutest(dispr_S)

#ordination by dominant_tree
phyloseq::plot_ordination(physeqSchorf, ordination_clr_S, type="samples", color="dominant_tree", shape="substrate") + 
  geom_point(size = 2) +
  coord_fixed(clr2S / clr1S) +
  stat_ellipse(aes(group = dominant_tree), linetype = 2) +
  scale_colour_manual(values = c("green","darkgreen"), name = "dominant tree species",
                      labels = c("Fagus sylvatica", "Pinus sylvestris")) +
  theme(legend.text = element_text(face = "italic", size = 12), legend.title = element_text(size = 12, face = "bold")) +
  theme(legend.position = "bottom",legend.direction = "vertical")+
  theme(legend.spacing = unit(2,"cm"))
ggsave(
  "PCA_Schorf.png",
  width = 8.3,
  height = 5.7,
  dpi = 800
)  

#ADONIS test
vegan::adonis(clr_dist_matrix_S ~ phyloseq::sample_data(ps_clr_S)$dominant_tree)
#Dispersion test and plot -> plot with distances to centroid
dispr_S <- vegan::betadisper(clr_dist_matrix_S, phyloseq::sample_data(ps_clr_S)$dominant_tree)
dispr_S
plot(dispr_S, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
#boxplot distance to centroid
boxplot(dispr_S, main ="", xlab = "")
#permutation test
permutest(dispr_S)

#ordination by ForMIclass -> no clustering
phyloseq::plot_ordination(physeqSchorf, ordination_clr_S, type="samples", color="ForMIclass") + 
  geom_point(size = 2) +
  coord_fixed(clr2S / clr1S) +
  stat_ellipse(aes(group = ForMIclass), linetype = 2)

##phylogenetic information in beta diversity - PCoA analysis
#generate distances
ord_unifrac <- ordinate(physeqSchorf, method = "PCoA", distance = "wunifrac") 
ord_unifrac_un <- ordinate(physeqSchorf, method = "PCoA", distance = "unifrac")
#Plot ordinations
a <- plot_ordination(physeqSchorf, ord_unifrac, color = "substrate") + geom_point(size = 2)
b <- plot_ordination(physeqSchorf, ord_unifrac_un, color = "substrate") + geom_point(size = 2)
cowplot::plot_grid(a, b, nrow = 1, ncol = 2, scale = .9, labels = c("Weighted", "Unweighted"))
# idea: create cowplot for substrate + dominant_tree

######################################################
############overlap analysis - Venn diagrams#################
######################################################

##########soil & bark per exploratory##################

ps_venn(physeqAlb, group = "dominant_tree", fraction = 0, weight = TRUE, relative = TRUE, plot = TRUE)
ps_venn(physeqSchorf, group = "dominant_tree", fraction = 0, weight = FALSE, relative = TRUE, plot = TRUE)

#####################Alb#######################

#plot soil Alb
ps_venn(physeqAlb_s, group = "dominant_tree", fraction = 0, weight = TRUE, relative = TRUE, plot = TRUE,
        fill = c("green","mediumseagreen"),labels = list(labels =c("Fagus sylvatica", "Picea abies"),cex=1.6,font=list(face=3)),
        quantities = list(type=c("percent"), labels = c("\n19% (3592)","\n9% (1198)","\n71% (1277)"), cex=1.6))
#png(filename="venn_Alb_soil.png", res=800, width = 1600, height = 900)

#list soil Alb
ps_venn(physeqAlb_s, group = "dominant_tree", fraction = 0, weight = FALSE, relative = TRUE, plot = FALSE)

#plot bark Alb
ps_venn(physeqAlb_b, group = "dominant_tree", fraction = 0, weight = TRUE, relative = TRUE, plot = TRUE,
        fill = c("green","mediumseagreen"), labels = list(labels =c("Fagus sylvatica", "Picea abies"),cex=1.6,font=list(face=3)),
        quantities = list(type=c("percent"), labels = c("\n11% (697)","\n7% (341)","\n82% (279)"), cex=1.6))
#list bark Alb
ps_venn(physeqAlb_b, group = "dominant_tree", fraction = 0, weight = FALSE, relative = TRUE, plot = FALSE)

#combined by venn_class
#plot venn_class Alb
ps_venn(physeqAlb, group = "venn_class", fraction = 0, weight = FALSE, relative = TRUE, plot = TRUE,
        main = "Venn Diagram - F.sylvatica & P.abies - Alb", fill = c("green","darkolivegreen2","darkgreen","darkolivegreen"))
ps_venn(physeqAlb, group = "venn_class", fraction = 0, weight = FALSE, relative = TRUE, plot = TRUE,
        main = "Venn Diagram - F.sylvatica & P.abies - Alb", fill = c("green","darkolivegreen2","darkgreen","darkolivegreen"))


###lists of venn_class Alb
shared_taxa_A <-ps_venn(physeqAlb, group = "venn_class", fraction = 0, weight = TRUE, relative = TRUE, plot = FALSE,
                      quantities = list(type=c("percent")))
fagus_bark_only_A <- data.frame("ASV_ID" = shared_taxa_A$F.sylvatica_bark)
fagus_soil_only_A <- data.frame("ASV_ID" = shared_taxa_A$F.sylvatica_soil)
fagus_soil_bark_A <- data.frame("ASV_ID" = shared_taxa_A$F.sylvatica_bark__F.sylvatica_soil)
fagus_specialists_A <- rbind(fagus_bark_only_A,fagus_soil_only_A,fagus_soil_bark_A)
fagus_specialists_A
picea_bark_only_A <- data.frame("ASV_ID"= shared_taxa_A$P.abies_bark)
picea_soil_only_A <- data.frame("ASV_ID"=shared_taxa_A$P.abies_soil)
picea_soil_bark_A <- data.frame("ASV_ID" =shared_taxa_A$F.sylvatica_soil__P.abies_bark__P.abies_soil)
picea_specialists_A <- rbind(picea_bark_only_A, picea_soil_only_A, picea_soil_bark_A)
overlap_all_A <- data.frame("ASV_ID" = shared_taxa_A$F.sylvatica_bark__F.sylvatica_soil__P.abies_bark__P.abies_soil)
tax_Alb <- data.frame(tax_table(physeqAlb))
tax_Alb
tax_Alb1 <- rownames_to_column(tax_Alb, var = "ASV_ID") %>% as_tibble()
tax_Alb1
generalists_Alb <- left_join(overlap_all_A,tax_Alb1, by = "ASV_ID")
generalists_Alb
write.csv(generalists_Alb, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/generalists_Alb.csv")
specialists_fagus_Alb <- left_join(fagus_specialists_A,tax_Alb1, by = "ASV_ID")
write.csv(specialists_fagus_Alb, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/specialists_fagus_Alb.csv")
specialists_picea_Alb <- left_join(picea_specialists_A, tax_Alb1, by = "ASV_ID")
write.csv(specialists_picea_Alb, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/specialists_picea_Alb.csv")

#count data + relative abundance
ps_venn(physeqAlb, group = "venn_class", fraction = 0, weight = TRUE, relative = TRUE, plot = TRUE,
        main = "Venn Diagram - F.sylvatica & P.abies - Alb", fill = c("green","darkolivegreen2","darkgreen","darkolivegreen"),
        quantities = list(type=c("percent")))

#count data + relative abundance (manual)
ps_venn(physeqAlb, group = "venn_class", fraction = 0, weight = TRUE, relative = TRUE, plot = TRUE,
        main = "Venn Diagram - F.sylvatica & P.abies - Alb \n weighted by abundance", fill = c("green","darkolivegreen2","darkgreen","darkolivegreen"),
        quantities = list(type=c("percent"), labels = c("4% (478)","7% (3405)","3% (256)","4% (1137)",
                                                        "3% (123)","8% (134)","<1% (8)","<1% (22)","17% (1076)","<1% (25)","12% (42)",
                                                        "7% (88)","6% (28)","2% (38)","26% (75)")))


#####################Schorfheide#######################

#plot soil Schorfheide
ps_venn(physeqSchorf_s, group = "dominant_tree", fraction = 0, weight = TRUE, relative = TRUE, plot = TRUE,
        main = "Venn Diagram - soil - Schorfheide", fill = c("green","darkgreen"))
ps_venn(physeqSchorf_s, group = "dominant_tree", fraction = 0, weight = TRUE, relative = TRUE, plot = TRUE,
        fill = c("green","mediumseagreen"), labels = list(labels =c("Fagus sylvatica", "Pinus sylvestris"),cex=1.6,font=list(face=3)),
        quantities = list(type=c("percent"), labels= c("\n8% (1412)", "\n4% (814)", "\n88% (885)"),cex=1.6))
#list soil Schorfheide
ps_venn(physeqSchorf_s, group = "dominant_tree", fraction = 0, weight = FALSE, relative = TRUE, plot = FALSE)

#plot bark Schorfheide
ps_venn(physeqSchorf_b, group = "dominant_tree", fraction = 0, weight = TRUE, relative = TRUE, plot = TRUE,
        fill = c("green","mediumseagreen"),  labels = list(labels =c("Fagus sylvatica", "Pinus sylvestris"),cex=1.6,font=list(face=3)),  
        quantities = list(type=c("percent"),labels = c("\n3% (478)", "\n7% (228)", "\n89% (202)"),cex =1.6))
ps_venn(physeqSchorf_b, group = "dominant_tree", fraction = 0.05, weight = TRUE, relative = TRUE, plot = TRUE,
        main = "Venn Diagram - bark - Schorfheide", fill = c("green","darkgreen"))
#list bark Schorfheide
ps_venn(physeqSchorf_b, group = "dominant_tree", fraction = 0, weight = FALSE, relative = TRUE, plot = FALSE)

#combined by venn_class
#plot venn_class Schorfheide
ps_venn(physeqSchorf, group = "venn_class", fraction = 0, weight = FALSE, relative = TRUE, plot = TRUE,
        main = "Venn Diagram - F.sylvatica & P.sylvestris - Schorfheide", fill = c("green","darkolivegreen2","darkgreen","darkolivegreen"))
ps_venn(physeqSchorf, group = "venn_class", fraction = 0, weight = TRUE, relative = TRUE, plot = TRUE,
        main = "Venn Diagram - F.sylvatica & P.sylvestris - Schorfheide", fill = c("green","darkolivegreen2","darkgreen","darkolivegreen"),
        quantities = list(type=c("percent")))

###lists of venn_class Alb
shared_taxa_S <-ps_venn(physeqSchorf, group = "venn_class", fraction = 0, weight = TRUE, relative = TRUE, plot = FALSE,
                        quantities = list(type=c("percent")))
shared_taxa_S
fagus_bark_only_S <- data.frame("ASV_ID" = shared_taxa_S$F.sylvatica_bark)
fagus_soil_only_S <- data.frame("ASV_ID" = shared_taxa_S$F.sylvatica_soil)
fagus_soil_bark_S <- data.frame("ASV_ID" = shared_taxa_S$F.sylvatica_bark__F.sylvatica_soil)
fagus_specialists_S <- rbind(fagus_bark_only_S,fagus_soil_only_S,fagus_soil_bark_S)
fagus_specialists_S
pinus_bark_only_S <- data.frame("ASV_ID"= shared_taxa_S$P.sylvestris_bark)
pinus_soil_only_S <- data.frame("ASV_ID"=shared_taxa_S$P.sylvestris_soil)
pinus_soil_bark_S <- data.frame("ASV_ID" =shared_taxa_S$P.sylvestris_bark__P.sylvestris_soil)
pinus_specialists_S <- rbind(pinus_bark_only_S, pinus_soil_only_S, pinus_soil_bark_S)
overlap_all_S <- data.frame("ASV_ID" = shared_taxa_S$F.sylvatica_bark__F.sylvatica_soil__P.sylvestris_bark__P.sylvestris_soil)
overlap_all_S
tax_Schorf <- data.frame(tax_table(physeqSchorf))
tax_Schorf
tax_Schorf1 <- rownames_to_column(tax_Schorf, var = "ASV_ID") %>% as_tibble()
tax_Schorf1
generalists_Schorf <- left_join(overlap_all_S,tax_Schorf1, by = "ASV_ID")
generalists_Schorf
write.csv(generalists_Schorf, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/generalists_Schorf.csv")
specialists_fagus_Schorf <- left_join(fagus_specialists_S,tax_Schorf1, by = "ASV_ID")
write.csv(specialists_fagus_Schorf, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/specialists_fagus_Schorf.csv")
specialists_pinus_Schorf <- left_join(pinus_specialists_S, tax_Schorf1, by = "ASV_ID")
write.csv(specialists_pinus_Schorf, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/specialists_pinus_Schorf.csv")


#count data + relative abundance (manual)
# plot layout: width 800, height 900
ps_venn(physeqSchorf, group = "venn_class", fraction = 0, weight = TRUE, relative = TRUE, plot = TRUE,
        main = "Venn Diagram - F.sylvatica & P.sylvestris - Schorfheide \n weighted by abundance", 
        fill = c("green","darkolivegreen2","darkgreen","darkolivegreen"),
        quantities = list(type=c("percent"), labels = c("1% (338)","3% (1312)","3% (175)","2% (778)",
                                                        "<1% (62)","10% (111)","3% (13)","<1% (14)","16% (737)","<1% (10)","5% (24)",
                                                        "6% (65)","3% (13)","3% (29)","46% (54)")))

#list venn_class Schorfheide
ps_venn(physeqSchorf, group = "venn_class", fraction = 0, weight = FALSE, relative = TRUE, plot = FALSE)


#####################Alb + Schorfheide#######################
#plot Alb+Schorfheide (all substrates)
ps_venn(physeq5, group = "tree_type", fraction = 0, weight = FALSE, relative = TRUE, plot = TRUE,
        main = "Venn Diagram - soil + bark - Alb + Schorfheide", fill = c("darkgreen","green"))





#######################################################################
#########statistical analysis of microbiome data######################
######################################################################

#install missing packages
.cran_packages <- c("cowplot", "picante", "HMP", "dendextend", "rms", "devtools")
.bioc_packages <- c("DESeq2", "microbiome", "metagenomeSeq", "ALDEx2")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(.bioc_packages)
devtools::install_github("adw96/breakaway")
devtools::install_github(repo = "UVic-omics/selbal")

#loading libraries
library("tidyverse")
library("DESeq2") 
library("microbiome")
library("vegan")
library("picante")
library("ALDEx2")
library("metagenomeSeq")
library("HMP")
library("dendextend")
library("selbal") #installation failed?
library("rms")
library("breakaway")
library("phyloseq")

#color palette setting
install.packages("paletteer")
library("paletteer")
palettes_d_names
paletteer_d("nord::aurora")
my_cols <- paletteer_c("viridis::inferno", n=18, direction = -1)
#"ggsci::default_igv"
#paletteer_c("viridis::inferno", n=18, direction = -1)
#"trekcolors::lcars_series"
# only 5 values in "nord::aurora"
# only 4 values "nord::frost"
box_cols <- c( "#5050FFFF", "#CE3D32FF", "#008099FF")

######visualizing relative abundance##########
#######alpha diversity##############

#e.g. stacked bar plots or faceted boxplots
#count phyla
table(phyloseq::tax_table(physeq4)[,"Phylum"])
#count orders
table(phyloseq::tax_table(physeq4)[,"Order"])
write.csv(table(phyloseq::tax_table(physeq4)[,"Order"]), "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/overview_orders.csv" )
#convert to relative abundance
ps_rel_abund = phyloseq::transform_sample_counts(physeq4, function(x){x / sum(x)})
phyloseq::otu_table(physeq4)[1:5, 1:5] #in this case (choice of display range) 0 counts
phyloseq::otu_table(ps_rel_abund)[1:5, 1:5] #logically same as above

####prepare data for barplots subset Alb + Schorfheide###
ps_rel_abund_AS = phyloseq::transform_sample_counts(physeq5, function(x){x / sum(x)})
phyloseq::otu_table(physeq5)[1:5, 1:5] #in this case (choice of display range) 0 counts
phyloseq::otu_table(ps_rel_abund_AS)[1:5, 1:5] #logically same as above
#by exploratory
phyloseq::plot_bar(ps_rel_abund_AS, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ exploratory, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
#by substrate
phyloseq::plot_bar(ps_rel_abund_AS, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ substrate, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
#by dominant tree
phyloseq::plot_bar(ps_rel_abund_AS, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ dominant_tree, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

################Schorfheide#######################

##order overview Schorfheide##
#count orders
table(phyloseq::tax_table(physeqSchorf)[,"Order"])
write.csv(table(phyloseq::tax_table(physeqSchorf)[,"Order"]), "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/overview_orders_Schorf.csv" )

#####prepare data for barplots subset Schorfheide####
ps_rel_abund_S = phyloseq::transform_sample_counts(physeqSchorf, function(x){x / sum(x)})
phyloseq::otu_table(physeqSchorf)[1:5, 1:5] #in this case (choice of display range) 0 counts
phyloseq::otu_table(ps_rel_abund_S)[1:5, 1:5] #logically same as above
####barplots relative abundance Schorfheide####
  #barplot Phyla by substrate##
phyloseq::plot_bar(ps_rel_abund_S, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ substrate, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  #barplot Phyla by Fagus vs. Pinus (by dominant tree)##
phyloseq::plot_bar(ps_rel_abund_S, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ dominant_tree, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  #TEST Order by Fagus vs. Pinus (by dominant tree)##
  #result: extension of ledgend makes plot unreadable
phyloseq::plot_bar(ps_rel_abund_S, fill = "Order") +
  geom_bar(aes(color = Order, fill = Order), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ dominant_tree, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

######subsets soil & bark Schorfheide###

##soil
#prepare data for barplots
ps_rel_abund_Ss = phyloseq::transform_sample_counts(physeqSchorf_s, function(x){x / sum(x)})
phyloseq::otu_table(physeqSchorf_s)[1:5, 1:5] #in this case (choice of display range) 0 counts
phyloseq::otu_table(ps_rel_abund_Ss)[1:5, 1:5] #logically same as above

#by dominant_tree subset soil
phyloseq::plot_bar(ps_rel_abund_Ss, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ dominant_tree, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

##bark
#prepare data for barplots
ps_rel_abund_Sb = phyloseq::transform_sample_counts(physeqSchorf_b, function(x){x / sum(x)})
phyloseq::otu_table(physeqSchorf_b)[1:5, 1:5] #in this case (choice of display range) 0 counts
phyloseq::otu_table(ps_rel_abund_Sb)[1:5, 1:5] #logically same as above

#by dominant_tree subset bark
phyloseq::plot_bar(ps_rel_abund_Sb, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ dominant_tree, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

######subsets per dominant_tree###

##Fagus Schorfheide
#prepare data for barplots
ps_rel_abund_SF = phyloseq::transform_sample_counts(physeqSchorf_F, function(x){x / sum(x)})
phyloseq::otu_table(physeqSchorf_F)[1:5, 1:5] #in this case (choice of display range) 0 counts
phyloseq::otu_table(ps_rel_abund_SF)[1:5, 1:5] #logically same as above

#by substrate
phyloseq::plot_bar(ps_rel_abund_SF, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ substrate, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

##Pinus Schorfheide
#prepare data for barplots
ps_rel_abund_SP = phyloseq::transform_sample_counts(physeqSchorf_P, function(x){x / sum(x)})
phyloseq::otu_table(physeqSchorf_P)[1:5, 1:5] #in this case (choice of display range) 0 counts
phyloseq::otu_table(ps_rel_abund_SP)[1:5, 1:5] #logically same as above

#by substrate
phyloseq::plot_bar(ps_rel_abund_SP, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ substrate, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

####
##boxplot abundance order level split by substrate Schorfheide
#####

#Schorfheide bark by dominant tree

# aggregate taxa on order level
physeq_order_Schorf_bark <- aggregate_taxa(physeqSchorf_b, level = 'Order')
# calculate relative abundance
ps_rel_abund_ord_Schorf_bark = phyloseq::transform_sample_counts(physeq_order_Schorf_bark, function(x){x / sum(x)})
# determine top 12
order_top12_Sb <- get_top_taxa(ps_rel_abund_ord_Schorf_bark, 12, relative = TRUE, other_label = "Others", discard_other = T)

# boxplot relative abundance dominant tree Schorfheide bark

phyloseq::psmelt(order_top12_Sb) %>%
  ggplot(data = ., aes(x = dominant_tree, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  scale_color_discrete(name = "Top 12 Orders - Schorfheide bark",
                       labels = c("Agaricales","Capnodiales","Chaetothyriales","Helotiales",
                                  "Hypocreales","Lecanorales","Orbiliales","Pleosporales",
                                  "Tremellales","Umbilicariales","Verrucariales","not assigned")) + 
  theme(legend.position = "bottom",legend.direction = "horizontal") +
  theme(axis.text.x= element_text(face = "italic", size = 12)) +
  scale_x_discrete(labels = c("Fagus sylvatica","Pinus sylvestris"))+
  theme( strip.text.x = element_blank() )+
  theme(legend.text = element_text(size = 12), legend.title = element_blank())+
  guides(color = guide_legend(override.aes = list(size=4)))+
  theme(legend.key.width = unit(2.5,"cm")) +
  labs(x = "", y = "relative abundance\n") +
  theme(axis.title.y = element_text(size = 12)) +
  facet_wrap(~ OTU, scales = "free") +
  ylim(0,0.6)
ggsave(
  "order_relab_Schorf_bark_top12.png",
  width = 13,
  height = 9.5,
  dpi = 800
)

# determine top 4
order_top4_Sb <- get_top_taxa(ps_rel_abund_ord_Schorf_bark, 4, relative = TRUE, other_label = "Others", discard_other = T)

# boxplot relative abundance dominant tree Schorfheide bark

phyloseq::psmelt(order_top4_Sb) %>%
  ggplot(data = ., aes(x = dominant_tree, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  scale_color_discrete(name = "Top 12 Orders - Schorfheide bark",
                      labels = c("Capnodiales","Chaetothyriales",
                                  "Lecanorales","not assigned")) + 
  theme(legend.position = "bottom",legend.direction = "horizontal") +
  theme(axis.text.x= element_text(face = "italic", size = 12)) +
  scale_x_discrete(labels = c("Fagus sylvatica","Pinus sylvestris"))+
  theme( strip.text.x = element_blank() )+
  theme(legend.text = element_text(size = 12), legend.title = element_blank())+
  guides(color = guide_legend(override.aes = list(size=4)))+
  theme(legend.key.width = unit(1.2,"cm")) +
  labs(x = "", y = "relative abundance\n") +
  theme(axis.title.y = element_text(size = 12)) +
  facet_wrap(~ OTU, scales = "free") +
  ylim(0,0.6)
ggsave(
  "order_relab_Schorf_bark_top4.png",
  width = 9.5,
  height = 7,
  dpi = 800
)

#average relative abundance top 4 orders Schorfheide bark

# subset of unassigned orders on bark Schorfheide
ps_rel_abund_Schorf_bark_order_na <- subset_taxa(physeqSchorf_b, is.na(Order))
# percentage of unassigned orders on bark Schorfheide
sum(sample_sums(ps_rel_abund_Schorf_bark_order_na))/sum(sample_sums(physeqSchorf_b))

# subset of Capnodiales on bark Schorfheide
ps_rel_abund_Schorf_bark_order_capno <- subset_taxa(physeqSchorf_b, Order == "o__Capnodiales")
# percentage of Capnodiales on bark Schorfheide
sum(sample_sums(ps_rel_abund_Schorf_bark_order_capno))/sum(sample_sums(physeqSchorf_b))

# subset of Lecanorales on bark Schorfheide
ps_rel_abund_Schorf_bark_order_lecano <- subset_taxa(physeqSchorf_b, Order == "o__Lecanorales")
# percentage of Lecanorales on bark Schorfheide
sum(sample_sums(ps_rel_abund_Schorf_bark_order_lecano))/sum(sample_sums(physeqSchorf_b))

# subset of Chaetothyriales on bark Schorfheide
ps_rel_abund_Schorf_bark_order_chaeto <- subset_taxa(physeqSchorf_b, Order == "o__Chaetothyriales")
# percentage of Chaetothyriales on bark Schorfheide
sum(sample_sums(ps_rel_abund_Schorf_bark_order_chaeto))/sum(sample_sums(physeqSchorf_b))

# subset of Trapeliales on bark Schorfheide (to compare to Swabian Alb)
ps_rel_abund_Schorf_bark_order_trapeli <- subset_taxa(physeqSchorf_b, Order == "o__Trapeliales")
# percentage of Trapeliales on bark Schorfheide
sum(sample_sums(ps_rel_abund_Schorf_bark_order_trapeli))/sum(sample_sums(physeqSchorf_b))

#Schorfheide soil by dominant tree

# aggregate taxa on order level
physeq_order_Schorf_soil <- aggregate_taxa(physeqSchorf_s, level = 'Order')
# calculate relative abundance
ps_rel_abund_ord_Schorf_soil = phyloseq::transform_sample_counts(physeq_order_Schorf_soil, function(x){x / sum(x)})

# determine top 12
order_top12_Ss <- get_top_taxa(ps_rel_abund_ord_Schorf_soil, 12, relative = TRUE, other_label = "Others", discard_other = T)

# boxplot relative abundance dominant tree Schorfheide soil

phyloseq::psmelt(order_top12_Ss) %>%
  ggplot(data = ., aes(x = dominant_tree, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  scale_color_discrete(name = "Top 12 Orders - Schorfheide soil",
                       labels = c("Agaricales","Boletales","Eurotiales","Filobasidiales",
                                  "Helotiales","Hypocreales","Mortierellales","Pezizales",
                                  "Russulales","Thelebolales","Tremellales","not assigned")) + 
  theme(legend.position = "bottom",legend.direction = "horizontal") +
  theme(axis.text.x= element_text(face = "italic", size = 12)) +
  scale_x_discrete(labels = c("Fagus sylvatica","Pinus sylvestris"))+
  theme( strip.text.x = element_blank() )+
  theme(legend.text = element_text(size = 12), legend.title = element_blank())+
  guides(color = guide_legend(override.aes = list(size=4)))+
  theme(legend.key.width = unit(2.5,"cm")) +
  labs(x = "", y = "relative abundance\n") +
  theme(axis.title.y = element_text(size = 12)) +
  facet_wrap(~ OTU, scales = "free") +
  ylim(0,0.5)
ggsave(
  "order_relab_Schorf_soil_top12.png",
  width = 13,
  height = 9.5,
  dpi = 800
)

# determine top 4
order_top4_Ss <- get_top_taxa(ps_rel_abund_ord_Schorf_soil, 4, relative = TRUE, other_label = "Others", discard_other = T)

# boxplot relative abundance dominant tree Schorfheide soil

phyloseq::psmelt(order_top4_Ss) %>%
  ggplot(data = ., aes(x = dominant_tree, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  scale_color_discrete(name = "Top 12 Orders - Schorfheide soil",
                       labels = c("Agaricales","Helotiales","Russulales",
                          "Tremellales","not assigned")) + 
  theme(legend.position = "bottom",legend.direction = "horizontal") +
  theme(axis.text.x= element_text(face = "italic", size = 12)) +
  scale_x_discrete(labels = c("Fagus sylvatica","Pinus sylvestris"))+
  theme( strip.text.x = element_blank() )+
  theme(legend.text = element_text(size = 12), legend.title = element_blank())+
  guides(color = guide_legend(override.aes = list(size=4)))+
  theme(legend.key.width = unit(1.2,"cm")) +
  labs(x = "", y = "relative abundance\n") +
  theme(axis.title.y = element_text(size = 12)) +
  facet_wrap(~ OTU, scales = "free") +
  ylim(0,0.5)
ggsave(
  "order_relab_Schorf_soil_top4.png",
  width = 9.5,
  height = 7,
  dpi = 800
)

#average relative abundance top 4 orders Schorfheide soil

# subset of Agaricales in soil Schorfheide
ps_rel_abund_Schorf_soil_order_agari <- subset_taxa(physeqSchorf_s, Order == "o__Agaricales")
# percentage of Agaricales in soil Schorfheide
sum(sample_sums(ps_rel_abund_Schorf_soil_order_agari))/sum(sample_sums(physeqSchorf_s))

# subset of Helotiales in soil Schorfheide
ps_rel_abund_Schorf_soil_order_heloti <- subset_taxa(physeqSchorf_s, Order == "o__Helotiales")
# percentage of Helotiales in soil Schorfheide
sum(sample_sums(ps_rel_abund_Schorf_soil_order_heloti))/sum(sample_sums(physeqSchorf_s))

# subset of Tremellales in soil Schorfheide
ps_rel_abund_Schorf_soil_order_tremella <- subset_taxa(physeqSchorf_s, Order == "o__Tremellales")
# percentage of Tremellales in soil Schorfheide
sum(sample_sums(ps_rel_abund_Schorf_soil_order_tremella))/sum(sample_sums(physeqSchorf_s))

# subset of Russulales in soil Schorfheide
ps_rel_abund_Schorf_soil_order_russula <- subset_taxa(physeqSchorf_s, Order == "o__Russulales")
# percentage of Russulales in soil Schorfheide
sum(sample_sums(ps_rel_abund_Schorf_soil_order_russula))/sum(sample_sums(physeqSchorf_s))

# subset of Mortierellales in soil Schorfheide (to compare to Swabian Alb)
ps_rel_abund_Schorf_soil_order_mortiere <- subset_taxa(physeqSchorf_s, Order == "o__Mortierellales")
# percentage of Mortierellales in soil Schorfheide
sum(sample_sums(ps_rel_abund_Schorf_soil_order_mortiere))/sum(sample_sums(physeqSchorf_s))


###heatmap Schorfheide top 50 taxa#### (image saved: heatmap_S_top50_order)
# determine top X
#top_50_taxa_S <- get_top_taxa(ps_rel_abund_S, 49, relative = TRUE, other_label = "Others") -> alternative approach
top_50_taxa_S <- prune_taxa(names(sort(taxa_sums(ps_rel_abund_S),TRUE)[1:50]), ps_rel_abund_S)
plot_heatmap(top_50_taxa_S, taxa.label="Order")

####boxplots of abundance on order level####

# aggregate taxa on order level
physeq_order_S <- aggregate_taxa(physeqSchorf, level = 'Order')
# calculate relative abundance
ps_rel_abund_ord_S = phyloseq::transform_sample_counts(physeq_order_S, function(x){x / sum(x)})
# determine top X
order_top12_S <- get_top_taxa(ps_rel_abund_ord_S, 11, relative = TRUE, other_label = "Others")

# boxplot relative abundance dominant tree
phyloseq::psmelt(order_top12_S) %>%
  ggplot(data = ., aes(x = dominant_tree, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")

################Alb#######################

##order overview Alb##
#count orders
table(phyloseq::tax_table(physeqAlb)[,"Order"])
write.csv(table(phyloseq::tax_table(physeqAlb)[,"Order"]), "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/overview_orders_Alb.csv" )

#####prepare data for barplots subset Alb####
ps_rel_abund_A = phyloseq::transform_sample_counts(physeqAlb, function(x){x / sum(x)})

phyloseq::otu_table(physeqAlb)[1:5, 1:5] #in this case (choice of display range) 0 counts
phyloseq::otu_table(ps_rel_abund_A)[1:5, 1:5] #logically same as above
####barplots relative abundance Schorfheide####
#barplot Phyla by substrate##
phyloseq::plot_bar(ps_rel_abund_A, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ substrate, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
#barplot Phyla by Fagus vs. Picea (by dominant tree)##
phyloseq::plot_bar(ps_rel_abund_A, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ dominant_tree, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
#TEST Order by Fagus vs. Picea (by dominant tree)##
#result: extension of ledgend makes plot unreadable
phyloseq::plot_bar(ps_rel_abund_A, fill = "Order") +
  geom_bar(aes(color = Order, fill = Order), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ dominant_tree, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

######subsets soil & bark Alb###

##soil
#prepare data for barplots
ps_rel_abund_As = phyloseq::transform_sample_counts(physeqAlb_s, function(x){x / sum(x)})
phyloseq::otu_table(physeqAlb_s)[1:5, 1:5] #in this case (choice of display range) 0 counts
phyloseq::otu_table(ps_rel_abund_As)[1:5, 1:5] #logically same as above

#by dominant_tree subset soil
phyloseq::plot_bar(ps_rel_abund_As, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ dominant_tree, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

##bark
#prepare data for barplots
ps_rel_abund_Ab = phyloseq::transform_sample_counts(physeqAlb_b, function(x){x / sum(x)})
phyloseq::otu_table(physeqAlb_b)[1:5, 1:5] #in this case (choice of display range) 0 counts
phyloseq::otu_table(ps_rel_abund_Ab)[1:5, 1:5] #logically same as above


#by dominant_tree subset bark
phyloseq::plot_bar(ps_rel_abund_Ab, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ dominant_tree, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

###boxplots of abundance on order level####

# aggregate taxa on order level
physeq_order_Alb <- aggregate_taxa(physeqAlb, level = 'Order')
# calculate relative abundance
ps_rel_abund_ord_Alb = phyloseq::transform_sample_counts(physeq_order_Alb, function(x){x / sum(x)})
# determine top X
order_top12_A <- get_top_taxa(ps_rel_abund_ord_Alb, 11, relative = TRUE, other_label = "Others")

# boxplot relative abundance dominant tree
phyloseq::psmelt(order_top12_A) %>%
  ggplot(data = ., aes(x = dominant_tree, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")

####
##boxplot abundance order level split by substrate Alb
#####

#Alb bark by dominant tree

# aggregate taxa on order level
physeq_order_Alb_bark <- aggregate_taxa(physeqAlb_b, level = 'Order')
# calculate relative abundance
ps_rel_abund_ord_Alb_bark = phyloseq::transform_sample_counts(physeq_order_Alb_bark, function(x){x / sum(x)})
# determine top 12
order_top12_Ab <- get_top_taxa(ps_rel_abund_ord_Alb_bark, 12, relative = TRUE, other_label = "Others", discard_other = T)

# boxplot relative abundance dominant tree Alb bark

phyloseq::psmelt(order_top12_Ab) %>%
  ggplot(data = ., aes(x = dominant_tree, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  scale_color_discrete(name = "Top 12 Orders - Swabian Alb bark",
                       labels = c("Calicales","Capnodiales","Chaetothyriales","Helotiales",
                                "Lecanorales","Myriangiales","Orbiliales","Pleosporales",
                                "Trapeliales","Tremellales","Verrucariales","not assigned")) +
  theme(legend.position = "bottom",legend.direction = "horizontal") +
  theme(axis.text.x= element_text(face = "italic", size = 12)) +
  scale_x_discrete(labels = c("Fagus sylvatica","Picea abies"))+
  theme( strip.text.x = element_blank() )+
  theme(legend.text = element_text(size = 12), legend.title = element_blank())+
  guides(color = guide_legend(override.aes = list(size=4)))+
  theme(legend.key.width = unit(2,"cm")) +
  labs(x = "", y = "relative abundance\n") +
  theme(axis.title.y = element_text(size = 12)) +
  facet_wrap(~ OTU, scales = "free") +
  ylim(0,0.6)
ggsave(
  "order_relab_Alb_bark_top12.png",
  width = 11,
  height = 9.5,
  dpi = 800
)

# determine top 4
order_top4_Ab <- get_top_taxa(ps_rel_abund_ord_Alb_bark, 4, relative = TRUE, other_label = "Others", discard_other = T)

# boxplot relative abundance dominant tree Alb bark

phyloseq::psmelt(order_top4_Ab) %>%
  ggplot(data = ., aes(x = dominant_tree, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  scale_color_discrete(name = "Top 4 Orders - Swabian Alb bark",
                       labels = c("Capnodiales","Lecanorales",
                                 "Trapeliales","not assigned")) +
  theme(legend.position = "bottom",legend.direction = "horizontal") +
  theme(axis.text.x= element_text(face = "italic", size = 12)) +
  scale_x_discrete(labels = c("Fagus sylvatica","Picea abies"))+
  theme( strip.text.x = element_blank() )+
  theme(legend.text = element_text(size = 12), legend.title = element_blank())+
  guides(color = guide_legend(override.aes = list(size=4)))+
  theme(legend.key.width = unit(1.2,"cm")) +
  labs(x = "", y = "relative abundance\n") +
  theme(axis.title.y = element_text(size = 12)) +
  facet_wrap(~ OTU, scales = "free") +
  ylim(0,0.6)
ggsave(
  "order_relab_Alb_bark_top4.png",
  width = 9.5,
  height = 7,
  dpi = 800
)

#average relative abundance top 4 orders Alb bark

# subset of unassigned orders on bark Alb
ps_rel_abund_Alb_bark_order_na <- subset_taxa(physeqAlb_b, is.na(Order))
# percentage of unassigned orders on bark Alb
sum(sample_sums(ps_rel_abund_Alb_bark_order_na))/sum(sample_sums(physeqAlb_b))

# subset of Capnodiales on bark Alb
ps_rel_abund_Alb_bark_order_capno <- subset_taxa(physeqAlb_b, Order == "o__Capnodiales")
# percentage of Capnodiales on bark Alb
sum(sample_sums(ps_rel_abund_Alb_bark_order_capno))/sum(sample_sums(physeqAlb_b))

# subset of Lecanorales on bark Alb
ps_rel_abund_Alb_bark_order_lecano <- subset_taxa(physeqAlb_b, Order == "o__Lecanorales")
# percentage of Lecanorales on bark Alb
sum(sample_sums(ps_rel_abund_Alb_bark_order_lecano))/sum(sample_sums(physeqAlb_b))

# subset of Trapeliales on bark Alb
ps_rel_abund_Alb_bark_order_trapeli <- subset_taxa(physeqAlb_b, Order == "o__Trapeliales")
# percentage of Trapeliales on bark Alb
sum(sample_sums(ps_rel_abund_Alb_bark_order_trapeli))/sum(sample_sums(physeqAlb_b))

# subset of Chaetothyrials on bark Alb (to compare to Schorfheide)
ps_rel_abund_Alb_bark_order_chaeto <- subset_taxa(physeqAlb_b, Order == "o__Chaetothyriales")
# percentage of Chaetothyriales on bark Alb
sum(sample_sums(ps_rel_abund_Alb_bark_order_chaeto))/sum(sample_sums(physeqAlb_b))

#Alb soil by dominant tree

# aggregate taxa on order level
physeq_order_Alb_soil <- aggregate_taxa(physeqAlb_s, level = 'Order')
# calculate relative abundance
ps_rel_abund_ord_Alb_soil = phyloseq::transform_sample_counts(physeq_order_Alb_soil, function(x){x / sum(x)})
# determine top 12
order_top12_As <- get_top_taxa(ps_rel_abund_ord_Alb_soil, 12, relative = TRUE, other_label = "Others", discard_other = T)

# boxplot relative abundance dominant tree Alb soil
phyloseq::psmelt(order_top12_As) %>%
  ggplot(data = ., aes(x = dominant_tree, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  scale_color_discrete(name = "Top 12 Orders - Swabian Alb bark",
                      labels = c("Agaricales","Atheliales","Cantharellales","Helotiales",
                                "Hypocreales","Mortierellales","Russulales","Sebacinales",
                               "Thelebolales","Thelephorales","Tremellales","not assigned")) + 
  theme(legend.position = "bottom",legend.direction = "horizontal") +
  theme(axis.text.x= element_text(face = "italic", size = 12)) +
  scale_x_discrete(labels = c("Fagus sylvatica","Picea abies"))+
  theme( strip.text.x = element_blank() )+
  theme(legend.text = element_text(size = 12), legend.title = element_blank())+
  guides(color = guide_legend(override.aes = list(size=4)))+
  theme(legend.key.width = unit(2,"cm")) +
  labs(x = "", y = "relative abundance\n") +
  theme(axis.title.y = element_text(size = 12)) +
  facet_wrap(~ OTU, scales = "free") +
  ylim(0,0.6)
ggsave(
  "order_relab_Alb_soil_top12.png",
  width = 11,
  height = 9.5,
  dpi = 800
)


# determine top 4
order_top4_As <- get_top_taxa(ps_rel_abund_ord_Alb_soil, 4, relative = TRUE, other_label = "Others", discard_other = T)

# boxplot relative abundance dominant tree Alb soil
phyloseq::psmelt(order_top4_As) %>%
  ggplot(data = ., aes(x = dominant_tree, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  scale_color_discrete(name = "Top 12 Orders - Swabian Alb bark",
                      labels = c("Agaricales","Helotiales",
                                "Mortierellales","Russulales")) + 
  theme(legend.position = "bottom",legend.direction = "horizontal") +
  theme(axis.text.x= element_text(face = "italic", size = 12)) +
  scale_x_discrete(labels = c("Fagus sylvatica","Picea abies"))+
  theme( strip.text.x = element_blank() )+
  theme(legend.text = element_text(size = 12), legend.title = element_blank())+
  guides(color = guide_legend(override.aes = list(size=4)))+
  theme(legend.key.width = unit(1.2,"cm")) +
  labs(x = "", y = "relative abundance\n") +
  theme(axis.title.y = element_text(size = 12)) +
  facet_wrap(~ OTU, scales = "free") +
  ylim(0,0.6)
ggsave(
  "order_relab_Alb_soil_top4.png",
  width = 9.5,
  height = 7,
  dpi = 800
)

#average relative abundance top 4 orders Alb soil

# subset of Agaricales in soil Alb
ps_rel_abund_Alb_soil_order_agari <- subset_taxa(physeqAlb_s, Order == "o__Agaricales")
# percentage of Agaricales in soil Alb
sum(sample_sums(ps_rel_abund_Alb_soil_order_agari))/sum(sample_sums(physeqAlb_s))

# subset of Helotiales in soil Alb
ps_rel_abund_Alb_soil_order_heloti <- subset_taxa(physeqAlb_s, Order == "o__Helotiales")
# percentage of Helotiales in soil Alb
sum(sample_sums(ps_rel_abund_Alb_soil_order_heloti))/sum(sample_sums(physeqAlb_s))

# subset of Mortierellales in soil Alb
ps_rel_abund_Alb_soil_order_mortiere <- subset_taxa(physeqAlb_s, Order == "o__Mortierellales")
# percentage of Mortierellales in soil Alb
sum(sample_sums(ps_rel_abund_Alb_soil_order_mortiere))/sum(sample_sums(physeqAlb_s))

# subset of Russulales in soil Alb
ps_rel_abund_Alb_soil_order_russula <- subset_taxa(physeqAlb_s, Order == "o__Russulales")
# percentage of Russulales in soil Alb
sum(sample_sums(ps_rel_abund_Alb_soil_order_russula))/sum(sample_sums(physeqAlb_s))

# subset of Tremellales in soil Alb (to compare to Schorfheide)
ps_rel_abund_Alb_soil_order_tremella <- subset_taxa(physeqAlb_s, Order == "o__Tremellales")
# percentage of Tremellales in soil Alb
sum(sample_sums(ps_rel_abund_Alb_soil_order_tremella))/sum(sample_sums(physeqAlb_s))



###heatmap Alb top 50 taxa####
# determine top X
#top_50_taxa_A <- get_top_taxa(ps_rel_abund_A, 49, relative = TRUE, other_label = "Others")
table(phyloseq::tax_table(ps_rel_abund_A)[,"Order"])
top_50_taxa_A <- prune_taxa(names(sort(taxa_sums(ps_rel_abund_A),TRUE)[1:50]), ps_rel_abund_A)
plot_heatmap(top_50_taxa_A, taxa.label="Order")

###
##### Phylum relative abundance####
#####faceted stacked bar plot with 4 subgroups###
######(grouped by venn class####
###

my_cols <- paletteer_d('ggsci::default_igv')

##Alb

#phyloseq::plot_bar(ps_rel_abund_A, fill = "Phylum") +
# geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
#scale_fill_manual(values = ggplot2::alpha(my_cols, 0.9), name = 'Phylum') +
#labs(x = "", y = "Relative Abundance\n") +
#facet_wrap(~ venn_class, scales = "free") +
#theme(panel.background = element_blank(),
#     axis.text.x=element_blank(),
#    axis.ticks.x=element_blank())

#phyloseq::plot_bar(ps_rel_abund_A, fill = "Phylum") +
 # geom_bar(aes(fill = Phylum), stat = "identity", position = "stack") +
  #scale_fill_manual(values = ggplot2::alpha(my_cols, 0.9), name = 'Phylum') +
  #labs(x = "", y = "Relative Abundance\n") +
  #facet_wrap(~ venn_class, scales = "free") +
  #theme(panel.background = element_blank(),
   #     axis.text.x=element_blank(),
    #    axis.ticks.x=element_blank())



### test 20220324####
# agglomerate taxa
glom_Alb <- tax_glom(ps_rel_abund_A, taxrank = 'Phylum')
# create dataframe from phyloseq object
dat_A <- psmelt(glom_Alb)
# convert Phylum to a character vector from a factor because R
dat_A$Phylum <- as.character(dat_A$Phylum)
# group dataframe by Phylum, calculate median rel. abundance
library(plyr)
medians_A <- ddply(dat_A, ~Phylum, function(x) c(median=median(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
other_A <- medians_A[medians_A$median <= 0.00005,]$Phylum

# change their name to "Remainder"
dat_A[dat_A$Phylum %in% other_A,]$Phylum <- 'Other'

# boxplot
ggplot(dat_A,
       aes(x=Phylum,
           y=Abundance)) + geom_boxplot() + coord_flip()

ggplot(dat_A,
       aes(x=Phylum,
           y=Abundance)) +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "stack")+
scale_fill_manual(values = pals::glasbey(n=18)) +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ venn_class, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
#wo sind die unknown, sample information fehlt
##test ende

#####phylum stacked barplots##

##20220609 stacked barplot with rel. abundance > 1 %

phy_abundfilt_A1 = filter_taxa(ps_rel_abund_A, function(x) sum(x) > .01, TRUE)

phyloseq::plot_bar(phy_abundfilt_A1, fill = "Phylum") +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "stack") +
  scale_fill_manual(values = pals::brewer.set1(n=9), 
                    labels =c("Ascomycota", "Basidiomycota", "Chytridiomycota",
                              "Glomeromycota", "Mortierellomycota", "Mucoromycota", 
                              "Olpidiomycota", "Rozellomycota", "not assigned")) +
  labs(x = "", y = "relative abundance\n(> 1%)") +
  facet_wrap(~ venn_class, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold"))+
  theme(legend.position = "bottom",legend.direction = "horizontal") +
  theme(strip.text.x = element_text(size = 14, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold")) +
  theme(axis.text = element_text(size = 14))

##20220726 stacked barplot with rel. abundance > 1 % + manual colors
phyloseq::plot_bar(phy_abundfilt_A1, fill = "Phylum") +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "stack") +
  scale_fill_manual(values = c("firebrick1", "deepskyblue3", "green3", "darkorchid3",
                               "darkorange1", "turquoise3", "black", "hotpink"), 
                    labels =c("Ascomycota", "Basidiomycota", "Chytridiomycota",
                              "Glomeromycota", "Mortierellomycota", "Mucoromycota", 
                              "Olpidiomycota", "Rozellomycota", "not assigned")) +
  labs(x = "", y = "relative abundance\n(> 1%)") +
  facet_wrap(~ venn_class, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20, face = "bold"))+
  theme(legend.position = "bottom",legend.direction = "horizontal") +
  theme(strip.text.x = element_text(size = 20, face = "bold")) +
  theme(axis.title = element_text(size = 20, face = "bold")) +
  theme(axis.text = element_text(size = 20))
ggsave(
  "phylum_relab_Alb_paper.png",
  width = 15,
  height = 10.64,
  dpi = 1200
)


#dataset without abundance filter

phyloseq::plot_bar(ps_rel_abund_A, fill = "Phylum") +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "stack") +
  scale_fill_manual(values = pals::glasbey(n=18), 
                    labels =c("Aphelidiomycota", "Ascomycota", "Basidiobolomycota", "Basidiomycota", "Calcarisporiellomycota",
                              "Chytridiomycota", "Entorrhizomycota", "Fungi phy Incertae sedis", "Glomeromycota", "Kickxellomycota",
                              "Monoblepharomycota", "Mortierellomycota", "Mucoromycota", "Neocallimastigomycota", "Olpidiomycota", "Rozellomycota",
                              "Zoopagomycota", "not assigned")) +
  labs(x = "", y = "relative abundance\n") +
  facet_wrap(~ venn_class, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"))+
  theme(legend.position = "bottom",legend.direction = "horizontal") +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 12, face = "bold")) +
  theme(axis.text = element_text(size = 12))
ggsave(
  "relabund_Alb_Phylum.png",
  width = 11.7,
  height = 8.3,
  dpi = 1200
)

# subset of unassigned Phyla on bark Alb
ps_rel_abund_Alb_bark_phylum_na <- subset_taxa(physeqAlb_b, is.na(Phylum) )

# percentage of unassigned Phyla on bark per sample Alb
sample_sums(ps_rel_abund_Alb_bark_phylum_na)/sample_sums(physeqAlb_b)

# percentage of unassigned Phyla on bark Alb
sum(sample_sums(ps_rel_abund_Alb_bark_phylum_na))/sum(sample_sums(physeqAlb_b))

# subset of unassigned Phyla on soil Alb
ps_rel_abund_Alb_soil_phylum_na <- subset_taxa(physeqAlb_s, is.na(Phylum) )

# percentage of unassigned Phyla on soil per sample Alb
sample_sums(ps_rel_abund_Alb_soil_phylum_na)/sample_sums(physeqAlb_s)

# percentage of unassigned Phyla on soil Alb
sum(sample_sums(ps_rel_abund_Alb_soil_phylum_na))/sum(sample_sums(physeqAlb_s))

# subset of Ascomycota on bark Alb
ps_rel_abund_Alb_bark_phylum_asco <- subset_taxa(physeqAlb_b, Phylum == "p__Ascomycota")

# percentage of Ascomycota on bark per sample Alb
sample_sums(ps_rel_abund_Alb_bark_phylum_asco)/sample_sums(physeqAlb_b)

# percentage of Ascomycota on bark Alb
sum(sample_sums(ps_rel_abund_Alb_bark_phylum_asco))/sum(sample_sums(physeqAlb_b))

# subset of Basidiomycota on soil Alb
ps_rel_abund_Alb_soil_phylum_basidio <- subset_taxa(physeqAlb_s, Phylum == "p__Basidiomycota")

# percentage of Basidiomycota on soil per sample Alb
sample_sums(ps_rel_abund_Alb_soil_phylum_basidio)/sample_sums(physeqAlb_s)

# percentage of Basidiomycota on soil Alb
sum(sample_sums(ps_rel_abund_Alb_soil_phylum_basidio))/sum(sample_sums(physeqAlb_s))


##Schorfheide
#phyloseq::plot_bar(ps_rel_abund_S, fill = "Phylum") +
# geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
#labs(x = "", y = "Relative Abundance\n") +
#facet_wrap(~ venn_class, scales = "free") +
#theme(panel.background = element_blank(),
#     axis.text.x=element_blank(),
#    axis.ticks.x=element_blank())

#phyloseq::plot_bar(ps_rel_abund_S, fill = "Phylum") +
 # geom_bar(aes(fill = Phylum), stat = "identity", position = "stack") +
  #scale_fill_manual(values = ggplot2::alpha(my_cols, 0.9), name = 'Phylum') +
  #labs(x = "", y = "Relative Abundance\n") +
  #facet_wrap(~ venn_class, scales = "free") +
  #theme(panel.background = element_blank(),
   #     axis.text.x=element_blank(),
    #    axis.ticks.x=element_blank())

##20220609 stacked barplot with rel. abundance > 1 %

phy_abundfilt_S1 = filter_taxa(ps_rel_abund_S, function(x) sum(x) > .01, TRUE)

phyloseq::plot_bar(phy_abundfilt_S1, fill = "Phylum") +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "stack") +
  scale_fill_manual(values = pals::brewer.set1(n=6), 
                    labels =c("Ascomycota", "Basidiomycota", "Mortierellomycota", "Mucoromycota",
                              "Rozellomycota", "not assigned")) +
  labs(x = "", y = "relative abundance\n(> 1%)") +
  facet_wrap(~ venn_class, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold"))+
  theme(legend.position = "bottom",legend.direction = "horizontal") +
  theme(strip.text.x = element_text(size = 14, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold")) +
  theme(axis.text = element_text(size = 14))

##20220726 stacked barplot with rel. abundance > 1 % + manual colors
phyloseq::plot_bar(phy_abundfilt_S1, fill = "Phylum") +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "stack") +
  scale_fill_manual(values = c("firebrick1", "deepskyblue3", "darkorange1", "turquoise3",
                               "hotpink"), 
                    labels =c("Ascomycota", "Basidiomycota", "Mortierellomycota", "Mucoromycota",
                              "Rozellomycota", "not assigned")) +
  labs(x = "", y = "relative abundance\n(> 1%)") +
  facet_wrap(~ venn_class, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20, face = "bold"))+
  theme(legend.position = "bottom",legend.direction = "horizontal") +
  theme(strip.text.x = element_text(size = 20, face = "bold")) +
  theme(axis.title = element_text(size = 20, face = "bold")) +
  theme(axis.text = element_text(size = 20))
ggsave(
  "phylum_relab_Schorf_paper.png",
  width = 11.7,
  height = 8.3,
  dpi = 1200
)

#dataset without abundance filter

phyloseq::plot_bar(ps_rel_abund_S, fill = "Phylum") +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "stack") +
  scale_fill_manual(values = pals::glasbey(n=18), 
                    labels =c("Aphelidiomycota", "Ascomycota", "Basidiobolomycota", "Basidiomycota", "Calcarisporiellomycota",
                              "Chytridiomycota", "Entorrhizomycota", "Fungi phy Incertae sedis", "Glomeromycota", "Kickxellomycota",
                              "Monoblepharomycota", "Mortierellomycota", "Mucoromycota", "Neocallimastigomycota", "Olpidiomycota", "Rozellomycota",
                              "Zoopagomycota", "not assigned")) +
  labs(x = "", y = "relative abundance\n") +
  facet_wrap(~ venn_class, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"))+
  theme(legend.position = "bottom",legend.direction = "horizontal") +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 12, face = "bold")) +
  theme(axis.text = element_text(size = 12))
ggsave(
  "relabund_Schorf_Phylum.png",
  width = 11.7,
  height = 8.3,
  dpi = 1200
)

# subset of unassigned Phyla on bark Schorfheide
ps_rel_abund_Schorf_bark_phylum_na <- subset_taxa(physeqSchorf_b, is.na(Phylum) )

# percentage of unassigned Phyla on bark per sample Schorfheide
sample_sums(ps_rel_abund_Schorf_bark_phylum_na)/sample_sums(physeqSchorf_b)

# percentage of unassigned Phyla on bark Schorfheide
sum(sample_sums(ps_rel_abund_Schorf_bark_phylum_na))/sum(sample_sums(physeqSchorf_b))

# subset of unassigned Phyla on soil Schorfheide
ps_rel_abund_Schorf_soil_phylum_na <- subset_taxa(physeqSchorf_s, is.na(Phylum) )

# percentage of unassigned Phyla on soil per sample Schorfheide
sample_sums(ps_rel_abund_Schorf_soil_phylum_na)/sample_sums(physeqSchorf_s)

# percentage of unassigned Phyla on soil Schorf
sum(sample_sums(ps_rel_abund_Schorf_soil_phylum_na))/sum(sample_sums(physeqSchorf_s))

# subset of Ascomycota on bark Schorf
ps_rel_abund_Schorf_bark_phylum_asco <- subset_taxa(physeqSchorf_b, Phylum == "p__Ascomycota")

# percentage of Ascomycota on bark per sample Schorfheide
sample_sums(ps_rel_abund_Schorf_bark_phylum_asco)/sample_sums(physeqSchorf_b)

# percentage of Ascomycota on bark Schorfheide
sum(sample_sums(ps_rel_abund_Schorf_bark_phylum_asco))/sum(sample_sums(physeqSchorf_b))

# subset of Basidiomycota on soil Schorfheide
ps_rel_abund_Schorf_soil_phylum_basidio <- subset_taxa(physeqSchorf_s, Phylum == "p__Basidiomycota")

# percentage of Basidiomycota on soil per sample Schorfheide
sample_sums(ps_rel_abund_Schorf_soil_phylum_basidio)/sample_sums(physeqSchorf_s)

# percentage of Basidiomycota on soil Schorfheide
sum(sample_sums(ps_rel_abund_Schorf_soil_phylum_basidio))/sum(sample_sums(physeqSchorf_s))




################all exploratories#######################

#barplots relative abundance
  #by substrate
phyloseq::plot_bar(ps_rel_abund, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ substrate, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  #by substrate with different color palette
phyloseq::plot_bar(ps_rel_abund, fill = "Phylum") +
  scale_fill_manual(values = ggplot2::alpha(my_cols, 0.9), name = "Phylum") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ substrate, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  #by exploratory
phyloseq::plot_bar(ps_rel_abund, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ exploratory, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  #by dominant_tree
phyloseq::plot_bar(ps_rel_abund, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ dominant_tree, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  #by tree_type
phyloseq::plot_bar(ps_rel_abund, fill = "Phylum") +
  scale_fill_manual(values = ggplot2::alpha(my_cols, 0.9), name = "Phylum") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ tree_type, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


#create boxplots of abundance by phylum!
#Agglomerate to phylum-level and rename
ps_phylum <- phyloseq::tax_glom(physeq4, "Phylum")
phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::otu_table(ps_phylum)[1:5, 1:5]
  # by substrate
phyloseq::psmelt(ps_phylum) %>%
  ggplot(data = ., aes(x = substrate, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")
  #by tree type
phyloseq::psmelt(ps_phylum) %>%
  ggplot(data = ., aes(x = tree_type, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")
  #by tree type relative abundance
phyloseq::psmelt(ps_rel_abund) %>%
  ggplot(data = ., aes(x = tree_type, y = "Relative Abundance\n")) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")
####boxplots of abundance by orders#####
ps_order <- phyloseq::tax_glom(physeq4, "Order")
phyloseq::taxa_names(ps_order) <- phyloseq::tax_table(ps_order)[, "Order"]
phyloseq::otu_table(ps_order)[1:5, 1:5]
phyloseq::psmelt(ps_order) %>%
  ggplot(data = ., aes(x = substrate, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")


#versuch mit anderer Variante:
physeq_order <- aggregate_taxa(physeq5, level = 'Order')
ps_rel_abund_ord = phyloseq::transform_sample_counts(physeq_order, function(x){x / sum(x)})
order_top25 <- get_top_taxa(ps_rel_abund_ord, 24, relative = TRUE, other_label = "Others")

ps_rel_abund_20 = phyloseq::transform_sample_counts(ps_order20, function(x){x / sum(x)})
phyloseq::otu_table(ps_order20)[1:5, 1:5] #in this case (choice of display range) 0 counts
phyloseq::otu_table(ps_rel_abund_20)[1:5, 1:5] #logically same as above
#barplot top 20 orders by substrate##
phyloseq::plot_bar(order_top25, fill = "Order") +
  geom_bar(aes(color = Order, fill = Order), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ substrate, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())



###hypothesis testing by HMP -> difference of phylum composition
# no sensible result from the test so far (p = 0)!!!!
#Subset groups
bark <- phyloseq::subset_samples(ps_phylum, substrate == "bark")
soil <- phyloseq::subset_samples(ps_phylum, substrate == "soil")
#Output OTU tables
bark_otu <- data.frame(phyloseq::otu_table(bark))
soil_otu <- data.frame(phyloseq::otu_table(soil))
#Group rare phyla (threshold 20 chosen)
bark_otu <- bark_otu %>%
  t(.) %>%
  as.data.frame(.) %>%
  mutate(Other = 	p__Aphelidiomycota + p__Basidiobolomycota + p__Calcarisporiellomycota + p__Entorrhizomycota + p__Fungi_phy_Incertae_sedis +
           p__Monoblepharomycota + p__Neocallimastigomycota + p__Olpidiomycota) %>%
  dplyr::select(-p__Aphelidiomycota, - p__Basidiobolomycota, - p__Calcarisporiellomycota, - p__Entorrhizomycota, - p__Fungi_phy_Incertae_sedis,
                - p__Monoblepharomycota, - p__Neocallimastigomycota, - p__Olpidiomycota)
soil_otu <- soil_otu %>%
  t(.) %>%
  as.data.frame(.) %>%
  mutate(Other = p__Aphelidiomycota + p__Basidiobolomycota + p__Calcarisporiellomycota + p__Entorrhizomycota + p__Fungi_phy_Incertae_sedis +
           p__Monoblepharomycota + p__Neocallimastigomycota + p__Olpidiomycota) %>%
  dplyr::select(-p__Aphelidiomycota, - p__Basidiobolomycota, - p__Calcarisporiellomycota, - p__Entorrhizomycota, - p__Fungi_phy_Incertae_sedis,
                - p__Monoblepharomycota, - p__Neocallimastigomycota, - p__Olpidiomycota)
#HMP test
group_data <- list(bark_otu, soil_otu)
(xdc <- HMP::Xdc.sevsample(group_data))

###hierarchical clustering by bray-curtis dissimilarity
#Extract OTU table and compute BC
ps_rel_otu <- data.frame(phyloseq::otu_table(ps_rel_abund))
ps_rel_otu <- t(ps_rel_otu)
bc_dist <- vegan::vegdist(ps_rel_otu, method = "bray")
as.matrix(bc_dist)[1:5, 1:5]
#Save as dendrogram
ward <- as.dendrogram(hclust(bc_dist, method = "ward.D2"))
#Provide color codes
meta <- data.frame(phyloseq::sample_data(ps_rel_abund))
colorCode <- c("bark" = "burlywood4", "soil" = "chocolate1")
labels_colors(ward) <- colorCode[meta$substrate][order.dendrogram(ward)]
#Plot
plot(ward)


###observed richness vs total reads
ggplot(data = data.frame("total_reads" =  phyloseq::sample_sums(physeq4),
                         "observed" = phyloseq::estimate_richness(physeq4, measures = "Observed")[, 1]),
       aes(x = total_reads, y = observed)) +
  geom_point() +
  geom_smooth(method="lm", se = FALSE) +
  labs(x = "\nTotal Reads", y = "Observed Richness\n")

###beta diversity####

##ordination with PCA -> CLR transformation -> values are not counts but dominance
(ps_clr <- microbiome::transform(physeq4, "clr"))  
phyloseq::otu_table(ps_clr)[1:5,1:5]
#PCA via phyloseq
ord_clr <- phyloseq::ordinate(ps_clr, "RDA")
#Plot scree plot
phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
#Examine eigenvalues and % prop. variance explained
head(ord_clr$CA$eig)
sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig))
#Scale axes and plot ordination
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)

#by substrate
phyloseq::plot_ordination(physeq4, ord_clr, type="samples", color="substrate") + 
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = substrate), linetype = 2)
#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "euclidean") 
#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$substrate)
#Dispersion test and plot
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(ps_clr)$substrate)
dispr
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
#boxplot distance to centroid
boxplot(dispr, main ="", xlab = "")
#permutation test
permutest(dispr)

#by tree_type
phyloseq::plot_ordination(physeq4, ord_clr, type="samples", color="tree_type") + 
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = substrate), linetype = 2)
#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "euclidean") 
#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$substrate)
#Dispersion test and plot
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(ps_clr)$tree_type)
dispr
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")

##phylogenetic information in beta diversity - PCoA analysis
#generate distances
ord_unifrac <- ordinate(ps_rare, method = "PCoA", distance = "wunifrac") 
ord_unifrac_un <- ordinate(ps_rare, method = "PCoA", distance = "unifrac")
#Plot ordinations
a <- plot_ordination(ps_rare, ord_unifrac, color = "substrate") + geom_point(size = 2)
b <- plot_ordination(ps_rare, ord_unifrac_un, color = "substrate") + geom_point(size = 2)
cowplot::plot_grid(a, b, nrow = 1, ncol = 2, scale = .9, labels = c("Weighted", "Unweighted"))

#######differential abundance testing#######
#Generate data.frame with OTUs and metadata
ps_wilcox <- data.frame(t(data.frame(phyloseq::otu_table(ps_clr))))
ps_wilcox$substrate <- phyloseq::sample_data(ps_clr)$substrate
#Define functions to pass to map
wilcox_model <- function(df){
  wilcox.test(abund ~ substrate, data = df)
}
wilcox_pval <- function(df){
  wilcox.test(abund ~ substrate, data = df)$p.value
}
#Create nested data frames by ASV and loop over each using map 
wilcox_results <- ps_wilcox %>%
  gather(key = ASV, value = abund, -substrate) %>%
  group_by(ASV) %>%
  nest() %>%
  mutate(wilcox_test = map(data, wilcox_model),
         p_value = map(data, wilcox_pval))                       
#Show results
head(wilcox_results)
head(wilcox_results$data[[1]])
wilcox_results$wilcox_test[[1]]
wilcox_results$p_value[[1]]
#unnesting
wilcox_results <- wilcox_results %>%
  dplyr::select(ASV, p_value) %>%
  unnest()
head(wilcox_results)
#Adding taxonomic labels
taxa_info <- data.frame(tax_table(ps_clr))
taxa_info <- taxa_info %>% rownames_to_column(var = "ASV")
#Computing FDR corrected p-values
wilcox_results <- wilcox_results %>%
  full_join(taxa_info) %>%
  arrange(p_value) %>%
  mutate(BH_FDR = p.adjust(p_value, "BH")) %>%
  filter(BH_FDR < 0.05) %>%
  dplyr::select(ASV, p_value, BH_FDR, everything())
#printing results
print.data.frame(wilcox_results)
##ANOVA-like differential expression (ALDEx2) #20220118 funktioniert nicht. verschiedene Fehler abhängig von $condition
#Run ALDEx2
aldex2_da <- ALDEx2::aldex(data.frame(phyloseq::otu_table(physeq4)), phyloseq::sample_data(physeq4)$tree_type, test="t", effect = TRUE, denom="iqlr")
#plot effect sizes
ALDEx2::aldex.plot(aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)
#extraction of PCs from PCA as subset for linear model
#Generate data.frame
clr_pcs <- data.frame(
  "pc1" = ord_clr$CA$u[,1],
  "pc2" = ord_clr$CA$u[,2],
  "pc3" = ord_clr$CA$u[,3],
  "substrate" = phyloseq::sample_data(ps_clr)$substrate
)
clr_pcs$substrate_num <- ifelse(clr_pcs$substrate == "bark", 0, 1) #soil =1 
head(clr_pcs)
#Specify a datadist object (for rms)
dd <- datadist(clr_pcs)
options(datadist = "dd")
#Plot the unconditional associations
a <- ggplot(clr_pcs, aes(x = pc1, y = substrate_num)) +
  Hmisc::histSpikeg(substrate_num ~ pc1, lowess = TRUE, data = clr_pcs) +
  labs(x = "\nPC1", y = "Pr(soil)\n")
b <- ggplot(clr_pcs, aes(x = pc2, y = substrate_num)) +
  Hmisc::histSpikeg(substrate_num ~ pc2, lowess = TRUE, data = clr_pcs) +
  labs(x = "\nPC2", y = "Pr(soil)\n")
c <- ggplot(clr_pcs, aes(x = pc3, y = substrate_num)) +
  Hmisc::histSpikeg(substrate_num ~ pc3, lowess = TRUE, data = clr_pcs) +
  labs(x = "\nPC3", y = "Pr(soil)\n")
cowplot::plot_grid(a, b, c, nrow = 2, ncol = 2, scale = .9, labels = "AUTO")
#Fit full model with splines (3 knots each)
m1 <- rms::lrm(substrate_num ~ rcs(pc1, 3) + rcs(pc2, 3) + rcs(pc3, 3), data = clr_pcs, x = TRUE, y = TRUE) #20220118 funktioniert nicht
#Grid search for penalties
pentrace(m1, list(simple = c(0, 1, 2), nonlinear = c(0, 100, 200)))


######################################################
############management intensity#########################
######################################################

################all exploratories#######################

#barplot management intesity by dominant tree##
phyloseq::plot_bar(physeq4, fill = "ForMIclass") +
  geom_bar(aes(color = ForMIclass, fill = ForMIclass), stat = "identity", position = "stack") +
  #labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ dominant_tree, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#barplot management intesity by tree_type##
phyloseq::plot_bar(physeq4, fill = "ForMIclass") +
  geom_bar(aes(color = ForMIclass, fill = ForMIclass), stat = "identity", position = "stack") +
  #labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ tree_type, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

######################################################
############network analysis#########################
###################################################### 

################Alb#######################

#######network analysis Physalia Metabarcoding course###

#by taxa
ig1 <- make_network(physeqAlb, "taxa", distance = "jaccard", max.dist = 0.95)
plot_network(ig1, physeqAlb, type="taxa", point_size = 5, label=NULL, color="Phylum", line_alpha = 0.05)

#network samples - substrate
ig2 <- make_network(physeqAlb, "samples", distance = "jaccard", max.dist = 0.95)
plot_network(ig2, physeqAlb, type="samples", point_size = 5, label=NULL, color="substrate", line_alpha = 0.05)

#network samples - dominant_tree
ig3 <- make_network(physeqAlb, "samples", distance = "jaccard", max.dist = 0.95)
plot_network(ig3, physeqAlb, type="samples", point_size = 5, label=NULL, color="dominant_tree", line_alpha = 0.05)

#######network analysis from example Lukas###

##subset Alb
#subset dataset to ASVs that occur in at least one percent of samples

phy_trans_A  = transform_sample_counts(physeqAlb, function(x) x / sum(x) )
phy_abundfilt_A = filter_taxa(phy_trans_A, function(x) sum(x) > .01, TRUE)
keep_names_A <- taxa_names(phy_abundfilt_A)

my_subset_A <- subset(otu_table(physeqAlb), rownames(otu_table(physeqAlb)) %in% keep_names_A)
phy_raw_abundfilt_A <- merge_phyloseq(my_subset_A, tax_table(physeqAlb), sample_data(physeqAlb))

#saveRDS(phy_raw_abundfilt_A, 'phy_fungi_raw_abundfilt.rds')
#check for location and execution on Malloy

# # this will be done on Malloy
# library(SpiecEasi)
# 
# phy_raw_abundfilt_A <- readRDS('phy_fungi_raw_abundfilt.rds')
# not required if execution on local machine
# 
pargs2 <- list(rep.num=50, seed=10010)
# not required if execution on local machine
# 
se_raw_abundfilt_A <- spiec.easi(phy_raw_abundfilt_A, method='mb', 
                                    lambda.min.ratio=1e-2, nlambda=70, 
                                    sel.criterion='bstars',  pulsar.select=TRUE,
                                    pulsar.params=pargs2)

#summary
se_raw_abundfilt_A$select$stars$summary

#get optimal parameters
getOptInd(se_raw_abundfilt_A)
# [1] 1
sum(getRefit(se_raw_abundfilt_A))/2


#network stability
getStability(se_raw_abundfilt_A)
 
# saveRDS(se_fun_raw_abundfilt, 'se_fun_raw_abundfilt.rds')
# not required if execution on local machine

# now back to local R for plotting the network

#se_fun_raw_abundfilt <- readRDS(here('se_fun_raw_abundfilt.rds'))
# not required if execution on local machine
fun.mb_A <- adj2igraph(getRefit(se_raw_abundfilt_A),  
                       vertex.attr=list(name=taxa_names(phy_raw_abundfilt_A)))

###
#plot the network
###

# convert to gephi 

fun.mb_A_gephi <- igraph.to.gexf(fun.mb_A)

write.gexf(fun.mb_A_gephi, output = here("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/fungi_A.gexf"))

#get the modules through gephi 

###
#get hub taxa, through Kleinbergs centrality
###


fun_between_A <- betweenness(fun.mb_A, directed = F)

fun_top5_between_A <- sort(fun_between_A, decreasing = T)[1:5]

names(fun_top5_between_A)

hub_taxa_fun_A <- subset_taxa(physeqAlb, taxa_names(physeqAlb) %in% names(fun_top5_between_A))
tax_table(hub_taxa_fun_A)

##subset Alb soil

#subset dataset to ASVs that occur in at least one percent of samples
phy_trans_As  = transform_sample_counts(physeqAlb_s, function(x) x / sum(x) )
phy_abundfilt_As = filter_taxa(phy_trans_As, function(x) sum(x) > .01, TRUE)
keep_names_As <- taxa_names(phy_abundfilt_As)
my_subset_As <- subset(otu_table(physeqAlb_s), rownames(otu_table(physeqAlb_s)) %in% keep_names_As)
phy_raw_abundfilt_As <- merge_phyloseq(my_subset_As, tax_table(physeqAlb_s), sample_data(physeqAlb_s))

pargs2 <- list(rep.num=50, seed=10010)

#SpiecEasi algorithm execution
se_raw_abundfilt_As <- spiec.easi(phy_raw_abundfilt_As, method='mb', 
                                  lambda.min.ratio=1e-2, nlambda=70, 
                                  sel.criterion='bstars',  pulsar.select=TRUE,
                                  pulsar.params=pargs2)

#summary
se_raw_abundfilt_As$select$stars$summary

#get optimal parameters
getOptInd(se_raw_abundfilt_As)
sum(getRefit(se_raw_abundfilt_As))/2

#network stability -> optimum close to 0.05
getStability(se_raw_abundfilt_As)

fun.mb_As <- adj2igraph(getRefit(se_raw_abundfilt_As),  
                        vertex.attr=list(name=taxa_names(phy_raw_abundfilt_As)))

##plot network
# convert to gephi 
fun.mb_As_gephi <- igraph.to.gexf(fun.mb_As)
write.gexf(fun.mb_As_gephi, output = here("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/fungi_As.gexf"))


##get hub taxa, through Kleinbergs centrality
fun_between_As <- betweenness(fun.mb_As, directed = F)
fun_top5_between_As <- sort(fun_between_As, decreasing = T)[1:5]
names(fun_top5_between_As)
hub_taxa_fun_As <- subset_taxa(physeqAlb_s, taxa_names(physeqAlb_s) %in% names(fun_top5_between_As))
tax_table(hub_taxa_fun_As)

##subset Alb bark

#subset dataset to ASVs that occur in at least one percent of samples
phy_trans_Ab  = transform_sample_counts(physeqAlb_b, function(x) x / sum(x) )
phy_abundfilt_Ab = filter_taxa(phy_trans_Ab, function(x) sum(x) > .01, TRUE)
keep_names_Ab <- taxa_names(phy_abundfilt_Ab)
my_subset_Ab <- subset(otu_table(physeqAlb_b), rownames(otu_table(physeqAlb_b)) %in% keep_names_Ab)
phy_raw_abundfilt_Ab <- merge_phyloseq(my_subset_Ab, tax_table(physeqAlb_b), sample_data(physeqAlb_b))

pargs2 <- list(rep.num=50, seed=10010)

#SpiecEasi algorithm execution
se_raw_abundfilt_Ab <- spiec.easi(phy_raw_abundfilt_Ab, method='mb', 
                                  lambda.min.ratio=1e-2, nlambda=70, 
                                  sel.criterion='bstars',  pulsar.select=TRUE,
                                  pulsar.params=pargs2)

#summary
se_raw_abundfilt_Ab$select$stars$summary

#get optimal parameters
getOptInd(se_raw_abundfilt_Ab)
sum(getRefit(se_raw_abundfilt_Ab))/2

#network stability -> optimum close to 0.05
getStability(se_raw_abundfilt_Ab)

fun.mb_Ab <- adj2igraph(getRefit(se_raw_abundfilt_Ab),  
                        vertex.attr=list(name=taxa_names(phy_raw_abundfilt_Ab)))

##plot network
# convert to gephi 
fun.mb_Ab_gephi <- igraph.to.gexf(fun.mb_Ab)
write.gexf(fun.mb_Ab_gephi, output = here("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/fungi_Ab.gexf"))


##get hub taxa, through Kleinbergs centrality
fun_between_Ab <- betweenness(fun.mb_Ab, directed = F)
fun_top5_between_Ab <- sort(fun_between_Ab, decreasing = T)[1:5]
names(fun_top5_between_Ab)
hub_taxa_fun_Ab <- subset_taxa(physeqAlb_b, taxa_names(physeqAlb_b) %in% names(fun_top5_between_Ab))
tax_table(hub_taxa_fun_Ab)


##subset Alb bark Fagus

#subset dataset to ASVs that occur in at least one percent of samples
phy_trans_AbF  = transform_sample_counts(physeqAlb_bF, function(x) x / sum(x) )
phy_abundfilt_AbF = filter_taxa(phy_trans_AbF, function(x) sum(x) > .01, TRUE)
keep_names_AbF <- taxa_names(phy_abundfilt_AbF)
my_subset_AbF <- subset(otu_table(physeqAlb_bF), rownames(otu_table(physeqAlb_bF)) %in% keep_names_AbF)
phy_raw_abundfilt_AbF <- merge_phyloseq(my_subset_AbF, tax_table(physeqAlb_bF), sample_data(physeqAlb_bF))

pargs2 <- list(rep.num=50, seed=10010)

#SpiecEasi algorithm execution
se_raw_abundfilt_AbF <- spiec.easi(phy_raw_abundfilt_AbF, method='mb', 
                                   lambda.min.ratio=1e-2, nlambda=70, 
                                   sel.criterion='bstars',  pulsar.select=TRUE,
                                   pulsar.params=pargs2)

#summary
se_raw_abundfilt_AbF$select$stars$summary

#get optimal parameters
getOptInd(se_raw_abundfilt_AbF)
sum(getRefit(se_raw_abundfilt_AbF))/2

#network stability -> optimum close to 0.05
getStability(se_raw_abundfilt_AbF)

fun.mb_AbF <- adj2igraph(getRefit(se_raw_abundfilt_AbF),  
                         vertex.attr=list(name=taxa_names(phy_raw_abundfilt_AbF)))

##plot network
# convert to gephi 
fun.mb_AbF_gephi <- igraph.to.gexf(fun.mb_AbF)
write.gexf(fun.mb_AbF_gephi, output = here("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/fungi_AbF.gexf"))


##get hub taxa, through Kleinbergs centrality
fun_between_AbF <- betweenness(fun.mb_AbF, directed = F)
fun_top5_between_AbF <- sort(fun_between_AbF, decreasing = T)[1:5]
names(fun_top5_between_AbF)
hub_taxa_fun_AbF <- subset_taxa(physeqAlb_bF, taxa_names(physeqAlb_bF) %in% names(fun_top5_between_AbF))
tax_table(hub_taxa_fun_AbF)
hubtax_AbF <- data.frame(tax_table(hub_taxa_fun_AbF))
hubtax_AbF <- rownames_to_column(hubtax_AbF, var = "ASV_ID") %>% as_tibble()
hubtax_AbF
write.csv(hubtax_AbF, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/hubtax_AbF.csv")


##subset Alb soil Fagus

#subset dataset to ASVs that occur in at least one percent of samples
phy_trans_AsF  = transform_sample_counts(physeqAlb_sF, function(x) x / sum(x) )
phy_abundfilt_AsF = filter_taxa(phy_trans_AsF, function(x) sum(x) > .01, TRUE)
keep_names_AsF <- taxa_names(phy_abundfilt_AsF)
my_subset_AsF <- subset(otu_table(physeqAlb_sF), rownames(otu_table(physeqAlb_sF)) %in% keep_names_AsF)
phy_raw_abundfilt_AsF <- merge_phyloseq(my_subset_AsF, tax_table(physeqAlb_sF), sample_data(physeqAlb_sF))

pargs2 <- list(rep.num=50, seed=10010)

#SpiecEasi algorithm execution
se_raw_abundfilt_AsF <- spiec.easi(phy_raw_abundfilt_AsF, method='mb', 
                                   lambda.min.ratio=1e-2, nlambda=70, 
                                   sel.criterion='bstars',  pulsar.select=TRUE,
                                   pulsar.params=pargs2)

#summary
se_raw_abundfilt_AsF$select$stars$summary

#get optimal parameters
getOptInd(se_raw_abundfilt_AsF)
sum(getRefit(se_raw_abundfilt_AsF))/2

#network stability -> optimum close to 0.05
getStability(se_raw_abundfilt_AsF)

fun.mb_AsF <- adj2igraph(getRefit(se_raw_abundfilt_AsF),  
                         vertex.attr=list(name=taxa_names(phy_raw_abundfilt_AsF)))

##plot network
# convert to gephi 
fun.mb_AsF_gephi <- igraph.to.gexf(fun.mb_AsF)
write.gexf(fun.mb_AsF_gephi, output = here("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/fungi_AsF.gexf"))


##get hub taxa, through Kleinbergs centrality
fun_between_AsF <- betweenness(fun.mb_AsF, directed = F)
fun_top5_between_AsF <- sort(fun_between_AsF, decreasing = T)[1:5]
names(fun_top5_between_AsF)
hub_taxa_fun_AsF <- subset_taxa(physeqAlb_sF, taxa_names(physeqAlb_sF) %in% names(fun_top5_between_AsF))
tax_table(hub_taxa_fun_AsF)
hubtax_AsF <- data.frame(tax_table(hub_taxa_fun_AsF))
hubtax_AsF <- rownames_to_column(hubtax_AsF, var = "ASV_ID") %>% as_tibble()
hubtax_AsF
write.csv(hubtax_AsF, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/hubtax_AsF.csv")

##subset Alb bark Picea

#subset dataset to ASVs that occur in at least one percent of samples
phy_trans_AbP  = transform_sample_counts(physeqAlb_bP, function(x) x / sum(x) )
phy_abundfilt_AbP = filter_taxa(phy_trans_AbP, function(x) sum(x) > .01, TRUE)
keep_names_AbP <- taxa_names(phy_abundfilt_AbP)
my_subset_AbP <- subset(otu_table(physeqAlb_bP), rownames(otu_table(physeqAlb_bP)) %in% keep_names_AbP)
phy_raw_abundfilt_AbP <- merge_phyloseq(my_subset_AbP, tax_table(physeqAlb_bP), sample_data(physeqAlb_bP))

pargs2 <- list(rep.num=50, seed=10010)

#SpiecEasi algorithm execution
se_raw_abundfilt_AbP <- spiec.easi(phy_raw_abundfilt_AbP, method='mb', 
                                   lambda.min.ratio=1e-2, nlambda=70, 
                                   sel.criterion='bstars',  pulsar.select=TRUE,
                                   pulsar.params=pargs2)

#summary
se_raw_abundfilt_AbP$select$stars$summary

#get optimal parameters
getOptInd(se_raw_abundfilt_AbP)
sum(getRefit(se_raw_abundfilt_AbP))/2

#network stability -> optimum close to 0.05
getStability(se_raw_abundfilt_AbP)

fun.mb_AbP <- adj2igraph(getRefit(se_raw_abundfilt_AbP),  
                         vertex.attr=list(name=taxa_names(phy_raw_abundfilt_AbP)))

##plot network
# convert to gephi 
fun.mb_AbP_gephi <- igraph.to.gexf(fun.mb_AbP)
write.gexf(fun.mb_AbP_gephi, output = here("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/fungi_AbP.gexf"))


##get hub taxa, through Kleinbergs centrality
fun_between_AbP <- betweenness(fun.mb_AbP, directed = F)
fun_top5_between_AbP <- sort(fun_between_AbP, decreasing = T)[1:5]
names(fun_top5_between_AbP)
hub_taxa_fun_AbP <- subset_taxa(physeqAlb_bP, taxa_names(physeqAlb_bP) %in% names(fun_top5_between_AbP))
tax_table(hub_taxa_fun_AbP)
hubtax_AbP <- data.frame(tax_table(hub_taxa_fun_AbP))
hubtax_AbP <- rownames_to_column(hubtax_AbP, var = "ASV_ID") %>% as_tibble()
hubtax_AbP
write.csv(hubtax_AbP, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/hubtax_AbP.csv")

##subset Alb soil Picea

#subset dataset to ASVs that occur in at least one percent of samples
phy_trans_AsP  = transform_sample_counts(physeqAlb_sP, function(x) x / sum(x) )
phy_abundfilt_AsP = filter_taxa(phy_trans_AsP, function(x) sum(x) > .01, TRUE)
keep_names_AsP <- taxa_names(phy_abundfilt_AsP)
my_subset_AsP <- subset(otu_table(physeqAlb_sP), rownames(otu_table(physeqAlb_sP)) %in% keep_names_AsP)
phy_raw_abundfilt_AsP <- merge_phyloseq(my_subset_AsP, tax_table(physeqAlb_sP), sample_data(physeqAlb_sP))

pargs2 <- list(rep.num=50, seed=10010)

#SpiecEasi algorithm execution
se_raw_abundfilt_AsP <- spiec.easi(phy_raw_abundfilt_AsP, method='mb', 
                                   lambda.min.ratio=1e-2, nlambda=70, 
                                   sel.criterion='bstars',  pulsar.select=TRUE,
                                   pulsar.params=pargs2)

#summary
se_raw_abundfilt_AsP$select$stars$summary

#get optimal parameters
getOptInd(se_raw_abundfilt_AsP)
sum(getRefit(se_raw_abundfilt_AsP))/2

#network stability -> optimum close to 0.05
getStability(se_raw_abundfilt_AsP)

fun.mb_AsP <- adj2igraph(getRefit(se_raw_abundfilt_AsP),  
                         vertex.attr=list(name=taxa_names(phy_raw_abundfilt_AsP)))

##plot network
# convert to gephi 
fun.mb_AsP_gephi <- igraph.to.gexf(fun.mb_AsP)
write.gexf(fun.mb_AsP_gephi, output = here("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/fungi_AsP.gexf"))


##get hub taxa, through Kleinbergs centrality
fun_between_AsP <- betweenness(fun.mb_AsP, directed = F)
fun_top5_between_AsP <- sort(fun_between_AsP, decreasing = T)[1:5]
names(fun_top5_between_AsP)
hub_taxa_fun_AsP <- subset_taxa(physeqAlb_sP, taxa_names(physeqAlb_sP) %in% names(fun_top5_between_AsP))
tax_table(hub_taxa_fun_AsP)
hubtax_AsP <- data.frame(tax_table(hub_taxa_fun_AsP))
hubtax_AsP <- rownames_to_column(hubtax_AsP, var = "ASV_ID") %>% as_tibble()
hubtax_AsP
write.csv(hubtax_AsP, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/hubtax_AsP.csv")



################Schorfheide#######################

#######network analysis SpiecEasi###

##subset Schorfheide

#subset dataset to ASVs that occur in at least one percent of samples
phy_trans_S  = transform_sample_counts(physeqSchorf, function(x) x / sum(x) )
phy_abundfilt_S = filter_taxa(phy_trans_S, function(x) sum(x) > .01, TRUE)
keep_names_S <- taxa_names(phy_abundfilt_S)
my_subset_S <- subset(otu_table(physeqSchorf), rownames(otu_table(physeqSchorf)) %in% keep_names_S)
phy_raw_abundfilt_S <- merge_phyloseq(my_subset_S, tax_table(physeqSchorf), sample_data(physeqSchorf))

#pargs2 <- list(rep.num=50, seed=10010) (executed for Alb)

#SpiecEasi algorithm execution
se_raw_abundfilt_S <- spiec.easi(phy_raw_abundfilt_S, method='mb', 
                                 lambda.min.ratio=1e-2, nlambda=70, 
                                 sel.criterion='bstars',  pulsar.select=TRUE,
                                 pulsar.params=pargs2)

#summary
se_raw_abundfilt_S$select$stars$summary

#get optimal parameters
getOptInd(se_raw_abundfilt_S)
sum(getRefit(se_raw_abundfilt_S))/2

#network stability -> optimum close to 0.05
getStability(se_raw_abundfilt_S)

fun.mb_S <- adj2igraph(getRefit(se_raw_abundfilt_S),  
                       vertex.attr=list(name=taxa_names(phy_raw_abundfilt_S)))

##plot network
# convert to gephi 
fun.mb_S_gephi <- igraph.to.gexf(fun.mb_S)
write.gexf(fun.mb_S_gephi, output = here("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/fungi_S.gexf"))


##get hub taxa, through Kleinbergs centrality
fun_between_S <- betweenness(fun.mb_S, directed = F)
fun_top5_between_S <- sort(fun_between_S, decreasing = T)[1:5]
names(fun_top5_between_S)
hub_taxa_fun_S <- subset_taxa(physeqSchorf, taxa_names(physeqSchorf) %in% names(fun_top5_between_S))
tax_table(hub_taxa_fun_S)

##subset Schorfheide soil

#subset dataset to ASVs that occur in at least one percent of samples
phy_trans_Ss  = transform_sample_counts(physeqSchorf_s, function(x) x / sum(x) )
phy_abundfilt_Ss = filter_taxa(phy_trans_Ss, function(x) sum(x) > .01, TRUE)
keep_names_Ss <- taxa_names(phy_abundfilt_Ss)
my_subset_Ss <- subset(otu_table(physeqSchorf_s), rownames(otu_table(physeqSchorf_s)) %in% keep_names_Ss)
phy_raw_abundfilt_Ss <- merge_phyloseq(my_subset_Ss, tax_table(physeqSchorf_s), sample_data(physeqSchorf_s))

pargs2 <- list(rep.num=50, seed=10010)

#SpiecEasi algorithm execution
se_raw_abundfilt_Ss <- spiec.easi(phy_raw_abundfilt_Ss, method='mb', 
                                  lambda.min.ratio=1e-2, nlambda=70, 
                                  sel.criterion='bstars',  pulsar.select=TRUE,
                                  pulsar.params=pargs2)

#summary
se_raw_abundfilt_Ss$select$stars$summary

#get optimal parameters
getOptInd(se_raw_abundfilt_Ss)
sum(getRefit(se_raw_abundfilt_Ss))/2

#network stability -> optimum close to 0.05
getStability(se_raw_abundfilt_Ss)

fun.mb_Ss <- adj2igraph(getRefit(se_raw_abundfilt_Ss),  
                        vertex.attr=list(name=taxa_names(phy_raw_abundfilt_Ss)))

##plot network
# convert to gephi 
fun.mb_Ss_gephi <- igraph.to.gexf(fun.mb_Ss)
write.gexf(fun.mb_Ss_gephi, output = here("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/fungi_Ss.gexf"))


##get hub taxa, through Kleinbergs centrality
fun_between_Ss <- betweenness(fun.mb_Ss, directed = F)
fun_top5_between_Ss <- sort(fun_between_Ss, decreasing = T)[1:5]
names(fun_top5_between_Ss)
hub_taxa_fun_Ss <- subset_taxa(physeqSchorf_s, taxa_names(physeqSchorf_s) %in% names(fun_top5_between_Ss))
tax_table(hub_taxa_fun_Ss)


##subset Schorfheide bark

#subset dataset to ASVs that occur in at least one percent of samples
phy_trans_Sb  = transform_sample_counts(physeqSchorf_b, function(x) x / sum(x) )
phy_abundfilt_Sb = filter_taxa(phy_trans_Sb, function(x) sum(x) > .01, TRUE)
keep_names_Sb <- taxa_names(phy_abundfilt_Sb)
my_subset_Sb <- subset(otu_table(physeqSchorf_b), rownames(otu_table(physeqSchorf_b)) %in% keep_names_Sb)
phy_raw_abundfilt_Sb <- merge_phyloseq(my_subset_Sb, tax_table(physeqSchorf_b), sample_data(physeqSchorf_b))

pargs2 <- list(rep.num=50, seed=10010)

#SpiecEasi algorithm execution
se_raw_abundfilt_Sb <- spiec.easi(phy_raw_abundfilt_Sb, method='mb', 
                                  lambda.min.ratio=1e-2, nlambda=70, 
                                  sel.criterion='bstars',  pulsar.select=TRUE,
                                  pulsar.params=pargs2)

#summary
se_raw_abundfilt_Sb$select$stars$summary

#get optimal parameters
getOptInd(se_raw_abundfilt_Sb)
sum(getRefit(se_raw_abundfilt_Sb))/2

#network stability -> optimum close to 0.05
getStability(se_raw_abundfilt_Sb)

fun.mb_Sb <- adj2igraph(getRefit(se_raw_abundfilt_Sb),  
                        vertex.attr=list(name=taxa_names(phy_raw_abundfilt_Sb)))

##plot network
# convert to gephi 
fun.mb_Sb_gephi <- igraph.to.gexf(fun.mb_Sb)
write.gexf(fun.mb_Sb_gephi, output = here("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/fungi_Sb.gexf"))


##get hub taxa, through Kleinbergs centrality
fun_between_Sb <- betweenness(fun.mb_Sb, directed = F)
fun_top5_between_Sb <- sort(fun_between_Sb, decreasing = T)[1:5]
names(fun_top5_between_Sb)
hub_taxa_fun_Sb <- subset_taxa(physeqSchorf_b, taxa_names(physeqSchorf_b) %in% names(fun_top5_between_Sb))
tax_table(hub_taxa_fun_Sb)

##subset Schorfheide bark Fagus

#subset dataset to ASVs that occur in at least one percent of samples
phy_trans_SbF  = transform_sample_counts(physeqSchorf_bF, function(x) x / sum(x) )
phy_abundfilt_SbF = filter_taxa(phy_trans_SbF, function(x) sum(x) > .01, TRUE)
keep_names_SbF <- taxa_names(phy_abundfilt_SbF)
my_subset_SbF <- subset(otu_table(physeqSchorf_bF), rownames(otu_table(physeqSchorf_bF)) %in% keep_names_SbF)
phy_raw_abundfilt_SbF <- merge_phyloseq(my_subset_SbF, tax_table(physeqSchorf_bF), sample_data(physeqSchorf_bF))

pargs2 <- list(rep.num=50, seed=10010)

#SpiecEasi algorithm execution
se_raw_abundfilt_SbF <- spiec.easi(phy_raw_abundfilt_SbF, method='mb', 
                                   lambda.min.ratio=1e-2, nlambda=70, 
                                   sel.criterion='bstars',  pulsar.select=TRUE,
                                   pulsar.params=pargs2)

#summary
se_raw_abundfilt_SbF$select$stars$summary

#get optimal parameters
getOptInd(se_raw_abundfilt_SbF)
sum(getRefit(se_raw_abundfilt_SbF))/2

#network stability -> optimum close to 0.05
getStability(se_raw_abundfilt_SbF)

fun.mb_SbF <- adj2igraph(getRefit(se_raw_abundfilt_SbF),  
                         vertex.attr=list(name=taxa_names(phy_raw_abundfilt_SbF)))

##plot network
# convert to gephi 
fun.mb_SbF_gephi <- igraph.to.gexf(fun.mb_SbF)
write.gexf(fun.mb_SbF_gephi, output = here("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/fungi_SbF.gexf"))


##get hub taxa, through Kleinbergs centrality
fun_between_SbF <- betweenness(fun.mb_SbF, directed = F)
fun_top5_between_SbF <- sort(fun_between_SbF, decreasing = T)[1:5]
names(fun_top5_between_SbF)
hub_taxa_fun_SbF <- subset_taxa(physeqSchorf_bF, taxa_names(physeqSchorf_bF) %in% names(fun_top5_between_SbF))
tax_table(hub_taxa_fun_SbF)
hubtax_SbF <- data.frame(tax_table(hub_taxa_fun_SbF))
hubtax_SbF <- rownames_to_column(hubtax_SbF, var = "ASV_ID") %>% as_tibble()
hubtax_SbF
write.csv(hubtax_SbF, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/hubtax_SbF.csv")


##subset Schorfheide soil Fagus

#subset dataset to ASVs that occur in at least one percent of samples
phy_trans_SsF  = transform_sample_counts(physeqSchorf_sF, function(x) x / sum(x) )
phy_abundfilt_SsF = filter_taxa(phy_trans_SsF, function(x) sum(x) > .01, TRUE)
keep_names_SsF <- taxa_names(phy_abundfilt_SsF)
my_subset_SsF <- subset(otu_table(physeqSchorf_sF), rownames(otu_table(physeqSchorf_sF)) %in% keep_names_SsF)
phy_raw_abundfilt_SsF <- merge_phyloseq(my_subset_SsF, tax_table(physeqSchorf_sF), sample_data(physeqSchorf_sF))

otu_table(physeqSchorf_bF)

pargs2 <- list(rep.num=50, seed=10010)

#SpiecEasi algorithm execution
se_raw_abundfilt_SsF <- spiec.easi(phy_raw_abundfilt_SsF, method='mb', 
                                   lambda.min.ratio=1e-2, nlambda=70, 
                                   sel.criterion='bstars',  pulsar.select=TRUE,
                                   pulsar.params=pargs2)

#summary
se_raw_abundfilt_SsF$select$stars$summary

#get optimal parameters
getOptInd(se_raw_abundfilt_SsF)
sum(getRefit(se_raw_abundfilt_SsF))/2

#network stability -> optimum close to 0.05
getStability(se_raw_abundfilt_SsF)

fun.mb_SsF <- adj2igraph(getRefit(se_raw_abundfilt_SsF),  
                         vertex.attr=list(name=taxa_names(phy_raw_abundfilt_SsF)))

##plot network
# convert to gephi 
fun.mb_SsF_gephi <- igraph.to.gexf(fun.mb_SsF)
write.gexf(fun.mb_SsF_gephi, output = here("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/fungi_SsF.gexf"))


##get hub taxa, through Kleinbergs centrality
fun_between_SsF <- betweenness(fun.mb_SsF, directed = F)
fun_top5_between_SsF <- sort(fun_between_SsF, decreasing = T)[1:5]
names(fun_top5_between_SsF)
hub_taxa_fun_SsF <- subset_taxa(physeqSchorf_sF, taxa_names(physeqSchorf_sF) %in% names(fun_top5_between_SsF))
tax_table(hub_taxa_fun_SsF)
hubtax_SsF <- data.frame(tax_table(hub_taxa_fun_SsF))
hubtax_SsF <- rownames_to_column(hubtax_SsF, var = "ASV_ID") %>% as_tibble()
hubtax_SsF
write.csv(hubtax_SsF, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/hubtax_SsF.csv")


##subset Schorfheide bark Pinus

#subset dataset to ASVs that occur in at least one percent of samples
phy_trans_SbP  = transform_sample_counts(physeqSchorf_bP, function(x) x / sum(x) )
phy_abundfilt_SbP = filter_taxa(phy_trans_SbP, function(x) sum(x) > .01, TRUE)
keep_names_SbP <- taxa_names(phy_abundfilt_SbP)
my_subset_SbP <- subset(otu_table(physeqSchorf_bP), rownames(otu_table(physeqSchorf_bP)) %in% keep_names_SbP)
phy_raw_abundfilt_SbP <- merge_phyloseq(my_subset_SbP, tax_table(physeqSchorf_bP), sample_data(physeqSchorf_bP))

pargs2 <- list(rep.num=50, seed=10010)

#SpiecEasi algorithm execution
se_raw_abundfilt_SbP <- spiec.easi(phy_raw_abundfilt_SbP, method='mb', 
                                   lambda.min.ratio=1e-2, nlambda=70, 
                                   sel.criterion='bstars',  pulsar.select=TRUE,
                                   pulsar.params=pargs2)

#summary
se_raw_abundfilt_SbP$select$stars$summary

#get optimal parameters
getOptInd(se_raw_abundfilt_SbP)
sum(getRefit(se_raw_abundfilt_SbP))/2

#network stability -> optimum close to 0.05
getStability(se_raw_abundfilt_SbP)

fun.mb_SbP <- adj2igraph(getRefit(se_raw_abundfilt_SbP),  
                         vertex.attr=list(name=taxa_names(phy_raw_abundfilt_SbP)))

##plot network
# convert to gephi 
fun.mb_SbP_gephi <- igraph.to.gexf(fun.mb_SbP)
write.gexf(fun.mb_SbP_gephi, output = here("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/fungi_SbP.gexf"))


##get hub taxa, through Kleinbergs centrality
fun_between_SbP <- betweenness(fun.mb_SbP, directed = F)
fun_top5_between_SbP <- sort(fun_between_SbP, decreasing = T)[1:5]
names(fun_top5_between_SbP)
hub_taxa_fun_SbP <- subset_taxa(physeqSchorf_bP, taxa_names(physeqSchorf_bP) %in% names(fun_top5_between_SbP))
tax_table(hub_taxa_fun_SbP)
hubtax_SbP <- data.frame(tax_table(hub_taxa_fun_SbP))
hubtax_SbP <- rownames_to_column(hubtax_SbP, var = "ASV_ID") %>% as_tibble()
hubtax_SbP
write.csv(hubtax_SbP, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/hubtax_SbP.csv")


##subset Schorfheide soil Pinus

#subset dataset to ASVs that occur in at least one percent of samples
phy_trans_SsP  = transform_sample_counts(physeqSchorf_sP, function(x) x / sum(x) )
phy_abundfilt_SsP = filter_taxa(phy_trans_SsP, function(x) sum(x) > .01, TRUE)
keep_names_SsP <- taxa_names(phy_abundfilt_SsP)
my_subset_SsP <- subset(otu_table(physeqSchorf_sP), rownames(otu_table(physeqSchorf_sP)) %in% keep_names_SsP)
phy_raw_abundfilt_SsP <- merge_phyloseq(my_subset_SsP, tax_table(physeqSchorf_sP), sample_data(physeqSchorf_sP))

pargs2 <- list(rep.num=50, seed=10010)

#SpiecEasi algorithm execution
se_raw_abundfilt_SsP <- spiec.easi(phy_raw_abundfilt_SsP, method='mb', 
                                   lambda.min.ratio=1e-2, nlambda=70, 
                                   sel.criterion='bstars',  pulsar.select=TRUE,
                                   pulsar.params=pargs2)

#summary
se_raw_abundfilt_SsP$select$stars$summary

#get optimal parameters
getOptInd(se_raw_abundfilt_SsP)
sum(getRefit(se_raw_abundfilt_SsP))/2

#network stability -> optimum close to 0.05
getStability(se_raw_abundfilt_SsP)

fun.mb_SsP <- adj2igraph(getRefit(se_raw_abundfilt_SsP),  
                         vertex.attr=list(name=taxa_names(phy_raw_abundfilt_SsP)))

##plot network
# convert to gephi 
fun.mb_SsP_gephi <- igraph.to.gexf(fun.mb_SsP)
write.gexf(fun.mb_SsP_gephi, output = here("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/fungi_SsP.gexf"))


##get hub taxa, through Kleinbergs centrality
fun_between_SsP <- betweenness(fun.mb_SsP, directed = F)
fun_top5_between_SsP <- sort(fun_between_SsP, decreasing = T)[1:5]
names(fun_top5_between_SsP)
hub_taxa_fun_SsP <- subset_taxa(physeqSchorf_sP, taxa_names(physeqSchorf_sP) %in% names(fun_top5_between_SsP))
tax_table(hub_taxa_fun_SsP)
hubtax_SsP <- data.frame(tax_table(hub_taxa_fun_SsP))
hubtax_SsP <- rownames_to_column(hubtax_SsP, var = "ASV_ID") %>% as_tibble()
hubtax_SsP
write.csv(hubtax_SsP, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/hubtax_SsP.csv")



#######network analysis Physalia Metabarcoding course###

#by taxa
ig4 <- make_network(physeqSchorf, "taxa", distance = "jaccard", max.dist = 0.95)
plot_network(ig4, physeqSchorf, type="taxa", point_size = 5, label=NULL, color="Phylum", line_alpha = 0.05)

#network samples - substrate
ig5 <- make_network(physeqSchorf, "samples", distance = "jaccard", max.dist = 0.95)
plot_network(ig5, physeqSchorf, type="samples", point_size = 5, label=NULL, color="substrate", line_alpha = 0.05)

#network samples - dominant_tree
ig6 <- make_network(physeqSchorf, "samples", distance = "jaccard", max.dist = 0.95)
plot_network(ig6, physeqSchorf, type="samples", point_size = 5, label=NULL, color="dominant_tree", line_alpha = 0.05)


############approach to load modularity classes to physeq object (Lukas)########

######Schorfheide#####

## check what modules are based on
# Read in the data on the modules obtained from Gephi. 
modules_fun_Schorf <- read.csv2(("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/nw_modules_schorf.csv"))
# Transform ASV information to dataframe to make handling easier.
otu_tab_fun_S <- as.data.frame(otu_table(phy_abundfilt_S))
# Calculate the total read count of each ASV.
otu_tab_fun_S$total <- rowSums(otu_tab_fun_S)
# Keep only the total read count.
otu_tab_fun_S <- otu_tab_fun_S %>%
  select(.,total) %>%
  rownames_to_column('Label')
# Merge module table and read count information.
otu_modules_S <- inner_join(otu_tab_fun_S, modules_fun_Schorf)
# Initiate an empty dataframe were we can store the total read counts per module.
otu_mod_fun_S <- data.frame(Label = otu_modules_S$Label,
                            Mod0 = rep(0, 596),
                            Mod1 = rep(0, 596),
                            Mod2 = rep(0, 596),
                            Mod3 = rep(0, 596),
                            Mod4 = rep(0, 596))

# If the ASV belongs to a module then write its read count, if not write zero.
otu_mod_fun_S$Mod0 <- if_else(otu_modules_S$modularity_class == 0, otu_modules_S$total, 0)
otu_mod_fun_S$Mod1 <- if_else(otu_modules_S$modularity_class == 1, otu_modules_S$total, 0)
otu_mod_fun_S$Mod2 <- if_else(otu_modules_S$modularity_class == 2, otu_modules_S$total, 0)
otu_mod_fun_S$Mod3 <- if_else(otu_modules_S$modularity_class == 3, otu_modules_S$total, 0)
otu_mod_fun_S$Mod4 <- if_else(otu_modules_S$modularity_class == 4, otu_modules_S$total, 0)


# set the ASV name back to the rownames.
rownames(otu_mod_fun_S) <- otu_mod_fun_S$Label
otu_mod_fun_S$Label <- NULL
# Create the phyloseq object to easily obtain the taxonomic information.
otu_mat_mod_fun_S <- as.matrix(otu_mod_fun_S)
tax_mat_mod_fun_S <- as.matrix(tax_table(phy_abundfilt_S))
OTU_mod_fun_S <- otu_table(otu_mat_mod_fun_S, taxa_are_rows = T)
TAX_mod_fun_S <- tax_table(tax_mat_mod_fun_S)
phy_modules_fun_S <- phyloseq(OTU_mod_fun_S, TAX_mod_fun_S, sampledataSchorf)
# Find out which modules contain more than 10 ASVs.
module_asv_count_fun_S <- colSums(otu_mat_mod_fun_S != 0)
which(module_asv_count_fun_S > 10)
sort(module_asv_count_fun_S)


####
# Find out which taxa are the top ones in the modules with more than ten ASVs.
####
# Subset by Module name.
phy_funS_mod0 <- prune_samples(sample_names(phy_modules_fun_S) == 'Mod0', phy_modules_fun_S)
phy_funS_mod1 <- prune_samples(sample_names(phy_modules_fun_S) == 'Mod1', phy_modules_fun_S)
phy_funS_mod2 <- prune_samples(sample_names(phy_modules_fun_S) == 'Mod2', phy_modules_fun_S)
phy_funS_mod3 <- prune_samples(sample_names(phy_modules_fun_S) == 'Mod3', phy_modules_fun_S)
phy_funS_mod4 <- prune_samples(sample_names(phy_modules_fun_S) == 'Mod4', phy_modules_fun_S)

# Remove any taxa without reads.
phy_funS_mod0_clean <- prune_taxa(taxa_sums(phy_funS_mod0) > 0, phy_funS_mod0)
phy_funS_mod1_clean <- prune_taxa(taxa_sums(phy_funS_mod1) > 0, phy_funS_mod1)
phy_funS_mod2_clean <- prune_taxa(taxa_sums(phy_funS_mod2) > 0, phy_funS_mod2)
phy_funS_mod3_clean <- prune_taxa(taxa_sums(phy_funS_mod3) > 0, phy_funS_mod3)
phy_funS_mod4_clean <- prune_taxa(taxa_sums(phy_funS_mod4) > 0, phy_funS_mod4)

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 0.
top_funS_mod0 <- get_top_taxa(phy_funS_mod0_clean, 1, discard_other = T)
tax_table(top_funS_mod0)
abundances(top_funS_mod0) / sum(abundances(phy_funS_mod0_clean))
# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 1.
top_funS_mod1 <- get_top_taxa(phy_funS_mod1_clean, 1, discard_other = T)
tax_table(top_funS_mod1)
abundances(top_funS_mod1) / sum(abundances(phy_funS_mod1_clean))
# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 2.
top_funS_mod2 <- get_top_taxa(phy_funS_mod2_clean, 1, discard_other = T)
tax_table(top_funS_mod2)
abundances(top_funS_mod2) / sum(abundances(phy_funS_mod2_clean))
# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 3.
top_funS_mod3 <- get_top_taxa(phy_funS_mod3_clean, 1, discard_other = T)
tax_table(top_funS_mod3)
abundances(top_funS_mod3) / sum(abundances(phy_funS_mod3_clean))
# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 4.
top_funS_mod4 <- get_top_taxa(phy_funS_mod4_clean, 1, discard_other = T)
tax_table(top_funS_mod4)
abundances(top_funS_mod4) / sum(abundances(phy_funS_mod4_clean))

phy_modules_Schorf <- phyloseq(OTU_mod_fun_S,TAX,sampledataSchorf)

######20220624 Schorfheide with new network#####

## check what modules are based on
# Read in the data on the modules obtained from Gephi. 
modules_fun_Schorf <- read.csv2(("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/nw_modules_schorf_2.csv"))
# Transform ASV information to dataframe to make handling easier.
otu_tab_fun_S <- as.data.frame(otu_table(phy_abundfilt_S))
# Calculate the total read count of each ASV.
otu_tab_fun_S$total <- rowSums(otu_tab_fun_S)
# Keep only the total read count.
otu_tab_fun_S <- otu_tab_fun_S %>%
  select(.,total) %>%
  rownames_to_column('Label')
# Merge module table and read count information.
otu_modules_S <- inner_join(otu_tab_fun_S, modules_fun_Schorf)
# Initiate an empty dataframe were we can store the total read counts per module.
otu_mod_fun_S <- data.frame(Label = otu_modules_S$Label,
                            Mod0 = rep(0, 596),
                            Mod1 = rep(0, 596),
                            Mod2 = rep(0, 596),
                            Mod3 = rep(0, 596))
                           # Mod4 = rep(0, 596))

# If the ASV belongs to a module then write its read count, if not write zero.
otu_mod_fun_S$Mod0 <- if_else(otu_modules_S$modularity_class == 0, otu_modules_S$total, 0)
otu_mod_fun_S$Mod1 <- if_else(otu_modules_S$modularity_class == 1, otu_modules_S$total, 0)
otu_mod_fun_S$Mod2 <- if_else(otu_modules_S$modularity_class == 2, otu_modules_S$total, 0)
otu_mod_fun_S$Mod3 <- if_else(otu_modules_S$modularity_class == 3, otu_modules_S$total, 0)
#otu_mod_fun_S$Mod4 <- if_else(otu_modules_S$modularity_class == 4, otu_modules_S$total, 0)


# set the ASV name back to the rownames.
rownames(otu_mod_fun_S) <- otu_mod_fun_S$Label
otu_mod_fun_S$Label <- NULL
# Create the phyloseq object to easily obtain the taxonomic information.
otu_mat_mod_fun_S <- as.matrix(otu_mod_fun_S)
tax_mat_mod_fun_S <- as.matrix(tax_table(phy_abundfilt_S))
OTU_mod_fun_S <- otu_table(otu_mat_mod_fun_S, taxa_are_rows = T)
TAX_mod_fun_S <- tax_table(tax_mat_mod_fun_S)
phy_modules_fun_S <- phyloseq(OTU_mod_fun_S, TAX_mod_fun_S, sampledataSchorf)
# Find out which modules contain more than 10 ASVs.
module_asv_count_fun_S <- colSums(otu_mat_mod_fun_S != 0)
which(module_asv_count_fun_S > 10)
sort(module_asv_count_fun_S)


####
# Find out which taxa are the top ones in the modules with more than ten ASVs.
####
# Subset by Module name.
phy_funS_mod0 <- prune_samples(sample_names(phy_modules_fun_S) == 'Mod0', phy_modules_fun_S)
phy_funS_mod1 <- prune_samples(sample_names(phy_modules_fun_S) == 'Mod1', phy_modules_fun_S)
phy_funS_mod2 <- prune_samples(sample_names(phy_modules_fun_S) == 'Mod2', phy_modules_fun_S)
phy_funS_mod3 <- prune_samples(sample_names(phy_modules_fun_S) == 'Mod3', phy_modules_fun_S)
#phy_funS_mod4 <- prune_samples(sample_names(phy_modules_fun_S) == 'Mod4', phy_modules_fun_S)

# Remove any taxa without reads.
phy_funS_mod0_clean <- prune_taxa(taxa_sums(phy_funS_mod0) > 0, phy_funS_mod0)
phy_funS_mod1_clean <- prune_taxa(taxa_sums(phy_funS_mod1) > 0, phy_funS_mod1)
phy_funS_mod2_clean <- prune_taxa(taxa_sums(phy_funS_mod2) > 0, phy_funS_mod2)
phy_funS_mod3_clean <- prune_taxa(taxa_sums(phy_funS_mod3) > 0, phy_funS_mod3)
#phy_funS_mod4_clean <- prune_taxa(taxa_sums(phy_funS_mod4) > 0, phy_funS_mod4)

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 0.
top_funS_mod0 <- get_top_taxa(phy_funS_mod0_clean, 1, discard_other = T)
tax_table(top_funS_mod0)
abundances(top_funS_mod0) / sum(abundances(phy_funS_mod0_clean))
# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 1.
top_funS_mod1 <- get_top_taxa(phy_funS_mod1_clean, 1, discard_other = T)
tax_table(top_funS_mod1)
abundances(top_funS_mod1) / sum(abundances(phy_funS_mod1_clean))
# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 2.
top_funS_mod2 <- get_top_taxa(phy_funS_mod2_clean, 1, discard_other = T)
tax_table(top_funS_mod2)
abundances(top_funS_mod2) / sum(abundances(phy_funS_mod2_clean))
# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 3.
top_funS_mod3 <- get_top_taxa(phy_funS_mod3_clean, 1, discard_other = T)
tax_table(top_funS_mod3)
abundances(top_funS_mod3) / sum(abundances(phy_funS_mod3_clean))
# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 4.
#top_funS_mod4 <- get_top_taxa(phy_funS_mod4_clean, 1, discard_other = T)
#tax_table(top_funS_mod4)
#abundances(top_funS_mod4) / sum(abundances(phy_funS_mod4_clean))

phy_modules_Schorf <- phyloseq(OTU_mod_fun_S,TAX,sampledataSchorf)

#####Alb

## check what modules are based on
# Read in the data on the modules obtained from Gephi. 
modules_fun_Alb <- read.csv2(("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/nw_modules_alb.csv"))
# Transform ASV information to dataframe to make handling easier.
otu_tab_fun_A <- as.data.frame(otu_table(phy_abundfilt_A))
# Calculate the total read count of each ASV.
otu_tab_fun_A$total <- rowSums(otu_tab_fun_A)
# Keep only the total read count.
otu_tab_fun_A <- otu_tab_fun_A %>%
  select(.,total) %>%
  rownames_to_column('Label')
# Merge module table and read count information.
otu_modules_A <- inner_join(otu_tab_fun_A, modules_fun_Alb)
# Initiate an empty dataframe were we can store the total read counts per module.
otu_mod_fun_A <- data.frame(Label = otu_modules_A$Label,
                            Mod0 = rep(0, 968),
                            Mod1 = rep(0, 968),
                            Mod2 = rep(0, 968),
                            Mod3 = rep(0, 968))
#Mod4 = rep(0, 968)

# If the ASV belongs to a module then write its read count, if not write zero.
otu_mod_fun_A$Mod0 <- if_else(otu_modules_A$modularity_class == 0, otu_modules_A$total, 0)
otu_mod_fun_A$Mod1 <- if_else(otu_modules_A$modularity_class == 1, otu_modules_A$total, 0)
otu_mod_fun_A$Mod2 <- if_else(otu_modules_A$modularity_class == 2, otu_modules_A$total, 0)
otu_mod_fun_A$Mod3 <- if_else(otu_modules_A$modularity_class == 3, otu_modules_A$total, 0)
#otu_mod_fun_A$Mod4 <- if_else(otu_modules_A$modularity_class == 4, otu_modules_A$total, 0)


# set the ASV name back to the rownames.
rownames(otu_mod_fun_A) <- otu_mod_fun_A$Label
otu_mod_fun_A$Label <- NULL
# Create the phyloseq object to easily obtain the taxonomic information.
otu_mat_mod_fun_A <- as.matrix(otu_mod_fun_A)
otu_mat_mod_fun_A
tax_mat_mod_fun_A <- as.matrix(tax_table(phy_abundfilt_A))
OTU_mod_fun_A <- otu_table(otu_mat_mod_fun_A, taxa_are_rows = T)
OTU_mod_fun_A
TAX_mod_fun_A <- tax_table(tax_mat_mod_fun_A)
phy_modules_fun_A <- phyloseq(OTU_mod_fun_A, TAX_mod_fun_A)
# Find out which modules contain more than 10 ASVs.
module_asv_count_fun_A <- colSums(otu_mat_mod_fun_A != 0)
which(module_asv_count_fun_A > 10)
sort(module_asv_count_fun_A)

####
# Find out which taxa are the top ones in the modules with more than ten ASVs.
####
# Subset by Module name.
phy_funA_mod0 <- prune_samples(sample_names(phy_modules_fun_A) == 'Mod0', phy_modules_fun_A)
phy_funA_mod1 <- prune_samples(sample_names(phy_modules_fun_A) == 'Mod1', phy_modules_fun_A)
phy_funA_mod2 <- prune_samples(sample_names(phy_modules_fun_A) == 'Mod2', phy_modules_fun_A)
phy_funA_mod3 <- prune_samples(sample_names(phy_modules_fun_A) == 'Mod3', phy_modules_fun_A)
#phy_funA_mod4 <- prune_samples(sample_names(phy_modules_fun_A) == 'Mod4', phy_modules_fun_A)

# Remove any taxa without reads.
phy_funA_mod0_clean <- prune_taxa(taxa_sums(phy_funA_mod0) > 0, phy_funA_mod0)
phy_funA_mod1_clean <- prune_taxa(taxa_sums(phy_funA_mod1) > 0, phy_funA_mod1)
phy_funA_mod2_clean <- prune_taxa(taxa_sums(phy_funA_mod2) > 0, phy_funA_mod2)
phy_funA_mod3_clean <- prune_taxa(taxa_sums(phy_funA_mod3) > 0, phy_funA_mod3)
#phy_funA_mod4_clean <- prune_taxa(taxa_sums(phy_funA_mod4) > 0, phy_funA_mod4)

# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 0.
top_funA_mod0 <- get_top_taxa(phy_funA_mod0_clean, 1, discard_other = T)
tax_table(top_funA_mod0)
abundances(top_funA_mod0) / sum(abundances(phy_funA_mod0_clean))
# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 1.
top_funA_mod1 <- get_top_taxa(phy_funA_mod1_clean, 1, discard_other = T)
tax_table(top_funA_mod1)
abundances(top_funA_mod1) / sum(abundances(phy_funA_mod1_clean))
# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 2.
top_funA_mod2 <- get_top_taxa(phy_funA_mod2_clean, 1, discard_other = T)
tax_table(top_funA_mod2)
abundances(top_funA_mod2) / sum(abundances(phy_funA_mod2_clean))
# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 3.
top_funA_mod3 <- get_top_taxa(phy_funA_mod3_clean, 1, discard_other = T)
tax_table(top_funA_mod3)
abundances(top_funA_mod3) / sum(abundances(phy_funA_mod3_clean))
# Retrieve the top taxa, its taxonomy and relative abundance for each organismal group within Module 4.
#top_funA_mod4 <- get_top_taxa(phy_funA_mod4_clean, 1, discard_other = T)
#tax_table(top_funA_mod4)
#abundances(top_funA_mod4) / sum(abundances(phy_funA_mod4_clean))

############determine modules in networks for Alb and Schorfheide########

####Alb####

taxmat_mod_A <- read.csv2("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/tax_mod_analysis_Alb.csv")
taxmat_mod_A
taxmat_mod_A = as.matrix(taxmat_mod_A, nrow = nrow(taxmat_mod_A),ncol = 8) 
taxmat_mod_A
#taxmat_mod_A = as.data.frame(taxmat_mod_A, nrow = nrow(taxmat_mod_A), ncol = 9)
colnames(taxmat_mod_A) <- c("ASV","Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species","Module")
rownames(taxmat_mod_A) <- taxmat_mod_A[,1]
#rownames(taxmax_mod_A) <- taxmat_mod_A %>% remove_rownames %>% column_to_rownames(var="ASV") #geht nur bei data frame
taxmat_mod_A <- taxmat_mod_A[,-1]
TAX_mod_A = tax_table(taxmat_mod_A)
TAX_mod_A

physeqAlb_mod = phyloseq(ASV, TAX_mod_A, sampledataAlb)
ps_rel_abund_A_mod = phyloseq::transform_sample_counts(physeqAlb_mod, function(x){x / sum(x)})

#plot_bar(ps_rel_abund_A_mod, x="Module")+
 # facet_wrap(~ venn_class, scales = "free") 

#plot_bar(ps_rel_abund_A_mod, x="Module", fill = "Module")+
  #facet_wrap(~ venn_class, scales = "free") +
  #scale_colour_manual(values = c("coral1","green3","deepskyblue","magenta2"))

plot_bar(ps_rel_abund_A_mod, x="Module", fill = "Module")+
  facet_wrap(~ venn_class, scales = "free") +
  theme(strip.text.x = element_text(size = 20, face = "italic"))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 20)) +
  theme(legend.position = "bottom",legend.direction = "horizontal") +
  theme(axis.title.y = element_text(size = 20))+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y = element_text(size = 16))+
  scale_colour_manual(values = c("coral1","green3","deepskyblue","magenta2"))+
  geom_bar(stat="identity")+
  labs(x = "", y = "proportion of total ASVs in %") +
  ylim(0, 30)
ggsave(
  "modules_Alb.png",
  width = 11.7,
  height = 8.3,
  dpi = 800
)

#20220728 use 4 module network and adjust colors --> no new names assigned to objects!
taxmat_mod_A <- read.csv2("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/tax_mod_analysis_Alb_2.csv")
taxmat_mod_A
taxmat_mod_A = as.matrix(taxmat_mod_A, nrow = nrow(taxmat_mod_A),ncol = 8) 
taxmat_mod_A
#taxmat_mod_A = as.data.frame(taxmat_mod_A, nrow = nrow(taxmat_mod_A), ncol = 9)
colnames(taxmat_mod_A) <- c("ASV","Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species","Module")
rownames(taxmat_mod_A) <- taxmat_mod_A[,1]
#rownames(taxmax_mod_A) <- taxmat_mod_A %>% remove_rownames %>% column_to_rownames(var="ASV") #geht nur bei data frame
taxmat_mod_A <- taxmat_mod_A[,-1]
TAX_mod_A = tax_table(taxmat_mod_A)
TAX_mod_A

physeqAlb_mod = phyloseq(ASV, TAX_mod_A, sampledataAlb)
ps_rel_abund_A_mod = phyloseq::transform_sample_counts(physeqAlb_mod, function(x){x / sum(x)})

plot_bar(ps_rel_abund_A_mod, x="Module", fill = "Module")+
  facet_wrap(~ venn_class, scales = "free") +
  theme(strip.text.x = element_text(size = 20, face = "italic"))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 20)) +
  theme(legend.position = "bottom",legend.direction = "horizontal") +
  theme(axis.title.y = element_text(size = 20))+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y = element_text(size = 16))+
  scale_fill_manual(values = c("coral2","green3","mediumorchid","steelblue1"))+
  geom_bar(stat="identity")+
  labs(x = "", y = "proportion of reads in %") +
  ylim(0, 30)

####Schorfheide####

#write.csv(tax_mat_mod_fun_S,"C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/tax_mat_S.csv")

taxmat_mod_S <- read.csv2(("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/tax_mod_analysis_Schorf.csv"))
taxmat_mod_S
taxmat_mod_S = as.matrix(taxmat_mod_S, nrow = nrow(taxmat_mod_S),ncol = 8) 
taxmat_mod_S
colnames(taxmat_mod_S) <- c("ASV","Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species","Module")
rownames(taxmat_mod_S) <- taxmat_mod_S[,1]
taxmat_mod_S <- taxmat_mod_S[,-1]
taxmat_mod_S
TAX_mod_S = tax_table(taxmat_mod_S)
TAX_mod_S

physeqSchorf_mod = phyloseq(ASV, TAX_mod_S, sampledataSchorf)
ps_rel_abund_S_mod = phyloseq::transform_sample_counts(physeqSchorf_mod, function(x){x / sum(x)})

plot_bar(ps_rel_abund_S_mod, x="Module", fill = "Module")+
  facet_wrap(~ venn_class, scales = "free") +
  theme(strip.text.x = element_text(size = 20, face = "italic"))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 20)) +
  theme(legend.position = "bottom",legend.direction = "horizontal") +
  theme(axis.title.y = element_text(size = 20))+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y = element_text(size = 16))+
  scale_colour_manual(values = c("coral1","green3","deepskyblue","magenta2"))+
  geom_bar(stat="identity")+
  labs(x = "", y = "proportion of total ASVs in %") +
  ylim(0, 20)
ggsave(
  "modules_Schorf.png",
  width = 11.7,
  height = 8.3,
  dpi = 800
)

#20220728 use 4 module network and adjust colors --> no new names assigned to objects!
taxmat_mod_S <- read.csv2(("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/tax_mod_analysis_Schorf_2.csv"))
taxmat_mod_S
taxmat_mod_S = as.matrix(taxmat_mod_S, nrow = nrow(taxmat_mod_S),ncol = 8) 
taxmat_mod_S
colnames(taxmat_mod_S) <- c("ASV","Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species","Module")
rownames(taxmat_mod_S) <- taxmat_mod_S[,1]
taxmat_mod_S <- taxmat_mod_S[,-1]
taxmat_mod_S
TAX_mod_S = tax_table(taxmat_mod_S)
TAX_mod_S

physeqSchorf_mod = phyloseq(ASV, TAX_mod_S, sampledataSchorf)
ps_rel_abund_S_mod = phyloseq::transform_sample_counts(physeqSchorf_mod, function(x){x / sum(x)})

plot_bar(ps_rel_abund_S_mod, x="Module", fill = "Module")+
  facet_wrap(~ venn_class, scales = "free") +
  theme(strip.text.x = element_text(size = 20, face = "italic"))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 20)) +
  theme(legend.position = "bottom",legend.direction = "horizontal") +
  theme(axis.title.y = element_text(size = 20))+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y = element_text(size = 16))+
  scale_fill_manual(values = c("coral2","green3","mediumorchid","steelblue1"))+
  geom_bar(stat="identity")+
  labs(x = "", y = "proportion of reads in %") +
  ylim(0, 18)



#########
############networks for trees Alb and Schorfheide########
###########

###Alb#####

###Alb Fagus_sylvatica

#subset dataset to ASVs that occur in at least one percent of samples
phy_trans_A_Fagus  = transform_sample_counts(physeqAlb_F, function(x) x / sum(x) )
phy_abundfilt_A_Fagus = filter_taxa(phy_trans_A_Fagus, function(x) sum(x) > .01, TRUE)
keep_names_A_Fagus <- taxa_names(phy_abundfilt_A_Fagus)
my_subset_A_Fagus <- subset(otu_table(physeqAlb_F), rownames(otu_table(physeqAlb_F)) %in% keep_names_A_Fagus)
phy_raw_abundfilt_A_Fagus <- merge_phyloseq(my_subset_A_Fagus, tax_table(physeqAlb_F), sample_data(physeqAlb_F))

pargs2 <- list(rep.num=50, seed=10010)

#SpiecEasi algorithm execution
se_raw_abundfilt_A_Fagus <- spiec.easi(phy_raw_abundfilt_A_Fagus, method='mb', 
                                       lambda.min.ratio=1e-2, nlambda=70, 
                                       sel.criterion='bstars',  pulsar.select=TRUE,
                                       pulsar.params=pargs2)

#summary
se_raw_abundfilt_A_Fagus$select$stars$summary

#get optimal parameters
getOptInd(se_raw_abundfilt_A_Fagus)
sum(getRefit(se_raw_abundfilt_A_Fagus))/2

#network stability -> optimum close to 0.05
getStability(se_raw_abundfilt_A_Fagus)

fun.mb_A_Fagus <- adj2igraph(getRefit(se_raw_abundfilt_A_Fagus),  
                             vertex.attr=list(name=taxa_names(phy_raw_abundfilt_A_Fagus)))

##plot network
# convert to gephi 
fun.mb_A_Fagus_gephi <- igraph.to.gexf(fun.mb_A_Fagus)
write.gexf(fun.mb_A_Fagus_gephi, output = here("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/fungi_A_Fagus.gexf"))


##get hub taxa, through Kleinbergs centrality
fun_between_A_Fagus <- betweenness(fun.mb_A_Fagus, directed = F)
fun_top5_between_A_Fagus <- sort(fun_between_A_Fagus, decreasing = T)[1:5]
names(fun_top5_between_A_Fagus)
hub_taxa_fun_A_Fagus <- subset_taxa(physeqAlb_F, taxa_names(physeqAlb_F) %in% names(fun_top5_between_A_Fagus))
tax_table(hub_taxa_fun_A_Fagus)
hubtax_AlbF <- data.frame(tax_table(hub_taxa_fun_A_Fagus))
hubtax_AlbF <- rownames_to_column(hubtax_AlbF, var = "ASV_ID") %>% as_tibble()
hubtax_AlbF

###Alb Picea_abies

#subset dataset to ASVs that occur in at least one percent of samples
phy_trans_A_Picea  = transform_sample_counts(physeqAlb_P, function(x) x / sum(x) )
phy_abundfilt_A_Picea = filter_taxa(phy_trans_A_Picea, function(x) sum(x) > .01, TRUE)
keep_names_A_Picea <- taxa_names(phy_abundfilt_A_Picea)
my_subset_A_Picea <- subset(otu_table(physeqAlb_P), rownames(otu_table(physeqAlb_P)) %in% keep_names_A_Picea)
phy_raw_abundfilt_A_Picea <- merge_phyloseq(my_subset_A_Picea, tax_table(physeqAlb_P), sample_data(physeqAlb_P))

pargs2 <- list(rep.num=50, seed=10010)

#SpiecEasi algorithm execution
se_raw_abundfilt_A_Picea <- spiec.easi(phy_raw_abundfilt_A_Picea, method='mb', 
                                       lambda.min.ratio=1e-2, nlambda=70, 
                                       sel.criterion='bstars',  pulsar.select=TRUE,
                                       pulsar.params=pargs2)

#summary
se_raw_abundfilt_A_Picea$select$stars$summary

#get optimal parameters
getOptInd(se_raw_abundfilt_A_Picea)
sum(getRefit(se_raw_abundfilt_A_Picea))/2

#network stability -> optimum close to 0.05
getStability(se_raw_abundfilt_A_Picea)

fun.mb_A_Picea <- adj2igraph(getRefit(se_raw_abundfilt_A_Picea),  
                             vertex.attr=list(name=taxa_names(phy_raw_abundfilt_A_Picea)))

##plot network
# convert to gephi 
fun.mb_A_Picea_gephi <- igraph.to.gexf(fun.mb_A_Picea)
write.gexf(fun.mb_A_Picea_gephi, output = here("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/fungi_A_Picea.gexf"))


##get hub taxa, through Kleinbergs centrality
fun_between_A_Picea <- betweenness(fun.mb_A_Picea, directed = F)
fun_top5_between_A_Picea <- sort(fun_between_A_Picea, decreasing = T)[1:5]
names(fun_top5_between_A_Picea)
hub_taxa_fun_A_Picea <- subset_taxa(physeqAlb_P, taxa_names(physeqAlb_P) %in% names(fun_top5_between_A_Picea))
tax_table(hub_taxa_fun_A_Picea)


###Schorfheide Fagus_sylvatica

#subset dataset to ASVs that occur in at least one percent of samples
phy_trans_S_Fagus  = transform_sample_counts(physeqSchorf_F, function(x) x / sum(x) )
phy_abundfilt_S_Fagus = filter_taxa(phy_trans_S_Fagus, function(x) sum(x) > .01, TRUE)
keep_names_S_Fagus <- taxa_names(phy_abundfilt_S_Fagus)
my_subset_S_Fagus <- subset(otu_table(physeqSchorf_F), rownames(otu_table(physeqSchorf_F)) %in% keep_names_S_Fagus)
phy_raw_abundfilt_S_Fagus <- merge_phyloseq(my_subset_S_Fagus, tax_table(physeqSchorf_F), sample_data(physeqSchorf_F))

pargs2 <- list(rep.num=50, seed=10010)

#SpiecEasi algorithm execution
se_raw_abundfilt_S_Fagus <- spiec.easi(phy_raw_abundfilt_S_Fagus, method='mb', 
                                       lambda.min.ratio=1e-2, nlambda=70, 
                                       sel.criterion='bstars',  pulsar.select=TRUE,
                                       pulsar.params=pargs2)

#summary
se_raw_abundfilt_S_Fagus$select$stars$summary

#get optimal parameters
getOptInd(se_raw_abundfilt_S_Fagus)
sum(getRefit(se_raw_abundfilt_S_Fagus))/2

#network stability -> optimum close to 0.05
getStability(se_raw_abundfilt_S_Fagus)

fun.mb_S_Fagus <- adj2igraph(getRefit(se_raw_abundfilt_S_Fagus),  
                             vertex.attr=list(name=taxa_names(phy_raw_abundfilt_S_Fagus)))

##plot network
# convert to gephi 
fun.mb_S_Fagus_gephi <- igraph.to.gexf(fun.mb_S_Fagus)
write.gexf(fun.mb_S_Fagus_gephi, output = here("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/fungi_S_Fagus.gexf"))


##get hub taxa, through Kleinbergs centrality
fun_between_S_Fagus <- betweenness(fun.mb_S_Fagus, directed = F)
fun_top5_between_S_Fagus <- sort(fun_between_S_Fagus, decreasing = T)[1:5]
names(fun_top5_between_S_Fagus)
hub_taxa_fun_S_Fagus <- subset_taxa(physeqSchorf_F, taxa_names(physeqSchorf_F) %in% names(fun_top5_between_S_Fagus))
tax_table(hub_taxa_fun_S_Fagus)


###Schorfheide Pinus_sylvestris

#subset dataset to ASVs that occur in at least one percent of samples
phy_trans_S_Pinus  = transform_sample_counts(physeqSchorf_P, function(x) x / sum(x) )
phy_abundfilt_S_Pinus = filter_taxa(phy_trans_S_Pinus, function(x) sum(x) > .01, TRUE)
keep_names_S_Pinus <- taxa_names(phy_abundfilt_S_Pinus)
my_subset_S_Pinus <- subset(otu_table(physeqSchorf_P), rownames(otu_table(physeqSchorf_P)) %in% keep_names_S_Pinus)
phy_raw_abundfilt_S_Pinus <- merge_phyloseq(my_subset_S_Pinus, tax_table(physeqSchorf_P), sample_data(physeqSchorf_P))

pargs2 <- list(rep.num=50, seed=10010)

#SpiecEasi algorithm execution
se_raw_abundfilt_S_Pinus <- spiec.easi(phy_raw_abundfilt_S_Pinus, method='mb', 
                                       lambda.min.ratio=1e-2, nlambda=70, 
                                       sel.criterion='bstars',  pulsar.select=TRUE,
                                       pulsar.params=pargs2)

#summary
se_raw_abundfilt_S_Pinus$select$stars$summary

#get optimal parameters
getOptInd(se_raw_abundfilt_S_Pinus)
sum(getRefit(se_raw_abundfilt_S_Pinus))/2

#network stability -> optimum close to 0.05
getStability(se_raw_abundfilt_S_Pinus)

fun.mb_S_Pinus <- adj2igraph(getRefit(se_raw_abundfilt_S_Pinus),  
                             vertex.attr=list(name=taxa_names(phy_raw_abundfilt_S_Pinus)))

##plot network
# convert to gephi 
fun.mb_S_Pinus_gephi <- igraph.to.gexf(fun.mb_S_Pinus)
write.gexf(fun.mb_S_Pinus_gephi, output = here("C:/Users/behof/Desktop/MSc Umweltwissenschaften/Master Thesis/data_analysis/fungi_S_Pinus.gexf"))


##get hub taxa, through Kleinbergs centrality
fun_between_S_Pinus <- betweenness(fun.mb_S_Pinus, directed = F)
fun_top5_between_S_Pinus <- sort(fun_between_S_Pinus, decreasing = T)[1:5]
names(fun_top5_between_S_Pinus)
hub_taxa_fun_S_Pinus <- subset_taxa(physeqSchorf_P, taxa_names(physeqSchorf_P) %in% names(fun_top5_between_S_Pinus))
tax_table(hub_taxa_fun_S_Pinus)

#20220705 approach for hub taxa from influential package
# Reconstructing the graph
#set.seed(70)
#My_graph <-  igraph::fun.mb_S(n = 50, m = 120, directed = TRUE)
# did not work

# Calculating the IVI values
My_graph_IVI <- ivi(fun.mb_S, directed = FALSE)

# Visualizing the graph based on IVI values
My_graph_IVI_Vis <- cent_network.vis(graph = fun.mb_S,
                                     cent.metric = My_graph_IVI,
                                     directed = FALSE,
                                     plot.title = "IVI-based Network",
                                     legend.title = "IVI value")

My_graph_IVI_Vis

#SIRIR model
# Reconstructing the graph
#My_graph <-  sif2igraph(Path = "Sample_SIF.sif")       

# Extracting the vertices
GraphVertices <- V(fun.mb_S)        

# Calculation of influence rank
Influence.Ranks <- sirir(graph = fun.mb_S,     
                         vertices = GraphVertices, 
                         beta = 0.5, gamma = 1, no.sim = 10, seed = 1234)
Influence.Ranks
write.csv(Influence.Ranks, "C:/Users/behof/Desktop/MSc Umweltwissenschaften/Paper_Masterthesis/hubtax_S_sirir.csv")

#comparison ranks Kleinberg centrality
##get hub taxa, through Kleinbergs centrality
#fun_between_S <- betweenness(fun.mb_S, directed = F)
fun_top596_between_S <- sort(fun_between_S, decreasing = T)[1:596]
write.csv(fun_top596_between_S,"C:/Users/behof/Desktop/MSc Umweltwissenschaften/Paper_Masterthesis/fun_top596_between_S.csv")

