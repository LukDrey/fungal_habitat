#################################################################
##                          Section 1                          ##
##                       Package Loading                       ##
#################################################################

library(phyloseq)

library(ggplot2)

library(ggvenn)

library(MicEco)

library(tidyverse)

library(ape)

library(rstatix)

library(ggpubr)

library(ggdist)

library(vegan)

library(paletteer)

library(ranacapa)

library(devtools)

library(fantaxtic)

library(SpiecEasi)

library(rgexf)

library(igraph)

library(here)

library(xlsx)

library(dplyr)

library(rlang)

library(metamicrobiomeR) 

library(remotes)

library(microbiome)

library(gghalves)

library(viridis)

library(pals)

library(influential)

library(ade4)

##################################################################
##                          Section 2                           ##
##         Data Loading, Data Cleaning and Initialization       ##
##################################################################

##---------------------------------------------------------------
##                        Taxonomy Table                        -
##---------------------------------------------------------------

# Load in taxonomy table. 
tax_fungi <- base::readRDS(here::here("fungi_tax_euk.rds")) %>% 
  base::as.data.frame() %>%
  tibble::rownames_to_column('sequence') %>%
  dplyr::rename(sequence_fungi = sequence)

# Load the fungal reads.
fungi_seqs_fasta <- Biostrings::readDNAStringSet(here::here('ASVs_fungi.fa'))

# Make a dataframe of the sequences and their ASV ID. 
seq_name_fungi <- base::names(fungi_seqs_fasta)
sequence_fungi <- base::paste(fungi_seqs_fasta)
fungi_rep_seqs <- base::data.frame(seq_name_fungi, sequence_fungi)

# Join the taxonomy table and the representative sequences
tax_clean_fungi <- dplyr::left_join(tax_fungi, fungi_rep_seqs, by = 'sequence_fungi')

# Split the taxonomy into different columns of taxonomic levels.
fungi_tax_fin <- tidyr::separate(tax_clean_fungi, Kingdom, c(NA, 'Kingdom') , sep = '__') %>% 
  tidyr::separate(Phylum, c(NA, 'Phylum') , sep = '__') %>% 
  tidyr::separate(Class, c(NA, 'Class') , sep = '__') %>% 
  tidyr::separate(Order, c(NA, 'Order') , sep = '__') %>% 
  tidyr::separate(Family, c(NA, 'Family') , sep = '__') %>% 
  tidyr::separate(Genus, c(NA, 'Genus') , sep = '__') %>% 
  tidyr::separate(Species, c(NA, 'Species') , sep = '__')

# Keep only Fungi 
fungi_tax_fin <- fungi_tax_fin %>% 
  dplyr::filter(Kingdom == "Fungi")

# Rename the ASV_ID column. 
fungi_tax_fin <- dplyr::rename(fungi_tax_fin, ASV_ID = seq_name_fungi)


# Set rownames.
base::row.names(fungi_tax_fin) <- fungi_tax_fin$ASV_ID

fungi_tax_fin$sequence_fungi <- NULL
fungi_tax_fin$ASV_ID <- NULL

##----------------------------------------------------------------
##                           Metadata                            -
##----------------------------------------------------------------

# Load in Metadata. 
metadata <- utils::read.csv(here::here('sample_data_3_AHS.csv'), sep = ';')

# Set the sample name as the rowname for the phyloseq creation.
base::row.names(metadata) <- metadata$sample
metadata$sample <- NULL

##----------------------------------------------------------------
##                          ASV table                            -
##----------------------------------------------------------------

# Load in ASV table previously curated using the LULU algorithm (https://doi.org/10.1038/s41467-017-01312-x).
# For code on how it was employed see https://github.com/LukDrey/beech_micro_communities. 
ASV_table_fungi_cur <- base::readRDS(here::here("ASV_table_fungi_cur.rds"))

# Keep only samples that do represent real tree swabs. Cut Controls. 
asv_fungi <- ASV_table_fungi_cur$curated_table %>% 
  dplyr::select(all_of(base::rownames(metadata)))

##---------------------------------------------------------------
##               Create the Phyloseq Objects                     -
##---------------------------------------------------------------

# First create a phyloseq object with all samples included.
# Create the matrices needed for the phyloseq functions.  
asv_mat <- base::as.matrix(asv_fungi)
taxmat <- base::as.matrix(fungi_tax_fin) 
  
# Create the ps object. 
ASV <- phyloseq::otu_table(asv_mat, taxa_are_rows = TRUE)
TAX <- phyloseq::tax_table(taxmat)
sampledata <- phyloseq::sample_data(metadata)

full_physeq <- phyloseq::phyloseq(ASV, TAX, sampledata)

# Filter out samples that were sampled on trees that are not from our three target species. 

filtered_physeq <- phyloseq::subset_samples(full_physeq, 
                                            dominant_tree %in% c("Fagus_sylvatica",
                                                                 "Pinus_sylvestris", 
                                                                 "Picea_abies")) %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

# Split the phyloseq object into the two exploratories.

physeq_alb <- phyloseq::subset_samples(filtered_physeq, exploratory == "Alb") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_sch <- phyloseq::subset_samples(filtered_physeq, exploratory == "Schorfheide") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

# Split the phyloseq object into bark and soil samples. 
physeq_bark <- phyloseq::subset_samples(filtered_physeq, substrate == "bark") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_soil <- phyloseq::subset_samples(filtered_physeq, substrate == "soil") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

# Split the phyloseq object into the two exploratories & the substrate.
# Swabian Alb 
physeq_alb_bark <- phyloseq::subset_samples(physeq_alb, substrate == "bark") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_alb_soil <- phyloseq::subset_samples(physeq_alb, substrate == "soil") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

# Schorfheide-Chorin
physeq_sch_bark <- phyloseq::subset_samples(physeq_sch, substrate == "bark") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_sch_soil <- phyloseq::subset_samples(physeq_sch, substrate == "soil") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

# Split the phyloseqobject by exploratories & the tree species. 
# Swabian Alb
physeq_alb_fagus <- phyloseq::subset_samples(physeq_alb, dominant_tree == "Fagus_sylvatica") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_alb_picea <- phyloseq::subset_samples(physeq_alb, dominant_tree == "Picea_abies") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

# Schorfheide Chorin
physeq_sch_fagus <- phyloseq::subset_samples(physeq_sch, dominant_tree == "Fagus_sylvatica") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_sch_pinus <- phyloseq::subset_samples(physeq_sch, dominant_tree == "Pinus_sylvestris") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

# Split the phyloseq object into the two exploratories & the substrate & the tree species.
# Swabian Alb
physeq_alb_bark_fagus <- phyloseq::subset_samples(physeq_alb_bark, dominant_tree == "Fagus_sylvatica") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_alb_soil_fagus <- phyloseq::subset_samples(physeq_alb_soil, dominant_tree == "Fagus_sylvatica") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_alb_bark_picea <- phyloseq::subset_samples(physeq_alb_bark, dominant_tree == "Picea_abies") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_alb_soil_picea <- phyloseq::subset_samples(physeq_alb_soil, dominant_tree == "Picea_abies") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

# Schorfheide-Chorin
physeq_sch_bark_fagus <- phyloseq::subset_samples(physeq_sch_bark, dominant_tree == "Fagus_sylvatica") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_sch_soil_fagus <- phyloseq::subset_samples(physeq_sch_soil, dominant_tree == "Fagus_sylvatica") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_sch_bark_pinus <- phyloseq::subset_samples(physeq_sch_bark, dominant_tree == "Pinus_sylvestris") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_sch_soil_pinus <- phyloseq::subset_samples(physeq_sch_soil, dominant_tree == "Pinus_sylvestris") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)


#################################################################
##                          Section 3                          ##
##                  Analyse the library sizes                  ##
#################################################################

################Swabian Alb###########################

# Create a sample_data column containing a column that corresponds to tree species and substrate.
physeq_alb_curve <- physeq_alb
phyloseq::sample_data(physeq_alb_curve) <- phyloseq::sample_data(physeq_alb) %>% 
  base::data.frame() %>%  
  dplyr::mutate(tree_substrate = base::paste(dominant_tree, substrate, sep = "-"))

# Create rarefaction curve and color the lines by substrate (bark/soil) and host tree species. 
rare_alb <- ranacapa::ggrare(physeq_alb_curve, step = 50, color = "tree_substrate", se = FALSE)

################Schorfheide###########################

# Create a sample_data column containing a column that corresponds to tree species and substrate.
physeq_sch_curve <- physeq_sch
phyloseq::sample_data(physeq_sch_curve) <- phyloseq::sample_data(physeq_sch) %>% 
  base::data.frame() %>%  
  dplyr::mutate(tree_substrate = base::paste(dominant_tree, substrate, sep = "-"))

# Create rarefaction curve and color the lines by substrate (bark/soil) and host tree species. 
rare_sch <- ranacapa::ggrare(physeq_sch_curve, step = 50, color = "tree_substrate", se = FALSE)


#################################################################
##                          Section 4                          ##
##                   Alpha diversity analyses                  ##
#################################################################

################Swabian Alb###########################

# Calculate the Shannon Diversity.
shannon_alb <- phyloseq::estimate_richness(physeq_alb, split = T, measures = 'Shannon') %>% 
  tibble::rownames_to_column()

# Extract the sample data to append it.
sample_data_alb <- base::data.frame(phyloseq::sample_data(physeq_alb)) %>% 
  tibble::rownames_to_column()

# Append sample_data to get the information on tree and substrate. 
div_data_alb <- dplyr::left_join(shannon_alb, sample_data_alb)

# Wilcoxon Rank Sum Test to test differences in shannon diversity between substrate and tree groups.
  
# bark vs soil
bark_alb = div_data_alb[div_data_alb$substrate == "bark",]
soil_alb = div_data_alb[div_data_alb$substrate == "soil",]
stats::wilcox.test(bark_alb$Shannon, soil_alb$Shannon)

# Differences between tree species within the substrate.
# bark Fagus vs. Picea
bark_fagus_alb = bark_alb[bark_alb$dominant_tree == "Fagus_sylvatica",]
bark_picea_alb = bark_alb[bark_alb$dominant_tree == "Picea_abies",]
stats::wilcox.test(bark_fagus_alb$Shannon, bark_picea_alb$Shannon)

# soil Fagus vs. Picea
soil_fagus_alb = soil_alb[soil_alb$dominant_tree == "Fagus_sylvatica",]
soil_picea_alb = soil_alb[soil_alb$dominant_tree == "Picea_abies",]
stats::wilcox.test(soil_fagus_alb$Shannon, soil_picea_alb$Shannon)

# Full Fagus vs. Picea Swabian Alb not subset by substrate
fagus_alb = div_data_alb[div_data_alb$dominant_tree == "Fagus_sylvatica",]
picea_alb = div_data_alb[div_data_alb$dominant_tree == "Picea_abies",]
stats::wilcox.test(fagus_alb$Shannon, picea_alb$Shannon)

################Schorfheide-Chorin###########################

# Calculate the Shannon Diversity.
shannon_sch <- phyloseq::estimate_richness(physeq_sch, split = T, measures = 'Shannon') %>% 
  tibble::rownames_to_column()

# Extract the sample data to append it.
sample_data_sch <- base::data.frame(phyloseq::sample_data(physeq_sch)) %>% 
  tibble::rownames_to_column()

# Append sample_data to get the information on tree and substrate. 
div_data_sch <- dplyr::left_join(shannon_sch, sample_data_sch)

# Wilcoxon Rank Sum Test to test differences in shannon diversity between substrate and tree groups.

# bark vs soil
bark_sch = div_data_sch[div_data_sch$substrate == "bark",]
soil_sch = div_data_sch[div_data_sch$substrate == "soil",]
stats::wilcox.test(bark_sch$Shannon, soil_sch$Shannon)

# Differences between tree species within the substrate.
# bark Fagus vs. Picea
bark_fagus_sch = bark_sch[bark_sch$dominant_tree == "Fagus_sylvatica",]
bark_pinus_sch = bark_sch[bark_sch$dominant_tree == "Pinus_sylvestris",]
stats::wilcox.test(bark_fagus_sch$Shannon, bark_pinus_sch$Shannon)

# soil Fagus vs. Picea
soil_fagus_sch = soil_sch[soil_sch$dominant_tree == "Fagus_sylvatica",]
soil_pinus_sch = soil_sch[soil_sch$dominant_tree == "Pinus_sylvestris",]
stats::wilcox.test(soil_fagus_sch$Shannon, soil_pinus_sch$Shannon)

# Full Fagus vs. Picea Swabian sch not subset by substrate
fagus_sch = div_data_sch[div_data_sch$dominant_tree == "Fagus_sylvatica",]
pinus_sch = div_data_sch[div_data_sch$dominant_tree == "Pinus_sylvestris",]
stats::wilcox.test(fagus_sch$Shannon, pinus_sch$Shannon)

#################################################################
##                          Section 4                          ##
##                    Beta diversity analyses                  ##
#################################################################

#########################NMDS-Ordination############################
# Ordinate the Swabian Alb phyloseq using an NMDS with Bray-Curtis distance.
nmds_alb <- phyloseq::ordinate(physeq_alb, method = "NMDS", distance = "bray")

# Plot the ordination. 
ordination_alb <- phyloseq::plot_ordination(physeq_alb, nmds_alb, type="samples", color="dominant_tree", shape="substrate") + 
  ggplot2::geom_point(size = 4) +
  ggplot2::stat_ellipse(ggplot2::aes(group = dominant_tree), linetype = 2) +
  ggplot2::scale_colour_manual(values = c("green","darkgreen"), name = "dominant tree species",
                               labels = c("Fagus sylvatica", "Picea abies")) +
  ggplot2::labs(subtitle = "(A) Swabian Alb") +
  ggplot2::theme(legend.position = "none",
        axis.text = ggplot2::element_text(size = 15),
        axis.title = ggplot2::element_text(size = 15))

ordination_alb

# Ordinate the Schorfheide-Chorin phyloseq using an NMDS with Bray-Curtis distance.
nmds_sch <- phyloseq::ordinate(physeq_sch, method = "NMDS", distance = "bray")

# Plot the ordination.
ordination_sch <- phyloseq::plot_ordination(physeq_sch, nmds_sch, type="samples", color="dominant_tree", shape="substrate") + 
  ggplot2::geom_point(size = 4) +
  ggplot2::stat_ellipse(ggplot2::aes(group = dominant_tree), linetype = 2) +
  ggplot2::scale_colour_manual(values = c("green","darkolivegreen4"), name = "dominant tree species",
                               labels = c("Fagus sylvatica", "Pinus sylvestris")) +
  ggplot2::labs(subtitle = "(B) Schorfheide-Chorin") +
  ggplot2::theme(legend.position = "none",
        axis.text = ggplot2::element_text(size = 15), 
        axis.title = ggplot2::element_text(size = 15))

ordination_sch

# Run the ordination on the full dataset to grab a nice legend. 
nmds_full <- phyloseq::ordinate(filtered_physeq, method = "NMDS", distance = "bray")

ordination_full <- phyloseq::plot_ordination(filtered_physeq, nmds_full, type="samples", color="dominant_tree", shape="substrate") + 
  ggplot2::geom_point(size = 4) +
  ggplot2::stat_ellipse(ggplot2::aes(group = dominant_tree), linetype = 2) +
  ggplot2::scale_colour_manual(values = c("green", "darkgreen", "darkolivegreen4"), name = "tree species",
                               labels = c("Fagus sylvatica","Picea abies", "Pinus sylvestris")) +
  ggplot2::theme(legend.text = ggplot2::element_text(face = "italic", size = 15),
        legend.title = ggplot2::element_text(size = 15, face = "bold"),
        legend.position = "right",
        legend.direction = "vertical",
        legend.spacing = ggplot2::unit(2,"cm"),
        axis.text = ggplot2::element_text(size = 15), 
        axis.title = ggplot2::element_text(size = 15))

ordination_full

# Extract the legend and store it as a ggplot object. 
ordination_legend <- ggpubr::get_legend(ordination_full)

# Combine the figures into one plot. 
ordination_final <- ggpubr::ggarrange(ordination_alb, ordination_sch,
                                      ncol = 1, nrow = 2,
                                      legend = "right", legend.grob = ordination_legend)

#################################################################
##                          Section 5                          ##
##                    Variance partitioning                    ##
#################################################################

# Function to extract an otu_table from a phyloseq object in the right format for vegan::varpart.
veganotu = function(physeq) {
  require("vegan")
  OTU = phyloseq::otu_table(physeq)
  if (phyloseq::taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}

################Swabian Alb###########################

# Extract OTU table Alb.
otu_alb <- veganotu(physeq_alb)

# Extract the sample data.
data_alb <- base::data.frame(phyloseq::sample_data(physeq_alb))

# Get the variation with vegan::varpart.
varp_alb <- vegan::varpart(otu_alb, ~ substrate, ~ dominant_tree, data = data_alb)

indfract_adj_r_alb <- varp_alb$part$indfract$Adj.R.squared
indfract_adj_r_alb

fract_adj_r_alb <- varp_alb$part$fract$Adj.R.squared
fract_adj_r_alb

### Schorfheide ###

# Extract OTU table Schorfheide.
otu_sch <- veganotu(physeq_sch)

# Extract sample data.
data_sch <- base::data.frame(phyloseq::sample_data(physeq_sch))

# Get the variation with vegan::varpart.
varp_sch <- vegan::varpart(otu_sch, ~ substrate, ~ dominant_tree, data = data_sch)

indfract_adj_r_sch <- varp_sch$part$indfract$Adj.R.squared
indfract_adj_r_sch

fract_adj_r_sch <- varp_sch$part$fract$Adj.R.squared
fract_adj_r_sch

div_labels <- c("\u03B1-Diversity", "\u03B2-Diversity")

ggpubr::ggbarplot(variance_full, x = "q_lev", y = "variance", 
                  fill = "variable", color = "variable", palette = alpha(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                                                                           "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                                                                         alpha = 0.9)) +
  ggplot2::facet_grid(organism ~ div_lev, space="free", scales="free", switch = "y") +
  ggplot2::scale_x_discrete(position = "top") +
  ggplot2::geom_hline(aes(yintercept = 0), linetype = "dashed") +
  ggplot2::coord_flip() + 
  ggplot2::ylab("Explained Variance") + 
  ggplot2::theme(text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.position = "right",
                 axis.title.y = element_blank(),
                 axis.line.y = element_blank(), 
                 axis.ticks.y = element_blank())


#################################################################
##                          Section 6                          ##
##                        Venn diagrams                        ##
#################################################################

#####################Alb#######################
###
# Venn diagrams by host tree species. 
###

# Check the absolute number of tree per category.
MicEco::ps_venn(physeq_alb_soil, group = "dominant_tree", fraction = 0, weight = F, relative = F, plot = TRUE)

# Check the relative number of tree per category.
MicEco::ps_venn(physeq_alb_soil, group = "dominant_tree", fraction = 0, weight = T, relative = T, plot = TRUE)

# Plot the final plot. 
venn_soil_alb <- MicEco::ps_venn(physeq_alb_soil, group = "dominant_tree", 
                                 fraction = 0, weight = F, relative = TRUE, plot = TRUE,
                                 fill = c("green","mediumseagreen"),
                                 labels = list(labels =c("Fagus sylvatica", "Picea abies"),
                                               cex=1.6, font=list(face=3)),
                                 quantities = list(type=c("percent"), 
                                                   labels = c("\n19% (3153)","\n9% (1025)","\n71% (1199)"), cex=1.6))
venn_soil_alb

#Check the absolute number of tree per category.
MicEco::ps_venn(physeq_alb_bark, group = "dominant_tree", fraction = 0, weight = F, relative = F, plot = TRUE)

#Check the relative number of tree per category.
MicEco::ps_venn(physeq_alb_bark, group = "dominant_tree", fraction = 0, weight = T, relative = T, plot = TRUE)

# Plot the final plot. 
venn_bark_alb <- MicEco::ps_venn(physeq_alb_bark, group = "dominant_tree",
                                 fraction = 0, weight = F, relative = F, plot = TRUE,
                                 fill = c("green","mediumseagreen"),
                                 labels = list(labels =c("Fagus sylvatica", "Picea abies"),
                                               cex=1.6,font=list(face=3)),
                                 quantities = list(type=c("percent"),
                                                   labels = c("\n13% (602)","\n9% (299)","\n78% (228)"), cex=1.6))
venn_bark_alb

###
# Venn diagrams by substrate.  
###

# Check the absolute number of tree per category.
MicEco::ps_venn(physeq_alb_fagus, group = "substrate", fraction = 0, weight = F, relative = F, plot = TRUE)

# Check the relative number of tree per category.
MicEco::ps_venn(physeq_alb_fagus, group = "substrate", fraction = 0, weight = T, relative = T, plot = TRUE)

# Plot the final plot. 
venn_fagus_alb <- MicEco::ps_venn(physeq_alb_fagus, group = "substrate",
                                  fraction = 0, weight = F, relative = F, plot = TRUE,
                                  fill = c("green","mediumseagreen"),
                                  labels = list(labels =c("Fagus bark", "Fagus soil"),
                                                cex=1.6,font=list(face=3)),
                                  quantities = list(type=c("percent"),
                                                    labels = c("\n15% (520)","\n29% (4042)","\n56% (310)"), cex=1.6))
venn_fagus_alb

# Check the absolute number of tree per category.
MicEco::ps_venn(physeq_alb_picea, group = "substrate", fraction = 0, weight = F, relative = F, plot = TRUE)

# Check the relative number of tree per category.
MicEco::ps_venn(physeq_alb_picea, group = "substrate", fraction = 0, weight = T, relative = T, plot = TRUE)

# Plot the final plot. 
venn_picea_alb <- MicEco::ps_venn(physeq_alb_picea, group = "substrate",
                                  fraction = 0, weight = F, relative = F, plot = TRUE,
                                  fill = c("green","mediumseagreen"),
                                  labels = list(labels =c("Picea bark", "Picea soil"),
                                                cex=1.6,font=list(face=3)),
                                  quantities = list(type=c("percent"),
                                                    labels = c("\n28% (378)","\n37% (2075)","\n35% (149)"), cex=1.6))
venn_picea_alb

#####################Schorfheide#######################
###
# Venn diagrams by host tree species. 
###

# Check the absolute number of tree per category.
MicEco::ps_venn(physeq_sch_soil, group = "dominant_tree", fraction = 0, weight = F, relative = F, plot = TRUE)

# Check the relative number of tree per category.
MicEco::ps_venn(physeq_sch_soil, group = "dominant_tree", fraction = 0, weight = T, relative = T, plot = TRUE)

# Plot the final plot. 
venn_soil_sch <- MicEco::ps_venn(physeq_sch_soil, group = "dominant_tree", 
                                 fraction = 0, weight = F, relative = TRUE, plot = TRUE,
                                 fill = c("green","mediumseagreen"),
                                 labels = list(labels =c("Fagus sylvatica", "Pinus sylvestris"),
                                               cex=1.6, font=list(face=3)),
                                 quantities = list(type=c("percent"), 
                                                   labels = c("\n8% (1282)","\n4% (734)","\n88% (838)"), cex=1.6))
venn_soil_sch

#Check the absolute number of tree per category.
MicEco::ps_venn(physeq_sch_bark, group = "dominant_tree", fraction = 0, weight = F, relative = F, plot = TRUE)

#Check the relative number of tree per category.
MicEco::ps_venn(physeq_sch_bark, group = "dominant_tree", fraction = 0, weight = T, relative = T, plot = TRUE)

# Plot the final plot. 
venn_bark_sch <- MicEco::ps_venn(physeq_sch_bark, group = "dominant_tree",
                                 fraction = 0, weight = F, relative = F, plot = TRUE,
                                 fill = c("green","mediumseagreen"),
                                 labels = list(labels =c("Fagus sylvatica", "Pinus sylvestris"),
                                               cex=1.6,font=list(face=3)),
                                 quantities = list(type=c("percent"),
                                                   labels = c("\n4% (411)","\n8% (198)","\n88% (162)"), cex=1.6))
venn_bark_sch

###
# Venn diagrams by substrate.  
###

# Check the absolute number of tree per category.
MicEco::ps_venn(physeq_sch_fagus, group = "substrate", fraction = 0, weight = F, relative = F, plot = TRUE)

# Check the relative number of tree per category.
MicEco::ps_venn(physeq_sch_fagus, group = "substrate", fraction = 0, weight = T, relative = T, plot = TRUE)

# Plot the final plot. 
venn_fagus_sch <- MicEco::ps_venn(physeq_sch_fagus, group = "substrate",
                                  fraction = 0, weight = F, relative = F, plot = TRUE,
                                  fill = c("green","mediumseagreen"),
                                  labels = list(labels =c("Fagus bark", "Fagus soil"),
                                                cex=1.6,font=list(face=3)),
                                  quantities = list(type=c("percent"),
                                                    labels = c("\n10% (379)","\n24% (1926)","\n66% (194)"), cex=1.6))
venn_fagus_sch

# Check the absolute number of tree per category.
MicEco::ps_venn(physeq_sch_pinus, group = "substrate", fraction = 0, weight = F, relative = F, plot = TRUE)

# Check the relative number of tree per category.
MicEco::ps_venn(physeq_sch_pinus, group = "substrate", fraction = 0, weight = T, relative = T, plot = TRUE)

# Plot the final plot. 
venn_pinus_sch <- MicEco::ps_venn(physeq_sch_pinus, group = "substrate",
                                  fraction = 0, weight = F, relative = F, plot = TRUE,
                                  fill = c("green","mediumseagreen"),
                                  labels = list(labels =c("Pinus bark", "Pinus soil"),
                                                cex=1.6,font=list(face=3)),
                                  quantities = list(type=c("percent"),
                                                    labels = c("\n22% (261)","\n29% (1473)","\n49% (99)"), cex=1.6))
venn_pinus_sch

#################################################################
##                          Section 7                          ##
##              Community Composition Barplots                 ##
#################################################################
my_cols <- paletteer::paletteer_d('ggsci::default_igv')

############## Swabian Alb ###############
physeq_alb_barplot <- physeq_alb
phyloseq::sample_data(physeq_alb_barplot) <- phyloseq::sample_data(physeq_alb) %>% 
  base::data.frame() %>%  
  dplyr::mutate(tree_substrate = base::paste(dominant_tree, substrate, sep = "-"))


# Subset the phyloseq object to the top 24 orders and put the rest in 
# a category "Others", based on relative abundance. 
phy_alb_ord_top25 <- fantaxtic::top_taxa(physeq_alb_barplot,
                                         tax_level = 'Order',
                                         n_taxa =  24,
                                         by_proportion = TRUE,
                                         merged_label = "Other",
                                         include_na_taxa = T) 
phy_alb_ord_top25_named <- fantaxtic::name_na_taxa(phy_alb_ord_top25$ps_obj, include_rank = T)

# Transform the subset dataset to compositional (relative) abundances.
phy_alb_ord_top25_named_plot <-  phy_alb_ord_top25_named %>%
  microbiome::aggregate_taxa(level = "Order") %>%  
  microbiome::transform(transform = "compositional") 

# Extract the names of the Orders.
taxa_names(phy_alb_ord_top25_named_plot) <- tax_table(phy_alb_ord_top25_named_plot)[, 4]

# Sort the taxa names alphabetically. 
taxa_names_alb_ord <- sort(taxa_names(phy_alb_ord_top25_named_plot))

# To get our desired plotting order and group names we need to change 
# the exploratory names and order them as factors.
sampledata_fungi <- data.frame(sample_data(phy_alb_ord_top25_named_plot))
sampledata_fungi <- sampledata_fungi %>% 
  mutate_at("tree_substrate",funs(str_replace(., "Fagus_sylvatica-bark", "F. sylvatica bark"))) %>% 
  mutate_at("tree_substrate",funs(str_replace(., "Fagus_sylvatica-soil", "F. sylvatica soil"))) %>% 
  mutate_at("tree_substrate",funs(str_replace(., "Picea_abies-bark", "P. abies bark"))) %>% 
  mutate_at("tree_substrate",funs(str_replace(., "Picea_abies-soil", "P. abies soil"))) 
sampledata_fungi$tree_substrate <- factor(sampledata_fungi$tree_substrate, 
                                       levels = c("F. sylvatica bark", "P. abies bark",
                                                  "F. sylvatica soil", "P. abies soil"))  

sample_data(phy_alb_ord_top25_named_plot) <- sample_data(sampledata_fungi)

# Custom plotting to make a nice stacked barplot. 
alb_ord_plots <- phy_alb_ord_top25_named_plot %>%
  microbiome::plot_composition(group_by =  'tree_substrate', otu.sort = taxa_names_alb_ord) +
  scale_fill_manual(values = ggplot2::alpha(my_cols, 0.9), name = 'Order') +
  guides(fill = guide_legend(title.position = 'top')) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = 'black', size = 0.5),
        axis.text.x =  element_blank(),
        axis.text.y =  element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 10),
        legend.position = 'bottom', 
        plot.title = element_text(vjust = -4, hjust = 0.03), 
        legend.text = element_text(colour = 'black', size = 7),
        legend.title =  element_text(size = 10),
        legend.key.size = unit(2.5, 'mm'),
        axis.ticks.length.x = unit(-0.2, "cm"), 
        legend.box.spacing = unit(-4, 'mm'),
        legend.background = element_rect(fill = 'transparent'),
        text = element_text(colour = 'black')) + 
  xlab('Sample') +
  ylab('Relative Abundance')   + 
  labs( subtitle = 'Swabian Alb')
alb_ord_plots


############## Schorfheide-Chorin ###############
physeq_sch_barplot <- physeq_sch
phyloseq::sample_data(physeq_sch_barplot) <- phyloseq::sample_data(physeq_sch) %>% 
  base::data.frame() %>%  
  dplyr::mutate(tree_substrate = base::paste(dominant_tree, substrate, sep = "-"))


# Subset the phyloseq object to the top 24 orders and put the rest in 
# a category "Others", based on relative abundance. 
phy_sch_ord_top25 <- fantaxtic::top_taxa(physeq_sch_barplot,
                                         tax_level = 'Order',
                                         n_taxa =  24,
                                         by_proportion = TRUE,
                                         merged_label = "Other",
                                         include_na_taxa = T) 
phy_sch_ord_top25_named <- fantaxtic::name_na_taxa(phy_sch_ord_top25$ps_obj, include_rank = T)

# Transform the subset dataset to compositional (relative) abundances.
phy_sch_ord_top25_named_plot <-  phy_sch_ord_top25_named %>%
  microbiome::aggregate_taxa(level = "Order") %>%  
  microbiome::transform(transform = "compositional") 

# Extract the names of the Orders.
taxa_names(phy_sch_ord_top25_named_plot) <- tax_table(phy_sch_ord_top25_named_plot)[, 4]

# Sort the taxa names alphabetically. 
taxa_names_sch_ord <- sort(taxa_names(phy_sch_ord_top25_named_plot))

# To get our desired plotting order and group names we need to change 
# the exploratory names and order them as factors.
sampledata_fungi <- data.frame(sample_data(phy_sch_ord_top25_named_plot))
sampledata_fungi <- sampledata_fungi %>% 
  mutate_at("tree_substrate",funs(str_replace(., "Fagus_sylvatica-bark", "F. sylvatica bark"))) %>% 
  mutate_at("tree_substrate",funs(str_replace(., "Fagus_sylvatica-soil", "F. sylvatica soil"))) %>% 
  mutate_at("tree_substrate",funs(str_replace(., "Pinus_sylvestris-bark", "P. sylvestris bark"))) %>% 
  mutate_at("tree_substrate",funs(str_replace(., "Pinus_sylvestris-soil", "P. sylvestris soil"))) 
sampledata_fungi$tree_substrate <- factor(sampledata_fungi$tree_substrate, 
                                          levels = c("F. sylvatica bark", "P. sylvestris bark",
                                                     "F. sylvatica soil", "P. sylvestris soil"))  

sample_data(phy_sch_ord_top25_named_plot) <- sample_data(sampledata_fungi)

# Custom plotting to make a nice stacked barplot. 
sch_ord_plots <- phy_sch_ord_top25_named_plot %>%
  microbiome::plot_composition(group_by =  'tree_substrate', otu.sort = taxa_names_sch_ord) +
  scale_fill_manual(values = ggplot2::alpha(my_cols, 0.9), name = 'Order') +
  guides(fill = guide_legend(title.position = 'top')) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = 'black', size = 0.5),
        axis.text.x =  element_blank(),
        axis.text.y =  element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 10),
        legend.position = 'bottom', 
        plot.title = element_text(vjust = -4, hjust = 0.03), 
        legend.text = element_text(colour = 'black', size = 7),
        legend.title =  element_text(size = 10),
        legend.key.size = unit(2.5, 'mm'),
        axis.ticks.length.x = unit(-0.2, "cm"), 
        legend.box.spacing = unit(-4, 'mm'),
        legend.background = element_rect(fill = 'transparent'),
        text = element_text(colour = 'black')) + 
  xlab('Sample') +
  ylab('Relative Abundance')   + 
  labs( subtitle = 'Schorfheide-Chorin')
sch_ord_plots














































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

