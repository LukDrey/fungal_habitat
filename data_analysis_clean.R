#################################################################
##                          Section 1                          ##
##                       Package Loading                       ##
#################################################################

library(phyloseq)

library(ggplot2)

library(MicEco)

library(tidyverse)

library(ggpubr)

library(vegan)

library(paletteer)

library(ranacapa)

library(fantaxtic)

library(SpiecEasi)

library(rgexf)

library(igraph)

library(here)

library(microbiome)

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
  phyloseq::subset_samples(exploratory %in% c("Alb", "Schorfheide")) %>%
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
rare_sch <- ranacapa::ggrare(physeq_sch_curve, step = 50,
                             color = "tree_substrate", se = FALSE) 

combined_rare_curves <- ggpubr::ggarrange(rare_alb, rare_sch,
                                          ncol = 2, nrow = 1)
combined_rare_curves

ggsave('combined_rare_curves.tiff', device = 'tiff',
       combined_rare_curves, width = 400, height = 240,
       units = 'mm', dpi = 300)  
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
##                          Section 5                          ##
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
full_ordination_ps <- filtered_physeq
sample_data(full_ordination_ps) <- data.frame(sample_data(full_ordination_ps)) %>% 
  dplyr::rename(habitat = substrate)

nmds_full <- phyloseq::ordinate(full_ordination_ps, method = "NMDS", distance = "bray")

ordination_full <- phyloseq::plot_ordination(full_ordination_ps, nmds_full,
                                             type="samples", color="dominant_tree", shape="habitat") + 
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
ordination_final

ggsave('ordination_final.tiff', device = 'tiff',
       ordination_final, width = 400, height = 240,
       units = 'mm', dpi = 300)

#################################################################
##                          Section 6                          ##
##                    Variance partitioning                    ##
#################################################################
##---------------------------------------------------------------
##                          Section 6.1                         -
##                        Alpha Diversity                       -
##---------------------------------------------------------------

################Swabian Alb###########################

# Remind us of the table we created earlier cotaining our metadata and shannon diversity.
head(div_data_alb)

# Create linear models lm of our alpha diversity measure explained by tree and substrate. 
lm_full_alb <- lm(Shannon ~ dominant_tree + substrate, 
                  data = div_data_alb, na.action = "na.fail")

base::summary(lm_full_alb)  
base::plot(effects::allEffects(lm_full_alb))

# Check model fit. Looks okay. 
base::plot(stats::residuals(lm_full_alb), stats::fitted(lm_full_alb))
stats::qqnorm(stats::residuals(lm_full_alb))
stats::qqline(stats::residuals(lm_full_alb))

# Now create linear models of the explanatory variables alone.  
lm_substrate_alb <- lm(Shannon ~ substrate, 
                  data = div_data_alb, na.action = "na.fail")

lm_tree_alb <- lm(Shannon ~ dominant_tree, 
                  data = div_data_alb, na.action = "na.fail")

# Now partition the variance. 
alb_variance_lm <- modEvA::varPart(A = summary(lm_substrate_alb)$r.squared,
                                 B = summary(lm_tree_alb)$r.squared,
                                 AB = summary(lm_full_alb)$r.squared,
                                 A.name = "Substrate",
                                 B.name = "Tree") %>% 
  mutate(across(Proportion, round, 3)) %>% 
  tibble::rownames_to_column() %>%  
  tibble::add_column(variable = c("Habitat", "Dominant tree", "Overlap", "Unexplained")) %>% 
  dplyr::select(variable, Proportion)
alb_variance_lm

################Schorfheide-Chorin###########################

# Remind us of the table we created earlier cotaining our metadata and shannon diversity.
head(div_data_sch)

# Create linear models lm of our alpha diversity measure explained by tree and substrate. 
lm_full_sch <- lm(Shannon ~ dominant_tree + substrate, 
                  data = div_data_sch, na.action = "na.fail")

base::summary(lm_full_sch)  
base::plot(effects::allEffects(lm_full_sch))

# Check model fit. Looks okay. 
base::plot(stats::residuals(lm_full_sch), stats::fitted(lm_full_sch))
stats::qqnorm(stats::residuals(lm_full_sch))
stats::qqline(stats::residuals(lm_full_sch))

# Now create linear models of the explanatory variables alone.  
lm_substrate_sch <- lm(Shannon ~ substrate, 
                       data = div_data_sch, na.action = "na.fail")

lm_tree_sch <- lm(Shannon ~ dominant_tree, 
                  data = div_data_sch, na.action = "na.fail")

# Now partition the variance. 
sch_variance_lm <- modEvA::varPart(A = summary(lm_substrate_sch)$r.squared,
                                   B = summary(lm_tree_sch)$r.squared,
                                   AB = summary(lm_full_sch)$r.squared,
                                   A.name = "Substrate",
                                   B.name = "Tree") %>% 
  mutate(across(Proportion, round, 3)) %>% 
  tibble::rownames_to_column() %>%  
  tibble::add_column(variable = c("Habitat", "Dominant tree", "Overlap", "Unexplained")) %>% 
  dplyr::select(variable, Proportion)
sch_variance_lm

##---------------------------------------------------------------
##                          Section 6.2                         -
##                         Beta Diversity                       -
##---------------------------------------------------------------
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

indfract_adj_r_alb <- varp_alb$part$indfract %>% 
  tibble::rownames_to_column()  %>%  
  tibble::add_column(variable = c("Habitat", "Overlap", "Dominant tree", "Unexplained")) %>% 
  dplyr::select(variable, Adj.R.squared) %>% 
  dplyr::rename(Proportion = Adj.R.squared) %>% 
  dplyr::mutate(across(Proportion, .fns = round, 3))
  
indfract_adj_r_alb

### Schorfheide ###

# Extract OTU table Schorfheide.
otu_sch <- veganotu(physeq_sch)

# Extract sample data.
data_sch <- base::data.frame(phyloseq::sample_data(physeq_sch))

# Get the variation with vegan::varpart.
varp_sch <- vegan::varpart(otu_sch, ~ substrate, ~ dominant_tree, data = data_sch)

indfract_adj_r_sch <- varp_sch$part$indfract %>% 
  tibble::rownames_to_column()  %>%  
  tibble::add_column(variable = c("Habitat", "Overlap", "Dominant tree", "Unexplained")) %>% 
  dplyr::select(variable, Adj.R.squared) %>% 
  dplyr::rename(Proportion = Adj.R.squared) %>% 
  dplyr::mutate(across(Proportion, .fns = round, 3))

indfract_adj_r_sch

##---------------------------------------------------------------
##                          Section 6.3                         -
##                            Plotting                          -
##---------------------------------------------------------------

# Combine alpha into one dataframe 

alb_variance_alpha <- alb_variance_lm  %>% 
  cbind(div_lev = rep("\u03B1-Diversity", nrow(alb_variance_lm))) %>% 
  cbind(exploratory = rep("Swabian Alb", nrow(alb_variance_lm)))

sch_variance_alpha <- sch_variance_lm  %>% 
  cbind(div_lev = rep("\u03B1-Diversity", nrow(sch_variance_lm))) %>% 
  cbind(exploratory = rep("Schorfheide-Chorin", nrow(sch_variance_lm)))

variance_lm <- rbind(alb_variance_alpha, sch_variance_alpha) %>% 
  dplyr::rename(variance = Proportion) %>% 
  dplyr::mutate(variance = (variance * 100))

# Combine beta into one dataframe 

alb_variance_nmds <- indfract_adj_r_alb  %>% 
  cbind(div_lev = rep("\u03B2-Diversity", nrow(indfract_adj_r_alb))) %>% 
  cbind(exploratory = rep("Swabian Alb", nrow(indfract_adj_r_alb)))

sch_variance_nmds <- indfract_adj_r_sch  %>% 
  cbind(div_lev = rep("\u03B2-Diversity", nrow(indfract_adj_r_sch))) %>% 
  cbind(exploratory = rep("Schorfheide-Chorin", nrow(indfract_adj_r_sch)))

variance_beta <- rbind(alb_variance_nmds, sch_variance_nmds) %>% 
  dplyr::rename(variance = Proportion) %>% 
  dplyr::mutate(variance = (variance * 100))

variance_full <- rbind(variance_lm, variance_beta)

# Create faceted multipanel plot with ggpubr and ggplot2. 

div_labels <- c("\u03B1-Diversity", "\u03B2-Diversity")

variance_barplot <- ggpubr::ggbarplot(variance_full, x = "exploratory", y = "variance", 
                  fill = "variable", color = "variable",
                  palette = ggplot2::alpha(c("#009E73",  "#E69F00",  
                                             "#56B4E9", "#999999"))) +
  ggplot2::facet_grid(exploratory ~ div_lev,
                      space="free", scales="free", switch = "y") +
  ggplot2::scale_x_discrete(position = "top") +
  ggplot2::geom_hline(aes(yintercept = 0), linetype = "dashed") +
  ggplot2::coord_flip() + 
  ggplot2::ylab("Explained Variance") + 
  ggplot2::theme(text = element_text(size = 15),
                 legend.title = element_blank(),
                 legend.position = c(1.05,0.5),
                 axis.title.y = element_blank(),
                 axis.line.y = element_blank(), 
                 axis.ticks.y = element_blank()) 
variance_barplot

ggsave('variance_barplot.tiff', device = 'tiff',
       variance_barplot, width = 400, height = 240,
       units = 'mm', dpi = 300)

#################################################################
##                          Section 7                          ##
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

# Create the final figure.

final_venn_diagrams <- ggpubr::ggarrange(venn_bark_alb, venn_soil_alb,
                                         venn_bark_sch, venn_soil_sch,
                                         venn_fagus_alb, venn_picea_alb,
                                         venn_fagus_sch, venn_pinus_sch,
                                         ncol = 2, nrow = 4, 
                                         labels = "AUTO")
final_venn_diagrams

ggsave('final_venn_diagrams.tiff', device = 'tiff',
       final_venn_diagrams, width = 400, height = 600,
       units = 'mm', dpi = 300)

#################################################################
##                          Section 8                          ##
##              Community Composition Barplots                 ##
#################################################################

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
phyloseq::taxa_names(phy_alb_ord_top25_named_plot) <- phyloseq::tax_table(phy_alb_ord_top25_named_plot)[, 4]

# Sort the taxa names alphabetically. 
taxa_names_alb_ord <- base::sort(phyloseq::taxa_names(phy_alb_ord_top25_named_plot))

# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
sampledata_alb <- base::data.frame(phyloseq::sample_data(phy_alb_ord_top25_named_plot))
sampledata_alb <- sampledata_alb %>% 
  mutate(across("tree_substrate", stringr::str_replace, "Fagus_sylvatica-bark", "F. sylvatica bark")) %>% 
  mutate(across("tree_substrate", stringr::str_replace, "Fagus_sylvatica-soil", "F. sylvatica soil")) %>% 
  mutate(across("tree_substrate", stringr::str_replace, "Picea_abies-bark", "P. abies bark")) %>% 
  mutate(across("tree_substrate", stringr::str_replace, "Picea_abies-soil", "P. abies soil")) 
sampledata_alb$tree_substrate <- factor(sampledata_alb$tree_substrate, 
                                       levels = c("F. sylvatica bark", "P. abies bark",
                                                  "F. sylvatica soil", "P. abies soil"))  

phyloseq::sample_data(phy_alb_ord_top25_named_plot) <- phyloseq::sample_data(sampledata_alb)


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
phyloseq::taxa_names(phy_sch_ord_top25_named_plot) <- phyloseq::tax_table(phy_sch_ord_top25_named_plot)[, 4]

# Sort the taxa names alphabetically. 
taxa_names_sch_ord <- sort(phyloseq::taxa_names(phy_sch_ord_top25_named_plot))

# To get our desired plotting order and group names we need to change 
# the exploratory names and order them as factors.
sampledata_sch <- data.frame(phyloseq::sample_data(phy_sch_ord_top25_named_plot))
sampledata_sch <- sampledata_sch %>% 
  mutate(across("tree_substrate", stringr::str_replace, "Fagus_sylvatica-bark", "F. sylvatica bark")) %>% 
  mutate(across("tree_substrate", stringr::str_replace, "Fagus_sylvatica-soil", "F. sylvatica soil")) %>% 
  mutate(across("tree_substrate", stringr::str_replace, "Pinus_sylvestris-bark", "P. sylvestris bark")) %>% 
  mutate(across("tree_substrate", stringr::str_replace, "Pinus_sylvestris-soil", "P. sylvestris soil"))
          
sampledata_sch$tree_substrate <- factor(sampledata_sch$tree_substrate, 
                                        levels = c("F. sylvatica bark", "P. sylvestris bark",
                                                   "F. sylvatica soil", "P. sylvestris soil"))  

phyloseq::sample_data(phy_sch_ord_top25_named_plot) <- phyloseq::sample_data(sampledata_sch)


#################################################################
##                          Section 8.2                        ##
##              Create great looking plots                     ##
#################################################################

my_cols <- Polychrome::createPalette(33, seedcolors = c("#ff0000", "#00ff00", "#0000ff"))

unique_orders <- sort(unique(c(taxa_names_alb_ord, taxa_names_sch_ord))) 

custom_sort <- c("Agaricales", "Archaeorhizomycetales",
                 "Atheliales", "Boletales",
                 "Caliciales", "Cantharellales",
                 "Capnodiales", "Chaetothyriales",
                 "Eurotiales", "Filobasidiales",
                 "Helotiales", "Hypocreales",
                 "Lecanorales", "Mortierellales",
                 "Mycosphaerellales", "Mytilinidales",
                 "Orbiliales", "Pezizales", "Phaeothecales",
                 "Pleosporales", "Russulales",
                 "Sebacinales", "Thelebolales",
                 "Thelephorales", "Trapeliales", 
                 "Tremellales", "Umbelopsidales", "Verrucariales",
                 "Unknown Ascomycota (Phylum)", "Unknown Dothideomycetes (Class)",
                 "Unknown Fungi (Kingdom)", "Unknown Rozellomycota (Phylum)",
                              "Other")

full_cols <- data.frame(order = custom_sort, color = my_cols)

alb_cols <- full_cols %>% 
  dplyr::filter(order %in% taxa_names_alb_ord)

sch_cols <- full_cols %>% 
  dplyr::filter(order %in% taxa_names_sch_ord)


scale_fill_alb <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = setNames(alb_cols$color, alb_cols$order), 
    ...
  )
}

scale_fill_sch <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = setNames(sch_cols$color, sch_cols$order), 
    ...
  )
}

############## Swabian Alb ###############

# Custom plotting to make a nice stacked barplot. 
alb_ord_soil_plots <- phyloseq::subset_samples(phy_alb_ord_top25_named_plot, substrate == "soil") %>%
  microbiome::plot_composition(group_by =  'tree_substrate', otu.sort = alb_cols$order) +
  scale_fill_alb() +
  guides(fill = guide_legend(title.position = 'top', nrow = 3)) +
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
        text = element_text(colour = 'black', size = 20),
        strip.text = element_text(face = "italic")) + 
  xlab('Sample') +
  ylab('Relative Abundance') + 
  labs(subtitle = "(B)")  
alb_ord_soil_plots

# Custom plotting to make a nice stacked barplot. 
alb_ord_bark_plots <- phyloseq::subset_samples(phy_alb_ord_top25_named_plot, substrate == "bark") %>%
  microbiome::plot_composition(group_by =  'tree_substrate', otu.sort = alb_cols$order) +
  scale_fill_alb() +
  guides(fill = guide_legend(title.position = 'top', ncol = 10)) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = 'black', size = 0.5),
        axis.text.x =  element_blank(),
        axis.text.y =  element_text(colour = "black", size = 10),
        axis.title.x = element_text(colour = "black", size = 10),
        axis.title.y = element_text(colour = "black", size = 10),
        legend.position = 'bottom', 
        plot.title = element_text(vjust = -4, hjust = 0.03), 
        legend.text = element_text(colour = 'black', size = 7),
        legend.title =  element_text(size = 10),
        legend.key.size = unit(2.5, 'mm'),
        axis.ticks.length.x = unit(-0.2, "cm"), 
        legend.box.spacing = unit(-4, 'mm'),
        legend.background = element_rect(fill = 'transparent'),
        text = element_text(colour = 'black', size = 20),
        strip.text = element_text(face = "italic")) + 
  xlab('Sample') +
  ylab('Relative Abundance') + 
  labs(title = 'Swabian Alb', subtitle = "(A)")
alb_ord_bark_plots

############## Schorfheide-Chorin ###############

# Custom plotting to make a nice stacked barplot. 
sch_ord_soil_plots <- phyloseq::subset_samples(phy_sch_ord_top25_named_plot, substrate == "soil") %>% 
  microbiome::plot_composition(group_by =  'tree_substrate', otu.sort = sch_cols$order) +
  scale_fill_sch() +
  guides(fill = guide_legend(title.position = 'top', ncol = 10)) +
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
        text = element_text(colour = 'black', size = 20),
        strip.text = element_text(face = "italic")) + 
  xlab('Sample') +
  ylab('Relative Abundance') + 
  labs(subtitle = "(B)")   
sch_ord_soil_plots

sch_ord_bark_plots <- phyloseq::subset_samples(phy_sch_ord_top25_named_plot, substrate == "bark") %>% 
  microbiome::plot_composition(group_by =  'tree_substrate', otu.sort = sch_cols$order) +
  scale_fill_sch() +
  guides(fill = guide_legend(title.position = 'top', ncol = 10)) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = 'black', size = 0.5),
        axis.text.x =  element_blank(),
        axis.text.y =  element_text(colour = "black", size = 10),
        axis.title.x = element_text(colour = "black", size = 10),
        axis.title.y = element_text(colour = "black", size = 10),
        legend.position = 'bottom', 
        plot.title = element_text(vjust = -4, hjust = 0.03), 
        legend.text = element_text(colour = 'black', size = 7),
        legend.title =  element_text(size = 10),
        legend.key.size = unit(2.5, 'mm'),
        axis.ticks.length.x = unit(-0.2, "cm"), 
        legend.box.spacing = unit(-4, 'mm'),
        legend.background = element_rect(fill = 'transparent'),
        text = element_text(colour = 'black', size = 20),
        strip.text = element_text(face = "italic")) + 
  xlab('Sample') +
  ylab('Relative Abundance') + 
  labs(title = 'Schorfheide-Chorin', subtitle = "(A)") 
sch_ord_bark_plots

######## Create final arranged plot #######

final_alb_community_barplot <- ggpubr::ggarrange(alb_ord_bark_plots, alb_ord_soil_plots,
                                                 ncol = 1, nrow = 2, hjust = 5,
                                                 legend = "bottom", common.legend = TRUE)
final_alb_community_barplot

ggsave('final_alb_community_barplot.tiff', device = 'tiff',
       final_alb_community_barplot, width = 400, height = 300,
       units = 'mm', dpi = 300)

final_sch_community_barplot <- ggpubr::ggarrange(sch_ord_bark_plots, sch_ord_soil_plots,
                                                 ncol = 1, nrow = 2,
                                                 legend = "bottom", common.legend = TRUE)
final_sch_community_barplot

ggsave('final_sch_community_barplot.tiff', device = 'tiff',
       final_sch_community_barplot, width = 400, height = 300,
       units = 'mm', dpi = 300)


combined_community_barplot <- ggpubr::ggarrange(final_alb_community_barplot,
                                                final_sch_community_barplot, 
                                                ncol = 1, nrow = 2)
combined_community_barplot

ggsave('combined_community_barplots.tiff', device = 'tiff',
       combined_community_barplot, width = 300, height = 600,
       units = 'mm', dpi = 300)

#################################################################
##                          Section 9                          ##
##                     Co-Occurence Networks                   ##
#################################################################

################Swabian Alb#######################

##subset Alb
#subset dataset to ASVs that occur in at least one percent of samples

phy_trans_A  <- phyloseq::transform_sample_counts(physeq_alb, function(x) x / sum(x) )
phy_abundfilt_A <- phyloseq::filter_taxa(phy_trans_A, function(x) sum(x) > .01, TRUE)
keep_names_A <- taxa_names(phy_abundfilt_A)

my_subset_A <- subset(otu_table(physeq_alb), rownames(otu_table(physeq_alb)) %in% keep_names_A)
phy_raw_abundfilt_A <- merge_phyloseq(my_subset_A, tax_table(physeq_alb), sample_data(physeq_alb))

#saveRDS(phy_alb_raw_abundfilt_A, 'phy_alb_raw_abundfilt.rds')

# # The following part was done on a massive RAM server for faster computation.
# library(SpiecEasi)
# 
# # Read in the abundance filtered phyloseq object. 
# phy_raw_abundfilt_A <- readRDS('phy_alb_raw_abundfilt.rds')
# 
# # Set up the parameters for the pulsar package 
# pargs <- list(rep.num=50, seed=10010)
# 
# # Run the spieceasi algorithm  
# se_raw_abundfilt_A <- SpiecEasi::spiec.easi(phy_raw_abundfilt_A, method='mb', 
#                                  lambda.min.ratio=1e-2, nlambda=70, 
#                                  sel.criterion='bstars',  pulsar.select=TRUE,
#                                  pulsar.params=pargs)
# 
# # Obtain the stability parameters.
# SpiecEasi::getStability(se_raw_abundfilt_A)
# 
# # Check if we have a none-zero network. 
# sum(SpiecEasi::getRefit(se_raw_abundfilt_A))/2
# 
# # Save the network.
# saveRDS(se_raw_abundfilt_A, 'se_raw_abundfilt_A.rds')

# Now we move back to our local computer.
# Read in the network again.
se_raw_abundfilt_A <- base::readRDS(here::here('se_raw_abundfilt_A.rds'))

# Convert it to an igraph object.
fun.mb_A <- SpiecEasi::adj2igraph(SpiecEasi::getRefit(se_raw_abundfilt_A),  
                       vertex.attr=list(name = phyloseq::taxa_names(phy_raw_abundfilt_A)))

# To obtain modules and get a nice plot we need to move to the OSS Gephi.
# Convert the igraph object to gephi format.
fun.mb_A_gephi <- rgexf::igraph.to.gexf(fun.mb_A)

# Write the object. 
rgexf::write.gexf(fun.mb_A_gephi, output = here::here("alb_network.gexf"))

# Get the figures and modules from Gephi. See the paper for details on the 
# algorithms we used inside Gephi. 

# Obtain the hub taxa based on Kleinbergs betweenness centrality.
fun_between_A <- igraph::betweenness(fun.mb_A, directed = F)

# Sort them by their values and subset to the top 5.
fun_top5_between_A <- base::sort(fun_between_A, decreasing = T)[1:5]

# Get the names.
base::names(fun_top5_between_A)

# Obtain the taxonomy from the phyloseq object. 
hub_taxa_fun_A <- phyloseq::subset_taxa(physeq_alb,
                                        phyloseq::taxa_names(physeq_alb)
                                        %in% base::names(fun_top5_between_A))
phyloseq::tax_table(hub_taxa_fun_A)


################Schorfheide-Chorin#######################

# Subset dataset to ASVs that occur in at least one percent of samples

phy_trans_S  <- phyloseq::transform_sample_counts(physeq_sch, 
                                                  function(x) x / sum(x))
phy_abundfilt_S <- phyloseq::filter_taxa(phy_trans_S, 
                                         function(x) sum(x) > .01, TRUE)
keep_names_S <- taxa_names(phy_abundfilt_S)

my_subset_S <- subset(otu_table(physeq_sch), 
                      rownames(otu_table(physeq_sch)) %in% keep_names_S)
phy_raw_abundfilt_S <- merge_phyloseq(my_subset_S, 
                                      tax_table(physeq_sch), 
                                      sample_data(physeq_sch))

# saveRDS(phy_raw_abundfilt_S, 'phy_sch_raw_abundfilt.rds')
# 
# # The following part was done on a massive RAM server for faster computation.
# # Read in the abundance filtered phyloseq object. 
# phy_raw_abundfilt_S <- readRDS('phy_sch_raw_abundfilt.rds')
# 
# # Set up the parameters for the pulsar package 
# pargs <- list(rep.num=50, seed=10010)
# 
# # Run the spieceasi algorithm  
# se_raw_abundfilt_S <- SpiecEasi::spiec.easi(phy_raw_abundfilt_S, method='mb', 
#                                             lambda.min.ratio=1e-2, nlambda=100, 
#                                             sel.criterion='bstars',  pulsar.select=TRUE,
#                                             pulsar.params=pargs)
# 
# # Obtain the stability parameters.
# SpiecEasi::getStability(se_raw_abundfilt_S)
# 
# # Check if we have a none-zero network. 
# sum(SpiecEasi::getRefit(se_raw_abundfilt_S))/2
# 
# # Save the network.
# saveRDS(se_raw_abundfilt_S, 'se_raw_abundfilt_S.rds')

# Now we move back to our local computer.
# Read in the network again.
se_raw_abundfilt_S <- base::readRDS(here::here('se_raw_abundfilt_S.rds'))

# Convert it to an igraph object.
fun.mb_S <- SpiecEasi::adj2igraph(SpiecEasi::getRefit(se_raw_abundfilt_S),  
                                  vertex.attr=list(name = phyloseq::taxa_names(phy_raw_abundfilt_S)))

# To obtain modules and get a nice plot we need to move to the OSS Gephi.
# Convert the igraph object to gephi format.
fun.mb_S_gephi <- rgexf::igraph.to.gexf(fun.mb_S)

# Write the object. 
rgexf::write.gexf(fun.mb_S_gephi, output = here::here("sch_network.gexf"))

# Get the figures and modules from Gephi. See the paper for details on the 
# algorithms we used inside Gephi. 

# Obtain the hub taxa based on Kleinbergs betweenness centrality.
fun_between_S <- igraph::betweenness(fun.mb_S, directed = F)

# Sort them by their values and subset to the top 5.
fun_top5_between_S <- base::sort(fun_between_S, decreasing = T)[1:5]
# Get the names.
names(fun_top5_between_S)

# Obtain the taxonomy from the phyloseq object. 
hub_taxa_fun_S <- phyloseq::subset_taxa(physeq_sch, 
                                        phyloseq::taxa_names(physeq_sch) 
                                        %in% base::names(fun_top5_between_S))
phyloseq::tax_table(hub_taxa_fun_S)

##----------------------------------------------------------------
##                          Section 9.1                          -
##          Relative abundances of habitats in network           -
##----------------------------------------------------------------

################Swabian Alb#######################
# Load in the modules obtained from Gephi. 
modules_alb <- utils::read.csv(here::here("modules_alb.csv"))

# Append this to the phyloseq object from which the networks were created. 
# Append it to the taxonomy.
tax_alb <- base::data.frame(phyloseq::tax_table(phy_raw_abundfilt_A)) %>% 
  tibble::rownames_to_column(var = "Label") %>% 
  dplyr::left_join(., modules_alb)

# Necessary steps to get it into the phyloseq object for later plotting.  
row.names(tax_alb) <- tax_alb$Label
tax_alb$Label <- NULL
taxmat_alb <- as.matrix(tax_alb)

# Create a new phyloseq object in order to not overwrite the old one.
physeq_alb_modules <- phy_raw_abundfilt_A

# Update the taxonomy table.
tax_table(physeq_alb_modules) <- tax_table(taxmat_alb)

# Update the sample_data to include tree_substrate variable.
phyloseq::sample_data(physeq_alb_modules) <- phyloseq::sample_data(sampledata_alb)

# Transform the ASV counts to relative abundance. 
physeq_alb_modules <- phyloseq::transform_sample_counts(physeq_alb_modules,
                                                        function(x){x / sum(x)})

alb_rel_abund_modules_barplot <- phyloseq::plot_bar(physeq_alb_modules, x = "modularity_class", fill = "modularity_class")+
  facet_wrap(~ tree_substrate, scales = "free", ncol = 1) +
  guides(fill = "none") +
  theme(strip.text.x = element_text(size = 20, face = "italic"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.position = "bottom",legend.direction = "horizontal",
        axis.title.y = element_text(size = 20),
        axis.text.x=element_blank(),
        axis.text.y = element_text(size = 16)) +
  scale_fill_manual(values = c("#51FEFF","#FE8334","#00D61C","#FF61EA")) +
  scale_x_discrete(limits = c("1", "3", "0", "2")) +
  geom_bar(stat="identity") +
  labs(x = "", y = "proportion of reads in %") +
  ylim(0, 30) 
alb_rel_abund_modules_barplot

ggpubr::ggexport(alb_rel_abund_modules_barplot, filename = "alb_modules_plot.tiff",
                 width = 900, height = 1800, 
                 res = 300)

################Schorfheide-Chorin#######################

# Load in the modules obtained from Gephi. 
modules_sch <- utils::read.csv(here::here("modules_sch.csv"))

# Append this to the phyloseq object from which the networks were created. 
# Append it to the taxonomy.
tax_sch <- base::data.frame(phyloseq::tax_table(phy_raw_abundfilt_S)) %>% 
  tibble::rownames_to_column(var = "Label") %>% 
  dplyr::left_join(., modules_sch)

# Necessary steps to get it into the phyloseq object for later plotting.  
row.names(tax_sch) <- tax_sch$Label
tax_sch$Label <- NULL
taxmat_sch <- as.matrix(tax_sch)

# Create a new phyloseq object in order to not overwrite the old one.
physeq_sch_modules <- phy_raw_abundfilt_S

# Update the taxonomy table.
tax_table(physeq_sch_modules) <- tax_table(taxmat_sch)

# Update the sample_data to include tree_substrate variable.
phyloseq::sample_data(physeq_sch_modules) <- phyloseq::sample_data(sampledata_sch)

# Transform the ASV counts to relative abundance. 
physeq_sch_modules <- phyloseq::transform_sample_counts(physeq_sch_modules,
                                                        function(x){x / sum(x)})

sch_rel_abund_modules_barplot <- plot_bar(physeq_sch_modules, x="modularity_class", fill = "modularity_class") +
  facet_wrap(~ tree_substrate, scales = "free", ncol = 1) +
  guides(fill = "none") +
  theme(strip.text.x = element_text(size = 20, face = "italic"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.position = "bottom",legend.direction = "horizontal",
        axis.title.y = element_text(size = 20),
        axis.text.x=element_blank(),
        axis.text.y = element_text(size = 16),
        legend.key = element_blank()) +
  scale_fill_manual(values = c("#FE8334", "#FF61EA", "#51FEFF","#00D61C")) +
  geom_bar(stat="identity") +
  labs(x = "", y = "proportion of reads in %") +
  ylim(0, 18)
sch_rel_abund_modules_barplot

ggpubr::ggexport(sch_rel_abund_modules_barplot, filename = "sch_modules_plot.tiff",
                 width = 900, height = 1800, 
                 res = 300)

#################################################################
##                          Section 10                         ##
##                    Miscellaneous Numbers                    ##
#################################################################

# Total number of ASVs
asv_num_total <- phyloseq::ntaxa(filtered_physeq)
asv_num_total

# Number of ASVs per study region. 
asv_num_alb <- phyloseq::ntaxa(physeq_alb)
asv_num_alb

asv_num_sch <- phyloseq::ntaxa(physeq_sch)
asv_num_sch

# Number of ASVs per substrate. 
asv_num_soil <- phyloseq::ntaxa(physeq_soil)
asv_num_soil

asv_num_bark <- phyloseq::ntaxa(physeq_bark)
asv_num_bark

# Number & percentage of ASVs not assignable at all taxonomic rank.
no_phylum_num <- sum(is.na(data.frame(phyloseq::tax_table(filtered_physeq))$Phylum))
no_phylum_num
no_phylum_num/asv_num_total * 100

no_class_num <- sum(is.na(data.frame(phyloseq::tax_table(filtered_physeq))$Class))
no_class_num
no_class_num/asv_num_total * 100

no_order_num <- sum(is.na(data.frame(phyloseq::tax_table(filtered_physeq))$Order))
no_order_num
no_order_num/asv_num_total * 100

no_family_num <- sum(is.na(data.frame(phyloseq::tax_table(filtered_physeq))$Family))
no_family_num
no_family_num/asv_num_total * 100

no_genus_num <- sum(is.na(data.frame(phyloseq::tax_table(filtered_physeq))$Genus))
no_genus_num
no_genus_num/asv_num_total * 100

no_species_num <- sum(is.na(data.frame(phyloseq::tax_table(filtered_physeq))$Species))
no_species_num
no_species_num/asv_num_total * 100

# Overlap between soil and bark in the regions. 
# Check the relative number of ASVs for the Swabian Alb.
alb_substrate_overlap <- MicEco::ps_venn(physeq_alb, group = "substrate", fraction = 0,
                weight = T, relative = T, plot = T)
alb_substrate_overlap

# Check the relative number of ASVs for Schorfheide-Chorin.
sch_substrate_overlap <- MicEco::ps_venn(physeq_sch, group = "substrate", fraction = 0,
                                         weight = T, relative = T, plot = T)
sch_substrate_overlap

# Average overlap between soil and bark  
(alb_substrate_overlap$data$original.values[3] + 
    sch_substrate_overlap$data$original.values[3]) / 2

# Overlap between tree species in the regions. 
# Check the relative number of ASVs for the Swabian Alb.
alb_tree_overlap <- MicEco::ps_venn(physeq_alb, group = "dominant_tree", fraction = 0,
                                         weight = T, relative = T, plot = T)
alb_tree_overlap

# Check the relative number of ASVs for Schorfheide-Chorin.
sch_tree_overlap <- MicEco::ps_venn(physeq_sch, group = "dominant_tree", fraction = 0,
                                         weight = T, relative = T, plot = T)
sch_tree_overlap

# Average overlap between soil and bark  
(alb_tree_overlap$data$original.values[3] + 
    sch_tree_overlap$data$original.values[3]) / 2

# Average number of reads per region. 
asv_num_reads_alb <- mean(phyloseq::sample_sums(physeq_alb))
asv_num_reads_alb

asv_num_reads_sch <- mean(phyloseq::sample_sums(physeq_sch))
asv_num_reads_sch

# Number of ASVs per region-tree species combination. 
asv_num_alb_fagus <- phyloseq::ntaxa(physeq_alb_fagus)
asv_num_alb_fagus

asv_num_alb_picea <- phyloseq::ntaxa(physeq_alb_picea)
asv_num_alb_picea

asv_num_sch_fagus <- phyloseq::ntaxa(physeq_sch_fagus)
asv_num_sch_fagus

asv_num_sch_pinus <- phyloseq::ntaxa(physeq_sch_pinus)
asv_num_sch_pinus

# Overlap between tree species in the regions and substrate. 
# Check the relative number of ASVs for the Swabian Alb.
alb_fagus_soil_overlap <- MicEco::ps_venn(physeq_alb_fagus, group = "substrate", fraction = 0,
                                    weight = T, relative = T, plot = T)
alb_fagus_soil_overlap

# Check the relative number of ASVs for the Swabian Alb.
alb_picea_soil_overlap <- MicEco::ps_venn(physeq_alb_picea, group = "substrate", fraction = 0,
                                          weight = T, relative = T, plot = T)
alb_picea_soil_overlap

# Check the relative number of ASVs for Schorfheide-Chorin.
sch_tree_overlap <- MicEco::ps_venn(physeq_sch, group = "dominant_tree", fraction = 0,
                                    weight = T, relative = T, plot = T)
sch_tree_overlap

# Average overlap between soil and bark  
(alb_tree_overlap$data$original.values[3] + 
    sch_tree_overlap$data$original.values[3]) / 2


# Number of ASVs entering the networks per region. 
asv_num_alb_network <- phyloseq::ntaxa(phy_abundfilt_A)
asv_num_alb_network

asv_num_sch_network <- phyloseq::ntaxa(phy_abundfilt_S)
asv_num_sch_network

# Percentage of Lecanoromycetes (of total #)  
asv_num_soil_lecanoromycetes <- phyloseq::ntaxa(
  phyloseq::subset_taxa(physeq_soil, Class == "Lecanoromycetes"))
asv_num_soil_lecanoromycetes
asv_num_soil_lecanoromycetes/asv_num_soil * 100

asv_num_bark_lecanoromycetes <- phyloseq::ntaxa(
  phyloseq::subset_taxa(physeq_bark, Class == "Lecanoromycetes"))
asv_num_bark_lecanoromycetes
asv_num_bark_lecanoromycetes/asv_num_bark * 100

# Percentage of Lecanoromycetes (relative abundance)
phy_class_bark <-  physeq_bark %>%
  microbiome::aggregate_taxa(level = "Class")  

# Calculate the total number of reads
total_reads_bark <- sum(rowSums(otu_table(phy_class_bark)))
total_reads_bark

# Calculate the total number of reads for the lecanoromycetes  
total_reads_lecanoromycetes_bark <- rowSums(otu_table(phyloseq::subset_taxa(phy_class_bark, Class == "Lecanoromycetes")))
total_reads_lecanoromycetes_bark

# Divide the two total to get the relative abundance
rel_abund_lecanoromycetes_bark <- total_reads_lecanoromycetes_bark / total_reads_bark * 100
rel_abund_lecanoromycetes_bark

phy_class_soil <-  physeq_soil %>%
  microbiome::aggregate_taxa(level = "Class")  

# Calculate the total number of reads
total_reads_soil <- sum(rowSums(otu_table(phy_class_soil)))
total_reads_soil

# Calculate the total number of reads for the lecanoromycetes  
total_reads_lecanoromycetes_soil <- rowSums(otu_table(phyloseq::subset_taxa(phy_class_soil, Class == "Lecanoromycetes")))
total_reads_lecanoromycetes_soil

# Divide the two total to get the relative abundance
rel_abund_lecanoromycetes_soil <- total_reads_lecanoromycetes_soil / total_reads_soil * 100
rel_abund_lecanoromycetes_soil

# Overlap between trees in the full substrate dataset. 
soil_tree_overlap <- MicEco::ps_venn(physeq_soil, group = "dominant_tree", fraction = 0,
                                    weight = T, relative = T, plot = T)
soil_tree_overlap

bark_tree_overlap <- MicEco::ps_venn(physeq_bark, group = "dominant_tree", fraction = 0,
                                     weight = T, relative = T, plot = T)
bark_tree_overlap

#################################################################
##                          Section 11                         ##
##                    Plotting alpha diversity                 ##
#################################################################

################Substrate differences#######################

      ################ASV Counts#######################

# Create a joint table of ASV numbers for the substrates and differences. 

alb_substrate_list <- MicEco::ps_venn(physeq_alb, group = "substrate",
                fraction = 0, weight = F, relative = F, plot = F)

asv_num_alb_soil <- phyloseq::ntaxa(physeq_alb_soil)
asv_num_alb_bark <- phyloseq::ntaxa(physeq_alb_bark)

alb_substrate_df <- base::data.frame(exploratory = "Swabian Alb",
                                     value = c(length(alb_substrate_list$soil), 
                                               length(alb_substrate_list$bark),
                                               length(alb_substrate_list$bark__soil),
                                               asv_num_alb_soil,
                                               asv_num_alb_bark),
                                     variable = c("soil only",
                                                  "bark only",
                                                  "shared",
                                                  "soil total",
                                                  "bark total"),
                                     stack = c(rep("yes", 3),
                                               rep("no1", 1),
                                               rep("no2", 1)))

sch_substrate_list <- MicEco::ps_venn(physeq_sch, group = "substrate",
                                      fraction = 0, weight = F, relative = F, plot = F)

asv_num_sch_soil <- phyloseq::ntaxa(physeq_sch_soil)
asv_num_sch_bark <- phyloseq::ntaxa(physeq_sch_bark)

sch_substrate_df <- base::data.frame(exploratory = "Schorfheide-Chorin",
                                     value = c(length(sch_substrate_list$soil), 
                                                  length(sch_substrate_list$bark),
                                                  length(sch_substrate_list$bark__soil),
                                                  asv_num_sch_soil,
                                                  asv_num_sch_bark),
                                     variable = c("soil only",
                                                  "bark only",
                                                  "shared",
                                                  "soil total",
                                                  "bark total"),
                                     stack = c(rep("yes", 3),
                                               rep("no1", 1),
                                               rep("no2", 1)))

# Combine the regional dataframes. 
full_substrate_df <- rbind(sch_substrate_df, alb_substrate_df)

full_substrate_df$exploratory <- base::factor(full_substrate_df$exploratory, 
                                              levels = c("Swabian Alb",
                                                         "Schorfheide-Chorin"))

full_substrate_df$variable <- base::factor(full_substrate_df$variable,
                                           levels = c("shared",
                                                      "bark only",
                                                      "soil only",
                                                      "bark total",
                                                      "soil total"))
# Create barplot.
asv_num_substrate_barplot <- ggpubr::ggbarplot(full_substrate_df, x = "stack", y = "value", 
                  fill = "variable", color = "variable", width = 0.9,
                  palette = c("#009E73",  "#E69F00", "#56B4E9",
                              ggplot2::alpha("#E69F00", 0.5),
                              ggplot2::alpha("#56B4E9", 0.5))) +
  ggplot2::facet_grid("exploratory",
                      space="free", scales="free", switch = "y") +
  ggplot2::scale_x_discrete(position = "top", limits = c("no2", "no1", "yes")) +
  ggplot2::coord_flip() + 
  ggplot2::ylab("Number of ASVs") + 
  ggplot2::labs(subtitle = "(A)") + 
  ggplot2::theme(text = element_text(size = 15),
                 legend.title = element_blank(),
                 axis.title.y = element_blank(),
                 axis.line.y = element_blank(), 
                 axis.ticks.y = element_blank(), 
                 axis.text.y = element_blank(), 
                 legend.position = "bottom") 
asv_num_substrate_barplot

      ################Relative abundance#######################

alb_bark_only_otu <- subset(otu_table(physeq_alb), 
                    rownames(otu_table(physeq_alb)) %in% alb_substrate_list$bark)

alb_soil_only_otu <- base::subset(otu_table(physeq_alb), 
                            rownames(otu_table(physeq_alb)) %in% alb_substrate_list$soil)

alb_shared_otu <- subset(otu_table(physeq_alb), 
                            rownames(otu_table(physeq_alb)) %in% alb_substrate_list$bark__soil)



rel_abund_bark_only_alb <- round(sum(colSums(alb_bark_only_otu)) /
                                   sum(taxa_sums(physeq_alb)) * 100, 2)
rel_abund_bark_only_alb

rel_abund_soil_only_alb <- round(sum(colSums(alb_soil_only_otu)) /
                                   sum(taxa_sums(physeq_alb)) * 100, 2)
rel_abund_soil_only_alb

rel_abund_shared_alb <- round(sum(colSums(alb_shared_otu)) /
                                sum(taxa_sums(physeq_alb)) * 100, 2)
rel_abund_shared_alb

rel_abund_bark_total_alb <- round(sum(taxa_sums(physeq_alb_bark)) / sum(taxa_sums(physeq_alb)) * 100, 2)
rel_abund_bark_total_alb

rel_abund_soil_total_alb <- round(sum(taxa_sums(physeq_alb_soil)) / sum(taxa_sums(physeq_alb)) * 100, 2)
rel_abund_soil_total_alb

rel_abund_alb_substrate_df <- base::data.frame(exploratory = "Swabian Alb",
                                     value = c(rel_abund_soil_only_alb, 
                                               rel_abund_bark_only_alb,
                                               rel_abund_shared_alb,
                                               rel_abund_soil_total_alb,
                                               rel_abund_bark_total_alb),
                                     variable = c("soil only",
                                                  "bark only",
                                                  "shared",
                                                  "soil total",
                                                  "bark total"),
                                     stack = c(rep("yes", 3),
                                               rep("no", 2)))

sch_bark_only_otu <- subset(otu_table(physeq_sch), 
                            rownames(otu_table(physeq_sch)) %in% sch_substrate_list$bark)

sch_soil_only_otu <- base::subset(otu_table(physeq_sch), 
                                  rownames(otu_table(physeq_sch)) %in% sch_substrate_list$soil)

sch_shared_otu <- subset(otu_table(physeq_sch), 
                         rownames(otu_table(physeq_sch)) %in% sch_substrate_list$bark__soil)



rel_abund_bark_only_sch <- round(sum(colSums(sch_bark_only_otu)) /
                                   sum(taxa_sums(physeq_sch)) * 100, 2)
rel_abund_bark_only_sch

rel_abund_soil_only_sch <- round(sum(colSums(sch_soil_only_otu)) /
                                   sum(taxa_sums(physeq_sch)) * 100, 2)
rel_abund_soil_only_sch

rel_abund_shared_sch <- round(sum(colSums(sch_shared_otu)) /
                                sum(taxa_sums(physeq_sch)) * 100, 2)
rel_abund_shared_sch

rel_abund_bark_total_sch <- round(sum(taxa_sums(physeq_sch_bark)) / sum(taxa_sums(physeq_sch)) * 100, 2)
rel_abund_bark_total_sch

rel_abund_soil_total_sch <- round(sum(taxa_sums(physeq_sch_soil)) / sum(taxa_sums(physeq_sch)) * 100, 2)
rel_abund_soil_total_sch

rel_abund_sch_substrate_df <- base::data.frame(exploratory = "Schorfheide-Chorin",
                                               value = c(rel_abund_soil_only_sch, 
                                                         rel_abund_bark_only_sch,
                                                         rel_abund_shared_sch,
                                                         rel_abund_soil_total_sch,
                                                         rel_abund_bark_total_sch),
                                               variable = c("soil only",
                                                            "bark only",
                                                            "shared",
                                                            "soil total",
                                                            "bark total"),
                                               stack = c(rep("yes", 3),
                                                         rep("no", 2)))

# Combine the regional dataframes. 
rel_abund_full_substrate_df <- rbind(rel_abund_sch_substrate_df,
                                     rel_abund_alb_substrate_df)

rel_abund_full_substrate_df$exploratory <- base::factor(rel_abund_full_substrate_df$exploratory, 
                                              levels = c("Swabian Alb",
                                                         "Schorfheide-Chorin"))

rel_abund_full_substrate_df$variable <- base::factor(rel_abund_full_substrate_df$variable,
                                           levels = c("shared",
                                                      "bark only",
                                                      "soil only",
                                                      "bark total",
                                                      "soil total"))
# Create barplot.
rel_abund_substrate_barplot <- ggpubr::ggbarplot(rel_abund_full_substrate_df, x = "stack", y = "value", 
                  fill = "variable", color = "variable", width = 0.9,
                  palette = c("#009E73",
                              "#E69F00",
                              "#56B4E9",
                              ggplot2::alpha("#E69F00", 0.5),
                              ggplot2::alpha("#56B4E9", 0.5))) +
  ggplot2::facet_grid("exploratory",
                      space="free", scales="free", switch = "y") +
  ggplot2::scale_x_discrete(position = "top", limits = c("no", "yes")) +
  ggplot2::coord_flip() + 
  ggplot2::ylab("Relative abundance (%)") + 
  ggplot2::labs(subtitle = "(B)") + 
  ggplot2::theme(text = element_text(size = 15),
                 legend.title = element_blank(),
                 axis.title.y = element_blank(),
                 axis.line.y = element_blank(), 
                 axis.ticks.y = element_blank(), 
                 axis.text.y = element_blank(), 
                 legend.position = "bottom") 
rel_abund_substrate_barplot

################Tree differences#######################

################ASV Counts#######################

# Create a joint table of ASV numbers for the substrates and differences. 

alb_tree_list_bark <- MicEco::ps_venn(physeq_alb_bark, group = "dominant_tree",
                                      fraction = 0, weight = F, relative = F, plot = F)

asv_num_alb_fagus_bark <- phyloseq::ntaxa(physeq_alb_bark_fagus)
asv_num_alb_picea_bark <- phyloseq::ntaxa(physeq_alb_bark_picea)

alb_tree_df_bark <- base::data.frame(exploratory = "Swabian Alb \n bark",
                                     value = c(length(alb_tree_list_bark$Fagus_sylvatica), 
                                               length(alb_tree_list_bark$Picea_abies),
                                               length(alb_tree_list_bark$Fagus_sylvatica__Picea_abies),
                                               asv_num_alb_fagus_bark,
                                               asv_num_alb_picea_bark),
                                     variable = c("Fagus only",
                                                  "Picea only",
                                                  "shared",
                                                  "Fagus total",
                                                  "Picea total"),
                                     stack = c(rep("yes", 3),
                                               rep("no1", 1),
                                               rep("no2", 1)))


alb_tree_list_soil <- MicEco::ps_venn(physeq_alb_soil, group = "dominant_tree",
                                      fraction = 0, weight = F, relative = F, plot = F)

asv_num_alb_fagus_soil <- phyloseq::ntaxa(physeq_alb_soil_fagus)
asv_num_alb_picea_soil <- phyloseq::ntaxa(physeq_alb_soil_picea)

alb_tree_df_soil <- base::data.frame(exploratory = "Swabian Alb \n soil",
                                     value = c(length(alb_tree_list_soil$Fagus_sylvatica), 
                                               length(alb_tree_list_soil$Picea_abies),
                                               length(alb_tree_list_soil$Fagus_sylvatica__Picea_abies),
                                               asv_num_alb_fagus_soil,
                                               asv_num_alb_picea_soil),
                                     variable = c("Fagus only",
                                                  "Picea only",
                                                  "shared",
                                                  "Fagus total",
                                                  "Picea total"),
                                     stack = c(rep("yes", 3),
                                               rep("no1", 1),
                                               rep("no2", 1)))

sch_tree_list_bark <- MicEco::ps_venn(physeq_sch_bark, group = "dominant_tree",
                                      fraction = 0, weight = F, relative = F, plot = F)

asv_num_sch_fagus_bark <- phyloseq::ntaxa(physeq_sch_bark_fagus)
asv_num_sch_pinus_bark <- phyloseq::ntaxa(physeq_sch_bark_pinus)

sch_tree_df_bark <- base::data.frame(exploratory = "Schorfheide- Chorin \n bark",
                                     value = c(length(sch_tree_list_bark$Fagus_sylvatica), 
                                               length(sch_tree_list_bark$Pinus_sylvestris),
                                               length(sch_tree_list_bark$Fagus_sylvatica__Picea_abies),
                                               asv_num_sch_fagus_bark,
                                               asv_num_sch_pinus_bark),
                                     variable = c("Fagus only",
                                                  "Pinus only",
                                                  "shared",
                                                  "Fagus total",
                                                  "Pinus total"),
                                     stack = c(rep("yes", 3),
                                               rep("no1", 1),
                                               rep("no2", 1)))


sch_tree_list_soil <- MicEco::ps_venn(physeq_sch_soil, group = "dominant_tree",
                                      fraction = 0, weight = F, relative = F, plot = F)

asv_num_sch_fagus_soil <- phyloseq::ntaxa(physeq_sch_soil_fagus)
asv_num_sch_pinus_soil <- phyloseq::ntaxa(physeq_sch_soil_pinus)

sch_tree_df_soil <- base::data.frame(exploratory = "Schorfheide- Chorin \n soil",
                                     value = c(length(sch_tree_list_soil$Fagus_sylvatica), 
                                               length(sch_tree_list_soil$Pinus_sylvestris),
                                               length(sch_tree_list_soil$Fagus_sylvatica__Picea_abies),
                                               asv_num_sch_fagus_soil,
                                               asv_num_sch_pinus_soil),
                                     variable = c("Fagus only",
                                                  "Pinus only",
                                                  "shared",
                                                  "Fagus total",
                                                  "Pinus total"),
                                     stack = c(rep("yes", 3),
                                               rep("no1", 1),
                                               rep("no2", 1)))

# Combine the regional dataframes. 
full_tree_df <- rbind(sch_tree_df_soil, sch_tree_df_bark, 
                      alb_tree_df_soil, alb_tree_df_bark)

full_tree_df$exploratory <- base::factor(full_tree_df$exploratory, 
                                              levels = c("Swabian Alb \n bark",
                                                         "Swabian Alb \n soil",
                                                         "Schorfheide- Chorin \n bark",
                                                         "Schorfheide- Chorin \n soil"))

full_tree_df$variable <- base::factor(full_tree_df$variable,
                                           levels = c("shared",
                                                      "Fagus only",
                                                      "Picea only",
                                                      "Fagus total",
                                                      "Picea total",
                                                      "Pinus only",
                                                      "Pinus total"))
# Create barplot.
asv_num_tree_barplot <- ggpubr::ggbarplot(full_tree_df, x = "stack", y = "value", 
                  fill = "variable", color = "variable", width = 0.9,
                  palette = c("#009E73",
                              "#E69F00",
                              "#56B4E9",
                              ggplot2::alpha("#E69F00", 0.5),
                              ggplot2::alpha("#56B4E9", 0.5),
                              "#26547C",
                              ggplot2::alpha("#26547C", 0.5))) +
  ggplot2::facet_grid("exploratory",
                      space="free", scales="free", switch = "y",
                      labeller = label_wrap_gen(width = 15)) +
  ggplot2::scale_x_discrete(position = "top", limits = c("no2", "no1", "yes")) +
  ggplot2::coord_flip() + 
  ggplot2::ylab("Number of ASVs") + 
  ggplot2::labs(subtitle = "(A)") + 
  guides(colour = guide_legend(nrow = 1)) + 
  ggplot2::theme(text = element_text(size = 15),
                 legend.title = element_blank(),
                 axis.title.y = element_blank(),
                 axis.line.y = element_blank(), 
                 axis.ticks.y = element_blank(), 
                 axis.text.y = element_blank(), 
                 legend.position = "bottom") 
asv_num_tree_barplot

################Relative abundance#######################

alb_bark_fagus_only_otu <- subset(otu_table(physeq_alb_bark), 
                            rownames(otu_table(physeq_alb_bark)) %in% alb_tree_list_bark$Fagus_sylvatica)

alb_bark_picea_only_otu <- base::subset(otu_table(physeq_alb_bark), 
                                  rownames(otu_table(physeq_alb_bark)) %in% alb_tree_list_bark$Picea_abies)

alb_bark_shared_otu <- subset(otu_table(physeq_alb_bark), 
                         rownames(otu_table(physeq_alb_bark)) %in% alb_tree_list_bark$Fagus_sylvatica__Picea_abies)

rel_abund_bark_fagus_only_alb <- round(sum(colSums(alb_bark_fagus_only_otu)) /
                                   sum(taxa_sums(physeq_alb_bark)) * 100, 2)
rel_abund_bark_fagus_only_alb

rel_abund_bark_picea_only_alb <- round(sum(colSums(alb_bark_picea_only_otu)) /
                                   sum(taxa_sums(physeq_alb_bark)) * 100, 2)
rel_abund_bark_picea_only_alb

rel_abund_bark_shared_alb <- round(sum(colSums(alb_bark_shared_otu)) /
                                sum(taxa_sums(physeq_alb_bark)) * 100, 2)
rel_abund_bark_shared_alb

rel_abund_fagus_total_alb_bark <- round(sum(taxa_sums(physeq_alb_bark_fagus)) /
                                     sum(taxa_sums(physeq_alb_bark)) * 100, 2)
rel_abund_fagus_total_alb_bark

rel_abund_picea_total_alb_bark <- round(sum(taxa_sums(physeq_alb_bark_picea)) /
                                     sum(taxa_sums(physeq_alb_bark)) * 100, 2)
rel_abund_picea_total_alb_bark

rel_abund_alb_tree_df_bark <- base::data.frame(exploratory = "Swabian Alb bark",
                                               value = c(rel_abund_bark_fagus_only_alb, 
                                                         rel_abund_bark_picea_only_alb,
                                                         rel_abund_bark_shared_alb,
                                                         rel_abund_fagus_total_alb_bark,
                                                         rel_abund_picea_total_alb_bark),
                                               variable = c("Fagus only",
                                                            "Picea only",
                                                            "shared",
                                                            "Fagus total",
                                                            "Picea total"),
                                               stack = c(rep("yes", 3),
                                                         rep("no", 2)))

alb_soil_fagus_only_otu <- subset(otu_table(physeq_alb_soil), 
                                  rownames(otu_table(physeq_alb_soil)) %in% alb_tree_list_soil$Fagus_sylvatica)

alb_soil_picea_only_otu <- base::subset(otu_table(physeq_alb_soil), 
                                        rownames(otu_table(physeq_alb_soil)) %in% alb_tree_list_soil$Picea_abies)

alb_soil_shared_otu <- subset(otu_table(physeq_alb_soil), 
                              rownames(otu_table(physeq_alb_soil)) %in% alb_tree_list_soil$Fagus_sylvatica__Picea_abies)

rel_abund_soil_fagus_only_alb <- round(sum(colSums(alb_soil_fagus_only_otu)) /
                                         sum(taxa_sums(physeq_alb_soil)) * 100, 2)
rel_abund_soil_fagus_only_alb

rel_abund_soil_picea_only_alb <- round(sum(colSums(alb_soil_picea_only_otu)) /
                                         sum(taxa_sums(physeq_alb_soil)) * 100, 2)
rel_abund_soil_picea_only_alb

rel_abund_soil_shared_alb <- round(sum(colSums(alb_soil_shared_otu)) /
                                     sum(taxa_sums(physeq_alb_soil)) * 100, 2)
rel_abund_soil_shared_alb

rel_abund_fagus_total_alb_soil <- round(sum(taxa_sums(physeq_alb_soil_fagus)) /
                                          sum(taxa_sums(physeq_alb_soil)) * 100, 2)
rel_abund_fagus_total_alb_soil

rel_abund_picea_total_alb_soil <- round(sum(taxa_sums(physeq_alb_soil_picea)) /
                                          sum(taxa_sums(physeq_alb_soil)) * 100, 2)
rel_abund_picea_total_alb_soil

rel_abund_alb_tree_df_soil <- base::data.frame(exploratory = "Swabian Alb soil",
                                               value = c(rel_abund_soil_fagus_only_alb, 
                                                         rel_abund_soil_picea_only_alb,
                                                         rel_abund_soil_shared_alb,
                                                         rel_abund_fagus_total_alb_soil,
                                                         rel_abund_picea_total_alb_soil),
                                               variable = c("Fagus only",
                                                            "Picea only",
                                                            "shared",
                                                            "Fagus total",
                                                            "Picea total"),
                                               stack = c(rep("yes", 3),
                                                         rep("no", 2)))

sch_bark_fagus_only_otu <- subset(otu_table(physeq_sch_bark), 
                                  rownames(otu_table(physeq_sch_bark)) %in% sch_tree_list_bark$Fagus_sylvatica)

sch_bark_pinus_only_otu <- base::subset(otu_table(physeq_sch_bark), 
                                        rownames(otu_table(physeq_sch_bark)) %in% sch_tree_list_bark$Pinus_sylvestris)

sch_bark_shared_otu <- subset(otu_table(physeq_sch_bark), 
                              rownames(otu_table(physeq_sch_bark)) %in% sch_tree_list_bark$Fagus_sylvatica__Pinus_sylvestris)

rel_abund_bark_fagus_only_sch <- round(sum(colSums(sch_bark_fagus_only_otu)) /
                                         sum(taxa_sums(physeq_sch_bark)) * 100, 2)
rel_abund_bark_fagus_only_sch

rel_abund_bark_pinus_only_sch <- round(sum(colSums(sch_bark_pinus_only_otu)) /
                                         sum(taxa_sums(physeq_sch_bark)) * 100, 2)
rel_abund_bark_pinus_only_sch

rel_abund_bark_shared_sch <- round(sum(colSums(sch_bark_shared_otu)) /
                                     sum(taxa_sums(physeq_sch_bark)) * 100, 2)
rel_abund_bark_shared_sch

rel_abund_fagus_total_sch_bark <- round(sum(taxa_sums(physeq_sch_bark_fagus)) /
                                          sum(taxa_sums(physeq_sch_bark)) * 100, 2)
rel_abund_fagus_total_sch_bark

rel_abund_pinus_total_sch_bark <- round(sum(taxa_sums(physeq_sch_bark_pinus)) /
                                          sum(taxa_sums(physeq_sch_bark)) * 100, 2)
rel_abund_pinus_total_sch_bark

rel_abund_sch_tree_df_bark <- base::data.frame(exploratory = "Schorfheide- Chorin bark",
                                               value = c(rel_abund_bark_fagus_only_sch, 
                                                         rel_abund_bark_pinus_only_sch,
                                                         rel_abund_bark_shared_sch,
                                                         rel_abund_fagus_total_sch_bark,
                                                         rel_abund_pinus_total_sch_bark),
                                               variable = c("Fagus only",
                                                            "Pinus only",
                                                            "shared",
                                                            "Fagus total",
                                                            "Pinus total"),
                                               stack = c(rep("yes", 3),
                                                         rep("no", 2)))

sch_soil_fagus_only_otu <- subset(otu_table(physeq_sch_soil), 
                                  rownames(otu_table(physeq_sch_soil)) %in% sch_tree_list_soil$Fagus_sylvatica)

sch_soil_pinus_only_otu <- base::subset(otu_table(physeq_sch_soil), 
                                        rownames(otu_table(physeq_sch_soil)) %in% sch_tree_list_soil$Pinus_sylvestris)

sch_soil_shared_otu <- subset(otu_table(physeq_sch_soil), 
                              rownames(otu_table(physeq_sch_soil)) %in% sch_tree_list_soil$Fagus_sylvatica__Pinus_sylvestris)

rel_abund_soil_fagus_only_sch <- round(sum(colSums(sch_soil_fagus_only_otu)) /
                                         sum(taxa_sums(physeq_sch_soil)) * 100, 2)
rel_abund_soil_fagus_only_sch

rel_abund_soil_pinus_only_sch <- round(sum(colSums(sch_soil_pinus_only_otu)) /
                                         sum(taxa_sums(physeq_sch_soil)) * 100, 2)
rel_abund_soil_pinus_only_sch

rel_abund_soil_shared_sch <- round(sum(colSums(sch_soil_shared_otu)) /
                                     sum(taxa_sums(physeq_sch_soil)) * 100, 2)
rel_abund_soil_shared_sch

rel_abund_fagus_total_sch_soil <- round(sum(taxa_sums(physeq_sch_soil_fagus)) /
                                          sum(taxa_sums(physeq_sch_soil)) * 100, 2)
rel_abund_fagus_total_sch_soil

rel_abund_pinus_total_sch_soil <- round(sum(taxa_sums(physeq_sch_soil_pinus)) /
                                          sum(taxa_sums(physeq_sch_soil)) * 100, 2)
rel_abund_pinus_total_sch_soil

rel_abund_sch_tree_df_soil <- base::data.frame(exploratory = "Schorfheide- Chorin soil",
                                               value = c(rel_abund_soil_fagus_only_sch, 
                                                         rel_abund_soil_pinus_only_sch,
                                                         rel_abund_soil_shared_sch,
                                                         rel_abund_fagus_total_sch_soil,
                                                         rel_abund_pinus_total_sch_soil),
                                               variable = c("Fagus only",
                                                            "Pinus only",
                                                            "shared",
                                                            "Fagus total",
                                                            "Pinus total"),
                                               stack = c(rep("yes", 3),
                                                         rep("no", 2)))


# Combine the regional dataframes. 
rel_abund_full_tree_df <- rbind(rel_abund_sch_tree_df_soil,
                                rel_abund_sch_tree_df_bark,
                                rel_abund_alb_tree_df_soil,
                                rel_abund_alb_tree_df_bark)

rel_abund_full_tree_df$exploratory <- base::factor(rel_abund_full_tree_df$exploratory, 
                                                   levels = c("Swabian Alb bark",
                                                              "Swabian Alb soil",
                                                              "Schorfheide- Chorin bark",
                                                              "Schorfheide- Chorin soil"))

rel_abund_full_tree_df$variable <- base::factor(rel_abund_full_tree_df$variable,
                                      levels = c("shared",
                                                 "Fagus only",
                                                 "Picea only",
                                                 "Fagus total",
                                                 "Picea total",
                                                 "Pinus only",
                                                 "Pinus total"))
# Create barplot.
rel_abund_tree_barplot <- ggpubr::ggbarplot(rel_abund_full_tree_df, x = "stack", y = "value", 
                  fill = "variable", color = "variable", width = 0.9,
                  palette = c("#009E73",
                              "#E69F00",
                              "#56B4E9",
                              ggplot2::alpha("#E69F00", 0.5),
                              ggplot2::alpha("#56B4E9", 0.5),
                              "#26547C",
                              ggplot2::alpha("#26547C", 0.5))) +
  ggplot2::facet_grid("exploratory",
                      space="free", scales="free", switch = "y",
                      labeller = label_wrap_gen(width = 15)) +
  ggplot2::scale_x_discrete(position = "top", limits = c("no", "yes")) +
  ggplot2::coord_flip() + 
  ggplot2::ylab("Relative abundance (%)") + 
  ggplot2::labs(subtitle = "(B)") + 
  guides(colour = guide_legend(nrow = 1)) + 
  ggplot2::theme(text = element_text(size = 15),
                 legend.title = element_blank(),
                 axis.title.y = element_blank(),
                 axis.line.y = element_blank(), 
                 axis.ticks.y = element_blank(), 
                 axis.text.y = element_blank(), 
                 legend.position = "bottom") 
rel_abund_tree_barplot


# Create the two final figures.

substrate_differences <- ggpubr::ggarrange(asv_num_substrate_barplot, rel_abund_substrate_barplot,
                                      ncol = 2, nrow = 1,
                                      legend = "bottom", common.legend = T)
substrate_differences

ggsave('substrate_differences.tiff', device = 'tiff',
       substrate_differences, width = 400, height = 240,
       units = 'mm', dpi = 300)

tree_differences <- ggpubr::ggarrange(asv_num_tree_barplot, rel_abund_tree_barplot,
                                           ncol = 2, nrow = 1,
                                           legend = "bottom", common.legend = T)
tree_differences

ggsave('tree_differences.tiff', device = 'tiff',
       tree_differences, width = 400, height = 240,
       units = 'mm', dpi = 300)
