############### Phylogeny Construction - Egyptian Seed Traits ######################
# This code is used to construct phylogenetic trees for seed traits in Egyptian plants
# Code developed by Garland Xie
####################################################################################

# Install packages and libraries (if you don't have them!) --------------------------

# Load up the easypackages 
if(!require(easypackages)){
  install.packages("easypackages")
  library(easypackages)
}

# check to see if you have these packages. If not, install them
packages("ape", "here", "brranching", "ggtree", "Hmisc", "magrittr", "phytools")

# use bioconductor to install ggtree
source("https://bioconductor.org/biocLite.R")

# Load libraries ---------------------------------------------------------------------
# ape: loading newick files
# here: constructing relative file paths
# brranching: constructing phylomatic trees
# magittr: pipe function
libraries("ape", "here", "brranching", "ggtree", "Hmisc", "phytools")

# Load datasets ----------------------------------------------------------------------
seed_clean_df <- readRDS(here("Raw Datasets/Species", "seed_clean"))

# Phylogeny construction -------------------------------------------------------------

# Used Zanne et al. 2014 as a backbone phylogeny; it has nearly every extant plant family
# Used APG to constrain orders and families (APG 2009)
# Taxonomic names were updated according to Qian et al. 2015 (see link below)

link = "https://raw.githubusercontent.com/jinyizju/S.PhyloMaker/master/PhytoPhylo"

# 103 species directly match with zanne2014 (65%)
# 55 species were placed on either genus or family level (35%)
# No species were excluded!
tree_phy <- phylomatic(seed_clean_df$spp, treeuri = link, get = 'POST')

# Weird: phylomatic functon gives out a unique phylo class object
# Doesn't work well with ggtree plotting functions
# Rewrite the tree as a newick file 
write.tree(tree_phy, here("Raw Datasets/Phylogeny", "phylo.new"))
tree_phy <- read.tree(here("Raw Datasets/Phylogeny", "phylo.new"))

# Adjust string format
tree_phy$tip.label <- capitalize(tree_phy$tip.label)

# Update string formats of species in seed_clean_df to match with tree_phy
seed_clean_df$spp <- gsub(" ", "_", seed_clean_df$spp)

# Check to see if taxon names in trait dataset is consistent with phylogeny
if(!all(seed_clean_df$spp %in% tree_phy$tip.label)) {
  
  # Re-update taxa in seed_clean_df
  seed_clean_df$spp[seed_clean_df$spp == "Cyperus_rotundus_var._rotundus"] <- "Cyperus_rotundus.var.rotundus"
  seed_clean_df$spp[seed_clean_df$spp == "Salsola_imbricata.supsp._imbricata"] <- "Salsola_imbricata.supsp.imbricata"
  
  # Re-update taxa in tree_phy
  tree_phy$tip.label[tree_phy$tip.label == "Cyperus_rotundusvar.rotundus"] <- "Cyperus_rotundus.var.rotundus"
}

# Check to see if tree is ultrametric
plot.phylo(tree_phy, show.tip.label = F) # looks ultrametric

# Multiple polytomous branch lengths in current tree
# Resolve them based on order of tree for certain phylo signal testss
if(!is.ultrametric(tree_phy)) {
  tree_phy <- force.ultrametric(tree_phy)
  tree_phy <- multi2di(tree_phy, random = F)
}

# adjust string format to match trait dataset with phylogeny
seed_clean_df <- gsub(" ", "_", seed_clean_df)

# Save the phylogenetic tree as an R object ---------------------------------------------------------------------------
saveRDS(tree_phy, here("R Objects", "tree_phy"))
saveRDS(seed_clean_df, here("Raw Datasets", "seed_clean"))

# Done! Move to 2-Phylo_Signals_Reham.R
      
    
tree_phy$tip.label[order(tree_phy$tip.label)] == seed_clean_df$spp[order(seed_clean_df$spp)]
