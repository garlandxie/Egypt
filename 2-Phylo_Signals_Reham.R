############### Phylogenetic Signals - Egyptian Seed Traits ########################
# This code is used to analyze phylogenetic tests for seed traits in Egyptian plants
# Code developed by Garland Xie
####################################################################################

# Install packages -------------------------------------------------------------------

# Load up "easypackages" - allows verification of multiple libraries
if(!require(easypackages)){
  install.packages("easypackages")
  library(easypackages)
}

# check to see if you have these packages. If not, install them
packages("ape", "here", "phytools", "caper", "Hmisc", "geiger")
 
# Load libraries ---------------------------------------------------------------------
# ape: loading newick files
# here: constructing relative file paths
libraries("ape", "here", "phytools", "caper", "Hmisc", "geiger")

# Load datasets -----------------------------------------------------------------------

# Check to see if you're on the correct folder
here() 

# ultrametric phylogenetic tree 
tree_phy <- readRDS(here("R Objects", "tree_phy"))

# trait dataset 
seed_clean_df <- readRDS(here("Raw Datasets/Species", "seed_clean"))

# Phylogenetic signal tests - Continuous --------------------------------------------

#update rownames for seed_clean_df using species
rownames(seed_clean_df) <- seed_clean_df$spp

# create vectors for each trait
seed_mass <- seed_clean_df$mass
seed_len <- seed_clean_df$len
seed_width <- seed_clean_df$width
seed_thick <- seed_clean_df$thick
seed_rain <- seed_clean_df$rain
seed_size <- seed_clean_df$size

# add species as names
names(seed_mass) <- seed_clean_df$spp
names(seed_len) <- seed_clean_df$spp
names(seed_thick) <- seed_clean_df$spp
names(seed_width) <- seed_clean_df$spp
names(seed_rain) <- seed_clean_df$spp
names(seed_size) <- seed_clean_df$spp

# Blomberg's K (continuous traits)
phylosig(tree_phy, seed_mass, method = "K", nsim = 999, test = T)
phylosig(tree_phy, seed_len, method = "K", nsim = 999, test = T)
phylosig(tree_phy, seed_width, method = "K", nsim = 999, test = T)
phylosig(tree_phy, seed_thick, method = "K", nsim = 999, test = T)

# Pagel's lambda (continuous traits) - uses LR tests for significance
phylosig(tree_phy, seed_mass, method = "lambda", test = T)
phylosig(tree_phy, seed_len, method = "lambda", test = T)
phylosig(tree_phy, seed_width, method = "lambda", test = T)
phylosig(tree_phy, seed_thick, method = "lambda", test = T)

# D-Statistic (binary traits) ---------------------------------------------------------

# duplicated labels and nodes in phylogeny; remove node labels
tree_phy$node.label <- NULL
# calculate d-statistic
phylo.d(data = seed_clean_df, phy = tree_phy, binvar = germ, names.col = spp)

# Pagel's lambda (discrete traits) -----------------------------------------------------

phy_size1 <-fitDiscrete(tree_phy, seed_size, transform = "lambda")
star_tree<- rescale(tree_phy, "lambda", 0)
phy_size0 <-fitDiscrete(star_tree, seed_size, transform = "lambda")
LR_size <-2*(phy_size1$opt$lnL - phy_size0$opt$lnL) 
pchisq(LR_size, df = 1 ,lower.tail = F)

phy_rain1 <-fitDiscrete(tree_phy, seed_rain, transform = "lambda")
star_tree<- rescale(tree_phy, "lambda", 0)
phy_rain0 <-fitDiscrete(star_tree, seed_rain, transform = "lambda")
LR_rain <-2*(phy_rain1$opt$lnL - phy_rain0$opt$lnL) 
pchisq(LR_rain, d = 1 ,lower.tail = F)





  