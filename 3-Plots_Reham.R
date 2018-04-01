############### Plots - Egyptian Seed Traits ######################
# This code is used to create plots
# Code developed by Garland Xie
####################################################################

#### Install packages (if you don't have them!) ####

# Load up the easypackages 
if(!require(easypackages)){
  install.packages("easypackages")
  library(easypackages)
}

# check to see if you have these packages. If not, install them
packages("ggtree")

# load libraries
# here: constructing relative file paths
# taxize: accessing TNRS database
# stringr: for counting words in a string
# taxonlookup: accessing genus and family names 
libraries("ggtree")

### Plots - Phylogenetic Trees #####

# get genus-level info 
groupInfo <- split(tree_phy$tip.label, gsub("_\\w+", "", tree_phy$tip.label))

# update genus-level info in tree_phy
tree_phy <- groupOTU(tree_phy, groupInfo)

# plot tree
ggtree(tree_phy, aes(color = group),layout = "circular") + 
  geom_treescale(x = 350, y = 1, linesize = 2) +
  geom_nodepoint(color = "black", size = 1) + 
  geom_tiplab(size = 1, aes(angle=angle), col = "black") 