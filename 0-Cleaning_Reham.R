############### Cleaning - Egyptian Seed Traits ######################
# This code is used to analyze phylogenetic tests for seed traits in Egyptian plants
# Code developed by Garland Xie
###############################################################################

# Install packages (if you don't have them!) ----------------------------------------------

# Load up the easypackages 
if(!require(easypackages)){
  install.packages("easypackages")
  library(easypackages)
}

# check to see if you have these packages. If not, install them
packages("taxize", "here", "tidyverse")

# load libraries
# here: constructing relative file paths
# taxize: accessing TNRS database
# stringr: for counting words in a string
# taxonlookup: accessing genus and family names 
libraries("here", "taxize", "tidyverse")

# Import data files --------------------------------------------------------------------------------
# folder name: Egypyt_Seed_Conservation"
# use relative file paths

# here() starts at C:/Users/garla/Google Drive/Research/Research Projects/Egypyt_Seed_Conservation
here()

# seed traits for each sampled species
seed_raw_df <- read_csv(here("Raw Datasets/Species", "Reham_seed_data.csv"))

# Cleaning -----------------------------------------------------------------------------------------

# keep raw data, but create new object for clean dataset
seed_clean_df <- seed_raw_df

# Update column names to make them shorter
colnames(seed_clean_df) <- c("spp", 'mass', "len", "width", "thick", "germ", "rain", "size")

# Remove empty values from seed_clean_df$Species
seed_clean_df <- seed_clean_df[seed_clean_df$spp != "", ]

# Remove commas from seed_clean_df$spp
seed_clean_df$spp <- gsub(",", "", seed_clean_df$spp) 

# Re-update synonyms (or typo mistakes)
seed_clean_df$spp[seed_clean_df$spp == "Psuedogynaphallium lutea-album"] <- "Helichrysum luteoalbum"
seed_clean_df$spp[seed_clean_df$spp == "Veronica-anagalis aquatica"] <- "Veronica anagallis-aquatica" 
seed_clean_df$spp[seed_clean_df$spp == "Eclypta alba"] <- "Eclipta alba"
seed_clean_df$spp[seed_clean_df$spp == "Ipomea carnea"] <- "Ipomoea carnea"
seed_clean_df$spp[seed_clean_df$spp == "Nympheae lotus"] <- "Nymphaea lotus"
seed_clean_df$spp[seed_clean_df$spp == "Nympheae nouchali"] <- "Nymphaea nouchali"
seed_clean_df$spp[seed_clean_df$spp == "Rorippa Palustris"] <- "Rorippa palustris"

# update taxonomic names for seed_clean_df for phylo signal tests 
seed_clean_df <- seed_clean_df[!(seed_clean_df$spp == "Amaranthus hyperidus v. hyperidus"), ]

# look for accepted taxonomic names from multiple databases
tnrs_df <- tnrs(query = seed_clean_df$spp)

# find number of words per string in matched names
wc <- sapply(tnrs_df$matchedname, function(x) str_count(x, '\\w+'))

# replace genus-only characters with species-level submitted names
tnrs_df$matchedname[wc == 1] <- tnrs_df$submittedname[wc == 1]

# remove any missing values
seed_clean_df <- seed_clean_df[complete.cases(seed_clean_df), ]

# check to see if taxon names are consistent in each dataset; if not, update seed_clean_df
# double-check this function later on; it seems  like re-update taxa should be outside of bracket.. 
if (!all(seed_clean_df$spp %in% tnrs_df$matchedname)) {
  
  # re-organize dataframes based on a condition
  seed_clean_df <- seed_clean_df[order(seed_clean_df$spp), ]
  tnrs_df <- tnrs_df[order(tnrs_df$matchedname), ]
  
  # re-update taxa
  seed_clean_df$spp <- tnrs_df$matchedname
}


# Save data ---------------------------------------------------------------------------------------------------------------
saveRDS(seed_clean_df, here("R Objects", "seed_clean"))

# Done! Move to 1-Phylo-Construct_Reham.R


