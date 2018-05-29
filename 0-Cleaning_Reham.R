############### Cleaning - Egyptian Seed Traits ####################################
# This code is used to analyze phylogenetic tests for seed traits in Egyptian plants
# Code developed by Garland Xie
####################################################################################

# Install packages (if you don't have them!) ----------------------------------

# Load up the easypackages 
if(!require(easypackages)){
  install.packages("easypackages")
  library(easypackages)
}

# check to see if you have these packages. If not, install them
packages("taxize", "here", "tidyverse", "readxl", "rebus", "magrittr")

# load libraries
# here: constructing relative file paths
# taxize: accessing TNRS database to update taxon names
# tidyverse: for making clean data
# readxl: importing excel documents
# rebus: human-readable custom pattern matching 
# magrittr: compound assignment operator 
libraries("here", "taxize", "tidyverse", "readxl", "rebus", "magrittr")

# Import data files -----------------------------------------------------------

# seed traits for each sampled species 
seed_raw <- read_xlsx(here("data/raw", "habitat_traits.xlsx"), sheet = "Habitats&traits")

# Cleaning -----------------------------------------------------------------------------------------

seed_clean <- seed_raw %>% 
  
  # rename columns for conciseness
  rename(spp = "Species",
         habitat = "Habitat",
         seed_mass = "seed mass(gm/100seed)", 
         seed_length = "seed length(mm", 
         seed_width = "seed width(mm)",
         seed_thick = "seed thickness(mm)", 
         seed_germ = "seed germination", 
         seed_rain = "Seed rain",
         seed_size = "Seed size") %>%
  
  # fix up typo mistakes (and possible taxonomic synonyms)
  # use case_when for multiple ifelse statements
  mutate(spp = case_when(spp == "Pseudognaphalium luteo-album" ~ "Helichrysum luteoalbum",
                         spp == "Veronica-anagalis aquatica" ~ "Veronica anagallis-aquatica",
                         spp == "Eclypta alba" ~ "Eclipta alba", 
                         spp == "Ipomea carnea" ~ "Ipomoea carnea",
                         spp == "Nympheae lotus" ~ "Nymphaea lotus", 
                         spp == "Rorippa Palustris" ~ "Rorippa palustris",
                         TRUE ~ spp)) 

# look for accepted taxonomic names using Taxonomic Name Resolution Database
tnrs_df <- tnrs(query = pull(seed_clean, spp))

## replace submitted names with accepted names 
# first, grab accepted names and submitted names 
tnrs_accepted <- pull(tnrs_df, acceptedname)
tnrs_submitted <- pull(tnrs_df, submittedname)

# them, create a logical vector that finds rows with any characters as TRUE; otherwise FALSE
contains_accepted_names <- str_detect(tnrs_accepted, pattern = wrd(1) %R% " " %R% wrd(1))

# then use ifelse statement: for rows that are TRUE, replace with the accepted names
# for rows that are FALSE, replace it with the submitted name 
# code works well except for 
seed_clean %<>% mutate(spp = ifelse(contains_accepted_names,
                                    tnrs_accepted, 
                                    tnrs_submitted)) 

# then just grab the genus and species as the first two string characters
seed_clean %<>% mutate(spp = word(spp, start = 1, end = 2, sep = fixed(" ")))

# Save data ---------------------------------------------------------------------------------------------------------------
saveRDS(seed_clean, here("data/working", "seed_clean.rds"))

# Done! Move to 1-Phylo-Construct_Reham.R

