##-------------------------------------------------------------------------
## Replace names of sequences with complete names for 
## BEAST analyses (GNUMBER/type/date)
## with type=0 for DS strains and type=1 for MDR strains
## 2021-06-21 Etthel Windels
##-------------------------------------------------------------------------



# Load libraries ----------------------------------------------------------

library(dplyr)
library(stringr)
library(tibble)
library(seqinr)


# Read files --------------------------------------------------------------

metadata <- read.csv('path_to_data', sep=',', dec='.', header=T)
alignment <- seqinr::read.fasta('path_to_fasta', seqtype='DNA', forceDNAtolower=FALSE)
metadata$MDR <- ifelse(metadata$dr_category=='pan-susceptible', 0, 1)


# Change date format ------------------------------------------------------

for (j in 1:length(metadata$sample_collection_date)){
  if (is.na(metadata$sample_collection_date[j])){
    metadata$new_date[j] <- NA
  }
  else{
    metadata$new_date[j]=str_glue(str_split_fixed(metadata$sample_collection_date[j],'\\.',3)[,1],'-',str_split_fixed(metadata$sample_collection_date[j],'\\.',3)[,2],'-20',str_split_fixed(metadata$sample_collection_date[j],'\\.',3)[,3])
  }
}


# Adjust metadata ---------------------------------------------------------

metadata <- 
  metadata %>%
  add_column(new_id = paste(metadata$g_number, metadata$MDR, metadata$new_date, sep = "/")) %>% 
  filter(g_number %in% names(alignment))
alignment <- alignment[names(alignment) %in% metadata$g_number]


# Adjust alignment names --------------------------------------------------

gnumber <- names(alignment)
full_names <- metadata$new_id[match(gnumber, metadata$g_number)]
print(full_names)
if (any(is.na(full_names))) {
  stop("Not all strains have full names in metadata.")
}
names(alignment) <- full_names


# Save output file --------------------------------------------------------

seqinr::write.fasta(alignment,names=names(alignment),file.out="output_path")
