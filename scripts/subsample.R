##-------------------------------------------------------------------------
## Take subsets of sequences for multitype birth-death analyses,
## considering 2 types: DS and MDR 
## (MDR also includes strains with compensatory mutations)
## 2021-05-21 Etthel Windels
##-------------------------------------------------------------------------



# Load libraries ----------------------------------------------------------

library(seqinr)
library(stringr)
library(dplyr)


# Read files --------------------------------------------------------------

metadata <- 
  read.csv('path_to_data', sep=',', dec='.', header=T) %>%
  .[!is.na(.$sample_collection_date),]  # only include samples for which sampling date is known
full_alignment <- 
  seqinr::read.fasta('path_to_fasta', seqtype='DNA', forceDNAtolower=FALSE) %>% # alignment of all sequences, with Gnumbers as sequence names
  .[names(.) %in% metadata$g_number] # only keep sequences with corresponding entry in filtered metadata table
metadata$MDR <- ifelse(metadata$dr_category=='pan-susceptible', 'DS', 'MDR')


# Function to generate random subset --------------------------------------

subsample <- function(fasta, metadata, n, lineage, seed){
  
  set.seed(seed)
  meta_lin <- metadata[metadata$mtb_lineage==lineage,]
  
  # Sample types with equal frequencies
  id <- c(sample(meta_lin$g_number[meta_lin$MDR=='DS'], round(n/2), replace=F),
          sample(meta_lin$g_number[meta_lin$MDR=='MDR'], round(n/2), replace=F))
  
  # Create FASTA file
  fasta_new <- fasta[names(fasta) %in% id]
  write.fasta(fasta_new,names=names(fasta_new),file.out=str_glue(lineage,'_',n,'_',seed,'.fasta'))
}


# Generate 3 subsets per lineage ------------------------------------------

subsample(fasta=full_alignment, metadata=metadata, n=200, lineage='L2', seed=123) 
subsample(fasta=full_alignment, metadata=metadata, n=200, lineage='L4', seed=123) 
subsample(fasta=full_alignment, metadata=metadata, n=200, lineage='L2', seed=456) 
subsample(fasta=full_alignment, metadata=metadata, n=200, lineage='L4', seed=456) 
subsample(fasta=full_alignment, metadata=metadata, n=200, lineage='L2', seed=789) 
subsample(fasta=full_alignment, metadata=metadata, n=200, lineage='L4', seed=789) 