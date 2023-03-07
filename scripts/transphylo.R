##-------------------------------------------------------------------------
## Run TransPhylo on predefined clades, dated with BEAST2
## 2021-11-25 Etthel Windels
##-------------------------------------------------------------------------



# Load libraries ----------------------------------------------------------

library(ape)
library(stringr)
library(TransPhylo)
library(dplyr)
library(tibble)


# Set working directory ---------------------------------------------------

setwd('path_to_files')


# Read metadata -----------------------------------------------------------

metadata <- read.csv('path_to_data', sep=',', dec='.', header=T)


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

metadata$transmitter <- rep(NA, dim(metadata)[1]) 
metadata$MDR <- ifelse(metadata$dr_category=='pan-susceptible', 0, 1)
metadata <- 
  metadata %>%
  add_column(new_id = paste(metadata$g_number, metadata$MDR, metadata$new_date, sep = "/"))


# Read dated trees --------------------------------------------------------

L2MCCtree1_nx <- ape::read.nexus("L2clade1_MCCtree")  
L2MCCtree2_nx <- ape::read.nexus("L2clade2_MCCtree")  
L2MCCtree3_nx <- ape::read.nexus("L2clade3_MCCtree")  
L2MCCtree4_nx <- ape::read.nexus("L2clade4_MCCtree")  
L2MCCtree5_nx <- ape::read.nexus("L2clade5_MCCtree")  
L4MCCtree1_nx <- ape::read.nexus("L4clade1_MCCtree")  
L4MCCtree2_nx <- ape::read.nexus("L4clade2_MCCtree")  
L4MCCtree3_nx <- ape::read.nexus("L4clade3_MCCtree")  
L4MCCtree4_nx <- ape::read.nexus("L4clade4_MCCtree")  


# Remove multifurcations - use 1 day as minimum branch length -------------

rm_multif <- function(tree){
  tree$tip.label <- str_replace_all(tree$tip.label,'_','/')
  tree.rm.multif <- multi2di(tree)
  tree.rm.multif$edge.length <- pmax(tree.rm.multif$edge.length,1/365)
  return(tree.rm.multif)
}

L2MCCtree1_nx <- rm_multif(L2MCCtree1_nx)
L2MCCtree2_nx <- rm_multif(L2MCCtree2_nx)
L2MCCtree3_nx <- rm_multif(L2MCCtree3_nx)
L2MCCtree4_nx <- rm_multif(L2MCCtree4_nx)
L2MCCtree5_nx <- rm_multif(L2MCCtree5_nx)
L4MCCtree1_nx <- rm_multif(L4MCCtree1_nx)
L4MCCtree2_nx <- rm_multif(L4MCCtree2_nx)
L4MCCtree3_nx <- rm_multif(L4MCCtree3_nx)
L4MCCtree4_nx <- rm_multif(L4MCCtree4_nx)


# Get last dates ----------------------------------------------------------------

get_last_date <- function(lineage, tree){
  if (lineage=='L2'){
    datesL2=read.table('dates_L2_iqtree.txt')
    dates <- datesL2[datesL2[,1] %in% tree$tip.label,]
    date.last <- max(dates[,2])
  }
  else if (lineage=='L4'){
    datesL4=read.table('dates_L4_iqtree.txt')
    dates <- datesL4[datesL4[,1] %in% tree$tip.label,]
    date.last <- max(dates[,2])
  }
  return(date.last)
}


# Run TransPhylo on all L2 trees ------------------------------------------

L2trees <- c(L2MCCtree1_nx,L2MCCtree2_nx,L2MCCtree3_nx,L2MCCtree4_nx,L2MCCtree5_nx)

for (j in L2trees){
  tree <- j
  date.last <- get_last_date("L2",tree)
  input_tree <- ptreeFromPhylo(tree,dateLastSample=date.last)
  
  # Starting values and settings
  
  neg <- 10/365       # Ne*g
  w.shape <- 1.3      # alpha parameter of generation time distribution (gamma distribution)
  w.scale <- 1        # beta parameter of generation time distribution (gamma distribution)
  ws.shape <- 1.3     # alpha parameter of sampling time distribution (gamma distribution)
  ws.scale <- 1       # beta parameter of sampling time distribution (gamma distribution)
  pi <- 0.5           # sampling density
  off.r <- 1
  off.p <- 0.5
  Re = off.r*off.p/(1-off.p); Re
  
  set.seed(1234)
  iters <- 1e5; thin <- 10
  
  # Infer transmission trees
  
  output_transphylo <-  inferTTree(input_tree, w.shape, w.scale, ws.shape, ws.scale, startPi=pi, startNeg = neg, startOff.r = off.r, startOff.p = off.p, 
                                   mcmcIterations = iters, thinning = thin, dateT = 2017, updatePi = TRUE, updateOff.r = TRUE, updateOff.p = FALSE, updateNeg = TRUE,
                                   verbose=T)

  # Infer transmitter status
  
  for (i in tree$tip.label){
    offspring_dist <- getOffspringDist(output_transphylo,  burnin = 0.2, show.plot = F, k=i)
    prob_transmitter <- sum(offspring_dist>0)/length(offspring_dist)
    metadata$transmitter[which(metadata$new_id==i)] <- ifelse(prob_transmitter>0.9, 'yes','no')
  }
}


# Run TransPhylo on all L4 trees ------------------------------------------

L4trees <- c(L4MCCtree1_nx,L4MCCtree2_nx,L4MCCtree3_nx,L4MCCtree4_nx)

for (j in L4trees){
  tree <- j
  date.last <- get_last_date("L4",tree)
  input_tree <- ptreeFromPhylo(tree,dateLastSample=date.last)
  
  # Starting values and settings
  
  neg <- 10/365       # Ne*g
  w.shape <- 1.3      # alpha parameter of generation time distribution (gamma distribution)
  w.scale <- 1        # beta parameter of generation time distribution (gamma distribution)
  ws.shape <- 1.3     # alpha parameter of sampling time distribution (gamma distribution)
  ws.scale <- 1       # beta parameter of sampling time distribution (gamma distribution)
  pi <- 0.5           # sampling density
  off.r <- 1
  off.p <- 0.5
  Re = off.r*off.p/(1-off.p); Re
  
  set.seed(1234)
  iters <- 1e5; thin <- 10
  
  # Infer transmission trees
  
  output_transphylo <-  inferTTree(input_tree, w.shape, w.scale, ws.shape, ws.scale, startPi=pi, startNeg = neg, startOff.r = off.r, startOff.p = off.p, 
                                   mcmcIterations = iters, thinning = thin, dateT = 2017, updatePi = TRUE, updateOff.r = TRUE, updateOff.p = FALSE, updateNeg = TRUE,
                                   verbose=T)
  
  # Infer transmitter status
  
  for (i in tree$tip.label){
    offspring_dist <- getOffspringDist(output_transphylo,  burnin = 0.2, show.plot = F, k=i)
    prob_transmitter <- sum(offspring_dist>0)/length(offspring_dist)
    metadata$transmitter[which(metadata$new_id==i)] <- ifelse(prob_transmitter>0.9, 'yes','no')
  }
}


# Save output file --------------------------------------------------------

write.csv(metadata, file="output_path")