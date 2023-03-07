##-------------------------------------------------------------------------
## Logistic regression analyses with correction for phylogenetic dependence
## to test which variables are associated to transmitter status 
## as inferred by TransPhylo
## 2021-11-26 Etthel Windels
##-------------------------------------------------------------------------



# Load libraries ----------------------------------------------------------

library(ape)
library(dplyr)
library(phylolm)
library(gdata)


# Read metadata -----------------------------------------------------------

metadata <- read.table("path_to_data", sep=',', dec='.', header=T)


# Read maximum likelihood tree --------------------------------------------

tree <- ape::read.tree('path_to_MLtree')
metadata <- metadata[metadata$g_number %in% tree$tip.label,]
metadata <- metadata %>% slice(match(tree$tip.label, g_number))
row.names(metadata) <- metadata$g_number


# Transform variables -----------------------------------------------------

for (i in 1:length(metadata$transmitter)){
  if (!is.na(metadata$transmitter[i])){
    if (metadata$transmitter[i]=='yes'){
        metadata$transmitter[i] <- 1
      }
    else if (metadata$transmitter[i]=='no'){
       metadata$transmitter[i] <- 0
    }
  }
}

metadata$transmitter <- as.numeric(metadata$transmitter)
metadata$mtb_lineage <- as.factor(metadata$mtb_lineage)
metadata$MDR <- ifelse(metadata$dr_category=='pan-susceptible',0,1)


# Select L2 and L4 isolates with inferred transmitter status --------------

metadata_transm <- metadata[!is.na(metadata$transmitter),]
metadata_transm <- metadata_transm[metadata_transm$mtb_lineage %in% c("L2","L4"),]
metadata_transm$mtb_lineage <- drop.levels(metadata_transm$mtb_lineage)


# GLM analyses ------------------------------------------------------------

# reference = L2
metadata_transm$mtb_lineage <- relevel(metadata_transm$mtb_lineage,"L2","L4")
levels(metadata_transm$mtb_lineage)
glm.fit <- glm(transmitter ~ mtb_lineage*MDR  + age + sex + hiv_status, data=metadata_transm, family = "binomial"(link=logit))
summary(phyloglm(transmitter ~  mtb_lineage*MDR  + age + sex + hiv_status, data=metadata_transm, tree, method='logistic_IG10', start.beta=coefficients(glm.fit), btol=10))

# reference = L4
metadata_transm$mtb_lineage <- relevel(metadata_transm$mtb_lineage,"L4","L2")
levels(metadata_transm$mtb_lineage)
glm.fit <- glm(transmitter ~ mtb_lineage*MDR  + age + sex + hiv_status, data=metadata_transm, family = "binomial"(link=logit))
summary(phyloglm(transmitter ~  mtb_lineage*MDR  + age + sex + hiv_status, data=metadata_transm, tree, method='logistic_IG10', start.beta=coefficients(glm.fit), btol=10))


# Get sample sizes per group ----------------------------------------------

metadata_transm <- metadata_transm[!is.na(metadata_transm$transmitter) & 
                                   !is.na(metadata_transm$mtb_lineage) & 
                                   !is.na(metadata_transm$MDR) & 
                                   !is.na(metadata_transm$age) & 
                                   !is.na(metadata_transm$sex) & 
                                   !is.na(metadata_transm$hiv_status),]
table(metadata_transm[metadata_transm$mtb_lineage=='L2' & metadata_transm$MDR==0,]$transmitter)
table(metadata_transm[metadata_transm$mtb_lineage=='L2' & metadata_transm$MDR==1,]$transmitter)
table(metadata_transm[metadata_transm$mtb_lineage=='L4' & metadata_transm$MDR==0,]$transmitter)
table(metadata_transm[metadata_transm$mtb_lineage=='L4' & metadata_transm$MDR==1,]$transmitter)
table(metadata_transm[metadata_transm$sex=='Female',]$transmitter)
table(metadata_transm[metadata_transm$sex=='Male',]$transmitter)
table(metadata_transm[metadata_transm$hiv_status=='HIV-',]$transmitter)
table(metadata_transm[metadata_transm$hiv_status=='HIV+',]$transmitter)
