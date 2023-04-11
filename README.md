# Phylodynamic estimation of the relative transmission fitness of MDR-TB
This repository contains the code and metadata associated with the phylodynamic analyses performed in
Loiseau, Windels et al. The relative transmission fitness of multidrug-resistant *Mycobacterium tuberculosis* in a drug resistance hotspot. *Nature Communications* 14, 1988 (2023). https://doi.org/10.1038/s41467-023-37719-y

These analyses aimed at estimating the transmission fitness (effective reproductive number and transmission rate) of multidrug-resistant *M. tuberculosis* strains, relative to their drug-susceptible counterparts. Using genome sequences sampled from TB patients in the country of Georgia (2014-2016), we fitted different multitype birth-death models in BEAST2 and inferred transmission trees using TransPhylo in R.

The raw sequencing data used in this study are available under project accession numbers [PRJEB39561](https://www.ebi.ac.uk/ena/browser/view/PRJEB39561) and [PRJEB50582](https://www.ebi.ac.uk/ena/browser/view/PRJEB50582).

- The folder `analyses/` contains the XML files used for the structured birth-death analyses in BEAST2.
- The folder `data/` contains the metadata file, including sample accession numbers.
- The folder `scripts/` contains R scripts used for pre-processing of the alignments, post-processing of the BEAST2 output, transmission tree inference (TransPhylo), and post-processing of the TransPhylo output.

