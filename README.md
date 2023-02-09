# Habitat and tree species identity shape aboveground and belowground fungal communities in Central European forests 
 
This repository contains the data and pre-processing pipeline, as well as the R script for the analysis accompanying our paper: 

> Hofmann B, Dreyling L, Dal Grande F, Otte J & Schmitt I 2023. Habitat and tree species identity shape aboveground and belowground fungal communities in Central European forests. Frontiers in Microbiology doi: [10.3389/fmicb.2023.1067906](https://www.frontiersin.org/articles/10.3389/fmicb.2023.1067906/abstract).

## Contacts

**Lukas Dreyling**  
Doctoral Candidate  
[E-Mail](mailto:lukas.dreyling@senckenberg.de)  

**Imke Schmitt**  
Principal Investigator  
[E-Mail](mailto:imke.schmitt@senckenberg.de)  

**Francesco Dal Grande**  
Principal Investigator  
[E-Mail](mailto:francesco.dalgrande@unipd.it)  

## Contents

1. [Pre-Processing Pipeline](01_processing_pipeline.txt)
2. [Data](02_Data.zip)
3. [Analysis Script](03_data_analysis.R)

Additionally you need to download the raw reads [here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA819266).  

If you lack the computing power to process the raw reads, the resulting ASV tables, FASTA files, taxonomy table and metadata are located [here](02_Data.zip).  

## Before starting

### You will need to have the following software installed.

#### Pre-Processing 
* fastQC http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
* Cutadapt https://cutadapt.readthedocs.io/
* R https://www.r-project.org/
    - dada2 https://benjjneb.github.io/dada2/index.html
    - ShortRead https://kasperdanielhansen.github.io/genbioconductor/html/ShortRead.html
    - Biostrings https://bioconductor.org/packages/release/bioc/html/Biostrings.html
* BLASTn https://www.ncbi.nlm.nih.gov/books/NBK279690/

#### Analysis
* R https://www.r-project.org/
* Rstudio https://www.rstudio.com/
  - here https://here.r-lib.org/
  - tidyverse https://www.tidyverse.org/
  - decontam https://benjjneb.github.io/decontam/
  - phyloseq https://joey711.github.io/phyloseq/
  - LULU https://github.com/tobiasgf/lulu
  - Biostrings https://bioconductor.org/packages/release/bioc/html/Biostrings.html
  - microbiome https://github.com/microbiome/microbiome
  - fantaxtic https://github.com/gmteunisse/Fantaxtic
  - igraph https://igraph.org/
  - SpiecEasi https://github.com/zdk123/SpiecEasi
  - ranacapa https://github.com/gauravsk/ranacapa
  - rgexf https://github.com/gvegayon/rgexf
  - vegan https://rdrr.io/cran/vegan/man/vegan-package.html
  - MicEco https://github.com/Russel88/MicEco
