#imports #linter

### MAIN COPY - DO NOT TOUCH
CES <- function(){
    #cancer effect size
    print("help")
}

getPDB <- function(PDBcode){
    filepath <- ""
    print(filepath)
    return (filepath)
}

dataParsing <- function(PDB,CES){
    Cordinates <- c(1:24) #This vector should be as big the amount of entries we have
return (Cordinates)
}

## uhhh cancereffecetsizer tutorial 
options(timeout = 600)
install.packages(c("remotes", "BiocManager"))
remotes::install_github("Townsend-Lab-Yale/cancereffectsizeR",
 dependencies = TRUE, repos = BiocManager::repositories())
dir.create("CES_tutorial")
setwd("CES_tutorial")
#loading CES and data.table package
library(cancereffectsizeR)
library(data.table)
remotes::install_github("nvelden/NGLVieweR")
library(NGLVieweR)

# downloading new package readr
install.packages("readr")
library(readr)
#downloading variant data from cBioPortal 
### UPDATE CHANGE FILE TO PULL A DOWNLOADED LOCAL FILE 
TCGA_maf_file <- "TCGA-LIHC.maf.gz"
if (!file.exists(TCGA_maf_file)) {
get_TCGA_project_MAF(project = "LIHC", file = TCGA_maf_file)
}
tcga_clinical <- fread(system.file())
## trying with CES functions
get_ces_signature_set(refset = )

### pulling MAF file fro cBIOportal; NOTE: after path to, use your actual path to downloaded MAF file
### complete later, to makke path generic
maf_data <-load_maf("path/to/..")

### pulling MAF file from local directory
maf_data <- load_maf("TCGA-HCC.maf.gz")