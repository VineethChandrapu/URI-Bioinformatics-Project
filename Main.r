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
remotes::install_github("Townsend-Lab-Yale/cancereffectsizeR", dependencies = TRUE, repos = BiocManager::repositories())


### pulling MAF file fro cBIOportal; NOTE: after path to, use your actual path to downloaded MAF file
maf_data <-load_maf("path/to/..")