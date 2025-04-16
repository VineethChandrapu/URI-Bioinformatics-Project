#imports #linter
library(rPDBapi)
library(cancereffectsizeR)
library(data.table)
library(ggplot2)


CES <- function(pdbcode){
####THIS CHUNK JUST RUNS CES TUTORIAL ON LOCAL DATA####
my_maf<- preload_maf("mc3.v0.2.8.PUBLIC.maf", refset = "ces.refset.hg19")

cesa<- CESAnalysis(refset = "ces.refset.hg19")

cesa <- load_maf(cesa, maf = my_maf)



}

getPDB <- function(){
    #downloads pdb file to local folder
    code <- readline(prompt = "enter pdb code: ")
    pdb_file <- get_pdb_file(pdb_id = code, filetype = "pdb", save = TRUE, path = "PDBfiles")
    return (code)
}

dataParsing <- function(PDB,CES){
    Cordinates <- c(1:24) #This vector should be as big the amount of entries we have
    
    return (Cordinates)
}
#"6V9C"
#pdb_code <- getPDB()
CES(pdb_code)


