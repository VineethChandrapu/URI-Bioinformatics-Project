#imports #linter
library(rPDBapi)
library(cancereffectsizeR)
library(data.table)
library(ggplot2)


CES <- function(pdbcode){
    #cancer effect size
    #geting maf files based of TCGA database here we are selecting Liver Hepatocellular Carcinoma
    tcga_maf_file <- "TCGA-lihc.maf.gz"
    if (!file.exists(tcga_maf_file)) {
        get_TCGA_project_MAF(project = "TCGA-lihc", filename = "TCGA-lihc.maf.gz")
    }
    maf <- preload_maf(maf = tcga_maf_file, refset = "ces.refset.hg38")
    cesa <- CESAnalysis(refset = "ces.refset.hg38")
    cesa <- load_maf(cesa = cesa, maf = maf)
    signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "LIHC", treatment_naive = TRUE)
    cesa <- trinuc_mutation_rates(
    cesa = cesa,
    signature_exclusions = signature_exclusions
    )
    cesa <- gene_mutation_rates(cesa, covariates = ces.refset.hg38$covariates$liver)
    cesa <- ces_variant(cesa, run_name = "example")
    plot_effects(effects = cesa$selection$example, color_by = "#DB382D", topn = 20)
    mut_effects <- mutational_signature_effects(cesa, cesa$selection$example)
    plot_signature_effects(mut_effects, viridis_option = "F", num_sig_groups = 5)



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


