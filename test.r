

# Main function
PullProteinStructure <- function(gene_name) {
  cat("Looking up human UniProt ID for gene:", gene_name, "\n")
 
  tryCatch({
    uniprot_id <- get_human_uniprot_id(gene_name)
    cat("Found human UniProt ID:", uniprot_id, "\n")
   
    pdb_file <- download_alphafold_pdb(uniprot_id, gene_name)
    cat("Process completed. File saved as:", pdb_file, "\n")
    return(pdb_file)
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    return(NULL)
  })
}


# Get gene name from command line argument or use default
gene_name <- "6V9C" 


# Execute main function
PullProteinStructure(gene_name)
