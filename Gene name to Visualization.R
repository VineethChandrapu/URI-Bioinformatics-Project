# Install required packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("bio3d", quietly = TRUE))
  BiocManager::install("bio3d")
if (!requireNamespace("r3dmol", quietly = TRUE))
  install.packages("r3dmol")
if (!requireNamespace("rgl", quietly = TRUE))
  install.packages("rgl")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")

# Load libraries
library(bio3d)
library(r3dmol)
library(rgl)
library(biomaRt)

# Function to get UniProt ID for a human gene
get_uniprot_id <- function(gene_symbol) {
  # Connect to the Ensembl database
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  
  # Define the attributes you want to retrieve
  attributes <- c("hgnc_symbol", "uniprot_gn_id", "uniprot_gn_symbol")
  
  # Define the filter
  filters <- "hgnc_symbol"
  
  # Query the database
  results <- getBM(attributes = attributes, 
                   filters = filters, 
                   values = gene_symbol, 
                   mart = ensembl)  
  # Return the results
  return(results)
}

# Replace "TP53" with any human gene symbol you're interested in
uniprot_data <- get_uniprot_id("TP53")
print(uniprot_data)

download_alphafold_pdb <- function(uniprot_id, output_dir = ".") {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Define the AlphaFold Database URL
  alphafold_url <- paste0("https://alphafold.ebi.ac.uk/files/AF-", uniprot_id, "-F1-model_v4.pdb")
  
  # Define the output file path
  output_file <- file.path(output_dir, paste0("AF-", uniprot_id, ".pdb"))
  
  # Try to download the file
  tryCatch({
    message(paste0("Downloading AlphaFold structure for UniProt ID: ", uniprot_id))
    download.file(url = alphafold_url, destfile = output_file, mode = "wb")
    message(paste0("Downloaded to: ", output_file))
    return(output_file)
  }, error = function(e) {
    message(paste0("Error downloading structure for ", uniprot_id, ": ", e$message))
    return(NULL)
  })
}

# Function that combines UniProt lookup and AlphaFold download
get_and_download_structure <- function(gene_symbol, output_dir = ".") {
  # Get UniProt ID (using the previous function)
  uniprot_data <- get_uniprot_id(gene_symbol)
  
  # Check if any results were found
  if (nrow(uniprot_data) == 0) {
    message(paste0("No UniProt ID found for gene: ", gene_symbol))
    return(NULL)
  }
  
  # Get the first UniProt ID
  first_uniprot_id <- uniprot_data$uniprot_gn_id[1]
  
  # If the first entry is empty, try to find a non-empty one
  if (first_uniprot_id == "" || is.na(first_uniprot_id)) {
    non_empty_ids <- uniprot_data$uniprot_gn_id[uniprot_data$uniprot_gn_id != "" & !is.na(uniprot_data$uniprot_gn_id)]
    if (length(non_empty_ids) > 0) {
      first_uniprot_id <- non_empty_ids[1]
    } else {
      message(paste0("No valid UniProt ID found for gene: ", gene_symbol))
      return(NULL)
    }
  }
  
  message(paste0("Using UniProt ID: ", first_uniprot_id))
  
  # Download AlphaFold structure
  pdb_file <- download_alphafold_pdb(first_uniprot_id, output_dir)
  
  return(pdb_file)
}

# Replace "TP53" with any human gene symbol you're interested in
# Specify an output directory (default is current directory)
pdb_file <- get_and_download_structure("BRCA1", output_dir = "alphafold_structures")

# If you want to view the structure (requires the bio3d package)
if (!is.null(pdb_file) && requireNamespace("bio3d", quietly = TRUE)) {
  library(bio3d)
  pdb <- read.pdb(pdb_file)
  print(pdb)
}

# Using r3dmol (interactive 3D visualization in RStudio Viewer pane)
visualize_with_r3dmol <- function(pdb_file) {
  # Create an r3dmol viewer
  r3dmol() %>%
    m_add_model(data = pdb_file, format = "pdb") %>%
    m_set_style(style = list(cartoon = list(color = "spectrum"))) %>%
    m_zoom_to() %>%
    m_spin()
}
viz_r3dmol <- visualize_with_r3dmol(pdb_file)
viz_r3dmol 

