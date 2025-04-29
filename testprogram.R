
# PDB server connection required - testing excluded
install.packages("bio3d", dependencies=TRUE)
library(bio3d)


##--- Distance Matrix Plot
pdb <- read.pdb( "6V9C" )
k <- dm(pdb,inds="calpha")
filled.contour(k, nlevels = 10)

## NOTE: FOLLOWING EXAMPLE NEEDS MUSCLE INSTALLED
if(check.utility("muscle")) {
  
  ##--- DDM: Difference Distance Matrix
  # Downlaod and align two PDB files
  pdbs <- pdbaln( get.pdb( c( "4q21", "521p"), path = tempdir() ), outfile = tempfile() )
  
  # Get distance matrix
  a <- dm.xyz(pdbs$xyz[1,])
  b <- dm.xyz(pdbs$xyz[2,])
  a
  b

  
  #imports #linter
  install.packages("rPDBapi")
  library(rPDBapi)
pdb_file <- get_pdb_file(pdb_id = "6V9C", filetype = "pdb", save = TRUE, path = "PDBFiles")
