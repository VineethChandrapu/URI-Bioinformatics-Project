#installing r packages this will take a while
install.packages(c("remotes", "BiocManager"))
remotes::install_github("Townsend-Lab-Yale/cancereffectsizeR", dependencies = TRUE, repos = BiocManager::repositories())
install.packages("rPDBapi", repos = "http://cran.us.r-project.org")
BiocManager::install(version = "3.20") # or newer, when available
remotes::install_github("Townsend-Lab-Yale/ces.refset.hg38@*release")
