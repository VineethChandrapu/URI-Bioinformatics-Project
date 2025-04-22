# ---------- Libraries ----------
if(!require(dplyr)){install.packages("dplyr")}
if(!require(dplyr)){install.packages("tidyverse")}
if(!require(dplyr)){install.packages("purrr")}
library(tidyverse)
library(dplyr)
library(purrr)
library(stringr)

# ---------- Read In and Organize Data to Be Used ----------
# Create a fake data frame for testing
df <- data.frame(
  AA <- c("AA1", "AA1", "AA1", "AA2", "AA2", "AA3", "AA3", "AA3"),
  Atom <- c("a1", "a2", "a3", "a1", "a2", "a1", "a2", "a3"),
  Mut <- c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE), # I don't know how mutations are indicated right now
  x <- c(3, 4, 5, 4, 7, 2, 8, 5),
  y <- c(7, 9, 3, 7, 5, 7, 2, 4),
  z <- c(2, 8, 1, 6, 9, 3, 1, 7)
)
colnames(df) <- c("AminoAcid", "Atom", "Mutation?", "X", "Y", "Z")
df

# Dataframe of mutated points only:
mut_df <- df[df$Mut == TRUE,]


# ---------- Part 1: Construct a Distance Matrix of ALL Amino Acids ----------
# Input: dataframe of all amino acids, atoms, and coordinates
# Output: matrix with minimum pairwise distance between each amino acid 

# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Make Matrices for each Amino Acid: (DOESN'T WORK RN)
#for (AA in df$AminoAcid)
#  {
#    coords_df <- subset(df, AminoAcid == AA, select=c('X', 'Y', 'Z'))
#
#    coords_mat <- as.matrix(coords_df)
#    rownames(coords_mat) <- c(df$Atom) #REPLACE WITH ROW NAMES VECTOR
#    coords_mat
#    
#    dist_mat <- dist(coords_mat, method = "euclidean", diag = TRUE, upper = FALSE)
#    dist_mat
#}
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# Find distance between every 2 amino acids:

# Group the data {use: grouped_by_AA}
  location_df <- df[,-3]
  grouped_by_AA <- df %>%
    group_by(AminoAcid) %>%
    nest()
  grouped_by_AA

# Euclidean function {use: euclidean_distance[i, j] }
    euclidean_distance <- function(point1, point2) {
      sqrt(sum((point1 - point2)^2))
    }

# How to use just one dataframe:
nrow(grouped_by_AA$data[[1]])

# Outputting wrong
df_lengths <- grouped_by_AA %>%
  for(i in .)
  {
    nrow(grouped_by_AA[[i]]) # Get the number of rows (length of col1) for each dataframe
  }
max(df_lengths)


# Looping through each dataframe to get distance matrices
for(i in 2:nrow(grouped_by_AA$data[[i]]))
  {
  for (j in 2:nrow(grouped_by_AA$data[[i+1]]))
    {
    distances <- matrix(NA, nrow = nrow(grouped_by_AA$data[[i-1]]), ncol = ncol(grouped_by_AA$data[[i]]))
    distances[i, j] <- euclidean_distance(
      c(grouped_by_AA$data[[i-1]]$X[i], grouped_by_AA$data[[i-1]]$Y[i], grouped_by_AA$data[[i-1]]$Z[i]),
      c(grouped_by_AA$data[[i]]$X[j], grouped_by_AA$data[[i]]$Y[j], grouped_by_AA$data[[i]]$Z[j])
    )
  #assign(i, data.frame(split(dat[[i]], rep(letters[1:5], each = 3))))
  print(distances)
  }
}



