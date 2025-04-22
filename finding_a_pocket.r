install.packages("tidyverse")
library(tidyverse)
library(dplyr)
library(purrr)

# ---------- Read In and Organize Data to Be Used ----------
# Make a Dataframe:
# REPLACE WITH ACTUAL DATA (needs to be read into data frame)
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

# Dataframe of mutation points
# mut_df <- df[df$Mut == TRUE,]


# ---------- Part 1: Construct a Distance Matrix of ALL Amino Acids ----------
# Input: dataframe of all amino acids, atoms, and coordinates
# Output: matrix with minimum pairwise distance between each amino acid 
# Diagonal values should be 0


# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Make Matrices for each Amino Acid: (DOESN'T WORK RN)
for (AA in df$AminoAcid)
  {
    coords_df <- subset(df, AminoAcid == AA, select=c('X', 'Y', 'Z'))

    coords_mat <- as.matrix(coords_df)
    rownames(coords_mat) <- c(df$Atom) #REPLACE WITH ROW NAMES VECTOR
    coords_mat
    
    dist_mat <- dist(coords_mat, method = "euclidean", diag = TRUE, upper = FALSE)
    dist_mat
}
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


# Find distance between 2 amino acids:
  # Group the data:
location_df <- df[,-3]
grouped_by_AA <- df %>%
  group_by(AminoAcid) %>%
  nest()
grouped_by_AA


  # Construct pairwise matrices for all amino acids:
euclidean_distance <- function(point1, point2) {
  sqrt(sum((point1 - point2)^2))
}

just_AA2 <- grouped_by_AA$data[[2]]
print(just_AA2)

for(i in 2:grouped_by_AA$data) {
  distances <- matrix(NA, nrow = nrow(grouped_by_AA$data[[i-1]]), ncol = nrow(grouped_by_AA$data[[i]]))
  distances[i, j] <- euclidean_distance(
    c(grouped_by_AA$data[[i-1]]$X[i], grouped_by_AA$data[[i-1]]$Y[i], grouped_by_AA$data[[i-1]]$Z[i]),
    c(grouped_by_AA$data[[i]]$X[j], grouped_by_AA$data[[i]]$Y[j], grouped_by_AA$data[[i]]$Z[j])
  )
  assign(i, data.frame(split(dat[[i]], rep(letters[1:5], each = 3))))
  print(distances)
}


# Construct pairwise matrices for all amino acids:
just_AA2 <- grouped_by_AA$data[[2]]
print(just_AA2)

nrow(grouped_by_AA$data[[1]])

for(i in grouped_by_AA$data) {
  matrix(sapply(1:nrow(grouped_by_AA$data[[1]]), function(x) euclidean_distance(mat1[x,], mat2[x,])), ncol = 5, byrow = FALSE)
}


# --------- GROSS VVV ----------





# ---------- Part 2:  ----------


# ---------- Part 3: Significance Matrix -----------
# Run upper triangle with statistical significance against "closeness" 


# ---------- Part 4: Reconstructing a Dataframe ----------
sig_df <- data.frame()