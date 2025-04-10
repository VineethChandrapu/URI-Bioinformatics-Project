# ---------- Read In and Organize Data to Be Used ----------
# Make a Dataframe:
# REPLACE WITH ACTUAL DATA (needs to be read into data frame)
df <- data.frame(
  AA <- c("AA1", "AA1", "AA3"),
  Atom <- c("a1", "a2", "a1"),
  Mut <- c(TRUE, FALSE, TRUE), # I don't know how mutations are indicated right now
  x <- c(3, 4, 5),
  y <- c(7, 9, 3),
  z <- c(2, 8, 1)
)
colnames(df) <- c("AminoAcid", "Atom", "Mutation?", "X", "Y", "Z")
df

# Dataframe of mutation points
mut_df <- df[df$Mut == TRUE,]
mut_df














# --------- GROSS VVV ----------


# ---------- Part 1: Construct a Distance Matrix of ALL Amino Acids ----------
# Input: dataframe of x,y,z values, labels of molecules
# Output: matrix with distance between each molecule
    # Diagonal values should be 0

# Make a Matrix:
coords_mat <- as.matrix(coords_df)
rownames(coords_mat) <- c("Atom1", "Atom2", "Atom3") #REPLACE WITH ROW NAMES VECTOR
coords_mat

dist_mat <- dist(coords_mat, method = "euclidean", diag = TRUE, upper = FALSE)
dist_mat

# Greatest distance between ATOMS (do we want this to be a metric for "closeness"?)
max(dist_mat)


# ---------- Part 2:  ----------


# ---------- Part 3: Significance Matrix -----------
# Run upper triangle with statistical significance against "closeness" 


# ---------- Part 4: Reconstructing a Dataframe ----------
sig_df <- data.frame()