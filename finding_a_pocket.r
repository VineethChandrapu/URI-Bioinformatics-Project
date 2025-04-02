# Part 1: Accept a dataframe and construct a matrix where both axes are all atoms in an amino acid and values are distance between them
    # Diagonal values should be 0
# Take minimum of these values as "necessary closeness of atoms"

# Part 2: Compile distance between atoms of all amino acids in a protein in a matrix of all_amino_acids x all_amino_acids

# Part 3: Run upper triangle with statistical significance against "closeness" determined in Part 1
# IF multiple mutation sites are close, this is a pocket