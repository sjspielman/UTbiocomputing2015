# Homework 4, question 2.
# Clean up this script such that it contains *3* functions (remember that you can use both positional and keyword arguments!)
# Your final script should produce the same output as this script.

import random

# Amino acid definitions. These lists should be kept at the top of the python script and used as global variables.
amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
polar       = ["Q", "N", "H", "S", "T", "Y", "C", "M", "W"]
nonpolar    = ["A", "I", "L", "F", "V", "P", "G"]
charged     = ["R", "K", "D", "E"]


######################### BELOW THIS LINE IS WHAT YOU SHOULD BE CHANGING ############################

# Generate a random protein sequence of a user-specified length, and print the sequence to screen.
size = raw_input("Enter the desired protein sequence length:\n")
protein = random.sample(amino_acids, int(size))
protein_string = "".join(protein)
print "The randomly generated protein sequence of length", size, "is", protein_string

# Determine the percentage of polar, nonpolar, and charged amino acids in the protein sequence, and print the result to screen.
total_polar = 0
for p in polar:
    total_polar += protein_string.count(p)
percent_polar = float(total_polar) / float(len(protein_string)) * 100
print "The percentage of polar residues is", percent_polar

total_nonpolar = 0
for np in nonpolar:
    total_nonpolar += protein_string.count(np)
percent_nonpolar = float(total_nonpolar) / float(len(protein_string)) * 100
print "The percentage of nonpolar residues is", percent_nonpolar


total_charged = 0
for c in charged:
    total_charged += protein_string.count(c)
percent_charged = float(total_charged) / float(len(protein_string)) * 100
print "The percentage of charged residues is", percent_charged


# From the original protein sequence, create a "subset" protein for each of the three groups (polar, nonpolar, charged) in which all non-polar/nonpolar/charged are replaced by "-"
# Print the resulting sequence to screen.

# Polar sequence. Turn nonpolar and charged residues into gaps
polar_only = ""
for residue in protein_string:
    if residue not in polar:
        polar_only += "-"
    else:
        polar_only += residue
print "The polar-only sequence is", polar_only
        
        
# Nonpolar sequence. Turn polar and charged residues into gaps
nonpolar_only = ""
for residue in protein_string:
    if residue not in nonpolar:
        nonpolar_only += "-"
    else:
        nonpolar_only += residue
print "The nonpolar-only sequence is", nonpolar_only


# Charged sequence. Turn polar and nonpolar residues into gaps
charged_only = ""
for residue in protein_string:
    if residue not in charged:
        charged_only += "-"
    else:
        charged_only += residue
print "The charged-only sequence is", charged_only














