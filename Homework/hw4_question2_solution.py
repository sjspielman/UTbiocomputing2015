# SOLUTION. Contact stephanie.spielman@gmail.com with any questions!
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

def generate_seq():
    ''' Prompt users for an integer, which is used to make a random protein sequence of that size.
        Return the randomly-generated protein sequence as a string.
    '''
    size = raw_input("Enter the desired protein sequence length:\n")
    protein = random.sample(amino_acids*500, int(size)) # hack! enjoy :)
    protein_string = "".join(protein)
    print "\n\nThe randomly generated protein sequence of length", size, "is", protein_string
    return protein_string
    
def compute_biochem(type, protein):
    ''' Compute the percentage of a certain type of amino acid (polar, nonpolar, or charged).
        Argument "type" is a string which must be either polar, nonpolar, or charged.
        Argument "protein" is a protein string on which we make calculations.
        Does not return any values.
    '''
    total = 0.
    if type == "polar":
        a = polar
    elif type == "nonpolar":
        a = nonpolar
    elif type == "charged":
        a = charged
    for entry in a:
        total += protein.count(entry)
    percent = (total / len(protein)) * 100
    print "The percentage of", type, "residues is", percent

def extract_biochem(type, protein):
    ''' 
        Extract from the provided "protein" string the sub-string of amino acids which belong to the specified "type" (polar, nonpolar, charged).
        All non-polar/nonpolar/charged are replaced by "-"
        Does not return any values.
    '''
    sub_protein = ""
    if type == "polar":
        a = polar
    elif type == "nonpolar":
        a = nonpolar
    elif type == "charged":
        a = charged
    for residue in protein:
        if residue not in a:
            sub_protein += "-"
        else:
            sub_protein +=  residue
    print "The", type, "-only sequence is", sub_protein
            


protein = generate_seq()
print
compute_biochem("polar", protein)
print
compute_biochem("nonpolar", protein)
print
compute_biochem("charged", protein)
print
extract_biochem("polar", protein)
print
extract_biochem("nonpolar", protein)
print
extract_biochem("charged", protein)
print






