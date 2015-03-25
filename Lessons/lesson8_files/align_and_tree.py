import sys
import os
import shutil

##### Input argument parsing and checking #####
usage = "\n\nUsage: python <directory> <alphabet> <seed>.\n\n<directory> refers to where files are stored \n<alphabet> should be either dna or protein (case insensitive) \n<seed> is a random integer to make RAXML results reproducible."

assert( len(sys.argv) == 4 ), usage

direc    = sys.argv[1]
alphabet = sys.argv[2].lower()
seed = sys.argv[3] # keep as string!

# Make sure direc ends with a slash, because we'll need to use it as part of a path!
if not direc.endswith("/"):
    direc += "/"
# Make sure alphabet is dna or protein
assert( alphabet == "dna" or alphabet == "protein" ), "Must specify either dna or protein as the alphabet."
#################################################



# Grab all the file names we need to process
files = os.listdir(direc)

# Perform alignment, phylogeny on each file
for file in files:
    # Make sure it's a file we want!
    if file.endswith(".fasta"):
        
        print "Processing", file
        
        # Determine a "job name" to name the output alignment and tree
        name = file.split(".fasta")[0]
        alignfile = direc + name + "_aligned.fasta"
        treefile  = direc + name + ".tree"
        
        # Build alignment
        command1 = "mafft --quiet " + direc+file + " > " + alignfile
        align = os.system(command1)
        assert( align == 0 ), "Alignment didn't run!"
        
        # Build phylogeny, and set up the RAxML command based on alphabet
        # Note that we send to out.txt just because raxml is noisy.
        if alphabet == "dna":
            command2 = "raxmlHPC-SSE3 -m GTRGAMMA -s " + alignfile + " -n tree -p " + seed + " > out.txt"
        elif alphabet == "protein":
            command2 = "raxmlHPC-SSE3 -m PROTGAMMALG -s " + alignfile + " -n tree -p " + seed + " > out.txt" 
        tree = os.system(command2)
        assert(tree == 0), "RAxML didn't run!"
        
        # Rename and save final tree, and remove raxml vomit
        shutil.move("RAxML_result.tree", direc + name + ".tree")
        os.system("rm RAxML_*") # To remove with a regular expression, use os.system (os.remove can't handle regex)
        os.remove("out.txt")
        
        
    

