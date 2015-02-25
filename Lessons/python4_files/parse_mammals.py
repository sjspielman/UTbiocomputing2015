import re

file = open("mammal_dat.nex", "r")
lines = file.readlines()
file.close()

# re.search will return None if the pattern is not matched
for line in lines:
    find_dim = re.search("^\s+Dimensions ntax=(\d+) nchar=(\d+)", line)
    if find_dim:
        print "found!!", line
        # Use .group(n) to *capture* what you placed in parentheses.
        ntaxa = find_dim.group(1) # these are strings, but can recast
        nchar = int(find_dim.group(2)) #recasted as integer
        break # stop looping over the file if we've found what we're looking for



# re.findall will return an empty list if nothing is found
for line in lines:
    # Capture the opossum sequence
    find_opossum = re.search("^Opossum\s+(.+)$", line)
    if find_opossum:
        sequence = find_opossum.group(1)
        
        # Find all ANNC motifs in the opossum sequence
        AnnC_motif = re.findall("A\w\wC", sequence)
        #print AnnC_motif
        
        # For fun put them in a dictionary!
        AnnC_motif_dict = {}
        for entry in AnnC_motif:
            if entry in AnnC_motif_dict:
                AnnC_motif_dict[entry] += 1
            else:
                AnnC_motif_dict[entry] = 1
        print AnnC_motif_dict
        break # no need to keep looping!
        
        
# re.split will return a list split on the specified regex. If no match, returns the whole string as a list (see why???)
species_dict = {}
for line in lines:
    is_species = re.search("^[A-Z].+\s+.+$", line)
    if is_species:
        l = re.split("\s+", line) # each line has a different number of spaces between the taxon name and sequence, so we can just split on a space!!
        species_dict[ l[0] ] = l[1]
print species_dict




        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
    
    
    
    
    
    
    
    
    