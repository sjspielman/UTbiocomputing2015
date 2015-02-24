import re


file = "hyphy_output_long.txt"

f = open(file, "r")
lines = f.readlines()
f.close()

# Collect values of interest
values = []
for line in lines:
    # Does this line contain site info?
    is_it_site = re.search("^Site", line)
    if is_it_site:
        newline = re.split("\s+", line)
        values.append( newline[4] ) # note that this is a STRING!!!
        
    # Does this line contain the scaling factor I want?
    find_raw_scaling = re.search("^Raw scaling factor:(.+)$", line)
    #find_raw_scaling = re.search("^Raw scaling factor:\d\.\d+$", line) # problem line since not general enough!!
    if find_raw_scaling:
        scale_factor = float( find_raw_scaling.group(1) )
        
#print scale_factor
#print values

# Two examples of saving to file
with open("dnds_outfile.txt", "w") as outf:
    
    for dnds in values:
        outf.write(str(dnds) + "\n") # why the \n?

with open("dnds_outfile2.txt", "w") as outf:
    outf.write("dnds\n") # header!
    outf.write("\n".join(values))

