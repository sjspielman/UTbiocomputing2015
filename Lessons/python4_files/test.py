import csv
import random

def compute_at_content(seq, digits = 5, percent = False):
    """ 
        Compute AT content from a provided sequence. By default, result is rounded to 5 significant digits.
        Positional arguments: 
            1. seq:  DNA sequence from which GC-content should be calculated
        Keyword arguments:
            1. digits: number of significant digits for results (default: 5)
            2. percent: return value as a percent rather than decimal (default: False)
    """
    a = seq.lower().count("a")
    t = seq.lower().count("t")
    at = round( float(a + t) / len(seq), digits)
    if percent:
        return at*100
    else:
        return at
        
        
        
def compute_gc_content(seq, digits = 5, percent = False):
    """ 
        Compute GC content from a provided sequence. By default, result is rounded to 5 significant digits.
        Positional arguments: 
            1. seq:  DNA sequence from which GC-content should be calculated
        Keyword arguments:
            1. digits: number of significant digits for results (default: 5)
            2. percent: return value as a percent rather than decimal (default: False)
    """
    g = seq.upper().count("G")
    c = seq.upper().count("C")
    gc = float(g + c) / len(seq)
    if percent:
        return round(gc * 100, digits)
    else:
        return round(gc, digits)
        
        
        
'''
infile = "crabs.csv"
with open(infile, 'rU') as inf:
    reader = csv.reader(inf)
    for row in reader:
        print row
    inf.seek(0)
    for row in reader:
        print "NEW", row

infile2 = "crabs_tdf.txt"
inf2 = open(infile2, "rU")
reader2 = csv.reader(inf2, delimiter = '\t')
for row in reader2:
    print row
inf2.seek(0)
reader_dict = csv.DictReader(inf2, delimiter = '\t')
for row in reader_dict:
    print row
inf2.close()
'''

# Create and write a dataset
dna = ["A", "C", "G", "T"]*1000
random.shuffle(dna)
with open("out.txt", "w") as outf:
    my_writer = csv.writer(outf, delimiter="\t")
    my_writer.writerow(["id", "at_content", "gc_content"])
    for i in range(20):
        seq = "".join( random.sample(dna, 10) )
        at = compute_at_content(seq)
        gc = compute_gc_content(seq)
        my_writer.writerow([seq, at, gc])

    
# Another way to write the dataset
with open("out2.txt", "w") as outf:
    outf.write("id,at_content,gc_content\n")
    for i in range(20):
        sequence = "".join( random.sample(dna, 10) )
        at = compute_at_content(sequence)
        gc = compute_gc_content(sequence)
        outf.write( str(i) + "," + str(at) + "," + str(gc) + '\n' )
    
    
    
    









