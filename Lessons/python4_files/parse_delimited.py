import csv

infile_csv = "AbilAhuG_uniprot_blastx.csv"
infile_tdf = "AbilAhuG_uniprot_blastx.txt"


with open(infile_csv, 'r') as inf:
    reader = csv.reader(inf)
    for row in reader:
        print row
    inf.seek(0) # required to loop another time. 
    for row in reader:
        print "NEW", row

# tab-delimited
inf2 = open(infile_tdf, 'r')
reader2 = csv.reader(inf2, delimiter = '\t')
for row in reader2:
    print row

# dictionaries instead of lists
inf2.seek(0)
reader_dict = csv.DictReader(inf2, delimiter = '\t')
for row in reader_dict:
    print row # each row is a dictionary!
inf2.close()