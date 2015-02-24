import csv

infile_csv = "AbilAhuG_uniprot_blastx.csv"
infile_tdf = "AbilAhuG_uniprot_blastx.txt"


with open(infile_csv, 'r') as inf:
    reader = csv.reader(inf)
    for row in reader:
        print row
    #inf.seek(0)
    for row in reader:
        print "NEW", row

infile2 = "crabs_tdf.txt"
inf2 = open(infile_tdf, 'r')
reader2 = csv.reader(inf2, delimiter = '\t')
for row in reader2:
    print row
inf2.seek(0)
reader_dict = csv.DictReader(inf2, delimiter = '\t')
for row in reader_dict:
    print row
inf2.close()