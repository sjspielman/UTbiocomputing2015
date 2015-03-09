from Bio import Entrez, SeqIO
Entrez.email = 'stephanie.spielman@gmail.com' # but please replace with your email!

def fetch_ncbi_record(id, db, savefile = None):
    ''' 
        Fetch an NCBI record with given id from given database (db). 
    '''
    
    # Can you guess why I have a while loop?
    no_response = True
    while no_response:
        try:
            # The rettype argument indicates that we want a GenBank file, and the retmode indicates that we want it returned as a plain text file
            fetch = Entrez.efetch(id=id, db=db, rettype="gb", retmode="text") 
            no_response = False
        except:
            pass
        record = SeqIO.read(fetch, "gb")
    
    if savefile:
        SeqIO.write(record, savefile, "genbank")
    
    return record
    

# http://www.ncbi.nlm.nih.gov/protein/NP_000549.1
hemoglobin_protein = fetch_ncbi_record("NP_000549.1", "protein", savefile = "hemo_prot.gb")

# Useful attributes! [use dir() to see available methods and attributes]
#print dir(hemoglobin_protein)
#print hemoglobin_protein.seq
#print hemoglobin_protein.description

annotations = hemoglobin_protein.annotations
#print annotations.keys()
#print annotations["db_source"]



# Capture the features attribute (which is a dictionary!)
# Each feature has three attributes: type, location, and qualifiers
# Generally, we're interested in qualifiers, which is itself also a dictionary!
features = hemoglobin_protein.features

# Examine all the entries
for entry in features:
    print entry

# Examine what's inside the CDS feature (as an example)
for entry in features:
    if entry.type == "CDS":
        print entry.qualifiers

# Grab information from the CDS feature
for entry in features:
    if entry.type == "CDS":
        for qual in entry.qualifiers:
            print qual, entry.qualifiers[ qual ]
        
        # Grab the CDS id
        codedby = entry.qualifiers['coded_by'][0]
        cds_id = codedby.split(":")[0]
        print "And here's the CDS id:", cds_id



# Now, we can fetch down this protein record's corresponding nucleotide record and grab whatever information we like from it.
hemoglobin_nuc = fetch_ncbi_record(cds_id, "nucleotide", savefile = "hemo_nuc.gb")

    