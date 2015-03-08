# Python VI: Biopython and NumPy

Both of these are very big and very useful packages.
Epic biopython [how-to guide](http://biopython.org/DIST/docs/tutorial/Tutorial.html). Here, ctl+F is your friend!

NumPy tutorials are all over the place, and the package has really good documentation.

This lesson introduces some basic usage of each of these modules. Note that there's also something called SciPy which is great for other scientific things, but we won't cover here. Just note that it's commonly packaged w/ NumPy.
If you know anyone in MatLab, yell at them and tell them to switch to NumPy+SciPy. It's better.




## Biopython

#### Reading, writing, and parsing sequence files

```python
from Bio import SeqIO, AlignIO
# There are two functions for each: read and parse. read reads a single sequence (SeqIO) or alignment (AligIO), and parse reads multiple

# Read in a file of un-aligned sequences. Turn into list with list() (otherwise remains generator, which is fine, but you can't index)
seqs = list( SeqIO.read("sequences.fasta", "fasta") )
align = AlignIO.read("alignment.phy", "phylip-relaxed")
# there's also one of these: SeqIO.to_dict. I've never used it, but go crazy.

# useful attributes include seq, id, description. Remember though, that seq is a sequence object, so to muck with it you have to convert to string.

```


```python
# Change file formats
temp = AlignIO.read("file.fasta", "fasta")
AlignIO.write(temp, "newfile.phy", "phylip") #object to write, filename, format 
```


#### BioPython Entrez

Biopython has built-in tools for collecting from Entrez and SwissProt/ExPASy ("Parsing sequences from the net")

```python
from Bio import ExPASy, SwissProt
from Bio import Entrez
Entrez.email = stephanie.spielman@gmail.com # But please replace with your email!


>>> accessions = ["O23729", "O23730", "O23731"]
>>> records = []

>>> for accession in accessions:
...     handle = ExPASy.get_sprot_raw(accession)
        try:...     
            record = SwissProt.read(handle)
...         records.append(record)
        except ValueException:
            print "WARNING: Accession %s not found" % accession 


# A function I wrote a while back which is awesome
def fetch_ncbi_record(id, db, ncbi_dir):
    ''' Fetch an NCBI record if not already fetched and saved. '''
    
    if not os.path.exists(ncbi_dir + id + ".txt"):
        # This is a useful construction for querying internet in general. Sometimes servers don't respond, and you don't want to kill your program if that happens! You just want to keep at it until you ping it successfully.
        no_response = True
        while no_response:
            try:
                fetch = Entrez.efetch(db=db, rettype="gb", retmode="text", id=id)
                no_response = False
            except:
                pass
        record = SeqIO.read(fetch, "gb")
    else:
        record = SeqIO.read(ncbi_dir + id + ".txt", "gb")
    return record
    
# Horribly disgusting. Best strategy is to figure it out for your files.
feat = prot_record.features
    cds_id = ''
    for entry in feat:
        if entry.type=='CDS':
            codedby = entry.qualifiers['coded_by'][0]
            cds_id = codedby.split(":")[0]
            assert(CDS_MATCH.match(cds_id)), "No coding sequence id was found."
```


#get_prodoc_entry
#To download a Prosite documentation record in HTML format
#get_prosite_entry
#To download a Prosite record in HTML format
#get_prosite_raw
#To download a Prosite or Prosite documentation record in raw format
#get_sprot_raw
#To download a Swiss-Prot record in raw format
#sprot_search_ful
#To search for a Swiss-Prot record
#sprot_search_de
#To search for a Swiss-Prot record

## NumPy

#### Creating and manipulating arrays

#### Indexing and slicing

#### Reading and writing data











