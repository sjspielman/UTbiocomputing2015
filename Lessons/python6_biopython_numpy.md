# Python VI: Biopython and NumPy

The Biopython and NumPy packages are both extremely useful libraries. An epic biopython [how-to guide](http://biopython.org/DIST/docs/tutorial/Tutorial.html) describes nearly all the capabilities (with examples!) of things to do with Biopython (ctl+F is your friend here!). The NumPy package has excellent [documentation](http://docs.scipy.org/doc/numpy/reference/), and the internet is chock-full of tutorials, so google away!


## Biopython

#### Reading and parsing sequence files

Biopython is an ideal tool for reading and writing sequence data. Biopython has two main modules for this purpose: `SeqIO` (sequence input-output) and `AlignIO` (alignment input-output). Each of these modules has two primary (although there are others!) functions: `read` and `parse`. The `read` function will read in a file with a single sequence or alignment, and the `parse` function will read in a file with multiple sequences or multiple sequence alignments. Each function takes two arguments: the file name and the sequence format (e.g. fasta, phylip, etc.)

```python
>>> from Bio import SeqIO, AlignIO

>>> # Read in a file of un-aligned sequences. Turn into list with list() (otherwise remains generator, which is fine, but you can't index)
>>> seqs = list( SeqIO.parse("seqs.fasta", "fasta") )

>>> # Read in an alignment file
>>> align = AlignIO.read("pb2.phy", "phylip-relaxed")
```
Biopython parses sequence files into, you guessed it, Biopython objects! Each sequence has several attributes (which you can examine with `dir()`), but the most important ones are `.seq`, `.id`, and `.description`.

```python
>>> # Biopython parses the file into a list of SeqRecord objects
>>> seqs = list( SeqIO.parse("seqs.fasta", "fasta") )
>>> print seqs 
[SeqRecord(seq=Seq('AGCTAGATCGATGC', SingleLetterAlphabet()), id='1', name='1', description='1', dbxrefs=[]), SeqRecord(seq=Seq('ATCGATACA', SingleLetterAlphabet()), id='2', name='2', description='2', dbxrefs=[]), SeqRecord(seq=Seq('ATACGAATAGCCTATACGTAGCATGCATGGGCTATAATTTTTT', SingleLetterAlphabet()), id='3', name='3', description='3', dbxrefs=[])]

>>> # Loop over sequences view important attributes, which can be converted to strings
>>> for record in seqs:
...   print "Record id:", record.id
...   print "Record sequence:", str(record.seq)

Record id: 1
The sequence record: AGCTAGATCGATGC
Record id: 2
The sequence record: ATCGATACA
Record id: 3
The sequence record: ATACGAATAGCCTATACGTAGCATGCATGGGCTATAATTTTTT
```

#### Manipulating and Writing Biopython objects

Working with sequence files can be a bit tricky since all sequences are really Biopython SeqRecord objects. Most sequence processing is done with the sequences converted to *strings*, but to interface with Biopython,

```python
>>> from Bio.Seq import Seq
>>> from Bio.SeqRecord import SeqRecord
>>> from Bio.Alphabet import *

>>> # Creating a Seq and SeqRecord object
>>> my_sequence = "ACGTACCGTTTTGGAACTTCC"

>>> # Recast to Seq object. Two arguments are the sequence and the alphabet (the latter is optional!!)
>>> my_biopython_seq = Seq(my_sequence, generic_dna)
>>> translated_seq = my_biopython_seq.translate() #You can convert this back to a string as desired
Seq('TYRFGTS', ExtendedIUPACProtein())

>>> # Create a SeqRecord object
>>> my_biopython_seqrec = SeqRecord(my_biopython_seq, id = "seq1") # Can add other arguments as desired
>>> my_biopython_seqrec
SeqRecord(seq=Seq('ACGTACCGTTTTGGAACTTCC', DNAAlphabet()), id='seq1', name='<unknown name>', description='<unknown description>', dbxrefs=[])

>>> # Biopython can write SeqRecord objects to files. Arguments are SeqRecord object(s), file name, sequence format
>>> SeqIO.write(my_biopython_seqrec, "newseq.fasta", "fasta") 

>>> # Write multiple sequences to file by providing SeqIO with a list of SeqRecord objects
>>> rec1 = SeqRecord(my_biopython_seq, id = "seq1")
>>> rec2 = SeqRecord(my_biopython_seq, id = "seq2")
>>> rec3 = SeqRecord(my_biopython_seq, id = "seq3")
>>> lots_of_records = [rec1, rec2, rec3]
>>> SeqIO.write(lots_of_records, "lots_of_records.phy", "phylip")
```

Note that `AlignIO.write()` works exactly as `SeqIO.write()`. 

#### Converting file formats
Biopython makes it very straight-forward to convert between sequence file formats - simply read in a file and write it out in the new format!

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











