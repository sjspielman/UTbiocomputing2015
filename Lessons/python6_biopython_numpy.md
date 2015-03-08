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


#### Biopython Entrez

Biopython has excellent built-in tools for collecting from Entrez and SwissProt/ExPASy (see section "Parsing sequences from the net" in the how-to guide linked above). An example of fetching and parsing NCBI sequence data is available in the script [entrez.py](python6_files/entrez.py).

## NumPy

The NumPy module is primarily useful for dealing with arrays and matrices. If you are familiar with MatLab, NumPy is (more or a less) a python substitute for much of MatLab's functionality. (Sidenote - there is a separate module called SciPy, often bundled with NumPy, which is also quite useful. Check out the docs!)

#### Creating and manipulating arrays

Unlike python lists, NumPy arrays are more similar to arrays in C/C++. All entries in an array must be of the same type, and the size of the array is predetermined (you cannot append straight-forwardly!). However, these arrays are much more powerful. They take up far less memory and are great for performing vectorized functions.

```python
>>> import numpy as np

>>> # Array creation routines
>>> z = np.zeros(10)
>>> oh = np.ones(10)
>>> my_array = np.array( [5, 6, 7, 8, 9.5])
>>> my_array / 5
array([ 1. ,  1.2,  1.4,  1.6,  1.9]) # Note how they are all treated as floats, since the 9.5 was included and all array elements are of the same type!

>>> # Array masking returns a *boolean* array
>>> z == 0
array([ True,  True,  True,  True,  True,  True,  True,  True,  True,  True], dtype=bool)





#### Indexing and slicing

#### Reading and writing data











