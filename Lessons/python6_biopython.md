# Python VI: Biopython

An epic biopython [how-to guide](http://biopython.org/DIST/docs/tutorial/Tutorial.html) describes nearly all the capabilities (with examples!) of things to do with Biopython (ctrl+F is your friend here!). While this library has lots of functionality, it is primarily useful for dealing with sequence data.

## Reading and parsing sequence files

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

## Manipulating and Writing Biopython objects

Working with sequence files can be a bit tricky since all sequences are really Biopython SeqRecord objects. Most sequence processing is done with the sequences converted to *strings*, but to interface with Biopython,

```python
>>> from Bio.Seq import Seq
>>> from Bio.SeqRecord import SeqRecord
>>> from Bio.Alphabet import *

>>> # Creating a Seq and SeqRecord object
>>> my_sequence = "ACGTACCGTTTTGGAACTTCC"

>>> # Recast to Seq object. Two arguments are the sequence and the alphabet (the latter is optional!!)
>>> my_biopython_seq = Seq(my_sequence, generic_dna)

>>> # Some exciting Biopython methods! All results are biopython objects which can be re-cast into strings with str()
>>> my_biopython_seq.complement()
Seq('TGCATGGCAAAACCTTGAAGG', DNAAlphabet())
>>> my_biopython_seq.reverse_complement()
Seq('GGAAGTTCCAAAACGGTACGT', DNAAlphabet())
>>> my_biopython_seq.transcribe()
Seq('ACGUACCGUUUUGGAACUUCC', RNAAlphabet())
>>> my_biopython_seq.translate()
Seq('TYRFGTS', ExtendedIUPACProtein())

>>> # Like strings, biopython seq objects are *immutable*
>>> my_biopython_seq[2] = "A"
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: 'Seq' object does not support item assignment


>>> # Create a SeqRecord object
>>> my_biopython_seqrec = SeqRecord(my_biopython_seq, id = "seq1") # Can add other arguments as desired
>>> my_biopython_seqrec
SeqRecord(seq=Seq('ACGTACCGTTTTGGAACTTCC', DNAAlphabet()), id='seq1', name='<unknown name>', description='<unknown description>', dbxrefs=[])



>>> # Biopython can write SeqRecord objects to files. Arguments are SeqRecord object(s), file name, sequence format
>>> SeqIO.write(my_biopython_seqrec, "newseq.fasta", "fasta") 

>>> # Write multiple sequences to file by providing SeqIO with a list of SeqRecord objects
>>> # Note that AlignIO.write works just like this!
>>> rec1 = SeqRecord(my_biopython_seq, id = "seq1")
>>> rec2 = SeqRecord(my_biopython_seq, id = "seq2")
>>> rec3 = SeqRecord(my_biopython_seq, id = "seq3")
>>> lots_of_records = [rec1, rec2, rec3]
>>> SeqIO.write(lots_of_records, "lots_of_records.phy", "phylip")
```


## Converting file formats
Biopython makes it very straight-forward to convert between sequence file formats - simply read in a file and write it out in the new format, or use the handy .convert() method. Again, these methods work with both `AlignIO()` and `SeqIO()`.

```python
>>> # Change file formats
>>> temp = AlignIO.read("file.fasta", "fasta")
>>> AlignIO.write(temp, "newfile.phy", "phylip") #object to write, filename, format 

>>> # Alternatively...
>>> AlignIO.convert("file.fasta", "fasta", "newfile.phy", "phylip")
```

## Interacting with sequence alignments

```python
>>> from Bio import AlignIO
>>> aln = AlignIO.read("pb2.fasta", "fasta")
>>> # Use len() to determine alignment size
>>> number_sequences = len(aln)
>>> print number_sequences
400
>>> number_columns = len(aln[0])
>>> print number_columns
2277

>>> # Extract alignment positions
>>> row5_column10_annoying = aln[5].seq[10] 
>>> row5_column10_easier   = aln[5,10]

>>> # Extract alignment columns
>>> col4 = aln[:,4]
'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'

>>> Extract alignment regions
>>> aln_chunk = aln[5:10, 9:13] # row indices 5-10 and column indices 9-13
>>> print aln_chunk
SingleLetterAlphabet() alignment with 5 rows and 4 columns
ATAA Av_ABC66491
ATAA Av_CAF33010
ATAA Av_ABI84837
ATAA Av_ABI85028
ATAA Av_ABB19754
```

Creating new sequence alignments is a bit tricky - you need to create multiple SeqRecord objects and merge them together into a MultipleSequenceAlignment object. We won't do that here, but you can look it up on the Biopython tutorial linked above.

## Useful function!

Here's a quick function that concatenates the description of identical sequences and prints them to a new file:
```python
# dictionary format will be  sequence:[description1, description2, ...]
def sequence_cleaner(fasta_file):
	sequences={}
 	for seq_record in SeqIO.parse(fasta_file, "fasta"):
 		sequence=str(seq_record.seq).upper()
    	if sequence not in sequences:
    		sequences[sequence]=seq_record.description 
    	else:
    		sequences[sequence]+="\t"+seq_record.description

	output_file=open("clean_"+fasta_file,"w")
	for sequence in sequences:
		output_file.write(">"+sequences[sequence]+"\n"+sequence+"\n")
	output_file.close()

sequence_cleaner(infile)
```

## Biopython Entrez

Biopython has excellent built-in tools for collecting from Entrez and SwissProt/ExPASy (see section "Parsing sequences from the net" in the how-to guide linked above). An example of fetching and parsing NCBI sequence data is available in the script [entrez.py](python6_files/entrez.py).






