#CHEATSHEET PYTHON 6: Biopython


##Sequence Input/Output
To read and write sequence or alignment data, enter the lines `from Bio import SeqIO` or `from Bio import AlignIO`, respectively.

Command  |  Description | Example
----------|-------------|----------
`SeqIO.read()` | read a single sequence from a file | `SeqIO.read("seq.fasta", "fasta")
`SeqIO.parse()` | read many sequences from a file | `SeqIO.parse("seqs.fasta", "fasta")
`AlignIO.read()` | read a single alignment from a file | `SeqIO.read("align.fasta", "fasta")
`AlignIO.parse()` | read many alignments from a file | `SeqIO.parse("aligns.phy", "phylip")
`SeqIO.write()` | write sequence(s) to a file | `SeqIO.write(seq_record(s), "seq_output.fasta", "fasta")
`AlignIO.write()` | write alignment(s) to a file | `AlignIO.write(alignment, "align_output.fasta", "fasta")

##Biopython objects
To use and/or manipulate these objects, enter the lines `from Bio.Seq import Seq`, `from Bio.SeqRecord import SeqRecord`, `from Bio.Alphabet import *`

Object  |  Useful attributes/methods
--------|------------------
`SeqRecord` | `.id`, `.seq`, `.description`
`Seq` | `.transcribe()`, `.translate()`, `.complement()`, and more..

##Collect NCBI Data
To download from NCBI, enter the line `from Bio import Entrez`. To avoid annoying warnings, also enter the command, `Entrez.email = "myemail@email.com"` (replaced with your email) before fetching!

Command | Description | Example
--------|-------------|---------
`Entrez.efetch()` | Fetch a record from the NCBI database (the type of id is up to you!) | `Entrez.efetch(id = "NP_000549.1", db = "protein", remode = "text", retype = "gb")` <br> `Entrez.efetch(id = "88758587", db = "nucleotide", remode = "text", retype = "gb")`

