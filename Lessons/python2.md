#Python 2: Loops!!



First, more useful functions

- `dir()` shows all methods you can use with a variable

```python
>>> name="Rebecca"
>>> dir(name)
['__add__', '__class__', '__contains__', '__delattr__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__getitem__', '__getnewargs__', '__getslice__', '__gt__', '__hash__', '__init__', '__le__', '__len__', '__lt__', '__mod__', '__mul__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__rmod__', '__rmul__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '_formatter_field_name_split', '_formatter_parser', 
'capitalize', 'center', 'count', 'decode', 'encode', 'endswith', 'expandtabs', 'find', 'format', 'index', 'isalnum', 'isalpha', 'isdigit', 'islower', 'isspace', 'istitle', 'isupper', 'join', 'ljust', 'lower', 'lstrip', 'partition', 'replace', 'rfind', 'rindex', 'rjust', 'rpartition', 'rsplit', 'rstrip', 'split', 'splitlines', 'startswith', 'strip', 'swapcase', 'title', 'translate', 'upper', 'zfill']

>>> mydict={'mole': 4, 'human': 7, 'mule': 5}
>>> dir(mydict)
['__class__', '__cmp__', '__contains__', '__delattr__', '__delitem__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__getitem__', '__gt__', '__hash__', '__init__', '__iter__', '__le__', '__len__', '__lt__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__setitem__', '__sizeof__', '__str__', '__subclasshook__', 
'clear', 'copy', 'fromkeys', 'get', 'has_key', 'items', 'iteritems', 'iterkeys', 'itervalues', 'keys', 'pop', 'popitem', 'setdefault', 'update', 'values', 'viewitems', 'viewkeys', 'viewvalues']

```

- using `%` notation in print statements makes for better output

	%d: insert number here <br>
	%f: insert float here <br>
	%s: insert string here <br>

```python
>>> print "My name is %s." % (name)
My name is Rebecca.

>>> print "Two divided by 5 is %f" %(2.0/5)
Two divided by 5 is 0.400000
>>> print "Two divided by 5 is %0.1f" %(2.0/5)
Two divided by 5 is 0.4
>>> print "Two divided by 5 is %d" %(2.0/5)
Two divided by 5 is 0
```

- `random` module is useful for generating fake data

```python
>>> import random
>>> dna = ["A","C","T","G"]
>>> dna=dna*250
>>> random.shuffle(dna)
>>> dna
['A', 'A', 'A', 'A', 'T', 'C', 'T', 'A', 'T', 'G', 'G', 'C', 'A', 'C', 'C', 'G', 'C', 'A', 'C', 'T', 'A', 'A', 'T', 'G', 'C', 'G', 'A', 'G', 'G', 'T', 'C', 'G', 'A', 'C', 'A', 'G', 'A', 'G', 'T', 'A', 'A', 'G', 'A', 'C', 'T', 'T', 'T', 'G', 'G', 'T', 'G', 'T', 'G', 'T', 'C', 'G', 'T', 'G', 'G', 'A', 'C', 'T', 'T', 'G', 'T', 'C', 'C', 'A', 'T', 'G', 'C', 'T', 'G', 'T', 'A', 'G', 'G', 'A', 'T', 'A', 'T', 'G', 'T', 'G', 'A', 'T', 'C', 'G', 'T', 'T', 'G', 'C', 'A', 'C', 'G', 'A', 'G', 'C', 'T', 'A', 'G', 'T', 'G', 'T', 'G', 'A', 'C', 'T', 'C', 'T', 'C', 'C', 'A', 'A', 'G', 'G', 'T', 'A', 'G', 'T', 'C', 'G', 'A', 'A', 'C', 'T', 'G', 'C', 'T', 'G', 'C', 'T', 'T', 'A', 'G', 'G', 'T', 'T', 'G', 'C', 'T', 'T', 'A', 'T', 'G', 'G', 'G', 'A', 'G', 'A', 'C', 'G', 'G', 'C', 'T', 'T', 'C', 'C', 'G', 'G', 'T', 'C', 'G', 'T', 'T', 'T', 'T', 'A', 'T', 'G', 'G', 'T', 'T', 'G', 'A', 'C', 'T', 'G', 'G', 'A', 'C', 'A', 'G', 'A', 'A', 'C', 'A', 'A', 'C', 'A', 'C', 'T', 'G', 'C', 'T', 'C', 'C', 'A', 'G', 'C', 'C', 'C', 'C', 'A', 'C', 'A', 'G', 'A', 'T', 'T', 'A', 'A', 'C', 'T', 'G', 'A', 'T', 'C', 'T', 'G', 'C', 'A', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'G', 'G', 'C', 'A', 'C', 'T', 'C', 'C', 'G', 'C', 'T', 'G', 'C', 'T', 'T', 'T', 'T', 'A', 'A', 'G', 'G', 'C', 'G', 'C', 'G', 'A', 'A', 'C', 'A', 'T', 'A', 'C', 'G', 'C', 'T', 'G', 'G', 'T', 'G', 'C', 'T', 'A', 'A', 'G', 'C', 'C', 'C', 'G', 'C', 'T', 'A', 'G', 'T', 'T', 'C', 'A', 'G', 'C', 'A', 'C', 'C', 'A', 'A', 'T', 'C', 'G', 'A', 'C', 'T', 'G', 'T', 'T', 'G', 'C', 'A', 'G', 'A', 'T', 'A', 'A', 'G', 'T', 'T', 'C', 'A', 'T', 'C', 'T', 'G', 'G', 'T', 'C', 'A', 'A', 'A', 'G', 'A', 'T', 'C', 'A', 'T', 'T', 'C', 'A', 'C', 'A', 'C', 'A', 'G', 'T', 'A', 'A', 'G', 'A', 'T', 'A', 'T', 'T', 'G', 'A', 'C', 'C', 'G', 'G', 'G', 'A', 'G', 'C', 'C', 'T', 'A', 'G', 'A', 'T', 'C', 'T', 'C', 'C', 'A', 'G', 'A', 'T', 'G', 'A', 'A', 'C', 'G', 'C', 'T', 'G', 'T', 'G', 'G', 'C', 'G', 'C', 'A', 'C', 'C', 'C', 'T', 'C', 'G', 'A', 'T', 'T', 'C', 'G', 'C', 'G', 'G', 'G', 'G', 'C', 'G', 'A', 'T', 'C', 'C', 'A', 'G', 'A', 'A', 'G', 'C', 'C', 'A', 'T', 'G', 'A', 'G', 'G', 'A', 'A', 'T', 'A', 'T', 'A', 'G', 'T', 'T', 'A', 'T', 'C', 'T', 'C', 'G', 'T', 'C', 'T', 'G', 'C', 'C', 'C', 'G', 'T', 'A', 'T', 'A', 'C', 'C', 'G', 'C', 'G', 'C', 'T', 'C', 'C', 'G', 'T', 'A', 'C', 'C', 'A', 'A', 'T', 'G', 'C', 'A', 'C', 'C', 'T', 'T', 'A', 'A', 'T', 'A', 'T', 'C', 'G', 'A', 'G', 'G', 'C', 'C', 'G', 'G', 'A', 'C', 'A', 'C', 'G', 'A', 'G', 'A', 'C', 'G', 'T', 'A', 'C', 'T', 'T', 'C', 'A', 'T', 'T', 'T', 'G', 'A', 'C', 'C', 'A', 'C', 'G', 'C', 'G', 'G', 'T', 'A', 'G', 'C', 'A', 'G', 'G', 'T', 'T', 'C', 'G', 'T', 'G', 'G', 'C', 'A', 'T', 'A', 'T', 'G', 'A', 'G', 'C', 'C', 'A', 'G', 'A', 'T', 'A', 'G', 'G', 'A', 'G', 'T', 'G', 'T', 'G', 'G', 'C', 'T', 'A', 'T', 'C', 'C', 'C', 'C', 'A', 'A', 'G', 'G', 'T', 'C', 'C', 'T', 'A', 'C', 'A', 'G', 'G', 'T', 'C', 'A', 'C', 'T', 'A', 'G', 'A', 'C', 'T', 'T', 'A', 'G', 'C', 'T', 'A', 'G', 'T', 'G', 'T', 'T', 'G', 'C', 'A', 'G', 'T', 'T', 'G', 'G', 'G', 'T', 'G', 'G', 'G', 'C', 'A', 'C', 'A', 'A', 'T', 'G', 'A', 'C', 'A', 'C', 'C', 'C', 'G', 'G', 'C', 'A', 'G', 'G', 'A', 'G', 'A', 'G', 'C', 'A', 'T', 'T', 'T', 'G', 'G', 'C', 'C', 'G', 'C', 'G', 'T', 'G', 'T', 'T', 'T', 'C', 'G', 'A', 'C', 'A', 'T', 'T', 'C', 'C', 'C', 'T', 'T', 'T', 'A', 'A', 'T', 'G', 'A', 'A', 'A', 'G', 'A', 'C', 'T', 'C', 'G', 'C', 'C', 'A', 'G', 'C', 'T', 'C', 'T', 'A', 'A', 'T', 'A', 'C', 'C', 'A', 'C', 'T', 'C', 'A', 'A', 'C', 'T', 'T', 'T', 'T', 'A', 'C', 'C', 'T', 'C', 'T', 'T', 'T', 'G', 'G', 'A', 'A', 'G', 'A', 'C', 'T', 'A', 'T', 'T', 'A', 'T', 'T', 'A', 'A', 'A', 'G', 'A', 'T', 'G', 'C', 'C', 'T', 'A', 'C', 'T', 'C', 'A', 'A', 'G', 'G', 'T', 'G', 'G', 'A', 'T', 'T', 'T', 'C', 'G', 'A', 'C', 'A', 'C', 'G', 'G', 'T', 'G', 'C', 'G', 'G', 'T', 'C', 'A', 'T', 'G', 'T', 'A', 'A', 'T', 'C', 'T', 'C', 'G', 'A', 'C', 'A', 'C', 'C', 'A', 'T', 'G', 'G', 'G', 'T', 'T', 'G', 'C', 'A', 'G', 'T', 'T', 'G', 'C', 'G', 'T', 'T', 'C', 'A', 'G', 'T', 'A', 'A', 'C', 'T', 'C', 'A', 'A', 'G', 'A', 'C', 'T', 'C', 'T', 'G', 'C', 'C', 'A', 'T', 'T', 'C', 'G', 'A', 'C', 'T', 'G', 'G', 'A', 'C', 'A', 'C', 'A', 'G', 'A', 'T', 'A', 'A', 'G', 'A', 'T', 'G', 'G', 'T', 'C', 'A', 'T', 'A', 'C', 'G', 'C', 'T', 'C', 'G', 'G', 'T', 'A', 'C', 'C', 'G', 'A', 'T', 'G', 'A', 'G', 'T', 'T', 'C', 'C', 'T', 'A', 'A', 'T', 'C', 'T', 'A', 'T', 'C', 'G', 'G', 'A', 'A', 'C', 'T', 'A', 'A', 'A', 'G', 'C', 'T', 'C', 'G', 'A', 'A', 'C', 'G', 'T', 'A', 'T', 'C', 'G', 'C', 'A', 'T', 'A', 'C', 'C', 'A', 'C', 'C', 'A', 'T', 'G', 'C', 'C', 'G', 'T', 'T', 'G', 'G', 'T', 'G', 'C', 'T', 'A', 'A', 'C', 'G', 'A', 'A', 'C', 'T', 'T', 'C', 'T', 'C', 'G', 'G', 'C', 'G', 'C', 'C', 'T', 'G', 'A', 'T', 'A', 'G', 'A', 'C', 'A', 'G', 'G', 'T', 'T', 'G', 'G', 'A', 'A', 'T', 'A', 'C', 'T', 'G', 'C', 'T', 'C', 'G', 'G', 'C', 'G', 'A', 'A', 'T', 'A', 'T', 'C', 'A', 'T', 'C', 'A', 'A', 'T', 'A', 'C', 'A', 'A', 'T', 'A', 'A', 'A', 'C', 'C', 'G', 'T', 'G', 'C', 'G', 'T', 'T', 'C']
>>> "".join(dna)
'AAAATCTATGGCACCGCACTAATGCGAGGTCGACAGAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGATTAACTGATCTGCAACGTACGTGGCACTCCGCTGCTTTTAAGGCGCGAACATACGCTGGTGCTAAGCCCGCTAGTTCAGCACCAATCGACTGTTGCAGATAAGTTCATCTGGTCAAAGATCATTCACACAGTAAGATATTGACCGGGAGCCTAGATCTCCAGATGAACGCTGTGGCGCACCCTCGATTCGCGGGGCGATCCAGAAGCCATGAGGAATATAGTTATCTCGTCTGCCCGTATACCGCGCTCCGTACCAATGCACCTTAATATCGAGGCCGGACACGAGACGTACTTCATTTGACCACGCGGTAGCAGGTTCGTGGCATATGAGCCAGATAGGAGTGTGGCTATCCCCAAGGTCCTACAGGTCACTAGACTTAGCTAGTGTTGCAGTTGGGTGGGCACAATGACACCCGGCAGGAGAGCATTTGGCCGCGTGTTTCGACATTCCCTTTAATGAAAGACTCGCCAGCTCTAATACCACTCAACTTTTACCTCTTTGGAAGACTATTATTAAAGATGCCTACTCAAGGTGGATTTCGACACGGTGCGGTCATGTAATCTCGACACCATGGGTTGCAGTTGCGTTCAGTAACTCAAGACTCTGCCATTCGACTGGACACAGATAAGATGGTCATACGCTCGGTACCGATGAGTTCCTAATCTATCGGAACTAAAGCTCGAACGTATCGCATACCACCATGCCGTTGGTGCTAACGAACTTCTCGGCGCCTGATAGACAGGTTGGAATACTGCTCGGCGAATATCATCAATACAATAAACCGTGCGTTC'
>>> 
```

#Decisions and conditions: For, If, While


##For loops

A for loop is a way to iterate through a list, range, dictionary, string, etc. Use it for whenever you are repeating something.

```
for item in thing:
	do this command
	do that command
	#return to for statement and move to next item
continue with main commands
```

Ok, let's test this out with the print statement. I use this statement in for loops often to check for errors and make sure that everything is working as I would expect.

```python
>>> dna='AGCT'
>>> for item in dna:
...     print item
... 
A
G
C
T
```

Let's try something a bit more complicated. Here we can do a repeat of a mathematical operation on items in a list.

```python
>>>sequences=['AGTCTA','AGTCAGTCAGTCAGT','ACTAGCTAGCTA','ACGTCAGTATCGTATTTTA','ACAGTCAGTGATCA','AGT','AGCTAGCTAGCTACGATGCTAGCTAGC']
>>> for i in sequences:
...     length=len(i)
...     Gcontent=i.count('G')
...     Ccontent=i.count('C')
...     Tcontent=i.count('T')
...     Acontent=i.count('A')
...     GCcontent=(Gcontent+Ccontent)/float(length)
...     print "GC content of %s is %f" %(i,GCcontent)
... 
GC content of AGTCTA is 0.333333
GC content of AGTCAGTCAGTCAGT is 0.466667
GC content of ACTAGCTAGCTA is 0.416667
GC content of ACGTCAGTATCGTATTTTA is 0.315789
GC content of ACAGTCAGTGATCA is 0.428571
GC content of AGT is 0.333333
GC content of AGCTAGCTAGCTACGATGCTAGCTAGC is 0.518519
```

But, if you see in this code, there is a method that is repeated four times for different value, so we should also use a for loop in this case.

```python
>>>GCdict={} #create an empty dictionary that we will fill with GCcontent values
>>>for seq in sequences:
...     seqdict={} #create another empty dictionary to *temporarily* store ACGT counts
...     length=len(seq)
...     for nuc in ['A','G','C','T']: #for each nucleotide
...             seqdict[nuc]=seq.count(nuc) #count the number of nucleotides in the sequence
...     GCcontent=(seqdict['G']+seqdict['C'])/float(length) #what happens if you don't do float()?
...		GCdict[seq]=GCcontent
...     print "GC content of %s is %f" %(seq,GCcontent)
... 
GC content of AGTCTA is 0.333333
GC content of AGTCAGTCAGTCAGT is 0.466667
GC content of ACTAGCTAGCTA is 0.416667
GC content of ACGTCAGTATCGTATTTTA is 0.315789
GC content of ACAGTCAGTGATCA is 0.428571
GC content of AGT is 0.333333
GC content of AGCTAGCTAGCTACGATGCTAGCTAGC is 0.518519
```

Now we can check what is in our dictionary.

```python
>>> for key in GCdict:
...     print '%s: %f' %(key,GCdict[key])
... 
ACAGTCAGTGATCA: 0.428571
ACGTCAGTATCGTATTTTA: 0.315789
ACTAGCTAGCTA: 0.416667
AGTCAGTCAGTCAGT: 0.466667
AGCTAGCTAGCTACGATGCTAGCTAGC: 0.518519
AGTCTA: 0.333333
AGT: 0.333333

```

A quick note: you can do for loops in bash too!

```
#for file in a list of text files in this directory, print filename
for i in *.txt; do echo $i; done
for i in *.fasta; cat i >> newseqs.fasta; done
```

And in many many other languages: http://rosettacode.org/wiki/Foreach


##if statements

- IF

If statements are used to test if your variable is true for some logical condition in order to determine what the program should do next. These are super handy and are essentially how you get the computer to make decisions for you.

```
if logical condition == True: #line leading to an indented block ends with a colon
	do this command  #indenting blocks of code indicates how the program flows
	do this command
continue with main commands
```

```python 
>>> numbers=range(20)
>>> numbers
[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
>>> if sum(numbers) >= 20:
...     print "Sum is greater than 20!"
... 
Sum is greater than 20!

>>> if sum(numbers) >= 20:
...     print "%d is greater than 20!" % (sum(numbers))
... 
190 is greater than 20!
>>> 

```

- IF ELSE

Adding an else to your statement allows the computer to decide what to do based on the conditional.

```
if logical condition == True: 
	do this command  
	do this command
else:
	do this command
```


```python
>>> dna='AGTCTGTAGTCTATAGA'
>>> if (dna.count('G') + dna.count('C'))/len(dna) > 0.5:
...     print "GC content higher than 50%"
... else:
...     print "GC content lower than 50%"
... 
GC content lower than 50%
```

You can do nested if/else statements.

```
if logical condition == True: 
	do this command  
	do this command
else:
	if logical condition == True: 
		do this command  
		do this command
	else:
		do this command
```

```python
>>> for item in newseqs:
...     if 'TTTTTT' in item:
...             print 'TTTTTT is in %s' %item
...     else:
...             if 'TTTTT' in item:
...                     print 'TTTTT is in %s' %item
...             else:
...                     if 'TTTT' in item:
...                             print 'TTTT is in %s' %item
...                     else:
...                             if 'TTT' in item:
...                                     print 'TTT is in %s' %item
...                             else:
...                                     if 'TT' in item:
...                                             print 'TT is in %s' %item
...                                     else:
...                                             print 'no T repeats are in %s' %item
... 
TTTT is in GAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGATTAACTGA
TTTT is in CGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGAC
no T repeats are in GCACTCCGCTG
TTTT is in TGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAA
TTTT is in TGCGAGGTCGACAGAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGATTAACTGATCTGCAACGTACGT
TTTT is in TAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAG
TTT is in CCGCACTAATGCGAGGTCGACAGAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATA
TT is in TGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCT
TTTT is in ACACTGCTCCAGCCCCACAGATTAACTGATCTGCAACGTACGTGGCACTCCGCTGCTTTTAAG
TTTT is in CGACAGAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGT
...
```
But, wait, wasn't that a terribly long way to do something? Yes it was! Here is a better way:

- IF ELIF

The if/elif statement can be used when you would otherwise use multiple nested if/else statements. They are performed *in order*.


```
if logical condition == True: 
	do this command  
	do this command
elif other logical condition == True:
	do this command
elif other logical condition == True:
	do this command
```

```python
>>> for item in newseqs:
...     if 'TTTTTT' in item:
...             print 'TTTTTT is in %s' %item
...     elif 'TTTTT' in item:
...             print 'TTTTT is in %s' %item
...     elif 'TTTT' in item:
...             print 'TTTT is in %s' %item
...     elif 'TTT' in item:
...             print 'TTT is in %s' %item
...     elif 'TT' in item:
...             print 'TT is in %s' %item
...     else:
...             print 'no T repeats are in %s' %item
... 
TTTT is in GAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGATTAACTGA
TTTT is in CGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGAC
no T repeats are in GCACTCCGCTG
TTTT is in TGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAA
TTTT is in TGCGAGGTCGACAGAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGATTAACTGATCTGCAACGTACGT
TTTT is in TAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAG
...
```

- IF and FOR together!

Using conditionals and for loops together makes your code very powerful.

```python
>>> newseqs=[]
>>> for i in range(100):
...     int1=random.randint(0,250)
...     int2=random.randint(0,250)
...     if int1 > int2:
...             newseqs.append(newdna[int2:int1])
...     else:
...             newseqs.append(newdna[int1:int2])
... 
>>> newseqs
['GAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGATTAACTGA', 'CGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGAC', 'GCACTCCGCTG', 'TGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAA', 'TGCGAGGTCGACAGAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGATTAACTGATCTGCAACGTACGT', 'TAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAG', 'CCGCACTAATGCGAGGTCGACAGAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATA', 'TGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCT', 'ACACTGCTCCAGCCCCACAGATTAACTGATCTGCAACGTACGTGGCACTCCGCTGCTTTTAAG', 'CGACAGAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGT', 'TTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTG', 'TTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGATTAACTGATCTGCAACGTACGTGGCACTCC', 'AAATCTATGGCACCGCACTAATGCGAGGTCGACAGAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAAC', 'AGTAAGACTTTGGTG', 'CATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCC', 'GTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAA', 'ACTCCGCTGCTTTTAAG', 'TTTTA', 'GAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCC', 'TTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGG', 'GAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGC', 'TGGACAGAACAACACTGCTCCAGCCCCACAGATTAA', 'GCACTA', 'TTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGATTAACTGATCTGCAACGTAC', 'GATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCC', 'GGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGAT', 'GTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGA', 'CTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGATTAACTGATCT', 'TATGGCACCGCACTAATGCGAGGTCGACAGAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGG', 'GCACTAATGCGAGGTCGACAGAGTAAGA', 'GTAGG', 'CTAATGCGAGGTCGACAGAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGAT', 'CGCACTAATGCGA', 'ATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGG', 'CTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGA', 'GTCGTGGACT', 'GTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTAT', 'GTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGATTA', 'ATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGATTAACTGATCTGCAACGTACG', 'CTGGACAGAACAACACTGCTCCAGCCCCACAGATTAACTGATCTGCAACGTACGTGGCAC', 'CGCACTAATGCGAGGTCGACAGAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGG', 'TAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGC', 'CGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGATTAACTGATCTGCAACGTACGTGGCACTCCGCTGCTTT', 'AGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACA', 'TGGACAGAACAACACTGCTCCAGCCCCACAGATTAACTGATCTGCAACGTACGTGGCA', 'TCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGATTAACTGATCTGCAACGTACGTGGCACTCCGCT', 'GCACCGCACTAATGCGAGGTCGACAGAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGATTAACTGATCTGCAACGTACGTGGCACTCCGCTGCTTTT', 'CGAGCTAGTGTGACTCTCCA', 'AAAATCTATGGC', '', 'TCTATGGCACCGCACTAATGCGAGGTCGACAGAGT', 'ACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTT', 'CCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCA', 'CTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGATTAACTGATCTGCA', 'TGCGAGGTCGACAGAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCA', 'TGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGAT', 'AATCTAT', 'CTGGACAGAACAACACTGCTCCAGCCCCACAGAT', 'TAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGATTAACTGATCTGCAACGTACG', 'AGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTC', 'AGTAAGA', 'GTTGCTTATGGGAGA', 'CAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGA', 'TCGTTGCACGAGCTAGTGTGACTCTCCA', 'CTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGA', 'GTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGC', 'GAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGATTAACTGATCTGCAACGTACGTGG', 'CATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGG', 'GACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCT', 'GACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCT', 'GTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGG', 'GAGGTCGACAGAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGATTAAC', 'TGGACAGAACAACACTGCTCCAGC', 'TATGGCACCGCACTAATGCGAGGTCGACAGAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAG', 'ATCTATGGCACCGCACTAATGCGAGGTCGACAGAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGA', '', 'TGCTTAGGTTGCTTATGGGAGAC', 'GGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGATTAA', 'CTAATGCGAGGTCGACAGAGTAAGACTTTGGTGTGTCG', 'CGCACTAATGCGAGGTCGACAGAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCG', 'TTGCACGAGCTAGTGTGACTCTCCAAGGTAGT', 'GTTGACTGGACAGAACAACACTGCTCCAGCCCCAC', 'GACAGAACAACACTGC', 'CCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCC', 'CGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACG', 'TGGCACCGCACTAATGCGAGGTCGACAGAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCG', 'CTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCT', 'GGACTTGTCCATGCTGTAGGATATGTGAT', 'AGATTAACTGATCT', 'ATCTATGGCACCGCACTAATGCGAGGTCGA', 'ATATGTGATCGTTGCACGAGCTAGTGTGACTC', 'CACTAATGCGAGGTCGACAGAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAGAACAACACTGCTCCAGCCCCACAGATTAACTGATCTGCAACGTACGTGG', 'GATTAACTGATCTGCAA', 'ATGGCACCGCACTAATGCGAGGTCGACAGAGTA', 'ACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGGCTTCCGGTCGTTTTATGGTTGACTGGACAG', 'ATGGCACCGCACTAATGCGAGGTCGACAGAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGT', 'TAGTGTGACTCTCCAAGGTAG', 'CACTAATGCGAGGTCGACAGAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGA', 'GCGAGGTCGACAGAGTAAGACTTTGGTGTGTCGTGGACTTGTCCATGCTGTAGGATATGTGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGC', 'TGATCGTTGCACGAGCTAGTGTGACTCTCCAAGGTAGTCGAACTGCTGCTTAGGTTGCTTATGGGAGACGG']

>>> for seq in newseqs:
...     if len(seq) < 10:
...             print seq
...             newseqs.remove(seq)
... 
TTTTA
GCACTA
GTAGG
AATCTAT
AGTAAGA

>>> len(newseqs)
93
```

- Important: You can make a dictionary without knowing/defining a priori what they keys are

```python
>>> dictionary={}
>>> seq='AGTAGTATTTGAGAGTTTAGATATAG'
>>> for letter in seq:
...     if letter in dictionary.keys():
...             dictionary[letter]+=1
...     else:
...             dictionary[letter]=1
... 
>>> dictionary
{'A': 9, 'T': 10, 'G': 7}
```



##While

While loops repeat until some condition is True. Be careful as they are prone to infinite loops. If you get caught in one, hit Ctrl+C.

```
while condition == True:
	do this command
	do that command
continue with normal commands
```

The following script
```python
import time #useful module for printing timestamps and counting time units during your scripts
i = 0
while i < 5:
        time.sleep(1) #wait one second before proceeding
        i+=1
        print i
print "Time is up!"
```
outputs:
```python
python timer.py
1
2
3
4
5
Time is up!
```

While loops can ensure that user input is in the right format. They are also useful when you want to manipulate part of the conditional that is being used within the loop.

Try using this program, which has a test for user input using a while loop:
```python
x=int(raw_input("Give me a number from 1 to 10: ")) #accepts user input, converts to an integer
while x not in range(0,11): #if input was not number between 1 and 10
        print "That's not a number between 1 and 10 :(" #prints a complaint
        x=int(raw_input("Give me a number from 1 to 10: ")) #asks for a new number
print "Thanks!" #thanks the user for being a nice user
```

Here we manipulate the conditional (x) within the loop. Could you do this with a for loop?
```python 
>>> x=2000
>>> while x > 40:
...     x=x/4
...     if x%2 == 0:
...             print "%d is even!" %x
...     else:
...             print "%d is not even!" %x
...             x+=5
... 
500 is even!
125 is not even!
32 is even!
```


##Comprehensions

Comprehensions are faster ways to execute loops. Your typical loop structure is shown below followed by the comprehension structure.

```
list1=[] #make empty list
for item in thing:
	list1.append(item)	

list1 = [item for item in thing]
```

examples:

```python
>>> mylist=[]
>>> for i in range(0,20):
...     mylist.append(i**3)
... 
>>> mylist
[0, 1, 8, 27, 64, 125, 216, 343, 512, 729, 1000, 1331, 1728, 2197, 2744, 3375, 4096, 4913, 5832, 6859]

>>> mylist=[i**3 for i in range(0,20)]
>>> mylist
[0, 1, 8, 27, 64, 125, 216, 343, 512, 729, 1000, 1331, 1728, 2197, 2744, 3375, 4096, 4913, 5832, 6859]
```

```python
>>> sequences=['AGTCTA','AGTCAGTCAGTCAGT','ACTAGCTAGCTA','ACGTCAGTATCGTATTTTA','ACAGTCAGTGATCA','AGT','AGCTAGCTAGCTACGATGCTAGCTAGC']
>>> mylist=[len(seq) for seq in sequences]
>>> mylist
[6, 15, 12, 19, 14, 3, 27]

>>> mylist=[(len(seq),seq) for seq in sequences]
>>> mylist
[(6, 'AGTCTA'), (15, 'AGTCAGTCAGTCAGT'), (12, 'ACTAGCTAGCTA'), (19, 'ACGTCAGTATCGTATTTTA'), (14, 'ACAGTCAGTGATCA'), (3, 'AGT'), (27, 'AGCTAGCTAGCTACGATGCTAGCTAGC')]

>>> dictionary=dict(mylist)
>>> dictionary
{3: 'AGT', 6: 'AGTCTA', 12: 'ACTAGCTAGCTA', 14: 'ACAGTCAGTGATCA', 15: 'AGTCAGTCAGTCAGT', 19: 'ACGTCAGTATCGTATTTTA', 27: 'AGCTAGCTAGCTACGATGCTAGCTAGC'}
>>> 

```










