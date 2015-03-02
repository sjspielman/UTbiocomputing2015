# Python 5: Code Hygiene

###Good coding practices:

		Modular coding, i.e. make each function do one thing and do it well. It’s easier to compose and test minimal, single-purpose functions.
	Test a function before you write another one!
	Comment a lot and add doc strings once you have a functioning function. Future you will be so happy.
	Use informative names for your variables. Function names should be verbs while instances should be nouns or noun phrases. In general, avoid short abbreviations.
	Once you're a seasoned python-er, keep your style up to standard by reading 'best practices' (eg [link](https://www.memonic.com/user/pneff/folder/python/id/1bufp))


###Ways to test your code:

	1. use a reduced file to test your code
	2. using print
	3. using assert with a boolean comparison
	4. try and except
	5. print a log file
	
	
###Once your code works:

	1. Clean it up! Delete or comment out old tests
	2. Add doc strings
	3. Test on more files.
	
	
###What about unit testing?
	
	Unless you're developing your own software, use other types of tests. We won't cover this today.



###Examples

	1. If you're working with large files, use a reduced file to test your code. 

```
head -1000 largeDataFile.csv > test.csv 
```

	If you plan to automate your analyses, test each function with only one file at first.

``` python
function_name('file1.txt')

#rather than:
for file in directory_list:
	function_name(file)
```

	2. As you're writing functions, use print to check that your output is what you expect.
	
```python
import os
from Bio.PDB.PDBParser import PDBParser
parser = PDBParser()

def make_filelist(directory, file_ending):

	"""This function makes a list of files with a certain ending in a certain directory.
	"""
	files=os.listdir(directory)
	files_pdb=[]
	for file in files:
		if file.endswith(file_ending):
			files_pdb.append(file)
	print files
	print files_pdb
	return files_pdb

make_filelist('.', '.pdb')
```
	
	If you're working with integers or floats or are getting a TypeError from python, you can use the type() function to check that the variable is in the format you expect.
	If you're having trouble indexing certain values, copy a subset into the python console and check your syntax.

	3. Python's assert function allows you to test a comparison and exits the script if not true.
	You can have it print an error message upon exiting.

```python
assert a == b, "Error: comparison %s == %s is false" %(value1,value2)
```

	Example: 
	
```python

num_list=[3,52,6,'b',2,463,'a']

def sum_num(num_list):
	totalSum=0
	for item in num_list:
		assert type(item)==int, "Uh oh, item '%s' not the right type! Exiting now." %item
		totalSum+=item
	print totalSum
			
sum_num(num_list)


```

	4. python's try and except allows you to trigger error messages for specific types of errors without killing the script.
	
```python
num_list=[3,52,6,'b',2,463,'a']

def sum_num(num_list):
	totalSum=0
	for item in num_list:
		try:
			totalSum+=item
		except TypeError:
			print "Warning, item '%s' not the right type! Continuing to next item" %item
	print totalSum
			
sum_num(num_list)


```

	
	5. Print stdout to a log file when using your computer to check for errors. Everything you 'print' will be concatenated to the logfile. 

```
python script.py > logfile
```
	
	However, you will need to print directly to a logfile if running on a cluster (in this case you will
	submit your file to the cluster rather than running python in the shell). Either open the logfile as 
	a global variable at the top of the file or open and close it every time you write to it.
	
```python
import time
logfile=open('jobname.log','a') 

def function_name(arguments):

	start=time.clock()
	code doing things
	and other things
	print >>logfile,"%s was parsed in %fs." %(input_file,(time.clock() - start))

function_name(arguments)



import time

def function_name(arguments):

	start=time.clock()
	code doing things
	and other things
	with open('jobname.log','a') as f:
		f.write("%s was parsed in %fs." %(input_file,(time.clock() - start)))

function_name(arguments)
```


###Now we will go through what modular coding and testing should look like with an example.
See the files for [python5](python5_files.zip)



	
