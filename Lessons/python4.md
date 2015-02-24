# Python IV: File Input/Output

## Opening and closing files
To interact with files in python, you must open the file and assign it to a variable - this is now a python file object. All operations on the contents of a file must be done using this variable (not the file name itself!). Once all operations are finished, the file must be closed. Importantly, your computer's operating system limits the number of files which can be open at once (type the command `ulimit -n` to the command line to see how many), so it is very important to always close files when you are finished.

There are two basic ways to open and close files, as follows. Note that these two chunks of code are equivalent in their output.
```python
# Open and close with open() and close()
file_handle = open("important_file.txt", "r") # Open as read-only
contents = file_handle.read() # Read contents of file
file_handle.close() # Close file when finished (important!!)

# Open with the *with* control flow command. The file automatically closes outside the scope of the with.
with open("important_file.txt", "r") as file_handle:
  contents = file_handle.read()
```

The `open()` function takes two arguments: i) the *name* (including full or relative path!) of the file to open, and ii) the *mode* in which the file should be read. Modes include read-only (`"r"`), write-only (`"w"`), append (`"a"`), or read+write (`"rw"` for PCs and `"r+"` for Mac/Linux). Writing and appending are similar to the bash operators `>` and `>>`; write will overwrite the file, but append will add to the bottom of an existing file. Note that if you open a file for writing (or appending), that file doesn't need to exist already, and a new file will be created with the specified name. Alternatively, if you attempt to open a file that does not exist for reading, you'll receive an error message.

```python
# Open a file for writing, and write to it
file_handle = open("file_to_fill.txt", "w") # Open file for writing
file_handle.write("I'm writing a sentence to this file!")
file_handle.close()

# Open a file for appending, and append text to it
file_handle = open("file_to_fill.txt", "a")
file_handle.write("I'm writing another line to this existing file!")
file_handle.close()
```


## Looping over file contents

Typically when interacting with files, we want to examine and perform computations with the file content. To do this, we need to read in the file contents after opening the file. There are two basic file methods, `.readlines()` and `.read()`, that exist in native Python for this purpose. 

```python
# Read file line-by-line with .readlines() [Note that this is slow for *MASSIVE* files!]
with open("file.txt", "r") as f:
    lines = f.readlines()
    for line in lines:
        print line
    # To loop another time, you must direct the readlines iterator to start from the top! The same would go for .read().
    f.seek(0)
    for line in lines:
        print line
        
# Read file as a whole, and then split into lines
with open("file.txt", "r") as f:
    contents = f.read()
    contents_list = contents.split("\n") # Split file contents on newline character
    for line in contents_list:
        print line
```

## The csv module

The **csv module** is a useful python module for parsing comma-separated files, tab-delimited files, or generally any files with delimited *fields*. The csv module provides functions for parsing a file which has already been opened using `open()`.
```python
import csv

# Read a standard csv file.
f = open("file.csv", "r")
reader = csv.reader(f)
for row in reader:
  print row
f.close()

# Read a file with a different delimiter (e.g. tab)
f = open("file.txt", "r")
reader = csv.reader(f, delimiter = '\t') # Specify the delimiter if it is not a comma!!
for row in reader:
  print row
f.close()
```

## File parsing

Several examples of file parsing are available in [project4_files](project4_files/).
The file [parse_delimited.py](project4_files/parse_delimited.py) contains examples for parsing and extracting information from csv and tab-delimited files. It parses the data files [AbilAhuG\_uniprot\_blastx.txt](project4_files/AbilAhuG_uniprot_blastx.txt) and [AbilAhuG\_uniprot\_blastx.csv](project4_files/AbilAhuG_uniprot_blastx.csv) (note that these are the same file, except one is tab-delimited and one is comma-delimited).
The file [parse_custom.py](project4_files/parse_custom.py) contains an example parser for a file with non-canonical formatting, [hyphy_output_long.txt](project4_files/hyphy_output_long.txt').

An integral aspect of this file parsing is the **re module**. This python module provides functions for using **regular expressions**. Regular expressions are essentially flexible pattern-matching symbols (see the lesson [cheatsheet](../../Cheatsheets/Cheatsheet_Python4.md) for some commonly-used regex's. The re module, and indeed regular expressions in general, are extremely powerful and endlessly useful. Note that the re module has many, many more available functions associated with it (see the [re module documentation](https://docs.python.org/2/library/re.html)) beyond what is discussed here, but the examples shown below are (in your teacher's opinion..) among the most useful.


```python
import re
# All functions shown below are executed with the general format, re.<function>(regex_pattern, string_to_look_for_regex_in)

# re.search function

# re.findall function

# re.split function

# re.sub function

```


















