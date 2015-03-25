# Merging Python and UNIX

## Automating file renaming

First, a quick trick for renaming lots of files! Let's imagine you have a bunch of files named, 
+ env1.csv
+ env2.csv
+ env3.csv
+ env4.csv
+ env5.csv
+ ... 

But, perhaps you're not satisfied with them being called "env", and you'd prefer them to be called environment1.csv, etc. Now, let's say there are hundreds of these files. Renaming these files individually would be hugely time-consuming, so.. regular expressions to the rescue! All you'll need is a text editor which can search/replace using regular expressions (good options are TextWrangler or Sublime Text). Follow these steps (or similar, depending on your specific data files):

**1.** From the command line, type `ls env*` to list all files named as env1.csv, env2.csv ...

**2.** Copy/paste this list of files into a text editor which can accomodate regular expressions.

**3.** You'll want each file name to have its own line, with no empty lines in between. You can clean up whitespace with the following two regular expression search/replace pairs (in order):

Search | Replace | What it does
-------|---------|-------------
`\s*`  | `\n`    | Replace all whitespace with a newline character
`\n+`  | `\n`    | Replace all "double newlines" with a single newline (no empty space between lines)

**4.** Finally, perform a search/replace with the following regular expressions:
<br>Search: `^env(\d+\.csv)`
<br>Replace: `mv env\1 environment\1`

**5.** Copy and paste the new search/replaced text file contents back to the command line, and viola! All files will have been renamed to environment1.csv, etc.


<br><br>
## The python `sys` module

The `sys` module interacts with the python interpreter. This module primarily comes with certain **variables** (rather than functions) which are particularly useful, two of which are described below.

### Passing command-line arguments `sys.argv`
Often, you'll want to pass input arguments to a script. All input arguments are stored in the variable `sys.argv` (note that you must import the `sys` module!). For example, you might have a script called `calc_dna.py` which performs certain calculations on a sequence data file, but each time you run the script, you might want to process a different file. One option to pass the file name as an input argument: `python calc_dna.py inputfile.fasta`, where `inputfile.fasta` is the single command-line argument. If you've loaded the `sys` module, all command line arguments will be stored in `sys.argv`:

Assume the following code is included in the script calc_dna.py. 
```python
import sys
print sys.argv
```
From the command line you run `python calc_dna.py inputfile.fasta` and get back:
`['calc_dna.py', 'inputfile.fasta']`


Notice that the first entry in `sys.argv` is the name of the script. After this come all command line arguments! In addition, all `sys.argv` entries will be **strings**. So remember that if you want to use an input argument as a number, you must convert it to a float or integer.

Generally, you should save the input arguments to a new variable inside the script:
```python
import sys
infile = sys.argv[1]
```

Error checking is often very useful here; you might have a script which **requires** two command line arguments. For instance, let's say you have a script (which does something..) which takes a file name and a number as its two arguments. To ensure that this always happens, help yourself out with assertion statements:

```python
import sys
# sys.argv must be of length 3 (script name, inputfile, number)
# print a usage statement if arguments not provided or provided incorrectly
assert( len(sys.argv) == 3 ), "Usage: python my_script.py <inputfile> <number>"
infile = sys.argv[1]
number = int( sys.argv[2] ) # remember to convert from string to int, as needed!
```


### Editting the path with `sys.path`

You can view everything in python's path by printing the contents of the list `sys.path`. This variable will you tell which directories on your computer that the current python interpretter is able to access. To edit this path in place, for instance by adding a directory to the path, simply use `.append()`:
```python
import sys
sys.path.append("/path/to/directory/that/python/should/know/about/")
```

<br><br>
## The python `os` and `shutil` modules

The `os` and `shutil` modules are useful for interacting with your computer's operating system (typically UNIX). With these modules, you can run commands from your python script which are analogous to UNIX commands like `cd` and `pwd`. 
<br>Some examples:

Module | Command  |  Description | Unix equivalent | Example
-------|----------|--------------|-----------------|--------
`os` | `os.listdir`| List all items in a given directory | `ls` | `os.system("/directory/of/interest/")`
`os` | `os.remove` | Remove a file | `rm` | `os.remove("i_hate_this_file.txt")`
`os` | `os.rmdir` | Remove a directory | `rm -r`| `os.rmdir("/i/hate/this/directory/")`
`os` | `os.mkdir`  | Create a new directory | `mkdir` |`os.mkdir("/path/to/brand/new/directory/")`
`os` | `os.mkdirs`  | Create many new directories | `mkdir`|`os.mkdir("/path/to/a/brand/new/directory/", "/path/to/another/brand/new/directory/")`
`os` | `os.chdir`  | Change directory where python is running | `cd` | `os.chdir("/another/directory/where/i/want/to/be/")`
`shutil` | `shutil.copy` | Copy a file | `cp` | `shutil.copy("old_file.txt", "new_file.txt")`
`shutil` | `shutil.move` | Move a file | `mv` | `shutil.move("old_file.txt", "new_file.txt")`


### Running external commands with `os`

You will often want to use Python scripting to automate analyses which use external programs or softwares. You can actually call these programs directly from your python script using the function `os.system()`. This function takes a single argument: the command you want to run (as a string). Anything that you could type into the command line can be given to `os.system`!

```python
# Create a multiple sequence alignment in MAFFT from python
import os
# Define input, output files
infile = "unaligned.fasta"
outfile = "aligned.fasta"
command = "mafft " + infile + " > " + outfile
os.system(command)
```

We can also incorporate `sys` to provide the input/output file names as arguments!
```python
# Create a multiple sequence alignment in MAFFT from python
import os
import sys

# Check and save input arguments
assert( len(sys.argv) == 3 ), "Usage: python align.py <inputfile> <outputfile>"
infile = sys.argv[1]
outfile = sys.argv[2]

# Run the alignment
command = "mafft " + infile + " > " + outfile
os.system(command)
```

Finally, you can check that the command has run properly by saving the output of `os.system()` (basically, save it to a variable). In UNIX, a returned value of 1 means an error occurred, but a returned value of 0 means everything went fine. Therefore, we want to make sure that `os.system()` returns a value of 0, by editting the last few lines:

```python
command = "mafft " + infile + " > " + outfile
aligned_properly = os.system(command)
assert(aligned_properly == 0), "MAFFT didn't work!"
```

# Part 2: UNIX/Bash

<br><br>
### `sed` is useful for quick and recursive replacements using Regex

* General pseudocode: `sed [-E] command/regex/replacement/optionalflag filename > newfile`
* My favorite pseudocode: `sed -E s/OLD/NEW/`
* Mac users must include `-E` to access regular expressions
* `sed` does not understand `\t` and `\n`, see below

```bash
# replace first instance of XX with YY for each line
sed s/XX/YY/ filename > newfile
```

```bash
# replace all instances of XX with YY 
# - `g` flag means 'global' and searches for all instances of the pattern
sed s/XX/YY/g filename > newfile
```

```bash
# replace all instances of XX with YY and of AA with ZZ
# - `-e` flag lets you execute multiple sed commands at once
sed -e s/XX/YY/g -e s/AA/ZZ/g filename > newfile
```

```bash
# keep only letter and space characters ([a-zA-Z' ']*) that come before a different type of character in each line
# - must escape all () using `\`, ie: \([regex]\)
# - must put whitespace and replacement ('\1\2') in quotes, or else it is interpreted as a separate command
sed -E s/\([a-zA-Z' ']*\)\(.*\)/'\1'/ example 
```

* To insert tabs (`\t`) you'll have to hit `ctrl + v`, then `Tab` while in the terminal environment
* For newline characters (`\n`), you have to code it directly into the line with `\ + enter`, for example:

```bash
sed -E s/\([0-9]\)/'\1\
'/ example 

```

### `sort` sorts the lines in a file by numbers then lowercase letters then uppercase letters
Command | Meaning | Example
----------|--------|---------
-b | ignore leading blanks | `sort -b filename > filename.sorted`
-r | reverse | `sort -r filename`
-k POS1 | sort by field/character indicated by POS1 | sort by field 2: `sort -k 2 filename` <br> sort by second character in field 2: `sort -k 2.2 filename`
-k POS1,POS2 | sort based on the characters from POS1 to POS2 | sort by characters in fields 2 and 3: `sort -k 2,3 filename` <br> sort starting with second character in field 2 up to and including field 3: `sort -k 2.2,3 filename`


### `uniq` counts the number of unique ADJACENT lines in a file; use `sort` beforehand to ready the file for `uniq`
Command | Meaning | Example
-------|--------|---------
-c | prefixes lines with the number of times they occur | `uniq -c filename`
-d | prints only repeated lines | `uniq -d filename`
-u | prints only unique lines | `uniq -u filename > filename.unique`
-f N | skips N number of lines | `uniq -f 30 filename`

A combination of grep, sort, and uniq is great for looking at raw sequence reads

```bash
# unzip file (`bunzip2`) and print contents to stdout rather than to a new file (`-c`)
# select any lines that start with A, C, T, or G (`grep ^[ACTG]`)
# keep lines that do not contain F, H, J, or B (`grep -v [FHJB]`)
# print the first 10000 lines, sort them alphabetically, show unique lines with counts (`head -10000 | sort | uniq -c`)
bunzip2 -c AzetekiG_SRR957179.fastq.bz2 | grep ^[ACTG] | grep -v [FHJB] | head -10000 | sort | uniq -c
```


### `awk` is useful for quick subsetting of tab- or csv-delimited datasets
* General pseudocode: `cat filename | awk '{ command }'`
* My favorite pseudocode: `cat filename | awk -Fdelimiter '{ print($linenumber,$otherlinenumber) }'`
* `awk`, unlike `sed`, _does_ understand `\n` and `\t`

```bash
# example using HW csv file
# print columns 2 and 3
cat WEEK_06_python5_HW.csv | awk -F, '{ print ($2,$3) }'
```

```bash
# print columns 2 and 3, add a `tab` between the items
cat WEEK_06_python5_HW.csv | awk -F, '{ print ($2"\t"$3) }'
```

```bash
# print columns 1 and 2, adding in a `tab` between, a new field "newline", and a new line at the end
cat WEEK_06_python5_HW.csv | awk -F, '{print ($1"\t"$2 "\tnewline\n")}'
```

```bash
# more complicated, a for loop that prints each item 
cat WEEK_06_python5_HW.csv | awk -F, '{for (i=1;i<=4;i++) {print ($i)}}'
# for each line, from items 1 until 4: `for (i=1;i<=4;i++)`
# print each item: `{print $i}`
```

```bash
# much more complicated, a for loop that prints a 10-nucleotide sequence that overlaps by 1 nucleotide along the entire sequence
# you should recall this first part that grabs the sequence data from an illumina output file
bunzip2 -c AzetekiG_SRR957179.fastq.bz2 | grep ^[ACTG] | grep -v [FHJB] | head -10000 | sort | uniq | \
awk '{for (i=1;i<=length($1)-10;i++) {print substr($1,i,10)}}'
# for 1 to the length of the line until 10 before the end: `for (i=1;i<=length($1)-10;i++)`
# print from the entire line (`$1`) a substring from i to 10: `{print substr($1,i,10)}`
```
