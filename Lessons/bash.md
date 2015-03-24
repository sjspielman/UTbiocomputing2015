# Merging python and unix, part II

## Bash is great for getting a quick look at your data and for simple regex replacements

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
# use Regex to keep only the characters ([a-zA-Z' ']*) that come before a number, ie not a space or letter, in each line
# - must escape all () using `\`, ie: \([regex]\)
# - must put whitespace and replacement (\1\2) in quotes, or else it is interpreted as a separate command
sed -E s/\([a-zA-Z' ']*\)\(.*\)/'\1'/ example 
```

* To insert tabs (`\t`) you'll have to hit `ctrl + v`, then `Tab` while in the terminal environment
* For newline characters (`\n`), you have to code it directly into the line, for example:

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
* awk, unlike sed, _does_ understand `\n` and `\t`

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
# print columns 1 and 2, adding in a tab between, a new field "newline", and a new line at the end
cat WEEK_06_python5_HW.csv | awk -F, '{print $1"\t"$2 "\tnewline\n"}'
```

```bash
# more complicated, a for loop that prints each item 
cat WEEK_06_python5_HW.csv | awk -F, '{for (i=1;i<=4;i++) {print $i}}'
# for each line, from items 1 until 4: `for (i=1;i<=4;i++)`
# print each item: `{print $i}`
```

```bash
# much more complicated, a for loop that prints a 10-nucleotide sequence that overlaps \
 by 1 nucleotide along the entire sequence
# you should recall this first part that grabs the sequence data from an illumina output file
bunzip2 -c AzetekiG_SRR957179.fastq.bz2 | grep ^[ACTG] | grep -v [FHJB] | head -10000 | sort | uniq | \
awk '{for (i=1;i<=length($1)-10;i++) {print substr($1,i,10)}}'
# for 1 to the length of the line until 10 before the end: `for (i=1;i<=length($1)-10;i++)`
# print from the entire line (`$1`) a substring from i to 10: `{print substr($1,i,10)}`
```
