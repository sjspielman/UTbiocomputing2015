#CHEATSHEET PYTHON 3: File Input/Output

##Functions and Methods

Command  |  Description | Example
----------|-------------|----------
`open` | Open a file as a python object |`{python} my_file = open("file.txt", "r")`
`close`| Close a file | `{python} my_file.close()
`seek` | Place the file iterator at a certain line. *Needed* to loop over a file contents multiple times! | `{python} my_file.seek(0)`

##csv Module functions and methods
Command  |  Description | Example
----------|-------------|----------
`csv.reader` | Setup a csv reader on an **opened** file handle | `{python} my_file = open("file.csv", "r")`<br>`{python} reader = csv.reader(my_file)`<br> Note that an argument `delimiter=...` may be provided for other separators, like `\t` for tabs.

##File modes
r
w
a
r+ - mac/linux rw
rw - pc rw
rU - read unicode characters



##Regular Expressions

Symbol | Meaning
-------|---------
`\s`   | space character
`\t`   | tab character
`\n`   | newline character (Mac and Linux! PCs may or may recognize this)
`\r`   | newline character (PCs! Mac and Linux may or may recnogize this)
`.`    | wildcard
`\d`   | digit (numbers only!)
`\w`   | letter or number (case insensitive)
`+`    | Symbol to append after a regular expression indicating "one or more of these" <br> E.g., `\d+` means match 1 or more numbers
`*`    | Symbol to append after a regular expression indicating "zero or more of these" <br> E.g., `\d*` means match 0 or more numbers
