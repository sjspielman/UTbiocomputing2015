# Navigating your computer: basic UNIX.
<br><br>

**Directory**
  - A *folder* on your computer which contains files. UNIX filesystems are organized as hierarchical directories. 
  - Forward slashes divide levels in the nested hierarchy of directories, e.g. `/top_level_directory/second_level_directory`
  - The directory at the top of this hierarchy is called the **root** directory and is denoted simply as `/`. 

**Path**
  - The *address* to a directory or file on your computer. There are, generally, two types of paths:
    1. **Absolute/full path** represents the path of a given directory/file beginning at the root directory.
    2. **Relative path** represents the path of a given directory/file relative to the working/current directory.

    For example, say you have a file "my\_favorite\_file.txt" located in the directory `/Users/myname/Desktop/my_directory`.
    - The **full path** to this file  is `/Users/myname/Desktop/my_directory/my_favorite_file.txt`.  
    - The **relative path** to this file depends on where you are on the computer. 
     - If you are calling this file from Desktop, the relative path would be `my_directory/my_favorite_file.txt`
     - If you are in `/Users/myname/`, the relative path becomes `Desktop/my_directory/my_favorite_file.txt`.
    
   Remember - Whenever you call the full path, you can reach the file from anywhere on your computer. Relative paths will change based on your current location.

**Roaming around**

Command | Translation | Examples
--------|-------------|---------
`cd` | **c**hange **d**irectory | `cd /absolute/path/of/the/directory/` <br> Go to the home directory by typing simply  `cd` or `cd ~` <br> Go up (back) a directory by typing `cd ..`
`pwd` | **p**rint **w**orking **d**irectory | `pwd`
`mkdir` | **m**ake **dir**ectory | `mkdir newDirectory` creates newDirectory in your current directory <br> Make a directory one level up with `mkdir ../newDirectory`
`cp` | **c**o**p**y | `cp file.txt newfile.txt` (and file.txt will still exist!)
`mv` | **m**o**v**e | `mv file.txt newfile.txt` (but file.txt will *no longer* exist!)
`rm` | **r**e**m**ove | `rm file.txt` removes file.txt <br> `rm -r directoryname/` removes the directory and all files within
`ls` | **l**i**s**t | `ls *.txt` lists all .txt files in current directory <br> `ls -a` lists all files including hidden ones in the current directory <br> `ls -l` lists all files in current directory including file sizes and timestamps <br> `ls -lh` does the same but changes file size format to be **h**uman-readable <br> `ls ../` lists files in the directory above the current one
`man` | **man**ual | `man ls` opens the manual for command `ls` (use `q` to escape page)
`grep` | **g**lobal **r**egular <br> **e**xpression **p**arser |  `grep ">" seqs.fasta` pulls out all sequence names in a fasta file <br> `grep -c ">" seqs.fasta` counts the number of those sequences <br> 
`cat` | con**cat**enate | `cat seqs.fasta` prints the contents of seqs.fasta to the screen (ie stdout)
`head` | **head** | `head seqs.fasta` prints the first 10 lines of the file <br> `head -n 3 seqs.fasta` prints first 3 lines
`tail` | **tail** | `tail seqs.fasta` prints the last 10 lines of the file <br> `tail -n 3 seqs.fasta` prints last 3 lines
`wc` | **w**ord **c**ount | `wc filename.txt` shows the number of new lines, number of words, and number of characters <br> `wc -l filename.txt` shows only the number of new lines <br> `wc -c filename.txt` shows only the number of characters

**Getting Comfortable**
  - tab completion
  - using `*` as a wildcard
  - up arrow can be used to call the last command
  - Ctrl + C kills current process
  - Ctrl + L (or `clear`) clears screen
  - relative paths:
    - `.` = here
    - `..` = one level up 
  - absolute paths:
    - `~` = home
    - `/` = root
    - `/usr/bin` = program file location
  - `>` redirects stdout to a file, *overwriting* file if it already exists
  - `>>` redirects stdout to a file, *appending* to the end of file if it already exists
  - `|` redirects stdout to become stdin for next command







