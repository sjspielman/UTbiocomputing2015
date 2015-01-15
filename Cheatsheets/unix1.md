# Navigating your computer: basic UNIX.
<br><br>

## Directories and Paths

**1. Directory**
  - A *folder* on your computer which contains files. UNIX filesystems are organized as hierarchical directories. 
  - Forward slashes divide levels in the nested hierarchy of directories, e.g. `/top_level_directory/second_level_directory`
  - The directory at the top of this hierarchy is called the **root** directory and is denoted simply as `/`. 

**2. Path**
  - The *address* to a directory or file on your computer. There are, generally, two types of paths:
    1. **Absolute/full path** represents the path of a given directory/file beginning at the root directory.
    2. **Relative path** represents the path of a given directory/file relative to the working/current directory.

    For example, say you have a file "my\_favorite\_file.txt" located in the directory `/Users/myname/Desktop/my_directory`.
    - The **full path** to this file  is `/Users/myname/Desktop/my_directory/my_favorite_file.txt`.  
    - The **relative path** to this file depends on where you are on the computer. 
     - If you are calling this file from Desktop, the relative path would be `my_directory/my_favorite_file.txt`
     - If you are in `/Users/myname/`, the relative path becomes `Desktop/my_directory/my_favorite_file.txt`.
    
   Remember - Whenever you call the full path, you can reach the file from anywhere on your computer. Relative paths will change based on your current location.

**3. Roaming around**

Command | Translation | Examples
--------|-------------|---------
`cd` | **c**hange **d**irectory | `cd /absolute/path/of/the/directory/` <br> Go to the home directory by typing simply  `cd` or `cd ~` <br> Go up (back) a directory by typing `cd ..`
`pwd` | **p**rint **w**orking **d**irectory | `pwd`
`mkdir` | **m**ake **dir**ectory | `mkdir newDirectory` creates newDirectory in your current directory <br> Make a directory one level up with `mkdir ../newDirectory`
`cp` | **c**o**p**y | `cp file.txt newfile.txt` (and file.txt will still exist!)
`mv` | **m**o**v**e | `mv file.txt newfile.txt` (but file.txt will *no longer* exist!)




- `~` stands for your **home** directory.
  - `..` means "go up one directory."
