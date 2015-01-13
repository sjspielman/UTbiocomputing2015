# Navigating your computer: basic UNIX.



## Directories and Paths

**1. Directory**
  - A *folder* on your computer which contains files. UNIX filesystems are organized as hierarchical directories.
  - The directory at the top of this hierarchy is called the **root** directory and is denoted "/". 

**2. Path**
  - The *address* to a directory or file on your computer. There are, generally, two types of paths:
    1. **Absolute/full path** represents the path of a given folder/file *starting from root*.
    2. **Relative path** represents the path of a given folder/file *relative to the working/current directory*.
    
    
    For example, say you have a file "my\_favorite\_file.txt" located in the directory `/Users/myname/Desktop/my_directory.`
    It's *full path* is `/Users/myname/Desktop/my_directory/my_favorite_file.txt`.
    It's *relative path* depends on where you are on the computer. If you are calling this file from Desktop, then it's relative path is `my_directory/my_favorite_file.txt`, which tells you how to get to the file from the current directory.

**3. Path shortcuts**
  - `~` stands for your **home** directory.
  - `..` means "go up one directory."