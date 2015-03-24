#CHEATSHEET LESSON 8: Merging Python and Bash


## Python `sys` module

Module | Command  |  Description | Example
-------|----------|--------------|---------
`sys` | `sys.argv` | Contains a list of command-line arguments (all strings!) |`save_me = sys.argv[1]`
`sys` | `sys.path` | Contains a list of directories in the python path |`sys.path.append("another_path_i_want_to_access")`

<br><br>

## Python `os` and `shutil` modules
Module | Command  |  Description | Unix equivalent | Example
-------|----------|--------------|-----------------|--------
`os` | `os.system` | Run a UNIX command from python | anything you would type into a command line | `os.system("pwd")`
`os` | `os.listdir`| List all items in a given directory | `ls` | `os.system("/directory/of/interest/")`
`os` | `os.remove` | Remove a file | `rm` | `os.remove("i_hate_this_file.txt")`
`os` | `os.rmdir` | Remove a directory | `rm -r`| `os.rmdir("/i/hate/this/directory/")`
`os` | `os.mkdir`  | Create a new directory | `mkdir` |`os.mkdir("/path/to/brand/new/directory/")`
`os` | `os.mkdirs`  | Create many new directories | `mkdir`|`os.mkdir("/path/to/a/brand/new/directory/", "/path/to/another/brand/new/directory/")`
`os` | `os.chdir`  | Change directory where python is running | `cd` | `os.chdir("/another/directory/where/i/want/to/be/")`
`shutil` | `shutil.copy` | Copy a file | `cp` | `shutil.copy("old_file.txt", "new_file.txt")`
`shutil` | `shutil.move` | Move a file | `mv` | `shutil.move("old_file.txt", "new_file.txt")`



# BASH section

