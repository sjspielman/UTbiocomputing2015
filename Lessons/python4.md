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

The `open()` function takes two arguments: i) the name of the file to open, and ii) the *mode* in which the file should be read. Modes include read-only (`"r"`), write-only (`"w"`), append (`"a"`), or read+write (`"rw"` for PCs and `"r+"` for Mac/Linux). Writing and appending are similar to the bash operators `>` and `>>`; write will overwrite the file, but append will add to the bottom of an existing file.

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








