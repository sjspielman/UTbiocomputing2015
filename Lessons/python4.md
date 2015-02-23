# Python IV: File Input/Output

## Opening and closing files

```python
# Open files and save as a file handle. All future operations are on this handle (not the file name itself!!)
file_handle = open("important_file.txt", "r") # Open as read-only
contents = file_handle.read() # Read contents of file
file_handle.close() # Close file when finished (important!!)
```









