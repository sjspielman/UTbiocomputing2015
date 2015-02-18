#CHEATSHEET PYTHON 3: Functions

##Functions 

Command  |  Description 
----------|-------------
`continue` | Proceed immediately to the next iteration of a loop (for/while)
`break`    | Exit immediately from a loop (for/while)

##Defining functions
```python
def divide_key(x, y, as_float = True, show_remainder = False):
  """ Divide two numbers, x and y. Takes optional arguments "as_float" (default True) and "show_remainder" (default False). """
    if show_remainder:
      print "The remainder is %d" %(x%y)
    if as_float:
        return float(x) / float(y)
    else:
        return x / y 
```


##Useful TextWrangler Commands
Note that you can open a files directly from the command line using the command "`open <file_name>`". The "`open`" command opens the given file in the default application for that file's extension. 
Alternatively, type "`open -t <file_name>`" to force your default text editor to open the file.
*You should set TextWrangler as the default application for `.py` files, and you should also set TextWrangler to be your default text editor.*

Command  |  Description 
----------|-------------
`⌘ + /` | Highlight some text and enter this command to comment out the text
`⌘ + <arrow_key>` | Navigate to the top/bottom of document with up/down arrows. <br> Navigate to beginning/end of line with left/right arrows.
`⌘ + click` | Highlight a column of text 






