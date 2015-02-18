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


##Useful Commands

*You should set TextWrangler as the default application for `.py` files, and you should also set TextWrangler to be your default text editor.* 
If you know the command for a PC for these short cuts, please let us know and we'll add it to the list!

Command        | OS: program         |  Description          
---------------|-----------------|------------
`open <file_name>` | mac: terminal | opens the given file in the default application for that file's extension 
`open -t <file_name>` | mac: terminal | force your default text editor to open the file
`⌘ + /` | mac: TextWrangler | Highlight some text and enter this command to comment out the text; does not work for .txt files
`⌘ + <arrow_key>` | mac: TextWrangler | Navigate to the top/bottom of document with up/down arrows. <br> Navigate to beginning/end of line with left/right arrows
`option + click` | mac: terminal  | Jump to where you click on the command line
`option + <arrow key>` | mac: terminal  | Jumps complete words
`option + click&drag` | mac: TextWrangler  | Highlight a column of text
`alt + click&drag` | PC: Notepad++ | Highlight a column of text






