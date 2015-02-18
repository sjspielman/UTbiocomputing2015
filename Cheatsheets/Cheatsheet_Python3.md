#CHEATSHEET PYTHON 3: Functions

##Functions 

Command  |  Description 
----------|-------------
`continue` | Proceed immediately to the next iteration of a loop (for/while)
`break`    | Exit immediately from a loop (for/while)

##Defining functions
Use the following syntax:
```python
def function_name(input_variable1, input_variable2, optional_argument = True/False):
  """Doc string to remind you of function use and input. The more detailed the better!"""
    for loops, if else, other statements
    more code
    return new_variable #this is what your function outputs
```

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

Command        | computer: program         |  Description          
---------------|---------------------------|------------
`open <file_name>` | mac: terminal | opens the given file in the default application for that file's extension 
`open -t <file_name>` | mac: terminal | force your default text editor to open the file
`option + click` | mac: terminal  | Jump to where you click on the command line
`option + <arrow key>` | mac: terminal  | Jumps complete words
`⌘ + /` | mac: TextWrangler | Highlight some lines and enter this command to uncomment/comment the text; does not work for .txt files
`⌘ + <arrow_key>` | mac: TextWrangler | Navigate to the top/bottom of document with up/down arrows. <br> Navigate to beginning/end of line with left/right arrows
`option + click&drag` | mac: TextWrangler  | Highlight a column of text
`alt + click&drag` | PC: Notepad++ | Highlight a column of text






