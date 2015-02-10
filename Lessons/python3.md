# Python I


## Functions
Functions are an integral part of programming. They are self-contained pieces of code which provide instructions for a given task. Functions are *called* with certain input, *execute* the code specified, and then *return* specified value(s) from the calculations performed. The presence of a function in a program does not guarantee that it will run - functions must be explicitly called to run. 

Functions should be an integral component of your scripts/program, allowing for modular, repeatable, and readable programs. *If/when you perform a certain calculation multiple times throughout your program, or you find yourself copy/pasting code into many locations in a program, it belongs in a function!* 


#### Syntax

Functions are defined using a `def` statement and their contents are indicated with whitespace (as with if, for, while):
```python
# Basic anatomy of a function which takes two arguments and returns a single value, x
def name_of_function(argument1, argument2):
    # code
    # more code
    return x
```
The input arguments (argument1, argument2) *do not* correspond to actual variables used in your code, and they only exist within the environment of the function. They are simply part of the generic function instructions and can be named however you like.

For example, we can create a simple function that adds two numbers together. First, think about the actual task of addition itself: you are given two numbers, you add them together, and then you report the result.

```python
# Simple addition function. 
def addition(x, y):
    sum = x + y
    return sum
    
# Another way to write the same function
def addition2(x, y):
    return x + y

>>> # We can now use the function with arbitrary input arguments. 
>>> addition(5, 6)
 11

>>> a = 5
>>> b = -6
>>> addition(a, b)
 -1

>>> # The input variables remained unchanged
>>> print a
 5

>>> print b
 6

>>> # Assign the output of a function to a variable
>>> c = addition(a,b)
>>> print c
 -1
```




#### Scope

Python more or less goes in order

#### Positional, keyword, and optional arguments

#### Modules

#### Docstrings













