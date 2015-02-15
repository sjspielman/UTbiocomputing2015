# Python III


## Functions
Functions are an integral part of programming. They are self-contained pieces of code which provide instructions for a given task. Functions are *called* with certain input, *execute* the code specified, and then *return* specified value(s) from the calculations performed. The presence of a function in a program does not guarantee that it will run - functions must be explicitly called to run. 

Functions should be an integral component of your scripts/program, allowing for modular, repeatable, and readable programs. Using functions dramatically lowers the chance of bugs occuring in your program and cleans up a lot of clutter.  **If/when you perform a certain calculation multiple times throughout your program, or you find yourself copy/pasting code into many locations in a program, that code belongs in a function!** 


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
>>> addition(85, 6)
 91

>>> a = 5
>>> b = -6
>>> addition(a, b)
 -1

>>> # The input variables remained unchanged
>>> print a, b
 5, -6

>>> # Variables defined inside the function only exist inside the function!!
>>> print x
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
NameError: name 'x' is not defined

>>> # Assign the output of a function to a variable
>>> c = addition(a,b)
>>> print c
 -1
 
>>> # Be sure to provide all arguments, or you'll get an error message!
addition(6) # only 1 argument provided
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: addition() takes exactly 2 arguments (1 given)
```

Functions can also return multiple values, stored as a *tuple*. This means that the values returned are separated by commas and constitute an *immutable container*. This is good news - as tuples can't be changed once created, the returned values are safe!
```python
# Function to divide two numbers and return the result and remainder
def divide_remain(x, y):
    div = x / y
    rem = x % y
    return div, rem

>>> # Run the function!
>>> divide_remain(77, 12)
(6, 5) # the order div, rem is preserved

>>> # Save the values returned into a single variable, and then use indexing to see the values
>>> a = divide_remain(77, 12)
>>> print a
(6, 5)
>>> type(a)
<type 'tuple'>
>>> print a[0]
6
>>> print a[1]
5

>>> # Save each returned value to its own variable
>>> a, b = divide_remain(77, 12)
>>> print a
6
>>> print type(a)
<type 'int'>

>>> print b
5
>>> print type(b)
<type 'int'>
```

#### Positional, keyword, and optional arguments

In the previous examples, all functions took a pre-specified number of arguments in a particular order. These are called *positional* arguments - the order in which you provided the arguments actually matters for how those arguments are used inside the functions. Other types of arguments allowing for more flexibility, however, are possible!

###### Keyword argument
If used in combo with positional arguments, all keyword arguments come at the *end*.
Example, 
```python

# Division function with two positional arguments
def divide_pos(x, y):
    return x / y 

>>> divide_pos(2,4)
2
>>> divide_pos(2,5) #womp womp! 
2

# Division with optional keyword argument, as_float. This function takes two positional arguments and one optional keyword argument.
def divide_key(x, y, as_float = True):
    if as_float:
        return float(x) / float(y)
    else:
        return x / y 

>>> divide_key(2,5)
2.5
>>> # When specifying the non-default value, you must provide the keyword
>>> divide_key(2,5, as_float = False)
2
```

A function can have as many keyword arguments as you want. When you call the function, you can specify these arguments in any order as long as they all come *after* all positional arguments. If you intend to use the default behavior of such an argument, no need to supply it! That's why the defaults are there.
```python
def divide(x, y, as_float = True, digits = 3, print_sentence = False, return_remainder = False):
    if as_float:
        div = round(float(x) / float(y), digits)
    else:
        div = round( x / y, digits)
    
    if print_sentence:
        if as_float:
            print "The result of %d / %d is %f." %(x, y, div)
        else:
            print "The result of %d / %d is %d." %(x, y, div)
    
    if return_remainder:
        return div, x%y
    else:
        return div
        
>>> a = divide(6, 767, print_sentence = True, digits = 10, return_remainder = True)
The result of 6 / 767 is 0.007823.
>>> print a # think about what type a will be and why! Also, think about why there are a different number of digits in the print statement and the final dividend returned.
(0.0078226858, 6)
```

#### Modules

Let's say you have several scripts which all perform similar tasks, and therefore require the same functions. One way to do this is simply to include your functions in every script. An alternative (and dare-I-say, better?) strategy is to create a stand-alone python script which contains only functions - this is a module! You can then import this module into the scripts which use these functions. This strategy will help ensure that you don't accidentally introduce bugs from copy/pasting the function, and more importantly allows you to change the function *only one time* as opposed to individually in each script where it's used (no matter how diligent you are, the latter strategy **will** introduce bugs!). 

For examples, see the scripts `my_math_module.py` (contains module functions) and `my_script.py` (actually runs stuff).



#### Docstrings
It is always (read: **always**) a good (read: **absolutely the most important**) idea to incorporate docstrings into your functions. Docstrings are essentially comments placed inside three quotation-mark bounds (""" words """) which explain the purpose, functionality, input arguments, and return values for your function. Docstrings are great because they explain to you and others looking at your code what exactly the function accomplishes, without the reader having to fully read and internalize all the code. Also, as a bonus, if you ever want to document your python code, there are awesome tools out there (like Sphinx) which will automatically create beautiful documentation from your python code using these docstrings. 
The docstrings are also shown whenever call `help()` on a given function.

Let's rewrite the `divide_remain()` function with docstrings included.
```python
def divide_remain(x, y):
    """ Function to compute the dividend and remainder between two numbers.
        Arguments: Two positional arguments, x and y. x is the numerator and y is the denominator.
        Returns: Dividend and remainder from x / y calculation.
    """
    div = x / y
    rem = x % y
    return div, rem

>>> help(divide_remain)
divide_remain(x, y)
    Function to compute the dividend and remainder between two numbers.
    Arguments: Two positional arguments, x and y. x is the numerator and y is the denominator.
    Returns: Dividend and remainder from x / y calculation.
```











