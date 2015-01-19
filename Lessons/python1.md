# Python I


## Object-oriented programming
Python is an object-oriented programming (OOP) language. Simply put, this means that data structures (e.g. variables) have certain properties, in particular **attributes** and **methods**. Attributes are usually descriptive features of the object, and methods are certain actions (functions) you can perform on the object itself. For example, consider a book. Books are physical objects with certain **attributes**, such as number of pages, number of words, dimensions (width, height, depth), cover art, etc. Similarly, there are several actions (**methods**) that we can do with books like reading, writing, throwing at people we hate, and sharing with people we love. A given realization of an object is called an **instance**. To continue with the book metaphor, "The Count of Monte Christo" and "Harry Potter" are each instances of book objects.

While these concepts might seem abstract, python syntax and behavior will make much more sense in light of the OOP paradigm upon which the language is based. In particular, attributes and methods are accessed in systematic ways, as follows:
```python
>>> # Assume mybook is a book instance

>>> # Call an attribute
>>> mybook.number_of_pages
>>> mybook.height
>>> mybook.width

>>> # Call a method
>>> mybook.read()
>>> mybook.write()
>>> mybook.throw_at_bad_people("Bob") # methods might take particular arguments, in this case who we're launching books at
```
As you see, attributes and methods are called after a . at the end of the instance name. Unlike attributes, however, end with parentheses. Sometimes, methods can take particular **arguments** which relate to their function, for instance when we threw the book at Bob.


## Operators

#### Mathematical operators
First, we have the standard operators for addition (+), subtraction (-), multiplication (*), and division (/). Additional mathematical operators include, 

1. Modulus operator `%`, which gives the remainder, e.g. `12 % 5` will result in 2.
2. Exponent operator, `**`, e.g. `5**2` will result in 25.

#### Logical operators
In addition to symbols used for basic calculations, python (like other languages!) has a series of symbols used to compare statements in a True/False context. Such statements are called **logical statements**, and the symbols we used to compare them are called **logical operators**. The most commonly used logical operators include, 

Symbol   | What it does | Example
---------|--------------|---------
  ==, is | Equals       |  `5 == 5` results in True <br> `9 is 9` results in True <br> `5 is 7` results in False
  !=, is not | Not equals       |  `5 != 5` results in False <br> `5 is not 7` results in True
 > | Greater than       |  `5 > 6` results in False <br> `11 > 6.23` results in True
 >= | Greater than or equal to |  `5 >= 4` results in True
< | Less than  |  `4 < 5` results in True
<= | Less than or equal to |  `4 <= 5` results in True

As the table above indicates, the symbol `==` is equivalent to using the word `is`, and similarly for the operation not equals, `!=` is equivalent to using the words `is not`. And above all, it is **very** important to remember that the equals logical comparison requires a **double equals sign** (`==`) - a single equals signs indicates variable assignment (see below).

We can also perform multiple comparisons at once, using the keywords `and` and/or `or`. 
For example, 
``` python
>>> # For "and," both logical statements must be true to result in True
>>> 5 == 5 and 6 == 6
True
>>> 5 > 7 and 8 < 10
False

>>> # For "or," at least one logical statement needs to be true to result in True
>>> 5 > 7 or 8 < 10
True

>>> # We can use the keyword "not" to negate a logical statement
>>> 3 < 11 and not 4 >= 7
True
```

## Variables
We assign values to a variable using the equals sign, `=`.
``` python
>>> a = 5
>>> # Check that the variable was correctly assigned using the "print" statement
>>> print a
5
```

All variables have a certain **type**. The variable type deterines what can be done with the variable itself. Standard python types include the following,

Variable Type   | Description | Casting
---------|--------------|---------
integer | whole number  | int()       
float   | decimal number | float()
string  | ordered, immutable character container | str()
list    | ordered, mutable container | list()
dictionary | unordered, mutable container | dict()
tuple | ordered, immutable container | tuple()

Remember, every time you create a variable, you are essentially creating an **instance** of that particular variable type. As a consequence, there are certain properties (attributes and methods!) associated with each type of variable.

### Integers and floats
Integers and floats are python's primary types for dealing with numbers. 

Integers are whole numbers only, but floats include decimal places. Whether a variable is an integer or float turns out to matter a lot - if you perform an operation with integers, the result will be an integer (even if the "real" answer is actually a float!)

```python

>>> a = 6
>>> type(a)
<type 'int'>
>>> # We can change the type of a variable by using casting
>>> a = float(a)
>>> type(a)
<type 'float'>

>>> # Caution! Variables can easily be reassigned!
>> a = 11
>>> print a
11

>>> # By adding a decimal point during assignment, we force the variable to be a float
>>> b = 6.
>>> type(b)
<type 'float'>

>>> # Be careful! If you perform operations with only integers, the result will always be an integer (rounding determines answer)
>>> x = 5
>>> y = 6
>>> x/y
0
>>> # One way to circumvent this issue is by casting the result as a float
>>> float(x / y)  
0.8333333333333334
>>> # Note that the above division will NOT change the casting of either x or y themselves
>>> x / y
0
>>> # Another solution would be to define either/both x or y as a float from the beginning
>>> x = 5.
>>> y = 6
>>> x / y
0.8333333333333334
```

### Strings
In python, a string is an **immutable** container of **characters**. By "immutable", we mean that it can't be changed - that is, once you create a string, you can't rewrite certain parts of it. It's an all-or-nothing thing. By characters, we basically mean "not numbers" - consequently, no mathematical operations can be done on strings. 
Strings are also **ordered** - python remembers the orders of values used in the list. This means we can **index** items in a list, again using brackets `[]`. **Importantly**, python indexes **starting at 0**, meaning the first item in a given string is the 0th entry.

```python
>>> # Assign strings using quotation marks
>>> name = "Stephanie"
>>> type(name)
<type 'str'>

>>> # Numbers can also be strings!
>>> age = "26"
>>> type(age)
<type 'str'>

>>> # Even though 26 is a number, we set up the variable age as a string, so no math can be performed with this variable
>>> age - 7
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: unsupported operand type(s) for -: 'str' and 'int'

>>> # But!! Since the value 26 is, in fact, a number we can recast age as an integer or float and then do maths with it!
>>> int(age) - 7
19
>>> # Again, age is still a string. We'd have to redefine the variable itself to make it an integer (or float) for good
>>> type(age)
<type 'str'>
>>> age = int(age)
>>> type(age)
<type 'int'>

>>> # We can't recast the variable name though, since letters just aren't numbers.
>>> name = float(name)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
ValueError: could not convert string to float: Stephanie

>>> # Strings are ordered, so we can index them
>>> name[5]
a
>>> # But strings are also immutable, so we can't edit strings in place.
>>> name[5] = "A"
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: 'str' object does not support item assignment

>>> # Index a range with x:y, where index x is included but index y is not.
>>> name[1:5]
'teph'

>>> # You can also index a range like x:y:z , where z is the step.
>>> name[1:8:2] # s TePhAnI e
'tpai'
```

Some useful string methods:

``` python
>>> my_string = "This is a really nice string for showing examples."

>>> # .upper() converts your string to upper-case
>>> my_string.upper()
'THIS IS A REALLY NICE STRING FOR SHOWING EXAMPLES.'

>>> # .lower() converts your string to lower-case
>>> my_string.lower()
'this is a really nice string for showing examples.'

>>> # .split() takes an *argument* to split the list on and creates a list containing each chunk
>>> my_string.split('i')
['Th', 's ', 's a really n', 'ce str', 'ng for show', 'ng examples.']

>>> # .strip() removes both leading/trailing whitespace, .rstrip() removes only trailing whitespace, and .lstrip() removes only leading whitespace
>>> # Instead of whitespace, you could also provide an argument to one of these functions to remove instead
>>> dna_string = "AAAAAGTCGAGGTAGCGAAAA"
>>> dna_string.strip("A")
'GTCGAGGTAGCG'
```


### Lists
Lists are defined using brackets (`[]`), and each list item can be any variable type. 

``` python
>>> # This list contains only integers
>>> numbers = [1,2,3,4,5]

>>> # This list contains integers and floats and strings, oh my!
>>> crazy_list = [5, 77.2, -9, "word", -1.32, "more words"]
```

Python lists are incredibly flexible. Like strings, they are ordered, so they support indexing. However, unlike strings, lists are **mutable**, meaning they can be changed! List items can be removed, changed, and new list items can even be added after they are defined. In other words, lists can be changed *in place without (re-)assigning anything*. 

``` python

>>> simple = [1,4,9,2,5,11]
>>> simple[4] # grab the 4th entry in this list
2
>>> simple[15] # what happens when the entry doesn't exist? you get an error message!
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
IndexError: list index out of range

>>> # Rewrite the list in place! You don't need to redefine simple - just modify an element using indexing.
>>> simple[2] = 888
simple
[1,4,888,2,5,11]

>>> # Some useful functions you can use on a list. These will NOT modify the list in any way.
>>> len(simple) # How many elements in a list?
6
>>> sorted(simple) # sort the items in ascending order
[1, 2, 4, 5, 11, 888]
>>> max(simple) # but probably only use this when the list has only numbers. (note that min() is also around!)
888
```

Some useful list methods include,

``` python
>>> # .append(value) adds value to the end of a list, ultimately modifying the list in place!
>>> simple.append(100.33)
simple
[1,4,888,2,5,11,100.33]
>>> len(simple)
7

>>> # .remove(value) removes all entries in a list corresponding to the provided argument, value
>>> simple.remove(11)
[1,4,888,2,5,100.33]

>>> # .insert(index, value) inserts the provided value into the specified index of the list
>>> simple.insert(0, 500)
[500,1,4,888,2,5,100.33]
```


### Dictionary
Dictionaries are defined using braces (`{}`), and they are essentially **unordered** lists of key:value pairs, and they are python's version of "associative arrays." Keys and values can be any type, although typically keys are either integers, floats, or strings. Dictionaries are incredibly useful for storing information; all keys must be unique, but values may be repeated.

```python
>>> taxonomy = {'gecko':'vertebrate', 'human':'vertebrate', 'squid':'mollusk', 'butterfly':'insect', 'oak tree': 'plant'}
```

Dictionaries are indexed using keys. As dictionaries are unordered, the particular order (e.g. gecko, vertebrate, squid..) is not preserved, but the key:value pairs are fixed.
```python
>>> taxonomy["gecko"]
'vertebrate'
>>> # Add a new key:value pair
>>> taxonomy["e. coli"] = "bactera"
>>> taxonomy
{'butterfly': 'insect', 'oak tree': 'plant', 'squid': 'mollusk', 'e. coli': 'bactera', 'human': 'vertebrate', 'gecko': 'vertebrate'}

>>> # the method .keys() pulls up all dictionary keys as a list
taxonomy.keys()
['butterfly', 'oak tree', 'squid', 'e. coli', 'human', 'gecko']

>>> # the method .values() pulls up all dictionary values as a list
['insect', 'plant', 'mollusk', 'bactera', 'vertebrate', 'vertebrate']


>>> # Values can be all kinds of things, even lists!
>>> meals = {"breakfast": ["coffee", "cereal", "banana"], "lunch": ["salad", "lemonade", "chicken fingers"], "dinner": ["steak", "asparagus", "beer", "more beer"], "dessert": ["ice cream", "chocolate sauce", "sprinkles"] }
```

### Tuples

Tuples are basically immutable lists, created with parentheses `()`. We won't use them much, but they're nice to be aware of! Aaand that's all we'll say for now.

### Useful functions

Some fun things we can do with these variables!

```python

>>> # "sep".join(list) will join items in list as a string together, with items separated by "sep"
>>> my_list = ["look", "this", "is", "a", "list"]
>>> " ".join(my_list)
'look this is a list'

>>> # Merge two lists together into a single list of paired tuples with the function zip()
>>> loci = ["locus1", "locus2", "locus3", "locus4"]
>>> snps = ["A/C", "T/G", "G/C", "A/G"]
>>> zip(loci, snps)
[('locus1', 'A/C'), ('locus2', 'T/G'), ('locus3', 'G/C'), ('locus4', 'A/G')]
>>> # Turn these zipped lists into a dictionary!
>>> dict(zip(loci, snps))
{'locus1': 'A/C', 'locus2': 'T/G', 'locus3': 'G/C', 'locus4': 'A/G'}

```














