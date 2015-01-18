# Python I

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

All variables have a certain **type**. The variable type deterines what can be done with the variable itself. Some basic types include, 

Variable Type   | Description | Casting
---------|--------------|---------
integer | whole number  | int()       
float   | decimal number | float()
string  | ordered, immutable character container | str()
list    | ordered, mutable container | list()
tuple   | ordered, immutable container | tuple()
dictionary | unordered, mutable container | dict()

Now, let's go into each variable type in depth.

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
```




