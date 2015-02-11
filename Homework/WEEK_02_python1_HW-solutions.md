### Python 1 Homework Exercises

**From the command line**, create a new directory for this homework called "python1-hw", enter this directory, and create a new file called "hw1.py". Open this file with your favorite text editor.
Complete the homework and all other exercises in this file. Be sure to use lots of comments! Every time you run a command, be sure to print it to stdout to make sure it works as you expect.


1. Define a variable called `mystring`, which contains a lengthy string of some kind (a random sentence, your address, song lyrics, a haiku, random letters, whatever).

   - Without redefining `mystring`, replace the first 4 occurrences of the letter "a" with the number 6. Use a single string method for this task.
       ```python
          mystring = "Do Koalas Prefer Chocolate Or Fruit, Generally Speaking?" # taxonomy!
          mystring.replace("a", "6", 4)
          mystring
          'Do Ko6l6s Prefer Chocol6te Or Fruit, Gener6lly Speaking?'
      ```
      
   - Redefine `mystring` such that all occurrences of some other letter are replaced with again a different letter (up to you!). Use a single string method for this task.
         ```python
          mystring.replace("r", "R")
          mystring
          'Do Koalas PRefeR Chocolate OR FRuit, GeneRally Speaking?'
      ```
      
   - Use indexing to replace the letter in the 5th position of `mystring` with the letter "X". How did that go?
         ```python
          mystring[4] = "X"
          Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
         TypeError: 'str' object does not support item assignment
      ```       
   *This doesn't work because strings are immutable!*

2. Define a list variable called `fruits` which contains the following 4 entries: banana, apple, grape, and plum.

   - Use the list method `.append()` to add "kiwi" to the end of the list.
   - Use the list method `.insert()` to add "pear" to the end of the list.
   - Use the list method `.insert()` to add "orange" to the **middle** of the list. Try to incorporate the `len()` function as part of your statement.
   - Use indexing to change the 2nd fruit in the list to "peach".

3. Create a variable called `rna` which contains some random RNA (A,C,G,and U) string of a good-sized length.
   - Use the `len()` function to determine the length of this variable.
   - Define a dictionary variable called `rna_dict`. Fill this dictionary with key:value pairs indicating how many A, C, G, and U are in the `rna` string. Try to use `.count()` method in your code!
   - Print a statement that says, "There are (some #) A nucleotides in this sequence." Do **not** use the letter "A" or the actual number as part of your print statement - use the dictionary itself in your printing.
   - Create another variable called `rna2` which contains the same sequence as `rna`, but with a few missing nucleotides ("N") tacked on. For this, do **not** type out (or copy/paste!) the `rna` contents, but rather use the `rna` variable itself and the `+` symbol to create `rna2`.
   - Create a new dictionary for the `rna2` nucleotide counts. Can you do this **without** retyping the same things as you typed to make the first dictionary?
   - BONUS! Use indexing to print the first key:value pair as a tuple. Did that work? No? Then hint! you'll need to use a dictionary method for this!
   
4. Python has many modules, or libraries, which can be imported into a script. Modules are basically pieces of code with lots of convenient functions that aren't normally available to you. For this exercise, we'll use the module "string", which contains lots of useful functions for dealing with string variables. In particular, we'll use the function `maketrans` and the string method `.translate()`.
   At the top of your script, add the line: `import string`. The `maketrans` function is very useful for "re-mapping" individual elements in strings. We'll use `maketrans` to create a translator map to determine the genetic complement of `rna`.
   
   ```python
    >>> # Import the python module "string"
    >>> import string
    >>> # Create a translator variable with the function maketrans. The function takes two string arguments, which should map 1:1.
    >>> # Note that, since this function is in the package string, we have to call it as "string.maketrans" (more on this notation in a few weeks!)
    >>> rna = "AACCUCUCAAGCGCAUCGAUCGA"
    >>> translate_to_complement = string.maketrans("ACGU", "UGCA") 
    
    >>> # Use the .translate() function in conjuction with the translator created with maketrans to determine the complement
    >>> complement = rna.translate( translate_to_complement )
    >>> print complement
    UUGGAGAGUUCGCGUAGCUAGCU
    >>> print rna
    AACCUCUCAAGCGCAUCGAUCGA

   ```
   
   - Now that you've seen an example, create a new translator with `maketrans` and play around with a new string to translate things!




