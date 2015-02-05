#CHEATSHEET PYTHON 2: LOOPS

##Functions 

Function  |  Description  | Example
----------|---------------|--------
`dir()` | shows all methods you can use with a variable | `dir(stringvariablename)` would show string methods
`range()` | makes a list of a certain range of numbers | `range(0,20)` makes a list from 0 to 19, 20 numbers where the first is inclusive and the last is exclusive <br> `range(0,20,2)` makes the same list but counts by 2
`input()` | asks the user for input, stores as variable | `input("Enter a phrase: ")` saves what user inputs as variable x


##Modules 

Module  |  Description  | Example
----------|---------------|--------
`random` | generating fake data | `random.random()` picks a value between 0 and 1 <br> `random.randint(0,10)` randomly picks a value between 0 and 10
`time` | anything to do with timing | `time.asctime( time.localtime(time.time()) )` prints current date and time


##Loops and simple examples

- FOR is for iterating over values in a list, string, file, range, etc
```														
for item in thing:
	do this command
	do that command
	#return to for statement and move to next item
```

Example: print 20 numbers (actually prints 0 to 19)
```python
for i in range(20):
	print i 
```

- IF is the basic decision making tool
```
if logical condition == True: #line leading to an indented block ends with a colon
	do this command  #indenting blocks of code indicates how the program flows
	do this command
continue with main commands
```

Example: check to see if a certain sequence is within a sequence
```python
sequence='AGGGTGTGTCCTGA'
if 'AGG' in sequence:
	print 'I found your sequence'
```

- IF ELSE is used for decision making in an either or context
```
if logical condition == True: 
	do this command  
	do this command
else:
	do this command
```

Example: Test whether sequence is RNA or DNA by seeing if there is a uracil
```python
seqs=['AUUGACAUCGAUCGA','AGACTGATCGATCTAG']
for seq in seqs:
	if 'U' in seq:
		print '%s is RNA' %seq
	else:
		print '%s is DNA' %seq
```


- IF ELIF is useful when you have more than one condition to check before deciding
```
if logical condition == True: 
	do this command  
	do this command
elif other logical condition == True:
	do this command
else:
	do this command
```

Example: take the above step one step further
```python
seqs=['AUUGACAUCGAUCGA','AGACTGATCGATCTAG','JIEONONE']
for seq in seqs:
	if 'U' in seq:
		print '%s is RNA' %seq
	elif 'T' in seq:
		print '%s is DNA' %seq
	else:
		print 'there are no Ts or Us in %s' %seq

- WHILE is useful for checking input types, and when your value in the conditional might change within the loop
```
while condition == True:
	do this command
	do that command
continue with normal commands
```

Example: check whether user input fits a certain criterion

```python
x=int(input("Give me a number from 1 to 10: ")) #accepts user input, converts to an integer
while x not in range(0,11): #if input was not number between 1 and 10
        print "That's not a number between 1 and 10 :(" #prints a complaint
        x=int(raw_input("Give me a number from 1 to 10: ")) #asks for a new number
print "Thanks!" #thanks the user for being a nice user
```


- COMPREHENSIONS are faster than for loops and good to use when you're sifting through a long list or range

```
list1=[] #make empty list
for item in thing:
	list1.append(item)	

list1 = [item for item in thing]
```

Example: make a list of the lengths of sequences in seqs list
```python
seqs=['AUUGACAUCGAUCGA','AGACTGATCGATCTAG']
l=[length(seq) for seq in seqs]
```

- FOR LOOP IN BASH

```for x in list; do command to variable x; done
```
Example: cat all fasta files in your current directory into a new file
```
for i in *.fasta; cat i >> newseqs.fasta; done 
```
