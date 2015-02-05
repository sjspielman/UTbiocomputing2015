#Python 2 Homework due Feb 18

1. Create your own python script that creates a genome with equal base frequencies of a *user-requested* size. Hint: use `input()`, `random`, and some of the text from the lesson.

2. Now rehash that code into a for loop that allows you to create 30 sequences of length 100. Save them in a list.

	- Using a for loop and string indexing, print the position of the first T for each sequence in a readable format (use %d with the print statement).
	- Now add a variable to save this number within the loop rather than printing it.
	- Use this variable to cut the part of the sequence that comes before that index out of each sequence. 
	- Repeat the last step but add an if statement so that it only cuts the sequence if the T is one of the first three nucleotides.
	
3. Write a for loop using your sequences that checks whether the sequence 'AAAA' is in them, and if it is, adds one to a counter. At the end of the loop, use the print statement with %d to state how many of your sequences had this subsequence.

