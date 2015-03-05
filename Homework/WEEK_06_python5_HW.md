## Homework for the Python 5 Testing module

Write a beautiful script with three functions. We will be using a csv file with a list of expenses per day during a field season.
```
Trip,Date,Price,Item,Paid with
quito,3-Sep,4.00,lunch,cash
quito,3-Sep,1.00,bus,cash
```

 1. The first function should read in the csv file [WEEK_06_python5_HW.csv](WEEK_06_python5_HW.csv), save it as a list, print the header and then remove it, and *return* the dataset as a list.

 2. The second function should compute the amount of money (column 3, Price) spent per trip (column 1, Trip). Check *(not manually!!)* that all of the Prices are floats, and fix them if they are not. Print to a logfile the total amount of money spent per trip.

 3. The third function should be a main() function that calls the other two.

Be sure to use print statements, assertions, and at least and one try-except clause as you write
in order to check your code. 


Bonus:

Write another function that calculates the amount of money spent per day and calculate an average value.
