# Introduction to R

R is a high-level programming language most typically used for statistical analysis and data manipulation.

**Why (when) to use R:**
* R is free and easy, many people use it which facilitates collaboration (though the same applies to Python)
* R has a huge (and constantly growing) base of peer-contributed packages and statistical routines
* R is an excellent platform for simulation and numerical routines
* R has a lot of useful syntax for data manipulation
* R has great graphics capabilities
* R interfaces very well with C++ (via Rcpp) if you need speed
* R markdown/Sweave/Shiny make it really easy to combine text, presentations, analysis, and output in a single file
	
**Why/when not to use R (in my experience):**
* R is not great for string manipulation (you can certainly do it but it's more cumbersome than in Python)
* R has little support for symbolic algebra (use SymPy or Sage or Mathematica)
* R is relatively slow (but see above about Rcpp)
* R is not meant to create stand-alone programs (ie. a graphical interface outside the console)

This lesson is not meant to be a comprehensive survey of syntax, or of statistical routines (which would take waaaay to long); but instead is meant to introduce some important functionalities and topics.

There is a helpful (if patronizing) interactive tutorial series called [Swirl](http://swirlstats.com/). I recommend it for beginners to programming (others will find it rather annoying).

In the following code snippets, anything following `#@` is output, and anything following `#` is a comment.

### Objects and assignment

R is an object oriented programming language. Assign values to objects with `=` or `<-`.
```r
x <- 1.3
fx1 <- cos(2*pi*x)
x <- 2.6
fx2 <- sin(2*pi*x)
y <- "yakkity yak" # text string
```
Here I am assigning some values to objects. `fx1` and `fx2` store the result of a trigonometric function applied to another number, which is stored in `x`. The format for calling a function is `function_name(function_arguments)`. I'll talk about functions in more detail later.

You can see the contents of an object by typing its name.
```r
fx1
	#@ -0.309017
y
	#@ "yakkity yak"
```

Objects have attributes. The most basic attribute of an object is its class(es).
```r
class(fx2)
	#@ "numeric"

class(y)
	#@ "character"
```

You will commonly see `numeric`, `character`, and `factor` classes. Think of factors like integers (ie. having an order), but where each integer has an arbitrary name.
```r
watterson <- factor(c("Calvin", "Hobbes", "Susie"))

watterson
	#@ [1] Calvin Hobbes Susie 
	#@ Levels: Calvin Hobbes Susie
	
as.numeric(watterson) ## the integer identity is revealed
	#@ [1] 1 2 3
```

Why use classes? Classes have intrinsic functions called *methods*; ie. functions which share the same name but which vary their behavior based on the class of the object that the function is applied to.
```r
## levels() shows levels of a factor
levels(watterson)

## but doesn't work with a character vector
levels(c("Calvin", Hobbes", "Susie"))
```

### Data structures

*Vectors* are concatenations of values of the same class. Concatenations are formed with the function `c()`.
```r
myVector <- c(fx1, fx2)
myVector
	#@ [1] -0.3090170 -0.5877853
```

Note that a concatenations of values of different classes will be coerced to to a common class.
```r
myVector <- c(y, fx1, fx2)
myVector 
	#@ [1] "yakkity yak" "-0.309016994374947" "-0.587785252292473"
# notice the parentheses, indicating that the numbers have been converted to strings
```

Access elements of a vector by square brackets `[]`.
```r
myVector[3]
	#@ "-0.587785252292473"
```
Here, we access the third element of `myVector`. I'll explain indexing in more detail later.

Vectors can have names associated with elements.
```r
names(myVector) <- c("chip", "and", "dale")
myVector
	#@          chip                  and                 dale 
	#@ "yakkity yak" "-0.309016994374947" "-0.587785252292473" 

myVector["chip"] # can access element by name
	#@          chip 
	#@ "yakkity yak" 
```

*Matrices* (and more generally *arrays*) are multi-dimensional generalizations of vectors. All the elements of an array all the same value, but in the shape of a rectangle, cube, etc.
```r
myMatrix <- matrix(c(1,2,3,4,5,6), nrow=2, ncol=3)
myMatrix
	#@      [,1] [,2] [,3]
	#@ [1,]    1    3    5
	#@ [2,]    2    4    6
	
myArray <- array(c(1,2,3,4,5,6,7,8), dim=c(2,2,2)) 
myArray
	#@ , , 1
	#@ 
	#@      [,1] [,2]
	#@ [1,]    1    3
	#@ [2,]    2    4
	#@ 
	#@ , , 2
	#@ 
	#@      [,1] [,2]
	#@ [1,]    5    7
	#@ [2,]    6    8
```

Access for matrices and arrays is similar vectors, but distinguish between dimensions with commas, `[,]` (ie. rows, columns)
```r
myMatrix[2,2]
	#@ [1] 4

myArray[1,2,2]
	#@ [1] 7
```

Leaving a dimension blank will extract all elements in that dimension.
```r
myMatrix[2,]
	#@ [1] 2 4 6

myArray[,,2]
	#@      [,1] [,2]
	#@ [1,]    5    7
	#@ [2,]    6    8
```

Matrices have linear algebra operators, for example for a given matrix **X** we might want to calculate the cross product **X'X**.
```r
t(myMatrix) %*% myMatrix
	#@      [,1] [,2] [,3]
	#@ [1,]    5   11   17
	#@ [2,]   11   25   39
	#@ [3,]   17   39   61
```

The most general class of data structure is a *list*.
```r
myList <- list()
```

Lists can contain objects of arbitrary class and length. For example, a list could contain a mixture of lists and single values. Assignment/access from lists can be by name or position of the element(s).
```r
## named assignment
myList$myValue <- fx1

## numerical assignment
myList[[2]] <- fx2

myList
	#@ $myValue
	#@ [1] -0.309017
	#@ 
	#@ [[2]]
	#@ [1] -0.5877853

## access by name (if present) or number
myList$myValue
	#@ [1] -0.309017

myList[["myValue"]]
	#@ [1] -0.309017

myList[[1]]
	#@ [1] -0.309017

```
Assigning to a name which does not already exist, assigns that name in the first available numerical position in the list.

A data frame is a list of vectors, where each vector has the same length, but can have different classes.
```r
myDataFrame <- data.frame(darn=c(1,2,3,4), geez=c(5,6,7,8))
myDataFrame
	#@   darn geez
	#@ 1    1    5
	#@ 2    2    6
	#@ 3    3    7
	#@ 4    4    8

myDataFrame$blarg <- c(1,2) ## when the length of a new entry is a factor of the row number, it will automatically expands to match length
myDataFrame
	#@   darn geez blarg
	#@ 1    1    5     1
	#@ 2    2    6     2
	#@ 3    3    7     1
	#@ 4    4    8     2

## when the length of a new entry is not a factor of the row number, it fails
myDataFrame$yonk <- c(1,2,3) 
	#@  Error in `$<-.data.frame`(`*tmp*`, "yonk", value = c(1, 2, 3)) : 
	#@   replacement has 3 rows, data has 4
``` 

Access individual elements of a data frame by name or double bracketed number (extracts column), by single bracket number (extracts element), or by brackets with comma (pulls out element in row, column)
```r
myDataFrame$darn
	#@ [1] 1 2 3 4

myDataFrame[[3]] ## third column
	#@ [1] 1 2 1 2

myDataFrame[1,2] ## first element of second column
	#@ [1] 5

colnames(myDataFrame) ## show column names
	#@ [1] "darn"  "geez"  "blarg"
```

View the structure of a data container with `str()`:
```r
str(myList)
	#@ List of 2
	#@  $ myValue: num -0.309
	#@  $        : num -0.588

str(myDataFrame)
	#@ 'data.frame':	4 obs. of  3 variables:
	#@  $ darn : num  1 2 3 4
	#@  $ geez : num  5 6 7 8
	#@  $ blarg: num  1 2 1 2
```

### Functions
All functions in R follow the same format: `function_name(function_argument_one, function_argument_two, etc.)`. 

For example, the functions `sin` and `cos` take a single argument: the value (or collection of values) to transform. The function `data.frame` takes an unlimited number of arguments, each of which is a column in the data frame.

You will often want to combine various functions and operations into a new function, which you can then repeatedly apply to data.
```r
## syntax for creating a function
myFunction <- function(arguments, more arguments){
	x <- whatever(crap, you, want)
	return(x)
}
```

For example, create a function to return the mean, standard error of the mean, and confidence interval of a vector:
```r
sumVec <- function(x, alpha){
	n <- length(x) # length of the vector
	se <- sd(x)/sqrt(n) # std error
	mn <- mean(x) # mean of x
	lo <- mn + qnorm(alpha/2)*se # low bound of confidence interval
	hi <- mn + qnorm(1-alpha/2)*se # hi bound of confidence interval
	return( c(mn=mn, se=se, lo=lo, hi=hi) )
}

## apply function to 100 random normal variates
set.seed(99)
normal_sample <- rnorm(100, 12, 2)
sumVec(normal_sample, 0.01)
	#@         mn         se         lo         hi 
	#@ 11.7919582  0.1801482 11.3279272 12.2559892
```

The arguments to a function can be named or not (in which case they are assumed to follow the order of arguments in the function definition).
```r
sumVec(x = normal_sample, alpha = 0.01)
	#@         mn         se         lo         hi 
	#@ 11.7919582  0.1801482 11.3279272 12.2559892
```

You can see the source code for a function by typing its name.
```r
## for example, the function to fit a linear regression model
lm
	#@ function (formula, data, subset, weights, na.action, method = "qr", 
	#@    model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE, 
	#@    contrasts = NULL, offset, ...) 
	#@ {
	#@    ret.x <- x
	#@    ret.y <- y
	#@    cl <- match.call()
	#@    mf <- match.call(expand.dots = FALSE)
	#@    m <- match(c("formula", "data", "subset", "weights", "na.action", 
	#@        "offset"), names(mf), 0L)
	#@    mf <- mf[c(1L, m)]
	#@ <snip>
```

### Packages

R has a vast array of packages. Packages live in a repository (GitHub is one such). [CRAN](http://cran.r-project.org/) is the primary repository for R packages. [Bioconductor](http://www.bioconductor.org/) also has many packages for bioinformatic analyses. There are also many in-progress packages on R-forge.

Install packages from CRAN with the command `install.packages(package_name)`. Load packages with `library(package_name)`
```r
install.packages("plyr", repos="http://cran.us.r-project.org")
install.packages("sos", repos="http://cran.us.r-project.org")
install.packages("ggplot2", repos="http://cran.us.r-project.org")
```
Another method is needed to install packages from Bioconductor (take a look at the website).

### Help

To see help regarding a specific function (note that the package that the function is in must be loaded!) do `?function_name`
```r
?apply
```

To search for a function across packages which are installed, do `??function_name`
```r
??coxph
```

To search using keywords for a specific type of function across the CRAN repository, use the function `findFn(keywords)` in the package `sos`.
```r
library(sos)
findFn("quantile regression")
findFn("repeat elements")
```

And of course, [stackoverflow](http://stackoverflow.com/questions/tagged/r).

### Debugging

Flag a function for debugging with `debug(function_name)`. When you run the function, it will pause before executing each line. To stop this behavior, use `undebug(function_name)`.
```r
sumVec( c(0:15), "a") # will fail because the confidence level needs to be numeric
	#@ Error in alpha/2 : non-numeric argument to binary operator

debug(sumVec)
sumVec( c(0:15), "a")
undebug(sumVec)
```
Use `n` to skip to the next line, `c` to complete a loop, and `Q` to quit debugging mode. Note that if the function calls other functions (and these are throwing the error), you'll also have to flag these for debugging.

Also, you can set R's error-handling behavior to enter debugging mode when something fails:
```r
option(error=recover) ## to set error recovery
sumVec( c(0:15), "a")
option(error=NULL) ## to stop error recovery
```

### Scope
It is important to briefly mention scope in the context of functions. When you save objects, they exist in your working environment. You can see the objects of your working environment by using `ls()`.
```r
ls()
```

When you run a function, it will first look for objects in among its inputs, and then in the working environment.
```r
y = 7

scopeTest <- function(x){
	return( x + y )
}

scopeTest(1)
	#@ [1] 8
```

### Reading in, writing out, and indexing data

First, we need to move into the working directory with our data. This is where files will get read in from, and written out to.
```r
## for example ... replace the following with a path specific to your computer, or use the R/RStudio GUI menu.
setwd("~/Dropbox/School/computing2015") 
```

It's most convenient to input data as a delimited table. Use the function `read.table(file, options)` to read in a delimited table as a data frame.
```r
# reads data from a comma-separated table with a header
Adoxo_count <- read.table("Adoxophyes_counts.csv", sep = ",", header = TRUE)
Adoxo_temp <- read.table("Adoxophyes_temp.csv", sep = ",", header = TRUE)
str(Adoxo_count)
	#@ 'data.frame':	2754 obs. of  3 variables:
	#@ $ day              : int  64 69 74 79 84 89 95 100 105 110 ...
	#@ $ year             : int  1961 1961 1961 1961 1961 1961 1961 1961 1961 1961 ...
	#@ $ Adoxophyes_honmai: int  0 0 1 3 0 0 4 0 0 0 ...
	
head(Adoxo_temp)
	#@   year day temperature
	#@ 1 1960   1       5.800
	#@ 2 1960   2       7.775
	#@ 3 1960   3       9.750
	#@ 4 1960   4      11.725
	#@ 5 1960   5      13.700
	#@ 6 1960   6       0.700
```
The argument `header = TRUE` indicates that the first row is column names.

The function `write.table(data, file, options)` takes `data` and writes it into `file`.
```r
# writes data out to a comma-separated table without the numeric row names
write.table(Adoxo_count, "Adoxo_count_copy.csv", sep = ",", row.names=F)  
```


##### Bracket notation for indexing
The following examples index a data frame, but the same applies to vectors/matrices/arrays.

Bracket notation is [integer position in 1st dimension, in 2nd dimension, in 3rd dimension, etc.] 

If one dimension is left blank, the entire row/column/whatever is returned. Note that R starts indexing at 1, not zero!
```r
Adoxo_temp[4,] # 4th row
	#@   year day temperature
	#@ 4 1960   4      11.725
```

A slice is a set of indices, for example the set \[ \{1,2,3,4,5\} \] is given by `1:5` in R, which can be used to extract elements 1 through 5:
```r
1:5
	#@ [1] 1 2 3 4 5

Adoxo_temp[4,1:2] # 4th row, elements from columns 1 through 2
	#@   year day
	#@ 4 1960   4
```

Slicing does not have to be contiguous: integers pieced together with `c()` are valid slices. Order matters, and elements can be duplicated.
```r
Adoxo_temp[4,c(1,3)] # 4th row, elements from columns 1 and 3
	#@   year temperature
	#@ 4 1960      11.725

subs <- c(2:1, 3, 519, 2) # 2nd through first element, 3rd element, 519th element, 2nd element again
Adoxo_temp[subs,] # replicates the rows listed in subs in the order given
	#@     year day temperature
	#@ 2   1960   2       7.775
	#@ 1   1960   1       5.800
	#@ 3   1960   3       9.750
	#@ 519 1961 153      17.300
	#@ 2.1 1960   2       7.775
```

Introducing a negative symbol before a vector of index numbers will drop these elements.
```r
## drop first column. 
## we use head because otherwise it would flood your console
head(Adoxo_temp[,-c(1)]) 
	#@   day temperature
	#@ 1   1       5.800
	#@ 2   2       7.775
	#@ 3   3       9.750
	#@ 4   4      11.725
	#@ 5   5      13.700
	#@ 6   6       0.700
``` 

Indexing is recursive.
```r
Adoxo_temp[-c(1),][2,] ## pulls out third row, because we dropped the first row, then asked for the second row of those that remain.
	#@   year day temperature
	#@ 3 1960   3        9.75
```


##### Logical indexing
We can evaluate logical statements using `==` (equals) `!=` (does not equal) `>=` (greater than or equal to) `|` (or) `&` (and).
```r
1==1
	#@ [1]  TRUE

c(1,2,3)==1
	#@ [1]  TRUE FALSE FALSE
```

If the logical evaluation is put in brackets, R will pull out those elements which are `TRUE`.
```r
Adoxo_head <- Adoxo_count[1:10,]
Adoxo_head
	#@    day year Adoxophyes_honmai
	#@ 1   64 1961                 0
	#@ 2   69 1961                 0
	#@ 3   74 1961                 1
	#@ 4   79 1961                 3
	#@ 5   84 1961                 0
	#@ 6   89 1961                 0
	#@ 7   95 1961                 4
	#@ 8  100 1961                 0
	#@ 9  105 1961                 0
	#@ 10 110 1961                 0

Adoxo_head$Adoxophyes_honmai == 3
	#@ [1] FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE

Adoxo_head[Adoxo_head$Adoxophyes_honmai == 3, ]
	#@   day year Adoxophyes_honmai
	#@ 4  79 1961                 3

Adoxo_head[Adoxo_head$Adoxophyes_honmai != 1 & Adoxo_head$Adoxophyes_honmai > 0, ]
	#@   day year Adoxophyes_honmai
	#@ 4  79 1961                 3
	#@ 7  95 1961                 4
```

The function `which()` returns the index numbers of those elements which match the logical evaluation,
```r
## note that the fourth and seventh elements have a value that is greater than 1
Adoxo_head$Adoxophyes_honmai
	#@ [1] 0 0 1 3 0 0 4 0 0 0
which(Adoxo_head$Adoxophyes_honmai > 1)
	#@ [1] 4 7
Adoxo_head$Adoxophyes_honmai[-which(Adoxo_head$Adoxophyes_honmai > 1)]
	#@ [1] 0 0 1 0 0 0 0 0
```

The function `match(x,y)` goes through each element of `x` and finds the index of the first value in `y` which matches that element. For example, we have two data frames of different length: the count data is only on some days of the year, while the temperature data is on all days of the year. We might want to extract the temperature values on those days which match the days in the count data.
```r
## to make things easier we'll only work with a single year.
## notice I use logical indexing to pull out only those data from the year 1961
Adoxo_count_1961 <- Adoxo_count[Adoxo_count$year == 1961, ]
Adoxo_temp_1961 <- Adoxo_temp[Adoxo_temp$year == 1961, ]

## just to make the example clearer, I'll scramble the order of the rows in the temperature data
set.seed(101)
Adoxo_temp_1961 <- Adoxo_temp_1961[sample(1:nrow(Adoxo_temp_1961)),]

## match returns a bunch of indices -- the number of which matches the number of rows in the count data
match(Adoxo_count_1961$day, Adoxo_temp_1961$day)
	#@ [1] 211  37  23 249 258 106 150 159 295 137 253 ... <snip>

## if I use these indices on the days in the temperature data ...
Adoxo_temp$day[match(Adoxo_count_1961$day, Adoxo_temp_1961$day)]
	#@ [1]  64  69  74  79  84  89  95 100 105 110 ... <snip>

## I get the same set of days as is given in the count data.
Adoxo_count_1961$day
	#@ [1]  64  69  74  79  84  89  95 100 105 110 ... <snip>
	
## Therefore I can add a temperature column to the count data very easily.
Adoxo_count_1961$temp <- Adoxo_temp_1961$temperature[match(Adoxo_count_1961$day, Adoxo_temp_1961$day)]
```

Note that missing values--denoted `NA` behave differently. You need to use `is.na()` for logical evaluations of missing values.
```r
NA==1
	#@ [1] NA

1==c(1,NA,4)
	#@ [1]  TRUE    NA FALSE

is.na(c(1,NA,4))
	#@ [1] FALSE  TRUE FALSE
```

Here are some other useful logical evaluations which operate on data containers.
```r
"monkey" %in% c("cow", "monkey", "Homo sapiens")
	#@ [1] TRUE

any(c(1,2,3,4,5) > 10)
	#@ [1] FALSE

all(c(1,2,3) < 5)
	#@ [1] TRUE
```

### Flow control and vectorization

Often we wish to repeat the same action, or apply it to multiple values/objects/whatever. Loops have the syntax `for(iterator in vector){ do stuff for each iterator }`.
```r
for(i in 1:5){ print(i) }
	#@ [1] 1
	#@ [1] 2
	#@ [1] 3
	#@ [1] 4
	#@ [1] 5
```

Note the curly brackets: as for functions, the curly brackets contain the operation you wish to perform at each iteration of the loop. `if` and `else` statements also require curly brackets:
```r
for(i in c(3,19,12,10)){
	if(i == 12) {
		print("whoops")
	} else {
		print(i)
	}
}
	#@ [1] 3
	#@ [1] 19
	#@ [1] "whoops"
	#@ [1] 10
```

Looping is extremely slow in R. Whenever possible, we wish to vectorize (perform an operation on all the elements of a data container simultaneously). Most operators and many functions in R are vectorized.
```r
c(1:10) + 5
	#@ [1]  6  7  8  9 10 11 12 13 14 15
c(1:10) + c(1:10) ## when objects have the same dimensions, the operations are done 'element-wise'
	#@ [1]  2  4  6  8 10 12 14 16 18 20
```

Many functions (not just mathematical operators) are vectorized:
```r
dnorm(Adoxo_temp$temperature[1:5], log=T) # evaluates standard normal pdf for each of the first 5 values in the temperature data
	#@ [1] -17.73894 -31.14425 -48.45019 -69.65675 -94.76394
```

The function `apply(data, margin, function, arguments)` is faster than looping over rows or columns, and automatically return the results in a convenient form. 
```r
## we put the temperature data into column form for the first few years
Adoxo_temp_wide <- data.frame(y1961 = Adoxo_temp$temp[Adoxo_temp$year == 1961], 
	y1962 = Adoxo_temp$temp[Adoxo_temp$year == 1962], 
	y1963 = Adoxo_temp$temp[Adoxo_temp$year == 1963])
head(Adoxo_temp_wide)
	#@   y1961 y1962  y1963
	#@ 1  1.30  2.40  5.000
	#@ 2  3.35  4.45  6.475
	#@ 3  5.40  6.50  7.950
	#@ 4  7.45  8.55  9.425
	#@ 5  9.50 10.60 10.900
	#@ 6  3.70  2.30 -0.500

## rather than looping like this ...
quants <- matrix(ncol=ncol(Adoxo_temp_wide), nrow=2)
for(i in 1:ncol(Adoxo_temp_wide)){
	quants[,i] <- quantile(Adoxo_temp_wide[,i], probs=c(0.025,0.975))
}
quants
	#@        [,1]  [,2]  [,3]
	#@ [1,]  1.720  2.21  0.71
	#@ [2,] 30.192 30.00 29.90

## use apply() function
quants <- apply(Adoxo_temp_wide, MARGIN=2, quantile, probs=c(0.025,0.975))
quants
	#@        y1961 y1962 y1963
	#@ 2.5%   1.720  2.21  0.71
	#@ 97.5% 30.192 30.00 29.90
```
The argument `MARGIN=1` specifies rows (think: rows, columns, etc.). The margin argument generalizes to any dimension. `apply()` works on data frames and vectors/arrays/matrices, and outputs an array.

`lapply(list, function, function arguments)` is the analogue of `apply()` for lists. It takes a list as input, and outputs a list.
```r
Adoxo_temp_list <- list(y1961 = Adoxo_temp$temp[Adoxo_temp$year == 1961], 
	y1962 = Adoxo_temp$temp[Adoxo_temp$year == 1962], 
	y1963 = Adoxo_temp$temp[Adoxo_temp$year == 1963])

# calculate mean for each element of the list
lapply(Adoxo_temp_list, mean)
	#@ $y1961
	#@ [1] 17.28863
	#@ 
	#@ $y1962
	#@ [1] 16.74288
	#@ 
	#@ $y1963
	#@ [1] 16.55151
```

### Split-apply-combine

Often we want to break our data into peices, do the same routine on each piece, and put the results back together in some form (like a table).

Examples:
* Break data apart and modify each piece in some way
* Break data apart and calculate summary statistics on each piece
* Break data apart and fit each piece to a model

The `apply()` function discussed earlier is applicable to this situation, but a more general solution is implemented in the package `plyr`.

Basic syntax: `**ply(data, how to split data up, function to apply to pieces of data)`

The stars are replaced with letters. The first star indicates the format of the data, the second star indicates the format of the output.

|Letter|Format|
|------|------|
|`a`|array|
|`d`|data frame|
|`l`|list|
|`_`|nothing|

So `aaply()` takes an array and returns an array, `alply()` takes an array and returns a list, `a_ply()` takes an array and returns nothing. The underscore `_` is used for the output format only; for example if we want to write out data to the harddrive but not to reassemble it back into an object. 

***Example of split-apply-combine:***

With the *Adoxophyes honmai* count data, we might want to split data up by year and calculate summary statistics. Some reasonable summary statistics might be the mean count, the ratio of zero counts to positive counts, and number of observations in the year. First, we write a function which will take a data frame with the *exact same structure* as the `Adoxo_count` data frame, and will calculate summary statistics on it.
```r
sum_func <- function(x){
	x <- x[!is.na(x$Adoxophyes_honmai),] # remove rows with missing values
	mn <- mean(x$Adoxophyes_honmai) # calculate mean
	zerosToPositives <- sum(x$Adoxophyes_honmai == 0)/sum(x$Adoxophyes_honmai > 0)
	num_obs <- nrow(x) # number of observations
	return( c(mean = mn, zero_ratio = zerosToPositives, number_observations = num_obs) )
}
```

If we apply `sum_func()` to the whole data frame, we get a summary for the whole data frame.
```r
sum_func(Adoxo_count)
	#@        mean          zero_ratio number_observations 
	#@ 105.8982495           0.1223905        2742.0000000
```

If we use `sum_func()` within `ddply()`:
```r
ddply(Adoxo_count, .(year), sum_func)
   	#@    year      mean zero_ratio number_observations
	#@ 1  1961  23.46296 0.28571429                  54
	#@ 2  1962  39.77778 0.38461538                  54
	#@ 3  1963  64.94444 0.17391304                  54
	#@ 4  1964  45.48148 0.20000000                  54
	#@ 5  1965  87.27778 0.10204082                  54
	#@ 6  1966  54.51852 0.05882353                  54
	#@ 7  1967  73.48148 0.10204082                  54
	#@ 8  1968 157.29630 0.01886792                  54
	#@ 9  1969  51.29630 0.12500000                  54
	#@ 10 1970 130.29630 0.17391304                  54
	#@ <snip>
```

Note that `sum_func()` takes a single argument `x`. In the context of `ddply`, `x` will be a data frame which is a subset of the entire data frame `Adoxo_count` -- a subset which represents a particular year. Also note that with data frames, the syntax `.(variable, another variable)` in the second argument to `**ply` indicates what variables to use to split up the data frame. If we had divided the time series into years and months, we could use `.(year, month)` to calculate summary statistics for each year-month combination in the data. 

The above is a very simple example; a more complicated application would be to fit a spline regression to the temperature data within each year, and return the fitted regression line.
```r
splineFit <- function(x){
	library(splines)
	fit <- lm(temperature ~ bs(day,8), data = x)
	return( data.frame(day=x$day, fit=fitted(fit) ) )
}

splineOut <- ddply(Adoxo_temp, .(year), splineFit)

head(splineOut)
  	#@   year day      fit
  	#@ 1 1960   1 9.179364
  	#@ 2 1960   2 8.901868
  	#@ 3 1960   3 8.641705
  	#@ 4 1960   4 8.398524
  	#@ 5 1960   5 8.171974
  	#@ 6 1960   6 7.961704
```

### Graphics

There are two main platforms for graphics in R: `base` and `grid`. Graphics are a huge topic with many associated packages. I'll only give a brief introduction to `ggplot2` which is probably the most popular package for R graphics.

`ggplot2` works with layers. For each layer, you need data, a description of the type of layer (points, lines, whatever), and a mapping from the data onto the graphics objects (like, what variable will go on the x-axis, a color scheme which describes grouping in the data, etc.). Then, for all layers, we need a description of the formatting of the figure as a whole (ie. font size, grid lines, etc.)

In `ggplot2` lingo:
The description of the type of layer is called a `geom`
The mapping of data to figure is called an aesthetic, or `aes`
The general formatting of the figure is called a theme

You can save each layer, or a combination of layers, as an object. If formatting options/mappings/data are specified for a layer, they will override the overall formatting options.

First, we start with a base layer, which will give the base data and base aesthetic to be used (unless overridden in a layer).
```r
Adoxo_temp$yearCat <- factor(Adoxo_temp$year) # treat year as categorical not continuous
myPlot <- ggplot(data = Adoxo_temp[Adoxo_temp$year > 1990,], aes(x = day))
```

The syntax `aes(x = day)` says to map `day` onto the x-axis.

Now, we add a point layer.
```r
myPlot <- myPlot + geom_point(aes(y=temperature, colour = yearCat), alpha=0.5, size = 1)
myPlot
```
The syntax `aes(y = temperature, colour = yearCat)` says to map `temperature` onto the y-axis, and to map `yearCat` onto colour (ie. colour-code years). `alpha` is an option that controls transparency. Note the syntax `geom_whatever()`

We add a line layer, using the results of the Fourier series in the `fourierOut` data table (note we have to specify that we are using different data in this layer).
```r
splineOut$yearCat <- factor(splineOut$year)
myPlot <- myPlot + geom_line(data=splineOut[splineOut$year > 1990,], aes(y = fit, colour = yearCat))
## plot with regression fit
myPlot
```

Currently, all the years are plotted together. We can split these into 'facets' (more `ggplot2` lingo):
```r
myPlot <- myPlot + facet_grid(~year)
## show plot split into facets
myPlot
```

We can change the overall formatting to a canned theme:
```r
myPlot <- myPlot + theme_minimal()
## plot with a less gray-heavy color scheme
myPlot
```

We can change the formatting beyond a canned theme, for example by adding a border around the panels.
```r
myPlot <- myPlot + theme(panel.border = element_rect(fill=NA))
## plot with border around facets
myPlot
```

Note that we could have done this in one line, adding all layers together.
```r
myPlot <- ggplot(data = Adoxo_temp[Adoxo_temp$year > 1990,], aes(x = day)) + geom_point(aes(y=temperature, colour = yearCat), alpha=0.5, size = 1) + geom_line(data=splineOut[splineOut$year > 1990,], aes(y = fit, colour = yearCat)) + facet_grid(~year) + theme_minimal() + theme(panel.border = element_rect(fill=NA))
```

Note: there is a function `qplot()` which provides a quick syntax for many basic plot types, without the layer-on-layer syntax. Many people learn using `qplot` because the syntax is closer to base R graphics. 
I don't recommend learning using this, because `qplot` is a less flexible syntax. Learn the long, arduous syntax first.

There are many nuances to `ggplot2`. You can easily create your own geoms, add custom graphical elements using package `gtable`, etc.  See the `ggplot2` online documentation.

