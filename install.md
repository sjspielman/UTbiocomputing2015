##Hello, new students!

There are a few things we need you to install prior to the course. Below, find your operating system and follow the instructions.

#All operating systems

We will use the anaconda distribution of Python. This version has some very nice support for mathematical and statistical programming. You can choose the version for your operating system here:  http://continuum.io/downloads

For Mac and PC, an installer will pop up, and guide you through the process. See below for Linux instructions.

##Mac

If you are using a Mac, you need a plain-text editor. This is used for editing code and scripts. We recommend TextWrangler:
http://www.barebones.com/products/textwrangler/

Under the Applications folder on your computer, please find the Utilities folder. In this folder, there is an application called 'Terminal'. Please double-click it. When it opens, type

```python
python
```

into the terminal. This should display the following, probably with some text:

```python
>>>
```

If you aren't able to open the terminal, or type 'python' into it, please e-mail me. Type Ctrl+D to exit from Python.

Lastly, Mac users will need a version of Git. It can be downloaded from here:
http://code.google.com/p/git-osx-installer 
Choose your operating system from the list and follow the instructions.

Once you've installed git, in your terminal, type:

```UNIX
git
```

As long as you don't see the following error:

```UNIX
-bash: git: command not found
```

You should be good to go!



##Windows

Windows users need a plain text editor for reading and editing code. We recommend using Notepad++: http://notepad-plus-plus.org/download/v6.4.3.html

Windows doesn't natively have the same terminal functionality of Mac or Linux, so you'll have to install it. Download and install git-bash from this site:

https://openhatch.org/missions/windows-setup/install-git-bash

Open the program by following the instructions here:

https://openhatch.org/missions/windows-setup/open-git-bash-prompt

When you open the program, call python by typing,

```Python

python

```

This should display the following, probably with some text:

```python
>>>
```

Ctrl + D will exit from Python.

Once you have exited Python, type:

```UNIX
git
```

If you do not receive any error messages, you should be set!


##Ubuntu



The native gEdit text editor should be fine for the code reading and writing we will do in this workshop. If you don't already have it, git can be installed in the software center.

Control + Alt + Tab will open a terminal window. Installing the Python Anaconda distribution on Linux is easy. Type 

```UNIX
./
```

and then drag the downloaded Anaconda file into the terminal. Press enter. Note that you may need to run the installer as root. This can be done like so:

```UNIX
sudo ./
```
followed by the dragged and dropped file. If you have permissions issues, you can fix these with 'chmod +x' followed by the filename


If you'd prefer to install git at the command line, you can do so by typing 

```UNIX
sudo apt-get install git
```

Once you've opened the terminal, type:

```UNIX
git
```

If you don't receive an error message saying the command was not found, you should be good to go.

Also, type:


```Python

python

```
This should display the following, probably with some text:

```python
>>>
```

Typing Ctrl+D will exit from python.


##Helpful Hints for all Operating Systems

If you do a new install of Python or git, it is helpful to close down your terminal program and restart it. This allows the terminal to update and use the new installs.

##Installing R

Later in the semester we will spend a few days with the R language. You can download [R] (http://cran.rstudio.com/) as well as a popular [IDE] (http://www.rstudio.com/ide/download/).
