import sys
import os
usage = "Usage: python inputfile outputfile"


# Input arguments!
n = 1
print sys.argv
print len(sys.argv)
#assert(len(sys.argv) == n), usage 

# We've seen this before
#print os.listdir("../")
# os.getcwd()
# os.chdir()
# os.mkdir() / os.mkdirs()
# os.remove()
# os.rmdir()

# shutil.copy(old, new)
# shutil.move(old, new)


# Run external commands with os.system
command = "mafft " + inputfile + " > " + outputfile
print command
a = os.system("mafft " + inputfile + " > " + outputfile)
assert(a == 0), "Didn't run properly :/"

