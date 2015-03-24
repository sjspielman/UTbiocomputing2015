from pyvolve import *
from random import randint
from Bio import AlignIO
tree = read_tree(file="tree.tre")
m = Model("nucleotide")
m.construct_model()

for i in range(100):
    print i
    save = []
    for n in range(10):
        p = Partition(size = randint(15, 30), models = m)
        Evolver(partitions = p, tree = tree, seqfile = "temp.fasta", infofile = None, ratefile = None)()
        s =AlignIO.read("temp.fasta", "fasta")
        save.append(str(s[0].seq))
        save.append(str(s[1].seq))
    with open("seqs" + str(i) + ".fasta", "w") as f:
        n = 0
        for entry in save:
            f.write(">taxon" + str(n) + "\n" + entry + "\n")
            n += 1
os.remove("temp.fasta")