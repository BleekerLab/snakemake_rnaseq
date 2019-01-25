import os
import sys

args = sys.argv
fa = open(args[1], "r")
name = ""
names = {}
combis = {}
for l in fa:
    if len(l) > 5:
        if l.startswith(">"):
            l = l.split(" ")
            if len(names)>0:
                names[name] = [str(len(names[name])), names[name]]
            name = ".".join(l[1].split(".")[:-1])
            no = l[0][1:]
            combis[no] = name
            names[name] = ""
        else:
            names[name] += l.rstrip()
names[name] = [str(len(names[name])), names[name]]
combis[no] = name
fa.close()

lijst = []
fu = open(args[2], "r")
for l in fu:
    l = l.rstrip().split("\t")
    if combis[l[0]] not in lijst:
        names[combis[l[0]]] += l[2:]
        lijst.append(combis[l[0]])
fu.close()

clusts = open("outFile.txt", "r")
clusters = {}
for line in clusts:
    line = line.rstrip()
    line = line.split("\t")
    if len(line) == 2:
        line.append("NaN")
    clusters[line[0]] = "\t".join(line[1:])
clusts.close()

inFile  = open(args[3], "r")
uitFile = open(args[4], "w")
for l in inFile:
    l = l.rstrip().split("\t")
    if "genes" in l[0]:
        l.append("\t".join([clusters["gene"], "length\tsequence\tlength of hit\te-value\tname and function"]))
    elif l[0] in clusters:
            l.append(clusters[l[0]] + "\t".join(names[l[0]]))
    else:
        if l[0] in names:
            l.append("\t".join(names[l[0]]))
    uitFile.write("\t".join(l)+"\n")
inFile.close()
uitFile.close()
