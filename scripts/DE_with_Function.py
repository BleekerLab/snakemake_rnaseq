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
inFile  = open(args[3], "r")
uitFile = open(args[4], "w")

for l in inFile:
    l = l.rstrip().split("\t")
    if "genes" in l[0]:
        l.append("length\tsequence\tlength of hit\te-value\tname and function")
    else:
        if l[0] in names:
            l.append("\t".join(names[l[0]]))
    uitFile.write("\t".join(l)+"\n")
inFile.close()
uitFile.close()
