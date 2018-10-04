import os
import sys

args = sys.argv

fa = open("../" + args[1], "r")
name = ""
names = {}
combis = {}
for l in fa:
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

fu = open("../" + args[2], "r")
for l in fu:
    l = l.rstrip().split("\t")
    names[combis[l[0]]] = l[1:]
fu.close()

inFile  = open("../" + args[3], "r")
uitFile = open("../" + args[4], "w")
print(names["STRG.447"])
for l in inFile:
    l = l.rstrip().split(",")
    if "Gene" in l[1]:
        print("check")
        l.append("length\tlength of hit\te-value\tname and function")
    else:
        if len(l) > 3:
            l.append("\t".join(names[l[1].strip("\"")]))
    uitFile.write("\t".join(l)+"\n")

