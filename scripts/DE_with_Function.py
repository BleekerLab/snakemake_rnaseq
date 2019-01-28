from optparse import OptionParser 

# get the command line arguments
parser = OptionParser(description="Script that combines the normalized and differential expression data as outputted by DESeq2.R, the clusters as outputted by plotscript.R and the results from the blast against the reference proteome. Resulting in a combined tab-delimeted data file")
parser.add_option('-r', '--result_file', 
                    type=str,
                    default="results/result.csv",
                    metavar="",
                    help = "path and name of of the input file, being the output file of the R-script DESeq2.R, default = results/result.csv")
parser.add_option('-f', '--fasta_file', 
                    type=str,
                    default="genome/stringtie_transcriptome.fasta",
                    metavar="",
                    help = "path and name of the stringtie transcriptome fasta file, a file containing the fastas of the de novo generated stringtie transcriptome, default = result/helperFile.csv")
parser.add_option('-b', '--blast_file', 
                    type=str,
                    default="results/stringtie_transcriptome_blast.txt",
                    metavar="",
                    help = "path and name of of the blast file, a tab delimited file containing the blast result for the de novo transcriptome, default = results/result.csv")
parser.add_option('-c', '--cluster_file', 
                    type=str,
                    default="result/clusters.csv",
                    metavar="",
                    help = "path and name of the clusters file, a tab delimited file containing cluster numbers and depending on the method of clustering correlation values or membership values., default = result/helperFile.csv")
parser.add_option('-o', '--output_file', 
                    type=str,
                    default="result/final.txt",
                    metavar="",
                    help = "path and name of the output file. A tab delimimited file containing normalised reads, differential expressions and (adjusted)p-values, number of clusters the genes partisioned to and the hypothetical function as found by a blast. default = result/final.txt")
parser.add_option('-p', '--working_dir', 
                    type=str,
                    default="temp/mapped/",
                    metavar="",
                    help = "path to bam files, to be removed from the sample names in the header of the final output file. default = temp/mapped/")


(options, args) = parser.parse_args()

# get combination of numbers and transcript name.
fa = open(options.fasta_file, "r")
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

# get 
lijst = []
fu = open(options.blast_file, "r")
for l in fu:
    l = l.rstrip().split("\t")
    if combis[l[0]] not in lijst:
        names[combis[l[0]]] += l[2:]
        lijst.append(combis[l[0]])
fu.close()

clusts = open(options.cluster_file, "r")
clusters = {}
for line in clusts:
    line = line.rstrip()
    line = line.split("\t")
    if len(line) == 2:
        line.append("NaN")
    clusters[line[0]] = "\t".join(line[1:])
clusts.close()

inFile  = open(options.result_file, "r")
uitFile = open(options.output_file, "w")
for l in inFile:
    l = l.rstrip().split("\t")
    if "genes" in l[0]:
        l = [x.replace(options.working_dir, "") for x in l]    # remove path from sample names
        l.append("\t".join([clusters["gene"], "length\tsequence\tlength of hit\te-value\tname and function"]))
    elif l[0] in clusters:
            l.append(clusters[l[0]] + "\t".join(names[l[0]]))
    else:
        if l[0] in names:
            l.append("\t".join(names[l[0]]))
    uitFile.write("\t".join(l)+"\n")
inFile.close()
uitFile.close()
