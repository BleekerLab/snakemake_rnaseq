from optparse import OptionParser 
import glob
import os

# get the command line arguments
parser = OptionParser(description="Script that combines the normalized and differential expression data as outputted by DESeq2.R, the clusters as outputted by plotscript.R and the results from the blast against the reference proteome. Resulting in a combined tab-delimeted data file")
parser.add_option('-r', '--result_file', 
                    type=str,
                    default="results/result.csv",
                    metavar="",
                    help = "path and name of of the input file, being the output file of the R-script DESeq2.R, default = results/result.csv")
parser.add_option('-a', '--anno_file', 
                    type=str,
                    default="annotations/Petunia_axillaris_v1.6.2_proteins.mapman_annotation.txt",
                    metavar="",
                    help = "path and name of the stringtie transcriptome fasta file, a file containing the fastas of the de novo generated stringtie transcriptome, default = result/helperFile.csv")
parser.add_option('-c', '--cluster_file', 
                    type=str,
                    default="results/clusters.txt",
                    metavar="",
                    help = "path and name of the clusters file, a tab delimited file containing cluster numbers and depending on the method of clustering correlation values or membership values., default = result/helperFile.csv")
parser.add_option('-o', '--output_file', 
                    type=str,
                    default="results/final.txt",
                    metavar="",
                    help = "path and name of the output file. A tab delimimited file containing normalised reads, differential expressions and (adjusted)p-values, number of clusters the genes partisioned to and the hypothetical function as found by a blast. default = result/final.txt")
parser.add_option('-p', '--working_dir', 
                    type=str,
                    default="temp/mapped/",
                    metavar="",
                    help = "path to bam files, to be removed from the sample names in the header of the final output file. default = temp/mapped/")


(options, args) = parser.parse_args()

clusts = open(options.cluster_file, "r")
clusters = {}
for line in clusts:
    line = line.rstrip()
    line = line.split("\t")
    if len(line) == 2:
        line.append("NaN")
    clusters[line[0].lower()] = "\t".join(line[1:])

clusts.close()

# if there is an annotation file available add annotaions.
if options.anno_file in glob.glob(options.anno_file):
    if glob.glob(options.anno_file)[0].endswith(".gz"):
        os.system(f"gunzip {options.anno_file}")
        annoFile = ".".join(options.anno_file.split(".")[:-1])
    else:
        annoFile = options.anno_file
    fa = open(annoFile, "r")
    annos = {}
    count = 0
    for l in fa:
        l = l.split("\t")
        l = [x.strip("'") for x in l]
        if len(l)>2:
            if len(l[2])>1:
                if l[2].lower() not in annos:
                    annos[l[2].lower()] = ".".join(l[1].split(".")[:3]).replace(" ", "_")
                else:
                    annos[l[2].lower()] += (";" + ".".join(l[1].split(".")[:3]).replace(" ", "_"))
    fa.close()
    
    inFile  = open(options.result_file, "r")
    uitFile = open(options.output_file, "w")
    for l in inFile:
        l = l.rstrip().split("\t")
        if "genes" in l[0]:
            l = [x.replace(options.working_dir, "") for x in l]    # remove path from sample names
            l.append("\t".join([clusters["gene"], "annotation"]))
        elif l[0].lower() in clusters:
            if l[0].lower() in annos:
                l.append(clusters[l[0].lower()] + "\t" + annos[l[0].lower()])
            else:
                l.append(clusters[l[0].lower()] + "\tno annotation available.")
        elif l[0].lower() in annos:
            l.append("NaN\tNaN\t" + annos[l[0].lower()])
        else:
            l.append("NaN\tNaN\tno annotation available.")
        uitFile.write("\t".join(l)+"\n")
    inFile.close()
    uitFile.close()

else:
    inFile  = open(options.result_file, "r")
    uitFile = open(options.output_file, "w")
    for l in inFile:
        l = l.rstrip().split("\t")
        if "genes" in l[0]:
            l = [x.replace(options.working_dir, "") for x in l]    # remove path from sample names
            l.append("\t".join([clusters["gene"]]))
        elif l[0].lower() in clusters:
            l.append(clusters[l[0].lower()])
        else:
            l.append("NaN\tNaN")
        uitFile.write("\t".join(l)+"\n")
    inFile.close()
    uitFile.close()