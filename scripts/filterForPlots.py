import numpy as np 
from optparse import OptionParser 

# get the command line arguments
parser = OptionParser(description="python script to create pvalue and log2(foldchange) filtered ")
parser.add_option('-i', '--input_file', 
                    type=str,
                    default="results/plotSelection.txt",
                    metavar="",
                    help = "path and name of of the input file, being the output file of the R-script DESeq2.R, default = results/plotSelection.txt")
parser.add_option('-f', '--helper_file', 
                    type=str,
                    default="result/helperFile.csv",
                    metavar="",
                    help = "path and name of the helper file, a tab delimited file containing one column samples and a column conditions, default = result/helperFile.csv")
parser.add_option('-o', '--output_file', 
                    type=str,
                    default="result/plotSelection.txt",
                    metavar="",
                    help = "path and name of the output file a tab delimimited file containing normalised reads of significantly differentially expressed genes of all the samples or averaged to the conditions, default = result/plotSelection.txt")
parser.add_option('-v', '--minimum_foldchange', 
                    type=float,
                    default=2,
                    metavar="",
                    help = "minimum log2(fold_change) in any combination of conditions; integer > 1, default = 2")
parser.add_option('-r', '--minimum_reads', 
                    type=int,
                    default=100,
                    metavar="", 
                    help = "minimum number of reads of all samples together; integer >= 0, default = 100")
parser.add_option('-p', '--maximum_pvalue', 
                    type=float,
                    default=0.05,
                    metavar="",
                    help = "maximum exepted adjusted pvalue; 0.0 to 1.0, default = 0.05")
parser.add_option('-a', '--average_samples', 
                    type=str,
                    default="yes",
                    metavar="",
                    help = "output needs to contain averages of conditions; yes or no, default = yes")

(options, args) = parser.parse_args()


# get a list of the conditions 
helper = open(options.helper_file, "r")
samples = {}
conditions = []
x=0
for l in helper:
    l = l.rstrip()
    if len(l) > 0:
        x += 1
        l = l.split(",")
        name = l[1]
        samples[x] = name
        if name not in conditions:
            conditions.append(name)
helper.close()

# function to average the samples within a condition
def getAverages(counts):
    global samples, conditions
    sets = {}
    averages = [counts[0]]
    for i in samples:
        if samples[i] in sets:
            sets[samples[i]].append(float(counts[i]))
        else:
            sets[samples[i]] = [float(counts[i])]
    for c in conditions:
        averages.append(str(np.mean(sets[c])))
    return(averages)

inputData = open(options.input_file, "r")
outAll    = open(options.output_file, "w")
totalBoth = 0
count     = 0
min_reads = options.minimum_reads
foldchange= options.minimum_foldchange
pvalue    = options.maximum_pvalue
avarage   = options.average_samples.upper()

for l in inputData:
    count += 1
    if l.startswith("genes") == False and len(l) > 5:
        l = l.replace("NA", "1")         # replace NA values in P-value to 1
        l = l.replace("#VALUE!", "0")    # replace NA values in foldchange to 0
        f = l.rstrip()                   # remove white space and line break at the end.
        f = f.split("\t")                # split string to list on tab
        if len(f) > x+1:
            fc = [float(f[i]) for i in range(x+2, len(f), 6)]       # list of faultchanges
            pv = [float(f[i]) for i in range(x+6, len(f), 6)]       # list of adjusted p-values
            counts = [float(f[i]) for i in range(1,x+1)]
            if sum(counts) > min_reads:             # total number of reads more then ~20 per sample
                for a1, a2 in zip(fc,pv):
                    if a1 > foldchange and a2 < pvalue and a2 == min(pv):
                        if avarage == "YES":
                            line = "\t".join(getAverages(f[:x+1]))
                        else:
                            line = "\t".join(f[:x+1])
                        outAll.write(line+"\n")
                        totalBoth += 1
                        break
                    elif a1 < -1*foldchange and a2 < pvalue and a2 == min(pv):
                        if avarage == "YES":
                            line = "\t".join(getAverages(f[:x+1]))
                        else:
                            line = "\t".join(f[:x+1])
                        outAll.write(line+"\n")
                        totalBoth += 1
                        break
        else:
            print("No fold changes and p-values were found. Try rerunning the DESeq2.R with a higher maximum fraction (to be adjusted in config.yaml)")
    else:
        if len(l) > 10:
            if avarage == "YES":
                line = "\t".join(conditions)
            else:
                line = "\t".join(l.split("\t")[1:x+1])
            outAll.write(l.split("\t")[0] + "\t" + line + "\n")

print(f"Total number of genes is: {str(count-1)}.")
print(f"Number of genes kept for plots is: {str(totalBoth)}.")

inputData.close() 
outAll.close()
