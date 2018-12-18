import os
import sys

tc      = sys.argv[1]      #path/name of transcriptome
outFile = sys.argv[2]      #path/name of output file
bams    = " ".join(sys.argv[3:])  #list of path/bam files

command = "featureCounts -t exon -a {0} -o {1} {2}".format(tc,outFile,bams) 
print (command)
os.system(command)