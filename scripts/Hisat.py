import sys
import os

args = sys.argv
os.system("cd ../")
num = len(args[2:])/3

for i in range(2,2+num):
    command = "hisat2 -x {} -1 {} -2 {} | samtools view -Sb -f 2 {}".format(args[1],args[i],args[i+num],args[i+2*num])
    print(command)
    #os.system(command)