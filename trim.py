import os

dataPath = "/home/maxime/Documents/ulug_depe/data/"
respath = "/home/maxime/Documents/ulug_depe/results/trimming/"
qcrespath = "/home/maxime/Documents/ulug_depe/results/trimqc/"

threads = 8

fqfiles = []
for afile in os.listdir(dataPath):
    if afile.endswith("R1.fastq.gz"):
        fqfiles.append(afile.split("_")[0])


for afile in fqfiles:
    basename = afile
    r1 = dataPath+basename+"_R1.fastq.gz"
    r2 = dataPath+basename+"_R2.fastq.gz"

    cmd = "AdapterRemoval --file1 "+r1+" --file2 "+r2+" --basename "+respath+basename+" --threads "+str(threads)+" --trimns --trimqualities --collapse"
    print("Trimming and cleaning for ", basename)
    print(cmd)
    os.system(cmd)

print("Quality Check with FASTQC")
os.system("fastqc "+respath+"*.collapsed --threads "+str(threads)+" --outdir "+qcrespath)

print("Generating report with MULTIQC")
os.system("multiqc "+qcrespath+ " --outdir "+qcrespath)
