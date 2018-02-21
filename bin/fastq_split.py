import argparse
import gzip
from collections import OrderedDict

def get_args():
    '''This function parses and return arguments passed in'''
    parser = argparse.ArgumentParser(
    prog='fastq-split',
    description=
    """
    Split a Paired-end {basename}.fastq file into {basename}.R1.fastq {basename}.R2.fastq
    Also works on gzipped fastq {basename}.fastq.gz
    """
    )
    parser.add_argument('infile', help="path to PE fastq")

    args = parser.parse_args()

    infile = args.infile


    return(infile)

def get_basename(file_name):
    if ("/") in file_name:
        basename = file_name.split("/")[-1].split(".")[0]
    else:
        basename = file_name.split(".")[0]
    return(basename)


def split_fastq (infile, basename):

    R1dict = OrderedDict()
    R2dict = OrderedDict()

    print("Splitting "+basename+ ".fastq into "+basename+".R1.fastq and "+basename+".R2.fastq")
    with open(infile, "r") as f:
        myflag = True
        for line in f:
            line = line.rstrip()
            if myflag == True:
                instrument = line.split()[0].split(":")[0]
                myflag = False

            if line.startswith(instrument):
                seqname = line
                try :
                    read = int(line.split()[1].split(":")[0])
                except IndexError as e:
                    print(line)
                    print(e)
                going = True
                continue
            elif line[0] !="@" and line[0] !="+" and read == 1 and going == True:
                R1dict[seqname] = [line]
                going = False
                continue
            elif line[0] !="@" and line[0] !="+" and read == 2 and going == True:
                R2dict[seqname] = [line]
                going = False
                continue
            else:
                if read == 1:
                    R1dict[seqname].append(line)
                elif read == 2:
                    R2dict[seqname].append(line)

    print("Writing "+basename+".R1.fastq")
    with open(basename+".R1.fastq", "w") as fq1:
        for akey in R1dict.keys():
            towrite1 = akey+"\n"+"\n".join(R1dict[akey])+"\n"
            fq1.write(towrite1)

    print("Writing "+basename+".R2.fastq")
    with gzip.open(basename+".R2.fastq", "w") as fq2:
        for akey in R2dict.keys():
            towrite2 = (akey+"\n"+"\n".join(R2dict[akey])+"\n")
            fq2.write(towrite2)

def split_fastq_gz(infile, basename):
    R1dict = OrderedDict()
    R2dict = OrderedDict()

    print("Splitting "+basename+ ".fastq.gz into "+basename+".R1.fastq.gz and "+basename+".R2.fastq.gz")
    with gzip.open(infile, "rb") as f:
        myflag = True
        for line in f:
            line = line.decode('utf-8')
            line = line.rstrip()
            if myflag == True:
                instrument = line.split()[0].split(":")[0]
                myflag = False

            if line.startswith(instrument):
                seqname = line
                try :
                    read = int(line.split()[1].split(":")[0])
                except IndexError as e:
                    print(line)
                    print(e)
                going = True
                continue
            elif line[0] !="@" and line[0] !="+" and read == 1 and going == True:
                R1dict[seqname] = [line]
                going = False
                continue
            elif line[0] !="@" and line[0] !="+" and read == 2 and going == True:
                R2dict[seqname] = [line]
                going = False
                continue
            else:
                if read == 1:
                    R1dict[seqname].append(line)
                elif read == 2:
                    R2dict[seqname].append(line)

    print("Writing "+basename+".R1.fastq.gz")
    with gzip.open(basename+".R1.fastq.gz", "wb") as fq1:
        for akey in R1dict.keys():
            towrite1 = akey+"\n"+"\n".join(R1dict[akey])+"\n"
            fq1.write(towrite1.encode('utf-8'))

    print("Writing "+basename+".R2.fastq.gz")
    with gzip.open(basename+".R2.fastq.gz", "wb") as fq2:
        for akey in R2dict.keys():
            towrite2 = (akey+"\n"+"\n".join(R2dict[akey])+"\n")
            fq2.write(towrite2.encode('utf-8'))




if __name__ == "__main__":
    infile = get_args()
    basename = get_basename(infile)
    if infile.endswith(".gz"):
        split_fastq_gz(infile, basename)
    else :
        split_fastq(infile, basename)
