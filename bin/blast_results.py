#!/usr/bin/env python3

import sys
import argparse

def get_args():
    '''This function parses and return arguments passed in'''
    parser = argparse.ArgumentParser(
    prog='Blast search from csv results',
    description='From a csv file, perform blast search to return LCA')
    parser.add_argument('infile', help="path to csv file")
    parser.add_argument(
    '-threads',
    default=4,
    help="Number of threads to use for blast")

    args = parser.parse_args()

    infile = args.infile
    threads = args.threads

    return(infile, threads)

def get_basename(samfile_name):
    if ("/") in samfile_name:
        basename = samfile_name.split("/")[-1].split(".")[0]
    else:
        basename = samfile_name.split(".")[0]
    return(basename)



def common_ancestor(ncbi_id):
    """
    Takes  ncbi_ids as input (ex: ['NC_012978.1','NC_0125678.1']), makes a call to JGI
    taxonomy API, and returns common_ancestor.

    INPUT:
        ncbi_id(list) ex: ['NC_012978.1','NC_0125678.1']
    OUPUT:
        common_ancestor(str) common_ancestor name
    """

    ncbi_id = ",".join(ncbi_id)
    request = "http://taxonomy.jgi-psf.org/sc/name/ancestor/c"+ncbi_id
    response = requests.get(request)
    answer = response.text
    specie = answer.split(":")[-1]
    return(common_ancestor)

if __name__ == "__main__":

    infile, threads = get_args()

    with open(infile, "r") as f:
        next(f)
        for line in f:
            basename = get_basename(infile)
            splitline = line.split(",")
            specie = splitline[0]
            ncbi_id = splitline[1]
            ncbi_id = ncbi_id.replace(" ","")
            sequence = splitline[4]
            with open(basename+"_"+ncbi_id+".fa","w") as fasta:
                fasta.write(">"+ncbi_id+"|"+specie+"\n")
                fasta.write(sequence+"\n")

            cmd = "blastn"

            break
