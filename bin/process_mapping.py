#!/usr/bin/env python3

import argparse
import requests
from pprint import pprint
import pickle

def get_args():
    '''This function parses and return arguments passed in'''
    parser = argparse.ArgumentParser(
    prog='Alignment to specie',
    description='From a sam file, returns specie names')
    parser.add_argument('mapsam', help="path to sam file")
    parser.add_argument(
    '-outdir',
    default="./",
    help="path to (existing) output directory")

    args = parser.parse_args()

    mysam = args.mapsam
    myoutdir = args.outdir

    return(mysam, myoutdir)

def request_to_specie(ncbi_id):
    """
    Takes a ncbi_id as input (ex: NC_012978.1), makes a call to JGI
    taxonomy API, and returns specie name.

    INPUT:
        ncbi_id(string) ex: NC_012978.1
    OUPUT:
        specie(str) specie name
    """

    request = "http://taxonomy.jgi-psf.org/sc/accession/"+ncbi_id
    response = requests.get(request)
    answer = response.text
    specie = answer.split(":")[-1]
    return(specie)

def get_basename(samfile_name):
    if ("/") in samfile_name:
        basename = samfile_name.split("/")[-1].split(".")[0]
    else:
        basename = samfile_name.split(".")[0]
    return(basename)


if __name__ == "__main__":
    mysam, myoutdir = get_args()

    matchdict = {}
    readdict = {}

    fileObject = open("/home/maxime/Documents/ulug_depe/scripts/organdiet/data/db/specie_taxonomy.pickle",'rb')
    specie_taxonomy = pickle.load(fileObject)
    id_threshold = 0.99
    basename = get_basename(mysam)
    with open(mysam, "r") as sam:
        with open(basename+"_species.csv","w") as fw:
            fw.write("specie, ncbi_id, identity_percentage, match_length, sequence\n")
            for line in sam:
                linestrip = line.rstrip()
                linesplit = linestrip.split("\t")
                for field in linesplit:
                    if ("XM:i:") in field:
                        mismatch = field.split(":")[2] # number of mismatches
                        mismatch = int(mismatch)
                        seq = linesplit[9]

                        seqlen = len(seq) # length of aligned read
                        identity = (seqlen-mismatch)/seqlen
                        dbmatch = linesplit[2] #match in database
                        readname = linesplit[0]
                        if identity >= id_threshold :
                            print(specie_taxonomy[dbmatch], identity, seqlen)
                            fw.write(specie_taxonomy[dbmatch]+", "+str(dbmatch)+", "+str(identity)+", "+str(seqlen)+", "+seq+"\n")

                    # if dbmatch not in matchdict.keys():
                    #     matchdict[dbmatch] = [identity]
                    # else :
                    #     matchdict[dbmatch].append(identity)
                    #
                    #     if readname not in readdict.keys():
                    #         readdict[readname] = [(dbmatch, specie_taxonomy[dbmatch], identity)]
                    #     else:
                    #         readdict[readname].append((dbmatch, identity))
    # pprint(readdict)
