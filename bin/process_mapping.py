#!/usr/bin/env python

import argparse
from pprint import pprint
# import pickle
# import requests

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
    parser.add_argument(
    '-minlen',
    default=28,
    help="Minimum length of match to report")
    parser.add_argument(
    '-idpercent',
    default=0.99,
    help="Minimum identity percentage of match to report")

    args = parser.parse_args()

    mysam = args.mapsam
    myoutdir = args.outdir
    minlen = args.minlen
    idpercent = args.idpercent


    return(mysam, myoutdir, minlen, idpercent)

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

def get_basename(samfile_name):
    if ("/") in samfile_name:
        basename = samfile_name.split("/")[-1].split(".")[0]
    else:
        basename = samfile_name.split(".")[0]
    return(basename)


if __name__ == "__main__":
    mysam, myoutdir, minlen, idpercent = get_args()

    matchdict = {}
    readdict = {}

    # fileObject = open("/home/maxime/Documents/ulug_depe/scripts/organdiet/data/db/specie_taxonomy.pickle",'rb')
    specie_taxonomy = pickle.load(fileObject)
    basename = get_basename(mysam)
    with open(mysam, "r") as sam:
        with open(basename+".best.aligned.fa","w") as fw:
        # with open(basename+"_species.fa","w") as fw:
            # fw.write("specie, ncbi_id, identity_percentage, match_length, sequence\n")
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
                        dbmatch = linesplit[2] #match in database, NCBI id
                        readname = linesplit[0]
                        if identity >= idpercent and seqlen > minlen:
                            # print(specie_taxonomy[dbmatch], identity, seqlen)
                            # fw.write(specie_taxonomy[dbmatch]+", "+str(dbmatch)+", "+str(identity)+", "+str(seqlen)+", "+seq+"\n")

                            if readname not in readdict.keys():
                                readdict[readname] = [seq]
                                fw.write(">"+readname+"\n"+seq+"\n")

    pprint(readdict)
