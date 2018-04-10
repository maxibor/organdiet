#!/usr/bin/env python

import argparse
from pprint import pprint


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

    basename = get_basename(mysam)
    with open(mysam, "r") as sam:
        with open(basename + ".best.aligned.fa", "w") as fw:
            for line in sam:
                linestrip = line.rstrip()
                linesplit = linestrip.split("\t")
                for field in linesplit:
                    if ("XM:i:") in field:
                        mismatch = field.split(":")[2]  # number of mismatches
                        mismatch = int(mismatch)
                        seq = linesplit[9]

                        seqlen = len(seq)  # length of aligned read
                        identity = (seqlen - mismatch) / seqlen
                        dbmatch = linesplit[2]  # match in database, NCBI id
                        readname = linesplit[0]
                        if identity >= idpercent and seqlen > minlen:
                            if readname not in readdict.keys():
                                readdict[readname] = [seq]
                                fw.write(">" + readname + "\n" + seq + "\n")

    pprint(readdict)
