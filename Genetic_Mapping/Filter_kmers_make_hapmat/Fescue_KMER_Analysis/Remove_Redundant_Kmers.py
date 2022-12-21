#!/usr/bin/env python

# Author email (darrianrtalamantes6@gmail.com)
# Author: Darrian Talamantes
# Objective: This will go through a hapmap file and remove Kmers/SNPs that are redundant in how they appear within
# all samples, keeping only 1 of them. These kmers/SNPs do not add any power when making a genetic map

import numpy as np
import pandas as pd
from os import walk
import argparse


def main():
    args = parse_args()

    # Run methods
    nonRedundantHapfile = removeRedundance(args.hapmatFile)
    saveHap(nonRedundantHapfile, args.hapmatFile, args.save)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-hap", "--hapmatFile", type=str, help="A hapmap file")
    parser.add_argument("-s", "--save", help="path and file name to save output array")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug
    return parser.parse_args()


def removeRedundance(hapmatFile):
    hapfile_OG = np.loadtxt(hapmatFile, dtype=str, skiprows=1, delimiter="\t")
    hapfile1 = np.delete(hapfile_OG, slice(0, 11), 1)
    hapfile1 = ["".join(item) for item in hapfile1.astype(str)]
    # # these next lines create dummy data.
    # hapfile1 = ["the", "of", "the", "peep", "eat", "the", "leap", "the", "feet", "feet", "feet", "read", "tell", "the"]
    # hapfileArray = np.array([[0],[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12],[13]])
    dupes = (list_duplicates_of(hapfile1))
    print(len(dupes))
    # #This is the time limiting step
    # # It may be helpfull to turn this into a dictionary with keys being the row number and
    # # searching those keys to pop em out of the list
    for i in reversed(dupes):
        hapfile_OG = np.delete(hapfile_OG, i, 0)
    return hapfile_OG


def list_duplicates_of(hapfile):
    # This creates a list of the positions of duplicates without listing the original in that list.
    oc_set = set()
    res = []
    for idx, kmer in enumerate(hapfile):
        if kmer not in oc_set:
            oc_set.add(kmer)
        else:
            res.append(idx)
    return res


def saveHap(finalHapFile, hapmapLoc, save):
    with open(hapmapLoc) as f:
        first_line = f.readline()
    np.savetxt(save, finalHapFile, delimiter="\t", fmt='%s')
    with open(save, 'r+') as f1:
        content = f1.read()
        f1.seek(0, 0)
        f1.write(first_line.rstrip('\r\n') + '\n' + content)


main()
