#!/usr/bin/env python

# Author email (darrianrtalamantes6@gmail.com)
# Author: Darrian Talamantes
# Objective: Use a file that contains all kmers after filtration to go through every individual kmer file and create a
# hapmat file.

import argparse
import numpy as np
import os as os
import pandas as pd


# Inputs, master kmer file, directory with parental kmers, directory with progeny kmers, corss im working on,
def main():
    # variable setting
    parentD = "/home/drt/Desktop/UGA/Wallace_Lab/Tall_Fescue_Grass/Kmer_genotyping/Genetic_Mapping/Data/example_Kmers" \
              "/parents/"
    progenyD = "/home/drt/Desktop/UGA/Wallace_Lab/Tall_Fescue_Grass/Kmer_genotyping/Genetic_Mapping/Data" \
               "/example_Kmers/Progeny/"
    cross = "314x312"
    masterKmersL = "/home/drt/Desktop/UGA/Wallace_Lab/Tall_Fescue_Grass/Kmer_genotyping/Genetic_Mapping/Data" \
                   "/example_Kmers/masterExample.txt"

    # this creates a list of files within the specified directory
    _, _, progenyfiles = next(os.walk(progenyD))
    print(progenyfiles)
    _, _, parentfiles = next(os.walk(parentD))
    print(parentfiles)

    # Create data for hapmat file
    masterKmers = importfile(masterKmersL)
    masterkeys = masterKmers.keys()
    masterkeys = list(masterkeys)
    masterKmers = np.array(masterkeys)
    hapmatArray = np.zeros(shape=(np.size(masterKmers), 11), dtype='U36')
    for x in range(np.size(masterKmers)):
        hapmatArray[x, 0] = masterKmers[x]
        hapmatArray[x, 1] = "A/C"
        hapmatArray[x, 2] = cross
        hapmatArray[x, 3] = x * 356  # Just making random uniq numbers
        hapmatArray[x, 4] = "+"
        hapmatArray[x, 5] = "NA"
        hapmatArray[x, 6] = "NA"
        hapmatArray[x, 7] = "NA"
        hapmatArray[x, 8] = "NA"
        hapmatArray[x, 9] = "NA"
        hapmatArray[x, 10] = "NA"

    # Fill in the first 11 columns of hapmat file using pandas data frame
    hapmat = pd.DataFrame(data=hapmatArray, index=None,
                          columns=['rs#', 'alleles', 'chrom', 'pos', 'strand', 'assembly#',
                                   'center', 'protLSID', 'assayLSIP', 'panelLSID', 'QCcode'],
                          dtype=str, copy=None)
    print(hapmat.shape)
    hapmat = append2hapmat(masterKmers, parentfiles, parentD, hapmat)
    hapmat = append2hapmat(masterKmers, progenyfiles, progenyD, hapmat)
    print(hapmat)


def append2hapmat(masterKmers, listOfiles, directoryOfiles, file2append):
    # This method iterates through the fed list of files and checks if the kmer from that file is in the master list
    # It creates a tiny array full of A's or C's that will then be added to the hapmat file
    for x in range(len(listOfiles)):
        tinymat = np.zeros(shape=(np.size(masterKmers)), dtype="U36")
        file = listOfiles[x]
        # print("reading file" + file)
        tempfile = importfile(directoryOfiles + "/" + file)
        for y in range(len(masterKmers)):
            if masterKmers[y] in tempfile:
                tinymat[y] = 'A'
                # print("Present")
            else:
                tinymat[y] = 'C'
                # print("Missing")
            file2append[file] = tinymat
    return file2append


def importfile(file):
    # method imports files
    d = {}
    two_col_file = open(file)
    for line in two_col_file:
        key, value = line.split()
        d[key] = value
    return d


main()
