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
    # variable inputs
    args = parse_args()
    parentD = args.parent_directory
    progenyD = args.progeny_directory
    cross = args.cross
    masterKmersL = args.master_kmer_file
    save = args.save_file

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
    # Fill in the first 11 columns of hapmat file using pandas data frame
    hapmat = pd.DataFrame(data=hapmatArray, index=None,
                          columns=['rs#', 'alleles', 'chrom', 'pos', 'strand', 'assembly#',
                                   'center', 'protLSID', 'assayLSIP', 'panelLSID', 'QCcode'],
                          dtype=str, copy=None)
    for x in range(len(parentfiles)):
        hapmat[parentfiles[x]] = ""
    for x in range(len(progenyfiles)):
        hapmat[progenyfiles[x]] = ""
    colnames = list(hapmat.columns)
    hapmat = hapmat.to_numpy(dtype='U113')
    precenceAbsance = np.zeros((len(masterkeys), len(progenyfiles) + len(parentfiles)), dtype=str)
    ACdata = append2hapmat(masterKmers, parentfiles, progenyfiles, parentD, progenyD, precenceAbsance)
    hapCol = 11
    FinalData = np.zeros(shape=(np.size(masterKmers) + 1, len(hapmat[0])), dtype='U113')
    # Puts in column names
    for x in range(len(colnames)):
        FinalData[0][x] = colnames[x]
    # Puts in metadata
    for x in range(11):
        for y in range(len(hapmat)):
            FinalData[y+1][x] = hapmat[y][x]
    # Puts in A's and C's
    for x in range(len(ACdata[0])):
        for y in range(len(ACdata)):
            FinalData[y + 1][hapCol] = ACdata[y][x]
        hapCol += 1
    np.savetxt(save, FinalData, delimiter='\t', fmt='%s')

    # This method will use the master kmer file to create a numpy array showing prescence abscence of kmers


def append2hapmat(masterKmers, listOfiles1, listOfiles2, directory1, directory2, precenceAbsance):
    fileNum = 0
    for file1 in listOfiles1:
        tempfile = importfile(directory1 + "/" + file1)
        for x in range(len(masterKmers)):
            kmer1 = masterKmers[x]
            if kmer1 in tempfile:
                precenceAbsance[x][fileNum] = 'A'
            else:
                precenceAbsance[x][fileNum] = 'C'
        fileNum += 1
    for file2 in listOfiles2:
        tempfile = importfile(directory2 + "/" + file2)
        for x in range(len(masterKmers)):
            kmer1 = masterKmers[x]
            if kmer1 in tempfile:
                precenceAbsance[x][fileNum] = 'A'
            else:
                precenceAbsance[x][fileNum] = 'C'
        fileNum += 1
        f = open("counter.txt", "w")
        f.write(str(fileNum))
        f.close()
    return precenceAbsance

#Imports a file into a dictionary
def importfile(file):
    # method imports files
    d = {}
    two_col_file = open(file)
    for line in two_col_file:
        key, value = line.split()
        d[key] = value
    return d

def parse_args():
    # Inputs, master kmer file, directory with parental kmers, directory with progeny kmers, corss im working on,
    # and save file path and name
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--master-kmer-file", help="File with all kmers being analyzed")
    parser.add_argument("-pd", "--parent-directory", help="directory with progeny kmers")
    parser.add_argument("-cd", "--progeny-directory", help="directory with progeny kmers")
    parser.add_argument("-s", "--save-file", help="path and file name to save hapmat file")
    parser.add_argument("-c", "--cross", help="current  biparental cross")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug
    return parser.parse_args()




main()
