#!/usr/bin/env python
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from os import walk
import argparse


def main():
    # args.seq_hits
    args = parse_args()
    filepath1 = args.kmer_count_directory
    _, _, filenames1 = next(walk(filepath1))


    # Make the upper and lower cutoffs options in the method below. They should be ints.     

    progeny_kmers = prescencecount(filenames1,filepath1, 2)
    f = open(args.kmer.save_file, "w")
    f.write(str(progeny_kmers))
    f.close()


# This method adds arguments into the program
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", "--kmer-count-directory", help="A directory only filled with kmer count files")
    parser.add_argument("-u", "--upper-bound", help="insert an int to set the upper bound")
    parser.add_argument("-l", "--lower-bound", help="insert an int to set the lower bound")
    parser.add_argument("-s", "--save-file", help="path and file name to save the filtered kmer file")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug
    return parser.parse_args()


# method imports files
def importfile(file):
    array = np.loadtxt(file, dtype=str)
    # Process this text file line by line.
    return array


# method counts how many files a kmer appers in in a directory.
def prescencecount(filenames,  filepath, cutoff):
    # The for loop below counts how many files a certain kmer appears in
    freq = {}
    for file in filenames:
        print(file)
        kmer_counts = importfile(filepath + file)
        # load the file one line at a time here.
        kmer_list = kmer_counts[:, 0]
        for kmer in kmer_list:
            if kmer in freq:
                freq[kmer] += 1
            else:
                freq[kmer] = 1
    # This is creating a second dictionary that will cut off any kmer found less than x times
    # THis part can be written to a file
    d = dict((k, v) for k, v in freq.items() if v >= cutoff)
    # print(d.items())
    return d


main()

