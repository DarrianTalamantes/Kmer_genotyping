#!/usr/bin/env python
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from os import walk
import argparse


def main():
    args = parse_args()
    # allows for outside argument to decide filepath
    filepath1 = args.kmer_count_directory
    _, _, filenames1 = next(walk(filepath1))


    # Make the upper and lower cutoffs options in the method below. They should be ints.

    if args.upper_bound is not None:
        prescencecount(filenames1, filepath1, args.upper_bound, False, args.save_file)
    elif args.lower_bound is not None:
        prescencecount(filenames1, filepath1, args.lower_bound, True, args.save_file)
    # f = open(args.save_file, "w")
    # f.write(str(progeny_kmers))
    # f.close()


# This method adds arguments into the program
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", "--kmer-count-directory", help="A directory only filled with kmer count files")
    parser.add_argument("-u", "--upper-bound", type=int, help="insert an int to set the upper bound")
    parser.add_argument("-l", "--lower-bound", type=int, help="insert an int to set the lower bound")
    parser.add_argument("-s", "--save-file", help="path and file name to save the filtered kmer file")
    parser.add_argument("--debug", default=False, action="store_true")
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
def prescencecount(filenames,  filepath, cutoff, method, savefile):
    # The for loop below counts how many files a certain kmer appears in
    freq = {}
    for file in filenames:
        print(file)
        kmer_counts = importfile(filepath + "/" + file)
        # load the file one line at a time here.
        # kmer_list = kmer_counts[:, 0] old code, delete if new works
        for kmer in kmer_counts[:, 0]:
            if kmer in freq:
                freq[kmer] += 1
            else:
                freq[kmer] = 1
    # This is creating a second dictionary that will cut off any kmer found less than x times
    # This part can be written to a file
    outfile = open(savefile, "w")
    if method:
        print("Went lower bound")
        for key, value in freq.items():
            if value >= cutoff:
                outfile.write(key + "\t" + str(value))
                outfile.write("\n")
                # freq = dict((k, v) for k, v in freq.items() if v >= cutoff) old code may be usefull
    if not method:
        print("Went upper bound")
        for key, value in freq.items():
            if value <= cutoff:
                outfile.write(key + "\t" + str(value))
                outfile.write("\n")
                # freq = dict((k, v) for k, v in freq.items() if v <= cutoff)
    outfile.close()
    # print(d.items())


main()

