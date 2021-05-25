from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import time
import errno
import os
from os import walk


def main():
    print("HELP!")
    # The code below this takes all the files in a directory and puts their names into a list
    filepath1 = "/home/drt83172/Documents/Tall_fescue/Kmer_analyses/Example/Data/KMER_Files/"
    _, _, filenames1 = next(walk(filepath1))
    filepath2 = "/home/drt83172/Documents/Tall_fescue/Kmer_analyses/Example/Data/KMER_Files_Pretend_parents/"
    _, _, filenames2 = next(walk(filepath2))

    # does prescence counts in both progeny and parental directories
    progeny_kmers = prescencecount(filenames1, 2, filepath1)
    parent_kmers = prescencecount(filenames2, 1, filepath2)

    # Next step is use the progeny file to filter the parent file or filter parent file
    parents_filtered = dict((key, parent_kmers[key]) for key in [k for k in progeny_kmers.keys() if k in parent_kmers])
    print(len(progeny_kmers))
    print(len(parents_filtered))
    # print(parents_filtered.values())

    # Then keep ones that only appear once in the new parent file
    parents_filtered_more = filter_for_rare(parents_filtered)
    print(len(parents_filtered_more))


# method imports files
def importfile(file):
    array = np.loadtxt(file, dtype=str)
    return array

# method counts how many files a kmer appers in in a directory.
def prescencecount(filenames, cutoff, filepath):
    # The for loop below counts how many files a certain kmer appears in
    freq = {}
    for file in filenames:
        print(file)
        kmer_counts = importfile(filepath + file)
        kmer_list = kmer_counts[:, 0]
        for kmer in kmer_list:
            if kmer in freq:
                freq[kmer] += 1
            else:
                freq[kmer] = 1
    # This is creating a second dictionary that will cut off any kmer found less than x times
    d = dict((k, v) for k, v in freq.items() if v >= cutoff)
    # print(d.items())
    return d


def filter_for_rare(dic_file):
    d = dict((k, v) for k, v in dic_file.items() if v <= 1)
    return d
main()
