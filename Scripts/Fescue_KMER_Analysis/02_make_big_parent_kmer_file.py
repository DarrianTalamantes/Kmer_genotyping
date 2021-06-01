#!/usr/bin/env python
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from os import walk



def main():

    filepath2 = "/scratch/drt83172/Wallace_lab/TallFescue/Data/Practice/Parent_KMERS/"
    _, _, filenames2 = next(walk(filepath2))
    parent_kmers = prescencecount(filenames2, 1, filepath2)
    f = open("/scratch/drt83172/Wallace_lab/TallFescue/Data/Practice/InterFiles/Kmers_Parent_filtered_step1.txt", "w")
    f.write(str(parent_kmers))
    f.close()

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


main()

