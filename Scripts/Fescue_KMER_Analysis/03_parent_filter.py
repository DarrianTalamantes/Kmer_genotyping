#!/usr/bin/env python
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from os import walk


def main():


    # Then keep ones that only appear once in the new parent file
    parents_kmers = importfile("/scratch/drt83172/Wallace_lab/TallFescue/Data/Practice/InterFiles/Kmers_Parent_filtered_step1.txt")
    parents_filtered = filter_for_rare(parents_kmers)
    print(len(parents_filtered))
    f = open("/scratch/drt83172/Wallace_lab/TallFescue/Data/Practice/InterFiles/Kmers_Parent_filtered_rare.txt", "w")
    f.write(str(parents_filtered))
    f.close()

# method imports files
def importfile(file):
    array = np.loadtxt(file, dtype=str)
    return array


def filter_for_rare(dic_file):
    d = dict((k, v) for k, v in dic_file.items() if v <= 1)
    return d


main()