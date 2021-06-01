#!/usr/bin/env python
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from os import walk


def main():

    parent_kmers = importfile("/scratch/drt83172/Wallace_lab/TallFescue/Data/Practice/InterFiles/Kmers_Parent_filtered_rare.txt")
    progeny_kmers = importfile("/scratch/drt83172/Wallace_lab/TallFescue/Data/Practice/InterFiles/Kmers_Progeny_filtered_step1.txt")


    # Next step is use the progeny file to filter the parent file or filter parent file
    parents_filtered_more = dict((key, parent_kmers[key]) for key in [k for k in progeny_kmers.keys() if k in parent_kmers])
    print(len(progeny_kmers))
    print(len(parents_filtered_more))
    # print(parents_filtered_more.values())

    np.savetxt("/scratch/drt83172/Wallace_lab/TallFescue/Data/Practice/Final_Kmers/Final_Kmers.txt", parents_filtered_more)
    f = open("/scratch/drt83172/Wallace_lab/TallFescue/Data/Practice/Final_Kmers/Final_Kmers.txt", "w")
    f.write(str(parents_filtered_more))
    f.close()

# method imports files
def importfile(file):
    d = {}
    two_col_file = open(file)
    for line in two_col_file:
        key, value = line.split()
    d[key] = value
    return d
main()