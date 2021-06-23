#!/usr/bin/env python
from matplotlib import pyplot as plt
import numpy as np
import argparse
import pandas as pd
from os import walk


def main():
    args = parse_args()
    progeny_kmers = importfile(args.filtered_progeny_file)


    # Next step is use the progeny file to filter the parent file
    d = {}
    two_col_file = open(args.filtered_parent_file)
    for line in two_col_file:
        key, value = line.split()
        d[key] = value
        # The code below keeps the values for d while filtering the keys for d using progeny_kmers
        parents_filtered_more = dict((key, d[key]) for key in [k for k in progeny_kmers.keys() if k in d])
    print(len(progeny_kmers))
    print(len(parents_filtered_more))
    # print(parents_filtered_more.values())

    outfile = open(args.save_file, "w")
    for key, value in parents_filtered_more.items():
        outfile.write(key + "\t" + str(value))
        outfile.write("\n")
    outfile.close()


# method imports files
def importfile(file):
    d = {}
    two_col_file = open(file)
    for line in two_col_file:
        key, value = line.split()
        d[key] = value
    return d

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--parent-file", help="File with kmers that will be filtered and saved")
    parser.add_argument("-c", "--filtered-progeny-file", help="File with kmers that will be used to filter")
    parser.add_argument("-s", "--save-file", help="path and file name to save the filtered kmer file")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug
    return parser.parse_args()


main()