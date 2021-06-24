#!/usr/bin/env python
from matplotlib import pyplot as plt
import numpy as np
import argparse
import pandas as pd
from os import walk


def main():
    args = parse_args()
    kmer_file = importfile(args.kmer_file)

    usefull_kmer_list = importfile(args.kmer_list)

    filtered_kmers = dict()
    # creates a set of the keys that are found in both these dictionaries
    kmer_file.intersection(usefull_kmer_list)
    for key in usefull_kmer_list:
        if key in kmer_file.keys():
            d = {key : kmer_file[key]}
            filtered_kmers.update(d)
        else:
            d = {key : kmer_file[key]}
            filtered_kmers.update(d)

    # The code below keeps the values for d while filtering the keys for d using progeny_kmers

    print(len(kmer_file))
    print(len(filtered_kmers))
    # print(parents_filtered_more.values())

    outfile = open(args.save_file, "w")
    for key, value in filtered_kmers.items():
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
    parser.add_argument("-i", "--kmer-file", help="File with kmers that will be filtered and saved")
    parser.add_argument("-l", "--kmer-list", help="Kmers to be used to filter the other file with")
    parser.add_argument("-s", "--save-file", help="path and file name to save the filtered kmer file")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug
    return parser.parse_args()


main()
