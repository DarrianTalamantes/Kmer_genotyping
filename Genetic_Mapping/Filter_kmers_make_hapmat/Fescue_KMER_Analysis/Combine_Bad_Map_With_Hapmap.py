#!/usr/bin/env python
# from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from os import walk
import argparse

# This code will combine the first map of a biparental cross I made with the hapmat file to
# add location data to the hapmat file. Beacuse I cant make map in sapelo I havent implemented this into sapelo.

def main():
    # args = parse_args()
    # changedis(args.hapmap, args.genetic_map,args.save_file)
    changedis("/home/drt83172/Documents/Tall_fescue/Kmer_Genotyping/Kmer_genotyping/Genetic_Mapping/Data/Genotype_Files/301x304_hapmapNonRedun.txt",
              "/home/drt83172/Documents/Tall_fescue/Kmer_Genotyping/Kmer_genotyping/Genetic_Mapping/Data/Maps/301x304_genotype_NonRedun_map.txt",
              "/home/drt83172/Documents/Tall_fescue/Kmer_Genotyping/Kmer_genotyping/Genetic_Mapping/Data/Maps/301_304_New_hapmap_NonRedun.txt")


# Will open genotype file and hapmap file. It will then compare the kmers and if they match replaces hapmap
# distance with genetic map distance
def changedis(hapmap, gp, save):
    hapfile = np.loadtxt(hapmap, dtype=str, skiprows=1)
    genmap = np.loadtxt(open(gp,'rt').readlines()[:-3], dtype=str, delimiter='\t', skiprows=10)
    # print(genmap[20][1]) # Iterating through [][1] is the genetic distances
    # print(hapfile[3][3]) # Iterating through [][3] is genetic distance
    # print(genmap.shape[0])
    # print(hapfile.shape[0])
    # I  round the positions to the nearest int here. May be a better way to keep order (Probs dont matter though)
    for gPos in range(genmap.shape[0]):
        kmer = genmap[gPos][0]
        for hPos in range(genmap.shape[0]):
            if kmer == hapfile[hPos][0]:
                dist = float(genmap[gPos][1])
                hapfile[hPos][3] = round(dist)
                print(kmer + "   " + hapfile[hPos][0])
                print(hapfile[hPos][3])
    # Grabs hapmap file first line, Saves numpy array to txt. prepends first line to saved file
    with open(hapmap) as f:
        first_line = f.readline()
    np.savetxt(save,hapfile, delimiter="\t", fmt='%s')
    with open(save, 'r+') as f1:
        content = f1.read()
        f1.seek(0, 0)
        f1.write(first_line.rstrip('\r\n') + '\n' + content)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-h", "--hapmap", help="The hapmap file")
    parser.add_argument("-m", "--genetic-map", type=int, help="the genetic map file")
    parser.add_argument("-s", "--save-file", help="path and file name to save the new hapmap file")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug
    return parser.parse_args()




main()








