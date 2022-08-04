#!/usr/bin/env python

# Author email (darrianrtalamantes6@gmail.com)
# Author: Darrian Talamantes
# Objective: Make a binomial distribution which will allow us to know how to filter the progeny files

import argparse
from os import walk
import numpy as np


# Inputs:
# The number of progeny in the Working_Kmers_Progeny folder "size"
# The number of trials "n" should be one
# The percent cahnce of success "p", we assume mendialian genetics so .5

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-tot", "--total-progeny", type=int, help="the total amount of progeny in a cross")
    parser.add_argument("-t", "--trials", type=int, help="insert the number of progeny you have here")
    parser.add_argument("-p", "--percent-success", type=float, help="this is your chance of success, here its always .5")
    parser.add_argument("-s", "--save-file", help="path and file name to save output array")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug
    return parser.parse_args()

def main():
    args = parse_args()

    a = np.random.binomial(n=args.total_progeny, p=args.percent_success, size=args.trials)
    std = np.std(a)
    x = np.mean(a)
    std2l = x - (2*std)
    std2u = (2*std) + x
    final_array = [round(std2l, 0), round(std2u, 0), x, std]
    np.savetxt(args.save_file, final_array, delimiter=',',fmt='%i')




main()
