#! python

# Goal of code is to insert uniq lines between every line to make a regular file into a fasta file

#packages
from contextlib import redirect_stdout
import re
import argparse




# Code start
def main():
    args = parse_args()
    fasta_creation(args.kmer,args.save)


def fasta_creation(kmer,save):    
    file = open(kmer, "r")
    kmer_file = file.readlines()
    x = 1

    with open(save, "w") as f:
        with redirect_stdout(f):
            for line in kmer_file:
                line2 = re.split(r'\t+', line)
                print(">" + line2[0] + "_" + str(x))
                print(line2[0])
                x += 1
    file.close()
    f.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-kmer", "--kmer", type=str, help="A kmer file")
    parser.add_argument("-s", "--save", help="path and file name to save output array")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug
    return parser.parse_args()
main()