import numpy as np
import argparse
import pandas as pd


def main():
    args = parse_args()
    parent_count = import_kmer_count_table(args.parent_counts)
    progeny_count = import_kmer_count_table(args.progeny_counts)
    scores = score_maker(parent_count, progeny_count)
    scores.to_csv(args.save_file)


def score_maker(parent_count, progeny_count):
    # Getting column and row names for later insertion
    parent_names = parent_count.index
    progeny_names = progeny_count.index
    # turning pandas into numpy
    parent_array = parent_count.to_numpy()
    progeny_array = progeny_count.to_numpy()
    # Making empty score array
    scoreCols = len(parent_array)
    scoreRows = len(progeny_array)
    score_array = np.zeros((scoreRows, scoreCols), dtype=int)
    print("rows", len(score_array), "columns", len(score_array[0]))
    # Doing the actual scoring
    for par in range(len(parent_array)):
        for prog in range(len(progeny_array)):
            for k_col in range(len(progeny_array[0])):
                if (progeny_array[prog][k_col] != 0) and (parent_array[par][k_col] != 0):
                    score_array[prog][par] += 1

    score_frame = pd.DataFrame(data=score_array, index=progeny_names, columns=parent_names)
    return score_frame


def import_kmer_count_table(file):
    data = pd.read_csv(file, sep=' ', header=0, na_values=" ")
    return data


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--parent-counts", help="Parent Kmer count table")
    parser.add_argument("-c", "--progeny-counts", help="Progeny Kmer count table")
    parser.add_argument("-s", "--save-file", help="where to save score file to")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug
    return parser.parse_args()


main()
