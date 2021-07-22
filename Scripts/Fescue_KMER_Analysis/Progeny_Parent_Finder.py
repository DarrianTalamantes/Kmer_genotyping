import numpy as np
import argparse
import pandas as pd
import subprocess


def main():
    parse_args()
    # # importing data
    centers = import_data("/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/R_Files/All_centers.txt")
    predicted = import_data("/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/R_Files/Predicted_Parents.txt")
    
    # # Step one use imported files to make a list of predicted parents
    found_parents_genetics = parent_finder(centers, predicted)
    found_parents_genetics.to_csv('/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/predicted_parents_genetics.csv')
    usable_parents_genetics = find_usable_parents(found_parents_genetics)
    usable_parents_genetics.to_csv('/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/usable_predicted_parents_genetics.csv')

    # # Confirms genetic data with maternal list and throws out those that dont work with both
    double_confirmed_parents = maternal_list_confirmer(found_parents_genetics)
    double_confirmed_parents.to_csv('/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/predicted_parents_double.csv')
    usable_double_confirmed = find_usable_parents(double_confirmed_parents)
    usable_double_confirmed.to_csv('/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/usable_predicted_parents_double.csv')

    # # adds maternal list to genetic data
    maternal_list_added = maternal_list_adder(found_parents_genetics)
    maternal_list_added.to_csv('/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/predicted_parents_mat_added.csv')
    usable_maternal_list_added = find_usable_parents(maternal_list_added)
    usable_maternal_list_added.to_csv('/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/usable_predicted_parents_mat_added.csv')

# # This method takes the maternal list and will add it in without confirming it with the genetics
def maternal_list_adder(genetic_list):
    genetic_list_array = genetic_list.to_numpy()
    progeny_names = genetic_list.index
    for row in range(len(genetic_list_array)):
        temp = progeny_names[row]
        parent = temp[0:3]
        for col in range(len(genetic_list_array[0])):
            if genetic_list_array[row][col] == int(parent):
                break
            elif genetic_list_array[row][col] == 0:
                genetic_list_array[row][col] = parent
                genetic_list_array[row][16] = 99
                break

    end_product = pd.DataFrame(genetic_list_array, index=genetic_list.index)
    print(end_product)
    return end_product


def maternal_list_confirmer(genetic_list):
    genetic_list_array = genetic_list.to_numpy()
    progeny_names = genetic_list.index
    for row in range(len(genetic_list_array)):
        temp = progeny_names[row]
        parent = temp[0:3]
        x = 0
        for col in range(len(genetic_list_array[0])):
            if genetic_list_array[row][col] == int(parent):
                x = 1
                break
        if x != 1:
            for col in range(len(genetic_list_array[0])):
                genetic_list_array[row][col] = 0

    end_product = pd.DataFrame(genetic_list_array, index=genetic_list.index)
    print(end_product)
    return end_product


# # This method takes an input file of the group centers and a file listing progeny sorted into groups 1 or 2 for
# # each parent. It outputs a dataframe that shows the predicted parents for each progeny based only on genetics.
def parent_finder(centers, predicted):
    centers_array = centers.to_numpy()
    predicted_array = predicted.to_numpy()

    # print(predicted.columns[0])
    # print(predicted.index[0])
    true_parents = np.zeros((len(predicted_array), len(predicted_array[0])), dtype=int)
    # # This loop adds our predicted parents
    for parent in range((len(centers.columns))):
        parent_name = centers.columns[parent]
        if parent_name == "302":
            # print("found a 302")
            parent_name = 303
        elif parent_name == "303":
            # print("found a 303")
            parent_name = 302
        if centers_array[0][parent] > centers_array[1][parent] and centers_array[0][parent] > centers_array[2][parent] and centers_array[0][parent] > centers_array[3][parent]:
            group = 1
        elif centers_array[1][parent] > centers_array[0][parent] and centers_array[1][parent] > centers_array[2][parent] and centers_array[1][parent] > centers_array[3][parent]:
            group = 2
        elif centers_array[2][parent] > centers_array[0][parent] and centers_array[2][parent] > centers_array[1][parent] and centers_array[2][parent] > centers_array[3][parent]:
            group = 3
        elif centers_array[3][parent] > centers_array[0][parent] and centers_array[3][parent] > centers_array[1][parent] and centers_array[3][parent] > centers_array[2][parent]:
            group = 4
        for progeny in range(len(predicted.index)):
            if predicted_array[progeny][parent] == group:
                # print(group, parent_name)
                # print(predicted.index[progeny])
                for i in range(0, len(predicted_array[0])):
                    if int(true_parents[progeny][i]) == 0 and int(true_parents[progeny][0]) != int(parent_name):
                            true_parents[progeny][i] = parent_name
                            break
    end_product = pd.DataFrame(true_parents, index=predicted.index)
    return end_product

    # print(progeny, predicted_array[progeny][parent])


# # uses the output of parent_finder method to expunge progeny who do not have exactly 2 parents assigned
def find_usable_parents(predicted_parents):
    not_enough_parents = set()
    too_many_parents = set()
    for col in range(1, 3):
        if col == 1:
            for row in range(len(predicted_parents)):
                if int(predicted_parents.iloc[row, col]) == 0:
                    not_enough_parents.add(predicted_parents.index[row])
        if col == 2:
            for row in range(len(predicted_parents)):
                if int(predicted_parents.iloc[row, col]) != 0:
                    too_many_parents.add(predicted_parents.index[row])
    usless_parents = not_enough_parents.symmetric_difference(too_many_parents)  # the amount usless here is correct
    usless_parents = list(usless_parents)
    for i in range(len(usless_parents)):
        predicted_parents = predicted_parents.drop(usless_parents[i])
    return predicted_parents


def import_data(file):
    data = pd.read_csv(file, sep=' ', header=0, na_values=" ")
    return data


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--identifiers", help="Parent Kmer count table")
    parser.add_argument("-c", "--centers", help="Progeny Kmer count table")
    parser.add_argument("-s", "--save-file", help="where to save score file to")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug
    return parser.parse_args()


main()
