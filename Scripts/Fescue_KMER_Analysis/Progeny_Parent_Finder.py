import numpy as np
import argparse
import pandas as pd
import subprocess

def main():
    # # use line below to test without redoing long step
    # best_parents = pd.read_csv('/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/predicted_parents_genetics.csv', index_col=0)
    # print(best_parents.iloc[2969])
    parse_args()

    # # Run R the first time to get All_centers.txt and Predicted_Parents.txt
    # # Args are random seed, files 1-4
    # subprocess.call(['Rscript', '/home/drt83172/Documents/Tall_fescue/Kmer_analyses/Scripts/Kmer_genotyping/Scripts/Score _analysis_auto.R', "10", "/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/R_Files/Score_table.csv", "/home/drt83172/Documents/Tall_fescue/half_key_parents.txt", "/home/drt83172/Documents/Tall_fescue/half_key_progeny.txt", "/home/drt83172/Documents/Tall_fescue/progeny_key.csv"])

    # # importing data for the first time
    # centers = import_data("/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/R_Files/All_centers.txt")
    # predicted = import_data("/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/R_Files/Predicted_Parents.txt")
    dead = import_data_csv("/home/drt83172/Documents/Tall_fescue/Plant_Info/Dead_Progeny.csv")

    # # Here we run the k-means grouping "resamples" amount of times and only keep ones that appear 90% or more of those times
    # resamples = 100
    # best_parents = resampling_kmeans(centers, predicted, resamples)
    # best_parents.to_csv('/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/predicted_parents_genetics.csv')
    # usable_parents_genetics = find_usable_parents(best_parents, dead)
    # usable_parents_genetics.to_csv('/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/usable_predicted_parents_genetics.csv')
    # print(len(usable_parents_genetics), "genetics only")

    # # # Confirms genetic data with maternal list and throws out those that dont work with both
    # best_parents = pd.read_csv('/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/predicted_parents_genetics.csv', index_col=0)
    # # print(best_parents.iloc[2969])
    # double_confirmed_parents = maternal_list_confirmer(best_parents)
    # double_confirmed_parents.to_csv('/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/predicted_parents_double.csv')
    # usable_double_confirmed = find_usable_parents(double_confirmed_parents, dead)
    # usable_double_confirmed.to_csv('/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/usable_predicted_parents_double.csv')
    # print(len(usable_double_confirmed), "double")
    #
    # # # adds maternal list to genetic data
    # best_parents = pd.read_csv('/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/predicted_parents_genetics.csv', index_col=0)
    # # print(best_parents.iloc[2969])
    # maternal_list_added = maternal_list_adder(best_parents)
    # maternal_list_added.to_csv('/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/predicted_parents_mat_added.csv')
    # usable_maternal_list_added = find_usable_parents(maternal_list_added, dead)
    # usable_maternal_list_added.to_csv('/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/usable_predicted_parents_mat_added.csv')
    # print(len(usable_maternal_list_added), "maternal list added")

    # # Uses the usable data frames created from find_usable_parents(double confirmed method to make a list of reciprical parents
    # # Semi unfinished
    # double_usable = pd.read_csv(
    #     "/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/usable_predicted_parents_double.csv", index_col=0)
    # reciprocal(double_usable)

    # # I use this to test how moving the cutoff value in the resampling method effects the result
    best_parents = pd.read_csv('/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/predicted_parents_genetics.csv',
                               index_col=0)
    resample = np.loadtxt("/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/consistancy.txt", dtype=float)
    beat_this = .85
    found_parents = quick_cutoff(resample, beat_this, best_parents)
    usable_parents_genetics = find_usable_parents(found_parents, dead)

def quick_cutoff(resample, beat_this, best_parents):
    resample_final = np.zeros((len(resample), len(resample[0])), dtype=int)
    # # This is where we make cut offs from the resample table and create the final resample table that we use after this.
    for row in range(len(resample)):
        for col in range(len(resample[0])):
            consistency = resample[row][col]
            if consistency >= beat_this:
                predicted_parent = column_finder(col, 2)
                for position in range(len(resample_final[0])):
                    position2 = resample_final[row][position]
                    if position2 == 0:
                        resample_final[row][position] = predicted_parent
                        break
    resample_final = pd.DataFrame(resample_final, index=best_parents.index,
                                  columns=best_parents.columns)
    return resample_final


# # Finding the reciprocal pairs of progeny. Must use double confirmed data to ensure we have data with maternal parents known
def reciprocal(predicted_parents):
    progeny_names = predicted_parents.index
    predicted_parents_array = predicted_parents.to_numpy()
    known_parent = np.zeros((len(progeny_names), 1), dtype='U100')

    for i in range(len(predicted_parents_array)):
        known_parent[i] = str(predicted_parents_array[i][0]) + "-" + str(predicted_parents_array[i][1])
    known_parent_set = set()
    for i in range(len(predicted_parents_array)):
        known_one = known_parent[i]
        known_one = known_one[-1]
        known_parent[i] = known_one
        known_parent_set.add(known_one)

    reciprocal_parents_set = set()
    for i in range(len(known_parent)):
        known_one = known_parent[i]
        known_one = known_one[-1]
        reciprocal = known_one[4:7] + known_one[3] + known_one[0:3]
        reciprocal_parents_set.add(reciprocal)

    for i in range(len(known_parent)):
        pair = known_parent[i]
        pair = pair[-1]
        if pair not in reciprocal_parents_set:
            reciprocal_parents_set.add(pair)

    return reciprocal_parents_set
    # setOfMarks = set(listOfMarks)


# # This method translates between column numbers and parent name
# # Input "True" if you are using a key first and anything else if using the value first
def column_finder(input, switcher):
    input = str(input)
    dic = {
        "301": "0",
        "302": "1",
        "303": "2",
        "304": "3",
        "305": "4",
        "306": "5",
        "307": "6",
        "308": "7",
        "310": "8",
        "312": "9",
        "313": "10",
        "314": "11",
        "315": "12",
        "316": "13",
        "318": "14",
        "319": "15",
        "320": "16",
    }
    if switcher == True:
        column = dic[input]
    else:
        column = list(dic.keys())[list(dic.values()).index(input)]

    return column


# # This method uses the parent finder 100 times and only keeps progeny parent pairs with statistical significance.
def resampling_kmeans(centers, predicted, resmaples):
    found_parents_genetics = parent_finder(centers, predicted)
    resample = np.zeros((len(found_parents_genetics), len(found_parents_genetics.columns)), dtype=int)
    beat_this = resmaples*.8
    for seed in range(resmaples):
        subprocess.call(['Rscript',
                         '/home/drt83172/Documents/Tall_fescue/Kmer_analyses/Scripts/Kmer_genotyping/Scripts/Score _analysis_auto.R',
                         str(seed), "/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/R_Files/Score_table.csv",
                         "/home/drt83172/Documents/Tall_fescue/half_key_parents.txt",
                         "/home/drt83172/Documents/Tall_fescue/half_key_progeny.txt",
                         "/home/drt83172/Documents/Tall_fescue/progeny_key.csv"])
        centers = import_data("/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/R_Files/All_centers.txt")
        predicted = import_data("/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/R_Files/Predicted_Parents.txt")
        found_parents_genetics = parent_finder(centers, predicted)
        found_parents_genetics_array = found_parents_genetics.to_numpy()
        for row in range(len(resample)):
            for col in range(len(resample[0])):
                predicted_parent = found_parents_genetics_array[row][col]
                if predicted_parent == 0:
                    break
                parent_column = column_finder(predicted_parent, True)
                parent_column = int(parent_column)
                resample[row][parent_column] += 1
    resample_final = np.zeros((len(found_parents_genetics), len(found_parents_genetics.columns)), dtype=int)
    # # This is where we make cut offs from the resample table and create the final resample table that we use after this.
    for row in range(len(resample)):
        for col in range(len(resample[0])):
            consistency = resample[row][col]
            if consistency >= beat_this:
                predicted_parent = column_finder(col, 2)
                for position in range(len(resample_final[0])):
                    position2 = resample_final[row][position]
                    if position2 == 0:
                        resample_final[row][position] = predicted_parent
                        break
    resample_final = pd.DataFrame(resample_final, index=found_parents_genetics.index, columns=found_parents_genetics.columns)
    np.savetxt("/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/consistancy.txt", resample, delimiter="\t")
    return resample_final


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
                genetic_list_array[row][col] = int(parent)
                break
    end_product = pd.DataFrame(genetic_list_array, index=genetic_list.index)
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
# # Also gets rid of all dead progeny
def find_usable_parents(predicted_parents, dead):
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
    print(len(too_many_parents), " have too many parents")
    print(len(not_enough_parents), " have too few parents")
    for i in range(len(usless_parents)):
        predicted_parents = predicted_parents.drop(usless_parents[i])

    # # Uses list of dead progeny and gets rid of them
    deads = 0
    deadlist = np.zeros((len(dead), 1), dtype='U100')
    for i in range(len(deadlist)):
        deadlist[i] = dead[i][0] + "-" + dead[i][1]
    for i in range(len(deadlist)):
        dead_one = deadlist[i]
        dead_one = dead_one[-1]
        if dead_one not in usless_parents:
            predicted_parents = predicted_parents.drop(dead_one)
            deads += 1
    print(deads, " got taken out for being dead af")
    return predicted_parents


def import_data(file):
    data = pd.read_csv(file, sep=' ', header=0, na_values=" ")
    return data

def import_data_csv(file):
    data = np.loadtxt(file, delimiter=",", dtype=str)
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
