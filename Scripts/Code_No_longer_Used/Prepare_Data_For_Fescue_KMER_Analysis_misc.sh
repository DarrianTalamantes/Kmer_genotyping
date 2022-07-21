#This is a bunch of one liners that create files that are used in other programs

#This is used to make a list of all KMERS that we can later use in python. It was used for practice.
cat E_coli_KMERS.txt KMERS1.txt KMERS2.txt | cut -f 1 | sort | uniq > KMER_List.txt
# This is used to make a key from the progeny file names ive been using to the way progeny are labelled in the feild.
cat progeny_parent.csv | cut -d "," -f 2,5 | cut -d "_" -f 1-17 | sed s'/L002/L002_KMERS_Rready.txt/'g > progeny_key.csv
# THis is used to turn the Usless progeny json file into something usefull
cat Usless_Parents.txt | sedg' | sed 's/]//g' | sed 's/"//g' | sed 's/,/\n/g' > Usless_list.txt






