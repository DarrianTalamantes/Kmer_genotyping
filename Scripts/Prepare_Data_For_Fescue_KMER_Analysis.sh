#This is used to make a list of all KMERS that we can later use in python
cat E_coli_KMERS.txt KMERS1.txt KMERS2.txt | cut -f 1 | sort | uniq > KMER_List.txt
# This is used to make a key from the progeny file names ive been using to the way progeny are labelled in the feild.
cat progeny_parent.csv | cut -d "," -f 2,5 | cut -d "_" -f 1-17 | sed s'/L002/L002_KMERS_Rready.txt/'g > progeny_key.csv






