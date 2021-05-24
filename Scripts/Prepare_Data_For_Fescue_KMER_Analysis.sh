#This is used to make a list of all KMERS that we can later use in python
cat E_coli_KMERS.txt KMERS1.txt KMERS2.txt | cut -f 1 | sort | uniq > KMER_List.txt
