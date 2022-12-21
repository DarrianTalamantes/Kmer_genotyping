#!/bin/bash
# Objective: This script will create many lists of progeny files that correspond to each set of parents 
# Author Darrian Talamantes

#variables 
parental_combos=/home/drt06/Documents/Tall_fescue/Plant_Info/Match_Parent_Data/All_combos.csv
predicted_parents=/home/drt06/Documents/Tall_fescue/Plant_Info/Match_Parent_Data/usable_predicted_parents_double.csv
progeny_key=/home/drt06/Documents/Tall_fescue/Plant_Info/Match_Parent_Data/progeny_key.csv


# The file usable_predicted_parents_double.csv was created uisng the meternal list to check the predicted genetic parents
# It got rid of any progeny that did not have exactly two parents and leaves dead plants in to have a stronger analysis.

>commands.txt
while read -r line
do 
a=$(echo $line | cut -d "," -f 1 )
b=$(echo $line | cut -d "," -f 2 | sed s'/ //g')
>${a}x${b}_Parents.txt
echo $a
echo $b
echo "cat $predicted_parents | cut -d "," -f 1-3 | egrep '${a},${b}|${b},${a}'>${a}x${b}_code.txt" >> commands.txt
echo "${a}_KMERS.txt">>${a}x${b}_Parents.txt
echo "${b}_KMERS.txt">>${a}x${b}_Parents.txt
done <$parental_combos

cat commands.txt | parallel --jobs 4 --progress

for line in $(ls | grep "code.txt")
do
echo "line is "$line 
arrIN=(${line//_/ })
cross=${arrIN[0]}
Rscript --vanilla match_all_progeny_to_parents_p2.R $progeny_key $line ${cross}_Progeny_Final.txt
done

# This entire file will be incorperated into Kmer_Filter_Sapelo2.sh
