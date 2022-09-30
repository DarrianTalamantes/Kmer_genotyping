#!/bin/bash
#HEADER FOR SUBMITTED SCRIPTS
#SBATCH --job-name=Progeny_Filtering
#SBATCH --partition=highmem_p 
#SBATCH  --nodes=1 
#SBATCH --ntasks-per-node=10
#SBATCH --time=150:00:00
#SBATCH --export=NONE
#SBATCH --mem=350gb
#SBATCH --mail-user=drt83172@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/drt83172/Wallace_lab/TallFescue/Scripts/OutFiles/%x_%j.out 
#SBATCH --error=/scratch/drt83172/Wallace_lab/TallFescue/Scripts/OutFiles/%x_%j.err 

# echo
# echo "Job ID: $PBS_JOBID"
# echo "Queue:  $PBS_QUEUE"
# echo "Cores:  $PBS_NP"
# echo "Nodes:  $(cat $PBS_NODEFILE | sort -u | tr '\n' ' ')"
# echo "mpirun: $(which mpirun)"
# echo

# cd $PBS_O_WORKDIR #to use the directory from which the job is submitted as the working directory (where to find input files or binaries)

#Peak memory consumption of thsi program was 756.11 GB
# #First three steps took 9 hours to complete.

#loading modules
module load Miniconda2/4.7.10
module load MSTmap/20210830-GCC-8.3.0
# activating conda
conda_env=KMER_Filter
. $(conda info --root)/etc/profile.d/conda.sh 
source activate $conda_env
module load MSTmap/20210830-GCC-8.3.0

# conda install -c anaconda numpy=1.19.2
# conda install -c conda-forge matplotlib=3.4.2
# conda install -c anaconda pandas=1.1.5
# conda install -c conda-forge argparse=1.4.0
# conda install -c conda-forge parallel=20210422
# conda install -c r r=3.6.0
# conda install -c r r-reshape=0.8.8 
# conda install -c r r-tidyverse=1.2.1
# conda install -c conda-forge r-vegan=2.5_7
 
# conda install -c bioconda mstmap 

# Making the correct directories 
# The next line is for practice 
# Home=/scratch/drt83172/Wallace_lab/TallFescue/Data/Practice

Home=/scratch/drt83172/Wallace_lab/TallFescue/Data
Scripts=/scratch/drt83172/Wallace_lab/TallFescue/Scripts
Progeny_KMERS=$Home/KMC_Progeny_Data
Parent_KMERS=$Home/KMC_Parent_Data
InterFiles=$Home/InterFiles
Final_Kmers=$Home/Final_Kmers
Parents_R_ready=$Home/Parents_R_ready
Progeny_R_ready=$Home/Progeny_R_ready
R_outputs=$Home/R_outputs
Score_table=$Home/Score_table
Working_Kmers_Parents=$Home/Working_Kmers_Parents 
Working_Kmers_Progeny=$Home/Working_Kmers_Progeny
Kmer_Lists=$Home/Kmer_Lists
Hapmat_FIles=$Home/Hapmat_Files
HetMap=$Scripts/Filter_Kmers/HetMap
Maps=$Home/Maps
Keys=$Home/Keys


if [ ! -e $Progeny_KMERS ] ; then mkdir $Progeny_KMERS; fi
if [ ! -e $InterFiles ] ; then mkdir $InterFiles; fi
if [ ! -e $Parent_KMERS ] ; then mkdir $Parent_KMERS; fi
if [ ! -e $Final_Kmers ] ; then mkdir $Final_Kmers; fi
if [ ! -e $Parent_stage_2 ] ; then mkdir $Parent_stage_2; fi
if [ ! -e $Parents_R_ready ] ; then mkdir $Parents_R_ready; fi
if [ ! -e $Progeny_R_ready ] ; then mkdir $Progeny_R_ready; fi
if [ ! -e $R_outputs ] ; then mkdir $R_outputs; fi
if [ ! -e $Score_table ] ; then mkdir $Score_table; fi
if [ ! -e $Working_Kmers_Parents ] ; then mkdir $Working_Kmers_Parents; fi
if [ ! -e $Working_Kmers_Progeny ] ; then mkdir $Working_Kmers_Progeny; fi
if [ ! -e $Kmer_Lists ] ; then mkdir $Kmer_Lists; fi
if [ ! -e $Hapmat_FIles ] ; then mkdir $Hapmat_FIles; fi
if [ ! -e $HetMap ] ; then mkdir $HetMap; fi
if [ ! -e $Maps ] ; then mkdir $Maps; fi
if [ ! -e $Keys ] ; then mkdir $Keys; fi


# # To make Kmers run the Merge_Kmer_making_x.sh script first.
# # Run program


# # Step 0
# # This step will create all lists of all predidcted progeny for every combination of parents
parental_combos=$Keys/All_combos.csv
predicted_parents=$Keys/predicted_parents_double_subset.csv  # This file is made from "Progeny_Parent_Finder.py"
progeny_key=$Keys/progeny_key.csv

rm $Keys/*_code.txt
rm $Keys/*_Parents.txt
while read -r line
do 
    a=$(echo $line | cut -d "," -f 1 )
    b=$(echo $line | cut -d "," -f 2 | sed s'/ //g')
    >${Keys}/${a}x${b}_Parents.txt
    echo "cat $predicted_parents | cut -d "," -f 1-3 | egrep '${a},${b}|${b},${a}'>$Keys/${a}x${b}_code.txt" >> commands.txt
    echo "${a}_KMERS.txt">>${Keys}/${a}x${b}_Parents.txt
    echo "${b}_KMERS.txt">>${Keys}/${a}x${b}_Parents.txt
done <$parental_combos
cat commands.txt | parallel --jobs 4 --progress


rm $Keys/*_Progeny_Final.txt
for line in $(ls $Keys| grep "code.txt")
do
    arrIN=(${line//_/ }) # makes the variable into an array that I sepetate by "_"
    cross=${arrIN[0]}  # calls the first section of the array
    Rscript --vanilla match_all_progeny_to_parents_p2.R $progeny_key $Keys/$line ${Keys}/${cross}_Progeny_Final.txt
    echo "Rscript --vanilla match_all_progeny_to_parents_p2.R $progeny_key $line ${Keys}/${cross}_Progeny_Final.txt"
done

echo "Step 0 complete"
# # Step .5
# # This step puts kmers that we will be working with into their own direcorty so we only work with them


# # This next part of the script will set 3 variables that are needed for the rest of the script (may be able to set these to static variables if needed
for line in $(ls $Keys | grep "_Progeny_Final.txt") #When creating these files in the previous step any that have 0 progeny are not made into files
do
    rm $Working_Kmers_Parents/*.txt
    rm $Working_Kmers_Progeny/*.txt
    arrIN=(${line//_/ }) 

    # # These are my special 3 variables
    cross=${arrIN[0]}
    Progeny_List=$line
    Parent_List=${cross}_Parents.txt
    echo
    echo "Current cross is $cross"


    # #Insert the parent list and the progeny lists here. Keeping this here to possibly be able to run the script peice by peice later.
    # # cross=312x314
    # # Progeny_List=${cross}_Progeny_Final.txt
    # # Parent_List=${cross}_Parents.txt

    >$InterFiles/DeleteMe01.txt
    for file in $(cat $Keys/${Progeny_List} | cut -d " " -f 2)
    do
        echo "cp $Progeny_KMERS/${file} $Working_Kmers_Progeny" >> $InterFiles/DeleteMe01.txt
    done
    cat $InterFiles/DeleteMe01.txt | parallel --jobs 3 --progress

    >$InterFiles/DeleteMe02.txt
    for file in $(cat $Keys/${Parent_List})
    do
        echo "cp  $Parent_KMERS/${file} $Working_Kmers_Parents" >> $InterFiles/DeleteMe02.txt
    done
    cat $InterFiles/DeleteMe02.txt | parallel --jobs 2 --progress







    echo "Step 0.5 complete"
    # # Step 1
    # # Uses a bionomial distribution assuming hardy weinburg equilibrium to find upper and lower limits on what would be 2 std's away from the mean 
    # # for how many times a kmer should appear thorught progeny files
    # # finds the the numbers at which a binomial distribution is 2 standard deviations away from the mean.
    tot_Prog=$(cat $Keys/${Progeny_List} | wc -l) # Makes each trial of the bionomial distibution have individuals equal to the amount of the cross
    trials=100 # This makes a bionomial distribution of 100 individual trials
    percent_success=.5 # This makes the distribution binomial
    python Binomial_Distributer.py -tot $tot_Prog -t $trials -p $percent_success -s $InterFiles/stdarray.txt
    lowlim=$(awk 'NR==1 {print; exit}' $InterFiles/stdarray.txt)
    uplim=$(awk 'NR==2 {print; exit}' $InterFiles/stdarray.txt)
    echo "python Binomial_Distributer.py -tot $tot_Prog -t $trials -p $percent_success -s $InterFiles/stdarray.txt"
    
    echo "Step 1 complete"
    # # Step 2
    # # Filters progeny kmers using the lower limit
    python upper_lower_filter.py -k $Working_Kmers_Progeny -l $lowlim -s $InterFiles/kmer_progeny_filterlower.txt
    echo "python upper_lower_filter.py -k $Working_Kmers_Progeny -l $lowlim -s $InterFiles/kmer_progeny_filterlower.txt"

    echo "Step 2 complete"
    # # Step 3
    # # Filters new set of progeny Kmers using the higher limit
    python upper_lower_filter.py -k $Working_Kmers_Progeny -u $uplim -s $InterFiles/kmer_progeny_filterupper.txt
    echo "python upper_lower_filter.py -k $Working_Kmers_Progeny -u $uplim -s $InterFiles/kmer_progeny_filterupper.txt"

    echo "Step 3 complete"
    # # Step 4
    # #  Ensure only kmers found in both progeny filter lists are in final progeny file
    python filter_one_dic_using_another.py -p $InterFiles/kmer_progeny_filterupper.txt -c $InterFiles/kmer_progeny_filterlower.txt -s $InterFiles/kmer_progeny_filteredfinal.txt
    echo "python filter_one_dic_using_another.py -p $InterFiles/kmer_progeny_filterupper.txt -c $InterFiles/kmer_progeny_filterlower.txt -s $InterFiles/kmer_progeny_filteredfinal.txt"
    
    echo "Step 4 complete"
    # # Step 5
    # # Filters parent Kmers to ensure they are only present in one parent
    
    FILE=$InterFiles/${cross}_parent_kmers_filtered.txt
    if [ -f "$FILE" ]; then
        echo "$FILE has already been made"
    else 
        echo "Creating $FILE ."
        python upper_lower_filter.py -k $Working_Kmers_Parents -u 1 -s $InterFiles/${cross}_parent_kmers_filtered.txt
        echo "python upper_lower_filter.py -k $Working_Kmers_Parents -u 1 -s $InterFiles/${cross}_parent_kmers_filtered.txt"
    fi

    echo "Step 5 complete"
    # #Step 6
    # # Makes master Kmer list by ensuring only kmers that appear in both progeny and parent lists show up 
    python filter_one_dic_using_another.py -p $InterFiles/${cross}_parent_kmers_filtered.txt -c $InterFiles/kmer_progeny_filteredfinal.txt -s $Kmer_Lists/${cross}.txt
    echo "python filter_one_dic_using_another.py -p $InterFiles/${cross}_parent_kmers_filtered.txt -c $InterFiles/kmer_progeny_filteredfinal.txt -s $Kmer_Lists/${cross}.txt"

    echo "Step 6 complete"
    # # Step 7
    # # Use the newly created master kmer file to create a hapmat file using      all kmer files
    python Making_Hapmat_FIle.py -pd $Working_Kmers_Parents -cd $Working_Kmers_Progeny -c $cross -m $Kmer_Lists/${cross}.txt -s $Hapmat_FIles/${cross}_hapmap.txt
    echo "python Making_Hapmat_FIle.py -pd $Working_Kmers_Parents -cd $Working_Kmers_Progeny -c $cross -m $Kmer_Lists/${cross}.txt -s $Hapmat_FIles/${cross}_hapmap.txt"

    echo "Step 7 complete"
    # # # Step 8 
    # # Removing Redundant Kmers
    python Remove_Redundant_Kmers.py -s $Hapmat_FIles/${cross}_hapmapNonRedun.txt -hap $Hapmat_FIles/${cross}_hapmap.txt
    echo "python Remove_Redundant_Kmers.py -s $Hapmat_FIles/${cross}_hapmapNonRedun.txt -hap $Hapmat_FIles/${cross}_hapmap.txt"

    echo "Step 8 complete"
    # # Step 9
    # # Create a genotype file from hapmap file that can be used with MSTmap
    # # First we create the header necessary for mstMap and then we append the genotype file
    number_of_loci=$(awk 'END { print NR - 1 }' $Hapmat_FIles/${cross}_hapmapNonRedun.txt )
    number_of_individual=$(head -n1 $Hapmat_FIles/${cross}_hapmapNonRedun.txt | cut -f 12- | sed 's/[^\t]//g' | wc -c)

    > $Hapmat_FIles/${cross}_genotype.txt
    echo "population_type DH" >> $Hapmat_FIles/${cross}_genotype.txt
    echo "population_name LG" >> $Hapmat_FIles/${cross}_genotype.txt
    echo "distance_function kosambi" >> $Hapmat_FIles/${cross}_genotype.txt
    echo "cut_off_p_value 2.0" >> $Hapmat_FIles/${cross}_genotype.txt
    echo "no_map_dist 15.0" >> $Hapmat_FIles/${cross}_genotype.txt
    echo "no_map_size 0" >> $Hapmat_FIles/${cross}_genotype.txt
    echo "missing_threshold 1.00" >> $Hapmat_FIles/${cross}_genotype.txt
    echo "estimation_before_clustering no" >> $Hapmat_FIles/${cross}_genotype.txt
    echo "detect_bad_data yes" >> $Hapmat_FIles/${cross}_genotype.txt
    echo "objective_function COUNT" >> $Hapmat_FIles/${cross}_genotype.txt
    echo "number_of_loci $number_of_loci" >> $Hapmat_FIles/${cross}_genotype.txt
    echo "number_of_individual $number_of_individual" >> $Hapmat_FIles/${cross}_genotype.txt
    printf "\n" >> $Hapmat_FIles/${cross}_genotype.txt

    cat  $Hapmat_FIles/${cross}_hapmapNonRedun.txt | cut -f 1,12- | sed 's/_/Zz/g' | sed 's/\.//g' |  sed 's/\tC\t/\tB\t/g' | sed 's/\/C/\/B/g' | sed 's/\tC/\tB/g'|  sed 's/rs#/locus_name/g' > $InterFiles/intergenotype.txt
    cat $InterFiles/intergenotype.txt >> $Hapmat_FIles/${cross}_genotype.txt

done
    
echo "Step 9 complete"
# # Step 10
# # Use MSTmap to make a genetic map. could not get this working on the cluster
# MSTMap.exe $Hapmat_FIles/${cross}_genotype.txt $Maps/$cross_map.txt














