#!/bin/bash
#HEADER FOR SUBMITTED SCRIPTS
#SBATCH --job-name=Progeny_Filtering
#SBATCH --partition=wallace_p 
#SBATCH  --nodes=1 
#SBATCH --ntasks-per-node=5
#SBATCH --time=60:00:00
#SBATCH --export=NONE
#SBATCH --mem=120gb
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


# # # To make Kmers run the Merge_Kmer_making_x.sh script first.
# # # Run program
# # # Step .5
# # # This step puts kmers that we will be working with into their own direcorty so we only work with them
# rm $Working_Kmers_Parents/*.txt
# rm $Working_Kmers_Progeny/*.txt
# 
# # Insert the parent list and the progeny lists here
Progeny_List=314x312_Progeny_Final.txt
Parent_List=314x312_Parents.txt
cross=314x312
# 
# >$InterFiles/DeleteMe01.txt
# for file in $(cat $Kmer_Lists/${Progeny_List})
# do
# echo "cp $Progeny_KMERS/${file} $Working_Kmers_Progeny" >> $InterFiles/DeleteMe01.txt
# done
# cat $InterFiles/DeleteMe01.txt | parallel --jobs 3 --progress
# 
# >$InterFiles/DeleteMe02.txt
# for file in $(cat $Kmer_Lists/${Parent_List})
# do
# echo "cp  $Parent_KMERS/${file} $Working_Kmers_Parents" >> $InterFiles/DeleteMe02.txt
# done
# cat $InterFiles/DeleteMe02.txt | parallel --jobs 2 --progress

# # Step 1
# # finds the the numbers at which a binomial distribution is 2 standard deviations away from the mean.
# tot_Prog=$(cat ${Kmer_Lists}/${Progeny_List} | wc -l) # Makes each trial of the bionomial distibution have 148 individuals
# trials=100 # This makes a bionomial distribution of 100 individuals
# percent_success=.5
# python Binomial_Distributer.py -tot $tot_Prog -t $trials -p $percent_success -s $InterFiles/stdarray.txt
# lowlim=$(awk 'NR==1 {print; exit}' $InterFiles/stdarray.txt)
# uplim=$(awk 'NR==2 {print; exit}' $InterFiles/stdarray.txt)
# echo "python Binomial_Distributer.py -tot $tot_Prog -t $trials -p $percent_success -s $InterFiles/stdarray.txt"
# 
# # # Step 2
# # # Filters progeny kmers using the lower limit
# python upper_lower_filter.py -k $Working_Kmers_Progeny -l $lowlim -s $InterFiles/kmer_progeny_filterlower.txt
# echo "python upper_lower_filter.py -k $Working_Kmers_Progeny -l $lowlim -s $InterFiles/kmer_progeny_filterlower.txt"
# 
# # # Step 3
# # # Filters new set of progeny Kmers using the higher limit
# python upper_lower_filter.py -k $Working_Kmers_Progeny -u $uplim -s $InterFiles/kmer_progeny_filterupper.txt
# echo "python upper_lower_filter.py -k $Working_Kmers_Progeny -u $uplim -s $InterFiles/kmer_progeny_filterupper.txt"

# # Step 4
#  #  Ensure only kmers found in both progeny filter lists are in final progeny file
# python filter_one_dic_using_another.py -p $InterFiles/kmer_progeny_filterupper.txt -c $InterFiles/kmer_progeny_filterlower.txt -s $InterFiles/kmer_progeny_filteredfinal.txt
# echo "python filter_one_dic_using_another.py -p $InterFiles/kmer_progeny_filterupper.txt -c $InterFiles/kmer_progeny_filterlower.txt -s $InterFiles/kmer_progeny_filteredfinal.txt"

# # Step 5
# # Filters parent Kmers to ensure they are only present in one parent
# python upper_lower_filter.py -k $Working_Kmers_Parents -u 1 -s $InterFiles/parent_kmers_filtered.txt
# echo "python upper_lower_filter.py -k $Working_Kmers_Parents -u 1 -s $InterFiles/parent_kmers_filtered.txt"
# 
# #Step 6
# # Makes master Kmer list by ensuring only kmers that appear in both progeny and parent lists show up 
# python filter_one_dic_using_another.py -p $InterFiles/parent_kmers_filtered.txt -c $InterFiles/kmer_progeny_filteredfinal.txt -s $Kmer_Lists/${cross}.txt
# echo "python filter_one_dic_using_another.py -p $InterFiles/parent_kmers_filtered.txt -c $InterFiles/kmer_progeny_filteredfinal.txt -s $Kmer_Lists/${cross}.txt"

# # Step 7
# # Use the newly created master kmer file to create a hapmat file using all kmer files
# python Making_Hapmat_FIle.py -pd $Working_Kmers_Parents -cd $Working_Kmers_Progeny -c $cross -m $Kmer_Lists/${cross}.txt -s $Hapmat_FIles/${cross}_hapmap.txt
# echo "python Making_Hapmat_FIle.py -pd $Working_Kmers_Parents -cd $Working_Kmers_Progeny -c $cross -m $Kmer_Lists/${cross}.txt -s $Hapmat_FIles/${cross}_hapmap.txt"

# # # Step 8
# # Create a genotype file from hapmap file that can be used with MSTmap
# # First we create the header necessary for mstMap and then we append the genotype file
# number_of_loci=$(awk 'END { print NR - 1 }' $Hapmat_FIles/${cross}_hapmap.txt )
# number_of_individual=$(head -n1 $Hapmat_FIles/${cross}_hapmap.txt | cut -f 12- | sed 's/[^\t]//g' | wc -c)
# 
# 
# > $Hapmat_FIles/${cross}_genotype.txt
# echo "population_type DH" >> $Hapmat_FIles/${cross}_genotype.txt
# echo "population_name LG" >> $Hapmat_FIles/${cross}_genotype.txt
# echo "distance_function kosambi" >> $Hapmat_FIles/${cross}_genotype.txt
# echo "cut_off_p_value 2.0" >> $Hapmat_FIles/${cross}_genotype.txt
# echo "no_map_dist 15.0" >> $Hapmat_FIles/${cross}_genotype.txt
# echo "no_map_size 0" >> $Hapmat_FIles/${cross}_genotype.txt
# echo "missing_threshold 1.00" >> $Hapmat_FIles/${cross}_genotype.txt
# echo "estimation_before_clustering no" >> $Hapmat_FIles/${cross}_genotype.txt
# echo "detect_bad_data yes" >> $Hapmat_FIles/${cross}_genotype.txt
# echo "objective_function COUNT" >> $Hapmat_FIles/${cross}_genotype.txt
# echo "number_of_loci $number_of_loci" >> $Hapmat_FIles/${cross}_genotype.txt
# echo "number_of_individual $number_of_individual" >> $Hapmat_FIles/${cross}_genotype.txt
# printf "\n" >> $Hapmat_FIles/${cross}_genotype.txt
# 
# cat  $Hapmat_FIles/${cross}_hapmap.txt | cut -f 1,12- | sed 's/_//g' | sed 's/\.//g' | sed 's/C/B/g' | sed 's/rs#/locus_name/g' > $InterFiles/intergenotype.txt
# cat $InterFiles/intergenotype.txt >> $Hapmat_FIles/${cross}_genotype.txt

# # Step 9
# # Use MSTmap to make a genetic map 
MSTMap.exe $Hapmat_FIles/${cross}_genotype.txt $Maps/$cross_map.txt



# # Step 8
# #This step must be done out of sapelo2. Convert the hapmap file into a VCF file using Tassel

# # Step 9 
# #  Follow the HetMap Pipeline (The filtering steps here seem to do nothing for us due to the nature of our data so we may want to just run mstMap?











