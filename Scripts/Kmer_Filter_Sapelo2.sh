#!/bin/bash
#HEADER FOR SUBMITTED SCRIPTS
#SBATCH --job-name=Score_maker
#SBATCH --partition=wallace_p 
#SBATCH --ntasks=1
#SBATCH  --nodes=1 
#SBATCH  --cpus-per-task=10
#SBATCH --time=24:00:00
#SBATCH --export=NONE
#SBATCH --mem=100gb
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


# Making the correct directories 
Home=/scratch/drt83172/Wallace_lab/TallFescue/Data
Progeny_KMERS=$Home/KMC_Progeny_Data/KMERS
InterFiles=$Home/InterFiles
Parent_KMERS=$Home/KMC_Parent_Data/KMERS
Final_Kmers=$Home/Kmers_Final
Parent_stage_2=$Home/Parent_stage_2
Parents_R_ready=$Home/Parents_R_ready
Progeny_R_ready=$Home/Progeny_R_ready
R_outputs=$Home/R_outputs


if [ ! -e $Progeny_KMERS ] ; then mkdir $Progeny_KMERS; fi
if [ ! -e $InterFiles ] ; then mkdir $InterFiles; fi
if [ ! -e $Parent_KMERS ] ; then mkdir $Parent_KMERS; fi
if [ ! -e $Final_Kmers ] ; then mkdir $Final_Kmers; fi
if [ ! -e $Parent_stage_2 ] ; then mkdir $Parent_stage_2; fi
if [ ! -e $Parents_R_ready ] ; then mkdir $Parents_R_ready; fi
if [ ! -e $Progeny_R_ready ] ; then mkdir $Progeny_R_ready; fi
if [ ! -e $R_outputs ] ; then mkdir $R_outputs; fi


# # Below is for practice 
# Home=/scratch/drt83172/Wallace_lab/TallFescue/Data/Practice
# Progeny_KMERS=$Home/Progeny_KMERS
# InterFiles=$Home/InterFiles
# Parent_KMERS=$Home/Parent_KMERS
# Final_Kmers=$Home/Final_Kmers
# Parent_stage_2=$Home/Parent_stage_2
# Parents_R_ready=$Home/Parents_R_ready
# Progeny_R_ready=$Home/Progeny_R_ready
# R_outputs=$Home/R_outputs
# Score_table=$Home/Score_table
# 
# 
# if [ ! -e $Progeny_KMERS ] ; then mkdir $Progeny_KMERS; fi
# if [ ! -e $InterFiles ] ; then mkdir $InterFiles; fi
# if [ ! -e $Parent_KMERS ] ; then mkdir $Parent_KMERS; fi
# if [ ! -e $Final_Kmers ] ; then mkdir $Final_Kmers; fi
# if [ ! -e $Parent_stage_2 ] ; then mkdir $Parent_stage_2; fi
# if [ ! -e $Parents_R_ready ] ; then mkdir $Parents_R_ready; fi
# if [ ! -e $Progeny_R_ready ] ; then mkdir $Progeny_R_ready; fi
# if [ ! -e $R_outputs ] ; then mkdir $R_outputs; fi
# if [ ! -e $Score_table ] ; then mkdir $Score_table; fi


# # Run program
# 
# # Filters for progeny kmers that appear at least x times
# python upper_lower_filter.py -k $Progeny_KMERS -l 6 -s $InterFiles/kmer_progeny_filter1.txt
# # Need to filter one parent file at a time using script 2
# >$InterFiles/DeleteMe.txt
# for i in $(ls $Parent_KMERS | cut -d . -f 1)
# do
# echo "python filter_one_dic_using_another.py -p $Parent_KMERS/${i}.txt -c $InterFiles/kmer_progeny_filter1.txt -s $Parent_stage_2/${i}_Stage2.txt" >> $InterFiles/DeleteMe.txt
# done
# cat $InterFiles/DeleteMe.txt | parallel --jobs 3 --progress
# # use script one on second stage parent files
# python upper_lower_filter.py -k $Parent_stage_2 -u 1 -s $Final_Kmers/usefull_kmers.txt
#  
# # Use the usefull_kmers list to extract usefull kmers from all other files
# >$InterFiles/EraseMe.txt
# for file in $(ls $Parent_stage_2 | cut -d . -f 1)
# do
# echo "python Getting_kmers_ready_for_R.py -i $Parent_stage_2/${file}.txt -l $Final_Kmers/usefull_kmers.txt -s $Parents_R_ready/${file}_Rready.txt" >> $InterFiles/EraseMe.txt
# done
# cat $InterFiles/EraseMe.txt | parallel --jobs 4 --progress
# 
# >$InterFiles/WipeMe.txt
# for file in $(ls $Progeny_KMERS | cut -d . -f 1)
# do
# echo "python Getting_kmers_ready_for_R.py -i $Progeny_KMERS/${file}.txt -l $Final_Kmers/usefull_kmers.txt -s $Progeny_R_ready/${file}_Rready.txt" >> $InterFiles/WipeMe.txt
# done
# cat $InterFiles/WipeMe.txt | parallel --jobs 4 --progress
# # commands are a directory for R ready parents and then R ready progeny, then output for parents and output for progeny.
# # R script only works in R studio and possibly on local computer. Does not work in sapelo. SImply add the outputs into sapelo if wanting to run on sapelo2.
# Rscript --vanilla Kmer_analysis.R /home/drt83172/Documents/Tall_fescue/Usefull_Kmers/Parents /home/drt83172/Documents/Tall_fescue/Usefull_Kmers/Progeny $R_outputs/R_parents.txt $R_outputs/R_progeny.txt
# 
# python Score_Table_creation.py -p $R_outputs/R_parents.txt -c $R_outputs/R_progeny.txt -s $R_outputs/Score_table.csv

# # This next code is to change the key we were given by flex seq to get know parents and progeny maternal pairs 
cat Parent_progeny.csv | cut -d "," -f 2 | cut -d "-" -f 1 > half_key_parents.txt
cat Parent_progeny.csv | cut -d "," -f 5 | cut -d "_" -f 1-17 | sed s'/L002/L002_KMERS_Rready.txt/'g > key_progeny_half.txt

# #  This creates another key for progeny sample code to customer code 
cat Parent_progeny.csv | cut -d "," -f 1,2 > Sample_to_customer_code.csv
cat Parent_progeny.csv | cut -d "," -f 2,5 | sed s'/L002_R1_001.fastq.gz/KMERS_Rready.txt/'g > progeny_key.csv 


# # This is how I found the depth that I will later incorporate into my data
bcftools stats -S Progeny_Names_VCF.txt UGA_149001_FlexSeqResults.vcf.gz | grep "PSC" | cut -f 3,10

# # Run this R script next with the score table output and the key edited with the code above
# Score _analysis.R

# # Final python program
python Progeny_Parent_Finder.py -i -c -s