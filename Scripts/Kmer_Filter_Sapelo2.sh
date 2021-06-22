#!/bin/bash
#HEADER FOR SUBMITTED SCRIPTS
#SBATCH --job-name=Filtering
#SBATCH --partition=highmem_p 
#SBATCH --ntasks=1
#SBATCH  --nodes=1 
#SBATCH  --cpus-per-task=10
#SBATCH --time=168:00:00
#SBATCH --export=NONE
#SBATCH --mem=800gb
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

# Making the correct directories 
Home=/scratch/drt83172/Wallace_lab/TallFescue/Data
Progeny_KMERS=$Home/KMC_Progeny_Data/KMERS
InterFiles=$Home/InterFiles
Parent_KMERS=$Home/KMC_Parent_Data/KMERS
Final_Kmers=$Home/Kmers_Final
Parent_stage_2=$Home/Parent_stage_2
if [ ! -e $Progeny_KMERS ] ; then mkdir $Progeny_KMERS; fi
if [ ! -e $InterFiles ] ; then mkdir $InterFiles; fi
if [ ! -e $Parent_KMERS ] ; then mkdir $Parent_KMERS; fi
if [ ! -e $Final_Kmers ] ; then mkdir $Final_Kmers; fi
if [ ! -e $Parent_stage_2 ] ; then mkdir $Parent_stage_2; fi

### Below is for practice 
# Home=/scratch/drt83172/Wallace_lab/TallFescue/Data/Practice
# Progeny_KMERS=$Home/Progeny_KMERS
# InterFiles=$Home/InterFiles
# Parent_KMERS=$Home/Parent_KMERS
# Final_Kmers=$Home/Final_Kmers
# Parent_stage_2=$Home/Parent_stage_2
# if [ ! -e $Progeny_KMERS ] ; then mkdir $Progeny_KMERS; fi
# if [ ! -e $InterFiles ] ; then mkdir $InterFiles; fi
# if [ ! -e $Parent_KMERS ] ; then mkdir $Parent_KMERS; fi
# if [ ! -e $Final_Kmers ] ; then mkdir $Final_Kmers; fi
# if [ ! -e $Parent_stage_2 ] ; then mkdir $Parent_stage_2; fi

# Run program

# Filters for progeny kmers that appear at least x times
# python upper_lower_filter.py -k $Progeny_KMERS -l 6 -s $InterFiles/kmer_progeny_filter1.txt
# Need to filter one parent file at a time using script 2
>$InterFiles/DeleteMe.txt
for i in $(ls $Parent_KMERS | cut -d . -f 1)
do
echo "python filter_one_dic_using_another.py -p $Parent_KMERS/${i}.txt -c $InterFiles/kmer_progeny_filter1.txt -s $Parent_stage_2/${i}_Stage2.txt" >> $InterFiles/DeleteMe.txt
done
cat $InterFiles/DeleteMe.txt | parallel --jobs 3 --progress
# use script one on second stage parent files
python upper_lower_filter.py -k $Parent_stage_2 -u 1 -s $Final_Kmers/usefull_kmers.txt
 









