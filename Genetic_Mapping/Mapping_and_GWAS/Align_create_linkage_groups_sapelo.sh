#!/bin/bash

#HEADER FOR SUBMITTED SCRIPTS
#SBATCH --job-name=Kmer_alignment
#SBATCH --partition=batch_p 
#SBATCH  --nodes=1 
#SBATCH --ntasks-per-node=4
#SBATCH --time=30:00:00
#SBATCH --export=NONE
#SBATCH --mem=50gb
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


# Code by Darrian Talamantes 
# Objective: Align Kmer files to Lolium pernne genome using BWA. Then I will seperate each chromosome into new linkage groups.


# conda envirnment is called Alignment on local PC

############## Activating envirment ###############
conda_env=Alignment
. $(conda info --root)/etc/profile.d/conda.sh 
source activate $conda_env



############## Variables ###############
# directories
alignment="/scratch/drt83172/Wallace_lab/TallFescue/Data/Kmer_Alignments"
# files
ref_genome="/scratch/drt83172/Wallace_lab/TallFescue/Data/Refrence/Lolium_pernne/Loliumpernne_genome.fasta"


#creates a bwa database
#bwa index Loliumpernne_genome.fasta

# make KMC files into fastaformat
# idea is to put fasta header before every single kmer
# ex: kmerNumber_progName : >01_301-41-2

for file in $(ls /home/drt06/Documents/Tall_fescue/Kmer_Genotyping/Kmer_genotyping/Genetic_Mapping/Data/Kmer_Files)
do
a=$(echo $line | cut -d "." -f 1)
cat $file | cut -f 1 | sed a\$a



done


# bwa commands to get alignment correct

# bwa aln <genome.fasta> <kmers.fasta> > <result.sai>  
# bwa samse <genome.fasta> <result.sai> <kmers.fasta> -f <mapping.sam> 