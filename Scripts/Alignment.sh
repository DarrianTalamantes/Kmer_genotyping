#!/bin/bash
#HEADER FOR SUBMITTED SCRIPTS
#SBATCH --job-name=Aligner
#SBATCH --partition=wallace_p 
#SBATCH --ntasks=1
#SBATCH  --nodes=1 
#SBATCH  --cpus-per-task=5
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
# conda install -c bioconda bwa=7.17 
# conda install -c bioconda samtools=1.9

# Making the correct directories 
Home=/scratch/drt83172/Wallace_lab/TallFescue/Data
Raw_Progeny=$Home/RawProgenyData
Aligned_progeny=$Home/Aligned_progeny
Aligned_progeny_bam=$Home/Aligned_progeny_bam
Sorted_progeny_bam=$Home/Sorted_progeny_bam

if [ ! -e $Raw_Progeny ] ; then mkdir $Raw_Progeny; fi
if [ ! -e $Aligned_progeny ] ; then mkdir $Aligned_progeny; fi





# #Setting other variables
PROCS=5
refGenome=$Home/Refrence/scaffolds_UC307-2944.fasta

#################################### Start of Code ##########################################
# bwa index -a bwtsw $refGenome


for sample in $(ls $Raw_Progeny | head | cut -d "_" -f 1-17) 
do

# bwa mem -t $PROCS $refGenome $Raw_Progeny/${sample}_R1_001.fastq.gz $Raw_Progeny/${sample}_R2_001.fastq.gz > $Aligned_progeny/$sample.aligned.sam

samtools view -S -b $Aligned_progeny/$sample.aligned.sam > $Sorted_progeny_bam/$sample.bam
samtools sort $Sorted_progeny_bam/$sample.bam -o $Sorted_progeny_bam/$sample.sorted.bam

done














