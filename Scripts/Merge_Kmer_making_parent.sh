#!/bin/bash
#HEADER FOR SUBMITTED SCRIPTS
#SBATCH --job-name=KMER_Merger_Parent
#SBATCH --partition=batch 
#SBATCH --ntasks=1
#SBATCH  --nodes=1 
#SBATCH  --cpus-per-task=10
#SBATCH --time=72:00:00
#SBATCH --export=NONE
#SBATCH --mem=50gb
#SBATCH --mail-user=drt83172@uga.edus
#SBATCH --mail-type=END,FAIL
#SBATCH --output=%x_%j.out 
#SBATCH --error=%x_%j.err 

# echo
# echo "Job ID: $PBS_JOBID"
# echo "Queue:  $PBS_QUEUE"
# echo "Cores:  $PBS_NP"
# echo "Nodes:  $(cat $PBS_NODEFILE | sort -u | tr '\n' ' ')"
# echo "mpirun: $(which mpirun)"
# echo

# cd $PBS_O_WORKDIR #to use the directory from which the job is submitted as the working directory (where to find input files or binaries)

#loading modules
module load parallel/20200422-GCCcore-8.3.0
module load Miniconda2/4.7.10

# THis script uses BBtools to maerge the reades together and then will use KMC3 to find Kmers
# conda envirnment to run this  is called KMERS
# contains bbtools and kmc3 
conda_env=KMERS
. $(conda info --root)/etc/profile.d/conda.sh 
conda activate $conda_env

#Settig variables for directories 
Data=/scratch/drt83172/Wallace_lab/TallFescue/Data
InterFiles=$Data/InterFiles
RawParentData=$Data/RawParentData
RawParentData_Merged=$Data/RawParentData_Merged
KMC_Parent_Data=$Data/KMC_Parent_Data

#Making directories if they are not made yet 
if [ ! -e $Data ] ; then mkdir $Data; fi
if [ ! -e $InterFiles ] ; then mkdir $InterFiles; fi
if [ ! -e $RawParentData ] ; then mkdir $RawParentData; fi
if [ ! -e $RawParentData_Merged ] ; then mkdir $RawParentData_Merged; fi
if [ ! -e $KMC_Parent_Data ] ; then mkdir $KMC_Parent_Data; fi

# Merging the raw data together  
>$InterFiles/DeleteMe.txt
for line in $(ls $RawParentData | sort | uniq | cut -d "-" -f 1)
do
echo "zcat ${RawParentData}/${line}-2944_R1_.fq.gz ${RawParentData}/${line}-2944_R2_.fq.gz | gzip > ${RawParentData_Merged}/${line}.fq.gz" >> $InterFiles/DeleteMe.txt
done
cat $InterFiles/DeleteMe.txt | parallel --jobs 2 --progress

for line in $(ls $RawParentData | sort | uniq | cut -d "-" -f 1)
do
# K is kmer size, m is memory allocated, t is threads
kmc -k25 -m30 -fq -t10 ${RawParentData_Merged}/${line}.fq.gz ${KMC_Parent_Data}/${line}_KMER_data ${RawParentData_Merged}
#This takes the files generated from the program before this one and makes a table of all the kmers
kmc_tools transform ${KMC_Parent_Data}/${line}_KMER_data dump -s ${KMC_Parent_Data}/${line}_KMERS.txt
done
