#!/bin/bash
#HEADER FOR SUBMITTED SCRIPTS
#SBATCH --job-name=KMER_Pro
#SBATCH --partition=wallace_p
#SBATCH --ntasks=1
#SBATCH  --nodes=1 
#SBATCH  --cpus-per-task=10
#SBATCH --time=72:00:00
#SBATCH --export=NONE
#SBATCH --mem=60gb
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
RawProgenyData=$Data/RawProgenyData
RawProgenyData_Merged=$Data/RawProgenyData_Merged
KMC_Progeny_Data=$Data/KMC_Progeny_Data

#Making directories if they are not made yet 
if [ ! -e $Data ] ; then mkdir $Data; fi
if [ ! -e $InterFiles ] ; then mkdir $InterFiles; fi
if [ ! -e $RawProgenyData ] ; then mkdir $RawProgenyData; fi
if [ ! -e $RawProgenyData_Merged ] ; then mkdir $RawProgenyData_Merged; fi
if [ ! -e $KMC_Progeny_Data ] ; then mkdir $KMC_Progeny_Data; fi

# Merging the raw data together  
>$InterFiles/DeleteMeProgeny.txt
for line in $(ls $RawProgenyData | cut -d "_" -f 1-17 | sort | uniq)
do
echo "zcat ${RawProgenyData}/${line}_R1_001.fastq.gz ${RawProgenyData}/${line}_R2_001.fastq.gz | gzip > ${RawProgenyData_Merged}/${line}.fq.gz" >> $InterFiles/DeleteMeProgeny.txt
done
cat $InterFiles/DeleteMeProgeny.txt | parallel --jobs 10 --progress

for line in $(ls $RawProgenyData | cut -d "_" -f 1-17 | sort | uniq)
do
# K is kmer size, m is memory allocated, t is threads
kmc -k25 -m30 -fq -t10 ${RawProgenyData_Merged}/${line}.fq.gz ${KMC_Progeny_Data}/${line}_KMER_data ${RawProgenyData_Merged}
#This takes the files generated from the program before this one and makes a table of all the kmers
kmc_tools transform ${KMC_Progeny_Data}/${line}_KMER_data dump -s ${KMC_Progeny_Data}/${line}_KMERS.txt
done

