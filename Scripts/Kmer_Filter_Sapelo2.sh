#!/bin/bash
#HEADER FOR SUBMITTED SCRIPTS
#SBATCH --job-name=Filtering
#SBATCH --partition=wallace_p 
#SBATCH --ntasks=1
#SBATCH  --nodes=1 
#SBATCH  --cpus-per-task=10
#SBATCH --time=72:00:00
#SBATCH --export=NONE
#SBATCH --mem=120gb
#SBATCH --mail-user=drt83172@uga.edu
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
module load Miniconda2/4.7.10

# activating conda
conda_env=KMER_Filter
. $(conda info --root)/etc/profile.d/conda.sh 
source activate $conda_env

# conda install -c anaconda numpy=1.19.2
# conda install -c conda-forge matplotlib=3.4.2
# conda install -c anaconda pandas=1.1.5
#conda install -c conda-forge argparse=1.4.0
# Run program

# Filters for progeny kmers that appear at least x times
python 01_progeny_filter.py
# Makes one huuuuuuuuuge pressence file for parent data.
python 02_make_big_parent_kmer_file.py
# Filters parent kmers to ones that only appear once through all the files
python 03_parent_filter.py
# Finds filtered parent kmers that will match up with the ones in the filtered progeny files
python 04_filter_par_using_pro.py
