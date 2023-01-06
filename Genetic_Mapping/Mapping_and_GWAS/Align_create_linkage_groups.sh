#!/bin/bash


# Code by Darrian Talamantes 
# Objective: Align Kmer files to Lolium pernne genome using BWA. Then I will seperate each chromosome into new linkage groups.


# conda envirnment is called Alignment on local PC

############## Variables ###############
# directories
alignment="/scratch/drt83172/Wallace_lab/TallFescue/Data/Kmer_Alignments/"
kmers="/home/drt06/Documents/Tall_fescue/Kmer_Genotyping/Kmer_genotyping/Genetic_Mapping/Data/Kmer_Files/"
Intermediates="/home/drt06/Documents/Tall_fescue/Kmer_Genotyping/Kmer_genotyping/Genetic_Mapping/Data/Intermediates/"
# files
ref_genome="/home/drt06/Documents/Tall_fescue/Kmer_Genotyping/Kmer_genotyping/Genetic_Mapping/Data/Ref_Genomes/Lolium_perenne/Loliumpernne_genome.fasta"


#creates a bwa database
#bwa index Loliumpernne_genome.fasta

# make KMC files into fastaformat
# idea is to put fasta header before every single kmer
# ex: kmerNumber_progName : >01_301-41-2
for file in $(ls $kmers)
do
    arrIN=(${file//./ }) 
    python fasta_maker.py -kmer $kmers/$file -s $Intermediates/${arrIN[0]}.fasta
done

# Carry out BWA alignments between kmers and lolium pernne
 


# # bwa commands to get alignment correct
# for file in $(ls $Intermediates | grep ".fasta")
# do
# a=$(echo $file | cut -d "." -f 1)
# bwa aln $ref_genome $Intermediates${file} > $Intermediates${a}.sai
# bwa samse $ref_genome $Intermediates${a}.sai $Intermediates${file} > $Intermediates${a}.sam
# done

# # seperate files by chromosome to use those as the linkage groups
# for file in (ls $Intermediates | grep ".sam")

# bwa aln /home/drt06/Documents/Tall_fescue/Kmer_Genotyping/Kmer_genotyping/Genetic_Mapping/Data/Ref_Genomes/Lolium_perenne/Loliumpernne_genome.fasta poopy.txt > result.sai 
# bwa samse /home/drt06/Documents/Tall_fescue/Kmer_Genotyping/Kmer_genotyping/Genetic_Mapping/Data/Ref_Genomes/Lolium_perenne/Loliumpernne_genome.fasta result.sai poopy.txt > align.sam
# bwa aln ref.fa short_read.fq > aln_sa.sai
# bwa samse bwa samse ref.fa aln_sa.sai short_read.fq > aln-se.sam
 