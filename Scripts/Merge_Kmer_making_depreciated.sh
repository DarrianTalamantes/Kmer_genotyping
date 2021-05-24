#bin/bash
# THis script uses BBtools to maerge the reades together and then will use KMC3 to find Kmers
# conda envirnment to run this  is called KMERS
# contains kmc3 
# Variables Here
RawData="/home/drt83172/Documents/Tall_fescue/RawData"
RawData_Merged="/home/drt83172/Documents/Tall_fescue/RawData_Merged"
KMCData="/home/drt83172/Documents/Tall_fescue/KMC_Data"

# Merging the raw data together  
for line in $(ls $RawData | sort | uniq | cut -d "-" -f 1)
do
zcat ${RawData}/${line}-2944_R1_.fq.gz ${RawData}/${line}-2944_R2_.fq.gz > ${RawData_Merged}/${line}.fq
gzip ${RawData_Merged}/${line}.fq
done


for line in $(ls $RawData | sort | uniq | cut -d "-" -f 1)
do
# K is kmer size, m is memory allocated, t is threads
kmc -k25 -m30 -fq -t8 ${RawData_Merged}/${line}.fq.gz ${KMCData}/${line}_KMER_data ${RawData_Merged}
#This takes the files generated from the program before this one and makes a table of all the kmers
kmc_tools transform ${KMCData}/${line}_KMER_data dump -s ${KMCData}/${line}_KMERS.txt
done
