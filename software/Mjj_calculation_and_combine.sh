#!/bin/bash
#SBATCH --job-name=Mjj_calculation_and_combine
#SBATCH -A tdlong_lab
#SBATCH -p standard

# sbatch Mjj_calculation_and_combine.sh
# $1 = save_data folder path
# $2 = software folder path
# $3 = path to input file of haplotype calls from the entire genome
# $4 = path to subset list of individuals
# $5 = path to full list of individuals
# $6 = path to output file of relatedness values

save_data=$1
software=$2
all_genome_haplotypes=$3
subset_individuals_address=$4
full_individuals_address=$5
Mjj_values=$6



# SUM
Mjj_sum="$save_data/Mjj_sum.tsv"
module load python/3.8.0
# calculates the sum of the Mjj values of the crosses
summing_program="$software/Mjj_sum.py"
cat $all_genome_haplotypes | python $summing_program $subset_individuals_address $full_individuals_address > $Mjj_sum

# rm $all_genome_haplotypes




# AVERAGE
# calculate the average of all of the sum values
averaging_program="$software/Mjj_average.py"
cat $Mjj_sum | python $averaging_program $subset_individuals_address > $Mjj_values

# rm $Mjj_sum