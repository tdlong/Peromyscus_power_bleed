#!/bin/bash
#SBATCH --job-name=Mjj_calculation
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --array=1-100

# sbatch /share/adl/pnlong/mouseproject/software/Mjj_calculation.sh

folder="/share/adl/pnlong/mouseproject"

all_genome_haplotypes="$folder/save_data/haplotypes_all_genome.tsv" # file that we will be drawing the data from

Mjj_sums_folder="$folder/save_data/Mjj_sums"

mkdir $Mjj_sums_folder

module load python/3.8.0
# split up the file of all haplotypes into 100 different ones, where we will run the python program which calculates the sum of the Mjj values of the crosses of the given file
summing_program="$folder/software/Mjj_sum.py"
cat $all_genome_haplotypes | sed -n $SLURM_ARRAY_TASK_ID~100p | python $summing_program > $Mjj_sums_folder/Mjj_sum.$SLURM_ARRAY_TASK_ID.tsv