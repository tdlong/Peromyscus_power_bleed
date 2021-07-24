#!/bin/bash
#SBATCH --job-name=generate_Mjj
#SBATCH -A tdlong_lab
#SBATCH -p standard

# sh generate_Mjj.sh
# $1 = save_data folder path
# $2 = software folder path
# $3 = path to list of gzipped vcf filepaths
# $4 = path to subset list of individuals
# $5 = path to full list of individuals

save_data=${1}
software=${2}
gzipped_vcf_files_list=${3}
subset_individuals_address=${4}
full_individuals_address=${5}


mkdir -p $save_data
mkdir -p $software



# SETUP
all_genome_haplotypes="$save_data/haplotypes_all_genome.tsv"
sh $software/Mjj_setup.sh $gzipped_vcf_files_list $all_genome_haplotypes $software
echo Setup complete



# GENERATE RELATEDNESS VALUES
Mjj_values="$save_data/Mjj_values_all.tsv"
sh $software/Mjj_calculation_and_combine.sh $save_data $software $all_genome_haplotypes $subset_individuals_address $full_individuals_address $Mjj_values
echo Kinship Matrix ready
