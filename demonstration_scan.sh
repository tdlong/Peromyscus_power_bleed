#!/bin/bash
#SBATCH --job-name=bleeding_time_demo
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --array=0-23
#SBATCH --time=2-00:00:00

# sbatch /share/adl/pnlong/mouseproject/software/demonstration_scan.sh phenotype shuffle_phenotype

folder="/share/adl/pnlong/mouseproject"

demonstration_trait="$folder/demonstration_trait"
mkdir $demonstration_trait

if test $2 = "FALSE"; then
	# for normal scan (with normal Y values):
	results_folder="$demonstration_trait/$1"
elif test $2 = "TRUE"; then
	# for control scan (with shuffled Y values):
	results_folder="$demonstration_trait/$1/control"
else
	echo invalid shuffle_phenotype value.
fi

mkdir $results_folder
snp_based_out="$results_folder/snpbased"
mkdir $snp_based_out
hap_based_out="$results_folder/hapbased"
mkdir $hap_based_out


save_data="$folder/save_data"
scanner="$folder/software/demonstration_scan.R"


completed_tests="$results_folder/completed_tests.txt"
# create completed scans list if it hasn't been created yet
printf "" >> $completed_tests

gzipped_vcf_files=($(cat $save_data/gzipped_vcf_files.txt | tr -s '\n' ' '))
chromosomes=gzipped_vcf_files
for i in "${!gzipped_vcf_files[@]}"; do

	# GET CHROMOSOME NUMBER
	arr_chromosome=($(echo ${gzipped_vcf_files[$i]} | tr -s '/' ' '))
	chromosomes[$i]=${arr_chromosome[5]}

done

# RUN QTL SCAN
module load R/3.6.2
R --vanilla -f $scanner --args phenotype_address="$demonstration_trait/phenotype.tsv" phenotype=$1 shuffle_phenotype=$2 snp_table_address="$save_data/snp_table_${chromosomes[$SLURM_ARRAY_TASK_ID]}.tsv" hap_table_address="$save_data/hap_table_${chromosomes[$SLURM_ARRAY_TASK_ID]}.tsv" kinship_address="$save_data/Mjj_values_all_pnl.tsv" batch_size=100 rows_in_chunk=5000 drop_PC_less_than=0.05 LOD_threshold=2 snp_out=$snp_based_out hap_out=$hap_based_out completed_files=$completed_tests

