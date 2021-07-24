#!/bin/bash
#SBATCH --job-name=scan
#SBATCH -A tdlong_lab
#SBATCH -p standard

# sh scan.sh
# $1 = save_data folder path
# $2 = software folder path
# $3 = path to gzipped vcf file of the chromosome that scans are performed on
# $4 = path to subset list of individuals
# $5 = path to full list of individuals
# $6 = path to phenotype file
# $7 = name of column for the phenotype in the phenotype file
# $8 = results folder path
# $9 = path to file where the filepath of a scan is printed to when the scan is completed
# $10 = genetic model (either "single", "multiple", or "rare")

save_data=${1}
software=${2}
gzipped_vcf_file=${3}
subset_individuals_address=${4}
full_individuals_address=${5}
phenotype_file=${6}
phenotype_column_name=${7}
results=${8}
completed_files_list=${9}
genetic_model=${10}


mkdir -p $save_data
mkdir -p $software



# MAKE NECESSARY SNP AND HAP TABLES
arr_chromosome=($(echo $gzipped_vcf_file | tr -s '/' ' '))
chromosome=${arr_chromosome[5]}
snp_out="$save_data/snp_table_$chromosome.tsv"
hap_out="$save_data/hap_table_$chromosome.tsv"

sh $software/hapsnp_table.sh $gzipped_vcf_file $save_data $software $subset_individuals_address $full_individuals_address
echo Necessary SNP and Haplotype data complete



# MAKE OUTPUT DIRECTORIES
mkdir -p $results

snpbased_out="$results/snpbased"
mkdir $snpbased_out
hapbased_out="$results/hapbased"
mkdir $hapbased_out

printf "" >> $completed_files_list

Mjj_values="$save_data/Mjj_values_all.tsv"



# RUN SCAN R PROGRAM
program="$software/scan.R"
module load R/3.6.2
R --vanilla -f $program --args test_type=$genetic_model phenotype_address=$phenotype_file phenotype_name=$phenotype_column_name individuals_address=$subset_individuals_address snp_table_address=$snp_out hap_table_address=$hap_out kinship_address=$Mjj_values batch_size=100 rows_in_chunk=5000 drop_PC_less_than=0.05 LOD_threshold=2 snp_out=$snpbased_out hap_out=$hapbased_out completed_files=$completed_files_list
echo SNP and Haplotype scan complete
