#!/bin/bash
#SBATCH --job-name=generate_Y
#SBATCH -A tdlong_lab
#SBATCH -p standard

# sh generate_Y.sh
# $1 = save_data folder path
# $2 = software folder path
# $3 = path to list of gzipped vcf filepaths
# $4 = path to gzipped vcf file of the chromosome for which we want our causative SNP from
# $5 = position of causative snp
# $6 = path to subset list of individuals
# $7 = path to full list of individuals

save_data=${1}
software=${2}
gzipped_vcf_files_list=${3}
number_of_SNPs_for_gbv=100000
gzipped_vcf_file=${4}
csnp_locus=${5}
subset_individuals_address=${6}
full_individuals_address=${7}


mkdir -p $save_data
mkdir -p $software



# GENETIC BACKGROUND VARIANCE
genetic_background_variance_output="$save_data/genetic_background_variance.tsv"
sh $software/genetic_background_variance.sh $save_data $software $gzipped_vcf_files_list $number_of_SNPs_for_gbv $full_individuals_address $genetic_background_variance_output
echo Genetic Background Variance complete



# MAKE NECESSARY SNP AND HAP TABLES
arr_chromosome=($(echo $gzipped_vcf_file | tr -s '/' ' '))
chromosome=${arr_chromosome[5]}
snp_out="$save_data/snp_table_$chromosome.tsv"

sh $software/hapsnp_table.sh $gzipped_vcf_file $save_data $software $subset_individuals_address $full_individuals_address
echo Necessary SNP data complete



# CREATE PHENOTYPE FILE
mkdir $save_data/phenotypes
module load R/3.6.2
R --vanilla -f $software/generate_Y.R --args output_prefix="$save_data/phenotypes/Y_" causative_snp_address=$snp_out causative_snp_locus=$csnp_locus range=25000 number_of_snps_per_csnp=10 individuals_address=$subset_individuals_address genetic_background_variance_address=$genetic_background_variance_output
echo Phenotype generated


