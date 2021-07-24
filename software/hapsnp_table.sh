#!/bin/bash
#SBATCH --job-name=hapsnp_table
#SBATCH -A tdlong_lab
#SBATCH -p standard

# sbatch hapsnp_table.sh
# $1 = path to gzipped_vcf_file
# $2 = save_data folder path
# $3 = software folder path
# $4 = path to subset list of individuals
# $5 = path to full list of individuals

gzipped_vcf_file=$1
save_data=$2
software=$3
subset_individuals_address=$4
full_individuals_address=$5



# GET CHROMOSOME NUMBER
arr_chromosome=($(echo $gzipped_vcf_file | tr -s '/' ' '))
chromosome=${arr_chromosome[5]}


wrangler="$software/hapsnp_tables_setup_manipulate.py"
module load bcftools/1.10.2
module load htslib/1.10.2
module load python/3.8.0


# SNP TABLE
snp_out="$save_data/snp_table_$chromosome.tsv"
if test -f "$snp_out"; then
	echo $snp_out exists
else
	zcat $gzipped_vcf_file | bcftools query -f '%POS [%DS{0} ]\n' | python $wrangler "snp" $subset_individuals_address $full_individuals_address > $snp_out
	echo $chromosome snp_table complete.
fi


# HAP TABLE
hap_out="$save_data/hap_table_$chromosome.tsv"
if test -f "$hap_out"; then
	echo $hap_out exists
else
	zcat $gzipped_vcf_file | bcftools query -f '%POS [%HD{0} %HD{1} %HD{2} %HD{3} %HD{4} %HD{5} %HD{6} %HD{7} ]\n' | python $wrangler "hap" $subset_individuals_address $full_individuals_address > $hap_out
	echo $chromosome hap_table complete.
fi