#!/bin/bash
#SBATCH --job-name=Mjj_setup
#SBATCH -A tdlong_lab
#SBATCH -p standard

# sbatch Mjj_setup.sh
# $1 = path to list of gzipped vcf filepaths
# $2 = path to input file of haplotype calls from the entire genome
# $3 = software folder path

gzipped_vcf_files_list=$1
all_genome_haplotypes=$2
software=$3



gzipped_vcf_files=($(cat $gzipped_vcf_files_list | tr '\n' ' '))


# load modules for vcf file reading
module load bcftools/1.10.2
module load htslib/1.10.2


# combine the VCFs into one file that can be read by a python program later
rm $all_genome_haplotypes
for gzipped_vcf_file in ${gzipped_vcf_files[@]}; do

	temp_file="$all_genome_haplotypes.temp"
	zcat $gzipped_vcf_file | head -n 12 > $temp_file
	zcat $gzipped_vcf_file | tail -n +13 | awk 'NR % 1000 == 0' >> $temp_file

	# so we will only take 0.1% of each chromosome, seeing as any more would be A) overkill and B) take too long
	cat $temp_file | bcftools query -f '%CHROM %POS [%HD{0} %HD{1} %HD{2} %HD{3} %HD{4} %HD{5} %HD{6} %HD{7} ]\n' | tr -s " " "\t" >> $all_genome_haplotypes
	rm $temp_file
	echo $gzipped_vcf_file has been processed
done