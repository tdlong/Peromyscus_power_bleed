#!/bin/bash
#SBATCH --job-name=Mjj_setup
#SBATCH -A tdlong_lab        ## account to charge
#SBATCH -p standard          ## partition/queue name

# sbatch /share/adl/pnlong/mouseproject/software/Mjj_setup.sh

folder="/share/adl/pnlong/mouseproject"

all_genome_haplotypes="$folder/save_data/haplotypes_all_genome.tsv"
rm $all_genome_haplotypes

gzipped_vcf_files=($(cat $folder/save_data/gzipped_vcf_files.txt | tr '\n' ' '))

# load modules for vcf file reading
module load bcftools/1.10.2
module load htslib/1.10.2

# combine the VCFs into one file that can be read by a python program later
for gzipped_vcf_file in ${gzipped_vcf_files[@]}; do
	# so we will only take 0.1% of each chromosome, seeing as any more would be A) overkill and B) take too long
	zcat $gzipped_vcf_file | bcftools query -f '%CHROM %POS [%HD{0} %HD{1} %HD{2} %HD{3} %HD{4} %HD{5} %HD{6} %HD{7} ]\n' | awk 'NR % 1000 == 0' | tr -s " " "\t" >> $all_genome_haplotypes
	echo $gzipped_vcf_file has been processed
done