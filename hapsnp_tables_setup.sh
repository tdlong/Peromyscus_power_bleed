#!/bin/bash
#SBATCH --job-name=hapsnp_tables_setup
#SBATCH -A tdlong_lab        ## account to charge
#SBATCH -p standard          ## partition/queue name

# FOR CHROMOSOME 23:
# sbatch /share/adl/pnlong/mouseproject/software/hapsnp_tables_setup.sh "/share/adl/tdlong/mouse_GWAS/STI8/Chr23/stitch.Chr23.vcf.gz"
# FOR CHROMOSOME 19 (NEGATIVE CONTROL):
# sbatch /share/adl/pnlong/mouseproject/software/hapsnp_tables_setup.sh "/share/adl/tdlong/mouse_GWAS/STI8/Chr19/stitch.Chr19.vcf.gz"


# this program will take the input of the gzipped vcf file we want to read, and the final folder we want to store the outputs in
# and output a snp table and haplotype table

# extract what chromosome we are working with
arr_chromosome=($(echo $1 | tr -s '/' ' '))
chromosome=${arr_chromosome[5]}

# output files
save_data="/share/adl/pnlong/mouseproject/save_data"
snp_out="$save_data/snp_table_$chromosome.tsv"
hap_out="$save_data/hap_table_$chromosome.tsv"
rm $snp_out
rm $hap_out

# save location of python manipulator file for readability
wrangle="/share/adl/pnlong/mouseproject/software/hapsnp_tables_setup_manipulate.py"

# load modules for vcf file reading
module load bcftools/1.10.2
module load htslib/1.10.2
module load python/3.8.0

# SNP TABLE
zcat $1 | bcftools query -f '%POS [%DS{0} ]\n' | python $wrangle "snp" > $snp_out

# HAP TABLE
zcat $1 | bcftools query -f '%POS [%HD{0} %HD{1} %HD{2} %HD{3} %HD{4} %HD{5} %HD{6} %HD{7} ]\n' | python $wrangle "hap" > $hap_out