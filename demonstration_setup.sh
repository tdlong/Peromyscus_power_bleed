#!/bin/bash
#SBATCH --job-name=demonstration_setup
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --output=demonstration_setup.txt
#SBATCH --array=0-23
#SBATCH --time=2-00:00:00

# sbatch /share/adl/pnlong/mouseproject/software/demonstration_setup.sh

folder="/share/adl/pnlong/mouseproject"

save_data="$folder/save_data"
wrangler="$folder/software/hapsnp_tables_setup_manipulate.py"
scanner="$folder/software/demonstration_scan.R"


# GET GZIPPED_VCF_FILE THAT IM WORKING ON
gzipped_vcf_files=($(cat $save_data/gzipped_vcf_files.txt | tr -s '\n' ' '))
gzipped_vcf_file=${gzipped_vcf_files[$SLURM_ARRAY_TASK_ID]}

# GET CHROMOSOME NUMBER
arr_chromosome=($(echo $gzipped_vcf_file | tr -s '/' ' '))
chromosome=${arr_chromosome[5]}

# load modules for vcf file reading
module load bcftools/1.10.2
module load htslib/1.10.2
module load python/3.8.0

# SNP TABLE
snp_out="$save_data/snp_table_$chromosome"
if test -f "$snp_out.tsv"; then
    echo $snp_out.tsv exists.
else
	zcat $gzipped_vcf_file | bcftools query -f '%POS [%DS{0} ]\n' | python $wrangler "snp" > $snp_out.tsv
	echo $chromosome snp_table complete.
fi

# HAP TABLE
hap_out="$save_data/hap_table_$chromosome"
if test -f "$hap_out.tsv"; then
	echo $hap_out.tsv exists.
else
	zcat $gzipped_vcf_file | bcftools query -f '%POS [%HD{0} %HD{1} %HD{2} %HD{3} %HD{4} %HD{5} %HD{6} %HD{7} ]\n' | python $wrangler "hap"> $hap_out.tsv
	echo $chromosome hap_table complete.
fi
