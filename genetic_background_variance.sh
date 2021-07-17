#!/bin/bash
#SBATCH --job-name=genetic_background_variance
#SBATCH -A tdlong_lab        ## account to charge
#SBATCH -p standard          ## partition/queue name

# sbatch /share/adl/pnlong/mouseproject/software/genetic_background_variance.sh 100000

folder="/share/adl/pnlong/mouseproject"

all_genome_dosages="$folder/save_data/entire_genome_dosages.txt"
rm $all_genome_dosages

module load bcftools/1.10.2
module load htslib/1.10.2

gzipped_vcf_files=($(cat $folder/save_data/gzipped_vcf_files.txt | tr '\n' ' '))

# The amount of loci from each chromosome that we want to contribute to out genetic background variance calculation
# determined by dividing the amount of files we are processing by the desired total amount of loci for our genetic background variance calculation (the first argument)
N=$(expr $1 / ${#gzipped_vcf_files[@]}) # N=4166

for gzipped_vcf_file in ${gzipped_vcf_files[@]}; do
  zcat $gzipped_vcf_file | bcftools query -f '%CHROM %POS [%DS{0} ]\n' | shuf -n $N >> $all_genome_dosages
  echo $gzipped_vcf_file has been processed
done

# now to run our program on the file we just generated
program="$folder/software/genetic_background_variance.py"
output="$folder/save_data/genetic_background_variance.tsv"

module load python/3.8.0
cat $all_genome_dosages | python $program > $output
rm $all_genome_dosages