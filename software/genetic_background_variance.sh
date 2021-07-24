#!/bin/bash
#SBATCH --job-name=genetic_background_variance
#SBATCH -A tdlong_lab
#SBATCH -p standard

# sbatch genetic_background_variance.sh
# $1 = save_data folder path
# $2 = software folder path
# $3 = path to list of gzipped vcf filepaths
# $4 = number of SNPs in consideration for genetic background variance (100000)
# $5 = path to full list of individuals
# $6 = path to output for genetic background variance

save_data=$1
software=$2
gzipped_vcf_files_list=$3
number_of_SNPs_for_gbv=$4
full_individuals_address=$5
output=$6





if test -f "$output"; then # if the genetic background varaiance file exists, then
  echo $output exists

else # if the genetic background varaiance file does not yet exist, then

  all_genome_dosages="$save_data/entire_genome_dosages.txt"
  rm $all_genome_dosages

  module load bcftools/1.10.2
  module load htslib/1.10.2

  gzipped_vcf_files=($(cat $gzipped_vcf_files_list | tr '\n' ' '))

  # The amount of loci from each chromosome that we want to contribute to out genetic background variance calculation
  # determined by dividing the amount of files we are processing by the desired total amount of loci for our genetic background variance calculation (the first argument)
  N=$(expr $number_of_SNPs_for_gbv / ${#gzipped_vcf_files[@]}) # N=4166

  for gzipped_vcf_file in ${gzipped_vcf_files[@]}; do

    temp_file="$all_genome_dosages.temp"
    zcat $gzipped_vcf_file | head -n 12 > $temp_file
    zcat $gzipped_vcf_file | tail -n +13 | shuf -n $N >> $temp_file

    cat $temp_file | bcftools query -f '%CHROM %POS [%DS{0} ]\n' >> $all_genome_dosages
    rm $temp_file
    echo $gzipped_vcf_file has been processed

  done



  # now to run our genetic background variance program on the file we just generated
  program="$software/genetic_background_variance.py"

  module load python/3.8.0
  cat $all_genome_dosages | python $program $full_individuals_address > $output
  rm $all_genome_dosages

fi