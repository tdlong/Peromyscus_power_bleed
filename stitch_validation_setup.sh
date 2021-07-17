#!/bin/bash
#SBATCH --job-name=stitch_validation_setup
#SBATCH -A tdlong_lab
#SBATCH -p standard

# sbatch /share/adl/pnlong/mouseproject/software/stitch_validation_setup.sh
folder="/share/adl/pnlong/mouseproject"




# RNA
big2="/share/adl/tdlong/mouse_GWAS/Ana/big2.txt.gz"
rna_out="$folder/save_data/rnaseq.tsv"
zcat $big2 | tr -s " " "\t" | head -n 1 > $rna_out # get correct column names
zcat $big2 | tr -s " " "\t" | cut --complement -f 1 | tail -n +2 >> $rna_out # remove line number column, get rid of now wrong column names due to spacings




# DNA
stitch_genomewide="$folder/save_data/stitch_genomewide.tsv"
rm $stitch_genomewide
# header
echo -e CHROM\\tPOS\\tID\\tDNA_GENO > $stitch_genomewide

module load bcftools/1.10.2
module load htslib/1.10.2
module load python/3.8.0
gzipped_vcf_files=($(cat $folder/save_data/gzipped_vcf_files.txt | tr -s '\n' ' '))
for gzipped_vcf_file in "${gzipped_vcf_files[@]}"; do
	# body
	zcat $gzipped_vcf_file | bcftools query -f '%CHROM %POS [%DS{0} ]\n' | python $folder/software/stitch_validation_setup.py >> $stitch_genomewide
done

# zcat /share/adl/tdlong/mouse_GWAS/Ana/all_genos.txt.gz | tr -s " " "\t" | tail -c +2 > $folder/save_data/stitch_genomewide.tsv # remove first # character