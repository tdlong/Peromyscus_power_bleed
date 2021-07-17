#!/bin/bash
#SBATCH --job-name=hapsnp
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --array=1
#SBATCH --time=2-00:00:00

# rerun if not all the tests are done, hope for better nodes

# sbatch /share/adl/pnlong/mouseproject/software/hapsnp_calculation.sh test_type save_file snp_table hap_table kinship_matrix snp_out hap_out completed_tests_list
# $1 = test_type -> ("single", "multiple", "rare")
# $2 = address of save_file
# $3 = address of snp_table
# $4 = address of hap_table
# $5 = address of the kinship matrix data
# $6 = address of the folder of the output for the snp-based test
# $7 = address of the folder of the output for the haplotype-based test
# $8 = address of list of completed files

chromosome=$(echo $2 | tr -cd '0-9') # so it only extracts the 23 or 19 part from the save_file name

snp_based_out="$6/Chr${chromosome}_${1}_snpbased"
hap_based_out="$7/Chr${chromosome}_${1}_hapbased"

# create completed_files_list if it hasn't been created yet
printf "" >> $8

# determine the csnp_effect_percent from 
csnp_effect=$(echo "scale=3; $(echo $8 | tr -dc '0-9') / 1000" | bc)

program="/share/adl/pnlong/mouseproject/software/hapsnp_calculation.R"
module load R/3.6.2
R --vanilla -f $program --args test_type=$1 save_file=$2 input_snp=$SLURM_ARRAY_TASK_ID snp_table_address=$3 hap_table_address=$4 kinship_address=$5 batch_size=100 rows_in_chunk=5000 drop_PC_less_than=0.05 LOD_threshold=2 csnp_effect_percent=$csnp_effect snp_out=$snp_based_out hap_out=$hap_based_out completed_files=$8