#!/bin/bash
#SBATCH --job-name=hapsnp_tests
#SBATCH -A tdlong_lab
#SBATCH -p free
#SBATCH --output=hapsnp_tests_runner.txt

# FOR 2.5% EFFECT
# sbatch /share/adl/pnlong/mouseproject/software/hapsnp_tests.sh "/share/adl/pnlong/mouseproject/results025" "/share/adl/pnlong/mouseproject/save_data/Mjj_values_all_pnl.tsv"
# FOR 5% EFFECT
# sbatch /share/adl/pnlong/mouseproject/software/hapsnp_tests.sh "/share/adl/pnlong/mouseproject/results050" "/share/adl/pnlong/mouseproject/save_data/Mjj_values_all_pnl.tsv"

results_folder=$1

folder="/share/adl/pnlong/mouseproject"
program="$folder/software/hapsnp_calculation.sh"


# format of the sbatch command:
# sbatch /share/adl/pnlong/software/hapsnp_calculation.sh test_type save_file snp_table snp_out hap_out

snp_table_Chr23="$folder/save_data/snp_table_Chr23.tsv"
hap_table_Chr23="$folder/save_data/hap_table_Chr23.tsv"
kinship_matrix=$2
completed_files_list="$results_folder/completed_tests.txt"

########### CHROMOSOME 23
save23="$folder/save_data/hapsnp_Chr23.Rsave"

# SINGLE
mkdir $results_folder/Chr23_single_snpbased
mkdir $results_folder/Chr23_single_hapbased
sbatch $program "single" $save23 $snp_table_Chr23 $hap_table_Chr23 $kinship_matrix "$results_folder/Chr23_single_snpbased" "$results_folder/Chr23_single_hapbased" $completed_files_list

# MULTIPLE
mkdir $results_folder/Chr23_multiple_snpbased
mkdir $results_folder/Chr23_multiple_hapbased
sbatch $program "multiple" $save23 $snp_table_Chr23 $hap_table_Chr23 $kinship_matrix "$results_folder/Chr23_multiple_snpbased" "$results_folder/Chr23_multiple_hapbased" $completed_files_list

# RARE
mkdir $results_folder/Chr23_rare_snpbased
mkdir $results_folder/Chr23_rare_hapbased
sbatch $program "rare" $save23 $snp_table_Chr23 $hap_table_Chr23 $kinship_matrix "$results_folder/Chr23_rare_snpbased" "$results_folder/Chr23_rare_hapbased" $completed_files_list

########### CHROMOSOME 19 (NEGATIVE CONTROL)
save19="$folder/save_data/hapsnp_Chr19.Rsave"

# SINGLE
mkdir $results_folder/Chr19_single_snpbased
mkdir $results_folder/Chr19_single_hapbased
sbatch $program "single" $save19 $snp_table_Chr23 $hap_table_Chr23 $kinship_matrix "$results_folder/Chr19_single_snpbased" "$results_folder/Chr19_single_hapbased" $completed_files_list

# MULTIPLE
mkdir $results_folder/Chr19_multiple_snpbased
mkdir $results_folder/Chr19_multiple_hapbased
sbatch $program "multiple" $save19 $snp_table_Chr23 $hap_table_Chr23 $kinship_matrix "$results_folder/Chr19_multiple_snpbased" "$results_folder/Chr19_multiple_hapbased" $completed_files_list

# RARE
mkdir $results_folder/Chr19_rare_snpbased
mkdir $results_folder/Chr19_rare_hapbased
sbatch $program "rare" $save19 $snp_table_Chr23 $hap_table_Chr23 $kinship_matrix "$results_folder/Chr19_rare_snpbased" "$results_folder/Chr19_rare_hapbased" $completed_files_list
