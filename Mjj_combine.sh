#!/bin/bash
#SBATCH --job-name=Mjj_combine
#SBATCH -A tdlong_lab
#SBATCH -p standard

# sbatch /share/adl/pnlong/mouseproject/software/Mjj_combine.sh "/share/adl/pnlong/mouseproject/save_data/Mjj_sums"

folder="/share/adl/pnlong/mouseproject"
Mjj_sums=$1

#combine all of the sum files into one
Mjj_sums_combined="$Mjj_sums/Mjj_sums_combined.tsv"
cat $Mjj_sums/Mjj_sum.*.tsv > $Mjj_sums_combined

#calculate the average of all of the (now combined) sum values
module load python/3.8.0
averaging_program="$folder/software/Mjj_average.py"
cat $Mjj_sums_combined | python $averaging_program > $folder/save_data/Mjj_values_all_pnl.tsv

rm -rf $Mjj_sums