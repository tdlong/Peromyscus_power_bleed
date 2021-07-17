#!/bin/bash
#SBATCH --job-name=figures
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH -c 12
#SBATCH --output=figures.txt

# sbatch /share/adl/pnlong/mouseproject/software/figures.sh

folder="/share/adl/pnlong/mouseproject"


software="$folder/software"
savedata="$folder/save_data"
plotsvalidation="$folder/plots_validation"
resultsfolder="$folder/results050"
plotsfolder="$folder/plots050"
phenotypefolder="$folder/demonstration_trait/BT2"


test_chromosome=23
control_chromosome=19
input_snp_index=50


module load R/3.6.2
# Supplementary Figure 1
R --vanilla -f $software/dosage_association.R --args hap_address="$savedata/hap_table_Chr23.tsv" plots_folder=$plotsvalidation

# Figure 1
# Supplementary Figure 2
R --vanilla -f $software/stitch_validation.R --args dna_location="$savedata/stitch_genomewide.tsv" rna_location="$savedata/rnaseq.tsv" kinship_address="$savedata/Mjj_values_all_pnl.tsv" output_path=$plotsvalidation

# Figure 2
rm $resultsfolder/raw_LODs_all_scans.tsv
R --vanilla -f $software/summary_significance_threshold_qq.R --args test_save="$savedata/hapsnp_Chr23.Rsave" control_save="$savedata/hapsnp_Chr19.Rsave" results_folder=$resultsfolder plots_folder=$plotsfolder

# Supplementary Figure 3
R --vanilla -f $software/summary_plots.R --args input_snp_index=$input_snp_index chromosome=$control_chromosome results_folder=$resultsfolder plots_folder=$plotsfolder

### SUMMARY STATISTICS
R --vanilla -f $software/summary_statistics.R --args test_save="$savedata/hapsnp_Chr23.Rsave" control_save="$savedata/hapsnp_Chr19.Rsave" results_folder=$resultsfolder

# Figure 3
# Supplementary Figure 4
R --vanilla -f $software/summary_statistics_plots.R --args results_folder=$resultsfolder plots_folder=$plotsfolder

# Figure 4
R --vanilla -f $software/summary_plots.R --args input_snp_index=$input_snp_index chromosome=$test_chromosome results_folder=$resultsfolder plots_folder=$plotsfolder

module load R/4.0.2
# Figure 5
R --vanilla -f $software/demonstration_manhattan.R --args results_folder=$phenotypefolder plots_folder=$phenotypefolder shuffle_phenotype=FALSE 

# Supplementary Figure 5
R --vanilla -f $software/demonstration_manhattan.R --args results_folder="$phenotypefolder/control" plots_folder=$phenotypefolder shuffle_phenotype=TRUE

module load R/3.6.2
# Supplementary Figure 6
rm $phenotypefolder/raw_LODs_all_scans.tsv
R --vanilla -f $software/demonstration_qq.R --args demonstration_folder=$phenotypefolder control_folder="$phenotypefolder/control" results_folder=$phenotypefolder plots_folder=$phenotypefolder

# Figure 6
# Supplementary Figure 7
R --vanilla -f $software/demonstration_chromosome_scans.R --args results_folder=$phenotypefolder plots_folder=$phenotypefolder




# scp -r pnlong@hpc3.rcic.uci.edu:/share/adl/pnlong/mouseproject/results050/summary_statistics*.tsv ~/Desktop/mouseproject/results050
# scp -r pnlong@hpc3.rcic.uci.edu:/share/adl/pnlong/mouseproject/plots_validation/*.tiff ~/Desktop/mouseproject/plots_validation
# scp -r pnlong@hpc3.rcic.uci.edu:/share/adl/pnlong/mouseproject/plots050/*.tiff ~/Desktop/mouseproject/plots050
# scp -r pnlong@hpc3.rcic.uci.edu:/share/adl/pnlong/mouseproject/plots050/LOD_scores/*.tiff ~/Desktop/mouseproject/plots050/LOD_scores
# scp -r pnlong@hpc3.rcic.uci.edu:/share/adl/pnlong/mouseproject/plots050/correlation/*.tiff ~/Desktop/mouseproject/plots050/correlation
# scp -r pnlong@hpc3.rcic.uci.edu:/share/adl/pnlong/mouseproject/demonstration_trait/BT2/*.tiff ~/Desktop/mouseproject/demonstration_trait/BT2

