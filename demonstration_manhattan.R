# README
# This program will build a manhattan plot for the demonstration genome scans (both marker- and Haplotype).
# Phillip Long
# June 25, 2021
' .sh runner
####
module load R/4.0.2
# BT1
phenotypefolder="/share/adl/pnlong/mouseproject/demonstration_trait/BT1"
# BT2
phenotypefolder="/share/adl/pnlong/mouseproject/demonstration_trait/BT2"

srun -p standard -c 12 R --vanilla -f /share/adl/pnlong/mouseproject/software/demonstration_manhattan.R --args results_folder="$phenotypefolder/control" plots_folder=$phenotypefolder shuffle_phenotype=TRUE
srun -p standard -c 12 R --vanilla -f /share/adl/pnlong/mouseproject/software/demonstration_manhattan.R --args results_folder=$phenotypefolder plots_folder=$phenotypefolder shuffle_phenotype=FALSE 

####
'

# LOADS
#############
library(tidyverse)
library(magrittr)
# library(devtools)
# install_github("drveera/ggman")
library(ggman)
# install.packages("ggpubr")
library(ggpubr)

load(file = "/share/adl/pnlong/mouseproject/software/useful_functions.Rsave")
# FOR TESTING
# arguments <- c(results_folder="/share/adl/pnlong/mouseproject/demonstration_trait/BT1", plots_folder="/share/adl/pnlong/mouseproject/demonstration_trait/BT1", shuffle_phenotype="FALSE")
arguments <- get_args()
results_folder = as.character(arguments["results_folder"])
#############

manhattan_plotter <- function(test_type) {	
	# FILE LOCATIONS
	#############
	cat(paste("\n", test_type, sep = ""))
	
	files <- list.files(path = paste(results_folder, test_type, sep = "/"),
											recursive = TRUE, # recursive = TRUE lists every single file from every sub and sub-sub directory in results_folders
											all.files = FALSE, # eliminates any invisible files
											full.names = TRUE) %>% # returns absolute paths rather than relative paths
		str_sort(x = ., numeric = TRUE)
	
	#############
	
	# IMPORT FILES, CONCATENATE INTO ONE TABLE
	#############
	LOD_table <- tibble(
		CHROM = character(),
		TEST_SNP = numeric(),
		LOD_VALUE = numeric()
	)
	
	for (filename in files) {
		
		# import the file info
		file_info <- scan(file = filename, what = character(), nlines = 1, quiet = TRUE) %>% # read the first line of the file in
			extract(-1) %>% # since first value is a "#"
			strsplit(split = "=") %>%
			two_dimensional_vector_to_tibble(x = ., names = c("KEY", "VALUE")) %>%
			deframe()
		
		# simplify chromosome number to an index in the chromosomes list for better memory purposes
		test_snp_chromosome = as.character(file_info["tSNP_chromosome"]) %>%
			if_else(condition = (. == "1621"),
							true = "16+21",
							false = .)
		
		# import the data
		current_LOD_table <- read_tsv(file = filename,
																	col_names = TRUE,
																	col_types = list(col_double(), col_double(), col_double()),
																	comment = "#",
																	progress = FALSE) %>%
			correct_LOD_values(x = .) %>%
			filter(LOD_METHOD == 1) %>%
			mutate(CHROM = rep(test_snp_chromosome, times = nrow(.))) %>%
			select(CHROM, TEST_SNP, LOD_VALUE)
		
		# get the MSM and its position
		msm <- filter(current_LOD_table, LOD_VALUE == max(current_LOD_table$LOD_VALUE))
		cat(paste("\n**********************\n",
							test_snp_chromosome, ": ", msm$LOD_VALUE[[1]], "@", msm$TEST_SNP[[1]],
							"\n**********************\n"), sep = "")
		
		
		# bind current LOD table to main LOD table
		LOD_table <- bind_rows(LOD_table, current_LOD_table)
	}
	
	LOD_table %<>%
		mutate(TEST_SNP = TEST_SNP / 1e+6, # convert to Mb
					 CHROM = factor(x = .$CHROM,
					 							 levels = str_sort(x = unique(.$CHROM), numeric = TRUE))) %>% # convert chromosome indexes back to actual chromosome names
		arrange(CHROM, TEST_SNP)
	
	#############
	
	# MAKE THE MANATTAN PLOT
	#############
	sig_thresholds <- read_tsv(file = "/share/adl/pnlong/mouseproject/results050/summary_statistics.tsv",
														 col_names = TRUE,
														 col_types = cols(),
														 comment = "#") %>%
		select(calculation_method, sig_threshold) %>%
		distinct() %>%
		deframe()
	
	results_folder_suffix = str_remove(string = test_type, pattern = "based")
	
	sig_threshold = as.numeric(sig_thresholds[results_folder_suffix])
	
	if (results_folder_suffix == "snp") point_size = 0.01 else if (results_folder_suffix == "hap") point_size = 0.1
	
	cat("Making manhattan plot.\n")
	manhattan_plot <- ggman(gwas = LOD_table,
													snp = "TEST_SNP",
													bp = "TEST_SNP",
													chrom = "CHROM",
													pvalue = "LOD_VALUE",
													sigLine = sig_threshold,
													lineColour = "red",
													pointSize = point_size,
													logTransform = FALSE, # since already in LOD form
													xlabel = "Chromosome",
													ylabel = "LOD Score",
													legend.remove = TRUE,
													title = "")
	
	cat("Manhattan plot finished.\n**********************\n")
	
	remove(LOD_table)
	
	return(manhattan_plot)
	
}

snpmanhattan = manhattan_plotter("snpbased")
hapmanhattan = manhattan_plotter("hapbased")

manhattan <- ggarrange(snpmanhattan, hapmanhattan,
											 labels = c("Marker", "Haplotype"),
											 ncol = 1, nrow = 2,
											 align = "v")

ggsave(plot = manhattan,
			 dpi = 320,
			 path = arguments["plots_folder"], filename = if_else(condition = (as.logical(arguments["shuffle_phenotype"]) == FALSE),
			 																										 true = "manhattan.tiff",
			 																										 false = "manhattan_permuted.tiff"))

remove(snpmanhattan, hapmanhattan, manhattan, manhattan_plotter)
#############

# GET MSM DATA
#############
files <- c(list.files(path = paste(results_folder, "snpbased", sep = "/"),
										recursive = TRUE, # recursive = TRUE lists every single file from every sub and sub-sub directory in results_folders
										all.files = FALSE, # eliminates any invisible files
										full.names = TRUE), # returns absolute paths rather than relative paths
					 list.files(path = paste(results_folder, "hapbased", sep = "/"),
					 					 recursive = TRUE, # recursive = TRUE lists every single file from every sub and sub-sub directory in results_folders
					 					 all.files = FALSE, # eliminates any invisible files
					 					 full.names = TRUE)) %>%
	str_sort(x = ., numeric = TRUE)

MSM_table <- tibble(
	TEST_TYPE = character(),
	CHROM = character(),
	MSM = numeric(),
	MSM_LOD = numeric()
)

for (filename in files) {
	
	# import the file info
	file_info <- scan(file = filename, what = character(), nlines = 1, quiet = TRUE) %>% # read the first line of the file in
		extract(-1) %>% # since first value is a "#"
		strsplit(split = "=") %>%
		two_dimensional_vector_to_tibble(x = ., names = c("KEY", "VALUE")) %>%
		deframe()
	
	
	test_snp_chromosome = as.character(file_info["tSNP_chromosome"]) %>%
		if_else(condition = (. == "1621"),
						true = "16+21",
						false = .)
	test_type = as.character(file_info["calculation_method"]) %>%
		str_remove(string = ., pattern = "_based")
	
	
	# import the data
	current_LOD_table <- read_tsv(file = filename,
																col_names = TRUE,
																col_types = list(col_double(), col_double(), col_double()),
																comment = "#",
																progress = FALSE) %>%
		correct_LOD_values(x = .) %>%
		filter(LOD_METHOD == 1) %>%
		filter(LOD_VALUE == max(LOD_VALUE)) %>% # Get the MSM
		mutate(CHROM = rep(test_snp_chromosome, times = nrow(.)),
					 TEST_TYPE = rep(test_type, times = nrow(.))) %>%
		rename("MSM" = TEST_SNP,
					 "MSM_LOD" = LOD_VALUE) %>%
		select(TEST_TYPE, CHROM, MSM, MSM_LOD)
	
	# bind current LOD table to main LOD table
	MSM_table <- bind_rows(MSM_table, current_LOD_table)
}

write_tsv(x = MSM_table,
					path = paste(results_folder, "demonstration_msm.tsv", sep = "/"))
#############
