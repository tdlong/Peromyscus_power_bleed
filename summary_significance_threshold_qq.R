# README
# This program makes some histograms to determine a significance threshold
# Phillip Long
# November 26, 2020

' .sh runner
####
rm $resultsfolder/raw_LODs_all_scans.tsv
####
module load R/3.6.2
# 2.5% effect
resultsfolder="/share/adl/pnlong/mouseproject/results025"
plotsfolder="/share/adl/pnlong/mouseproject/plots025"
# 5.0% effect
resultsfolder="/share/adl/pnlong/mouseproject/results050"
plotsfolder="/share/adl/pnlong/mouseproject/plots050"

srun -p standard R --vanilla -f /share/adl/pnlong/mouseproject/software/summary_significance_threshold_qq.R --args test_save="/share/adl/pnlong/mouseproject/save_data/hapsnp_Chr23.Rsave" control_save="/share/adl/pnlong/mouseproject/save_data/hapsnp_Chr19.Rsave" results_folder=$resultsfolder plots_folder=$plotsfolder
####
'
# LOADS
#############
library(tidyverse)
library(magrittr)
load("/share/adl/pnlong/mouseproject/software/useful_functions.Rsave")
arguments <- get_args()
# for testing purposes:
# arguments <- c(test_save="/share/adl/pnlong/mouseproject/save_data/hapsnp_Chr23.Rsave", control_save="/share/adl/pnlong/mouseproject/save_data/hapsnp_Chr19.Rsave", results_folder="/share/adl/pnlong/mouseproject/results025", plots_folder="/share/adl/pnlong/mouseproject/plots025")
# arguments <- c(test_save="/share/adl/pnlong/mouseproject/save_data/hapsnp_Chr23.Rsave", control_save="/share/adl/pnlong/mouseproject/save_data/hapsnp_Chr19.Rsave", results_folder="/share/adl/pnlong/mouseproject/results050", plots_folder="/share/adl/pnlong/mouseproject/plots050")

files_to_keep_from_saves <- c("csnp_chromosome_number")

load(file = as.character(arguments["test_save"]),
		 envir = (Chr_test_env <- new.env())) # Chr23, store in its own environment
remove(list = ls(Chr_test_env)[!(ls(Chr_test_env) %in% files_to_keep_from_saves)], # get rid of objects not in list of objects we want to save
			 envir = Chr_test_env)

load(file = as.character(arguments["control_save"]),
		 envir = (Chr_control_env <- new.env())) # Chr19, store in its own environment
remove(list = ls(Chr_control_env)[!(ls(Chr_control_env) %in% files_to_keep_from_saves)], # get rid of objects not in list of objects we want to save
			 envir = Chr_control_env)

#############

# FILE LOCATIONS
#############
is_data_file <- function(files) {
	x <- vector(mode = "logical", length = length(files))
	
	for (i in seq_along(files)) {
		if ((length(unlist(strsplit(x = files[[i]], split = "\\/"))) == 8) & (endsWith(x = files[[i]], suffix = ".tsv"))) {
			x[[i]] = TRUE
		} else {
			x[[i]] = FALSE
		}
	}
	
	return(x)
}

results_folder = arguments["results_folder"]
files <- list.files(path = results_folder,
										recursive = TRUE, # recursive = TRUE lists every single file from every sub and sub-sub directory in results_folder
										all.files = FALSE, # eliminates any invisible files
										full.names = TRUE) %>% # returns absolute paths rather than relative paths
	str_sort(numeric = TRUE) %>%
	.[is_data_file(.)] # filter out so that there is the files we want
#############

# IMPORT FILES, CONCATENATE INTO ONE TABLE
#############
output_file = paste(results_folder, "raw_LODs_all_scans.tsv", sep = "/")

if (!file.exists(output_file)) {
	# write the header
	write_lines(x = paste("Calculation Method", "Test Type", paste("QTL on Chr. ", Chr_control_env$csnp_chromosome_number, sep = ""), paste("QTL on Chr. ", Chr_test_env$csnp_chromosome_number, sep = ""), sep = "\t"),
							path = output_file,
							sep = "\n",
							append = FALSE)
	
	test_types <- c("single", "multiple", "rare")
	calculation_methods <- c("snpbased", "hapbased")
	calc_model <- expand.grid(TEST_TYPE = test_types,
														CALCULATION_METHOD = calculation_methods) %>%
		tibble()
	
	
	make_quantile_percents <- function(base_interval_size = 0.001, final_interval_size) {
		quantile_vec <- seq(from = 0, to = 1, by = base_interval_size)
		
		interval_size = base_interval_size / 10
		while (interval_size >= final_interval_size) {
			# quantile_vec[[(length(quantile_vec) - 1)]]) is the second to last value
			current_vec <- seq(from = (quantile_vec[[(length(quantile_vec) - 1)]] + interval_size),
												 to = (1 - interval_size),
												 by = interval_size)
			
			quantile_vec <- c(quantile_vec, current_vec) %>%
				sort(decreasing = FALSE) %>%
				unique()
			
			interval_size = interval_size / 10
		}
		
		return(unique(quantile_vec))
	}
	
	current_calc_model_index = 0
	for (row_num in 1:nrow(calc_model)) {
		cat("CONTROL_MAXIMUM\tTEST_MAXIMUM\n")
		test_type = as.character(calc_model$TEST_TYPE[[row_num]])
		calculation_method = as.character(calc_model$CALCULATION_METHOD[[row_num]])
		current_files <- files[str_detect(string = files, pattern = test_type)
													 & str_detect(string = files, pattern = calculation_method)] %>% # should be 200 (100 per chromosome)
			str_sort(numeric = TRUE) %>%
			matrix(nrow = length(.) %/% 2, ncol = 2)
		
		current_files_line_count = wc_l(current_files[1, 1]) - 2
		
		final_interval = -8
		if (calculation_method == "hapbased") {
			current_files_line_count = current_files_line_count / 2 # only count the lines for LOD_METHOD == 1 (half of the lines are LOD_METHOD == 1)
			quantiles <- make_quantile_percents(base_interval_size = 0.001, final_interval_size = 10^(final_interval))
		} else if (calculation_method == "snpbased") {
			quantiles <- make_quantile_percents(base_interval_size = 0.001, final_interval_size = 10^(final_interval - 1)) # because 10x more values in SNPbased
		}
		
		final_length = current_files_line_count * (length(current_files) / 2) # /2 because of 2 chromosome values per row
		LOD_control <- vector(mode = "numeric", length = final_length)
		LOD_test <- vector(mode = "numeric", length = final_length)

		for (i in 1:nrow(current_files)) {
			
			filenames <- c("control" = current_files[i, 1], "test" = current_files[i, 2])
			
			current_data_control <- read_tsv(file = filenames["control"], col_names = TRUE, comment = "#", col_types = cols(), progress = FALSE) %>%
				correct_LOD_values(x = .) %>%
				filter(LOD_METHOD == 1) %>%
				pull(LOD_VALUE)
			
			current_data_test <- read_tsv(file = filenames["test"], col_names = TRUE, comment = "#", col_types = cols(), progress = FALSE) %>%
				correct_LOD_values(x = .) %>%
				filter(LOD_METHOD == 1) %>%
				pull(LOD_VALUE)
			
			cat(paste(max(current_data_control), "\t", max(current_data_test), sep = ""))
			if (max(current_data_control) >= max(current_data_test)) cat(" * * * * *")
			cat("\n") # print a newline
			
			
			LOD_control[(((i - 1) * current_files_line_count) + 1):(i * current_files_line_count)] <- current_data_control
			LOD_test[(((i - 1) * current_files_line_count) + 1):(i * current_files_line_count)] <- current_data_test
		}
		
		current_LOD_table <- tibble(
			CALCMETHOD = rep(x = calc_model$CALCULATION_METHOD[[row_num]], times = length(quantiles)),
			TESTTYPE = rep(x = calc_model$TEST_TYPE[[row_num]], times = length(quantiles)),
			CHRCONTROL = quantile(x = LOD_control, probs = quantiles),
			CHRTEST = quantile(x = LOD_test, probs = quantiles)
		)
		
		write_tsv(x = current_LOD_table,
							path = output_file,
							col_names = FALSE,
							append = TRUE)

		current_calc_model_index = current_calc_model_index + length(current_files)
		cat(paste("\n******************************************************************************\n",
							test_type, ", ", calculation_method, "\n",
							current_calc_model_index, " files scanned.",
							"\n******************************************************************************\n", sep = ""))
	}
	
}
#############

# MAKE QQPLOT
#############
LOD_table <- read_tsv(file = output_file, col_names = TRUE, progress = show_progress()) %>%
	mutate(`Calculation Method` = factor(x = .$`Calculation Method`, levels = c("hapbased", "snpbased")),
				 `Test Type` = factor(x = .$`Test Type`, levels = c("single", "multiple", "rare"))) # order with factors so that the faceting in ggplot is consistent with other plots

calculation_method <- c(hapbased = "Haplotype", snpbased = "Marker")
test_type <- c(single = "Single SNP", multiple = "10 SNPs", rare = "All SNPs")

qq_plot <- ggplot(data = LOD_table) +
	facet_grid(`Test Type` ~ `Calculation Method`, labeller = labeller(`Calculation Method` = calculation_method, `Test Type` = test_type)) +
	geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "red", alpha = .3) +
	geom_point(mapping = aes(x = `QTL on Chr. 19`, y = `QTL on Chr. 23`), size = 0.75, alpha = 0.25)

ggsave(plot = qq_plot,
			 dpi = 320,
			 path = arguments["plots_folder"],
			 filename = "qq_plot.tiff")
#############
