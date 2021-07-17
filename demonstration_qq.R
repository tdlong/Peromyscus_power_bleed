# README
# This program makes some histograms to determine a significance threshold
# Phillip Long
# November 26, 2020

' .sh runner
####
rm /share/adl/pnlong/mouseproject/demonstration_trait/BT1/raw_LODs_all_scans.tsv
rm /share/adl/pnlong/mouseproject/demonstration_trait/BT2/raw_LODs_all_scans.tsv
####
module load R/3.6.2
# BT1
phenotypefolder="/share/adl/pnlong/mouseproject/demonstration_trait/BT1"
# BT2
phenotypefolder="/share/adl/pnlong/mouseproject/demonstration_trait/BT2"

srun -p free R --vanilla -f /share/adl/pnlong/mouseproject/software/demonstration_qq.R --args demonstration_folder=$phenotypefolder control_folder="$phenotypefolder/control" results_folder=$phenotypefolder plots_folder=$phenotypefolder

####
'
# LOADS
#############
library(tidyverse)
library(magrittr)
load("/share/adl/pnlong/mouseproject/software/useful_functions.Rsave")
arguments <- get_args()
# for testing purposes:
# arguments <- c(demonstration_folder="/share/adl/pnlong/mouseproject/demonstration_trait/BT1", control_folder="/share/adl/pnlong/mouseproject/demonstration_trait/BT1/control", results_folder="/share/adl/pnlong/mouseproject/demonstration_trait/BT1", plots_folder="/share/adl/pnlong/mouseproject/demonstration_trait/BT1")
results_folder = arguments["results_folder"]
demonstration_folder = arguments["demonstration_folder"]
control_folder = arguments["control_folder"]
#############

# IMPORT FILES, CONCATENATE INTO ONE TABLE
#############
output_file = paste(results_folder, "raw_LODs_all_scans.tsv", sep = "/")

if (!file.exists(output_file)) {
	# WRITE THE HEADER
	#############
	write_lines(x = paste("IS_HAP", "DEMO", "CONTROL", sep = "\t"),
							path = output_file,
							sep = "\n",
							append = FALSE)
	#############
	
	
	# FILE LOCATIONS
	#############
	files_snpbased <- tibble(
		DEMO = str_sort(x = list.files(path = paste(demonstration_folder, "snpbased", sep = "/"), recursive = TRUE, all.files = FALSE, full.names = TRUE), numeric = TRUE),
		CONTROL = str_sort(x = list.files(path = paste(control_folder, "snpbased", sep = "/"), recursive = TRUE, all.files = FALSE, full.names = TRUE), numeric = TRUE)
	) %>%
		mutate(TEST_TYPE = rep("snp", times = nrow(.)))
	files_hapbased <- tibble(
		DEMO = str_sort(x = list.files(path = paste(demonstration_folder, "hapbased", sep = "/"), recursive = TRUE, all.files = FALSE, full.names = TRUE), numeric = TRUE),
		CONTROL = str_sort(x = list.files(path = paste(control_folder, "hapbased", sep = "/"), recursive = TRUE, all.files = FALSE, full.names = TRUE), numeric = TRUE)
	) %>%
		mutate(TEST_TYPE = rep("hap", times = nrow(.)))
	
	
	files <- bind_rows(files_snpbased, files_hapbased) %>%
		select(TEST_TYPE, DEMO, CONTROL) %>%
		mutate(IS_HAP = as.numeric(recode(.$TEST_TYPE, "snp" = FALSE, "hap" = TRUE)))
	#############
	
	quantiles_to_extract <- seq(from = 0.01, to = 0.99, by = 0.01)
	files_scanned_count = 0
	
	for (row_num in 1:nrow(files)) {
		
		current_files <- files[row_num, ]
		
		# LOAD IN LOD VALUES
		#############
		demo_LODs <- read_tsv(file = as.character(current_files$DEMO),
													comment = "#",
													col_names = TRUE,
													col_types = cols(col_double(), col_double(), col_double())) %>%
			correct_LOD_values(x = .) %>%
			filter(LOD_METHOD == 1) %>%
			pull(LOD_VALUE) %>%
			sort()
		control_LODs <- read_tsv(file = as.character(current_files$CONTROL),
														 comment = "#",
														 col_names = TRUE,
														 col_types = cols(col_double(), col_double(), col_double())) %>%
			correct_LOD_values(x = .) %>%
			filter(LOD_METHOD == 1) %>%
			pull(LOD_VALUE) %>%
			sort()
		
		files_scanned_count = files_scanned_count + 2
		#############
		
		# GET QUANTILES
		#############
		demo_values <- c(as.numeric(quantile(x = demo_LODs,
																				 probs = quantiles_to_extract,
																				 na.rm = TRUE)),
										 demo_LODs[demo_LODs >= as.numeric(quantile(x = demo_LODs, probs = 0.99, na.rm = TRUE))]
		)
		# extract quantiles 1-99%
		# then extract all values greater or equal to the 99th quantile
		
		control_values <- c(as.numeric(quantile(x = control_LODs,
																						probs = quantiles_to_extract,
																						na.rm = TRUE)),
												tail(x = control_LODs, n = (length(demo_values) - length(quantiles_to_extract)))
		)
		# this takes the top x largest LOD values in the controls,
		# x being the amount of values extracted from the demo_LODs that were in the 99th quantile
		
		#############
		
		# MAKE THE CURRENT TABLE TO BE WRITTEN TO THE OUTPUT FILE
		#############
		current_qq <- tibble(
			DEMO = demo_values,
			CONTROL = control_values
		) %>%
			mutate(IS_HAP = rep(x = as.numeric(current_files$IS_HAP), times = nrow(.))) %>%
			select(IS_HAP, DEMO, CONTROL)
		
		write_tsv(x = current_qq,
							path = output_file,
							col_names = FALSE,
							append = TRUE)
		#############
		
		cat(paste("\n******************************************************************************\n",
							files_scanned_count, " files scanned.",
							"\n******************************************************************************\n", sep = ""))
	}
	
}
#############

# MAKE HISTOGRAM
#############
LOD_table <- read_tsv(file = output_file,
											col_names = TRUE,
											progress = show_progress(),
											col_types = cols(col_character(), col_double(), col_double())) %>%
	group_by(IS_HAP) %>%
	mutate(IS_HAP = factor(x = IS_HAP, levels = c("1", "0")),
				 DEMO = sort(x = DEMO),
				 CONTROL = sort(x = CONTROL))

qq_plot <- ggplot(data = LOD_table) +
	facet_wrap(facets = vars(IS_HAP),
						 nrow = 1, ncol = 2,
						 labeller = labeller(IS_HAP = c("1" = "Haplotype", "0" = "Marker"))) +
	geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "red", alpha = .3) +
	geom_point(mapping = aes(x = CONTROL, y = DEMO), size = 0.75, alpha = 0.25) +
	labs(x = "Permuted Bleeding Time",
			 y = "Bleeding Time")

ggsave(plot = qq_plot,
			 dpi = 320,
			 path = arguments["plots_folder"],
			 filename = "qq_plot.tiff")
#############
