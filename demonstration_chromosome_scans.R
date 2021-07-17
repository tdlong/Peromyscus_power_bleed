# README
# This program will build position vs. LOD score plots for individual chromosomes from the demonstration scan with notable points.
# Phillip Long
# June 29, 2021
' .sh runner
####
module load R/3.6.2
# BT1
phenotypefolder="/share/adl/pnlong/mouseproject/demonstration_trait/BT1"
# BT2
phenotypefolder="/share/adl/pnlong/mouseproject/demonstration_trait/BT2"

srun -p standard -c 2 R --vanilla -f /share/adl/pnlong/mouseproject/software/demonstration_chromosome_scans.R --args results_folder=$phenotypefolder plots_folder=$phenotypefolder
####
'

# LOADS
#############
library(tidyverse)
library(magrittr)

load(file = "/share/adl/pnlong/mouseproject/software/useful_functions.Rsave")
# FOR TESTING
# arguments <- c(results_folder="/share/adl/pnlong/mouseproject/demonstration_trait/BT1", plots_folder="/share/adl/pnlong/mouseproject/demonstration_trait/BT1")
arguments <- get_args()
results_folder = as.character(arguments["results_folder"])
#############

# GET THE INFORMATION OF THE SCANS WE WANT TO LOOK AT
#############
chromosomes <- c("3", "13", "15", "22")
files_info <- chromosomes %>%
	expand.grid(TEST_TYPE = c("snp", "hap"),
							CHROM = .) %>%
	as_tibble() %>%
	mutate(TEST_TYPE_FORMAL = factor(x = recode(TEST_TYPE, "snp" = "Marker", "hap" = "Haplotype"),
																	 levels = c("Marker", "Haplotype")),
				 CHROMOSOME_FORMAL = factor(x = paste("Chr", CHROM, sep = ""),
				 													 levels = paste("Chr", str_sort(x = unique(.$CHROM), numeric = TRUE), sep = "")),
				 FILES = paste(results_folder, "/", TEST_TYPE, "based/Chr", CHROM, "_", TEST_TYPE, "based.tsv", sep = "")) %>%
	arrange(TEST_TYPE_FORMAL, CHROMOSOME_FORMAL)

files <- files_info$FILES

#############

# IMPORT FILES, CONCATENATE INTO ONE TABLE
#############
LOD_table <- tibble(
	TEST_TYPE = character(),
	CHROM = character(),
	TEST_SNP = numeric(),
	LOD_VALUE = numeric()
)

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
		mutate(CHROM = rep(test_snp_chromosome, times = nrow(.)),
					 TEST_TYPE = rep(test_type, times = nrow(.)))
	
	# bind current LOD table to main LOD table
	LOD_table <- bind_rows(LOD_table,
												 (current_LOD_table %>%
												  	select(TEST_TYPE, CHROM, TEST_SNP, LOD_VALUE)))
	
	MSM_table <- bind_rows(MSM_table,
												 (current_LOD_table %>%
												  	filter(LOD_VALUE == max(LOD_VALUE)) %>% # Get the MSM
												  	rename("MSM" = TEST_SNP,
												  				 "MSM_LOD" = LOD_VALUE) %>%
												  	select(TEST_TYPE, CHROM, MSM, MSM_LOD)))
	
}

LOD_table %<>%
	mutate(TEST_SNP = TEST_SNP / 1e+6, # convert to Mb
				 CHROM = factor(x = paste("Chr", .$CHROM, sep = ""),
				 							 levels = paste("Chr", str_sort(x = unique(.$CHROM), numeric = TRUE), sep = "")),
				 TEST_TYPE = factor(x = recode(TEST_TYPE, "snp" = "Marker", "hap" = "Haplotype"),
				 									 levels = c("Marker", "Haplotype"))) %>% # convert chromosome indexes back to actual chromosome names
	arrange(TEST_TYPE, CHROM, TEST_SNP)


#############

# MAKE THE FULL CHROMOSOMES PLOT
#############
sig_thresholds <- read_tsv(file = "/share/adl/pnlong/mouseproject/results050/summary_statistics.tsv",
													 col_names = TRUE,
													 col_types = cols(),
													 comment = "#") %>%
	select(calculation_method, sig_threshold) %>%
	distinct() %>%
	deframe()

sig_threshold_lines <- expand.grid(TEST_TYPE = unique(LOD_table$TEST_TYPE), # chooses all of the test_type values
																	 CHROM = unique(LOD_table$CHROM)) %>% # the chromosome of the causative SNP
	mutate(sig_threshold = recode(as.character(.$TEST_TYPE), "Marker" = sig_thresholds["snp"], "Haplotype" = sig_thresholds["hap"]))



# shading of the points
alpha_value_normal = 0.05
LOD_table %<>%
	mutate(ALPHA_SHADES = recode(.x = LOD_table$TEST_TYPE,
															 "Haplotype" = (alpha_value_normal * 10), # because there are 10 times more values in the Marker-based, make Haplotype-based 10 times darker so the shades match
															 "Marker" = alpha_value_normal)) %>%
	left_join(x = ., y = sig_threshold_lines,
						by = c("TEST_TYPE" = "TEST_TYPE", "CHROM" = "CHROM")) %>%
	unite(TEST_TYPE, CHROM, col = "TEST_TYPE_CHROM", sep = "-", remove = FALSE) %>%
	mutate(TEST_TYPE_CHROM = factor(x = .$TEST_TYPE_CHROM,
																	levels = paste(files_info$TEST_TYPE_FORMAL, files_info$CHROMOSOME_FORMAL, sep = "-")))

chromosome_scans_plot <- ggplot(data = LOD_table) +
	facet_wrap(facets = vars(TEST_TYPE_CHROM),
						 nrow = 2, ncol = length(chromosomes),
						 drop = TRUE,
						 scales = "free_x") +
	geom_hline(mapping = aes(yintercept = sig_threshold), # and different mapping argument, again
						 color = "red", size = .5, alpha = .6) +
	geom_point(mapping = aes(x = TEST_SNP, y = LOD_VALUE, alpha = ALPHA_SHADES),
						 size = .1, na.rm = TRUE) +
	scale_alpha_continuous(guide = "none") +
	labs(x = "Locus (Mb)",
			 y = "LOD Score")

ggsave(plot = chromosome_scans_plot,
			 dpi = 320,
			 path = arguments["plots_folder"], filename = "chromosome_scans.tiff")
#############

# MAKE REGION AROUND MSM PLOTS
#############
# GET THE CENTERS (BY CHROMOSOME) FOR THE PLOTS
MSM_table %<>%
	left_join(x = ., y = enframe(x = sig_thresholds, name = "TEST_TYPE", value = "THRESHOLD"), by = c("TEST_TYPE")) %>%
	mutate(DIFF = c(MSM_LOD - THRESHOLD)) %>%
	arrange(desc(DIFF))

centering_scans_info <- MSM_table[match(x = chromosomes, table = MSM_table$CHROM), ] %>%
	select(TEST_TYPE, CHROM) %>%
	mutate(TEST_TYPE_FORMAL = recode(TEST_TYPE, "snp" = "Marker", "hap" = "Haplotype"),
				 CHROMOSOME_FORMAL = paste("Chr", CHROM, sep = ""))

LOD_table_1 <- subset(LOD_table, FALSE)
for (chrom in as.character(unique(LOD_table$CHROM))) {
	
	LOD_table_current <- filter(LOD_table, CHROM == chrom)
	
	# do we want to center around the snp-based or hap-based MSM?
	centering_scan <- filter(centering_scans_info, CHROMOSOME_FORMAL == chrom)
	
	MSM = LOD_table_current %>%
		filter(as.character(TEST_TYPE) == as.character(centering_scan$TEST_TYPE_FORMAL)) %>% # choose whether we want the MSM from snpbased or hapbased
		filter(LOD_VALUE == max(LOD_VALUE)) %>%
		pull(TEST_SNP) %>%
		extract(1)
	
	LOD_table_current <- filter(LOD_table_current, abs(TEST_SNP - MSM) <= 5) # So 10 Mb in total region around the MSM
	
	
	LOD_table_1 <- bind_rows(LOD_table_1, LOD_table_current)
}

chromosome_scans_msm_window_plot <- ggplot(data = LOD_table_1) +
	facet_wrap(facets = vars(TEST_TYPE_CHROM),
						 nrow = 2, ncol = length(chromosomes),
						 drop = TRUE,
						 scales = "free_x") +
	geom_hline(mapping = aes(yintercept = sig_threshold), # and different mapping argument, again
						 color = "red", size = .5, alpha = .6) +
	geom_point(mapping = aes(x = TEST_SNP, y = LOD_VALUE, alpha = ALPHA_SHADES),
						 size = .5, na.rm = TRUE) +
	scale_alpha_continuous(guide = "none") +
	scale_x_continuous(labels = scales::number_format(accuracy = 1)) +
	theme(axis.text.x = element_text(size = 6)) +
	labs(x = "Locus (Mb)",
			 y = "LOD Score")

ggsave(plot = chromosome_scans_msm_window_plot,
			 dpi = 320,
			 path = arguments["plots_folder"], filename = "chromosome_scans_msm_window.tiff")

# write MSMs from "region around" plots to .tsv
MSM_table <- LOD_table_1 %>%
	group_by(TEST_TYPE, CHROM) %>%
	summarise(LOD = max(LOD_VALUE),
						MSM = 1e+06 * as.numeric(filter(., LOD_VALUE == LOD) %>% pull(TEST_SNP))) %>%
	ungroup() %>%
	mutate_if(is.factor, as.character) %>%
	mutate(CHROM = as.numeric(str_remove(string = .$CHROM, pattern = "Chr"))) %>%
	select(TEST_TYPE, CHROM, MSM, LOD) %>%
	arrange(CHROM) %>%
	pivot_longer(cols = c("MSM", "LOD"), names_to = "TYPE", values_to = "VALUE") %>%
	unite(TYPE, TEST_TYPE, col = "TYPE", sep = "-", remove = TRUE) %>%
	pivot_wider(names_from = TYPE, values_from = VALUE) %>%
	rename("Chr" = CHROM) %>%
	select(Chr, `MSM-Marker`, `MSM-Haplotype`, `LOD-Marker`, `LOD-Haplotype`)

write_tsv(x = MSM_table,
					path = paste(results_folder, "demonstration_msm_from_windows.tsv", sep = "/"))
#############
