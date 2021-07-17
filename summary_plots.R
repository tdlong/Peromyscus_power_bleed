# README
# This program makes an image with a plot for every test type in it, when given a causative SNP index.
# /share/adl/pnlong/mouseproject/results/LOD_plots
# Phillip Long
# November 26, 2020

' .sh runner
####
input_snp_index=50

test_chromosome=23
control_chromosome=19

module load R/3.6.2

# 2.5% effect
resultsfolder="/share/adl/pnlong/mouseproject/results025"
plotsfolder="/share/adl/pnlong/mouseproject/plots025"
srun -p standard R --vanilla -f /share/adl/pnlong/mouseproject/software/summary_plots.R --args input_snp_index=$input_snp_index chromosome=$test_chromosome results_folder=$resultsfolder plots_folder=$plotsfolder
srun -p standard R --vanilla -f /share/adl/pnlong/mouseproject/software/summary_plots.R --args input_snp_index=$input_snp_index chromosome=$control_chromosome results_folder=$resultsfolder plots_folder=$plotsfolder
	
# 5.0% effect
resultsfolder="/share/adl/pnlong/mouseproject/results050"
plotsfolder="/share/adl/pnlong/mouseproject/plots050"
srun -p standard R --vanilla -f /share/adl/pnlong/mouseproject/software/summary_plots.R --args input_snp_index=$input_snp_index chromosome=$test_chromosome results_folder=$resultsfolder plots_folder=$plotsfolder
srun -p standard R --vanilla -f /share/adl/pnlong/mouseproject/software/summary_plots.R --args input_snp_index=$input_snp_index chromosome=$control_chromosome results_folder=$resultsfolder plots_folder=$plotsfolder
####
'

# LOADS
#############
library(tidyverse)
library(magrittr)
load("/share/adl/pnlong/mouseproject/software/useful_functions.Rsave")
arguments <- get_args()
# for testing purposes:
# arguments <- c(input_snp_index="50", chromosome="23", results_folder="/share/adl/pnlong/mouseproject/results050", plots_folder="/share/adl/pnlong/mouseproject/plots050")
# arguments <- c(input_snp_index="50", chromosome="19", results_folder="/share/adl/pnlong/mouseproject/results050", plots_folder="/share/adl/pnlong/mouseproject/plots050")

#############

# FILE LOCATIONS
#############
results_folder = as.character(arguments["results_folder"])
csnp_chromosome = as.character(arguments["chromosome"])
all_files <- list.files(path = results_folder,
												recursive = TRUE, # recursive = TRUE lists every single file from every sub and sub-sub directory in results_folders
												all.files = FALSE, # eliminates any invisible files
												full.names = TRUE) # returns absolute paths rather than relative paths

input_snp_index = as.numeric(arguments["input_snp_index"])

files <- all_files[endsWith(x = all_files, suffix = paste(".", input_snp_index, ".tsv", sep = "")) & str_detect(string = all_files, pattern = paste("Chr", csnp_chromosome, sep = ""))]
remove(all_files)
#############

# IMPORT FILES, CONCATENATE INTO ONE TABLE
#############
LOD_table <- tibble(
	TEST_SNP = numeric(),
	LOD_METHOD = character(),
	LOD_VALUE = numeric(),
	CAUSATIVE_SNP = numeric(),
	CAUSATIVE_SNP_CHROM = character(),
	METHOD = character(),
	TEST_TYPE = character(),
	TEST_SNP_CHROM = character()
)

LOD_table_1 <- tibble(
	TEST_SNP = numeric(),
	LOD_VALUE = numeric(),
	CAUSATIVE_SNP = numeric(),
	CAUSATIVE_SNP_CHROM = character(),
	METHOD = character(),
	TEST_TYPE = character(),
	TEST_SNP_CHROM = character()
)

for (filename in files) {
	
	# import the file info
	file_info <- scan(file = filename, what = character(), nlines = 1, quiet = TRUE) %>% # read the first line of the file in
		extract(-1) %>% # since first value is a "#"
		strsplit(split = "=") %>%
		two_dimensional_vector_to_tibble(x = ., names = c("KEY", "VALUE")) %>%
		deframe()
	
	# the causative snp
	cSNP_locus = as.numeric(file_info["cSNP_locus"])
	
	# calculation_method + cSNP_chromosome
	cSNP_chromosome = as.character(file_info["cSNP_chromosome"]) %>%
		paste("Chr.", ., sep = " ")
	
	calculation_method_raw = as.character(file_info["calculation_method"])
	calculation_method = case_when( # nicely format
		calculation_method_raw == "hap_based" ~ "Haplotype",
		calculation_method_raw == "snp_based" ~ "Marker"
	)
	
	# test_type
	test_type_raw = as.character(file_info["test_type"])
	test_type = case_when(# nicely format
		test_type_raw == "single" ~ "Single SNP",
		test_type_raw == "multiple" ~ "10 SNPs",
		test_type_raw == "rare" ~ "All SNPs"
	)
	
	test_snp_chromosome = as.character(file_info["tSNP_chromosome"])
	
	# add type column
	current_LOD_table <- read_tsv(file = filename,
																col_names = TRUE,
																col_types = list(col_double(), col_character(), col_double()),
																comment = "#",
																progress = FALSE) %>%
		mutate(LOD_METHOD = recode(LOD_METHOD, "1" = "Approximation Method", "2" = "Exact Method")) %>%
		mutate(CAUSATIVE_SNP = rep(cSNP_locus, nrow(.)),
					 CAUSATIVE_SNP_CHROM = rep(cSNP_chromosome, nrow(.)),
					 METHOD = rep(calculation_method, nrow(.)),
					 TEST_TYPE = rep(test_type, nrow(.)),
					 TEST_SNP_CHROM = rep(test_snp_chromosome, nrow(.)))
	
	# bind current LOD table to main LOD table
	LOD_table <- bind_rows(LOD_table, current_LOD_table)
	
	
	
	current_LOD_table_1 <- read_tsv(file = filename,
																	col_names = TRUE,
																	col_types = list(col_double(), col_character(), col_double()),
																	comment = "#",
																	progress = FALSE) %>%
		correct_LOD_values(x = .) %>%
		filter(LOD_METHOD == 1) %>%
		select(-LOD_METHOD) %>%
		mutate(CAUSATIVE_SNP = rep(cSNP_locus, nrow(.)),
					 CAUSATIVE_SNP_CHROM = rep(cSNP_chromosome, nrow(.)),
					 METHOD = rep(calculation_method, nrow(.)),
					 TEST_TYPE = rep(test_type, nrow(.)),
					 TEST_SNP_CHROM = rep(test_snp_chromosome, nrow(.)))
	
	LOD_table_1 <- bind_rows(LOD_table_1, current_LOD_table_1)
}

tsnp_chromosome = unique(LOD_table$TEST_SNP_CHROM)[[1]] # in case theres more than one for some wierd reason...

LOD_table <- LOD_table %>%
	arrange(CAUSATIVE_SNP_CHROM, METHOD, TEST_TYPE, LOD_METHOD) %>%
	mutate_if(is.character, as.factor) %>%
	mutate(TEST_SNP = TEST_SNP / 1e+6,
				 CAUSATIVE_SNP = CAUSATIVE_SNP / 1e+6) %>% # convert both the TESTING and CAUSATIVE SNP values from bases to Mb (megabases)
	select(METHOD, TEST_TYPE, CAUSATIVE_SNP, everything()) %>%
	rename("Locus (Mb)" = TEST_SNP,
				 "LOD Calculation Method" = LOD_METHOD,
				 "LOD Score" = LOD_VALUE,
				 "Causative SNP" = CAUSATIVE_SNP,
				 "Causative SNP Chromosome" = CAUSATIVE_SNP_CHROM,
				 "Calculation Method" = METHOD,
				 "Test Type" = TEST_TYPE)

LOD_table_1 <- LOD_table_1 %>%
	arrange(CAUSATIVE_SNP_CHROM, METHOD, TEST_TYPE) %>%
	mutate_if(is.character, as.factor) %>%
	mutate(TEST_SNP = TEST_SNP / 1e+6,
				 CAUSATIVE_SNP = CAUSATIVE_SNP / 1e+6) %>% # convert both the TESTING and CAUSATIVE SNP values from bases to Mb (megabases)
	select(METHOD, TEST_TYPE, CAUSATIVE_SNP, everything()) %>%
	rename("Locus (Mb)" = TEST_SNP,
				 "LOD Score" = LOD_VALUE,
				 "Causative SNP" = CAUSATIVE_SNP,
				 "Causative SNP Chromosome" = CAUSATIVE_SNP_CHROM,
				 "Calculation Method" = METHOD,
				 "Test Type" = TEST_TYPE)
	
causative_snp = filter(LOD_table, `Causative SNP Chromosome` == "Chr. 23") %>%
	pull(`Causative SNP`) %>%
	extract(1)
#############

# BUILD THE METHOD 1 PLOT
#############
# plotting the csnp line based on if its the right csnp chromsome or not
csnp_lines <- expand.grid(`Test Type` = unique(LOD_table_1$`Test Type`), # chooses all of the test_type values
						`Causative SNP Chromosome` = paste("Chr. ", tsnp_chromosome, sep = ""), # the chromosome of the causative SNP
						`Calculation Method` = unique(LOD_table_1$`Calculation Method`),
						cSNP = causative_snp)

# plotting the significance threshold differently depending on the calculation method
sig_thresholds <- read_tsv(file = paste(results_folder, "summary_statistics.tsv", sep = "/"), col_names = TRUE, col_types = cols(), comment = "#") %>%
	select(calculation_method, sig_threshold) %>%
	mutate(calculation_method = recode(.$calculation_method, "hap" = "Haplotype", "snp" = "Marker")) %>%
	distinct() %>%
	deframe()
	
sig_threshold_lines <- expand.grid(`Test Type` = unique(LOD_table_1$`Test Type`), # chooses all of the test_type values
																	 `Causative SNP Chromosome` = unique(LOD_table_1$`Causative SNP Chromosome`), # the chromosome of the causative SNP
																	 `Calculation Method` = unique(LOD_table_1$`Calculation Method`)) %>%
	mutate(sig_threshold = recode(.$`Calculation Method`, !!!sig_thresholds))


# shading of the points
alpha_value_normal = 0.01
# order with faceting variables factors so that the faceting in ggplot is consistent with other plots
LOD_table_1 %<>%
	mutate(ALPHA_SHADES = recode(.x = .$`Calculation Method`,
															 "Haplotype" = (alpha_value_normal * 10), # because there are 10 times more values in the Marker, make Haplotype 10 times darker so the shades match
															 "Marker" = alpha_value_normal),
				 `Calculation Method` = factor(x = .$`Calculation Method`, levels = c("Haplotype", "Marker")),
				 `Test Type` = factor(x = .$`Test Type`, levels = c("Single SNP", "10 SNPs", "All SNPs")))

# # make the caption
# main_subtitle = paste("Causative Gene from Chromosome ", csnp_chromosome,
# 											(if(csnp_chromosome == tsnp_chromosome) paste(": ", (causative_snp * 1e+6), sep = "")), # if I am on chromosome 23
# 											sep = "")

main_plot <- ggplot(data = LOD_table_1) +
	facet_grid(`Test Type` ~ `Calculation Method`) +
	geom_vline(data = csnp_lines, # note the different data argument
						 mapping = aes(xintercept = cSNP), # and different mapping argument
						 color = "blue", size = .5, alpha = .6, linetype = "solid") +
	geom_hline(data = sig_threshold_lines, # note the different data argument, again
						 mapping = aes(yintercept = sig_threshold), # and different mapping argument, again
						 color = "red", size = .5, alpha = .6, linetype = "solid") +
	geom_point(data = LOD_table_1,
						 mapping = aes(x = `Locus (Mb)`, y = `LOD Score`, alpha = ALPHA_SHADES),
						 size = .1, na.rm = TRUE) +
	scale_alpha_continuous(guide = FALSE)

# save plot as tiff
ggsave(plot = main_plot,
			 dpi = 320,
			 filename = paste("Chr", csnp_chromosome, "_LOD_approximation.", input_snp_index, ".tiff", sep = ""),
			 path = paste(as.character(arguments["plots_folder"]), "LOD_scores", sep = "/"))
#############

# BUILD METHODS CORRELATION PLOT
#############
LOD_table_corr <- LOD_table %>%
	filter(`Calculation Method` == "Haplotype") %>% # since Haplotype are the only type where there are different LODs
	pivot_wider(names_from = `LOD Calculation Method`, values_from = `LOD Score`)

corr_plot <- ggplot(data = LOD_table_corr) +
	geom_abline(slope = 1, intercept = 0, color = "red") +
	geom_point(mapping = aes(x = `Approximation Method`, y = `Exact Method`), alpha = 0.5, na.rm = TRUE)
# labs(title = "LOD Score Calculation Method, Correlation Comparison",
# 		 subtitle = "Approximation versus Exact",
# 		 caption = "The exact LOD score calculation is only performed if the approximation method first returns\na score greater than 2. Hence, there are no exact values if the approximation value is less than 2.")

# save plot as tiff
ggsave(plot = corr_plot,
			 dpi = 320,
			 filename = paste("Chr", csnp_chromosome, "_LOD_correlation.", input_snp_index, ".tiff", sep = ""),
			 path = paste(as.character(arguments["plots_folder"]), "correlation", sep = "/"))
#############

