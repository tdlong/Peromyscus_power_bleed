# README
# This program makes some plots for summary statistics based off of
# /share/adl/pnlong/mouseproject/results
# Phillip Long
# February 20, 2021
' .sh runner
####
module load R/3.6.2
# 2.5% effect
resultsfolder="/share/adl/pnlong/mouseproject/results025"
plotsfolder="/share/adl/pnlong/mouseproject/plots025"
# 5.0% effect
resultsfolder="/share/adl/pnlong/mouseproject/results050"
plotsfolder="/share/adl/pnlong/mouseproject/plots050"

srun -p free R --vanilla -f /share/adl/pnlong/mouseproject/software/summary_statistics_plots.R --args results_folder=$resultsfolder plots_folder=$plotsfolder
####
'
# LOADS
#############
library(tidyverse)
library(magrittr)
load("/share/adl/pnlong/mouseproject/software/useful_functions.Rsave")
arguments <- get_args()
# arguments <- c(results_folder="/share/adl/pnlong/mouseproject/results050", plots_folder="/share/adl/pnlong/mouseproject/plots050")
# aruments <- c(results_folder="/share/adl/pnlong/mouseproject/results025", plots_folder="/share/adl/pnlong/mouseproject/plots025")
results_folder = arguments["results_folder"]
plots_folder = arguments["plots_folder"]
#############

# ADJUST SUMMARY TABLE
#############
summary_table_raw <- read_tsv(file = paste(results_folder, "summary_statistics.tsv", sep = "/"),
															col_names = TRUE)

tSNP_chromosome = unique(summary_table_raw$tSNP_chromosome)[[1]]
siglog <- summary_table_raw %>%
	select(calculation_method, sig_threshold) %>%
	distinct() %>%
	deframe()

summary_table <- summary_table_raw %>%
	arrange(desc(cSNP_chromosome), cSNP, test_type, calculation_method, sig_threshold) %>%
	mutate(cSNP_chromosome = as.factor(paste("Chr.", cSNP_chromosome, sep = " ")),
				 test_type = factor(x = recode(test_type, "single" = "Single SNP", "multiple" = "10 SNPs", "rare" = "All SNPs"), levels = c("Single SNP", "10 SNPs", "All SNPs")),
				 calculation_method = factor(x = recode(calculation_method, "snp" = "Marker", "hap" = "Haplotype"), levels = c("Haplotype", "Marker")), # order with factors so that the faceting in ggplot is consistent with other plots
				 sig_threshold = factor(x = paste(">", sig_threshold, sep = " "),
				 											 levels = paste(">", as.character(sort(siglog)), sep = " ")))

tSNP_chromosome = unique(summary_table$tSNP_chromosome)
#############

# SUMMARIZE THE SUMMARY STATISTICS
#############
percent_w_hits <- summary_table %>%
	group_by(cSNP_chromosome, test_type, calculation_method) %>%
	summarise(PERCENT_W_HITS = round(x = 100 * (sum(has_hit) / n()), digits = 0))

summary_table_summary <- summary_table %>%
	filter(has_hit == TRUE) %>%
	group_by(cSNP_chromosome, test_type, calculation_method) %>%
	summarise(AVG_HITS = round(x = mean(n_hits), digits = 0),
						AVG_DIST_MSM_MB = round(x = mean(avgdist_MSM) / 1e+06, digits = 2),
						MED_DIST_MSM_MB = round(x = median(avgdist_MSM) / 1e+06, digits = 2),
						AVG_SPAN_MB = round(x = mean(hit_span, na.rm = TRUE) / 1e+06, digits = 2),
						AVG_LOD_MSM = round(x = mean(MSM_LOD), digits = 2)) %>%
	left_join(x = ., y = percent_w_hits, by = c("cSNP_chromosome" = "cSNP_chromosome", "test_type" = "test_type", "calculation_method" = "calculation_method")) %>%
	rename_all(toupper) %>%
	select(CSNP_CHROMOSOME, TEST_TYPE, CALCULATION_METHOD, PERCENT_W_HITS, AVG_HITS, AVG_SPAN_MB, AVG_DIST_MSM_MB, MED_DIST_MSM_MB, AVG_LOD_MSM) %>%
	arrange(TEST_TYPE)

write_tsv(x = summary_table_summary,
					path = paste(results_folder, "summary_statistics_summary.tsv", sep = "/"),
					col_names = TRUE)



percent_w_hits_summary <- summary_table %>%
	group_by(cSNP_chromosome, calculation_method) %>%
	summarise(POWER = round(x = 100 * (sum(has_hit) / n()), digits = 0))

summary_table_summary_summary <- summary_table %>%
	filter(has_hit == TRUE) %>%
	group_by(cSNP_chromosome, calculation_method) %>%
	summarise(MEAN_HIT_COUNT = round(x = mean(n_hits), digits = 0),
						LOCALIZATION_MSM_AVGDIST = round(x = mean(avgdist_MSM) / 1e+06, digits = 2),
						LOCALIZATION_MSM_MEDDIST = round(x = median(avgdist_MSM) / 1e+06, digits = 2),
						LOCALIZATION_SPAN = round(x = mean(hit_span, na.rm = TRUE) / 1e+06, digits = 2),
						MEAN_MSM_LOD = round(x = mean(MSM_LOD), digits = 2)) %>%
	left_join(x = ., y = percent_w_hits_summary, by = c("cSNP_chromosome" = "cSNP_chromosome", "calculation_method" = "calculation_method")) %>%
	rename_all(toupper) %>%
	mutate(TEST_TYPE = rep("Average", times = nrow(.))) %>%
	select(CSNP_CHROMOSOME, TEST_TYPE, CALCULATION_METHOD, POWER, MEAN_HIT_COUNT, LOCALIZATION_SPAN, LOCALIZATION_MSM_AVGDIST, LOCALIZATION_MSM_MEDDIST, MEAN_MSM_LOD) %>%
	arrange(CALCULATION_METHOD)

write_tsv(x = summary_table_summary_summary,
					path = paste(results_folder, "summary_statistics_summary_summary.tsv", sep = "/"),
					col_names = TRUE)


# make the table into an image
# library(knitr)
# library(kableExtra)
# library(webshot)
# summary_table_summary %>%
# 	filter(CSNP_CHROMOSOME == tSNP_chromosome) %>% # so just Chromosome 23, because we know the statistics of Chromosome 19 Scans
# 	select(CSNP_CHROMOSOME, TEST_TYPE, CALCULATION_METHOD, PERCENT_W_HITS, AVG_HITS, AVG_SPAN, AVG_DIST_MSM, AVG_LOD_MSM) %>%
# 	rename("cSNP Chromosome" = CSNP_CHROMOSOME,
# 				 "Model" = TEST_TYPE,
# 				 "Test Type" = CALCULATION_METHOD,
# 				 "Power" = PERCENT_W_HITS,
# 				 "Mean 'hit' count" = AVG_HITS,
# 				 "Mean span" = AVG_SPAN,
# 				 "Mean MSM distance" = AVG_DIST_MSM,
# 				 "Mean MSM LOD Score" = AVG_LOD_MSM) %>%
# 	knitr::kable(x = ., padding = 5, digits = 3) %>% # create kable table from summary_table_summary
# 	kable_styling() %>%
# 	save_kable(x = ., file = paste(plots_folder, "summary_table_summary.png", sep = "/"))
#############

# LOD of MSM versus distance of MSM from cSNP
#############
msmlod_msmdist_plot <- ggplot(data = summary_table %>%
																filter(n_hits > 0) %>%
																mutate(avgdist_MSM = avgdist_MSM / 1e+06),
															mapping = aes(x = MSM_LOD, y = avgdist_MSM)) +
	facet_grid(~ calculation_method) +
	geom_point(mapping = aes(color = test_type)) +
	labs(x = "LOD Score of MSM",
			 y = "Distance of MSM from Causative Gene (Mb)") +
	scale_color_discrete(name = "Model")
	
ggsave(plot = msmlod_msmdist_plot,
			 dpi = 320,
			 path = plots_folder,
			 filename = "msmlod_msmdist.tiff")
#############

# OUTLIER MSM PLOTS
#############
outlier_files <- summary_table_raw %>%
	filter(n_hits > 0 & avgdist_MSM / 1e+06 > 7.5) %>%
	mutate(cSNP_chromosome = paste("Chr", .$cSNP_chromosome, sep = ""),
				 calculation_method = paste(.$calculation_method, "based", sep = ""),
				 ADDRESS = paste(results_folder, paste(cSNP_chromosome, test_type, calculation_method, sep = "_"), paste(paste(cSNP_chromosome, test_type, calculation_method, sep = "_"), snp_index, "tsv", sep = "."), sep = "/")) %>%
	pull(ADDRESS) %>%
	str_sort(x = ., numeric = TRUE)

LOD_table <- tibble(
	TEST_SNP = numeric(),
	LOD_VALUE = numeric(),
	CAUSATIVE_SNP = numeric(),
	CAUSATIVE_SNP_CHROM = character(),
	METHOD = character(),
	TEST_TYPE = character(),
	TEST_SNP_CHROM = character()
)

for (filename in outlier_files) {
	
	# import the file info
	file_info <- scan(file = filename, what = character(), nlines = 1, quiet = TRUE) %>% # read the first line of the file in
		extract(-1) %>% # since first value is a "#"
		strsplit(split = "=") %>%
		two_dimensional_vector_to_tibble(x = ., names = c("KEY", "VALUE")) %>%
		deframe()
	
	# import the data
	current_LOD_table <- read_tsv(file = filename,
																col_names = TRUE,
																col_types = list(col_double(), col_character(), col_double()),
																comment = "#",
																progress = FALSE) %>%
		correct_LOD_values(x = .) %>%
		filter(LOD_METHOD == 1) %>%
		select(-LOD_METHOD)
	
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
	
	z = nrow(current_LOD_table)
	# add type column
	current_LOD_table %<>%
		mutate(CAUSATIVE_SNP = rep(cSNP_locus, z),
					 CAUSATIVE_SNP_CHROM = rep(cSNP_chromosome, z),
					 METHOD = rep(calculation_method, z),
					 TEST_TYPE = rep(test_type, z),
					 TEST_SNP_CHROM = rep(test_snp_chromosome, z))
	
	# bind current LOD table to main LOD table
	LOD_table <- bind_rows(LOD_table, current_LOD_table)
}

tsnp_chromosome = unique(LOD_table$TEST_SNP_CHROM)[[1]] # in case theres more than one for some wierd reason...

LOD_table <- LOD_table %>%
	arrange(CAUSATIVE_SNP_CHROM, METHOD, TEST_TYPE) %>%
	mutate("Calculation Method-Test Type-Causative SNP" = as.factor(paste(.$METHOD, "-", .$TEST_TYPE, "\n", .$CAUSATIVE_SNP, sep = "")),
				 TEST_SNP = TEST_SNP / 1e+6,
				 CAUSATIVE_SNP = CAUSATIVE_SNP / 1e+6) %>% # convert both the TESTING and CAUSATIVE SNP values from bases to Mb (megabases)
	mutate_if(is.character, as.factor) %>%
	select(-TEST_SNP_CHROM, TEST_TYPE, CAUSATIVE_SNP, everything()) %>%
	rename("Locus (Mb)" = TEST_SNP,
				 "LOD Score" = LOD_VALUE,
				 "Causative SNP" = CAUSATIVE_SNP,
				 "Causative SNP Chromosome" = CAUSATIVE_SNP_CHROM,
				 "Calculation Method" = METHOD,
				 "Test Type" = TEST_TYPE)

# plotting the csnp line based on if its the right csnp chromsome or not
csnp_lines <- LOD_table %>%
	select(`Causative SNP Chromosome`, `Calculation Method-Test Type-Causative SNP`, `Causative SNP`) %>%
	distinct()

# plotting the significance threshold differently depending on the calculation method
sig_thresholds <- summary_table %>%
	select(calculation_method, sig_threshold) %>%
	distinct() %>%
	mutate(sig_threshold = as.numeric(str_remove(string = sig_threshold, pattern = "> "))) %>%
	deframe()

sig_threshold_lines <- expand.grid(`Calculation Method` = unique(LOD_table$`Calculation Method`),
																	 `Causative SNP Chromosome` = unique(LOD_table$`Causative SNP Chromosome`), # the chromosome of the causative SNP
																	 `Calculation Method-Test Type-Causative SNP` = unique(LOD_table$`Calculation Method-Test Type-Causative SNP`)) %>%
	mutate(sig_threshold = recode(as.character(.$`Calculation Method`), !!!sig_thresholds)) %>%
	as_tibble() %>%
	select(-`Calculation Method`)


# shading of the points
alpha_value_normal = 0.01
# order with faceting variables factors so that the faceting in ggplot is consistent with other plots
LOD_table %<>%
	mutate(alpha_shades = recode(.x = as.character(.$`Calculation Method`),
				 											"Haplotype" = (alpha_value_normal * 10), # because there are 10 times more values in the Marker, make Haplotype 10 times darker so the shades match
				 											"Marker" = alpha_value_normal))

manhattan_plot <- ggplot(data = LOD_table) +
	facet_wrap(facets = vars(`Calculation Method-Test Type-Causative SNP`),
						 nrow = 2, ncol = 2) +
	geom_vline(data = csnp_lines, # note the different data argument
						 mapping = aes(xintercept = `Causative SNP`), # and different mapping argument
						 color = "blue", size = .5, alpha = .6, linetype = "solid") +
	geom_hline(data = sig_threshold_lines, # note the different data argument, again
						 mapping = aes(yintercept = sig_threshold), # and different mapping argument, again
						 color = "red", size = .5, alpha = .6, linetype = "solid") +
	geom_point(data = LOD_table,
						 mapping = aes(x = `Locus (Mb)`, y = `LOD Score`, alpha = alpha_shades),
						 size = .1, na.rm = TRUE) +
	scale_alpha_continuous(guide = FALSE)

ggsave(plot = manhattan_plot,
			 dpi = 320,
			 filename = "outlier_msmdist.tiff",
			 path = plots_folder)

#############
