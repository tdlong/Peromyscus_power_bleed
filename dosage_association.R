# README
# This program determines dosage correlations for both the haplotype- and marker-based tests
# Phillip Long
# May 29, 2021
' .sh runner
####
module load R/3.6.2
srun -p standard -c 2 R --vanilla -f /share/adl/pnlong/mouseproject/software/dosage_association.R --args hap_address="/share/adl/pnlong/mouseproject/save_data/hap_table_Chr23.tsv" plots_folder="/share/adl/pnlong/mouseproject/plots_validation"
####
'

# LOADS
##############
library(tidyverse)
library(magrittr)
# install.packages("ggpubr")
library(ggpubr)
library(gridExtra)

load(file = "/share/adl/pnlong/mouseproject/software/useful_functions.Rsave")
arguments <- get_args()
# for testing purposes:
# arguments <- c(hap_address="/share/adl/pnlong/mouseproject/save_data/hap_table_Chr23.tsv", plots_folder="/share/adl/pnlong/mouseproject/plots_validation")
plots_folder = arguments["plots_folder"]
##############

# WRANGLE HAPLOTYPE AND SNP TABLES
##############
individuals <- read_lines(file = "/share/adl/pnlong/mouseproject/save_data/individuals.tsv", n_max = 1) %>% # read the first line of the haplotype table (the column names), since that includes the individual names
	strsplit(x = ., split = "\t") %>%
	unlist() %>%
	as.character()

coverages <- read_tsv(file = "/share/adl/pnlong/mouseproject/save_data/coverages.tsv", col_names = TRUE, col_types = cols(col_character(), col_number(), col_character())) %>%
	filter(TYPE == "gDNA" & INDV %in% individuals) %>%
	.[trunc(seq(from = 1, to = nrow(.), length.out = 6)), ] %>%
	select(INDV, MEAN_DEPTH) %>%
	arrange(MEAN_DEPTH) %>%
	mutate(MEAN_DEPTH = round(x = MEAN_DEPTH, digits = 1))

every_nth_line_pipe <- function(file, n = 10) {
	# begin at the second line (i = 1) because first line is column headers
	pipe_command = paste("tail -n +2 ", file, " | awk 'NR%", n, "==1'", sep = "")
	
	return(pipe(pipe_command))
}
##############

# QUANTILE PLOT
##############
# c("POS", unlist(map(.x = individuals, .f = ~ paste(paste("HAP", 1:8, sep = "_"), .x, sep = "-")))) generates column names for the hap table
hap_table_quantile <- read_tsv(file = every_nth_line_pipe(arguments["hap_address"], n = 10),
															 col_names = c("POS", unlist(map(.x = individuals, .f = ~ paste(paste("HAP", 1:8, sep = "_"), .x, sep = "-"))))) %>%
	select(POS, contains(coverages$INDV)) %>%
	pivot_longer(cols = names(.)[names(.) != "POS"], # all but the POS column
							 names_to = "DOSAGE_NUMBER-ID",
							 values_to = "DOSAGE_HOME") %>%
	separate(col = `DOSAGE_NUMBER-ID`, into = c("DOSAGE_NUMBER", "INDV"), sep = "-") %>%
	mutate(DOSAGE_NUMBER = as.numeric(str_extract(string = DOSAGE_NUMBER, pattern = "[0-9]+")),
				 INDV = factor(x = .$INDV, levels = coverages$INDV)) %>%
	group_by(INDV) %>%
	mutate(DOSAGE_NEIGHBOR = c(DOSAGE_HOME[-(1:8)], rep(x = NA_real_, times = 8))) %>% # the dosage of the neighboring locus
	ungroup() %>%
	filter(!is.na(DOSAGE_NEIGHBOR)) %>%
	mutate(ABS_DIFF = abs(DOSAGE_HOME - DOSAGE_NEIGHBOR)) %>%
	group_by(POS, INDV) %>%
	summarise(AVG_ABS_DIFF = sqrt(mean(ABS_DIFF ^ 2))) %>% # 0.50, 0.95,0.98,0.99,0.999
	group_by(INDV) %>%
	summarise("50%" = quantile(x = AVG_ABS_DIFF, probs = 0.50),
						"95%" = quantile(x = AVG_ABS_DIFF, probs = 0.95),
						"98%" = quantile(x = AVG_ABS_DIFF, probs = 0.98),
						"99%" = quantile(x = AVG_ABS_DIFF, probs = 0.99),
						"99.9%" = quantile(x = AVG_ABS_DIFF, probs = 0.999)) %>%
	pivot_longer(cols = c(`50%`, `95%`, `98%`, `99%`, `99.9%`), names_to = "QUANTILE", values_to = "QUANTILE_VALUE") %>%
	left_join(x = ., y = coverages, by = c("INDV" = "INDV")) %>%
	mutate(QUANTILE = factor(x = QUANTILE, levels = c("50%", "95%", "98%", "99%", "99.9%"))) %>%
	unite(INDV, MEAN_DEPTH, col = "INDV_COVERAGE", sep = " (") %>%
	mutate(INDV_COVERAGE = factor(x = paste(INDV_COVERAGE, "X)", sep = ""),
																levels = c(unite(data = coverages, INDV, MEAN_DEPTH, col = "INDV_COVERAGE", sep = " (") %>%
																	pull(INDV_COVERAGE) %>%
																	paste(., "X)", sep = ""))))

quantiles_plot <- ggplot(data = hap_table_quantile, mapping = aes(x = INDV_COVERAGE, y = QUANTILE_VALUE, color = QUANTILE, group = QUANTILE)) +
	geom_line() +
	geom_point() +
	scale_color_discrete(name = "Quantile") +
	scale_x_discrete(limits = hap_table_quantile$INDV_COVERAGE) +
	labs(x = "Individual",
			 y = "Mean Absolute Difference") +
	theme(axis.text.x = element_text(angle = 90))

ggsave(plot = quantiles_plot,
			 dpi = 320,
			 path = plots_folder, filename = "quantiles.tiff")
##############

# DOSAGE PLOT
##############
min_individual = coverages$INDV[[1]]
max_individual = coverages$INDV[[length(coverages$INDV)]]
point_size = 0.01
legend_point_size = 3
legend_name = "Haplotype Number"

hap_table_dosages <- read_tsv(file = every_nth_line_pipe(arguments["hap_address"], n = 10),
															col_names = c("POS", unlist(map(.x = individuals, .f = ~ paste(paste("HAP", 1:8, sep = "_"), .x, sep = "-"))))) %>%
	pivot_longer(cols = names(.)[names(.) != "POS"], # all but the POS column
							 names_to = "DOSAGE_NUMBER-ID",
							 values_to = "DOSAGE_HOME") %>%
	separate(col = `DOSAGE_NUMBER-ID`, into = c("DOSAGE_NUMBER", "INDV"), sep = "-") %>%
	mutate(POS = POS / 1e+06,
				 DOSAGE_NUMBER = factor(x = paste("#", str_remove(string = .$DOSAGE_NUMBER, pattern = "HAP_"), sep = ""), levels = paste("#", 1:8, sep = "")))


min_dosage_plot <- ggplot(data = filter(hap_table_dosages, INDV == min_individual)) +
	geom_point(mapping = aes(x = POS, y = DOSAGE_HOME, color = DOSAGE_NUMBER), size = point_size) +
	scale_color_discrete(name = legend_name) +
	guides(color = guide_legend(override.aes = list(size = legend_point_size))) +
	labs(x = "Locus (Mb)",
			 y = "Haplotype Dosage",
			 caption = paste("Individual #", min_individual,
			 								" (", filter(coverages, INDV == min_individual)$MEAN_DEPTH, "X)", sep = ""))
	

max_dosage_plot <- ggplot(data = filter(hap_table_dosages, INDV == max_individual)) +
	geom_point(mapping = aes(x = POS, y = DOSAGE_HOME, color = DOSAGE_NUMBER), size = point_size) +
	scale_color_discrete(name = legend_name) +
	guides(color = guide_legend(override.aes = list(size = legend_point_size))) +
	labs(x = "Locus (Mb)",
			 y = "Haplotype Dosage",
			 caption = paste("Individual #", max_individual,
			 								" (", filter(coverages, INDV == max_individual)$MEAN_DEPTH, "X)", sep = ""))

mean_dosage_plot <- ggplot(data = hap_table_dosages %>%
													 	group_by(POS, DOSAGE_NUMBER) %>%
													 	summarise(MEAN_DOSAGE = mean(DOSAGE_HOME))) +
	geom_point(mapping = aes(x = POS, y = MEAN_DOSAGE, color = DOSAGE_NUMBER), size = point_size) +
	scale_color_discrete(name = legend_name) +
	guides(color = guide_legend(override.aes = list(size = legend_point_size))) +
	labs(x = "Locus (Mb)",
			 y = "Haplotype Dosage",
			 caption = "Average over all individuals")


dosage_plots <- ggpubr::ggarrange(min_dosage_plot, max_dosage_plot, mean_dosage_plot,
																	ncol = 1, nrow = 3,
																	align = "v",
																	legend = "top",
																	common.legend = TRUE)
# save the plots
ggsave(plot = dosage_plots,
			 dpi = 320,
			 path = plots_folder, filename = "dosages.tiff")
	
##############

# CREATE FOUR-PANELED FIGURE FROM QUANTILES AND DOSAGE PLOTS
##############
dosages_quantiles_plot <- gridExtra::grid.arrange(dosage_plots, quantiles_plot,
																									ncol = 3, nrow = 1,
																									layout_matrix = rbind(c(1, 1, 2)))
# save the plots
in_height = 7
ggsave(plot = dosages_quantiles_plot,
			 units = "in", width = (in_height * 1.5), height = in_height,
			 dpi = 320,
			 path = plots_folder, filename = "dosages_quantiles.tiff")
##############

