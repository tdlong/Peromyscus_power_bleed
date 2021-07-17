# README
# This program will validate STITCH with RNAseq data collected from livers.
# Phillip Long
# May 1, 2020
' .sh runner
####
module load R/3.6.2
mkdir /share/adl/pnlong/mouseproject/plots_validation
srun -p standard -c 5 R --vanilla -f /share/adl/pnlong/mouseproject/software/stitch_validation.R --args dna_location="/share/adl/pnlong/mouseproject/save_data/stitch_genomewide.tsv" rna_location="/share/adl/pnlong/mouseproject/save_data/rnaseq.tsv" kinship_address="/share/adl/pnlong/mouseproject/save_data/Mjj_values_all_pnl.tsv" output_path="/share/adl/pnlong/mouseproject/plots_validation"
####
'

# LOADS
#############
library(tidyverse)
library(magrittr)
# install.packages("ggpubr")
library(ggpubr)
library(gridExtra)

load("/share/adl/pnlong/mouseproject/software/useful_functions.Rsave")
# arguments <- c(dna_location="/share/adl/pnlong/mouseproject/save_data/stitch_genomewide.tsv", rna_location="/share/adl/pnlong/mouseproject/save_data/rnaseq.tsv", kinship_address="/share/adl/pnlong/mouseproject/save_data/Mjj_values_all_pnl.tsv", output_path="/share/adl/pnlong/mouseproject/plots_validation")
arguments <- get_args()
output_path = as.character(arguments["output_path"])
dna_location = as.character(arguments["dna_location"]) # "/share/adl/pnlong/mouseproject/save_data/stitch_genomewide.tsv"
# location of STITCH imputation values ^^
rna_location = as.character(arguments["rna_location"]) # "/share/adl/pnlong/mouseproject/save_data/rnaseq.tsv"
# location of RNAseq values ^^
#############

# IMPORT GDNA AND RNASEQ TABLES AND MERGE
#############
# 6 individuals we chose with RNAseq spleen data
individuals_file = file("/share/adl/pnlong/mouseproject/save_data/rnaseq_individuals.tsv")
individuals <- readLines(con = individuals_file, n = -1) %>%
	strsplit(x = ., split = "\t") %>% unlist()
close(individuals_file)

# STITCH imputations data
dna <- read_tsv(file = dna_location, col_names = TRUE, col_types = cols(col_character(), col_double(), col_double(), col_double()))

# pivot_longer(cols = all_of(individuals), names_to = "ID", values_to = "DNA_GENO") %>%
# filter(TRT == "DF" & K == 8) %>% # since the only STITCH Data we care to validate is the DNA FULL and K-value of 8 is the same as K-value of 16, so I just had to choose one or the other
# select(-TRT, -K)
# # * "DF" = DNA FULL means I ran STITCH the normal way using roughly 300 mice for which we have low pass sequence data, or which the 6 animals above are the subset for which we happen to have spleen RNAseq data.

# RNAseq data
rna <- read_tsv(file = rna_location, col_names = TRUE) %>%
	select(-8, -10:-18) %>%
	rename("CHROM" = Chr, "POS" = Pos) %>%
	pivot_longer(cols = paste("R", individuals, sep = ""), names_to = "ID", values_to = "RNA_GENO") %>%
	mutate(ID = as.numeric(str_remove(string = ID, pattern = "R")),
				 CHROM = str_remove(string = CHROM, pattern = "Chr"))

# merge the RNAseq and STITCH imputations on shared SNPs
dna_rna <- inner_join(x = dna, y = rna, by = c("CHROM", "POS", "ID")) %>% # merge on shared SNPs between RNAseq and STITCH imputation dosages
	select(-CHROM) %>%
	mutate(DIFF = abs(DNA_GENO - RNA_GENO),
				 ID = as.character(ID)) # calculate absolute error

remove(dna, rna)
#############

# ERROR AS A FUNCTION OF COVERAGE PLOTS
#############
# get the subset of 297 individuals that we are running QTL detection on
individuals_subset <- read_lines(file = "/share/adl/pnlong/mouseproject/save_data/individuals.tsv", n_max = 1) %>% # read the first line of the haplotype table (the column names), since that includes the individual names
	strsplit(x = ., split = "\t") %>%
	unlist() %>%
	as.character()

# coverages for gDNA in 6 individuals
# filter to 297 individuals we chose or 6 RNAseq data individuals
coverages <- read_tsv(file = "/share/adl/pnlong/mouseproject/save_data/coverages.tsv", col_names = TRUE) %>%
	rename("ID" = INDV, "COVERAGE" = MEAN_DEPTH) %>%
	filter(TYPE == "gDNA" & (ID %in% individuals_subset | ID %in% individuals)) %>% # the RNAseq coverages are unimportant, since the RNAseq is merely used as a validator and not in the experiment, so filter out RNAseq coverage info
	mutate(COVERAGE = log10(COVERAGE))

# calculate percent of SNPs with error over 0.10 per individual
dna_rna_coverage <- dna_rna %>%
	group_by(ID) %>%
	summarise(ERROR_PERCENT = (100 * (sum(if_else(condition = (DIFF >= 0.10), TRUE, FALSE)) / n())),
						CORRELATION = cor(x = RNA_GENO, y = DNA_GENO)) %>%
	left_join(x = ., y = select(coverages, ID, COVERAGE), by = c("ID")) # will automatically filter out coverages not in individuals

# create shared x-axis limits for the two graphs, so they are aligned and easier to read
coverage_binwidths = 0.1
coverage_lower_limit = min(dna_rna_coverage$COVERAGE) - (2 * coverage_binwidths)
histo_outlier_count = sum((coverages$COVERAGE <= coverage_lower_limit)) # the number of observations that get cut off by the lower limit in the histogram
coverage_upper_limit = max(coverages$COVERAGE) + (2 * coverage_binwidths) # *2 to provide some buffer space and not cut anything off at the upper limit

coverage_label = "Raw Sequence Coverage (log10)"

# percent error as a function of coverage plot
coverage_vs_percent_error_plot <- ggplot(data = dna_rna_coverage, mapping = aes(x = COVERAGE, y = ERROR_PERCENT)) +
	geom_point(mapping = aes(color = ID), size = 3) + # geom_smooth(method = "lm", se = FALSE, color = "grey") + # add line of best fit
	ylim(0, NA) + # set lower limit to 0, but automatic upper limit
	xlim(coverage_lower_limit, coverage_upper_limit) + # set shared x-axis limits for the two graphs
	labs(x = coverage_label,
			 y = "% of SNPs with Absolute Error >0.10") +
	scale_color_discrete(name = "Deermouse ID") +
	theme(legend.position = "top")

# histogram showing coverage frequencies
coverage_histo <- ggplot(data = filter(coverages, ID %in% individuals_subset)) + # filter out the 6 individuals with RNAseq data that will not be used later on in the QTL Detection experiment (so just the 297 chosen ones)
	geom_histogram(mapping = aes(x = COVERAGE), binwidth = coverage_binwidths) +
	xlim(coverage_lower_limit, coverage_upper_limit) + # set shared x-axis limits for the two graphs
	labs(x = coverage_label,
			 y = "Number of Individuals")

print(paste(histo_outlier_count, " outliers with coverages < ", format(round(coverage_lower_limit, digits = 1), nsmall = 2), " are excluded. -0.5, 0, and 0.5 are 0.3X, 1X, and 3X respectively.", sep = ""))


coverage_plots <- ggpubr::ggarrange(coverage_vs_percent_error_plot, coverage_histo,
														ncol = 1, nrow = 2,
														align = "v",
														common.legend = TRUE, legend = "top")
ggsave(plot = coverage_plots,
			 dpi = 320,
			 path = output_path, filename = "coverages.tiff")
#############

# CUMULATIVE ERROR
#############
# cumulative error table
cumulative_error <- dna_rna %>%
	select(-c(POS, RNA_GENO, DNA_GENO)) %>%
	arrange(ID, DIFF) %>%
	mutate(PERCENT = 100 * c(unlist(map(.x = (group_by(.data = ., ID) %>% summarise(N = n()) %>% pull(N)),
																			.f = ~ c((1:.x)/.x)))))

# cumulative error lineplot
cumulative_error_plot <- ggplot(data = cumulative_error) +
	geom_line(mapping = aes(x = DIFF, y = PERCENT, color = ID)) +
	lims(x = c(0, 1.1), # just show errors from 0.0 to 1.1 absolute difference (since 99+% have errors less than 1.0)
			 y = c(75, 100)) + # just show 75% and up
	labs(x = "Absolute Error",
			 y = "% of SNPs") +
	scale_color_discrete(name = "Deermouse ID")
ggsave(plot = cumulative_error_plot,
			 dpi = 320,
			 path = output_path, filename = "cumulative_error.tiff")

#############

# MAKE PANELED PLOT FROM COVERAGES PLOTS AND CUMULATIVE ERROR PLOT
#############
coverages_cumulative_error_plots <- gridExtra::grid.arrange(coverage_vs_percent_error_plot, coverage_histo, (cumulative_error_plot + theme(legend.position = "none")),
																														nrow = 2, ncol = 3,
																														layout_matrix = rbind(c(1, 1, 3),
																																									c(2, 2, 3)))

in_width = 9
ggsave(plot = coverages_cumulative_error_plots,
			 units = "in", width = in_width, height = round(in_width * (2/3), digits = 2),
			 dpi = 320,
			 path = output_path, filename = "coverages_cumulative_error.tiff")

#############

# KINSHIP MATRIX PLOTS
#############
individuals <- read_lines(file = "/share/adl/pnlong/mouseproject/save_data/individuals.tsv", n_max = 1) %>% # read the first line of the haplotype table (the column names), since that includes the individual names
	strsplit(x = ., split = "\t") %>%
	unlist() %>%
	as.character()

# import full Mjj Values table set that will later be used in kinship matrix
Mjj_values_normal <- read_tsv(file = as.character(arguments["kinship_address"]),
															col_names = TRUE,
															col_types = list(col_character(), col_character(), col_double()))



# CODE TO GENERATE K-VALUE
# K-VALUE is a constant that will rescale Mjj values so that the average of the two parent-offsprings is 0.5
K_parents <- c("19998", "20036") # chosen because they have a known relationship with individual 20351
K_offspring <- ("20351")
K_individuals <- c(K_parents, K_offspring)
K_Mjj_values_subset <- Mjj_values_normal %>%
	filter(INDIVIDUAL.x %in% K_individuals, INDIVIDUAL.y %in% K_individuals) %>% # filter down to the K individuals
	filter(!(INDIVIDUAL.x %in% K_parents & INDIVIDUAL.y %in% K_parents)) # filter so that the Mjj values we have are ones that show a parent-offspring relationship, not a parent-parent relationship

# take the average of the two parent-offspring Mjj values
K_parentoffspring_average = K_Mjj_values_subset %>%
	pull(MJJ) %>%
	mean()

# calculate the K constant that will make the parent-offspring average be 0.5
# K = 0.5 / ((0.619 + 0.635) / 2) = 0.7973824
K = 0.5 / K_parentoffspring_average


Mjj_values_normal %<>%
	filter((INDIVIDUAL.x %in% individuals & INDIVIDUAL.y %in% individuals) | (INDIVIDUAL.x %in% K_individuals & INDIVIDUAL.y %in% K_individuals)) %>% # filter so there are only the individuals used in this program
	mutate_if(is.character, as.numeric) %>% # convert INDIVIDUAL columns to numeric so they can be arranged
	arrange(INDIVIDUAL.x, INDIVIDUAL.y) %>%
	mutate(MJJ = MJJ * K)


kinship_plot <- ggplot() +
	geom_histogram(data = Mjj_values_normal, mapping = aes(x = MJJ), binwidth = 0.001) +
	geom_segment(data = K_Mjj_values_subset, mapping = aes(x = MJJ * K, y = 200, xend = MJJ * K, yend = 20), arrow = arrow(), lineend = "round", linejoin = "mitre", size = 1, color = "red") +
	labs(x = "Relatedness",
			 y = "Count")

ggsave(plot = kinship_plot,
			 dpi = 320,
			 path = output_path, filename = "kinship_matrix.tiff")
#############

