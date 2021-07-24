# README
# This program creates a tsv of 3 phenotypes at a given locus:
# - Single SNP Causative
# - 10 SNPs Causative
# - All SNPs Causative
# In the output file, along with these phenotypes are individuals, the causative SNP, and the causative SNP chromosome.
# Phillip Long
# December 29, 2020

# INSTRUCTIONS TO RUN IN SHELL
##############
' .sh runner
module load R/3.6.2
# FOR CHROMOSOME 23:
srun -p free R --vanilla -f /share/adl/pnlong/mouseproject/software/generate_Y.R --args output_prefix="/share/adl/pnlong/mouseproject/save_data/phenotypes/Y_" causative_snp_address="/share/adl/pnlong/mouseproject/save_data/snp_table_Chr23.tsv" causative_snp_locus=100 range=25000 number_of_snps_per_csnp=10 individuals_address="/share/adl/pnlong/mouseproject/save_data/individuals.tsv" genetic_background_variance_address="/share/adl/pnlong/mouseproject/save_data/genetic_background_variance.tsv"

# FOR CHROMOSOME 19:
srun -p free R --vanilla -f /share/adl/pnlong/mouseproject/software/generate_Y.R --args output_prefix="/share/adl/pnlong/mouseproject/save_data/phenotypes/Y_" causative_snp_address="/share/adl/pnlong/mouseproject/save_data/snp_table_Chr19.tsv" causative_snp_locus=100 range=25000 number_of_snps_per_csnp=10 individuals_address="/share/adl/pnlong/mouseproject/save_data/individuals.tsv" genetic_background_variance_address="/share/adl/pnlong/mouseproject/save_data/genetic_background_variance.tsv"
'
##############

# LOADS
##############
library(tidyverse)
library(magrittr)

# a function to cleanly extract the arguments of an R script, returning a named vector where the names are argument names and values are argument values
get_args <- function() {
	args <- t(matrix(unlist(strsplit(x = commandArgs(trailingOnly = TRUE),
																	 split = "=")),
									 nrow = 2))
	
	vec <- args[,2]
	names(vec) <- args[,1]
	return(vec)
}

arguments <- get_args()
# for testing purposes:
# arguments <- c(output_prefix="/share/adl/pnlong/mouseproject/save_data/phenotypes/Y_", causative_snp_address="/share/adl/pnlong/mouseproject/save_data/snp_table_Chr19.tsv", causative_snp_locus="100", range="25000", number_of_snps_per_csnp="10", individuals_address="/share/adl/pnlong/mouseproject/save_data/individuals.tsv", genetic_background_variance_address="/share/adl/pnlong/mouseproject/save_data/genetic_background_variance.tsv")
##############

# LOAD IN INDIVIDUALS
##############
individuals_file = file(as.character(arguments["individuals_address"]))
individuals <- readLines(con = individuals_file, n = -1) %>%
	strsplit(x = ., split = "\t") %>% unlist()
close(individuals_file)
##############

# CAUSATIVE AND TEST SNP CHROMOSOME NUMBERS
##############
causative_snp_address = as.character(arguments["causative_snp_address"]) # the address of the SNP table from which the causative SNPs will be extracted
csnp_chromosome_number = strsplit(x = causative_snp_address, split = "\\.") %>% unlist() %>%
	extract(1) %>%
	strsplit(x = ., split = "_") %>% unlist() %>%
	extract(length(.)) %>%
	str_remove(string = ., pattern = "Chr")
##############

# FULL SNP TABLE
#############
# full snp table which we can pull markers from to be causative snps
full_snp_table <- read_tsv(file = causative_snp_address,
													 col_names = TRUE,
													 col_types = cols())
#############

# CAUSATIVE SNP INDICES
#############
# a function to determine if a given locus is valid, or if it contains errors (if values are all 0's or 2's)
is_valid_snp <- function(locus) {
	if (!(locus %in% full_snp_table$POS)) return(FALSE)
	
	snp_dosage <- unname(as.numeric(full_snp_table[(match(x = locus, table = full_snp_table$POS)[[1]]), -1]))
	if (sum(snp_dosage) == 0 | sum(snp_dosage) == (length(snp_dosage) * 2)) return(FALSE) else return(TRUE)
}

# a function that will find a new valid locus (if necessary) nearby, given a locus
choose_new_snp_if_necessary <- function(locus) {
	if (!(locus %in% full_snp_table$POS)) locus = full_snp_table$POS[[which.min(abs(full_snp_table$POS - locus))]]
	# if the locus is valid, then return the same locus, if not, run another iteration until a valid snp is found
	while(!is_valid_snp(locus)) {
		# selecting the new locus from the full_snp_table$POS because then that locus is also in the haplotype table
		# selecting the next locus basically
		locus = full_snp_table$POS[[(match(x = locus, table = full_snp_table$POS)[[1]]) + 1]]
	}
	return(locus) # return a position on the genome
}

causative_snp_loci <- choose_new_snp_if_necessary(locus = as.numeric(arguments["causative_snp_locus"]))

# now to start by extracting all loci within an inputted range of our causative snps
# this set (1-100 loci in causative_snp_loci) of sets (all loci within range of given locus) will be used for both the normal and rare caculations
range = as.numeric(arguments["range"])
within_range <- full_snp_table$POS[(full_snp_table$POS >= causative_snp_loci - range) & (full_snp_table$POS <= causative_snp_loci + range)] %>%
	.[Vectorize(is_valid_snp)(.) == TRUE] # make sure all of the loci are valid

# for the normal multiple causative snp test, just extract randomly from each set of loci in within_range
set.seed(0.001)
normal_multiple_snps_loci <- sort(sample(x = within_range, size = as.numeric(arguments["number_of_snps_per_csnp"])))
#############

# CAUSATIVE SNP TABLES
#############
# get the dosages of a given locus
get_dosages <- function(locus) {
	i = match(x = locus, table = full_snp_table$POS)[[1]] # [[1]] at the end just in case of duplicate entries
	return(unname(as.numeric(full_snp_table[i, -1])))
}

single_causative_snp_dosages <- get_dosages(locus = causative_snp_loci)

# take the 10 loci per causative SNP and sum them
multiple_causative_snps_dosages <- map(.x = normal_multiple_snps_loci, .f = get_dosages) %>%
	as.data.frame(col.names = normal_multiple_snps_loci) %>%
	mutate(SUM = rowSums(.)) %>% 
	pull(SUM)

# we will mimic the effect of rare by scaling down all the values very small of the range we have (within_range), then rescale again later
rare_causative_snps_dosages <- map(.x = within_range, .f = get_dosages) %>%
	map(.x = ., .f = ~ scale(.x)) %>%
	as.data.frame(col.names = within_range) %>%
	mutate(SUM = rowSums(.)) %>% 
	pull(SUM)
#############

# CREATE PHENOTYPE
#############
csnp_effect_percent = 0.05
ve_effect_percent = 0.50 # ve = environmental variance
set.seed(seed = causative_snp_loci) # will allow for a unique set of random numbers for the VE to be generated, but will allow these numbers to stay the same if the test is rerun

# background variance
# import gbv, which has genetic_background_variance for all individuals, but then filter so it only includes the individuals from individuals vector
# make the phenotype based off of the background variance
phenotype <- read_tsv(file = as.character(arguments["genetic_background_variance_address"]),
											col_names = TRUE) %>%
	filter(INDIVIDUAL %in% as.character(individuals)) %>%
	mutate(INDIVIDUAL = as.numeric(INDIVIDUAL)) %>% # convert individual to numeric so it can be reordered
	arrange(INDIVIDUAL) %>% # reorder least-to-greatest for individuals so it matches up with kinship matrix
	mutate(BG = (sqrt(ve_effect_percent - csnp_effect_percent) * c(scale(.$GBV_VALUE))), # rescale bg to 47.5% of 45% of total variance
				 VE = rnorm(n = nrow(.),
				 					 mean = 0,
				 					 sd = sqrt(ve_effect_percent)), # add environmental variance (which stays constant across all tests) at 50%
				 Y_single = c(scale((BG + VE + (sqrt(csnp_effect_percent) * c(scale(single_causative_snp_dosages)))))),
				 Y_multiple = c(scale((BG + VE + (sqrt(csnp_effect_percent) * c(scale(multiple_causative_snps_dosages)))))),
				 Y_rare = c(scale((BG + VE + (sqrt(csnp_effect_percent) * c(scale(rare_causative_snps_dosages)))))),
				 CSNP_CHROM = rep(csnp_chromosome_number, times = nrow(.)),
				 CSNP = rep(causative_snp_loci, times = nrow(.))) %>%
	select(CSNP_CHROM, CSNP, INDIVIDUAL, Y_single, Y_multiple, Y_rare)

write_tsv(x = phenotype,
					path = paste(as.character(arguments["output_prefix"]), "Chr", csnp_chromosome_number, ".", format(x = causative_snp_loci, scientific = FALSE, big.mark = ""), ".tsv", sep = ""),
					col_names = TRUE, # include column names
					append = FALSE)
#############