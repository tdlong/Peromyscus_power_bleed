# README
# This program will save() the CONSTANTS needed to calculate LOD scores with:
# - a single causative SNP
# - multiple causative SNPs
# - multiple rare causative SNPs (by mimicking the effect of rare haplotypes through scaling values so that they are very small)
#
# For each of these 3 tests, 2 methods are used to calculate LOD scores:
# - haplotype based
# - SNP dosage based
#
# Finally, these 6 tests (3 x 2) need to be done twice with causative SNPs on a different chromosome (a negative control):
# - Chromosome 23 (normal)
# - Chromosome 19 (negative control with the causative snps on Chromosome 19)
#
# Again, this program sets up the constants for these tests in a file that will be loaded for scanning in other R scripts.
# Phillip Long
# December 29, 2020

# INSTRUCTIONS TO RUN IN SHELL
##############
' .sh runner
module load R/3.6.2
# FOR CHROMOSOME 23:
srun -p free R --vanilla -f /share/adl/pnlong/mouseproject/software/hapsnp_calculation_save.R --args tsnp_chromosome_number=23 output="/share/adl/pnlong/mouseproject/save_data/hapsnp_Chr23.Rsave" causative_snp_address="/share/adl/pnlong/mouseproject/save_data/snp_table_Chr23.tsv" desired_number_of_causative_snps=100 range=25000 number_of_snps_per_csnp=10

# FOR CHROMOSOME 19:
srun -p free R --vanilla -f /share/adl/pnlong/mouseproject/software/hapsnp_calculation_save.R --args tsnp_chromosome_number=23 output="/share/adl/pnlong/mouseproject/save_data/hapsnp_Chr19.Rsave" causative_snp_address="/share/adl/pnlong/mouseproject/save_data/snp_table_Chr19.tsv" desired_number_of_causative_snps=100 range=25000 number_of_snps_per_csnp=10
'
##############

# LOADS
##############
library(tidyverse)
library(magrittr)
load(file = "/share/adl/pnlong/mouseproject/software/useful_functions.Rsave")
arguments <- get_args()
# for testing purposes:
# arguments <- c(tsnp_chromosome_number="23", output="/share/adl/pnlong/mouseproject/save_data/hapsnp_Chr19.Rsave", causative_snp_address="/share/adl/pnlong/mouseproject/save_data/snp_table_Chr19.tsv", desired_number_of_causative_snps="100", range="25000", number_of_snps_per_csnp="10")
##############

# CAUSATIVE AND TEST SNP CHROMOSOME NUMBERS
##############
causative_snp_address = as.character(arguments["causative_snp_address"]) # the address of the SNP table from which the causative SNPs will be extracted
csnp_chromosome_number = str_extract(string = causative_snp_address, pattern = "[[:digit:]]+") # csnp = causative SNP

tsnp_chromosome_number = as.character(arguments["tsnp_chromosome_number"]) # tsnp_chromosome = testing SNP chromosome (the chromosome the scan will be performed on)
##############

# FULL SNP TABLE
#############
# full snp table which we can pull markers from to be causative snps
full_snp_table <- read_tsv(file = causative_snp_address,
                           col_names = TRUE)
#############

# CAUSATIVE SNP INDICES
#############
# a function to determine if a given locus is valid, or if it contains errors (if values are all 0's or 2's)
is_valid_snp <- function(locus) {
  snp_dosage <- unname(as.numeric(full_snp_table[(match(x = locus, table = full_snp_table$POS)[[1]]), -1]))
  if (sum(snp_dosage) == 0 | sum(snp_dosage) == (length(snp_dosage) * 2)) return(FALSE) else return(TRUE)
}

# a function that will find a new valid locus (if necessary) nearby, given a locus
choose_new_snp_if_necessary <- function(locus) {
  # if the locus is valid, then return the same locus, if not, run another iteration until a valid snp is found
  while(!is_valid_snp(locus)) {
    # selecting the new locus from the full_snp_table$POS because then that locus is also in the haplotype table
    # selecting the next locus basically
    locus = full_snp_table$POS[[(match(x = locus, table = full_snp_table$POS)[[1]]) + 1]]
  }
  return(locus) # return a position on the genome
}

desired_number_of_causative_snps = as.numeric(arguments["desired_number_of_causative_snps"])

causative_snp_loci <- full_snp_table$POS[seq(from = round(nrow(full_snp_table) * 0.25), to = round(nrow(full_snp_table * 0.75)), length.out = desired_number_of_causative_snps)] %>%
  map(.f = choose_new_snp_if_necessary) %>% unlist()
# full_snp_table$POS assures that the causative snps are in both the haplotype table and snp table
# causative_snp_loci is a list of (generally) 100 SNPs that will be used as causative in the single-SNP causative model, the most basic of the three

# now to start by extracting all loci within an inputted kB of our causative snps
# this set (1-100 loci in causative_snp_loci) of sets (all loci within range of given locus) will be used for both the normal and rare caculations
kB = as.numeric(arguments["range"])
within_kB_range <- map(.x = causative_snp_loci,
                       .f = function(locus) full_snp_table$POS[(full_snp_table$POS >= locus - kB) & (full_snp_table$POS <= locus + kB)]) %>%
  map(.f = function(loci_list) { loci_list[Vectorize(is_valid_snp)(loci_list) == TRUE] }) # make sure all of the loci are valid
# withink_kB_range includes all of the loci, not just the ones in hap_table, since the actual position value doesnt matter

# for the normal multiple causative snp test, just extract randomly from each set of loci in within_Kb_range
normal_multiple_snps_loci <- map(.x = within_kB_range,
                                 .f = ~ sort(sample(x = .x, size = as.numeric(arguments["number_of_snps_per_csnp"]))))
#############

# CAUSATIVE SNP TABLES
#############
# get the dosages of a given locus
get_dosages <- function(locus) {
  i = match(x = locus, table = full_snp_table$POS)[[1]] # [[1]] at the end just in case of duplicate entries
  return(unname(as.numeric(full_snp_table[i, -1])))
}

single_causative_snp_dosages <- map(.x = causative_snp_loci,
                                    .f = get_dosages)

# take the 10 loci per causative SNP and sum them
multiple_causative_snps_dosages <- map(.x = normal_multiple_snps_loci,
                                       .f = function(loci_list) {
                                         map(.x = loci_list, .f = get_dosages) %>%
                                           as.data.frame(col.names = loci_list) %>%
                                           mutate(SUM = rowSums(.)) %>% 
                                           pull(SUM)
                                       })

# we will mimic the effect of rare by scaling down all the values very small of the range we have (within_kB_range), then rescale again later
rare_multiple_causative_snps_dosages <- map(.x = within_kB_range,
                                            .f = function(loci_list) {
                                              map(.x = loci_list, .f = get_dosages) %>%
                                                map(.f = ~ scale(.x)) %>%
                                                as.data.frame(col.names = loci_list) %>%
                                                mutate(SUM = rowSums(.)) %>% 
                                                pull(SUM)
                                            })
#############

save(list = c("csnp_chromosome_number", "tsnp_chromosome_number", "causative_snp_loci", "normal_multiple_snps_loci", "within_kB_range", "single_causative_snp_dosages", "multiple_causative_snps_dosages", "rare_multiple_causative_snps_dosages"),
     file = as.character(arguments["output"]))