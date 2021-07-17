# README
# This program calculates some summary statistics for every single test in
# /share/adl/pnlong/mouseproject/results/Chr*
# Phillip Long
# November 26, 2020
' .sh runner
####
# 2.5% effect
resultsfolder="/share/adl/pnlong/mouseproject/results025"
# 5.0% effect
resultsfolder="/share/adl/pnlong/mouseproject/results050"

module load R/3.6.2
srun -p free R --vanilla -f /share/adl/pnlong/mouseproject/software/summary_statistics.R --args test_save="/share/adl/pnlong/mouseproject/save_data/hapsnp_Chr23.Rsave" control_save="/share/adl/pnlong/mouseproject/save_data/hapsnp_Chr19.Rsave" results_folder=$resultsfolder
####
'

# LOADS
#############
library(tidyverse)
library(magrittr)
load("/share/adl/pnlong/mouseproject/software/useful_functions.Rsave")
arguments <- get_args()
# for testing purposes:
# arguments <- c(test_save="/share/adl/pnlong/mouseproject/save_data/hapsnp_Chr23.Rsave", control_save="/share/adl/pnlong/mouseproject/save_data/hapsnp_Chr19.Rsave", results_folder="/share/adl/pnlong/mouseproject/results050")

files_to_keep_from_saves <- c("individuals", "csnp_chromosome_number", "causative_snp_loci", "single_causative_snp_dosages", "multiple_causative_snps_dosages", "rare_multiple_causative_snps_dosages")
# c("csnp_chromosome_number", "tsnp_chromosome_number", "individuals", "loci_from_hap_table", "hap_table", "causative_snp_loci", "single_causative_snp_dosages", "multiple_causative_snps_dosages", "rare_multiple_causative_snps_dosages", "kinship_matrix", "background_variance")

load(file = as.character(arguments["test_save"]),
     envir = (Chr_test_env <- new.env())) # Chr23, store in its own environment
remove(list = ls(Chr_test_env)[!(ls(Chr_test_env) %in% files_to_keep_from_saves)], # get rid of objects not in list of objects we want to save
       envir = Chr_test_env)

load(file = as.character(arguments["control_save"]),
     envir = (Chr_control_env <- new.env())) # Chr19, store in its own environment
remove(list = ls(Chr_control_env)[!(ls(Chr_control_env) %in% files_to_keep_from_saves)], # get rid of objects not in list of objects we want to save
       envir = Chr_control_env)

# Access environment variables with:
# environment_name$variable_name
# (e.g. Chr_test_env$x or Chr_control_env$x)

individuals <- Chr_test_env$individuals # since individuals is the same for both
remove(individuals, envir = Chr_test_env)
remove(individuals, envir = Chr_control_env)
#############

# FILE LOCATIONS
#############
results_folder = as.character(arguments["results_folder"])
files <- list.files(path = results_folder,
                    recursive = TRUE, # recursive = TRUE lists every single file from every sub and sub-sub directory in results_folder
                    all.files = FALSE, # eliminates any invisible files
                    full.names = TRUE) %>% # returns absolute paths rather than relative paths
  str_sort(numeric = TRUE)
#############

# SIGNIFICANCE THRESHOLDS
#############
sig_thresholds <- c("snp" = 7.5, "hap" = 6.0) ################################################################################################################################################################################# Significance Thresholds

# FOR DEMONSTRATION SCANS LATER ON??????
# sig_thresholds <- c("snp" = 5.0, "hap" = 4.0)

#############

# DETERMINE COLUMN NAMES
#############
# name the vector
stat_names <- c("tSNP_chromosome", "cSNP_chromosome", "test_type", "calculation_method", "snp_index", "cSNP", "freqcSNP", "sig_threshold", "has_hit", "n_hits", "avgdist_hits", "hit_span", "avgdist_MSM", "MSM_LOD")

# GOALS OF SUMMARY TABLE:
# % of scans with hits
# average # of hits; hit
# average distance MSM from cSNP; hit
# average "span" of hits; hit
# average LOD of MSM

#############

# DEFINE SUMMARY STATISTIC FUNCTION
#############
summary_statistics <- function(filename) {
  
  # LOAD TABLE, EXTRACT CAUSATIVE SNP INDEX AND FIRST LINE INFORMATION
  #############
  file_info <- scan(file = filename, what = character(), nlines = 1, quiet = TRUE) %>% # read the first line of the file in
    extract(-1) %>% # since first value is a "#"
    strsplit(split = "=") %>%
    two_dimensional_vector_to_tibble(x = ., names = c("KEY", "VALUE")) %>%
    deframe() # converts to nameed vector
  
  LOD_table <- read_tsv(file = filename, col_names = TRUE, comment = "#", col_types = cols(), progress = FALSE) %>%
    correct_LOD_values(x = .) %>%
    filter(LOD_METHOD == 1) # select only method 1 (the approximation method)
  
  input_snp_index = as.numeric(unlist(strsplit(x = filename, split = "\\."))[[2]])
  #############
  
  # INITIALIZE MAIN VECTOR
  #############
  # create return vector
  main_vector <- vector(mode = "character", length = length(stat_names))
  # vector indexer
  index = 1
  #############
  
  # GET THE TABLE INFORMATION INTO THE VECTOR
  #############
  # c("tSNP_chromosome", "cSNP_chromosome", "test_type", "calculation_method", "cSNP") # first 5 columns
  
  # tSNP_chromosome
  main_vector[index] = as.character(file_info["tSNP_chromosome"])
  # cSNP_chromosome
  cSNP_chromosome = file_info["cSNP_chromosome"]
  main_vector[index + 1] = as.character(cSNP_chromosome)
  
  is_test_chr = case_when(
    cSNP_chromosome == as.character(Chr_test_env$csnp_chromosome_number) ~ TRUE,
    cSNP_chromosome == as.character(Chr_control_env$csnp_chromosome_number) ~ FALSE
  )
  if (is_test_chr) Chr_env <- Chr_test_env else Chr_env <- Chr_control_env
  
  # test_type
  test_type = file_info["test_type"]
  main_vector[index + 2] = test_type
  # calculation_method
  calculation_method = file_info["calculation_method"] %>%
    strsplit(split = "_") %>% # [[1]] [1] "hap" "based"
    unlist() %>% # > [1] "hap" "based"
    extract(1)
  main_vector[index + 3] = calculation_method
  
  # input snp index
  main_vector[index + 4] = as.character(input_snp_index)
  
  # causative snp loci
  causative_snp_locus = as.numeric(file_info["cSNP_locus"])
  main_vector[index + 5] = as.character(causative_snp_locus) # since the vector is a character vector at the moment

  index = index + 6
  #############
  
  # CAUSATIVE SNP FREQUENCY
  #############
  # get the dosages (depending on the test type)
  causative_snp_dosages <- case_when(
    test_type == "single" ~ Chr_env$single_causative_snp_dosages[[input_snp_index]],
    test_type == "multiple" ~ Chr_env$multiple_causative_snps_dosages[[input_snp_index]],
    test_type == "rare" ~ Chr_env$rare_multiple_causative_snps_dosages[[input_snp_index]]
  )

  # frequency is what percentage of the alleles express the causative snp
  frequency = sum(causative_snp_dosages) / (length(individuals) * 2) # *2 because 2 times the number of individuals is number of alleles
  
  main_vector[index] = as.character(frequency)
  index = index + 1
  #############
  
  # SIG THRESHOLD
  #############
  sig_threshold = as.numeric(sig_thresholds[calculation_method])
  
  # significance threshold being used
  main_vector[index] = as.character(sig_threshold)
  index = index + 1
  #############
  
  # DOES TEST HAVE A HIT
  #############
  hits <- LOD_table %>%
    filter(LOD_VALUE >= sig_threshold) %>%
    arrange(TEST_SNP)
  
  hits_count = nrow(hits)
  has_hits = (hits_count > 0)
  
  main_vector[index] = as.character(has_hits)
  main_vector[index + 1] = as.character(hits_count)
  main_vector[index + 2] = if_else(condition = is_test_chr,
                                   true = hits %>%
                                     mutate(DIST = abs(causative_snp_locus - TEST_SNP)) %>%
                                     pull(DIST) %>%
                                     mean(na.rm = TRUE) %>%
                                     as.character(),
                                   false = NA_character_)
  main_vector[index + 3] = if_else(condition = has_hits,
                                   true = as.character(max(hits$TEST_SNP) - min(hits$TEST_SNP)),
                                   false = NA_character_)
  index = index + 4
  #############
  
  
  
  # FIND MSM FIND DISTANCE BETWEEN MSM AND CAUSATIVE SNP AND AMOUNT OF MSM
  #############
  # get the max lod score (the MSM)
  maxlod = max(LOD_table$LOD_VALUE)
  
  # calculate the average of the absolute difference between the MSM(s) and cSNP
  msm_dist <- if_else(condition = is_test_chr,
                      true = LOD_table %>%
                        filter(LOD_VALUE == maxlod) %>%
                        pull(TEST_SNP) %>%
                        subtract(causative_snp_locus) %>%
                        abs() %>%
                        mean() %>%
                        as.character(),
                      false = NA_character_) # the distance is wrong and thus NA if the causative snp is outside of the test snp chromosome
  
  # add average absolute difference of significant markers to the main vector
  main_vector[index] = msm_dist
  main_vector[index + 1] = as.character(maxlod)
  
  index = index + 2
  #############
  
  return(main_vector)
}
#############

# RUN SUMMARY STATISTICS ON FILES AND WRITE TO A TSV
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

summary_values <- files[is_data_file(files)] %>%
  map(.x = ., .f = summary_statistics)

# convert to tibble
summary_table <- two_dimensional_vector_to_tibble(x = summary_values, names = stat_names) %>%
  type_convert() %>%
  arrange(desc(cSNP_chromosome), cSNP, test_type, calculation_method)

write_tsv(x = summary_table,
          path = paste(results_folder, "summary_statistics.tsv", sep = "/"),
          col_names = TRUE)
#############