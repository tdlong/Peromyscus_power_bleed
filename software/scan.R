# README
# This program will calculate LOD scores for the loci of a given chromosome
# using both a marker- and haplotype-based calculation method.
# Phillip Long
# December 29, 2020
####################################################


####################################################
# INGREDIENTS FOR BOTH

# LOADS
#############
library(tidyverse)
library(magrittr)
library(devtools)
library(lmtest)
# install_github("variani/lme4qtl")
library(lme4qtl)
# install_github("variani/matlm")
library(matlm)
# install_github("variani/wlm")
library(wlm)
library(Matrix)
library(LaF)

# a function to cleanly extract the arguments of an R script, returning a named vector where the names are argument names and values are argument values
get_args <- function() {
	args <- t(matrix(unlist(strsplit(x = commandArgs(trailingOnly = TRUE),
																	 split = "=")),
									 nrow = 2))
	
	vec <- args[,2]
	names(vec) <- args[,1]
	return(vec)
}
wc_l <- function(file) {
	if (!file.exists(file)) {
		return(0)
	} else {
		x <- file(file)
		count = length(readLines(con = x))
		close(x)
		
		return(count)
	}
}

# FOR TESTING
# arguments <- c(test_type="single", phenotype_address="", phenotype_name="Y_single", individuals_address="/share/adl/pnlong/mouseproject/save_data/individuals.tsv", snp_table_address="/share/adl/pnlong/mouseproject/save_data/snp_table_Chr23.tsv", hap_table_address="/share/adl/pnlong/mouseproject/save_data/hap_table_Chr23.tsv", kinship_address="/share/adl/pnlong/mouseproject/save_data/Mjj_values_all_pnl.tsv", batch_size="100", rows_in_chunk="5000", drop_PC_less_than="0.05", LOD_threshold="2", snp_out="/share/adl/pnlong/mouseproject/results/Chr23_single_snpbased/single_csnp_snpbased", hap_out="/share/adl/pnlong/mouseproject/results/Chr23_single_hapbased/single_csnp_hapbased", completed_files="/share/adl/pnlong/mouseproject/results/completed_tests.txt")
arguments <- get_args()
completed_tests_path = as.character(arguments["completed_files"])
#############

# LOAD IN INDIVIDUALS
##############
individuals_file = file(as.character(arguments["individuals_address"]))
individuals <- readLines(con = individuals_file, n = -1) %>%
	strsplit(x = ., split = "\t") %>% unlist()
close(individuals_file)
##############

# PHENOTYPE AND SCAN INFO
#############
phenotype <- read_tsv(file = as.character(arguments["phenotype_address"]), col_names = TRUE, col_types = cols(INDIVIDUAL = col_character())) %>%
	filter(INDIVIDUAL %in% individuals) %>%
	arrange(as.numeric(INDIVIDUAL)) # alphanumerically sort INDIVIDUAL IDs in ascending order

individuals <- unique(phenotype$INDIVIDUAL)

snp_table_address = as.character(arguments["snp_table_address"])
tsnp_chromosome_number = strsplit(x = snp_table_address, split = "\\.") %>% unlist() %>%
	extract(1) %>%
	strsplit(x = ., split = "_") %>% unlist() %>%
	extract(length(.)) %>%
	str_remove(string = ., pattern = "Chr")

test_type = trimws(tolower(as.character(arguments["test_type"]))) %>%
	ifelse(test = (. %in% c("single", "multiple", "rare")),
				 yes = .,
				 no = "\b \b")
if (("CSNP_CHROM" %in% colnames(phenotype)) & ("CSNP" %in% colnames(phenotype)) & (test_type %in% c("single", "multiple", "rare"))) {
	# This information will later be printed to the first line of each of the 2 output files
	# (another  value will be added later: hap-based or snp-based):
	# - Causative SNP Chromosome Number
	# - Test SNP Chromosome Number
	# - Test Type
	# - Calculation Method (to be added later, though)
	
	csnp_chromosome_number = as.character(unique(phenotype$CSNP_CHROM))
	csnp_chrom = paste("cSNP_chromosome", csnp_chromosome_number,
										 sep = "=")
	
	causative_snp = format(x = unique(phenotype$CSNP), scientific = FALSE, big.mark = "")
	csnp_position = paste("cSNP_locus", causative_snp,
												sep = "=")
	
	tsnp_chrom = paste("tSNP_chromosome", tsnp_chromosome_number,
										 sep = "=")

	test_type_text = paste("test_type", test_type,
												 sep = "=")
	
	
	scan_info <- paste(csnp_chrom, csnp_position, tsnp_chrom, test_type_text, sep = "\t")
} else {
	tsnp_chrom = paste("tSNP_chromosome", tsnp_chromosome_number,
										 sep = "=")
	
	csnp_chromosome_number = tsnp_chromosome_number
	causative_snp = "\b \b"
	
	scan_info <- paste(tsnp_chrom, sep = "\t")
}

phenotype <- phenotype %>%
	select(INDIVIDUAL, as.character(arguments["phenotype_name"])) %>%
	set_colnames(c("INDIVIDUAL", "Y"))
#############

# KINSHIP MATRIX
#############
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
	filter(INDIVIDUAL.x %in% K_individuals, INDIVIDUAL.y %in% K_individuals) # filter down to the K individuals

# take the average of the two parent-offspring Mjj values
K_parentoffspring_average = K_Mjj_values_subset %>%
	filter(!(INDIVIDUAL.x %in% K_parents & INDIVIDUAL.y %in% K_parents)) %>% # filter so that the Mjj values we have are ones that show a parent-offspring relationship, not a parent-parent relationship
	pull(MJJ) %>%
	mean()

# calculate the K constant that will make the parent-offspring average be 0.5
# K = 0.5 / ((0.619 + 0.635) / 2) = 0.7973824
K = 0.5 / K_parentoffspring_average



Mjj_values_normal %<>%
	filter(INDIVIDUAL.x %in% individuals, INDIVIDUAL.y %in% individuals) %>% # filter so there are only the individuals used in this program
	mutate_if(is.character, as.numeric) %>% # convert INDIVIDUAL columns to numeric so they can be arranged
	arrange(INDIVIDUAL.x, INDIVIDUAL.y) %>%
	mutate(MJJC = MJJ * K) %>%
	select(-MJJ)

# now, convert the table of kinship values (which is shaped like a triangle), into a matrix

Mjj_values_flip <- Mjj_values_normal %>%
	rename(INDIVIDUAL.x = INDIVIDUAL.y, INDIVIDUAL.y = INDIVIDUAL.x)

Mjj_values_diagonal <- tibble(
	INDIVIDUAL.x = as.numeric(individuals),
	INDIVIDUAL.y = as.numeric(individuals),
	MJJC = rep(x = 1.0, times = length(individuals))
)

Mjj_values <- bind_rows(Mjj_values_normal, Mjj_values_flip, Mjj_values_diagonal) %>%
	arrange(INDIVIDUAL.x, INDIVIDUAL.y)
remove(Mjj_values_normal, Mjj_values_flip, Mjj_values_diagonal)

kinship_matrix <- matrix(data = Mjj_values$MJJC, nrow = length(individuals), ncol = length(individuals), dimnames = list(individuals, individuals))

remove(Mjj_values) # to save memory
#############

####################################################



####################################################
# MARKER
snp_out_path = trimws(paste(as.character(arguments["snp_out"]), paste(paste("Chr", csnp_chromosome_number, "_", test_type, "_snpbased",  sep = ""), causative_snp, "tsv", sep = "."), sep = "/")) # address of the output file for the marker-based scan

if (!file.exists(snp_out_path)) { # if snp_out_path does not exist yet
	# GFIT NULL MODEL
	#############
	# fit the model without causative SNP
	gfit_null_snp <- lme4qtl::relmatLmer(formula = Y ~ (1 | INDIVIDUAL),
																			 data = phenotype,
																			 relmat = list(INDIVIDUAL = kinship_matrix),
																			 REML = TRUE)
	#############
	
	# VARIANCE/COVARIANCE MATRIX, DEAL WITH SMALL VALUES AND EVD COMPOSITION
	#############
	# extract something called the variance/covariance matrix associated with the null model
	V <- lme4qtl::varcov(gfit_null_snp, idvar = "INDIVIDUAL")
	remove(gfit_null_snp)
	# take very small values in the matrix and force them to equal zero
	V[abs(V) < 1e-10] <- 0
	
	# carry out something called a EVD composition of the var/covar matrix
	# this is like a Cholesky decomposition (2nd year college matrix algebra) but numerically more stable
	# this is basically an “expensive” transformation of the var/covar matrix to a “numerically better form” for downstream computations
	decomp <- wlm::decompose_varcov(V, method = "evd", output = "all")
	remove(V)
	#############
	
	# CREATE OUTPUT TIBBLE
	#############
	snp_out <- tibble(
		TEST_SNP = numeric(),
		LOD_VALUE = numeric()
	)
	#############
	
	# OPEN CONNECTION TO SNP_TABLE FILE
	#############
	# open connection using LaF package
	snp_table.laf <- LaF::laf_open(detect_dm_csv(filename = snp_table_address, sep = "\t", header = FALSE))
	#############
	
	# MAIN FOR LOOP
	#############
	snp_nrow = wc_l(snp_table_address) # wc_l is one of my own functions that takes the number of lines in a given file (similar to wc -l in linux)
	rows_in_chunk = as.numeric(arguments["rows_in_chunk"]) # numbers of rows I scan in a "chunk"
	batchsize = as.numeric(arguments["batch_size"]) # number of rows in a batch for the matlm::matlm function
	
	# scan the SNP table in chunks, since all in one go uses too much memory to import the entire SNP table
	for (start in seq(from = 2, to = snp_nrow, by = rows_in_chunk)) { # start at 2 since the line 1 is a header
		
		# SET UP CURRENT BLOCK
		#############
		# move cursor to the starting row
		LaF::goto(x = snp_table.laf,
							i = start)
		
		# set up the current block of snp_table
		snp_table <- LaF::next_block(x = snp_table.laf, nrows = rows_in_chunk) %>%
			as_tibble() %>%
			rename("POS" = 1) %>%
			mutate(POS = as.numeric(POS),
						 SUM = rowSums(x = .[, -1])) %>% # make sure rows are valid; start by summing up the values of the rows (except for the POS column)
			filter(SUM != 0, SUM != (length(individuals) * 2)) %>% # filter out invalid rows (i.e. rows whose values are all 0s or 2s)
			select(-SUM) # sum column no longer needed, only use is for filtering out invalid snps
		
		#############
		
		# PREDICTOR TIBBLE
		#############
		# transpose snp_table so it can be fed to matlm function
		predictor <- snp_table %>%
			select(all_of(individuals)) %>%
			t()
		#############
		
		# APPROXIMATE P VALUES AND CONVERT TO LOD
		#############
		# now this is the fast approximation...
		# note it expects the matrix of SNPs to test as columns not rows, hence the "t",
		# you could do this outside the function of course
		# you would typically fit the "fixed" effects under the null here (i.e., sex, weight, etc),
		# but our null has only random effects
		passoc_gls <- matlm::matlm(formula = Y ~ 1,
															 data = phenotype,
															 pred = predictor,
															 ids = INDIVIDUAL,
															 transform = decomp$transform,
															 batch_size = batchsize,
															 verbose = 2)
		# matlm::matlm outputs p-values.
		
		# It is so bleeding edge I am not even exactly sure what is outputted.
		# Likely p-values or natural log p-values,
		# so I would have to look at the table to convert to LOD scores (so we can compare with the other two methods).
		# I am hoping I can just look at the table and tell.
		
		# convert p-values to LOD
		LODs <- -log(passoc_gls$tab$pval) / log(10) # 0 values become NaN because log_10(0) = infinity
		remove(passoc_gls)
		# LOD score is approximately -log(p)/log(10).
		
		#############
		
		# CREATE LOCAL TIBBLE TO ADD TO OUTPUT TIBBLE
		#############
		# local tibble (i.e. it is local to this iteration of the for-loop)
		snp_current <- tibble(
			TEST_SNP = snp_table$POS,
			LOD_VALUE = LODs)
		#############
		
		# APPEND CURRENT TIBBLE WITH OUTPUT TIBBLE
		#############
		# append the local tibble to the output tibble
		snp_out <- bind_rows(snp_out, snp_current)
		remove(snp_current)
		#############
	}
	
	# EXPORT SNP OUT
	#############
	snp_out <- snp_out %>%
		mutate(LOD_METHOD = rep(1, times = nrow(.))) %>% # LOD_METHOD refers to the preciseness of the method used to calculate the LOD Score. this method isnt that precise and is more of an approximation, and the approximation method is LOD_METHOD=1
		select(TEST_SNP, LOD_METHOD, LOD_VALUE)
	
	# snp_based_comment gives important file information, like the statistical test or genetic model
	snp_based_comment = paste("#", scan_info, "calculation_method=snp_based", sep = "\t")
	# write first line of the output (which is a comment line), and overwrite any previous files
	write_lines(x = snp_based_comment,
							path = snp_out_path,
							sep = "\n", # since the comment is one string, the sep needs to be \n so the line ends with a newline
							append = FALSE) # so a new file is created
	
	# append actual data following the comment line
	write_tsv(x = snp_out,
						path = snp_out_path,
						col_names = TRUE, # include column names
						append = TRUE) # append to the file info
	
	
	# add this filename to the completed file list (so I know it is done)
	write_lines(x = snp_out_path,
							path = completed_tests_path,
							sep = "\n", # since the path is one string, the sep needs to be \n so the line ends with a newline
							append = TRUE) # so that it is added to the completed_tests_path file
	#############
}
if (file.exists(snp_out_path)) {
	snp_already_complete = TRUE
}
####################################################



####################################################
# HAPLOTYPE
hap_out_path = trimws(paste(as.character(arguments["hap_out"]), paste(paste("Chr", csnp_chromosome_number, "_", test_type, "_hapbased",  sep = ""), causative_snp, "tsv", sep = "."), sep = "/"))
hap_nrow = wc_l(hap_out_path) # count the number of lines already in the haplotype-based output

# GET USEFUL DATA FOR HAPLOTYPE TEST
#############
hap_table_address = as.character(arguments["hap_table_address"])

# use base-R's pipe() to extract the first line of the hap_table to get column names as a vector
hap_table_colnames_connection = pipe(paste("head -n 1", hap_table_address, sep = " "))
hap_table_colnames <- readLines(con = hap_table_colnames_connection, n = 1) %>%
	strsplit(x = ., split = "\t") %>% unlist()
close(con = hap_table_colnames_connection)
remove(hap_table_colnames_connection)

# extract the loci in the hap_table
loci_from_hap_table <- read_tsv(file = hap_table_address, col_names = TRUE, col_types = cols_only(POS = 'i')) %>%
	pull(POS)
#############

# run if the number of lines in the haplotype output is not the full amount it should be
if (hap_nrow < ((length(loci_from_hap_table) * 2) + 2)) { # multiply by 2 since there are 2 LOD values per locus in output table, add 2 to account for the comment and header, checking if the haptest is complete
	
	hap_based_comment = paste("#", scan_info, "calculation_method=hap_based", sep = "\t")
	
	# CHECK IF COMMENT AND HEADER NEEDS TO BE PRINTED
	#############
	if (hap_nrow == 0) { # if file is empty and file needs to be created with proper starting lines
		# write the hap_based_comment, which gives vital file information, like the causative SNP, the statistical test, and genetic model
		write_lines(x = hap_based_comment, path = hap_out_path, sep = "\n", append = FALSE)
		# write the column names
		write_lines(x = paste("TEST_SNP", "LOD_METHOD", "LOD_VALUE", sep = "\t"), path = hap_out_path, sep = "\n", append = TRUE)
		
		hap_nrow = 2 # because there are now two rows at the beginning
	}
	#############
	
	# GFIT NULL MODEL
	#############
	# use first 15 Principal Component Scores of the kinship matrix as an approximation
	# then I use the much much faster linear model
	# PCSs are automatically sorted from biggest to smallest
	reduced_matrix_pcs <- stats::princomp(covmat = kinship_matrix)
	yr <- (lm(formula = Y ~ (reduced_matrix_pcs$loadings[, 1:15]), data = phenotype))$resid
	remove(reduced_matrix_pcs)
	
	# fit the model with only the causative SNP
	gfit_null_hap <- lme4qtl::relmatLmer(formula = Y ~ (1 | INDIVIDUAL),
																			 data = phenotype,
																			 relmat = list(INDIVIDUAL = kinship_matrix),
																			 REML = FALSE)
	#############
	
	# FOR LOOP
	#############
	drop_PC_less_than = as.numeric(arguments["drop_PC_less_than"]) # drop principal component scores less than (generally 0.05, but is an argument that the program takes)
	LOD_threshold = as.numeric(arguments["LOD_threshold"]) # the threshold needed for a LOD score to surpass to run the second, more accurate, but slower test
	first_iter_value = (hap_nrow - 2) %/% 2 # the value of the first iteration of the for loop (determined by how complete the output table currently is)
	for (i in seq(from = first_iter_value, to = length(loci_from_hap_table), by = 1)) { # (hap_nrow - 2 + 1) because -2 to account for comment and header and +1 so that we start scanning on the one after we stopped at
		is_first_iter = (i == first_iter_value)
		
		# IF IT IS THE VERY FIRST ITERATION
		#############
		# WRITE OUT THIS SPECIFIC TEST INFO ON FIRST ITERATION
		if (is_first_iter) cat(paste(hap_based_comment, "\n", sep = "")) # print out test information on first info and print a newline character so it looks good (this is more for me, this is not printed to the output)
		
		# SKIP CURRENT ITERATION IF THE FINAL TEST_SNP OF CURRENT HAP_OUT WAS COMPLETE
		# if this is the first iteration, and the ending locus of the current hap out file already has 2 values (LOD_1 and LOD_2), skip the current iteration and move onto the next
		if (is_first_iter & (((hap_nrow - 2) / 2) %% 1 == 0)) next()
		#############
		
		# LOAD HAPLOTYPE TABLE LINES
		#############
		# use base-R's pipe() to extract lines from the hap_table one at a time
		hap_table_connection = pipe(paste("head -n", format(x = (i + 1), scientific = FALSE, big.mark = ""), hap_table_address, "| tail -n 1", sep = " "))
		hap_table <- readLines(con = hap_table_connection, n = 1) %>%
			strsplit(x = ., split = "\t") %>% unlist() %>%
			as.numeric() %>%
			as_tibble_row(.name_repair = "minimal") %>%
			set_colnames(hap_table_colnames) %>%
			select(POS, ends_with(individuals)) %>% # filter down to the columns of the individuals that we want
			pivot_longer(cols = names(.)[names(.) != "POS"],
									 names_to = "TYPE",
									 values_to = "VALUE") %>%
			separate(col = TYPE,
							 into = c("DOSE", "INDIVIDUAL"),
							 sep = "-") %>%
			pivot_wider(names_from = DOSE,
									values_from = VALUE) %>%
			pivot_wider(names_from = POS,
									values_from = paste("HAP", 1:8, sep = "_")) %>%
			select(-INDIVIDUAL)
		
		close(con = hap_table_connection)
		remove(hap_table_connection)
		
		#############
		
		# ACTUALLY CALCULATE THE TWO LOD SCORES
		#############
		
		# I take the haplotype scores and create something called pricipal component scores
		# these are a transformation of the raw haplotypes "rotated" so that they are uncorrelated
		# the variance due to these uncorrelated variables should be the same, but the modeling is easier
		# I drop columns of the 8 column rotated scores that account for less than 5% of total variation
		# this should make the model fitting more robust as well
		
		# principal component score (pcs)
		pcs <- stats::prcomp(x = hap_table)
		pcs_filter <- pcs$x[, ((pcs$sdev^2 / sum(pcs$sdev^2)) > drop_PC_less_than)]
		
		# so initially fit something really simple
		LOD_1 = -log(x = (stats::anova(lm(yr ~ pcs_filter))$`Pr(>F)`[1])) / log(x = 10) # convert p-values to LOD -log(p) / log(10) where p = anova(lm(yr ~ pcs_filter))$`Pr(>F)`[1]
		LOD_2 = "NA" # set the second LOD value to NA for now
		
		# if the first LOD score is >2 (which should happen about 1% of time since a 10^-2 = 0.01)
		# fit the full mixed model ... if too many haps are >2 ... we could set the threshold at 3, etc.
		# because the full mixed model takes significantly longer, we only fit it if the more approximate test (used to generate the first LOD score)
		# yields a value greater than a given threshold
		if (LOD_1 > LOD_threshold){
			mixed_model_gfit <- stats::update(object = gfit_null_hap,
																				formula = . ~ . + pcs_filter) # fit the mixed model from the null model generated earlier
			mixed_model_pcs <- stats::anova(mixed_model_gfit, gfit_null_hap)
			
			LOD_2 = -(stats::pchisq(q = mixed_model_pcs[2, 6],
															df = mixed_model_pcs[2, 7],
															lower.tail = FALSE,
															log.p = TRUE)) / log(10) # convert to LOD score
		}
		#############
		
		# WRITE OUT INFO AND DATA ACCORDING TO CERTAIN CONDITIONS
		#############
		
		# ALWAYS WRITE FIRST LINE (THE LOD_1 VALUE) UNLESS THE FINAL LOCUS OF CURRENT HAP_OUT WAS INCOMPLETE
		# paste(TEST_SNP, LOD_METHOD, LOD_VALUE, sep = "\t")
		# write the LOD_1
		if (!is_first_iter | (is_first_iter & (((hap_nrow - 2) / 2) %% 1 != 0.5))) { # if it is the first iteration and the value is 0.5, this is the only case where no line for LOD_1 should be printed, since it was already printed in the running of this program before
			write_lines(x = paste(loci_from_hap_table[[i]], 1, LOD_1, sep = "\t"),
									path = hap_out_path,
									sep = "\n",
									append = TRUE)
		}
		
		# ALWAYS WRITE THE SECOND LINE (THE LOD_2 VALUE)
		# write the LOD_2 line
		write_lines(x = paste(loci_from_hap_table[[i]], 2, LOD_2, sep = "\t"),
								path = hap_out_path,
								sep = "\n",
								append = TRUE)
		#############
		
		# PROGRESS UPDATE FOR ME
		#############
		# print out percent scanned every 5% to indicate progress
		if (i %in% as.integer(seq(from = 0.05 * length(loci_from_hap_table),
															to = length(loci_from_hap_table),
															length.out = 20))) {
			cat(paste(round(x = (i / length(loci_from_hap_table)) * 100, digits = 0), # percentage scanned
								"% scanned, currently scanning locus ",
								i, # how many loci have actually been scanned
								".\n",
								sep = ""))
		}																																																																									 
		#############
		
	}
	#############
	
	# ADD THE HAP_OUT FILE PATH TO THE COMPLETED FILES LIST AND REEVALUATE THE LINE COUNT OF HAP_NROW
	#############
	# this will only be run if the for loop for hapbased test is completed, since if the program gets cancelled while running the for loop, this never gets reached
	# add this filename to the completed file list so I know this scan is complete
	write_lines(x = hap_out_path,
							path = completed_tests_path,
							sep = "\n", # since the path is one string, the sep needs to be \n so the line ends with a newline
							append = TRUE) # so that it is added to the completed_tests_path file
	
	hap_nrow = wc_l(hap_out_path) # update the line count of the hap_out file
	#############
}
if (hap_nrow >= ((length(loci_from_hap_table) * 2) + 2)){
	hap_already_complete = TRUE
}
####################################################

if (snp_already_complete & hap_already_complete) {
	# so if both the snp and hap tests were completed in a previous running of this code or this running of the code
	cat("hapsnp complete.\nWBqVU4pz\n") # random sequence that makes filtering out completed slurm files easy via finding slurm files with this sequence in them
}
# grep -H "WBqVU4pz" /share/adl/pnlong/mouseproject/slurm_files/*.out | cut -d : -f 1 | xargs rm

# End of Program.
# [{-_-}]
# ZZZzz zz z...
