# README
# This program will calculate LOD scores for the loci of a given chromosome
# using both a haplotype- and SNP-based calculation method.
# Phillip Long
# June 16, 2021
' .sh runner
####
module load R/3.6.2
R --vanilla -f /share/adl/pnlong/mouseproject/software/demonstration_scan.R --args phenotype_address="/share/adl/pnlong/mouseproject/demonstration_trait/phenotype.tsv" phenotype="BT1" shuffle_phenotype="FALSE" snp_table_address="/share/adl/pnlong/mouseproject/save_data/snp_table_Chr23.tsv" hap_table_address="/share/adl/pnlong/mouseproject/save_data/hap_table_Chr23.tsv" kinship_address="/share/adl/pnlong/mouseproject/save_data/Mjj_values_all_pnl.tsv" batch_size=100 rows_in_chunk=5000 drop_PC_less_than=0.05 LOD_threshold=2 snp_out=$snp_based_out hap_out=$hap_based_out completed_files=$completed_files_list
####
'
####################################################


####################################################
# INGREDIENTS FOR BOTH

# LOADS
#############
library(tidyverse)
library(magrittr)
library(lmtest)
library(lme4qtl)
library(matlm) 
library(wlm)
library(Matrix)
library(LaF)

load(file = "/share/adl/pnlong/mouseproject/software/useful_functions.Rsave")
# FOR TESTING
# arguments <- c(phenotype_address="/share/adl/pnlong/mouseproject/demonstration_trait/phenotype.tsv", phenotype="BT1", shuffle_phenotype="FALSE", snp_table_address="/share/adl/pnlong/mouseproject/save_data/snp_table_Chr9.tsv", hap_table_address="/share/adl/pnlong/mouseproject/save_data/hap_table_Chr9.tsv", kinship_address="/share/adl/pnlong/mouseproject/save_data/Mjj_values_all_pnl.tsv", batch_size="100", rows_in_chunk="5000", drop_PC_less_than="0.05", LOD_threshold="2", snp_out="/share/adl/pnlong/mouseproject/demonstration_trait/BT1/snpbased", hap_out="/share/adl/pnlong/mouseproject/demonstration_trait/BT1/hapbased", completed_files="/share/adl/pnlong/mouseproject/demonstration_trait/BT1/completed_tests.txt")   
arguments <- get_args()
completed_tests_path = as.character(arguments["completed_files"])
#############

# INDIVIDUALS
##############
individuals_file = file("/share/adl/pnlong/mouseproject/save_data/individuals.tsv")
individuals <- readLines(con = individuals_file, n = -1) %>%
	strsplit(x = ., split = "\t") %>% unlist()
close(individuals_file)
##############

# PHENOTYPE (BG)
#############
y_column = as.character(arguments["phenotype"])
phenotype <- read_tsv(file = as.character(arguments["phenotype_address"]), col_names = TRUE, col_types = cols(col_number(), col_character(), col_number(), col_number())) %>%
	select(ID, y_column) %>%
	rename("INDIVIDUAL" = ID, "Y" = y_column) %>%
	filter(INDIVIDUAL %in% individuals) %>%
	arrange(INDIVIDUAL)

nonduplicate_individuals <- phenotype %>%
	group_by(INDIVIDUAL) %>%
	summarise(COUNT = n()) %>%
	filter(COUNT == 1) %>% # remove duplicated individuals
	pull(INDIVIDUAL)

# set the individuals list to the individuals in the phenotype, so that the duplicate individual(s) are not included
individuals <- nonduplicate_individuals
phenotype %<>%
	filter(INDIVIDUAL %in% nonduplicate_individuals)
remove(nonduplicate_individuals)

# # FOR SHUFFLED Y VALUES CONTROL
if (as.logical(arguments["shuffle_phenotype"]) == TRUE) {
	set.seed(0.001) # just to make it reproducible
	phenotype %<>%
		mutate(Y = sample(Y))
}
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

# SCAN CHROMSOME INFO
#############
snp_table_address = as.character(arguments["snp_table_address"])

tsnp_chrom = strsplit(x = snp_table_address, split = "/") %>% unlist() %>%
	extract(length(.)) %>%
	strsplit(x = ., split = "\\.") %>% unlist() %>%
	extract(1) %>%
	strsplit(x = ., split = "_") %>% unlist() %>%
	extract(length(.)) %>%
	str_remove(string = ., pattern = "Chr")

tsnp_chrom_label = paste("tSNP_chromosome", tsnp_chrom,
												 sep = "=")
#############

####################################################



####################################################
# MARKER
snp_out_path = paste(as.character(arguments["snp_out"]), "/Chr", tsnp_chrom, "_snpbased.tsv", sep = "") # address of the output file for the marker-based scan

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
	# open connection
	snp_table.laf <- LaF::laf_open(detect_dm_csv(filename = snp_table_address, sep = "\t", header = TRUE))
	#############
	
	# MAIN FOR LOOP
	#############
	
	# get SNP table column names
	snp_table_colnames <- names(snp_table.laf) %>%
		str_remove(string = ., pattern = "X")

	snp_nrow = wc_l(snp_table_address) # wc -l /share/adl/pnlong/data/SNP_table_Chr23_certain_individuals.tsv
	rows_in_chunk = as.numeric(arguments["rows_in_chunk"])
	batchsize = as.numeric(arguments["batch_size"])
	# scan the SNP table in chunks, since all in one go uses too much memory to import the entire SNP table
	for (start in seq(from = 1, to = snp_nrow, by = rows_in_chunk)) { # start at 2 since the line 1 is a header
		
		# SET UP CURRENT BLOCK
		#############
		# move cursor to the starting row
		LaF::goto(x = snp_table.laf,
							i = start)
		
		# set the current block of snp_table
		snp_table <- LaF::next_block(x = snp_table.laf, nrows = rows_in_chunk) %>%
			as_tibble() %>%
			set_colnames(snp_table_colnames) %>%
			select(POS, all_of(individuals)) %>%
			mutate(POS = as.numeric(POS),
						 SUM = rowSums(x = .[, -1])) %>% # make sure rows are valid; start by summing up the values of the rows (except for the POS column)
			filter(SUM != 0, SUM != (length(individuals) * 2)) %>% # filter out non valid rows
			select(-SUM) # sum column no longer needed, only use is for filtering
		
		#############
		
		# PREDICTOR TIBBLE
		#############
		# unnest and transpose snp_table so it can be fed to matlm function
		predictor <- snp_table %>%
			select(-POS) %>%
			t()
		#############
		
		# APPROXIMATE P VALUES AND CONVERT TO LOD
		#############
		# now this is the fast approximation…
		# note it expects the matrix of SNPs to test as columns not rows, hence the "t",
		# you could do this outside the function of course
		# I am not sure if the "1" below should be a "0", I think a one
		# you would typically fit the "fixed" effects under the null here (i.e., sex, weight, etc),
		# but our null has only random effects.
		passoc_gls <- matlm::matlm(formula = Y ~ 1,
															 data = phenotype,
															 pred = predictor,
															 ids = INDIVIDUAL,
															 transform = decomp$transform,
															 batch_size = batchsize,
															 verbose = 2)
		
		# It is so bleeding edge I am not even exactly sure what is outputted.
		# Likely p-values or natural log p-values,
		# so I would have to look at the table to convert to LOD scores (so we can compare with the other two methods).
		# I am hoping I can just look at the table and tell.
		LODs <- -log(passoc_gls$tab$pval) / log(10) # 0 values become NaN because log_10(0) = infinity
		remove(passoc_gls)
		# likely has a column called p-values, or ln_p of log10_p or something like that.
		# ln usually mean log base e, log10 is usually log base 10, but “log” can mean either,
		# so the worst case would be a column call log_p then we don’t know if it natural log or log bae 10 (as both are commonly used).
		# If I look at the first 10 observations of the tab compared to the first 10 LOD scores using the other two methods,
		# we can likely figure it out (as log10 values are just natural logs divided by ln(10)).
		# Looking at the code it may be outputting straight p-values.
		# If so the LOD score is approximately -log(p)/log(10).
		# The goal would be to figure out what is actually output, convert to a LOD, and compare to the other two {slower?} methods via an adopt plot.
		
		#############
		
		# CREATE TIBBLE TO ADD TO OUTPUT TIBBLE LATER
		#############
		snp_current <- tibble(
			TEST_SNP = snp_table$POS,
			LOD_VALUE = LODs)
		#############
		
		# APPEND CURRENT TIBBLE WITH OUTPUT TIBBLE
		#############
		snp_out <- bind_rows(snp_out, snp_current)
		remove(snp_current)
		#############
	}
	
	# EXPORT SNP OUT
	#############
	snp_out <- snp_out %>%
		mutate(LOD_METHOD = rep(1, nrow(.))) %>% # LOD_METHOD refers to the preciseness of the method used to calculate the LOD Score. this method isnt that precise and is more of an approximation
		select(TEST_SNP, LOD_METHOD, LOD_VALUE)
	
	snp_based_comment = paste("#", tsnp_chrom_label, "calculation_method=snp_based", sep = "\t")
	# write first line of the output (which is a comment line), and overwrite any previous files
	write_lines(x = snp_based_comment,
							path = snp_out_path,
							sep = "\n", # since the comment is one string, the sep needs to be \n so the line ends with a newline
							append = FALSE) # so a new file is created
	# append actual data following the comment line
	write_tsv(x = snp_out,
						path = snp_out_path,
						col_names = TRUE,
						append = TRUE) # append to the file info
	
	
	# add this filename to the completed file list
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
hap_out_path = paste(as.character(arguments["hap_out"]), "/Chr", tsnp_chrom, "_hapbased.tsv", sep = "")
hap_nrow = wc_l(hap_out_path)

# GET USEFUL DATA FOR HAPLOTYPE TEST
#############
hap_table_address = as.character(arguments["hap_table_address"])

hap_table_colnames_connection = pipe(paste("head -n 1", hap_table_address, sep = " "))
hap_table_colnames <- readLines(con = hap_table_colnames_connection, n = 1) %>%
	strsplit(x = ., split = "\t") %>% unlist()
close(con = hap_table_colnames_connection)
remove(hap_table_colnames_connection)

loci_from_hap_table <- read_tsv(file = hap_table_address, col_names = TRUE, col_types = cols_only(POS = 'i')) %>%
	pull(POS)
#############

# run if the number of lines in the haplotype output is not the full amount it should be
if (hap_nrow < ((length(loci_from_hap_table) * 2) + 2)) { # multiply by 2 since there are 2 LOD values per locus in hap_table, add 2 to account for the comment and header, checking if the haptest is complete
	
	hap_based_comment = paste("#", tsnp_chrom_label, "calculation_method=hap_based", sep = "\t")
	
	# CHECK IF COMMENT AND HEADER NEEDS TO BE PRINTED
	#############
	if (hap_nrow == 0) { # if file is empty and file needs to be created with proper starting lines
		write_lines(x = hap_based_comment,
								path = hap_out_path,
								sep = "\n",
								append = FALSE)
		write_lines(x = paste("TEST_SNP", "LOD_METHOD", "LOD_VALUE", sep = "\t"),
								path = hap_out_path,
								sep = "\n",
								append = TRUE)
		
		hap_nrow = 2 # because there are now two rows at the beginning
	}
	#############
	
	# GFIT NULL MODEL
	#############
	# use first 15 PCs of the kinship matrix as an approximation
	# then I use the much much faster linear model
	# PCs are automatically sorted from biggest to smallest
	reduced_matrix_pcs <- stats::princomp(covmat = kinship_matrix)
	yr <- (lm(formula = Y ~ (reduced_matrix_pcs$loadings[, 1:15]), data = phenotype))$resid
	remove(reduced_matrix_pcs)
	
	# fit the model with only the causative SNP
	# gfit_null_hap <- relmatLmer(formula = (phenotype$Y) ~ (1 | individuals),
	#                             relmat = list(individuals = kinship_matrix),
	#                             REML = FALSE)
	gfit_null_hap <- lme4qtl::relmatLmer(formula = Y ~ (1 | INDIVIDUAL),
															data = phenotype,
															relmat = list(INDIVIDUAL = kinship_matrix),
															REML = FALSE)
	#############
	
	# FOR LOOP
	#############
	drop_PC_less_than = as.numeric(arguments["drop_PC_less_than"])
	LOD_threshold = as.numeric(arguments["LOD_threshold"])
	first_iter_value = (hap_nrow - 2) %/% 2
	
	for (i in seq(from = first_iter_value, to = length(loci_from_hap_table), by = 1)) { # (hap_nrow - 2 + 1) because -2 to account for comment and header and +1 so that we start scanning on the one after we stopped at
		is_first_iter = (i == first_iter_value)
		
		# IF IT IS THE VERY FIRST ITERATION
		#############
		# WRITE OUT THIS SPECIFIC TEST INFO ON FIRST ITERATION
		if (is_first_iter) cat(paste(hap_based_comment, "\n", sep = "")) # print out test information on first info and print a newline character so it looks good (this is more for me, this is not printed to the output)
		
		# SKIP CURRENT ITERATION IF THE FINAL TEST_SNP OF CURRENT HAP_OUT WAS COMPLETE
		# if this is the first iteration, and the ending test_snp of the current hap out file already has 2 values, skip current iteration
		if (is_first_iter & (((hap_nrow - 2) / 2) %% 1 == 0)) next()
		#############
		
		# LOAD HAPLOTYPE TABLE LINES
		#############
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
		LOD_1 = -log(x = (stats::anova(lm(yr ~ pcs_filter))$`Pr(>F)`[1])) / log(x = 10)
		LOD_2 = "NA"
		
		# if the LOD score is >2 (which should happen about 1% of time since a 10^-2 = 0.01)
		# fit the full mixed model ... if too many haps are >2 ... we could set the threshold at 3, etc.
		if (LOD_1 > LOD_threshold){
			mixed_model_gfit <- stats::update(object = gfit_null_hap,
																				formula = . ~ . + pcs_filter)
			mixed_model_pcs <- stats::anova(mixed_model_gfit, gfit_null_hap)
			
			LOD_2 = -(stats::pchisq(q = mixed_model_pcs[2, 6],
															df = mixed_model_pcs[2, 7],
															lower.tail = FALSE,
															log.p = TRUE)) / log(10)
		}
		#############
		
		# WRITE OUT INFO AND DATA ACCORDING TO CERTAIN CONDITIONS
		#############
		
		# ALWAYS WRITE FIRST LINE (THE LOD_1 VALUE) UNLESS THE FINAL TEST_SNP OF CURRENT HAP_OUT WAS INCOMPLETE
		# paste(TEST_SNP, LOD_METHOD, LOD_VALUE, sep = "\t")
		# write the LOD_1
		if (!is_first_iter | (is_first_iter & (((hap_nrow - 2) / 2) %% 1 != 0.5))) { # if it is the first iteration and the value is 0.5, this is the only case where no line for LOD_1 should be printed, since it was already printed in the run before
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
		# print out values every once in a while to indicate progress
		if (i %in% vector_separate_at_indices(x = seq_along(loci_from_hap_table), n = 20)) cat(paste(round((i / length(loci_from_hap_table)) * 100), # percentage scanned
																																																 "% scanned, currently scanning locus ",
																																																 i, # how many loci have actually been scanned
																																																 ".\n",
																																																 sep = ""))
		#############
		
	}
	#############
	
	# ADD THE HAP_OUT FILE PATH TO THE COMPLETED FILES LIST AND REEVALUATE THE LINE COUNT OF HAP_NROW
	#############
	# this will only be run if the for loop for hapbased test is completed, since if the program gets cancelled while running the for loop, this never gets reached
	# add this filename to the completed file list
	write_lines(x = hap_out_path,
							path = completed_tests_path,
							sep = "\n", # since the path is one string, the sep needs to be \n so the line ends with a newline
							append = TRUE) # so that it is added to the completed_tests_path file
	
	hap_nrow = wc_l(hap_out_path)
	#############
}
if (hap_nrow >= ((length(loci_from_hap_table) * 2) + 2)){
	hap_already_complete = TRUE
}
####################################################

if (snp_already_complete & hap_already_complete) {
	# so if both the snp and hap tests were completed in a previous running of this code
	cat("scan complete.\nWBqVU4pz\n") # random sequence that makes filtering out completed slurm files easy via finding slurm files with this sequence in them
}
# grep -H "WBqVU4pz" /share/adl/pnlong/mouseproject/slurm_files/*.out | cut -d : -f 1 | xargs rm

# End of Program.
# [{-_-}]
# ZZZzz zz z...

