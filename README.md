# Peromyscus_power_bleed
Code associated with the peromyscus power / bleed-time paper (Long et al 2021)

## Prerequisites
Some files and directories needed include:
- a `save_data` directory, where important files like the kinship matrix or simulated phenotypes are stored
- a `software` directory, download [here](https://github.com/tdlong/Peromyscus_power_bleed/tree/main/software)
- a `gzipped_vcf_files_list` file, a list of gzipped filepaths to *.vcf* files containing STITCH data for a chromosome, one filepath per line (unzip *gzipped_vcf_files.txt.gz* in the *.tar* folder for an example)
- an `individuals_full` file, a tab-delimited list of all individual IDs used by STITCH in imputing genotypes (unzip and use *individuals/all_individuals.tsv.gz* in the *.tar* folder)
- an `individuals_subset` file, a tab-delimited list of the subset of 297 individual IDs used throughout the paper (unzip and use *individuals/individuals.tsv.gz* in the *.tar* folder)

## Kinship matrix
*`generate_Mjj.sh`* generates a kinship matrix, provided the following arguments:
1. path to `save_data` directory
2. path to `software` directory
3. path to **unzipped** `gzipped_vcf_files_list`
4. path to **unzipped** `individuals_subset`
5. path to **unzipped** `individuals_full`



## Simulating phenotypes
*`generate_Y.sh`* generates simulated phenotypes for 3 different genetic models at a given causative SNP, provided the following arguments:
1. path to `save_data` directory
2. path to `software` directory
3. path to **unzipped** `gzipped_vcf_files_list`
4. path to gzipped *.vcf* file of STITCH imputation data for the causative SNP's chromosome
5. position of causative SNP (note that if the given SNP is not in the gzipped *.vcf* file, the closest valid SNP is used)
6. path to **unzipped** `individuals_subset`
7. path to **unzipped** `individuals_full`



## Scans
*`scan.sh`* performs both a marker- and haplotype-based scan of a chromosome given a phenotype, provided the following arguments:
1. path to `save_data` directory
2. path to `software` directory
3. path to gzipped *.vcf* file of STITCH imputation data for the chromosome that will be scanned
4. path to **unzipped** `individuals_subset`
5. path to **unzipped** `individuals_full`
6. path to phenotype file (containing at least an **INDIVIDUAL** and phenotype column, often called **Y**)
7. name of the column with phenotype data in the phenotype file of the previous argument
8. path to results directory (subdirectories `snpbased` and `hapbased` will be created for the outputs of marker- and haplotype-based scans, respectively)
9. path to completed scans file where a scan's filepath is printed once that scan is complete
10. genetic model (either "single" [single SNP causative], "multiple" [10 SNPs causative], "rare" [all SNPs causative], or "unknown")



## Bleeding time data (and normalization)
```R
k5 <- read.table(file = "k5.txt")

# normalize Bleeding time
hist(k5$BleedTime)
oo = rank(k5$BleedTime)/(397 - 1)/(2 * 397)
k5$normBleed = qnorm(oo)
plot(x = k5$normBleed, y = k5$BleedTime)

# check for fixed factors affecting bleeding time
anova(lm(formula = normBleed ~ Sex + Age + Weight + timeDOB,
         data = k5))
k5$residNormBleed <- lm(formula = normBleed ~ Sex + Age + timeDOB,
                        data = k5)$resid
write.table(x = k5,
            file = "k5.txt")
```
```
Analysis of Variance Table

Response: normBleed
           Df  Sum Sq Mean Sq F value    Pr(>F)    
Sex         1   4.919   4.919  6.2166 0.0130673 *  
Age         1  10.073  10.073 12.7291 0.0004046 ***
Weight      1   1.385   1.385  1.7503 0.1866126    
timeDOB     1  32.254  32.254 40.7589 4.874e-10 ***
Residuals 392 310.202   0.791                      
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
