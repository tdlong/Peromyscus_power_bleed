# Peromyscus_power_bleed
Code associated with the peromyscus power / bleed-time paper (Long et al 2021)

## Genetic background variance
- *genetic_background_variance.py*
- *genetic_background_variance.sh*

## Kinship matrix
- *Mjj_average.py*
- *Mjj_calculation.sh*
- *Mjj_combine.sh*
- *Mjj_setup.sh*
- *Mjj_sum.py*

## Simulation scans
- *hapsnp_calculation.R*
- *hapsnp_calculation.sh*
- *hapsnp_calculation_save.R*
- *hapsnp_tables_setup.sh*
- *hapsnp_tables_setup_manipulate.py*
- *hapsnp_tests.sh*

## Bleeding time data (and normalization)
```R
k5=read.table("k5.txt")

# normalize Bleeding time
hist(k5$BleedTime)
oo = rank(k5$BleedTime)/397 - 1/(2*397)
k5$normBleed = qnorm(oo)
plot(k5$normBleed,k5$BleedTime)

# check for fixed factors affecting bleeding time
anova(lm(normBleed~Sex+Age+Weight+timeDOB,data=k5))
k5$residNormBleed = lm(normBleed~Sex+Age+timeDOB,data=k5)$resid
write.table(k5,"k5.txt")
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

## Bleeding time scans
- *demonstration_chromosome_scans.R*
- *demonstration_manhattan.R*
- *demonstration_qq.R*
- *demonstration_scan.R*
- *demonstration_scan.sh*
- *demonstration_setup.sh*
- *dosage_association.R*

## Figures and Tables
- *stitch_validation.R*
- *stitch_validation_setup.py*
- *stitch_validation_setup.sh*
- *summary_plots.R*
- *summary_significance_threshold_qq.R*
- *summary_statistics.R*
- *summary_statistics_plots.R*
- *figures.sh*

## Other
- *useful_functions.R*
