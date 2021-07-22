# Peromyscus_power_bleed
Code associated with the peromyscus power / bleed-time paper (Long et al 2021)

## Kinship matrix
*generate_Mjj.sh* generates a kinship matrix.

## Simulating phenotypes
*generate_Y.sh* generates simulated phenotypes for 3 different genetic models at a causative locus.

## Scans
*scan.sh* performs both a marker- and haplotype-based scan of a chromosome given a phenotype.

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
