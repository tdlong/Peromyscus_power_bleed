# Peromyscus_power_bleed
Code associated with the peromyscus power / bleed-time paper (Long et al 2021)

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

## next section
