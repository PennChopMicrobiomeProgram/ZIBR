# ZIBR (Zero-Inflated Beta Random Effect model)

<!-- badges: start -->
  [![R-CMD-check](https://github.com/PennChopMicrobiomeProgram/ZIBR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PennChopMicrobiomeProgram/ZIBR/actions/workflows/R-CMD-check.yaml)
  [![codecov](https://codecov.io/gh/PennChopMicrobiomeProgram/ZIBR/graph/badge.svg?token=6A7MIF2IPE)]( https://app.codecov.io/gh/PennChopMicrobiomeProgram/ZIBR)
  [![Published in Bioinformatics](https://img.shields.io/badge/Published%20in-Bioinformatics-blue)](https://academic.oup.com/bioinformatics/article-abstract/32/17/2611/2450750?rss=1)
  <!-- badges: end -->

## Introduction
The longitudinal microbiome compositional data are highly skewed, bounded in [0,1), and often sparse with many zeros. In addition, the observations from repeated measures are correlated. We propose a two-part zero-inflated Beta regression model with random effects (ZIBR) for testing the association between microbial abundance and clinical covariates for longitudinal microbiome data. The model includes a logistic component to model presence/absence of the microbe in samples and a Beta component to model non-zero microbial abundance and each component includes a random effect to take into account the correlation among repeated measurements on the same subject.

The details of the statistical model are as follows:
<img src="vignettes/zibr_stat_model.png" width="600" align="center">

The ZIBR model combines the logistic regression and Beta regression in one model. Each regression part includes random effects to account for correlations acorss time points. We call these two regressions in ZIBR model as logistic component and Beta component. These two components model two different aspects of the data. The logistic component models presence/absence of the microbe and Beta component models non-zero microbial abundance.

Accordingly, we can test three biologically relevant null hypotheses:  
- H0: α_j = 0.  This is to test the coefficients in the logistic component, if the covariates are associated with the bacterial taxon by affecting its presence or absence;  
- H0: β_j = 0.  This is to test the coefficients in the Beta component, if the taxon is associated with the covariates by showing different abundances;  
- H0: α_j = 0 and β_j = 0 for each covariate X_j and Z_j. This is to joinly test the coefficients in both logistic and Beta components, if the covariates affect the taxon both in terms of presence/absence and its non-zero abundance.  

## Installation
You can install our ZIBR package from CRAN

```r
install.packages("ZIBR")
```

Or get the dev version from GitHub

```r
#install.packages("devtools")
devtools::install_github("PennChopMicrobiomeProgram/ZIBR")
library(ZIBR)
```

## Basic Usage

```r
zibr_fit <- zibr(logistic_cov=logistic_cov,beta_cov=beta_cov,Y=Y,subject_ind=subject_ind,time_ind=time_ind)
```

- **logistic_cov**: covariates for the logistic component. Rows: samples. Columns: covariates.  
- **beta_cov**: covariates for the Beta component. Rows: samples. Columns: covariates.  
- **Y**: the response variable (i.e the bacterial relative abundance). It is a vector with values in [0,1).  
- **subject_ind**: the variable with subject IDs.   
- **time_ind**: the variable with time points.   

The ordering of the samples in the above matrix or vectors must be consistent.

The zibr function will return the following results

```r
zibr_fit
```

- **logistic_est_table**: the estimated coefficients for logistic component.  
- **logistic_s1_est**: the estimated standard deviation for the random effect in the logistic component.  
- **beta_est_table**: the estimated coefficients for Beta component.  
- **beta_s2_est**: the estimated standard deviation for the random effect in the Beta component.  
- **beta_v_est**: the estimated dispersion parameter in the Beta component.  
- **joint_p**: the p-values for jointly testing each covariate in both logistic and Beta component.  

## Troubleshooting

### Missing values
If there are missing values in certain time points, they can be imputed as following:
1. Calculate the mean or median of values from previous time point(s) and later time points(s). Use such values to replace the missing values.
2. Group the time point with missing values with other time points. For example, if you have T1, T2, T3 and T4 and T1 has missing values, you can group T1 and T2 as one time point.

After the missing values are imputed, the data can be fed into ZIBR.

### Other issues

If you have other problems with the package or features/fixes to suggest, please open an issue on the [GitHub issues page](https://github.com/PennChopMicrobiomeProgram/ZIBR/issues).

## Citation
Eric Z. Chen and Hongzhe Li (2016). A two-part mixed effect model for analyzing longitudinal microbiome data. Bioinformatics. [Link](https://academic.oup.com/bioinformatics/article-abstract/32/17/2611/2450750?rss=1)

## Contact
Maintained by Charlie Bushman (ctbushman\@gmail.com) and the Penn CHOP Microbiome Program.

## Updates

<!---
- Add likelihood to the output.
- variable name is missing in beta.est.table
-->
