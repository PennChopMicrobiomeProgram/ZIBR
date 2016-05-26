# ZIBR (Zero-Inflated Beta Random Effect model)

## Introduction
The longitudinal microbiome compositional data are highly skewed, bounded in [0,1), and often sparse with many zeros. In addition, the observations from repeated measures are correlated. We propose a two-part zero-inflated Beta regression model with random effects (ZIBR) for testing the association between microbial abundance and clinical covariates for longitudinal microbiome data. The model includes a logistic component to model presence/absence of the microbe in samples and a Beta component to model non-zero microbial abundance and each component includes a random effect to take into account the correlation among repeated measurements on the same subject.

The details of the statistical model are as follows:
<img src="inst/image/zibr.png" width="600" align="center">

The ZIBR model combines two regression models, logistic regression and beta regression. Each regression has random effects to account for correlations acorss time points. We call these two regressions in ZIBR model as logistic component and beta component. These two components model two different aspects of the data. The logistic component models presence/absence of the microbe and Beta component to model non-zero microbial abundance.

Accordingly, we can test three biologically relevant null hypotheses:  
1. **H0: α_j = 0** the covariates are associated with the bacterial taxon by affecting its presence or absence. This is to test the coefficients in the logistic component;  
2. **H0: β_j = 0**  the taxon is associated with the covariates by showing different abundances. This is to test the coefficients in the beta component;  
3. ** H0: α_j = 0 and β_j = 0 for each covariate X_j and Z_j** the covariates affect the taxon both in terms of presence/absence and its non-zero abundance. This is to joinly test the coefficients in the logistic and beta component. 

## Installation
You can install our ZIBR package from Github
```r
install.packages("devtools")
devtools::install_github("chvlyl/ZIBR")
library(ZIBR)
```

## Basic Usage

```r
zibr(logistic.cov=logistic.cov,beta.cov=beta.cov,Y=Y,subject.ind=subject.ind,time.ind=time.ind)
```

- **logistic.cov**: covariates for the logistic component. Rows: samples. Columns: covariates.  
- **beta.cov**: covariates for the beta component. Rows: samples. Columns: covariates.  
- **Y**: the response variable (i.e the bacterial relative abundance). It is a vector with values in [0,1).  
- **subject.ind**: the variable with subject IDs.   
- **time.ind**: the variable with time points.   

The ordering of the samples in the above matrix or vectors must be consistent. 

The zibr function will return the following results:
- **logistic.est.table**: the estimated coefficients for logistic component.  
- **logistic.s1.est**: the estimated standard deviation for the random effect in the logistic component.  
- **logistic.est.table**: the estimated coefficients for logistic component.  
- **beta.s2.est**: the estimated standard deviation for the random effect in the beta component.  
- **beta.v.est**: the estiamted dispersion parameter in the beta component.  


## Examples
The following function will simulate some data according to the zero-inflated beta random effect model. We specify the covariates in the logistic component (X) and covariates in the beta component to be the same (i.e set Z=X).

```r
sim <- simulate_zero_inflated_beta_random_effect_data(
    subject.n=100,time.n=5,
    X = as.matrix(c(rep(0,50*5),rep(1,50*5))),
    Z = as.matrix(c(rep(0,50*5),rep(1,50*5))),
    alpha = as.matrix(c(-0.5,1)),
    beta = as.matrix(c(-0.5,0.5)),
    s1 = 1,s2 = 0.8,
    v = 5,
    sim.seed=100)
```

The simulation function returns the the bacterial abundance (Y) simulated according to the above parameter settings. This function also returns the other variables such as X, Z, alpha, beta etc. that we just specified. It also returns two variables subject.ind and time.ind, which are subject IDs and time points for each subject.


We can run the zibr function to fit the zero-inflated beta random effect model on the simulated data.
```r
zibr.fit <- zibr(logistic.cov = sim$X, beta.cov = sim$Z, Y = sim$Y,
    subject.ind = sim$subject.ind,time.ind = sim$time.ind)
zibr.fit
```



## Citation
Eric Z. Chen and Hongzhe Li (2016). A two-part mixed effect model for analyzing longitudinal microbiome data. Bioinformatics. [Link](http://bioinformatics.oxfordjournals.org/content/early/2016/05/14/bioinformatics.btw308.short?rss=1)

## Contact
Feel free to contact me by chvlyl AT gmail.com
