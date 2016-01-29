# ZIBRE (Zero-Inflated Beta Random Effect model)

## Introduction
The longitudinal microbiome compositional data are highly skewed, bounded in [0,1), and often sparse with many zeros. In addition, the observations from repeated measures are correlated. We propose a two-part zero-inflated Beta regression model with random effects (ZIBRE) for testing the association between microbial abundance and clinical covariates for longitudinal microbiome data. The model includes a logistic component to model presence/absence of the microbe in samples and a Beta component to model non-zero microbial abundance and each component includes a random effect to take into account the correlation among repeated measurements on the same subject.

## Installation
You can install our ZIBRE package from Github
```r
install.packages("devtools")
devtools::install_github("chvlyl/ZIBRE")
library(ZIBRE)
```

## Basic Usage

```r
zibre(logistic.cov=logistic.cov,beta.cov=beta.cov,Y=Y,subject.ind=subject.ind,time.ind=time.ind)
```

## Citation
Eric Z. Chen and Hongzhe Li (2016). A two-part mixed effect model for analyzing longitudinal microbiome data
