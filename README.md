# PLNseq

PLNseq: A multivariate Poisson lognormal distribution for high-throughput correlated RNA-sequencing read count data

This R package conducts differential expression (DE) analysis using high throughput next-generation sequencing read count data generated from correlated samples. The marginal distribution of the read count is the compounding of the Poisson distribution and the lognormal distribution (`PLN' distribution for short), and the correlation between the read counts of each matched sample set is modeled by the multivariate lognormal distribution. This package provides estimates of rho (correlation coefficient matrix in multivariate lognormal distribution) and its standard error, sigma (standard deviation of lognormal distribution), log-fold change (defined as the difference between log-gene expression of matched samples), and p-value for detecting differentially expressed genes.

There are three main functions: LRtest1 (with a common correlation shared by all genes), LRtest2 (with gene cluster specific correlations), and PLN_ANOVA (with a rank-reduced ANOVA model).

The manual file is "PLNseq-manual.pdf". 

Installation of PLNseq in R:

> library(‘devtools’);

> install_github(‘zhanghfd/PLNseq’);
