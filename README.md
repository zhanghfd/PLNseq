# PLNseq

PLNseq: A multivariate Poisson lognormal distribution for high-throughput correlated RNA-sequencing read count data
This R package conducts differential expression (DE) analysis using high throughput  RNA-seq read count data generated from correlated samples. The marginal distribution of the read count is the compounding of the Poisson distribution and the log-normal distribution (‘PLN’ distribution for short), and the correlation between the read counts of the correlated samples is modeled by multivariate log-normal distribution with correlation coefficient parameters, which are assumed to be common for all genes. This R package provides an estimate of the correlation coefficient parameters (standard errors), fold change (defined as the difference between log-gene expression) and p-value for testing differentially expressed genes.

The manual file is "PLNseq-manual.pdf". 

Installation of PLNseq in R:

library(‘devtools’);

install_github(‘zhanghfd/PLNseq’);
