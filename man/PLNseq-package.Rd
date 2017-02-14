\name{PLNseq-package}
\alias{PLNseq-package}
\docType{package}
\title{
Differential expression analysis using matched read count data
}
\description{

This R package conducts differential expression (DE) analysis using high throughput next-generation sequencing read count data generated from correlated samples. The marginal distribution of the read count is the compounding of the Poisson distribution and the lognormal distribution (`PLN' distribution for short), and the correlation between the read counts of each matched sample set is modeled by the multivariate lognormal distribution with correlation coefficient matrix that is assumed to be common for all genes. This package provides estimates of rho (correlation coefficient matrix in multivariate lognormal distribution) and its standard error, sigma (standard deviation of lognormal distribution), log2-fold change (defined as the difference between log2-gene expression of matched samples) and p-value for detecting differentially expressed genes.

There are three main functions: LRtest1 (with a common correlation shared by all genes), LRtest2 (with gene cluster specific correlations), and PLN_ANOVA (with a rank-reduced ANOVA model).

}
\details{
\tabular{ll}{
Package: \tab PLNseq\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-06-18\cr
License: \tab GPL (>= 3) \cr
}
}
\author{
Hong Zhang

Maintainer: Hong Zhang <zhanghd@fudan.edu.cn>
}

\references{
Zhang, H., Xu, J., Jiang N., Hu, X., and Luo, Z. (2015). PLNseq: A multivariate Poisson lognormal distribution for high-throughput matched RNA-sequencing read count data. Statistics in Medicine 34: 1577-1589.
}
\keyword{ package }

%\seealso{
%}
%\examples{
%}
