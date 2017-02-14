\name{PLN-ANOVA}
\alias{PLN-ANOVA}

\title{
Differential expression analysis based on a rank-reduced ANOVA model.
}
\description{
This function returns the u estimates and DE analysis results including estimated log-fold changes and the correponding estimated standard deviations, and DE test p-values.
}
\usage{
    PLN_ANOVA(d,n.top=1e3)
}

\arguments{
    \item{d}{
    This is a PLNseq object.
    }
    \item{n.top}{
    The number of genes used to estimated u and v parameters.
    }
}

\value{
    \item{ANOVA}{A list containg `u', `log2FoldChange', and `sd.log2FoldChange', and `p.value'}
}


\examples{

    data(lung);
    count = lung[,c(2:4,8:10)];
    d = PLNobject(count,conditionNumber=2);
    d = sizeFactor(d,maxCount=2e3);

    ##  Not run:
    ##  d = commonSigma(d);  
    ##  d = genewiseSigma(d);
    ##  d = correlationCoefficient(d);
    ##  d = PLN_ANOVA(d,n.top=1e3);
}

\keyword{differential expression}
\keyword{rank-reduced ANOVA model}

