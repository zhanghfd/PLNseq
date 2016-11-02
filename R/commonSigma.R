commonSigma <-
function(d){
    
    dat0 = sweep(d$count,2,d$sample$sizeFactor,'/');

    I = ncol(dat0);
    J = nrow(dat0);
    alpha = matrix(rep(rowMeans(dat0),I),ncol=I);
    
    ##########################################
    # deviance-based estimates for sigma
    ##########################################
    fn = function(sig){
        phi.inv = 1/(exp(sig^2)-1);
        res = 2 * sum(dat0*log(dat0/alpha)-(dat0+phi.inv)*log((phi.inv+dat0)/(phi.inv+alpha)))-(I-1)*J;
        return(res);
    }

    res = uniroot(fn, c(0.1,0.8), tol=1e-6);
    d$commonSigma = res$root;

    return(d);

}
