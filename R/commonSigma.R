commonSigma <-
function(d){
    
    dat0 = sweep(d$count,2,d$sample$sizeFactor,'/');

    R = d$conditionNumber;
    I = ncol(dat0)/R;
    J = nrow(dat0);
    
    x = array(NA,c(J,R,I));
    for(r in 1:R){
      x[,r,] = dat0[,1:I+(r-1)*I];
    }
    
    alpha = apply(x,1:2,mean);
    
    ##########################################
    # deviance-based estimates for sigma
    ##########################################
    fn = function(sig){
        phi.inv = 1/(exp(sig^2)-1);
        res = 2*sum(x*log(sweep(x,1:2,alpha,'/'))-(x+phi.inv)*log(sweep(x+phi.inv,1:2,alpha+phi.inv,'/'))) - R*(I-1)*J;
        return(res);
    }

    res = uniroot(fn, c(0.1,0.8), tol=1e-6);
    d$commonSigma = res$root;

    return(d);

}
