LRtest2 <-
function(d,M,use.commonSigma=FALSE,id=NULL){

    count = d$count;
    R = d$conditionNumber;
    I = ncol(count)/R;
    J = nrow(count);
    
    dat = array(NA,c(J,R,I));
    for(r in 1:R){
      dat[,r,] = count[,1:I+(r-1)*I];
    }

    log.size = matrix(log(d$sample$sizeFactor),ncol=R); # I x R matrix
    
    est.mu.lrt = function(I,mu0,sigma0,z,x){
      
      sig.z = sigma0*z;
      
      fn1 = function(mu){
        loglik = matrix(0,M,I);
        for(i in 1:I){
          for(r in 1:R){
            log.ld = sig.z[,r] + mu[r] + log.size[i,r];
            loglik[,i] = loglik[,i] + x[r,i] * log.ld - exp(log.ld);
          }
        }
        loglik.max = apply(loglik,2,max);
        l = sum(loglik.max+log(colMeans(exp(sweep(loglik,2,loglik.max)))));
        return(-l);
      }
      
      res1 = optim(mu0,fn1,method='Nelder-Mead',control=list(maxit=1e4));
      
      muhat = res1$par;
      
      fn0 = function(mu){
        loglik = matrix(0,M,I);
        for(i in 1:I){
          for(r in 1:R){
            log.ld = sig.z[,r] + mu + log.size[i,r];
            loglik[,i] = loglik[,i] + x[r,i] * log.ld - exp(log.ld);
          }
        }
        loglik.max = apply(loglik,2,max);
        l = sum(loglik.max+log(colMeans(exp(sweep(loglik,2,loglik.max)))));
        return(-l);
      }
      res0 = optimize(fn0,lower=min(muhat),upper=max(muhat));
      
      lrt = 2*(res0$objective-res1$value);
      
      return(list(mu1=res1$par,mu0=res0$min,lrt=lrt));
      
    }
  

    mu.mar = d$mu;

    alpha = apply(dat,1:2,mean);
    mu.mar = log(alpha)-d$commonSigma^2/2;

    if(use.commonSigma){
        if(is.null(d$commonSigma)){
          stop("Run 'commonSigma' first.");
        }else{
          sigma.mar = rep(d$commonSigma,J);
        }
    }else{
      if(is.null(d$genewiseSigma)){
        stop("Run 'genewiseSigma' first.");
      }else{
        sigma.mar = d$genewiseSigma;
      }
    }

    #################################
    # likelihood ratio test statistic
    #################################

    if(is.null(id)){
        id = 1:J;
    }
    
    id = intersect(id,1:J);

    J0 = length(id);
    
    LR0 = matrix(NA,J0,R+2);
    lrt0 = rep(NA,J0);
    mu0 = matrix(NA,J0,R);

    j1 = 0;
    for(j in id){
        rho = matrix(d$rho[,d$cluster[j]],R);
        z = mvrnorm(n = M, mu=rep(0,R), Sigma=rho);
        x = dat[j,,];
        res = est.mu.lrt(I,mu.mar[j,],sigma.mar[j],z,x);
        j1 = j1 + 1;
        mu0[j1,] = res$mu1;
        lrt0[j1] = res$lrt;
    }

    LR0[,1:R] = sweep(mu0,1,mu0[,1]);
    LR0[,1+R] = lrt0;

    LR0[,2+R] = ifelse(lrt0 < 20, 1-pchisq(lrt0,R-1), -pchisq(lrt0,R-1,log.p=TRUE));

    if(is.null(d$LR)){
      d$LR = data.frame(matrix(NA,J,R+2));
      names(d$LR) = c(paste('logFC',1:R,sep=''),'LR statistic','p.value');
      d$LR[id,] = LR0;
    }else{
      d$LR[id,] = LR0;
    }
    
    return(d);

}
