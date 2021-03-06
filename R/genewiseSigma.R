genewiseSigma <-
  function(d,w=25){
    
    dat = sweep(d$count,2,d$sample$sizeFactor,'/');
    
    J = nrow(dat);
    R = d$conditionNumber;
    I = ncol(dat)/R;
    
    x = array(NA,c(J,R,I));
    for(r in 1:R){
      x[,r,] = dat[,1:I+(r-1)*I];
    }
    
    alpha = apply(x,1:2,mean);
    
    w1 = w/((I-1)*R);
    
    top.n = 3e2;
    sigmas = seq(0.01,1,.01);
    n.sigma = length(sigmas);

    fn = function(sig){
      phi.inv = 1/(exp(sig^2)-1);
      res = 0;
      for(r in 1:R){
        res = res + 2*apply(x[,r,]*log(sweep(x[,r,],1,alpha[,r],'/'))-(x[,r,]+phi.inv)*log(sweep(x[,r,]+phi.inv,1,alpha[,r]+phi.inv,'/')),1,sum) - (I-1);
      }
      return(res);
    }
    
    scores = matrix(NA, J, n.sigma);
    for(k in 1:n.sigma){
      scores[,k] = fn(sigmas[k]);
    }
    
    mN = rowMeans(dat);
    
    rk = round(rank(mN));
    
    index = order(mN);
    
    i1 = index[1:(top.n-1)];
    score1 = colMeans(scores[i1,]);
    ranked.mscore = matrix(NA,J,n.sigma);
    
    for(j in 1:(top.n/2)){
      ranked.mscore[j,] = score1;
    }
    
    for(j in (top.n/2+1):(J-top.n/2)){
      tmp = index[(j-top.n/2+1):(j+top.n/2-1)];
      ranked.mscore[j,] = colMeans(scores[tmp,]);
    }
    
    i2 = index[(J-top.n+2):J];
    score2 = colMeans(scores[i2,]);
    for(j in (J-top.n/2+1):J){
      ranked.mscore[j,] = score2;
    }
    
    
    sigma = rep(NA,J);
    for(j in 1:J){
      sigma[j] = sigmas[which.min(abs(w1*ranked.mscore[rk[j],]+scores[j,]))];
    }
    
    d$genewiseSigma = sigma;
    
    return(d);
  }
