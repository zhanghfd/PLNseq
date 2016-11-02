
correlationCoefficient <-
function(d){

    ####################
    # estimate rho
    ####################
    R = d$conditionNumber;
    I = ncol(d$count)/R;
    J = nrow(d$count);
    count = sweep(d$count,2,d$sample$sizeFactor,'/');
    
    dat = array(NA,c(J,R,I));
    for(r in 1:R){
        dat[,r,] = count[,1:I+(r-1)*I];
    }

    sig = mean(d$commonSigma);

    phi = (exp(sig^2)-1);

    rhos = array(NA,c(J,R,R));

    for(j in 1:J){
        x = t(dat[j,,]);
        v = apply(x,2,var);
        if(any(v==0)){
            rhos[j,,] = diag(R);
        }else{
            rhos[j,,] = cor(x);
        }
    }
    
    n = I;
    tmp = gamma((n-1)/2)/gamma(1/2)/gamma((n-2)/2);
    
    cor.new <- 
    function(r){
        n.r = length(r);
        res = rep(1/tmp,n.r);
        
        for(i in 1:n.r){
          r0 = r[i];
          if(!r0 %in% c(0,1,-1)){
            f.integral <- 
              function(t){
                return(t^{-1/2}*(1-t)^{(n-2)/2-1}/(1-t*(1-r0^2))^{1/2})
              }
            res[i] = integrate(f.integral,0,1)$value;
          }
        }
        res = tmp * res * r;
        return(res);
    }
    
    for(r1 in 1:(R-1)){
      for(r2 in r1:R){
        rhos[,r1,r2] = cor.new(rhos[,r1,r2]);
        rhos[,r2,r1] = rhos[,r1,r2];
      }
    }
#    rhos = rhos * (I-1)/(I-1.75);

    rhohat = matrix(1,R,R);
    s.rho = matrix(0,R,R);
    inv.alpha = 1/apply(dat,1:2,mean);
    for(r1 in 1:(R-1)){
        for(r2 in (r1+1):R){
            tmp = sqrt((phi+inv.alpha[,r1]) * (phi+inv.alpha[,r2])) * rhos[,r1,r2];
            rhohat[r1,r2] = log(1+mean(tmp))/(sig^2);
            s.rho[r1,r2] = sd(tmp)/sqrt(J)/sig^2/(1+mean(tmp));
        }
    }
    for(r1 in 2:R){
        for(r2 in 1:(r1-1)){
            rhohat[r1,r2] = rhohat[r2,r1];
            s.rho[r1,r2] = s.rho[r2,r1];
        }
    }

    d$rho = rhohat;
    d$rho.se = s.rho;

    return(d);

}
