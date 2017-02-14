
PLN_ANOVA = function(d,n.top=1e3){

  R = d$conditionNumber;  
  count = d$count;
  ct = sweep(count,2,d$sample$sizeFactor,'/');
  
  J = nrow(count);
  I = ncol(count)/R;
  y = array(NA,c(J,R,I));
  for(r in 1:R){
    y[,r,] = ct[,1:I+(r-1)*I]
  }
  
  condition = as.factor(rep(1:R,each=I));
  p.value = rep(NA,J);
  for(j in 1:J){
    x = ct[j,];
    fit = lm(x~condition);
    p.value[j] = anova(fit)[1,5];
  }
  top.id = order(p.value)[1:n.top];
  
  mu = apply(y,1:2,mean);
  etahat = log(mu);
  muhat = mean(etahat);
  
  alphahat = colMeans(etahat) - muhat;
  betahat = rowMeans(etahat) - muhat;
  z = etahat - muhat - outer(betahat,alphahat,'+');
  
  u1 = c(1,rep(-1/(R-1),R-1));
  u0 = u1 + 1;
  ep = 1e-6;
  z1 = z[top.id,];
  
  if(R==2){
    uhat = c(1,-1);
    tmp = as.numeric(z%*%(uhat));
    vhat = (tmp)/sum(uhat^2);
  }else{
    while(max(abs(u1-u0))>ep){
      u0 = u1;
      v1 = as.numeric(z1%*%(u0))/sum(u0^2);
      u1 = as.numeric(v1%*%z1)/sum(v1^2);
      u1 = u1/u1[1];
    }
    uhat = u1;
    vhat = as.numeric(z%*%(uhat))/sum(uhat^2);
  }
  
  
  sigma = d$genewiseSigma;  
  rho = d$rho;
  
  var.v = rep(0,J);
  for(r1 in 1:R){
    for(r2 in 1:R){
      if(r1==r2){
        var.v = var.v + uhat[r1]^2*(exp(sigma^2)-1+1/mu[,r1])/I;
      }else{
        var.v = var.v + uhat[r1]*uhat[r2]*(exp(rho[r1,r2]*sigma^2)-1)/I;
      }
    }
  }
  var.v = var.v/sum(uhat^2)^2;
  
  sd.v = sqrt(var.v);
  stat = (vhat/sd.v)^2;
  
  th = qchisq(1-1e-5,1);
  p.new = ifelse(stat<th,1-pchisq(stat,1),-pchisq(stat,1,log=TRUE));
  log2fold.change = outer(as.numeric(vhat),as.numeric(uhat-uhat[1]))/log(2);
  sd.log2foldchange = outer(as.numeric(sd.v),abs(as.numeric(uhat-uhat[1])))/log(2);
  d$ANOVA = list(u=uhat,log2FoldChange=log2fold.change,sd.log2FoldChange=sd.log2foldchange,p.value=p.new));
  
  return(d);
}
