sizeFactor <-
  function(d,maxCount=NA){
    count = d$count;
    ms = as.numeric(apply(as.matrix(count),2,'median'));
    m0 = exp(mean(log(ms)));
    
    d$sample$sizeFactor = ms/m0;
    names(d$sample)[4] = 'sizeFactor';
        
    row.max = apply(count,1,max);
    
    if(!is.na(maxCount)){
      count[row.max>maxCount,] = round(sweep(count[row.max>maxCount,],1,row.max[row.max>maxCount],'/')*maxCount);
    }
    
    count[count==0] = 1;
    
    d$count = count;
    
    rm(count,ms,m0);
    
    return(d);
  }
