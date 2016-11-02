PLNobject <-
  function(count,conditionNumber=2){
    
    nSample = ncol(count);
    
    if(nSample%%conditionNumber!=0){
      stop(paste('Number of conditions should be',conditionNumber,sep=' '));
    }else{
      I = nSample/conditionNumber;
    }
    
    if(is.null(names(count))){
      samplename = NULL;
      for(k in 1:conditionNumber){
        samplename = c(samplename,paste('C',k,'.Sample',1:I,sep=''));
      }
    }else{
      samplename = names(count);
    }
    
    sample = data.frame(samplename,as.numeric(colSums(count)),as.numeric(apply(count,2,median)));
    names(sample) = c('sampleName','totalCount','medianCount');
    
    d = list(count=as.matrix(count),conditionNumber=conditionNumber,sample=sample);

    rm(samplename,sample);
    return(d);
  }
