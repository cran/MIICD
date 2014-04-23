ANDA.coxph <-
function(formula , data , imax = 25 , k = 10 , th0 = 1e-3 )
  {
  #function parameters
  prep<-preproc.coxph( data , k  )
  data2<-prep$data2
  data1<-prep$data1
  t1<-prep$t1
  or<-prep$or
  I<-prep$I 
    
  mm<-model.matrix(formula , data)
  dim_beta<-ncol(mm)-1
  beta<-matrix(0,ncol=dim_beta,nrow=1)
  dn<-dimnames(mm)[[2]][-1]
  sigma<-matrix(rep(1e-3,dim_beta^2),ncol=dim_beta,dimnames=list(dn,dn))  
  
  beta_AN<-mvrnorm(n=k,mu=beta,Sigma=sigma)
  
  #Step1
  # generate sets for the algorithm initialization
  sets<-sapply(1:k , get.set , data)
  
  # derive Breslow estimate of the baseline hazard
  bbh1<-apply(sets  ,  2  ,  get.est  ,  data  ,  t1  ,  formula = ~treatment )
  bbh2<-sapply(bbh1 , function(x) x$hazard)
  bbh<-rowMeans(bbh2)
   
  #initial linear predictors
  Z<-apply( beta_AN , 1 , function(x)  mm[,-1]%*%as.matrix(x))
  surv<-data.frame(time = t1 , surv = exp(-bbh))
  #plot(surv$time,surv$surv,type='l')  
  
  ss1<-apply(data2 , 1 , function(x ) subset( surv , time >=  as.numeric(x['left']) & time <= as.numeric(x['right']) ) )
  tk2<-sapply(seq_len(nrow(data2)) ,function(X)  ss1[[X]]$time)
  
  
  #Step 2  
  # make matrix to get results
  cat('\n\nIterations\n' )
  beta_iter<-matrix(NA,ncol=imax,nrow=dim_beta,dimnames=list(dn,1:imax))
  sigma_iter<-matrix(NA,ncol=imax,nrow=dim_beta,dimnames=list(dn,1:imax))
  sigma1_iter<-matrix(NA,ncol=imax,nrow=dim_beta,dimnames=list(dn,1:imax))
  th_iter<-matrix(NA,ncol=imax,nrow=dim_beta,dimnames=list(dn,1:imax))
  
  #progression bar
  pb <- txtProgressBar(style = 2 , char = '.')
  i<-0
  
  repeat {
    i<-i+1 
    setTxtProgressBar(pb ,  i%%4/150 ) 
    if( i > imax ){setTxtProgressBar(pb ,  0.02 )
    break }
    
    #Get the samples from drown betas from the normal mixture
    samples<-t(sapply( seq_len(nrow(data2)) , function(X) {
    pk2<-sapply( 1:k ,function(x) c(0,diff(1-ss1[[X]]$surv)^exp(Z[I,][X,x]) ) )
    apply( pk2 , 2 ,function(x)
    if( sum(x) ) sample(  tk2[[X]] , size = 1 , prob = x ) 
    else  sample(  tk2[[X]] , size = 1 )  )}))
    
    samples2<-rbind(samples,data1)[or,]
    
    
    #Get estimates
    est_1<-apply(samples2  ,  2  ,  get.est  ,  data  ,  t1  , formula )
    #Get the beta from the successive sets
    betas<-matrix(unlist(sapply(est_1  ,  function(x) x$beta))  ,  nrow  =  dim_beta  ,  dimnames  =  list(dn  ,  1:k))
    #Get the variance 1: beta x t(beta)
    betas2<-array(sapply(1:k,function(x) betas[,x]%*%t(betas[,x])), dim = c(dim_beta ,  dim_beta , k) , dimnames = list(dn , dn ,1:k))
    #Get mean of beta over samples
    beta <- rowMeans(matrix(unlist(sapply(est_1  ,  function(x) x$beta))  ,  nrow  =  dim_beta  ,  dimnames  =  list(dn  ,  1:k)))
    #calculate (1/m sum beta)t(1/m sum beta)
    biv<-beta%*%t(beta)
    #Get the sigma in an d3 array
    sigma1<-array( sapply( est_1 , function(x) x$sigma)  ,  dim  =  c(dim_beta  ,  dim_beta  ,  k) )
    msigma1<-apply(sigma1,c(1,2),mean)
    #add the two terms to get the estimate of the first term of the variance variance (3d array)
    wiv1<-sigma1+betas2
    #average over datasets
    #wiv2<-matrix(rowMeans(matrix(apply(wiv1 , 3 , unlist) , ncol = k)) , ncol = dim_beta , dimnames = list(dn , dn))
    wiv2<-apply(wiv1,c(1,2),mean)
    
    #combine the two term and get an estimate of the variance
    sigmaf<-wiv2-biv
    #sample beta_j from a possibliy multivariate normal with mean means of the betas and sigma obtained before
    beta_AN<-mvrnorm(n=k,mu=beta,Sigma=sigmaf)
    #get new linear predictor
    Z<-apply( beta_AN , 1 , function(x)  mm[,-1]%*%as.matrix(x))
    #get the new estimate of the survival
    surv_<-rowMeans(sapply(est_1 , function(x) x$surv))
    surv <- data.frame(time = t1 , surv = surv_)
    #save estimate of beta and sigma for the prensent iteration
    beta_iter[,i]<-beta
    sigma_iter[,i]<-diag(sigmaf)
    sigma1_iter[,i]<-diag(msigma1) 
    #convergence criteria
    th_iter[,i]<-rowMeans(beta_iter,na.rm=T)
  }
  th_iter<-matrix(th_iter[!is.na(th_iter)],nrow=dim_beta)
  close(pb)
  #compute the mean of beta and sigma for at most the 25 last iterations without the 5 first results 
  
  ret<-list(beta=beta_iter,sr_sigma=sigma_iter,n_iter=i,conv=th_iter,sigma1=sigma1_iter)
  return(ret)
}
