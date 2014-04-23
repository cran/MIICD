PMDA.coxph <-
function(formula , data , imax = 25 , k = 10 , th0 = 1e-3 ){
  #function parameters
  
  prep<-preproc.coxph( data , k )
  data2<-prep$data2
  data1<-prep$data1
  t1<-prep$t1
  or<-prep$or
  I<-prep$I 
  
  mm<-model.matrix(formula , data)
  dim_beta<-ncol(mm)-1
  beta<-matrix(0,ncol=dim_beta,nrow=1)
  dn<-dimnames(mm)[[2]][-1]
  
  #Step1
  sets<-sapply(1:k , get.set , data)
  bbh1<-apply(sets  ,  2  ,  get.est  ,  data  ,  t1  ,  formula = ~treatment )
  bbh2<-sapply(bbh1 , function(x) x$hazard)
  bbh<-rowMeans(bbh2)
   
  #initial linear predictors
  Z<-mm[,-1]%*%t(beta)
  surv<-data.frame(time = t1 , surv = exp(-bbh))
  #plot(surv$time,surv$surv,type='s')
  
  ss1<-apply(data2 , 1 , function(x ) subset( surv , time >=  as.numeric(x['left']) & time <= as.numeric(x['right']) ) )
  tk2<-sapply(seq_len(nrow(data2)) ,function(X)  ss1[[X]]$time)
  
  #Step 2  
  cat('\nIterations\n' )
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

    Z<-matrix(rep(Z,k),ncol=k,byrow=F)
    #Get the samples from drown betas from the normal mixture
    samples<-t(sapply( seq_len(nrow(data2)) , function(X) {
    pk2<-sapply( 1:k ,function(x) c(0,diff(1-ss1[[X]]$surv)^exp(Z[I,][X,x]) ) )
    apply( pk2 , 2 ,function(x)
    if( sum(x) ) sample(  tk2[[X]] , size = 1 , prob = x ) 
    else  sample(  tk2[[X]] , size = 1 )  )}))
    
    samples2<-rbind(samples,data1)[or,]
    
    est_1<-apply(samples2  ,  2  ,  get.est  ,  data  ,  t1  , formula )
        
    betas<-matrix(unlist(sapply(est_1  ,  function(x) x$beta))  ,  nrow  =  dim_beta  ,  dimnames  =  list(dn  ,  1:k))
    beta <- rowMeans(matrix(unlist(sapply(est_1  ,  function(x) x$beta))  ,  nrow  =  dim_beta  ,  dimnames  =  list(dn  ,  1:k)))
    sigma<-array(sapply(est_1 , function(x) x$sigma)  ,  dim  =  c(dim_beta  ,  dim_beta  ,  k))
    msigma1<-apply(sigma,c(1,2),mean)
    
    surv_<-rowMeans(sapply(est_1 , function(x) x$surv))
    
    delta_sigma<-matrix( (1 + ( 1/k ))*rowSums( ( betas-beta )**2 / ( k-1 ) ),ncol=dim_beta ,dimnames=list('',dn))
    Z <- mm[,-1]%*%as.matrix(beta)
    
    sigma<-diag(as.matrix(msigma1))+diag(as.matrix(delta_sigma))
    surv <- data.frame(time = t1 , surv = surv_)
    
    beta_iter[,i]<-beta
    sigma_iter[,i]<-sigma
    sigma1_iter[,i]<-diag(msigma1)
    #convergence criteria
    th_iter[,i]<-rowMeans(beta_iter,na.rm=T)
  }

  th_iter<-matrix(th_iter[!is.na(th_iter)],nrow=dim_beta)
  close(pb)
  
    ret<-list(beta=beta_iter,sr_sigma=sigma_iter,n_iter=i,conv=th_iter, sigma1 = sigma1_iter)
  return(ret)
}
