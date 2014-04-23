MIICD.crreg <-
function( formula , data , imax = 25 , k = 10 , th0 = 1e-3 , status , trans , cens , keep , method = c('PMDA','ANDA') , model = c('FG','Cox'))
{
  
  cl <- match.call()
  cat('\nMultiple Imputation for Interval Censored Data for proportional hazards regression with competing risks\n')
  cat("\nCall:\n", paste(deparse(cl), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if(imax<10){
    imax<-10
    cat('\nWarning in' , deparse(cl) , ' :\nimax < 10, set to 10\n' )
    }
  cat('Interval-censored response of a competing risks model\n')
  x<-nrow(data)
  cat('No.Observation:', x , '\n')
  cat('Patern:\n')
  stat<-ifelse(data[,status]==cens,'unknown (right-censored)',as.character(data[,status]))
  type<-ifelse(data$right==data$left , 'exact' , NA )
  type<-ifelse(data$right!=data$left & data$right!=Inf , 'interval-censored' , type )
  type<-ifelse(data[,status]==cens  , 'right-censored' , type )
  print(table('Cause'=stat, type))
  
  method <- match.arg(method)
  model <- match.arg(model)
  
  if( method == 'PMDA' )
  res<-PMDA.crreg( formula , data , imax  , k , th0 ,  status , trans , cens , keep , model )
  else
  if( method == 'ANDA' )
  res<-ANDA.crreg( formula = formula , data , imax , k , th0 ,  status , trans , cens  , keep , model )

  nmean<-25
  mim<-min(res$n_iter-5,nmean)
  mean_beta  <- rowMeans(tail(res$beta , mim),na.rm=T)
  mean_sigma <- rowMeans(tail(res$sr_sigma , mim),na.rm=T)
  mean_sigma1 <- rowMeans(tail(res$sigma1 , mim),na.rm=T)
  #Compute the Pvalues
  pv<-1-pchisq((mean_beta/(mean_sigma^.5))**2 , df = 1)
  #print the results to the terminal
  

  df1<-data.frame('coef'=mean_beta,'exp(coef)'=exp(mean_beta),'se(coef)'=mean_sigma**.5,'z'= mean_beta/mean_sigma^.5,pv=pv,'.'= '')
  colnames(df1)<-c('coef','exp(coef)' , 'se(coef)' , 'z' , 'p', '')
  cat("Coefficients:\n")
  print(format(df1  ,  digits  =  max(3L ,  getOption("digits") - 3L))   ,  print.gap  =  3L  ,  quote  =  FALSE)
  cat('\n')
  cat(paste('Estimates computed using the last ',mim,' iterations over ',res$n_iter-1,'\n\n',sep=''))
    #return the results
  
  
  ret<-list('Mean beta' =  mean_beta ,  'Mean sigma'  =  mean_sigma  ,  'betas' = res$beta , call = cl , df = df1 , niter = res$n_iter , conv = res$conv , 'sigma1' = mean_sigma1)
  class(ret)<-'MIICD'
  return(ret)
}
