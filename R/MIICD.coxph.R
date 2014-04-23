MIICD.coxph <-
function(formula , data , imax = 25 , k = 10 , th0 = 1e-3 ,  method = c('PMDA','ANDA') )
{
  cl <- match.call()
  cat('\nMultiple Imputation for Interval Censored Data for Cox proportional hazards regression\n')
  cat("\nCall:\n", paste(deparse(cl), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if(imax<10){
    imax<-10
    cat('\nWarning in' , deparse(cl) , ' :\nimax < 10, set to 10\n' )
    }
  cat('Interval-censored Response of a Proportional Hazard model\n')
  x<-nrow(data)
  cat('No.Observation:', x , '\n')
  cat('Patern:\n')
  stat<-ifelse(data$right == Inf ,'right-censored', 'event of interest')
  type<-ifelse(data$right == data$left , 'exact' , NA )
  type<-ifelse(data$right!=data$left & data$right!=Inf , 'interval-censored' , type )
  type<-ifelse(data$right == Inf  , 'right-censored' , type )
  print(table('Cause'=stat, type))
  
  method<-match.arg(method)
  if(method=='PMDA')
  res<-PMDA.coxph( formula = formula , data = data , k = k , imax = imax , th0 = th0 )
  else
  if (method == 'ANDA')
  res<-ANDA.coxph(formula = formula , data = data , imax = imax , k = k , th0 = th0 )  
 
  nmean<-25
  mim<-min(res$n_iter-5,nmean)
  mean_beta  <- rowMeans(tail(res$beta , mim),na.rm=T)
  mean_sigma <- rowMeans(tail(res$sr_sigma , mim),na.rm=T)
  #Compute the Pvalues
  pv<-1-pchisq((mean_beta/(mean_sigma^.5))**2 , df = 1)
  #print the results to the terminal
  df1<-data.frame('coef'=mean_beta,'exp(coef)'=exp(mean_beta),'se(coef)'=mean_sigma**.5,'z'= mean_beta/mean_sigma^.5,pv=pv,'.'= '')
  colnames(df1)<-c('coef','exp(coef)' , 'se(coef)' , 'z' , 'p', '')
  cat("\nCoefficients:\n")
  print(format(df1  ,  digits  =  max(3L ,  getOption("digits") - 3L))   ,  print.gap  =  3L  ,  quote  =  FALSE)
  cat('\n')
  cat(paste('Estimates computed using the last ',mim,' iterations over ',res$n_iter-1,'\n\n',sep=''))
  #return the results
  ret<-list('Mean beta' =  mean_beta ,  'Mean sigma'  =  mean_sigma**.5  ,  'betas' = res$beta , 'call' = cl , df = df1 , niter = res$n_iter , conv = res$conv , 'sigma1' = res$sigma1 )
  class(ret)<-'MIICD'
  return(ret)
}
