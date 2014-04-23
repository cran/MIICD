get.est.FG <-
function( x , data , t1 , status , trans , cens , keep ,  formula ){
  
  mm<-model.matrix(formula , data)
  m1<-colMeans(mm)[-1]
  df1 <- data.frame( data , 'ns' = x )
  FGR1<-paste("FGR( Hist( time = ns , event = ",status,", cens.code = '", cens ,"' )", deparse(formula) ,", data = df1 , cause = '",trans,"' )" ,sep='' )
  FGR1
  FGR2<-eval(parse(text=FGR1)) 
  vc1<-vcov(FGR2$crrFit)
  coef1<-FGR2$crrFit$coef
  ci<-cuminc(ftime=df1$ns, fstatus=df1[,status],   rho=0, cencode = cens )
  p1<-paste(1,trans)
  H<- -log(1-ci[[p1]]$est)
    bz0<-sum(m1*coef1)
    H <- H * exp( -bz0 )
  
  ci2<-data.frame(time=ci[[p1]]$time,hazard=H)
  
  bh2 <- rbind( c( 0 , 0 ) , ci2 )
  ts <- sapply( t1 , function( x ) tail( which( bh2$time <= x ) , 1 ) )
  bh3 <- bh2[ ts , ]
  surv_1 <- exp( -bh3$hazard )
return( list( surv = surv_1 , sigma = vc1 , beta = coef1 , hazard=bh3$hazard) )
}
