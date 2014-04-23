get.est.cr.cox <-
function( x , data , t1 , status , trans , cens , keep , formula ){

  
  df1 <- data.frame( data , 'ns' = x + 1e-6 )
  df2 <- crprep( data = df1 , Tstop = 'ns' , status = status , trans = trans , cens = cens , keep = keep )
  s1 <- with( df2 , Surv( time = Tstart , time2 = Tstop , type = 'counting' , event = status==trans ) )
  head(df2)
  formula1 <- formula( paste( 's1' , deparse( formula ) ) )
  cph1 <- coxph( formula1 , data = df2 , weights = df2$weight.cens )

vc1 <- vcov( cph1 )
coef1 <- cph1$coefficients
bh <- basehaz( cph1 , centered = F )
 #plot(bh$time,exp(-bh$hazard))
  bh2 <- rbind( c( 0 , 0 ) , bh )
ts <- sapply( t1 , function( x ) tail( which( bh2$time <= x ) , 1 ) )
bh3 <- bh2[ ts , ]
bh3 <- bh2[ ts , ]
surv_1 <- exp( -bh3$hazard )
return( list( surv = surv_1 , sigma = vc1 , beta = coef1) )
}
