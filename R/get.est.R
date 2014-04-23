get.est <-
function(x  ,  data  ,  t1  ,  formula ){

  s1<-Surv(time = as.numeric(x) , event = data$right != Inf )
  #data.frame(s1,data)
  formula1<-formula(paste('s1' , deparse(formula)))
  cph1<-coxph(formula1 , data = data)
  vc1<-vcov(cph1)
  coef1<-cph1$coefficients
  bh<-basehaz(cph1 ,  centered  =  F)
  bh2<-rbind(c(0 , 0) , bh)
  ts<-sapply(t1 , function(x) tail(which(bh2$time<= x) , 1))
  bh3<-bh2[ts , ]
  surv_1<-exp(-bh3$hazard)
  return(list(surv = surv_1 , sigma = vc1 , beta = coef1 , hazard=bh3$hazard))

}
