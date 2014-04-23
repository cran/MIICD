preproc.crreg <-
function(  data = data  , k = k  , trans = trans , status = status , lt=3*nrow(data)){

  rownames(data)<-seq_len(nrow(data))
  I <- data[,status]==trans &  data[,'left'] != data[,'right' ] & data[,'right'] != Inf
  data2<-data[I,]
  t1<-sort(unique(c(data2$left , data2$right)))
  t1<-seq(from=0 , to = range(t1)[2], length=lt)
  dataE<-data[!rownames(data)%in%rownames(data2),]
  or<-order(c(as.numeric(rownames(data2)),as.numeric(rownames(dataE))))
  data1<-t(apply( dataE , 1 , function(x) as.numeric(rep(x[1],k))))

return(list(data2=data2,t1=t1,data1=data1,or=or,I=I))
}
