allcors=function(data,is.DV=NULL) {
  source("/Users/eric8226/Box Sync/Eric's R functions/anycor.R")
  allcors.mat=matrix(NA,ncol(data),ncol(data))
  rownames(allcors.mat)=colnames(data)
  colnames(allcors.mat)=colnames(data)
  all.tests=list()
  allcors.df=matrix(,0,6)
  colnames(allcors.df)=c("var1","var2","r2","p.val","test.type","beta")
  # force certain variables to be the DV
  if (is.null(is.DV)) {is.DV=rep(FALSE,ncol(data))}
  for (i in 1:(nrow(allcors.mat)-1)) {
    for (j in (1+i):ncol(allcors.mat)) {
      i<<-i
      j<<-j
      x=rownames(allcors.mat)[i]
      y=colnames(allcors.mat)[j]
      is.DV_ij=paste(c("","x")[1+is.DV[i]*1],
                     c("","y")[1+is.DV[j]*1],sep="")
      out=anycor(x,y,is.DV_ij,data)
      # save effect size, p-value, and test output
      allcors.df=rbind(allcors.df,c(x,y,out$r2,out$p.val,out$test.type,out$beta.val))
      all.tests[[paste(x,y,sep="_X_")]]=out$test.out
    }
  }
  allcors.df=as.data.frame(allcors.df)
  for (i in c(3,4,6)) {
    allcors.df[,i]=as.numeric(as.character(allcors.df[,i]))
  }
  return(list(allcors.df,all.tests))
}