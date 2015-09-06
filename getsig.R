getsig=function(allcors.df,p.crit,var.x,var.x.class=FALSE) {
  sig.indxs=allcors.df$p<p.crit
  # a few tests don't make sense and yield NaN (have verified these)
  sig.indxs[is.na(sig.indxs)]=FALSE
  allcors.df.sig=allcors.df[sig.indxs,]
  allcors.df.sig=allcors.df.sig[order(allcors.df.sig$p),]
  allcors.df.sig=allcors.df.sig[order(allcors.df.sig$r2,decreasing=T),]
  
  # find all significant correlations with a specific variable
  var.x=var.x
  sig.x=(var.x==allcors.df.sig[,1])|(var.x==allcors.df.sig[,2])
  
  # find all significant correlations with a class of variables
  if (var.x.class==T) {
    sig.x=grepl(var.x,allcors.df.sig[,1])|grepl(var.x,allcors.df.sig[,2])
  }
  
  show.sig=allcors.df.sig[sig.x,]
  show.sig[,1:2]=as.matrix(show.sig[,1:2])
  show.sig[!grepl(var.x,show.sig$var1),c(1,2)]=as.matrix(show.sig[!grepl(var.x,show.sig$var1),c(2,1)])
  show.sig=show.sig[order(show.sig$var2),]
  return(show.sig)
}
