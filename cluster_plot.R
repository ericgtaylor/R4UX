cluster.plot=function(cvar,DV,allvars.sub,plot.folder=paste(getwd(),"/",sep="")) {
  title=paste(DV," by ",unlist(strsplit(cvar,"\\."))[1],sep="")
  title=gsub("\\.","_",title)
  res=200
  png(paste(plot.folder,title,".png"),
      res*(5+log(max(1,ncol(allvars.sub)-4))),res*4,res=res)
  par(xpd=T)
  # add margin space depending on number of rows
  add.bmar=ceiling(nrow(x.sum.graph)/4)-2
  par(mar=c(6+add.bmar,4,4,2))
  xs=with(allvars,
          bargraph.CI(allvars[,cvar],allvars[,DV],xaxt="n",yaxt="n",
                      ylim=c(0,4),ylab=DV,cex.lab=.8,
                      main=title))
  axis(2,0:4,c("$0","$10","$100","$1k","$10k"),cex.axis=.8)
  # to add text labels for above-threshold features
  # text(xs$xvals,-.1,gsub("_","\n",names(xs$vals)),srt=90,cex=.4,adj=c(1,.5))
  # add clustering squares
  cnames=names(xs$vals)
  vnames=rownames(x.sum.graph)
  for (i in 1:length(cnames)) {
    cs=4
    rs=ceiling(nrow(x.sum.graph)/4)
    xl=xs$xvals[i]-.5
    bxs.per.row=c(rep(4,floor(length(vnames)/4)),length(vnames)%%4)
    var.means=apply(allvars.sub[allvars[,cvar]==cnames[i],],2,mean,na.rm=T)
    # reorder to match graph (top to bottom in clustering graph goes left to right in this graph)
    var.means=var.means[vnames]
    rect.cols.range=seq(1,.25,length.out=10)
    is.binary=max(allvars.sub,na.rm=T)==1
    rect.cols=cut(var.means,seq(c(1,-.0001)[1+is.binary*1], # use -1 to include 0s
                                max(allvars.sub,na.rm=T),length.out=11),
                  rect.cols.range)
    rect.cols=as.numeric(as.character(rect.cols))
    # record the x/y coordinate so use for plotting legend
    x1s=c()
    y1s=c()
    for (r in 1:rs) {
      bw=.25
      bh=bw*(par()$usr[4]/(1/bw))
      x1=xl+bw*(seq(bxs.per.row[r])-1);x1s=c(x1s,x1)
      x2=xl+bw*(seq(bxs.per.row[r]))
      y1=-.025/bh-bh*(r-1);y1s=c(y1s,y1)
      y2=-.025/bh-bh*r
      rect(x1,y1,x2,y2,col=grey(rect.cols[(1:bxs.per.row[r])+4*(r-1)]))
    }
  }
  text(xs$xvals,-2*.025/bh-bh*rs,paste("N=",table(allvars[,cvar]),sep=""),
       adj=c(.5,1),cex=.8)
  # make x-offset center of bars
  xoffset=mean(x1s[2:3])-mean(xs$xvals)
  # make y-offset twice as far down as N=... text
  yoffset.val=-2*.025/bh-bh*rs-3*bh
  yoffset=y1s[1]-yoffset.val
  xs.text=x1s-xoffset
  ys.text=y1s-yoffset
  h.stretch=mean(xs$xvals)/diff(range(xs.text))
  v.stretch=.33*diff(par()$usr[3:4])/diff(range(xs.text))
  xs.text=xs.text*h.stretch
  ys.text=ys.text*v.stretch
  xoffset=mean(xs.text[2:3])-mean(xs$xvals)
  yoffset=ys.text[1]-yoffset.val
  xs.text=xs.text-xoffset
  ys.text=ys.text-yoffset
  text(xs.text,rep(ys.text,bxs.per.row[bxs.per.row>0]),
       names(var.means),adj=c(.5,.5),font=2,cex=.65)
  
  dev.off()
}

# vertical histograms
# for (i in 1:length(cnames)) {
#   hist.out=hist(allvars$MRR.log[allvars$features.cluster==cnames[1]],
#                 breaks=c(-4,seq(-0,6,by=.25)),plot=F)
#   reg=par()$usr
#   xi=xs$xvals[i]
#   par(fig=c(c(xi-.5,xi+.5)/diff(reg[1:2]),
#             c(reg[3],reg[4])/diff(reg[3:4])),new=T)
#   par(mar=c(0,0,0,0))
#   barplot(hist.out$counts,horiz=T,add=T,width=.1)
# }

# define categories
# cats.mat=matrix(,nrow(allvars.features),0)
# x=allvars.features
# clusters=cluster.out$cluster.names.for.allvars
# cnames=sort(unique(clusters))
# tol=1
# for (i in 1:length(cnames)) {
#   bounds=apply(x[clusters==cnames[i],],2,quantile,na.rm=T)[c(2,4),]
#   ci=apply(!apply(x,1,in.range,bounds=bounds),2,sum)<=tol
#   ci[is.na(ci)]=FALSE
#   cats.mat=cbind(cats.mat,ci)
# }