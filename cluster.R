cluster=function(x,dis.method="binary",c.method="ward.D",
                 Nc.max=10,Nc.set=NULL,min.cluster.size=1,Nc.criterion="within.cluster.ss",
                 trim=T,cutoff=.33,cutoff.abs=NULL,ivars=NULL,
                 res=250,plot.folder=paste(getwd(),"/",sep=""),
                 binary.summary=F,save.out=T) {
  
  # if set the cluster size create new folder for that number
  if (!is.null(Nc.set)) {
    plot.folder.sub=paste(plot.folder,"Nc=",Nc.set,"/",sep="")
    if (!dir.exists(plot.folder.sub)) {
      dir.create(plot.folder.sub)
    }
    plot.folder=plot.folder.sub
  }
  
  # load necessary packages
  library(gclus)
  library(fpc) # http://www.rdocumentation.org/packages/fpc/functions/fpc-package
  library(NbClust) # http://www.jstatsoft.org/v61/i06/paper
  library(sp)
  library(rgeos)
  library(dplyr)
  
  # output
  dis.mat=dist(x,method=dis.method)
  clust.sol=hclust(dis.mat,method=c.method)
  clust.sol.ord=reorder(clust.sol,dis.mat)
  clust.ord=clust.sol.ord$order
  sort.by.var.means=order(apply(x,2,sum))
  
  # only do these things if you are not asking for a specific cluster,
  # in which case I assume you've already run the general Nc algorithm
  if (is.null(Nc.set)) {
    
    # -----------------------------
    # create dissimilarity matrix
    # -----------------------------
    #   if (nrow(x)>125) {
    #     samp.rows=sample(1:nrow(x),125)
    #   } else {
    #     samp.rows=1:nrow(x)
    #   }
    cmat=dmat.color(as.matrix(dis.mat)[clust.ord,clust.ord],
                    col=grey(seq(.1,.9,length.out=8)))
    # plot it
    # max character width/height is appx .1/.12 inch
    # to fit labels into .05 inch matrix row height, divide cell width by max height and make 75% smaller
    cell.width=15/nrow(x)
    cex=(cell.width/.12)*.75
    rowlabs=rownames(as.matrix(dis.mat))[clust.ord]
    rowlabs.inches=max(nchar(rowlabs))*.1*cex
    width.mat=nrow(cmat)*cell.width+rowlabs.inches+.1
    height.mat=nrow(cmat)*cell.width
    png(paste(plot.folder,"dismat ",c.method,".png",sep=""),
        res*width.mat,res*height.mat,res=res)
    par(bg="transparent")
    par(mai=c(.1,rowlabs.inches,.1,.1))
    par(xpd=T)
    plotcolors(cmat)
    text(0,1:nrow(cmat),rowlabs,cex=cex,adj=c(1,.5))
    dev.off()
    
    # -----------------------------
    # create dendogram
    # -----------------------------
    if (nrow(x)<250) {
      png(paste(plot.folder,"dendo ",c.method,".png",sep=""),
          res*nrow(x)/15,res*nrow(x)/45,res=res)
      par(bg="transparent")
      par(mar=c(1,0,1,0))
      plot(clust.sol.ord,cex=.2,xlab="",ylab="",xaxt="n",yaxt="n")
      dev.off()
    }
    
    # find best number of clusters based on wb.ratio
    # my favorites
    # get.stats=c("average.between","average.within","within.cluster.ss","avg.silwidth",
    #             "pearsongamma","dunn","dunn2","entropy","wb.ratio","ch","sindex","min.cluster.size")
    # but go with within.cluster.ss since it is easiest to interpret
    # also record min.cluster.size for filtering
    get.stats=c(Nc.criterion,"min.cluster.size")
    stats.mat=matrix(,0,length(get.stats))
    for (ci in 2:Nc.max) {
      clust.sol.cut=cutree(clust.sol,k=ci)
      stats=cluster.stats(dis.mat,clust.sol.cut)
      stats.mat=rbind(stats.mat,as.numeric(unlist(stats[get.stats])))
    }
    stats.mat=as.data.frame(as.matrix(stats.mat))
    colnames(stats.mat)=get.stats
    
    # ----------------------------------------------------------------
    # choose number of clusters using elbow point of average between-cluster dissimilarity
    # first take only points that satisfy minimum cluster size requirement
    # ----------------------------------------------------------------
    wcss=stats.mat$within.cluster.ss[stats.mat$min.cluster.size>=min.cluster.size]
    # define xs and ys
    ys=wcss
    xs=seq(from=min(ys),to=max(ys),length.out=length(ys))
    wcss.line.slope=(last(ys)-first(ys))/(last(xs)-first(xs))
    wcss.line.intercept=first(ys)-wcss.line.slope*first(xs)
    
    # find a, b, and c for ax+by+c=0
    # based on y=Mx+B, -Mx+y-B=0, so a=-M, b=1, and c=-B
    a=-wcss.line.slope
    b=1
    c=-wcss.line.intercept
    # get coordinates of closest point on line
    x0=(b*(xs-a*ys)-a*c)/(a^2+b^2)
    y0=(a*(-b*xs+a*ys)-b*c)/(a^2+b^2)
    
    y2=last(ys)
    y1=first(ys)
    x2=last(xs)
    x1=first(xs)
    dists=abs((y2-y1)*xs-(x2-x1)*ys+x2*y1-y2*x1)/
      sqrt((y2-y1)^2+(x2-x1)^2)
    max.pt=which.max(dists)
    Nc=max.pt+1
    
    png(paste(plot.folder,"optimal Nc ",c.method,".png",sep=""),
        5*res,5*res,res=res)
    par(bg="transparent")
    par(mar=c(4,4,1,1))
    plot(xs,ys,type="b",xaxt="n",yaxt="n",
         ylab="within cluster sums of squares",xlab="number of clusters")
    axis(1,xs,1:length(xs)+1,cex.axis=.8)
    axis(2,cex.axis=.8)
    segments(first(xs),first(ys),last(xs),last(ys))
    segments(xs[max.pt],ys[max.pt],x0[max.pt],y0[max.pt],col="red",lwd=2)
    dev.off()
  }
  
  # allow user to override optimal number of clusters
  if (!is.null(Nc.set)) {Nc=Nc.set}
  
  # recreate basic objects with new orders and clusters
  clust.sol.cut.ord=cutree(clust.sol,k=Nc)[clust.ord]
  clust.vals=unique(clust.sol.cut.ord)
  x.cord=x[clust.ord,sort.by.var.means]
  
  # ----------------------------------------------------------------
  # create table for discriminating clusters based on variables
  # ----------------------------------------------------------------
  # record means (or percentages) of variables for each cluster
  x.cmeans=matrix(,0,ncol(x.cord))
  for (i in clust.vals) {
    x.cmeans=rbind(x.cmeans,round(apply(x.cord[clust.sol.cut.ord==i,],2,mean),2))
  }
  colnames(x.cmeans)=colnames(x.cord)
  rownames(x.cmeans)=paste("cluster",clust.vals)
  
  var.discr=apply(x.cmeans,2,sd)
  var.freq=apply(x.cord,2,mean)
  var.freq.max=apply(x.cmeans,2,max)
  var.discr.df=data.frame(discr=round(var.discr,2),
                          freq=round(var.freq,2),
                          freq.max=round(var.freq.max,2))
  var.discr.df=var.discr.df[order(var.discr.df$disc,decreasing=T),]
  
  bar.width=.2
  cex=(bar.width/.12)*.5
  barlabs=rownames(var.discr.df)
  barlabs.inches=max(nchar(barlabs))*.1*cos(.25*pi)*cex
  width.bp=ncol(x)*bar.width+4*.2
  height.bp=3+barlabs.inches+.1
  png(paste(plot.folder,"feature variance across clusters.png",sep=""),
      res*width.bp,res*height.bp,res=res)
  par(mai=c(barlabs.inches,4*.2,.1,.1))
  par(xpd=T)
  xs=barplot(var.discr.df$discr,ylab="features by variance across clusters")
  text(xs,-.02*max(var.discr.df$discr),barlabs,
       cex=cex,srt=45,adj=c(1,.5))
  dev.off()
  
  # ---------------------------------------
  # plot all data with clusters overlayed
  # ---------------------------------------
  # record what clusters have sufficient non-0 or above-median value for each variable
  # for non-binary data, if cutoff.abs is NULL, use the variable mean as the cutoff
  if (dis.method!="binary" & is.null(cutoff.abs)) {
    cutoff=apply(x.cord,2,mean)
  } else {
    if (!is.null(cutoff.abs)) {
      cutoff=cutoff.abs
      # repeat the same value for each variable if only one is given
      cutoff=rep(cutoff.abs,ncol(x))[1:ncol(x)]
    }
  }
  abv_cutoff=x.cmeans>matrix(cutoff,nrow(x.cmeans),ncol(x.cmeans),byrow=T)
  trim.cols=apply(abv_cutoff,2,sum)==0
  # I think for non-binary data, 1+ clusters should pass the above cutoff for each variable
  # but just in case, make sure not to exclude any variables
  if (dis.method!="binary" | trim==F) {
    trim.cols=rep(FALSE,ncol(x.cord))
    names(trim.cols)=colnames(abv_cutoff)
  }
  x.cord=x.cord[,!trim.cols]
  
  # reorder the vars so that they group according to what clusters use them
  eval(parse(text=paste("all.combos=expand.grid(",
                        paste(rep("c(0,1),",length(clust.vals)-1),collapse=""),
                        "c(0,1))",sep="")))
  all.combos=all.combos[order(apply(all.combos,1,sum)),]
  vbc.long=t(abv_cutoff[,!trim.cols])*1
  vbc.new=matrix(,0,length(clust.vals))
  for (i in 1:nrow(all.combos)) {
    combo.matches.data=as.numeric(all.combos[i,])==as.matrix(t(vbc.long))
    combo.matches=apply(combo.matches.data,2,sum)==length(clust.vals)
    vbc.new=rbind(vbc.new,vbc.long[combo.matches,])
    if (sum(combo.matches)==1) {
      rownames(vbc.new)[nrow(vbc.new)]=rownames(vbc.long)[combo.matches]
    }
  }
  # reorganize x.cord based on vbc ordering
  vbc.ord=data.frame(var=rownames(vbc.new),
                     vbc.ord=1:nrow(vbc.new))
  x.cols=data.frame(var=colnames(x.cord),orig.ord=1:ncol(x.cord))
  x.cols=merge(x.cols,vbc.ord,by="var",all.x=T)
  x.cols=arrange(x.cols,orig.ord)
  x.cord=x.cord[,order(x.cols$vbc.ord)]
  # reorganize x.cmeans to match
  x.cmeans=x.cmeans[,colnames(x.cord)]
  
  # first, try to graph with square cells
  # determine plot width, which is 15'' plus necessary margin for rowlabels
  cell.width=15/nrow(x.cord)
  cex=(cell.width/.12)*.5
  row.cex=cex;col.cex=cex
  rowlabs=colnames(x.cord)
  rowlabs.inches=max(nchar(rowlabs))*.1*cex
  width.clust=nrow(x.cord)*cell.width+rowlabs.inches+.1
  # determine plot height
  if (is.null(ivars)) {
    if (is.null(rownames(x.cord))) {
      collabs=NULL
      collabs.inches=0
    } else {
      collabs=rownames(x.cord)
      collabs.inches=max(nchar(collabs))*.1*cex
    }
  } else {
    collabs=ivars
    collabs.strwidths=apply(nchar(as.matrix(collabs)),2,max,na.rm=T)
    collabs.inches=sum(collabs.strwidths)*.1*cex
  }
  height.clust=ncol(x.cord)*cell.width+.12+.1
  # if square cells result in width>2*height, then set 2/1 ratio
  # change rowlabels so they fit just within implied cell height
  if (width.clust/height.clust>4) {
    cell.height=(15/4)/ncol(x.cord)
    row.cex=(cell.height/.12)*.5
    rowlabs.inches=max(nchar(rowlabs))*.1*row.cex
    width.clust=nrow(x.cord)*cell.width+rowlabs.inches+.1
    col.cex=cex
    height.clust=15/2+collabs.inches+.1
  }
  png(paste(plot.folder,"clustered vars by role"," ",c.method,".png",sep=""),
      res*width.clust,res*height.clust,res=res)
  par(bg="transparent")
  par(mai=c(collabs.inches,rowlabs.inches,.1,.1))
  par(xpd=T)
  image(1:nrow(x.cord),1:ncol(x.cord),as.matrix(x.cord),
        col=grey(seq(1,0,length.out=10)),xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
  text(0,1:ncol(x.cord),colnames(x.cord),adj=c(1,.5),cex=row.cex)
  clust.vals=unique(clust.sol.cut.ord)
  rgb.clusts.mat=matrix(,0,3)
  for (j in 1:length(clust.vals)) {
    rgb.clusts=runif(3)
    rgb.clusts.mat=rbind(rgb.clusts.mat,rgb.clusts)
    rect(min(which(clust.sol.cut.ord==j))-.5,.5,
         max(which(clust.sol.cut.ord==j))+.5,ncol(x.cord)+.5,col=rgb(rgb.clusts[1],
                                                                     rgb.clusts[2],
                                                                     rgb.clusts[3],alpha=.5))
  }
  # add all item-level variables
  if (!is.null(ivars)) {
    # set how much space to add between variables
    ivars.strwidths=collabs.strwidths
    ivars.strwidths.mat=matrix(ivars.strwidths[1:(ncol(ivars)-1)],
                               ncol(ivars)-1,ncol(ivars)-1)
    ivars.strwidths.mat[lower.tri(ivars.strwidths.mat,diag=F)]=0
    ivars.strwidths.sums=apply(ivars.strwidths.mat,2,sum)
    ivars.sep=c(0,-ivars.strwidths.sums*cex*1.5)
    for (i in 1:3) {
      text(1:nrow(x),ivars.sep[i],ivars[clust.ord,i],
           adj=c(1,.5),srt=90,cex=cex)
    }  
  }
  dev.off()
  
  # -----------------------------------------------
  # create summary verion of variables by cluster
  # -----------------------------------------------
  if (binary.summary==T) {
    x.sum=t(vbc.new)
  } else {
    x.sum=x.cmeans[,!as.vector(trim.cols[colnames(t(vbc.new))])]
  }
  cell.width=.2
  cex=(cell.width/.12)*.5
  rowlabs=colnames(x.sum)
  rowlabs.inches=max(nchar(rowlabs))*.1*cex
  lmar=rowlabs.inches+.1
  width.clust.sum=nrow(x.sum)*cell.width*1.5+lmar
  collabs=paste("N=",table(clust.sol.cut.ord)[unique(clust.sol.cut.ord)],sep="")
  collabs.inches=max(nchar(collabs))*.1*cex
  bmar=collabs.inches+.1
  height.clust.sum=ncol(x.sum)*cell.width+bmar
  png(paste(plot.folder,"clustered vars by role sum"," ",c.method,".png",sep=""),
      res*width.clust.sum,res*height.clust.sum,res=res)
  par(bg="transparent")
  par(mai=c(bmar,lmar,.1,.1))
  par(xpd=T)
  image(1:nrow(x.sum),1:ncol(x.sum),x.sum,
        col=grey(seq(1,.25,length.out=10)),xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
  text(.5-.05*nrow(x.sum),1:ncol(x.sum),colnames(x.sum),adj=c(1,.5),cex=cex)
  for (j in 1:length(clust.vals)) {
    rect(which(clust.vals==j)-.5,.5,
         which(clust.vals==j)+.5,ncol(x[,!trim.cols])+.5,col=rgb(0,
                                                                 0,
                                                                 0,alpha=.1))
  }
  text(1:nrow(x.sum),.5-.05*ncol(x.sum),collabs,
       adj=c(1,.5),cex=cex,srt=45)
  dev.off()
  
  # create cluster names
  cluster.names=c()
  x.sum=t(vbc.new)
  for (i in 1:nrow(x.sum)) {
    cluster.names=c(cluster.names,
                    paste(colnames(x.sum)[x.sum[i,]==1],sep="",collapse="_"))  
  }
  if (length(cluster.names)!=length(unique(cluster.names))) {
    cluster.names=paste(1:length(cluster.names),cluster.names,sep="__")
  }
  clust.sol.cut=cutree(clust.sol,Nc)
  cluster.names=data.frame(cluster.names=cluster.names,
                           cluster.indexes=unique(clust.sol.cut[clust.ord]))
  # verify correct cluster names
  # round(apply(allvars[clust.sol.cut==1,
  #                     grepl("role",colnames(allvars))]==1,2,mean,na.rm=T))
  cluster.names.ord=select(arrange(cluster.names,cluster.indexes),cluster.names)
  clust.sol.cut.factor=as.factor(clust.sol.cut)
  levels(clust.sol.cut.factor)=as.matrix(cluster.names.ord)
  cluster.names.for.allvars=as.character(clust.sol.cut.factor)
  # verify correct cluster again
  # round(apply(allvars[clust.sol.cut.factor=="Owner",
  #                     grepl("role",colnames(allvars))]==1,2,mean,na.rm=T))
  
  # revert to correct x.sum for output
  if (binary.summary==T) {
    x.sum=t(vbc.new)
  } else {
    x.sum=x.cmeans[,!as.vector(trim.cols[colnames(t(vbc.new))])]
  }
  
  out=list(
    sort.by.var.means=sort.by.var.means,
    dis.mat=dis.mat,
    clust.sol=clust.sol,
    clust.sol.ord=clust.sol.ord,
    clust.ord=clust.ord,
    x.sum=x.sum,
    cluster.names.for.allvars=cluster.names.for.allvars
  )
  if (save.out==T) {
    return(out)
  }
}
