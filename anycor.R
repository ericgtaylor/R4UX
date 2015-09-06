# conduct test of x/y relationship for all combinations of variable classes
anycor=function(x,y,is.DV=NULL,data) {
  x.class=class(data[,x])
  y.class=class(data[,y])
  # factor x factor: run a Chi-squared test
  if (x.class=="factor" & y.class=="factor") {
    rc.table=table(data[,c(x,y)])
    test.out=chisq.test(rc.table)
    r2=test.out$stat/sum(rc.table)
    p.val=test.out$p.val
    test.type="chisq"
    beta.val=NA
    if (nrow(rc.table)==ncol(rc.table) & nrow(rc.table)==2) {
      props=rc.table/matrix(apply(rc.table,1,sum),2,2)
      effect.mag=props[2,2]-props[1,2]
    }
  }
  # factor x numeric
  if ((x.class=="factor")+(y.class=="factor")==1) {
    # ensure means are calculable
    means1=eval(parse(text=paste("with(data,tapply(",y,",",x,",mean,na.rm=T))",sep="")))
    means2=eval(parse(text=paste("with(data,tapply(",x,",",y,",mean,na.rm=T))",sep="")))
    if (sum(is.na(means1))==0 | sum(is.na(means2))==0) {
      if (is.DV=="" | is.DV=="xy") {
        # if no preference for DV, choose numeric
        # factor x numeric: run an anova
        DV=c(x,y)[1+(x.class=="factor")*1]
        IV=c(x,y)[1+(y.class=="factor")*1]
        test=lm(data[,DV]~data[,IV])
        test.out=anova(test)
        r2=test.out$Sum[1]/sum(test.out$Sum)
        p.val=test.out$Pr[1]
        test.type="anova"
        beta.val=coef(test)[2]
        eval(parse(text=paste("props=with(data,tapply(",y,",",x,",mean,na.rm=1))",sep="")))
        if (sum(unique(data[,x])%in%c("no","yes"))==2) {
          effect.mag=props["yes"]-props["no"]
        }
      }
      else {
        DV=c(x,y)[1+(is.DV=="y")*1]
        IV=c(y,x)[1+(is.DV=="y")*1]
        DV.levels=unique(data[,DV])
        DV.levels=DV.levels[!is.na(DV.levels)]
        # factor x numeric (but need factor as DV): logistic regression
        if (length(DV.levels)==2) {
          test=glm(data[,DV]~data[,IV],family="binomial")
          test.out=summary(test)
          r2=1-(test$deviance/-2)/(test$null.deviance/-2)
          p.val=coef(test.out)[2,4]
          test.type="logreg"
          beta.val=coef(test)[2]
          IV.vals=apply(data,2,unique)[colnames(data)==IV]
          IV.vals.range=range(as.numeric(unlist(IV.vals)),na.rm=1)
          if ((IV.vals.range[1]==1 & IV.vals.range[2]==6) & all.equal(sort(unique(data$y)),c("no","yes"))) {
            eval(parse(text=paste("props=with(data,tapply((",y,"=='yes')*1,",x,">3,mean,na.rm=1))",sep="")))
            effect.mag=props[names(props)==TRUE]-props[names(props)==FALSE]
          }
        }
        # factor x numeric (but need numeric as DV): regression
        else {
          test=lm(data[,DV]~data[,IV])
          test.out=anova(test)
          r2=test.out$Sum[1]/sum(test.out$Sum)
          p.val=test.out$Pr[1]
          test.type="anova"
          beta.val=coef(test)[2]
          if (all.equal(sort(unique(data$y)),c("no","yes"))) {
            props=eval(parse(text=paste("with(data,tapply(",y,",",x,",mean,na.rm=1))",sep="")))
            effect.mag=props["yes"]-props["no"]
          }
        }
      }
    }
    else {
      test.out=NA
      r2=NA
      p.val=NA
      test.type=NA
      beta.val=NA
    }
  }
  # numeric x numeric: regression
  if (!(x.class=="factor") & !(y.class=="factor")) {
    xy=matrix(c(data[,x],data[,y]),nrow(data),2)
    xy=xy[!is.na(apply(xy,1,sum)),]
    if (length(xy)>2) {
      if (length(xy)<40 | sum(apply(xy,2,var)==0)) {
        test.out=NA
        r2=NA
        p.val=NA
        test.type=NA
        beta.val=NA
      }
      else {
        test=lm(xy[,1]~xy[,2])
        test.out=summary(test)
        r2=cor(xy[,1],xy[,2])^2
        p.val=test.out$coefficients[2,4]
        test.type="linreg"
        beta.val=coef(test)[2]
      } 
    }
    else {
      test.out=NA
      r2=NA
      p.val=NA
      test.type=NA
      beta.val=NA
    }
  }
  
  # return effect size, p-value, and test.out
  return(list(r2=r2,
              p.val=p.val,
              test.out=test.out,
              test.type=test.type,
              beta.val=beta.val))
}