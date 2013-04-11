#PurBayes
library(rjags)

write.PB<-function(fn.jags,prior=NULL,het=FALSE,germ=FALSE){
  if(germ==TRUE){
    germ.string<-c("Z~dbin(q,M)",
      "q~dunif(0,1)")}else{
    germ.string<-""}
  if(is.null(prior)==TRUE)
    prior<-"dunif(0,1)"
  if(het==FALSE){ 
    model.string<-c("model{",germ.string,
      "for(i in 1:N.snv){",
      ifelse(germ==TRUE,"n.alt[i] ~ dbin(q*pur,n.tot[i])}","n.alt[i] ~ dbin(0.5*pur,n.tot[i])}"),
      paste("pur ~ ",prior,"}",sep=""))
      }else{
    model.string<-c("model{",
      germ.string,
      "for(j in 1:n.pop){",
      "lambda[j]~dunif(0,1)}",
      "lambda.srt<-sort(lambda)",
      "for(i in 1:N.snv){",
      "n.cat[i] ~ dcat(kappa)",
      ifelse(germ==TRUE,"n.alt[i] ~ dbin(q*lambda.srt[n.cat[i]],n.tot[i])}",
        "n.alt[i] ~ dbin(0.5*lambda.srt[n.cat[i]],n.tot[i])}"),
      "kappa ~ ddirich(alpha[])",
      "pur<-max(lambda)",
      "}"
      )
      }
	fileconn<-file(fn.jags)
  writeLines(model.string,fileconn)
  close(fileconn)
  }
  
PB.plot<-function(N,Y,out.PB){
  plot(N,Y,cex=0.75,pch=16,xlab="Total Reads",ylab="Mutant Allele Reads")
  n.pop<-out.PB$n.pop
  PB.post<-as.matrix(out.PB$PB.post)
  N.max<-max(N)
  if(n.pop==1){
    val.j<-quantile(as.matrix(PB.post),c(0.025,0.5,0.975))
    lines(c(0,N.max),c(0,N.max*0.5*val.j[2]))
    lines(c(0,N.max),c(0,N.max*0.5*val.j[1]),lty=2)
    lines(c(0,N.max),c(0,N.max*0.5*val.j[3]),lty=2)
    }else{
  for(i in 1:n.pop){
    pop.i<-ncol(PB.post)-1-n.pop+i
    val.j<-quantile(PB.post[,pop.i],c(0.025,0.5,0.975))
    lines(c(0,N.max),c(0,N.max*0.5*val.j[2]))
    lines(c(0,N.max),c(0,N.max*0.5*val.j[1]),lty=2)
    lines(c(0,N.max),c(0,N.max*0.5*val.j[3]),lty=2)
    }
  }
  }

PurBayes<-function(N,Y,M=NULL,Z=NULL,pop.max=5,prior=NULL,burn.in=50000,n.post=10000,fn.jags="PB.jags",plot=FALSE){
  n.pop<-1
  germ.dat<-ifelse(is.null(M)==F&&is.null(Z)==F,TRUE,FALSE)
  write.PB(fn.jags,prior,het=FALSE,germ=germ.dat)
  if(germ.dat==TRUE){
    pb.dat.old<-list("N.snv"=length(N),"n.tot"=N,"n.alt"=Y,"M"=sum(M),"Z"=sum(Z))}else{
    pb.dat.old<-list("N.snv"=length(N),"n.tot"=N,"n.alt"=Y)}
  
  pb.m.old<-jags.model(file=fn.jags,pb.dat.old,n.chains=2,n.adapt=1000)
  update(pb.m.old,burn.in)
  dic.old<-dic.samples(pb.m.old,1000,type="popt")
  pb.list<-list()
  dic.list<-list()
  pb.list[[1]]<-pb.m.old
  dic.list[[1]]<-dic.old
  
  pop.max.break<-0
  repeat{
    n.pop<-n.pop+1
    if(n.pop>pop.max){
      print("Warning: PurBayes has reached the defined pop.max, consider increasing this value")
      pop.max.break<-1
      break
      }
    write.PB(fn.jags,prior,het=TRUE,germ=germ.dat)
    if(germ.dat==TRUE){
      pb.dat<-list("N.snv"=length(N),"n.tot"=N,"n.alt"=Y,"M"=sum(M),"Z"=sum(Z),"n.pop"=n.pop,"alpha"=rep(1,n.pop))}else{
      pb.dat<-list("N.snv"=length(N),"n.tot"=N,"n.alt"=Y,"n.pop"=n.pop,"alpha"=rep(1,n.pop))}
    pb.m<-jags.model(file=fn.jags,pb.dat,n.chains=2,n.adapt=1000)
    update(pb.m,burn.in)
    dic.new<-dic.samples(pb.m,1000,type="popt")
    pb.list[[n.pop]]<-pb.m
    dic.list[[n.pop]]<-dic.new
    dic.check<-diffdic(dic.old,dic.new)
    dic.diff<-sum(dic.check)
  
    if(dic.diff<0)
      break
    pb.m.old<-pb.m
    dic.old<-dic.new
    }
  
  #Identify optimal model
  if(pop.max.break==1){
    n.pop.fin<-pop.max}else{
    
    model.check<-matrix(NA,nrow=n.pop,ncol=4)
    dic.sd<-sqrt(length(dic.check))*sd(dic.check)
    model.check[n.pop,]<-c(n.pop,sum(dic.new$deviance)+sum(dic.new$penalty),dic.diff,dic.sd)
    model.check[n.pop-1,]<-c(n.pop-1,sum(dic.old$deviance)+sum(dic.old$penalty),0,0)
    if(nrow(model.check)>2){
      for(i in 1:(n.pop-2)){
        dic.i<-dic.list[[i]]
        dic.check.i<-diffdic(dic.i,dic.old)
        dic.diff.i<-sum(dic.check.i)
        dic.sd.i<-sqrt(length(dic.check.i))*sd(dic.check.i)
        model.check[i,]<-c(i,sum(dic.i$deviance)+sum(dic.i$penalty),dic.diff.i,dic.sd.i)
        }
      n.pop.fin<-min(which(model.check[,3]<=model.check[,4]))
      }else{
      n.pop.fin<-1
      }
    }
  print(paste("PurBayes detected ",n.pop.fin," population(s) of variants",sep=""))
  pb.m.final<-pb.list[[n.pop.fin]]
  if(n.pop.fin==1){
    pb.post<-coda.samples(pb.m.final,"pur",n.post)
    }else{
    pb.post<-coda.samples(pb.m.final,c("pur","kappa","lambda.srt"),n.post)
    }
  PB.out<-list("n.pop"=n.pop.fin,"PB.post"=pb.post)
  if(plot==TRUE)
    PB.plot(N,Y,PB.out)
  return(PB.out)
  }