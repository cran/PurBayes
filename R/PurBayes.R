#PurBayes

write.PB<-function(fn.jags,prior=NULL,het=FALSE){
  if(is.null(prior)==TRUE)
    prior<-"dunif(0,1)"
  if(het==FALSE){ 
    model.string<-c("model{",
      "for(i in 1:N.snv){",
      "n.alt[i] ~ dbin(0.5*pur,n.tot[i])}",
      paste("pur ~ ",prior,"}",sep=""))
      }else{
    model.string<-c("model{",
      "for(j in 1:n.pop){",
      "lambda[j]~dunif(0,1)}",
      "lambda.srt<-sort(lambda)",
      "for(i in 1:N.snv){",
      "n.cat[i] ~ dcat(kappa)",
      "n.alt[i] ~ dbin(0.5*lambda.srt[n.cat[i]],n.tot[i])}",
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

PurBayes<-function(N,Y,pop.max=5,prior=NULL,burn.in=50000,n.post=10000,fn.jags="PB.jags",plot=FALSE){
  require(rjags)
  n.pop<-1
  write.PB(fn.jags,prior,het=FALSE)
  pb.dat.old<-list("N.snv"=length(N),"n.tot"=N,"n.alt"=Y)
  pb.m.old<-jags.model(file=fn.jags,pb.dat.old,n.chains=2,n.adapt=1000)
  update(pb.m.old,burn.in)
  dic.old<-dic.samples(pb.m.old,1000,type="popt")
  repeat{
    n.pop<-n.pop+1
    if(n.pop>pop.max){
      print("PurBayes has reached the defined pop.max, consider increasing this value")
      break
      }
    write.PB(fn.jags,prior,het=TRUE)
    pb.dat<-list("N.snv"=length(N),"n.tot"=N,"n.alt"=Y,"n.pop"=n.pop,"alpha"=rep(1,n.pop))
    pb.m<-jags.model(file=fn.jags,pb.dat,n.chains=2,n.adapt=1000)
    update(pb.m,burn.in)
    dic.new<-dic.samples(pb.m,1000,type="popt")
    dic.check<-diffdic(dic.old,dic.new)
    dic.diff<-sum(dic.check)
    dic.sd<-sqrt(length(dic.check))*sd(dic.check)
    if(dic.diff<0|dic.diff<dic.sd)
      break
    pb.m.old<-pb.m
    dic.old<-dic.new
    }
  n.pop.fin<-n.pop-1
  print(paste("PurBayes detected ",n.pop.fin," population(s) of variants",sep=""))
  if(n.pop.fin==1){
    pb.post<-coda.samples(pb.m.old,"pur",n.post)
    }else{
    pb.post<-coda.samples(pb.m.old,c("pur","kappa","lambda.srt"),n.post)
    }
  PB.out<-list("n.pop"=n.pop.fin,"PB.post"=pb.post)
  if(plot==TRUE)
    PB.plot(N,Y,PB.out)
  return(PB.out)
  }
