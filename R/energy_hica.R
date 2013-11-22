energy_hica <-
function(HICA.obj,maxcomp=1,nlevel=1,plot=TRUE){
 X=HICA.obj$X
 basis=HICA.obj$basis
 maxlev=dim(HICA.obj$aggregation)[1]
 n <- dim(X)[1]
 p <- dim(X)[2]
 if (maxcomp<1 || maxcomp > p){
  stop(paste("the number of components must to be between 1 and ",p,sep=""))
 }
 if (nlevel<1 || nlevel > maxlev){
  stop(paste("the number of levels must to be between 1 and ",maxlev,sep=""))
 }
 Xin=t(t(X)-colMeans(X))
 Etot=array(0,c(p,(p-1),maxcomp))
 Efin=matrix(0,maxcomp,(p-1))
 comp=matrix(0,maxcomp,(p-1))
 for (nn in 1:maxcomp){
  if (nn==1){
   for (i in (p-nlevel):(p-1)){
    for (j in 1:p){
     C1=t(cbind(basis[[i]][,j]))
     Etot[j,i,nn]=sum((Xin%*%t(C1)%*%solve(C1%*%t(C1)))^2)
    }
   }
   for (i in (p-nlevel):(p-1)){
    comp[nn,i]=which(rank(Etot[,i,nn])==p)
    Efin[nn,i]=Etot[comp[nn,i],i,nn] 
   }
  }
  if (nn>1){
   for (i in (p-nlevel):(p-1)){
    for (j in 1:p){
     if (sum(j==comp[,i])==0){
      C=NULL
      for (index in 1:(nn-1)){
       C=cbind(C,basis[[i]][,comp[index,i]])
      }
      C=t(cbind(C,basis[[i]][,j]))
      E=Xin%*%t(C)%*%solve(C%*%t(C))
      Etot[j,i,nn]=sum(E^2)
      for (r in 1:(nn-1)){
       for (c in (r+1):nn){
        Etot[j,i,nn]=Etot[j,i,nn]+2*sum(E[,r]*E[,c]*(C[r,]%*%C[c,]))
       }
      }
     }  
    }
   }
   for (i in (p-nlevel):(p-1)){
    comp[nn,i]=which(rank(Etot[,i,nn])==p)
    Efin[nn,i]=Etot[comp[nn,i],i,nn]
   }
  }
 }
 if (plot){
  vector.en=matrix(0,maxcomp,nlevel)
  for (i in 1:nlevel){
   for (j in 1:maxcomp){
    vector.en[j,i]=Efin[j,(p-nlevel+i-1)]/sum(Xin^2)
   }
  }
  leg=NULL
  for (k in 1:maxcomp){
   leg[k]=paste("K=",k,sep="")
  }
  dev.new()
  plot(((p-nlevel):(p-1)),vector.en[1,],type="b",col=1,lwd=2,ylim=c(0,(1+maxcomp/5)),xlab="Level",ylab="Energy",
  main="HICA energy",mgp=c(2.5,1,0))
  for (k in 2:maxcomp){
   points(((p-nlevel):(p-1)),vector.en[k,],type="b",col=k,lty=k,lwd=2)
  }
  legend((p-nlevel),(1+maxcomp/5),lty=1:maxcomp,col=1:maxcomp,lwd=3,legend=leg)
 }
 list(components=comp, energy=Efin, HICA.obj=HICA.obj)
}
