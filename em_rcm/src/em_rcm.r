library(mvtnorm)  #Multivariate Normal and t Distributions.
library(e1071) # fuzzy k means


########
### weighted rank of x
########
#dyn.load("src/wtdr.dll")
#
#fwtdr=function(x,w=1){
#  x=as.matrix(x)
#  l=dim(x)
#  n=l[1]
#  d=l[2]
#  if (length(w)==1){
#    w=rep(1,n)/n    
#  }
#  d=dim(x)
#  if(any(d[2]==1,is.vector(x))) return (rank(x))
#  tmp=.C("wtdr",as.double(x),as.integer(d),as.double(w),res=double(d[1]*d[2]))$res
#  tmp=matrix(tmp,ncol=d[2])
#  return(tmp)
#}

## Using pure R version of fwtdr
source("src/fwtdr.r")

#########
### End weighted rank of x
########

########
### Spatial Rank  Pseudo-Median index
########
fctindx=function(r){ # (rank) center index
  out=eucldisc(r)
  fctindx=which.min(out)
  fctindx
}

########
### End Spatial Rank  Pseudo-Median
########

########
### euclidean disitance, give a matirx n*d, return n-vector distance from center
########
eucldisc=function(x){
  u=dim(x)
  n=u[1]
  d=u[2]
  sqx=x*x
  eucldisc=apply(sqx,1,sum)
  eucldisc=sqrt(eucldisc)
  eucldisc=as.vector(eucldisc)
  eucldisc
}

########
### End euclidean disitance
########

#######
##  Use fwtdr with defaul w=1/(dim(x)[1]) return spatial rank and find the
##  Pseudo Spatial Median
#######
fwtdr.center=function(x,w=1){
    #w=1 will go into fwtdr and therefore initialize as w=c(1/n,...,1/n)
    wtdr=fwtdr(x,w)
    ctindx=fctindx(wtdr)
    mu=x[ctindx,]
    mu=as.double(mu)
    return(mu)
}

########
##  End fwtdr.center
########



#########
### E-step, return  P(Z_i=j|X,theta) as n vector  for Gussian Dist.
#########
fTij=function(x,mul,sgml,taul){
  #mul is c*d matrix
  #sgml is list
  #taul is vector
  n=dim(x)[1]
  d=dim(x)[2]
  c=dim(mul)[1]   #number of components

  fTij=matrix(0,nrow=n,ncol=c)
  dstmtx=dmvnorm(x,mul[1,],sgml[[1]])
  for (j in 2:c){
    dst=dmvnorm(x,mul[j,],sgml[[j]])
    dstmtx=cbind(dstmtx,dst)
  }
  dstmtx[dstmtx==0]=1 ## Strange! But just make really far 
                      ## outlier have the same density w.r.t different components
                      ## the 1 doesn't have to be 1, can be any constant.
  wtddstmtx=dstmtx%*%diag(taul)
  denom=apply(wtddstmtx,1,sum)  #return as a n*1 vector
  for (j in 1:c){
    fTij[,j]=wtddstmtx[,j]/denom  #return as n*1 vector
  }
  fTij
}
#########
### End E-step
#########

#########
### M-Step
#########

mix_rcm=function(x,mul,sgml,taul){
  #mul is c*d matrix
  #sgml is list
  #taul is vector
  n=dim(x)[1]
  d=dim(x)[2]
  c=dim(mul)[1]   #number of components
  #if (sum(taul)!= 1){print("taul not sum up 1")}
  if (c!=length(sgml)| c!=length(taul)){
    print("Number of component doesn't match!")
    break
  }
  Tij=fTij(x,mul,sgml,taul)
  
  write(paste("Tij=",head(Tij,1),"date=",date(),sep=" "),'error.log',append=T)

  for (j in 1:c){    # update mean and covariance
    wtdr=fwtdr(x,Tij[,j])
	
	write(paste("wtdr=",length(wtdr),"date=",date(),sep=" "),'error.log',append=T)
    
	ctindx=fctindx(wtdr)
	
	write(paste("j=",j,"ctindx=",ctindx,"date=",date(),sep=" "),'error.log',append=T)
	
    mul[j,]=x[ctindx,]

    ## Use weighted MRCM find Sigma (eigenvalu=mean abs deviation)
    s=sum(Tij[,j])
    wtdrcm=t(diag(Tij[,j]/s)%*%wtdr)%*%wtdr
    evec=eigen(wtdrcm)$vectors
    lambda=rep(0,d)

    ##Use Kai's weighted lambda
    jthcomp_indx=floor(n*(1-s/n))+1
    for (k in 1:d){
      med=t(mul[j,])%*%evec[,k]
      wtdprj=Tij[,j]*(x%*%evec[,k]-rep(med,n))
      abswtdprj=abs(wtdprj)
      wtdprj_jthcomp=wtdprj[which(abswtdprj>=sort(abswtdprj)[jthcomp_indx])]
      lambda[k]=median(abs(wtdprj_jthcomp))*1.4826 #1.4826 consistent coef.

      #wtdTij=Tij[,j]*(1-eucldisc(wtdr))
      #swtdTij=sum(wtdTij)
      #wtdlambda=t(wtdTij/swtdTij)%*%abs(x%*%evec[,k]-rep(t(mul[[j]])%*%evec[,k]))
      #lambda[k]=wtdlambda*1.5
    }
    ##End Use Kai's weightd lambda


#   ##Use Dr. Dang's weighted lambda
#   for (k in 1:d){
#      med=t(Tij[,j]/s)%*%x%*%evec[,k]
#      lambda[k]=t(Tij[,j]/s)%*%abs(x%*%evec[,k]-rep(med,n))*1.253 #1.253 consitent coef
#   }
#   ## End Use Dr. Dang's weighted lambda
    mrcm=evec%*%diag(lambda^2)%*%t(evec)


    sgml[[j]]=mrcm+0.00000001*diag(1,d)
  }
  taul=apply(Tij,2,mean)
  list(mul,sgml,taul)
}

########
### End M-step
########

########
###  em iteration;  criteria (loops<100) & (taul[1] stable)
########

em_rcm_GM=function(x,mul,sgml=0,taul=0){
  n=dim(x)[1]
  d=dim(x)[2]
  c=dim(mul)[1]

  ## default sgml & taul
  if (!is.list(sgml)){
    rm(sgml)
    sgm=diag(1,d)
    sgml=list(sgm)
    for (k in 2:c){
      sgml=c(sgml,list(sgm))
    }
  }
  if (length(taul)==1){
    rm(taul)
    taul=rep(1/c,c)
  }

  ##
  mull=list(mul)   # 2 levels
  sgmll=list(sgml) # 3 levels
  taull=list(taul) # 2 levels
  for (t in 1:20){
	print(paste('The ',t,' iteration',proc.time()[3]))
  
    mix=mix_rcm(x,mull[[t]],sgmll[[t]],taull[[t]])
    mull=c(mull,list(mix[[1]]))
    sgmll=c(sgmll,list(mix[[2]]))
    taull=c(taull,list(mix[[3]]))
    if (abs(taull[[t+1]][1]-taull[[t]][1])<0.0001)
      {break}
  }
  ## (t+1)=ending position
  list(mull[[t+1]],sgmll[[t+1]],taull[[t+1]],t+1,mull,sgmll,taull)

}