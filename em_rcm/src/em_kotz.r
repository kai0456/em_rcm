library(mvtnorm)  #Multivariate Normal and t Distributions.
library(e1071) # fuzzy k means


########
### weighted rank of x
########
#dyn.load("src/wtdr.dll")

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
#########
### End weighted rank of x
########




### Weighted Rank function written in R code
fwtdr=function(x,w){
  n=dim(x)[1]
  d=dim(x)[2]
  s=sum(w)
  fwtdr=sapply(1:n,fwtdr_ith,x=x,w=w,n=n,d=d)
  fwtdr=fwtdr/s
  fwtdr=t(fwtdr)
  return(fwtdr)
}
#########
fwtdr_ith=function(i,x,w,n,d){
  xmtx=matrix(x[i,],nrow=n,ncol=d,byrow=T)
  xmtx=xmtx-x
  rm(x)
  nomx=eucldisc(xmtx)
  invnomx=1/nomx
  rm(nomx)
  invnomx[!is.finite(invnomx)]=0
  invnomx=matrix(invnomx,n,d)
  fwtdr_ith=t(w)%*%(xmtx*invnomx) ## haven't divide by s
  return(fwtdr_ith)
}
##### Weighted Rank function written in R code

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

mix_kotz=function(x,mul,sgml,taul){
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
  
    for (j in 1:c){   # update median and covariance
  	
  	sgmj = sgml[[j]]
  	eve = eigen(sgmj)
  	evec = eve$vectors
  	evev = eve$values
  	y = x%*%evec%*%diag(evev^(-1/2))%*%t(evec) 	   # standardize
    wtdr=fwtdr(y,Tij[,j])
    ctindx=fctindx(wtdr)
    mul[j,]= t(x[ctindx,])
       
   ## do loop for covariance 
   	 s=sum(Tij[,j])    
     x_mius_mu=t(x)-x[ctindx,]   
     for (k in 1:100) {
      y_mius_ymu = t(y)-y[ctindx,]
     dis=eucldisc(t(y_mius_ymu))
     invbottom=dis^(-1)
     invbottom[!is.finite(invbottom)]=0    
     sgmj1 = x_mius_mu %*% diag(Tij[,j]/s*invbottom)%*% t(x_mius_mu)+0.00000001*diag(1,d)
    a =norm(sgmj1-sgmj)
     	 if (a<0.0001)
      {break}
      sgmj = sgmj1
      eve = eigen(sgmj)
  	evec = eve$vectors
  	evev = eve$values
  	y = x%*%evec%*%diag(evev^(-1/2))%*%t(evec)    	
     }  
    sgml[[j]] =sgmj1  
          
    
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

em_kotz_GM=function(x,mul,sgml=0,taul=0){
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
  for (t in 1:100){
    mix=mix_kotz(x,mull[[t]],sgmll[[t]],taull[[t]])
    mull=c(mull,list(mix[[1]]))
    sgmll=c(sgmll,list(mix[[2]]))
    taull=c(taull,list(mix[[3]]))
    #if (abs(taull[[t+1]][1]-taull[[t]][1])<0.0001)
    if (abs(taull[[t+1]][1]-taull[[t]][1])<0.001)
      {break}
  }
  ## (t+1)=ending position
  list(mull[[t+1]],sgmll[[t+1]],taull[[t+1]],t+1,mull,sgmll,taull)

}