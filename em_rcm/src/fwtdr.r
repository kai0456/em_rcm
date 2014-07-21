#####################
## Date:06/08/2014
## Functions:
##  Calculate weighted rank of x
#####################

vec_norm=function(x){
  x=as.vector(x)
  return(sqrt(sum(x^2)))
}

spatial_sign=function(x){
  x=as.vector(x)
  x_norm=vec_norm(x)
  if(x_norm == 0){
    return(rep(0,length(x)))
  } else {
    return(x/x_norm)
  }
}

fwtdr=function(dat,w=1){
  n=dim(dat)[1]
  p=dim(dat)[2]
  if (length(w)==1){
    w=rep(1/n,n)
  }
  dat=as.matrix(dat)
  

  spatialRank=matrix(NA,n,p)
  for (i in 1:n){
    x=matrix(dat[i,],n,p,byrow=TRUE)
    s=x-dat
    s_list=list()
    for (k in 1:n){
      s_list[[k]]=as.vector(s[k,])
    }
    spatialSign=sapply(s_list,spatial_sign)
    spatialSign=t(spatialSign)
    spatialRank[i,]=w%*%spatialSign
  }
  return(spatialRank)
}
