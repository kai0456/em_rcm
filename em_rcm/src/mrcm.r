source("src/em_rcm.r")

mrcm=function(x){
  x=as.matrix(x)
  temp=dim(x)
  n=temp[1]
  d=temp[2]
  r=fwtdr(x)
  rcm=t(r)%*%r
  evec=eigen(rcm)$vectors
  eval=vector()
  for(i in 1:d){
    val=x%*%evec[,i]
    eval[i]=(mad(val))^2
  }
  sgm_diag=diag(eval)
  return(evec%*%sgm_diag%*%t(evec))
}

