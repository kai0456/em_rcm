## This weighted Spatial Rank function is pretty inefficient
## But it helps to replace fwtdr function in em_rcm.r in case we cann't
## run the "wtdr.dll" under non-i386 R system.

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