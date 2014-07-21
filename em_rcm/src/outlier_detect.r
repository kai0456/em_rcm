
######
##  Dectect whether x is outlier from N(mu,sgm)_j's, given eps
######

outlier.detect=function(x,mul,sgml,taul=1,eps=0.01){
  n=dim(x)[1]
  d=dim(x)[2]
  if(!is.matrix(mul)) return(print("mul is not matrix!"))
  if(!is.list(sgml)) return(print("sgml is not list!"))
  c=dim(mul)[1]
  if(length(taul)==1){taul=rep(1,c)/c}
  
  if(length(sgml)!=c) return(print("Size of sgml is not right!"))
  res=sapply(1:c,MD,x=x,mul=mul,sgml=sgml)
  colnames(res)=seq(1:c)
  res_class=sapply(1:n,max_indx,x=res,taul=taul,eps=eps)
  return(res_class)
}

#########
##
#########


##########

MD=function(j,x,mul,sgml){
  return(mahalanobis(x,center=mul[j,],cov=sgml[[j]]))
}

#########

max_indx=function(i,x,taul=taul,eps){
  d=dim(x)[2]
  indx=which.min(t(x[i,]))
  p2err=t(taul)%*%pchisq(x[i,],d,lower.tail=FALSE)  
  if(p2err<eps) return("outlier") else return(indx)
}

#########

