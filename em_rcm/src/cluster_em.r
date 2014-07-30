source("src/em_rcm.r")
source("src/em_reg.r")

#Does not contain outlier detector.

cluster_em_rcm=function(x,mul,sgml=0,taul=0,outlier=TRUE,epsilon=.01){
  n=dim(x)[1]
  d=dim(x)[2]
  c=dim(mul)[1]
  
  x=as.matrix(x)
  if (outlier){
    em=em_rcm_GM(x,mul,sgml,taul)
    mul=em[[1]]
    sgml=em[[2]]
  
    res=integer(n)
    dsty=double(c)
    for (i in 1:n){
      for(j in 1:c){
        dsty[j]=dmvnorm(x[i,],mul[j,],sgml[[j]])
        belong=which.max(dsty)
        res[i]=belong
      }
    }
    return(list(cbind(x,res),mul,sgml))
   }
  
}

############

cluster_em_reg=function(x,mul,sgml=0,taul=0){
  n=dim(x)[1]
  d=dim(x)[2]
  c=dim(mul)[1]
  
  x=as.matrix(x)

  em=em_reg_GM(x,mul,sgml,taul)
  mul=em[[1]]
  sgml=em[[2]]

  res=integer(n)
  dsty=double(c)
  for (i in 1:n){
    for(j in 1:c){
      dsty[j]=dmvnorm(x[i,],mul[j,],sgml[[j]])
      belong=which.max(dsty)
      res[i]=belong
    }
  }
  list(cbind(x,res),mul,sgml)


}
