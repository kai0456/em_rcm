source("src/em_rcm.r")
source("src/em_reg.r")
source('src/em_kotz.r')

# Contain outlier detector with density tolerance eps=.01
###########
###########
cluster_em_rcm=function(x,mul,sgml=0,taul=0,outlier=TRUE,eps=0.01){
  n=dim(x)[1]
  d=dim(x)[2]
  c=dim(mul)[1]
  qtchi=(qchisq(1-eps,d))
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
  

  x=as.matrix(x)
  if (outlier){
    em=em_rcm_GM(x,mul,sgml,taul)
    mul=em[[1]]
    sgml=em[[2]]
    taul=em[[3]]
    res=sapply(1:n,ith_belong,x=x,mul=mul,sgml=sgml,c=c,eps=qtchi)
    return(list(as.data.frame(cbind(x,res)),mul,sgml,taul))
   } else{
    qtchi=Inf
    res=sapply(1:n,ith_belong,x=x,mul=mul,sgml=sgml,c=c,eps=qtchi)
    return(list(as.data.frame(cbind(x,res)),mul,sgml,taul))
   }
}
###############
###############

###############
###############
cluster_em_reg=function(x,mul,sgml=0,taul=0,outlier=TRUE,eps=0.01){
  n=dim(x)[1]
  d=dim(x)[2]
  c=dim(mul)[1]
  qtchi=(qchisq(1-eps,d))
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
  

  x=as.matrix(x)
  if (outlier){
    em=em_reg_GM(x,mul,sgml,taul)
    mul=em[[1]]
    sgml=em[[2]]
    taul=em[[3]]
    res=sapply(1:n,ith_belong,x=x,mul=mul,sgml=sgml,c=c,eps=qtchi)
    return(list(as.data.frame(cbind(x,res)),mul,sgml,taul))
   } else{
    qtchi=Inf
    res=sapply(1:n,ith_belong,x=x,mul=mul,sgml=sgml,c=c,eps=qtchi)
    return(list(as.data.frame(cbind(x,res)),mul,sgml,taul))
   }
}

##################
##################


###########
###########
cluster_em_kotz=function(x,mul,sgml=0,taul=0,outlier=TRUE,eps=0.01){
  n=dim(x)[1]
  d=dim(x)[2]
  c=dim(mul)[1]
  qtchi=(qchisq(1-eps,d))
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
  

  x=as.matrix(x)
  if (outlier){
    em=em_kotz_GM(x,mul,sgml,taul)
    mul=em[[1]]
    sgml=em[[2]]
    taul=em[[3]]
    res=sapply(1:n,ith_belong,x=x,mul=mul,sgml=sgml,c=c,eps=qtchi)
    return(list(as.data.frame(cbind(x,res)),mul,sgml,taul))
   } else{
    qtchi=Inf
    res=sapply(1:n,ith_belong,x=x,mul=mul,sgml=sgml,c=c,eps=qtchi)
    return(list(as.data.frame(cbind(x,res)),mul,sgml,taul))
   }
}
###############
###############




#######
## jth_row vector density
#######
jth_dsty=function(j,x,mul,sgml){
  #dsty=dmvnorm(x,mul[j,],sgml[[j]])
   dsty=mahalanobis(x,mul[j,],sgml[[j]])
}
#######
## End jth_row vector density
#######

#######
## belongs ith component out of c components
#######
ith_belong=function(i,x,mul,sgml,c,eps){
  dsty=sapply(1:c,jth_dsty,x=x[i,],mul=mul,sgml=sgml)
  if (min(dsty)>eps) {belong="outlier"} else {belong=which.min(dsty)}
  return(belong)
}

########
## End belongs ith component out of c components
########











