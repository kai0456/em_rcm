source("src/em_rcm.r")
source("src/em_reg.r")
source("src/em_kotz.r")
source('src/GM_score.r')

library(pROC)

# Compute AUC for outlier detectors given an outlier label 
# 1. compute density of ith observation w.r.t each component j
# 2. find the maximum of density (md) of ith observation among all j
# 3. roc(label,md)
###########
###########
auc_em_rcm=function(x,mul,sgml=0,taul=0,label){
	
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
  
  x=as.matrix(x)

    em=em_rcm_GM(x,mul,sgml,taul)
    mul=em[[1]]
    sgml=em[[2]]
    taul=em[[3]]
    md=sapply(1:n,ith_mdsty,x=x,mul=mul,sgml=sgml,c=c)
    roc_rcm = roc(label, md)
    return(roc_rcm$auc)   
}
###############
###############

###############
###############
auc_em_reg=function(x,mul,sgml=0,taul=0,label){
	
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
  
  x=as.matrix(x)

    em=em_reg_GM(x,mul,sgml,taul)
    mul=em[[1]]
    sgml=em[[2]]
    taul=em[[3]]
    md=sapply(1:n,ith_mdsty,x=x,mul=mul,sgml=sgml,c=c)
    roc_reg = roc(label, md)
    return(roc_reg$auc)   
}
#################################
##################


auc_em_kotz=function(x,mul,sgml=0,taul=0,label){
	
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
  
  x=as.matrix(x)

    em=em_kotz_GM(x,mul,sgml,taul)
    mul=em[[1]]
    sgml=em[[2]]
    taul=em[[3]]
    md=sapply(1:n,ith_mdsty,x=x,mul=mul,sgml=sgml,c=c)
    roc_kotz = roc(label, md)
    return(roc_kotz$auc)   
}
###############


#######
##  density of x w.r.t jth component
#######
jth_dsty=function(j,x,mul,sgml){
  dsty=dmvnorm(x,mul[j,],sgml[[j]])
  # dsty=mahalanobis(x,mul[j,],sgml[[j]])
}
#######
## End density of x wrt jth component
#######

#######
## maximum density of ith observation
#######
ith_mdsty=function(i,x,mul,sgml,c,eps){
  dsty=sapply(1:c,jth_dsty,x=x[i,],mul=mul,sgml=sgml)
  mdsty =max(dsty)
  return(mdsty)
}

########
## End maximum density of ith observation
########












