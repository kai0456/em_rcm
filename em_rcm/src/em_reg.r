
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

mix_reg=function(x,mul,sgml,taul){
  #mul is c*d matrix
  #sgml is list
  #taul is vector
  n=dim(x)[1]
  d=dim(x)[2]
  c=dim(mul)[1]   #number of components
  if (sum(taul)!= 1)
  {print("taul not sum up 1")}
  if (c!=length(sgml)| c!=length(taul)){
    print("Number of component doesn't match!")
    break
  }
  Tij=fTij(x,mul,sgml,taul)

  for (j in 1:c){    # update mean and covariance
    s=sum(Tij[,j])
    mul[j,]=t(x)%*%Tij[,j]/s
    mulj=as.vector(mul[j,])

    x_mius_mu_dn=t(x)-mulj
    sgml[[j]]=x_mius_mu_dn%*%diag(Tij[,j]/s)%*%t(x_mius_mu_dn)+0.00000001*diag(1,d)
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

em_reg_GM=function(x,mul,sgml=0,taul=0){
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
  mull=list(mul)   # 3 levels
  sgmll=list(sgml) # 3 levels
  taull=list(taul) # 2 levels
  for (t in 1:100){
    mix=mix_reg(x,mull[[t]],sgmll[[t]],taull[[t]])
    mull=c(mull,list(mix[[1]]))
    sgmll=c(sgmll,list(mix[[2]]))
    taull=c(taull,list(mix[[3]]))
    if (abs(taull[[t+1]][1]-taull[[t]][1])<0.0001)
      {break}
  }
  ## (t+1)=ending position
  list(mull[[t+1]],sgmll[[t+1]],taull[[t+1]],t+1,mull,sgmll,taull)

}