####################
## Date: 11/25/2013
## Author: Kai Yu
## Function: For Mixture of Gaussian model, return density value of a give point
####################


function=GM_score(x,mul,sgml,taul){
  x=as.vector(x)
  c=dim(mul)[1]
  mahalanobisDistance=vector()
  for (j in c(1:c)){
    mahalanobisDistance[j]=mahalanobis(x,center=mul[j,],cov=sgml[[j]])
  }
  score=t(taul)%*%mahalanobisDistance
  return(score)
}