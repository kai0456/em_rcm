library(e1071)

k_means=function(x,k){
  x=as.matrix(x)
  res=kmeans(x,k)
  tss=res$tot.withinss
  for (i in 1:30){  
    res1=kmeans(x,k)
    tss1=res1$tot.withinss
    if (tss<tss1){res=res}else{res=res1}
    tss=res$tot.withinss
  }
  #return(res$centers)
  structure(list(centers=res$centers,cluster=res$cluster))    
}

############## scratch experiment
#a=matrix(0,nrow=10,ncol=8)
#for (i in 1:10){
#  a[i,]=sort(k_means(train_mat,8)[[1]])-sort(k[2:10])
#  
#}