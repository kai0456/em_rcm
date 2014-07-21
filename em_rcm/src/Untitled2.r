fish=read.csv("Fish/allfishes.txt",header=F)

fish=fish[,c(1:4,7,10,12,16:20,30)]

k=c(128,297,172,42,36,53,39,60,76,86)

train_list=list()
test_list=list()
type1err=type2err=vector()
##type1err= Cond. Prob.(Not Outlier | Outlier)
##type2err= Cond. Prob.(Outlier | Not Outlier)


pdf(file="Fish/Rplot.pdf")
par(mfrow=c(2,2))
idx=seq(1:10)
for (i in 1:10){
  train_mat=subset(fish,V30!=i)
  test_mat=subset(fish,V30==i)

  train=lapply(idx[-i],function(j,mtx) subset(mtx,V30==j)[,-13],mtx=fish)
  ## train is a list of 9, i.e. train[[1]]=fish with class id=1 but without class
  ## id shown in last column by ([,-13])

  train_mat=train_mat[,-13]
  test_mat=test_mat[,-13]
  train_mat=as.matrix(train_mat)
  test_mat=as.matrix(test_mat)

  boxplot(train_mat,main=paste("test set=",i,sep=""))

  
}
 dev.off()