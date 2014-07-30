################
##plot em_rcm and em_reg on the same graph without outlier.
############

setwd("/home/ubuntu/My_Paper_Github/em_rcm")
source("src/em_rcm.r")
source("src/em_reg.r") 
source("usps/dimension_reduction.r")
library('ElemStatLearn')


x=list()
for (i in 1:7291){
	x[[i]]=dimension_reduction(zip.train[i,])
}

zip.train.reduced=matrix(unlist(x),nrow=7291,byrow=T)

train=zip.train.reduced[,c(2:65)]

fuzzy=kmeans(train,10)
mul=as.matrix(fuzzy$centers)


proc.time()
  
output=em_rcm_GM(train,mul)

proc.time()
