################
##plot em_rcm and em_reg on the same graph without outlier.
############

setwd("/home/ubuntu/My_Paper/em_rcm")
source("src/em_rcm.r")
source("src/em_reg.r") 
library('ElemStatLearn')

train=zip.train[,c(2:257)]

fuzzy=kmeans(train,10)
mul=as.matrix(fuzzy$centers)

train_small=train[c(1:1024),]

proc.time()
  
output=em_rcm_GM(train_small,mul)

proc.time()

## Test Github branch feature1