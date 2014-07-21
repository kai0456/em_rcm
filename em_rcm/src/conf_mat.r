#######
##Confusion matrix
#######
# final return by comparing at diagonal matrix

cm <- function (actual, predicted){
  actual=as.vector(actual)
  predicted=as.vector(predicted)
	# Produce a confusion matrix

	t<-table(predicted,actual)
	# there is a potential bug in here if columns are tied for ordering
	t[apply(t,2,function(c) which.max(c)),]
}