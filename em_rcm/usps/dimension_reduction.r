## Dimension Reduction for usps data

dimension_reduction_matrix_form=function(x){
	## Take the matrix form x to reduce dimension, e.g 256*256 -> 64*64
	x=as.matrix(x)
	if (dim(x)[1] != dim(x)[2]){print('Matrix x has differnt width and length')}
	colNum=dim(x)[1]
	rowNum=colNum/2
	
	## Reduce the dimesion of x by half
	## For example, leftMatrix_{8,16} %*% x %*% rightMutiplier_{16,8}
	leftMutiplier=matrix(0,nrow=rowNum,ncol=colNum)

	for (i in 1:rowNum){
		for (j in (2*i-1):(2*i)){
			leftMutiplier[i,j]=1/4
		}

	}

	rightMutiplier=t(leftMutiplier)*4
	
	x_reduced=leftMutiplier %*% x %*% rightMutiplier

}

##

dimension_reduction=function(x){
	x=as.vector(x)
	digit=x[1]
	print(paste('Digit ', digit, ' is taken to reduce dimension. :)'))
	
	n=length(x)-1
	if (n%%2 != 0){print('vector length is not a multiple of 2')}
	
	vec=x[c(2:length(x))]
	rowNum=colNum=sqrt(length(vec))
	x_matrix_form=matrix(vec,nrow=rowNum,byrow=FALSE)
	x_matrix_form=x_matrix_form[,sort(seq(1:dim(x_matrix_form)[2]),decreasing=T)]
	
	x_matrix_reduced=dimension_reduction_matrix_form(x_matrix_form)
	
	x_matrix_reduced=x_matrix_reduced[,sort(seq(1:dim(x_matrix_reduced)[2]),decreasing=T)]
	
	x_result=c(x[1],as.vector(x_matrix_reduced))
	
	return(x_result)	
}

##
zip4image=function (zip,line) 
{
	if (is.vector(zip)) {
		im <- zip
	} else {
		im <- zip[line, ]
	}
    
    print(paste("digit ", im[1], " taken"))
	rowNum=sqrt(length(zip)-1)
	if (rowNum%%2 != 0) {print(paste('Length of ',zip,' is not a multiple of 2'))}
    im <- im[-1]
    im <- t(matrix(im, rowNum, rowNum, byrow = TRUE))
    im <- im[, rowNum:1]
    return(im)
}


# image(im, col=gray(256:0/256), zlim=c(0,1), xlab="", ylab="" )
