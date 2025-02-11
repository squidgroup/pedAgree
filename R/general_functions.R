exp2lat<- function(mean,cov){
	mean <- as.matrix(mean)
	cov <- as.matrix(cov)
	
	mean_out <- rep(NA,nrow(cov),ncol(cov))
	for(i in 1:nrow(cov)) mean_out[i] <- log(mean[i]^2/sqrt(mean[i]^2+cov[i,i]))
	
	cov_out <- matrix(NA,nrow(cov),ncol(cov))
	for(i in 1:nrow(cov)){
		for(j in 1:ncol(cov)) cov_out[i,j] <- log(1 + cov[i,j]/(mean[i]*mean[j]))
	}
	return(list(mean=mean_out,cov=cov_out))
}

## function to get means and covariance matrix from latent (normal) to expected (lognormal) scale
lat2exp<- function(mean,cov){
	mean <- as.matrix(mean)
	cov <- as.matrix(cov)
	
	mean_out <- rep(NA,nrow(cov),ncol(cov))
	for(i in 1:nrow(cov)) mean_out[i] <- exp(mean[i]+ cov[i,i]/2)

	cov_out <- matrix(NA,nrow(cov),ncol(cov))
	for(i in 1:nrow(cov)){
		for(j in 1:ncol(cov)) cov_out[i,j] <- exp(mean[i]+mean[j] + (cov[i,i]+cov[j,j])/2)*(exp(cov[i,j])-1)
	}
	return(list(mean=mean_out,cov=cov_out))
}



is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

between <- function(x,min,max) x>=min & x<=max


mult_rownames <- function(..., names){

}

## fucntion to grab diagonal element from matrix
quick_diag <- function(matrix,element){
	matrix[element,element]
}

row_col_names <- function(matrix,names){
	rownames(matrix) <- colnames(matrix) <- names
	matrix	
}
