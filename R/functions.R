
exp2lat<- function(mean,cov){
	mean <- as.matrix(mean)
	cov <- as.matrix(cov)
	
	mean_out <- rep(NA,2)
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
	
	mean_out <- rep(NA,2)
	for(i in 1:nrow(cov)) mean_out[i] <- exp(mean[i]+ cov[i,i]/2)

	cov_out <- matrix(NA,nrow(cov),ncol(cov))
	for(i in 1:nrow(cov)){
		for(j in 1:ncol(cov)) cov_out[i,j] <- exp(mean[i]+mean[j] + (cov[i,i]+cov[j,j])/2)*(exp(cov[i,j])-1)
	}
	return(list(mean=mean_out,cov=cov_out))
}


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



transform_det <- function(det, env){
	# print(environment())
	# variance components

	components <- c("G", "PE", "MG", "ME", "C", "Y", "E")
	# vital rates that can vary 

	rates <- c("p_breed","fecundity","p_sire","p_retain","juv_surv","adult_surv")
	empty_mat <- matrix(0,ncol=length(rates),nrow=length(rates),dimnames=list(rates,rates))
# G <- PE <- MG <- ME <- C <- Y <- E <- empty_mat


	for(component in components){
		mat <- det[[component]]
		if(!is.null(mat)){
			if(!is.matrix(mat)) mat <- as.matrix(mat)
			rownames(mat)<-colnames(mat)<-det[["rates"]]
		}
	
		x <- empty_mat
		x[rownames(mat),colnames(mat)] <- mat
		assign(component,x,pos=env)
	}
}


