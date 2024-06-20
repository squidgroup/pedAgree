
## function to check parameters input to simulate_pedigree
general_check <- function(name,env, rate=TRUE, sex_specific=TRUE){
	x <- get(name,pos=env)

	if(!is.vector(x)) stop(paste(name,"is not a vector"))
	
	if(rate){
		if(!all(between(x,0,1))) stop(paste(name,"must be between 0 and 1"))	
	}else{
		if(!all(is.wholenumber(x) & x>0)) stop(paste(name,"must be a whole number and >0"))	
	}
	
	
	if(sex_specific){
		if(!length(x) %in% 1:2) stop(paste(name,"should be length 1 or 2"))
	
		assign(paste0(name,"_f"),x[1],pos=env)
		if(length(x)==2){	
			assign(paste0(name,"_m"),x[2],pos=env)
		}else{
			assign(paste0(name,"_m"),x[1],pos=env)
		}
	}else{
		if(!length(x) == 1) stop(paste(name,"should be length 1"))
	}
}




## function to fill out deterministic parameters

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


### function that pairs males and females
### priority for pairing is given to unpaired males, if specified
sample_males <- function(females,males,unpaired_males=NULL){

	if(is.null(unpaired_males)) unpaired_males <- males
	# if unpaired_males are not specified, assumes all males are unpaired

	## when more breeding females, then all males are paired, and then sampled for extra ones
	if(length(females)>=length(unpaired_males)){
		sample(c(unpaired_males, sample(males,length(females)-length(unpaired_males), replace = TRUE)), replace=FALSE)
	}else{
		sample(unpaired_males, length(females), replace=FALSE)
	}
}



#	Function to sampling how many of the offspring are sired by different males. It takes the probability p, and works out how many the paired male sired, and then gives the same probability to subsequent males. This means that extra pair males will be few, and have several offspring if p_sire if high
fill_sires <- function(n,p){
	n_remaining<-n
	sire <- 1
	sires <- rep(NA,length=n)
	while(n_remaining>0){
		n_sired <- (stats::rbinom(1,n_remaining,p))
		if(n_sired>0){	 
			sires[(n-sum(is.na(sires)) +1):(n-sum(is.na(sires)) + n_sired)] <- sire
			n_remaining <- sum(is.na(sires))
				 sire <- sire+1}
	}
	sires
}


