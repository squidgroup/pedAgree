# library(MCMCglmm)
## what do we output - how do we incorporate alive individuals that don't breed or multiple measurements
## output pedigree and data structure?

	# years = 5
	# n_females = 10
	# p_breed = 1
	# fecundity = 4
	# p_sire = 1
	# p_retain = 0.8
	# polgyny_rate = 0
	# juv_surv = 0.5
	# adult_surv = 0
	# immigration = 0
	# constant_pop = TRUE 
	# known_age_structure = FALSE

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


simulate_pedigree <- function(
	years = 5,
	n_females = 50,
	
	p_breed = 1, #probability that a female breeds
	fecundity = 4, #number of zygotes a female produces each year
	p_sire = 1, # probability that 'social' male sires all offspring 
	p_retain = 0, #probability that social partnership is retained
	# polgyny_rate = 0,
	juv_surv = 0.25,
	adult_surv = 0.5,

	det = NULL,#list(components=c("p_breed","fecundity"), G, PE, MG, ME, E),
	
	immigration = 0,
	constant_pop = TRUE ,
	known_age_structure = FALSE){

options(stringsAsFactors=FALSE)

	if(missing(det)) det <- list(components=c("p_breed","fecundity","p_sire","p_retain","juv_surv","adult_surv"), G, PE, MG, ME, E)

	fecundity_l <- log(fecundity) - (G["fecundity","fecundity"]+PE["fecundity","fecundity"]+E["fecundity","fecundity"])/2

## if using threshold, presumably can use pnorm/qnorm to get between observed and latent scale?
  juv_surv_l <- qnorm(juv_surv,0,G["juv_surv","juv_surv"]+PE["juv_surv","juv_surv"]+E["juv_surv","juv_surv"])


	# det_growth_rate <- (juv_surv * fecundity)/2 + adult_surv + immigration 

	# v_as <- adult_surv * (1-adult_surv)
	# v_js <- juv_surv * (1-juv_surv)
	# v_fec <- fecundity
	# v_rec <- fecundity^2*v_js + juv_surv^2*v_fec + v_js*v_fec

	# stoch_growth_rate <- det_growth_rate - (v_as + v_rec)/(2*n_females)

	# immigration <- 1- stoch_growth_rate


	# make pedigree for base population
	pedigree <- data.frame(
		animal = paste0("0_",1:(n_females*2)),
		dam = NA,
		sire = NA,
		sex = rep(c("F","M"),each=n_females),
		cohort=NA
	)

		### store breeding values in ped - one record per individual
	bv <- MASS::mvrnorm(n=n_females*2, mu=0, G)
	pe <- MASS::mvrnorm(n=n_females*2, mu=0, PE)

	pedigree$f_bv <- bv[,1]
	pedigree$js_bv <- bv[,2]
	pedigree$as_bv <- bv[,3]


	# make list that stores who is alive in each year
	dat <- list()
	dat[[1]] <-  data.frame(animal = pedigree$animal, sex = pedigree$sex, age=NA)
	## think about using rgeom here = should it be different for year 1
	## do we start with age unknown as that would be realistic to the sampling?
	# plot(table(rgeom(1000,0.5)+1))
	# dgeom(0:10,0.5)

# 	e <- MASS::mvrnorm(n=nrow(dat[[1]]), mu=0, E)


	pairs <- list()

	# year=1
	for(year in 1:years){
		## import individuals that are around as adults pre-breeding (dat[[year]])

		# probability of breeding just on females? assume that male breeding is dependent on females?

		# get vectors of females and males available to breed, accounting for the probability of females breeding
		females <- subset(dat[[year]],sex=="F")$animal
		breeding_females <- females[as.logical(rbinom(length(females),1,p_breed))]
		males <- subset(dat[[year]],sex=="M")$animal

		#work out number of pairs that can be formed 
		n_pair <- length(breeding_females)#min(length(breeding_females),length(males))
		# print(length(males))
		# print(length(breeding_females))



		
		### assign 'social' male
		if(year==1){
			## should this be with replacement?
			social_male<-sample(males,n_pair, replace=TRUE)	
		}else{
			social_male<-sapply(breeding_females,function(bf){
				if(bf %in% pairs[[year-1]]$female 
				& pairs[[year-1]]$male[match(bf,pairs[[year-1]]$female)] %in% males
				& rbinom(1,1,p_retain)==1 ){
				  pairs[[year-1]]$male[match(bf,pairs[[year-1]]$female)]
				}else{
					sample(males,1, replace=TRUE)
				}
			})
		}
		
		
# if female existed last year, who was male. else sample new male
		# if male alive, then rbinom(1,1,p_retain), 
		  # if retain then male, 
		  # else sample new male
		pairs[[year]] <- data.frame(
			female = breeding_females,
			male = social_male
			)		

		# number of offspring per female
		# n_juv <- rpois(n_pair,fecundity)
		n_juv <- rep(fecundity,n_pair)

		# exp(rnorm(n_pair,fecundity_l, )
	f_index <- match(breeding_females,pedigree$animal)
	f_m_index <- match(breeding_females,pedigree$dam)
	f_e <- rnorm(n_pair, 0, E[1,1])

f_exp <- exp(fecundity_l + 
	rowSums(cbind(pedigree[f_index,c("f_bv","f_pe"],pedigree[f_index,c("f_mg","f_me"]))
+ f_e)

		## make ped incorporating EPP and fecundity
		ped <- data.frame(
			animal=paste0(year,"_",1:sum(n_juv)),
			# fecundity
			dam=rep(breeding_females,n_juv),
			# EPP
			sire=c(lapply(1:n_pair,function(i){
					n_sired <- rbinom(1,n_juv[i],p_sire)
					c(rep(social_male[i], n_sired), sample(males,n_juv[i]-n_sired,replace=TRUE))
				}), recursive=TRUE),
			#equal sex ratio
			sex=sample(c("M","F"),sum(n_juv),replace=TRUE),
			cohort=year)
		pedigree <- rbind(pedigree,ped)

		# print(table(aggregate(animal~dam+cohort,pedigree, function(x)length(x))[,3]))
		# print(table(aggregate(animal~dam+cohort,pedigree, function(x)length(x))[,2]))

		## create individuals present in the next year
		next_year_ind <- if(constant_pop){
			rbind(
				### need to ensure equal sex ratio of recruits, otherwise population size fluctuations
				ped[sample(which(ped[,"sex"]=="F"), juv_surv*fecundity*n_females/2, replace=FALSE),c(1,4)],
				ped[sample(which(ped[,"sex"]=="M"), juv_surv*fecundity*n_females/2, replace=FALSE),c(1,4)],
				if(adult_surv>0){	
					cbind(animal=sample(females, adult_surv*n_females, replace=FALSE), sex="F")},
				if(adult_surv>0){	
					cbind(animal=sample(males, adult_surv*n_females, replace=FALSE), sex="M")}
			)
		}else{
			rbind(
				ped[as.logical(rbinom(nrow(ped),1,juv_surv)),c(1,4)],
				cbind(animal=females[as.logical(rbinom(length(females),1,adult_surv))], sex="F"),
				cbind(animal=males[as.logical(rbinom(length(males),1,adult_surv))], sex="M")
			)
		}
	  next_year_ind$age<-(year+1) - pedigree[match(next_year_ind$animal,pedigree$animal),"cohort"]
	  dat[[year+1]] <- next_year_ind
	  ### save fecundity, survival etc into dat, so makes analysis easier
	}

	return(list(pedigree=pedigree,data_str=do.call(rbind,dat)))

}
