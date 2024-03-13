# library(MCMCglmm)
## what do we output - how do we incorporate alive individuals that don't breed or multiple measurements
## output pedigree and data structure?

# rm(list=ls())
# source("/Users/joelpick/github/squidPed/R/functions.R")
# 	years = 5
# 	n_females = 100
# 	p_breed = 0.7
# 	fecundity = 4
# 	fixed_fecundity = FALSE
# 	p_sire = 1
# 	p_retain = 0.8
# 	polgyny_rate = 0
# 	juv_surv = 0.25
# 	adult_surv = 0.5
# 	det = list(rates=c("fecundity"), G=1)
# 	immigration = 0
# 	constant_pop = TRUE 
# 	known_age_structure = FALSE
#   verbose=TRUE



## first go - make genetic variation in JS and F

simulate_pedigree_var <- function(
	years = 5,
	n_females = 50,
	
	p_breed = 1, #probability that a female breeds
	fecundity = 4, #mean number of zygotes a female produces each year
	fixed_fecundity = FALSE,
	p_sire = 1, # probability that 'social' male sires all offspring 
	p_retain = 0, #probability that social partnership is retained
	# polgyny_rate = 0,
	juv_surv = 0.25,
	adult_surv = 0.5,

	det = NULL,#list(rates=c("p_breed","fecundity"), G, PE, MG, ME, Y, C, E), ## maybe also have year and cohort effects?
	
	immigration = 0,
	constant_pop = TRUE ,
	known_age_structure = FALSE,
	verbose=FALSE 
){

	options(stringsAsFactors=FALSE)

	## get environment, for transform det function
	Renv <- environment()
	# make all covariance matrices
	transform_det(det,env=Renv)

	P <- G + PE + MG + ME + C + Y + E

## latent scale mean fecundity
	fecundity_l <- log(fecundity) - (P["fecundity","fecundity"])/2

## latent scale mean juvenile survival
## if using threshold, presumably can use pnorm/qnorm to get between observed and latent scale?
  adult_surv_l <- stats::qnorm(adult_surv, 0, sqrt(P["adult_surv","adult_surv"]+1))

  juv_surv_l <- stats::qnorm(juv_surv, 0, sqrt(P["juv_surv","juv_surv"]+1))

  ## random effect means
	mus <- rep(0,ncol(G))

	## year and cohort effects
  year_effect <- MASS::mvrnorm(n=years, mu=mus, Y)
  cohort_effect <- MASS::mvrnorm(n=years, mu=mus, C)
  ## should cohort effects for previous years be simulated?

### GROWTH RATE
	det_growth_rate <- (juv_surv * fecundity)/2 + adult_surv + immigration 

	v_as <- adult_surv * (1-adult_surv)
	v_js <- juv_surv * (1-juv_surv)
	v_fec <- fecundity
	v_rec <- fecundity^2*v_js + juv_surv^2*v_fec + v_js*v_fec

	stoch_growth_rate <- det_growth_rate - (v_as + v_rec)/(2*n_females)

	# immigration <- 1- stoch_growth_rate

	cat("deterministic growth rate =", det_growth_rate,"\n")
	cat("stochastic growth rate =", stoch_growth_rate,"\n")


### BASE POPULATION
	# make pedigree for base population
	pedigree <- data.frame(
		animal = paste0("0_",1:(n_females*2)),
		dam = NA,
		sire = NA,
		sex = rep(c("F","M"),each=n_females),
		cohort=NA
	)

	# simulate breeding values (bv), permanent environment effects (pe), maternal genetic effects (mg) and maternal environment effects (me) for base population
	bv <- MASS::mvrnorm(n=n_females*2, mu=mus, G)
  
  ## first year also don't have known mothers, so MG and ME get added into their PE
	pe <- MASS::mvrnorm(n=n_females*2, mu=mus, PE + MG + ME)
  
  # these effects are the effect that the focal individual has on its offspring
	mg <- MASS::mvrnorm(n=n_females*2, mu=mus, MG)
	me <- MASS::mvrnorm(n=n_females*2, mu=mus, ME)

	rownames(bv)<-rownames(pe)<-rownames(mg)<-rownames(me)<-pedigree[,"animal"]


	# make list that stores who is alive in each year
	dat <- list()
	dat[[1]] <-  data.frame(animal = pedigree$animal, sex = pedigree$sex, age=NA, year=1)
	## think about using rgeom here = should it be different for year 1
	## do we start with age unknown as that would be realistic to the sampling?
	# plot(table(rgeom(1000,0.5)+1))
	# dgeom(0:10,0.5)

# 	e <- MASS::mvrnorm(n=nrow(dat[[1]]), mu=0, E)

	if(verbose) cat("initial sims done \n")

	pairs <- list()

	# year=1
	for(year in 1:years){

		# probability of breeding just on females? assume that male breeding is dependent on females?

		# get vectors of females and males available to breed (around as adults pre-breeding (dat[[year]])), accounting for the probability of females breeding
		females <- subset(dat[[year]],sex=="F")$animal
		breeding_females <- females[as.logical(rbinom(length(females),1,p_breed))]
		males <- subset(dat[[year]],sex=="M")$animal

		#work out number of pairs that can be formed 
		n_pair <- length(breeding_females)#min(length(breeding_females),length(males))
		# print(length(males))
		# print(length(breeding_females))


		### assign 'social' male
		if(year==1 || adult_surv==0){ # or adult_surv=0?
			## should this be with replacement?
			social_male<-sample(males,n_pair, replace=TRUE)	#TURE?
		}else{
			# if female existed last year, who was male. else sample new male
			# if male alive, then rbinom(1,1,p_retain), 
			# if retain then male, else sample new male
			social_male<-sapply(breeding_females,function(bf){
				if(bf %in% pairs[[year-1]]$female 
				& pairs[[year-1]]$male[match(bf,pairs[[year-1]]$female)] %in% males
				& stats::rbinom(1,1,p_retain)==1 ){
				  pairs[[year-1]]$male[match(bf,pairs[[year-1]]$female)]
				}else{
					sample(males,1, replace=TRUE)
				}
			})
		}
		
		pairs[[year]] <- data.frame(
			female = breeding_females,
			male = social_male
			)		

##### FECUNDITY

		#breeding female mother IDs
		bf_dams <- pedigree$dam[match(breeding_females,pedigree$animal)]
		#breeding female cohorts
		bf_cohort <- pedigree$cohort[match(breeding_females,pedigree$animal)]
		
		# residual variance
		if(year==1){
			## first years e = c+e, as dont know what c is (unless implement known ages)
			f_e <- stats::rnorm(n_pair, 0, sqrt(C["fecundity","fecundity"] + E["fecundity","fecundity"]))
			# as_e <- rnorm(n_pair, 0, sqrt(C["adult_surv","adult_surv"] + E["adult_surv","adult_surv"]))
			# for males and femaels
		}else{
			f_e <- stats::rnorm(n_pair, 0, sqrt(E["fecundity","fecundity"]))	
			# as_e <- rnorm(n_pair, 0, sqrt(E["adult_surv","adult_surv"]))	
		}		
		## save E in dat?

		## multivariate E? residual covariance in survival and fecundity?

		f_exp <- exp(fecundity_l + 
			Matrix::rowSums(cbind(
				bv[breeding_females, "fecundity"],
				pe[breeding_females, "fecundity"],
				### with overlapping gens some of these might be NA?
				if(year!=1) mg[match(bf_dams,rownames(mg)), "fecundity"],
				if(year!=1) me[match(bf_dams,rownames(mg)), "fecundity"],
				year_effect[year, "fecundity"],
				if(year!=1) cohort_effect[bf_cohort, "fecundity"]
				),na.rm=TRUE)
			+ f_e)

		# number of offspring per female
		n_juv <- if(fixed_fecundity) {
			rep(fecundity,n_pair)
		}else{
			stats::rpois(n_pair,f_exp)
			# rcomp(n_pair, mu = f_exp, nu = 2)

		}

	if(verbose) cat("year",year,"fecundity \n")

### Generate offspring, pedigree and genetic and permanent environmental values

		## make ped incorporating EPP and fecundity
		ped <- data.frame(
			animal=paste0(year,"_",1:sum(n_juv)),
			# fecundity
			dam=rep(breeding_females,n_juv),
			# EPP
			sire=c(lapply(1:n_pair,function(i){
					n_sired <- stats::rbinom(1,n_juv[i],p_sire)
					c(rep(social_male[i], n_sired), sample(males,n_juv[i]-n_sired,replace=TRUE))
				}), recursive=TRUE),
			#equal sex ratio
			sex=sample(c("M","F"),sum(n_juv),replace=TRUE),
			cohort=year)

		## calculate genetic values from parents plus Mendelian sampling deviation
		bv_add <- (bv[ped$dam,] + bv[ped$sire,])/2 + MASS::mvrnorm(n=nrow(ped), mu=mus, G/2)
		mg_add <- (mg[ped$dam,] + mg[ped$sire,])/2 + MASS::mvrnorm(n=nrow(ped), mu=mus, MG/2)

		## generate environment effects
		pe_add <- MASS::mvrnorm(n=nrow(ped), mu=mus, PE)
		me_add <- MASS::mvrnorm(n=nrow(ped), mu=mus, ME)

		rownames(bv_add)<-rownames(pe_add)<-rownames(mg_add)<-rownames(me_add)<-ped[,"animal"]

		pedigree <- rbind(pedigree,ped)
		bv <- rbind(bv,bv_add)
		pe <- rbind(pe,pe_add)
		mg <- rbind(mg,mg_add)
		me <- rbind(me,me_add)

	if(verbose) cat("year",year,"ped \n")

		# print(table(aggregate(animal~dam+cohort,pedigree, function(x)length(x))[,3]))
		# print(table(aggregate(animal~dam+cohort,pedigree, function(x)length(x))[,2]))

##### create individuals present in the next year


	##### Juvenile Survival
	
		# residual variance
		# same as PE
		# js_e <- rnorm(n_pair, 0, E["juv_surv","juv_surv"])	
			
		js_exp <- stats::pnorm(juv_surv_l + 
			Matrix::rowSums(cbind(
				bv[ped$animal, "juv_surv"]
				, pe[ped$animal, "juv_surv"]
				, mg[ped$dam, "juv_surv"]
				, me[ped$dam, "juv_surv"]
				, year_effect[year, "juv_surv"]
				# , cohort_effect[year, "juv_surv"]
				## there could be cohort covariance - juvenile survival and adult traits - bad conditions affect later life traits
				## there could be year covariance - years where 
			),na.rm=TRUE)
			# + js_e
			## don't think can have residual covariance - doesn't link to anything 
		)

	## check the conversion works 

		#breeding female mother IDs
		f_dams <- pedigree$dam[match(females,pedigree$animal)]
		m_dams <- pedigree$sire[match(males,pedigree$animal)]
		#breeding female cohorts
		f_cohort <- pedigree$cohort[match(females,pedigree$animal)]
		m_cohort <- pedigree$cohort[match(males,pedigree$animal)]
		

#females[as.logical(rbinom(length(females),1,adult_surv))]
		as_exp_f <- stats::pnorm(adult_surv_l + 
			Matrix::rowSums(cbind(
				bv[females, "adult_surv"],
				pe[females, "adult_surv"],
				if(year!=1) mg[match(f_dams,rownames(mg)), "adult_surv"],
				if(year!=1) me[match(f_dams,rownames(me)), "adult_surv"],
				year_effect[year, "adult_surv"],
				if(year!=1) cohort_effect[f_cohort, "adult_surv"]
			),na.rm=TRUE)
			# + as_e
		)
		as_exp_m <- stats::pnorm(adult_surv_l + 
			Matrix::rowSums(cbind(
				bv[males, "adult_surv"],
				pe[males, "adult_surv"],
				if(year!=1) mg[match(m_dams,rownames(mg)), "adult_surv"],
				if(year!=1) me[match(m_dams,rownames(me)), "adult_surv"],
				year_effect[year, "adult_surv"],
				if(year!=1) cohort_effect[m_cohort, "adult_surv"]
			),na.rm=TRUE)
			# + as_e
		)

		next_year_ind <- if(constant_pop){
			rbind(
				### need to ensure equal sex ratio of recruits, otherwise population size fluctuations
				# ped[sample(which(ped[,"sex"]=="F"), juv_surv*fecundity*n_females/2, replace=FALSE),c(1,4)],
				# ped[sample(which(ped[,"sex"]=="M"), juv_surv*fecundity*n_females/2, replace=FALSE),c(1,4)],
				ped[sample(which(ped[,"sex"]=="F"), juv_surv*fecundity*n_females/2, replace=FALSE, prob=js_exp[which(ped[,"sex"]=="F")]),c(1,4)],
				### i think this sampling is right - needs checking
				ped[sample(which(ped[,"sex"]=="M"), juv_surv*fecundity*n_females/2, replace=FALSE, prob=js_exp[which(ped[,"sex"]=="M")]),c(1,4)],
				if(adult_surv>0){	
					cbind(animal=sample(females, adult_surv*n_females, replace=FALSE), sex="F")},#, prob=as_exp_f
				if(adult_surv>0){	
					cbind(animal=sample(males, adult_surv*n_females, replace=FALSE), sex="M")}#, prob=as_exp_m
			)
		}else{
			rbind(
				## variation in JS here
				# ped[as.logical(rbinom(nrow(ped),1,juv_surv)),c(1,4)],
				ped[as.logical(rbinom(nrow(ped),1,js_exp)),c(1,4)],
				
				## variation in AS here
				if(adult_surv>0){	
				cbind(animal=females[as.logical(rbinom(length(females),1,as_exp_f))], sex="F")},
				if(adult_surv>0){	
				cbind(animal=males[as.logical(rbinom(length(males),1,as_exp_m))], sex="M")}
			)
		}
	  next_year_ind$age <- (year+1) - pedigree[match(next_year_ind$animal,pedigree$animal),"cohort"]
 	  next_year_ind$year <- year+1
	  dat[[year+1]] <- next_year_ind
	  ### save fecundity, survival etc into dat, so makes analysis easier
	if(verbose) cat("year",year,"next year \n")

		### update dat to include sruvvial and reproduction in the current year 
				dat[[year]]$survival <- as.numeric(dat[[year]]$animal %in% next_year_ind$animal)
			dat[[year]]$fecundity <- 0
			# dat[[year]]$recruits <- 0

			## female fecundity
			index_f<-which(dat[[year]]$animal%in%breeding_females)
			dat[[year]]$fecundity[index_f] <- n_juv[match(dat[[year]]$animal[index_f],breeding_females)]

			##male fecundity
			male_fecundity<-table(ped$sire)
			index_m <- which(dat[[year]]$animal%in%names(male_fecundity))
			dat[[year]]$fecundity[index_m] <- male_fecundity[match(dat[[year]]$animal[index_m],names(male_fecundity))]

			
			## if in last year, the next year in wont have any survival or reproudction
		if(year ==years){
			dat[[year+1]]$survival <- NA
			dat[[year+1]]$fecundity <- NA
		}

			#dam_fec <- table(ped$dam)
		

	}
	

	return(list(pedigree=pedigree,data_str=do.call(rbind,dat), y=year_effect, c=cohort_effect, a=bv,pe=pe,mg=mg,me=me))

}
